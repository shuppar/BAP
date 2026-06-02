#!/usr/bin/env python
"""
build_yaml.py — generate config/brain.yaml and config/placenta.yaml from sample_metadata.csv

Run from the Analysis/ directory:
    uv run python scripts/build_yaml.py

This script is idempotent — re-run it any time the metadata CSV changes.
Writes two YAML files (one per tissue) with:
  - sample manifest (id, paths, group, age, sex, pool, donor_id, library)
  - sex marker genes for the Y-chromosome inferred sex check (Phase 0)
  - QC thresholds
  - Random seed

Dev-mode subsetting lives in config/dev.yaml (hand-maintained), which pulls
sample records from brain.yaml via `samples_from:` and selects an explicit
subset by `subset.sample_ids`.
"""

from pathlib import Path
import sys

import pandas as pd
import yaml


# -----------------------------------------------------------------------------
# Paths and config — adjust if you move things around
# -----------------------------------------------------------------------------
METADATA_CSV = Path("sample_metadata.csv")
CONFIG_DIR = Path("config")


# -----------------------------------------------------------------------------
# Shared config blocks (same for brain and placenta)
# -----------------------------------------------------------------------------
SHARED_CONFIG = {
    "qc": {
        "pct_mt_max": 1.0,          # snRNA: should be near zero
        "pct_hemo_max": 5.0,        # critical for placenta
        "n_mads": 5,                # MAD-based bounds on n_genes & log-counts
        "min_counts": 500,          # hard UMI floor (debris below this)
        "min_genes": 200,           # hard gene floor (uninformative below this)
    },
    "sex_markers": {
        "y_linked": ["Ddx3y", "Uty", "Eif2s3y", "Kdm5d"],
        "x_linked": ["Xist"],
    },
    "composition": {
        # Phase 8a: min donors (samples) per group to attempt scCODA. Below 3 the
        # Bayesian variance estimate is essentially unconstrained. Lower to 2 to
        # attempt n=2 groups (e.g. P1 Late Stress) — those run flagged
        # 'unreliable_n<3' in composition_results.csv. CLI --min-donors overrides.
        "min_donors": 3,
    },
    "pathways": {
        # Phase 8c: GSEA on DE stats. Gene sets come from a single MSigDB export
        # TSV produced by scripts/fetch_genesets.R (msigdbr, mouse-native — the
        # decoupler get_resource(MSigDB, mouse) path is broken). Generate once:
        #   Rscript scripts/fetch_genesets.R --out refs/msigdb_mouse.tsv
        "geneset_tsv": "refs/msigdb_mouse.tsv",
        # MH=hallmark, M2=Reactome(curated), M5=GO:BP(ontology), M8=cell-type sigs
        # (M8 esp. useful for 7b subclusters: "does this match a known state?").
        "collections": ["MH", "M2", "M5", "M8"],
        "min_genes_per_set": 5,
        "use_builtin_stress_sets": False,   # optional niche supplement (off)
        "run_tf_activity": False,           # CollecTRI TF activity (needs network)
    },
    "random_seed": 42,
}


# -----------------------------------------------------------------------------
# Per-tissue reference / annotation config (Phase 7 + 7c).
#
# Emitted into each tissue YAML so it survives regeneration. Edit the YAML
# directly afterwards if you need to point at local model/atlas files — re-running
# build_yaml.py will reset these to the defaults below, so keep custom paths noted.
#
# celltypist_models: per-age CellTypist model (Phase 7). Brain P1 has a built-in
#   model; adult brain and placenta need a custom .pkl or fall back to markers.
# reference: scANVI label-transfer config (Phase 7c, 07c_label_transfer.py).
#   - ref_h5ad: path to a labeled reference AnnData (e.g. ABC Atlas subset).
#   - labels_key: .obs column in the reference holding cell type labels.
#   - region_key: .obs column holding region labels. Set to null for tissues
#     with NO spatial/regional reference (e.g. placenta) — regional claims are
#     then skipped entirely and only cell-type labels are transferred.
#   - region_concentration_threshold: a cell type earns a regional label only if
#     >= this fraction of its reference cells fall in a single region. Below it,
#     the region label is suppressed (no fuzzy regional claims). 0.8 = 80%.
# -----------------------------------------------------------------------------

REFERENCE_CONFIG = {
    "brain": {
        "celltypist_models": {
            "P1": "Developing_Mouse_Brain.pkl",   # only built-in mouse brain model
            # "4W":  "<custom Allen BCA .pkl>",    # no built-in adult model
            # "3mo": "<custom Allen BCA .pkl>",
        },
        "reference": {
            "ref_h5ad": None,                      # path to ABC Atlas reference subset
            "labels_key": "cell_type",             # cell type column in the reference
            "region_key": "region",                # region column → enables regional claims
            "region_concentration_threshold": 0.8,
        },
    },
    "placenta": {
        # No built-in CellTypist placenta model; supply a custom .pkl if available.
        "celltypist_models": {},
        "reference": {
            "ref_h5ad": None,                      # path to a placenta reference (if any)
            "labels_key": "cell_type",
            "region_key": None,                    # NO spatial reference → cell-type only
            "region_concentration_threshold": 0.8, # ignored when region_key is null
        },
    },
}


def build_sample_entry(row: pd.Series) -> dict:
    """Convert one CSV row into one sample dict for the YAML.

    The sex field is normalized: 'TBD' -> 'unknown' (will be inferred in Phase 0).
    """
    sex = row["sex_declared"]
    if sex == "TBD":
        sex = "unknown"
    return {
        "id": row["sample_id"],
        "donor_id": row["donor_id"],
        "h5": row["h5_path"],
        "raw_h5": row["raw_h5_path"],
        "age": row["age"],
        "group": row["group"],
        "sex": sex,
        "pool": row["pool"],
        "library": row["library"],
    }



# -----------------------------------------------------------------------------
# Declarative contrast spec (Phase 8: 8a composition, 8b DE, 8c pathway,
# 8e communication, 8g cross-age). The downstream engines iterate over this —
# adding a contrast = edit here + re-run build_yaml.py, not a code change.
#
# Verbatim from project doc §6. Brain has all contrasts; placenta keeps only the
# analyzable ones (no complete factorial — E12.5 has Early+Relaxed, E18.5 has
# Late+Relaxed; cross-age placenta is NOT comparable — see project doc §2).
#
# Per contrast:
#   design      : DESeq2-style formula (pseudobulk; pup = statistical unit)
#   group_by    : split data into independent analyses by this obs key
#   test        : the factor tested (or 'group_omnibus' for the F-test)
#   levels      : [test, reference] for pairwise (omit for omnibus)
#   flag        : primary | secondary | confounded_with_pool | underpowered_exploratory | derived
#   confound_warnings : optional per-level notes surfaced in output
# -----------------------------------------------------------------------------

BRAIN_CONTRASTS = {
    "early_vs_relaxed_per_age": {
        "description": "Early Stress vs Relaxed, within each age — primary",
        "design": "~ sex + pool + group", "group_by": "age", "test": "group",
        "levels": ["Early_Stress", "Relaxed"], "flag": "primary",
    },
    "late_vs_relaxed_per_age": {
        "description": "Late Stress vs Relaxed, within each age — primary (P1 pool-confounded)",
        "design": "~ sex + pool + group", "group_by": "age", "test": "group",
        "levels": ["Late_Stress", "Relaxed"], "flag": "primary",
        "confound_warnings": {
            "P1": "Late Stress at P1 is Pool3 only; Relaxed at P1 is Pool2 — pool vs group fully confounded.",
        },
    },
    "omnibus_3group_per_age": {
        "description": "F-test: do the three groups differ at each age?",
        "design": "~ sex + pool + group", "group_by": "age",
        "test": "group_omnibus", "flag": "primary",
    },
    "early_vs_late_per_age": {
        "description": "Early vs Late Stress — secondary, do the two timings differ?",
        "design": "~ sex + pool + group", "group_by": "age", "test": "group",
        "levels": ["Early_Stress", "Late_Stress"], "flag": "secondary",
    },
    "within_group_across_age": {
        "description": "Developmental trajectory within each group — pool-confounded with age",
        "design": "~ sex + age", "group_by": "group", "test": "age",
        "pairwise": [["P1", "4W"], ["4W", "3mo"], ["P1", "3mo"]],
        "flag": "confounded_with_pool",
    },
    "within_age_sex_stratified": {
        "description": "Stress contrasts within each sex × age — exploratory (n=2)",
        "design": "~ pool + group", "group_by": ["age", "sex"], "test": "group",
        "levels": ["Early_Stress", "Relaxed"], "flag": "underpowered_exploratory",
    },
    "group_x_age_interaction": {
        "description": "Does the stress effect change with age?",
        "design": "~ sex + pool + group * age", "test": "group:age",
        "flag": "underpowered_exploratory",
    },
    # Post-hoc set operations (run by 08g_cross_age.py, not a DE engine)
    "persistence_early": {
        "description": "Persistent/early-only/emergent/transient of Early-vs-Relaxed DEGs",
        "source_contrast": "early_vs_relaxed_per_age",
        "ages_required": ["P1", "4W", "3mo"], "flag": "derived",
    },
    "persistence_late": {
        "description": "Same for Late-vs-Relaxed DEGs (P1 carries confound flag)",
        "source_contrast": "late_vs_relaxed_per_age",
        "ages_required": ["P1", "4W", "3mo"], "flag": "derived",
    },
}

# Placenta: only analyzable contrasts (incomplete factorial — project doc §2).
# E12.5 = Early+Relaxed only; E18.5 = Late+Relaxed only; no cross-age.
PLACENTA_CONTRASTS = {
    "early_vs_relaxed_E12.5": {
        "description": "E12.5 placenta Early vs Relaxed — primary (mostly Pool3; 2 Relaxed in Pool4)",
        "design": "~ pool + group", "group_by": "age", "test": "group",
        "levels": ["Early_Stress", "Relaxed"], "flag": "primary",
        "subset": {"age": "E12.5"},
    },
    "late_vs_relaxed_E18.5": {
        "description": "E18.5 placenta Late vs Relaxed — primary (all Pool4, clean)",
        "design": "~ sex + group", "group_by": "age", "test": "group",
        "levels": ["Late_Stress", "Relaxed"], "flag": "primary",
        "subset": {"age": "E18.5"},
    },
}

CONTRASTS = {"brain": BRAIN_CONTRASTS, "placenta": PLACENTA_CONTRASTS}

STRESS_FOCUSED_CELL_TYPES = [
    "microglia", "oligodendrocyte_lineage", "excitatory_neurons",
    "inhibitory_neurons", "astrocytes",
]


def build_yaml_for_tissue(df: pd.DataFrame, tissue: str) -> dict:
    """Build the config dict for one tissue."""
    sub = df[df["tissue"] == tissue].copy()

    cfg = {
        "tissue": tissue,
        "group_reference": "Relaxed",       # Relaxed is baseline; +logFC = upregulated in stress
        "results_dir": f"results/{tissue}",
        "samples": [build_sample_entry(r) for _, r in sub.iterrows()],
        **SHARED_CONFIG,
    }
    # Per-tissue annotation + reference config (Phase 7 / 7c).
    ref = REFERENCE_CONFIG.get(tissue, {})
    if "celltypist_models" in ref:
        cfg["annotation"] = {"celltypist_models": ref["celltypist_models"]}
    if "reference" in ref:
        cfg["reference"] = ref["reference"]
    # Declarative contrast spec (Phase 8).
    if tissue in CONTRASTS:
        cfg["contrasts"] = CONTRASTS[tissue]
        cfg["stress_focused_cell_types"] = STRESS_FOCUSED_CELL_TYPES
    return cfg


def main():
    if not METADATA_CSV.is_file():
        sys.exit(f"ERROR: {METADATA_CSV} not found. Run from the Analysis/ directory.")

    df = pd.read_csv(METADATA_CSV)
    CONFIG_DIR.mkdir(exist_ok=True)

    for tissue in ["brain", "placenta"]:
        cfg = build_yaml_for_tissue(df, tissue)
        out = CONFIG_DIR / f"{tissue}.yaml"
        with out.open("w") as f:
            # default_flow_style=False -> block style (one key per line)
            # sort_keys=False -> preserve our intentional ordering
            yaml.safe_dump(cfg, f, default_flow_style=False, sort_keys=False, width=120)
        n = len(cfg["samples"])
        print(f"  Wrote {out}  ({n} samples)")

    # Summary
    print()
    print("Summary:")
    for tissue in ["brain", "placenta"]:
        sub = df[df["tissue"] == tissue]
        print(f"  {tissue} ({len(sub)} samples):")
        print(sub.groupby(["age", "group"]).size().to_string().replace("\n", "\n    "))
        print()


if __name__ == "__main__":
    main()
