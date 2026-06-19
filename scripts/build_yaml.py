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
        # Phase 8a (propeller): min donors/group to RUN a stratum (`min_donors`),
        # and the threshold at/above which it's 'ok' rather than 'low_n'
        # (`reliable_donors`). We RUN at n>=2 so P1 Late Stress and the sex-specific
        # strata are attempted, but FLAG any stratum with a group <3 as low_n.
        # CLI --min-donors / --reliable-donors override. (propeller replaced scCODA;
        # limma's empirical-Bayes moderation makes n=2 usable-but-weak, hence
        # run-and-flag rather than hard-skip.)
        "min_donors": 2,
        "reliable_donors": 3,
    },
    # Phase 8 sex stratification — ONE declarative dimension applied to EVERY
    # contrast in EVERY 8x stage (8a..8g), via _utils.iter_strata. 'combined' =
    # sex stays a covariate (pooled run); 'M'/'F' = subset to that sex (sex then
    # auto-drops from the design). This supersedes the standalone
    # within_age_sex_stratified contrast (which only covered Early-vs-Relaxed).
    # sex-specific strata are underpowered (~halved n) and get flagged low_n.
    "strata": {
        "sex": ["combined", "M", "F"],
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
    # Phase 4 (04_integration_prep.py) — HVG selection + cell cycle scoring.
    # n_hvg is overridden per tissue in build_yaml_for_tissue (3000 brain /
    # 2000 placenta). seurat_v3 flavor, batch_key=pool. HVG exclusion lists
    # (mito/ribo/hemo/sex + placenta pregnancy genes) are hardcoded in the
    # script, not config-driven.
    "integration": {
        "n_hvg": 3000,                      # tissue default; placenta overridden to 2000
        "batch_key": "pool",                # scVI batch_key (technical multiplex)
    },
    # Phase 5 (05_integration.py) — scVI model. accelerator/precision are
    # auto-selected by _utils.select_accelerator (gpu+bf16 on the Ada, cpu+fp32
    # otherwise), so they're not set here. condition_cell_cycle=false: cell
    # cycle is real biology (P1 progenitors, trophoblast) and a candidate stress
    # phenotype — we do NOT regress it out by default. Flip to true only if the
    # phase UMAP shows cycle driving artifactual cross-cell-type clusters.
    "scvi": {
        "n_latent": 30,
        "n_layers": 2,
        "max_epochs": 400,
        "batch_size": 1024,
        "early_stopping_patience": 30,
        "condition_cell_cycle": False,
    },
    # Phase 5 builds the neighbor graph + UMAP with these; Phase 6 REUSES that
    # graph for clustering (Option B — figure and clusters share one graph).
    # Tuned for 400-700K-cell atlases: n_neighbors=30 (cohesive clusters),
    # min_dist=0.3 + spread=1.2 (clean, well-separated islands without scatter).
    # `resolution` (Phase 6 Leiden) is set per tissue: placenta=2.0; brain omits
    # it to auto-select via the knee of the resolution-vs-nclusters curve.
    "clustering": {
        "n_neighbors": 30,
        "min_dist": 0.3,
        "spread": 1.2,
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
# celltypist_models: per-age CellTypist model (Phase 7). Brain 4W/3mo use the
#   ABC-trained pkls; P1 is null (annotated by scANVI instead — see scanvi_p1).
# scanvi_p1: P1 scANVI label-transfer config (Phase 7 subprocess run_scanvi_p1.py).
#   - ref_h5ad: Rosenberg P2-brain reference (prepare_rosenberg_reference.py).
#   - labels_key: .obs column in the reference holding the fine cell-type label.
#   - config_dir: where the rosenberg_*.csv derivation maps live (default config).
# -----------------------------------------------------------------------------

REFERENCE_CONFIG = {
    "brain": {
        # Three-tier CellTypist (Phase 7), per age:
        #   class    -> canonical per-(cluster x age) label (celltypist_class),
        #               assigned by majority vote of per-cell class predictions
        #               WITHIN each (Phase 6 Leiden cluster x age) group.
        #   subclass -> kept PER-CELL (celltypist_subclass); consumed at the
        #               subcluster level by 7b/7d.
        #   region   -> kept PER-CELL (celltypist_region); used at Phase 9 for
        #               regional matching to human datasets (e.g. mouse Isocortex
        #               cells vs human dlPFC). Per-cell from a model = context-
        #               aware, unlike abc_subclass_to_region.csv which buries
        #               cross-regional subclasses (cortical interneurons) in
        #               `multi_region`.
        # P1 is annotated by scANVI label transfer from Rosenberg 2018 (see
        # scanvi_p1 below), NOT CellTypist — Di Bella is cortex-only and
        # mislabeled ~42% of whole-brain P1 cells as erythrocyte. So P1's
        # CellTypist entry is null; Phase 7 routes P1 to the scANVI subprocess.
        # 4W/3mo share the ABC-trained pkls written by train_celltypist_brain.py.
        "celltypist_models": {
            "P1":  {"class": None,                                            # P1 -> scANVI (see scanvi_p1)
                    "subclass": None,
                    "region": None},
            "4W":  {"class": "refs/celltypist_brain_adult_class.pkl",         # ABC WMB-10Xv3 class (34)
                    "subclass": "refs/celltypist_brain_adult_subclass.pkl",   # ABC subclass (~334)
                    "region": "refs/celltypist_brain_adult_region.pkl"},      # ABC anatomical_division_label (~12)
            "3mo": {"class": "refs/celltypist_brain_adult_class.pkl",
                    "subclass": "refs/celltypist_brain_adult_subclass.pkl",
                    "region": "refs/celltypist_brain_adult_region.pkl"},
        },
        # P1 scANVI label transfer (Phase 7 subprocess: run_scanvi_p1.py).
        # Transfers the published Rosenberg fine label; Phase 7 derives
        # class/region/broad via the config/rosenberg_*.csv maps written by
        # prepare_rosenberg_reference.py.
        "scanvi_p1": {
            "ref_h5ad":   "refs/rosenberg_p2brain_reference.h5ad",
            "labels_key": "rosenberg_fine",
            "config_dir": "config",
        },
        # Deterministic class->broad map. Phase 7 now derives celltypist_broad
        # for ALL ages: ABC classes (4W/3mo) via class_to_broad_csv, Rosenberg
        # classes (P1) via config/rosenberg_class_to_broad.csv (auto-found via
        # scanvi_p1.config_dir). subclass->region CSV is informational only.
        "class_to_broad_csv":     "refs/abc_class_to_broad.csv",
        "subclass_to_region_csv": "refs/abc_subclass_to_region.csv",
    },
    "placenta": {
        # No built-in CellTypist placenta model. Primary track is STAMP
        # reference-based Spearman correlation (Liu et al. 2024 eLife) —
        # see scripts/build_placenta_reference.py for the .h5 reference.
        "celltypist_models": {},
        "stamp": {
            "reference":           "refs/stamp/stamp_ref_allcells.h5",
            "gap_threshold":       0.05,   # per-cell low-confidence threshold
            "purity_threshold":    0.5,    # per-cluster low-purity threshold
            "min_cluster_size":    50,     # clusters smaller -> 'under_populated'
        },
    },
}



def load_markers_yaml(tissue, config_dir="config"):
    """Merge config/{tissue}_markers.yaml into the annotation block if present.

    Expected structure in the marker file:
        annotation:
          markers:
            "Cell type A": [GeneA, GeneB, ...]
    Returns the markers dict (empty if file absent). Keeps marker source files
    out of build_yaml.py itself — they live in config/ and survive regen.
    """
    from pathlib import Path as _P
    import yaml as _yaml
    mp = _P(config_dir) / f"{tissue}_markers.yaml"
    if not mp.exists():
        return {}
    doc = _yaml.safe_load(mp.read_text())
    return (doc.get("annotation", {}) or {}).get("markers", {}) or {}


def build_sample_entry(row: pd.Series) -> dict:
    """Convert one CSV row into one sample dict for the YAML.

    Sex field precedence: `assigned_sex` (written by Phase 0 from Y-score
    + Xist ambiguity rule) overrides `sex_declared`. `assigned_sex` is the
    source of truth for downstream covariates — see INSTRUCTIONS.md.
    Falls back to sex_declared (with TBD → unknown) if assigned_sex is
    missing (e.g. before Phase 0 has run).
    """
    assigned = row.get("assigned_sex") if "assigned_sex" in row.index else None
    if pd.notna(assigned) and str(assigned).strip() != "":
        sex = str(assigned)
    else:
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
    # NOTE: within_age_sex_stratified was REMOVED — superseded by the declarative
    # sex strata (combined/M/F applied to EVERY contrast via _utils.iter_strata).
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
    import copy
    sub = df[df["tissue"] == tissue].copy()

    # Deep-copy SHARED_CONFIG so per-tissue overrides below don't leak into the
    # other tissue (nested dicts like integration/clustering are mutated).
    shared = copy.deepcopy(SHARED_CONFIG)

    cfg = {
        "tissue": tissue,
        "group_reference": "Relaxed",       # Relaxed is baseline; +logFC = upregulated in stress
        "results_dir": f"results/{tissue}",
        "samples": [build_sample_entry(r) for _, r in sub.iterrows()],
        **shared,
    }

    # --- Tissue-specific overrides ---
    if tissue == "placenta":
        # Placenta uses fewer HVGs (less cell-type diversity than whole brain);
        # pregnancy genes already excluded in Phase 4.
        cfg["integration"]["n_hvg"] = 2000
        # Placenta Leiden resolution locked at 2.0 (resolves trophoblast subtypes);
        # brain omits `resolution` to auto-select via the knee curve.
        cfg["clustering"]["resolution"] = 2.0
    # Per-tissue annotation + reference config (Phase 7 / 7c).
    ref = REFERENCE_CONFIG.get(tissue, {})
    if "celltypist_models" in ref:
        annotation = {"celltypist_models": ref["celltypist_models"]}
        # Carry the deterministic mapping CSVs + P1 scANVI block if declared.
        for k in ("class_to_broad_csv", "subclass_to_region_csv", "scanvi_p1"):
            if k in ref:
                annotation[k] = ref[k]
        cfg["annotation"] = annotation
    # Placenta STAMP primary track — flatten the stamp block into annotation
    # under the keys that 07_annotation.py reads (stamp_reference, etc.).
    if "stamp" in ref:
        stamp = ref["stamp"]
        ann = cfg.setdefault("annotation", {})
        ann["stamp_reference"]        = stamp.get("reference")
        ann["stamp_gap_threshold"]    = stamp.get("gap_threshold", 0.05)
        ann["stamp_purity_threshold"] = stamp.get("purity_threshold", 0.5)
        ann["stamp_min_cluster_size"] = stamp.get("min_cluster_size", 50)
    # NOTE: the old top-level `reference:` block (Phase 7c scANVI) is removed —
    # 7c is deleted; P1 scANVI now lives under annotation.scanvi_p1 and runs
    # inside Phase 7 via the run_scanvi_p1.py subprocess.
    # Merge curated marker sets from config/{tissue}_markers.yaml if present.
    # Survives regen — marker source files live in config/, never overwritten.
    markers = load_markers_yaml(tissue)
    if markers:
        cfg.setdefault("annotation", {})["markers"] = markers
        print(f"  [{tissue}] merged {len(markers)} marker sets from config/{tissue}_markers.yaml")
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
