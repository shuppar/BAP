#!/usr/bin/env python
"""
08e_comms.py — Phase 8e: cell-cell communication (liana-py). COMPUTE ONLY.

Writes complete, source-of-truth CSVs. All plotting lives in
`08e_comms_summary.py`, which reads ONLY these CSVs (mirrors the 8b/8c
compute/summary split). This script does no plotting and imports no plot module.

Three arms (toggle with --arms; each writes only its own CSV(s) so one arm can
be re-run without touching the others):

  baseline    — rank_aggregate per group×age (pooled cells). DESCRIPTIVE
                landscape only, NOT a stress test (pooled cells = no donor
                structure). Writes 08e_lr_baseline.csv (+ specificity_fdr if
                permutations are on) and 08e_sender_receiver.csv.
  differential — df_to_lr on 8b Wald stats (PRIMARY inferential arm). Inherits
                8b's animal-as-unit pseudobulk design; shares signal with 8b by
                construction (it re-expresses 8b DE in LR space, a transcriptional
                proxy for signalling change — frame as such). Feeds 8f view 4.
                Writes 08e_lr_differential.csv.
  perdonor    — rank_aggregate per donor → Mann-Whitney U across donors
                (independent-design corroboration; animal is the unit). Expected
                mostly null / low_n at n≈4. Writes 08e_lr_per_donor.csv (+
                08e_lr_quantified.csv with per-test p/FDR).

Covers all three group comparisons: ES-v-Rel, LS-v-Rel, ES-v-LS.

CANONICAL CELL-TYPE KEY (see snRNAseq_project_summary.md ⚠️ note):
  brain    -> celltypist_broad   (8 types at level=='whole')
  placenta -> celltype_majority  (~20 types)
The differential arm filters the 8b master CSV to level=='whole' and joins on
this key. Contrast names differ by tissue (brain *_per_age with age in
group_level; placenta age-baked *_E12.5/*_E18.5) — resolved by family prefix.

No lognorm layer is persisted in the h5ad (raw counts in .X) -> add_lognorm()
after load.

Usage:
  # smoke test (no permutations — fast; no specificity p-values)
  uv run python scripts/08e_comms.py --config config/placenta.yaml --n-perms 0
  # production
  uv run python scripts/08e_comms.py --config config/brain.yaml \\
      --n-perms 1000 --n-jobs 8
  # re-run a single arm (e.g. after the 8b CSV changes)
  uv run python scripts/08e_comms.py --config config/brain.yaml \\
      --arms differential --n-perms 1000

Deps: uv add liana
"""

import argparse
import os
import sys
import warnings
from itertools import combinations
from pathlib import Path

# Pin per-process thread pools to 1 BEFORE importing numpy/scanpy. The per-donor
# arm fans out across processes; without this, 8 workers × many BLAS/OpenMP
# threads oversubscribe the box. Must be set before the first numpy import.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from _utils import load_config, add_lognorm, phase_table_dir, parallel_map

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

PHASE = "08e_communication"   # keep folder name stable for downstream consumers

# Canonical whole-level cell-type key per tissue (see ⚠️ note in project summary).
CANONICAL_KEY = {"brain": "celltypist_broad", "placenta": "celltype_majority"}
# Fallback search order if tissue not in CANONICAL_KEY or key absent.
LABEL_KEY_PRIORITY = [
    "subcluster_name", "subcluster",          # subcluster objects (deferred scope)
    "celltypist_broad", "celltype_majority",
    "manual_annotation", "celltypist_class",
]


# ============================================================================
# CLI
# ============================================================================

def parse_args():
    p = argparse.ArgumentParser(description="Phase 8e compute: cell-cell communication")
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--arms", default="baseline,differential,perdonor",
                   help="Comma-separated subset of {baseline,differential,perdonor}. "
                        "Each arm writes only its own CSV(s).")
    p.add_argument("--n-perms", type=int, default=1000,
                   help="LIANA permutations (0 = skip -> no specificity_rank/p-values; "
                        "1000 = production/CellPhoneDB standard).")
    p.add_argument("--min-cells", type=int, default=5,
                   help="Min cells per cell type (default 5).")
    p.add_argument("--expr-prop", type=float, default=0.1,
                   help="Min expression proportion per LR gene (default 0.1).")
    p.add_argument("--n-jobs", type=int, default=8,
                   help="Parallel workers for the per-donor arm (default 8). "
                        "Each LIANA call is RAM-hungry with permutations on.")
    p.add_argument("--celltype-key", default=None,
                   help="Override the obs cell-type column (default: canonical per tissue).")
    p.add_argument("--subcluster", default=None,
                   help="Run on an 08c subcluster object instead of the main object "
                        "(slug = h5ad basename, e.g. 'opc_oligodendrocytes'). "
                        "Scope note: subcluster loop is deferred; main objects only.")
    p.add_argument("--node-scheme", default="broad", choices=["broad", "subtype"],
                   help="'broad' = canonical key (celltypist_broad / celltype_majority). "
                        "'subtype' = build comms_subtype: focal subcluster substates "
                        "(>= --min-node-cells) as nodes, smaller ones + non-focal cells "
                        "kept at parent broad label. Exploratory. Main object only.")
    p.add_argument("--levels", default="whole",
                   help="Comma-separated subset of {whole,regional}. 'regional' adds one "
                        "baseline slice per celltypist_region (brain only). Baseline arm only.")
    p.add_argument("--min-node-cells", type=int, default=300,
                   help="--node-scheme subtype: a substate needs >= this many cells to be "
                        "its own node; smaller substates collapse to the parent broad label.")
    p.add_argument("--head", type=int, default=0, metavar="N",
                   help="ULTRA-FAST smoke: keep only the first N cells of the object "
                        "right after load (before subtype build + per-donor slicing), "
                        "and force n_perms=0. Validates code paths in ~1-2 min; "
                        "numbers MEANINGLESS. Try --head 5000.")
    p.add_argument("--smoke", type=int, default=0, metavar="N",
                   help="Smoke mode: subsample to N cells per (group×age) and force "
                        "n_perms=0. Validates code paths + CSV schemas in seconds; "
                        "numbers are MEANINGLESS. Try --smoke 500.")
    return p.parse_args()


# ============================================================================
# Loading / key resolution
# ============================================================================

def resolve_celltype_key(adata, tissue, explicit=None, subcluster=None):
    if explicit:
        if explicit not in adata.obs.columns:
            sys.exit(f"ERROR: --celltype-key '{explicit}' not in adata.obs.")
        print(f"  Cell-type key (explicit): '{explicit}'")
        return explicit
    # Subcluster objects key on the fine subtype label, NOT the main-object
    # canonical (which holds the coarse parent label, e.g. 'Immune').
    if subcluster:
        for k in ("subcluster_name", "subcluster"):
            if k in adata.obs.columns:
                print(f"  Cell-type key (subcluster): '{k}'")
                return k
        sys.exit("ERROR: subcluster run but no 'subcluster_name'/'subcluster' column.")
    key = CANONICAL_KEY.get(tissue)
    if key and key in adata.obs.columns:
        print(f"  Cell-type key (canonical for {tissue}): '{key}'")
        return key
    for k in LABEL_KEY_PRIORITY:
        if k in adata.obs.columns:
            print(f"  Cell-type key (fallback): '{k}'")
            return k
    sys.exit(f"ERROR: no cell-type label column found. Looked for: {LABEL_KEY_PRIORITY}")


def check_var_names_are_symbols(adata):
    if any(str(v).startswith("ENSMUSG") for v in list(adata.var_names[:5])):
        sys.exit("ERROR: var_names look like Ensembl IDs — liana needs gene symbols.")


def load_annotated_h5ad(cfg, subcluster_slug=None):
    base = Path(cfg["results_dir"]) / "h5ad"
    if subcluster_slug:
        p = base / "08c_subclustered" / f"{subcluster_slug}.h5ad"
        if not p.is_file():
            sys.exit(f"ERROR: subcluster h5ad not found: {p}")
        print(f"  Loading subcluster: {p}")
        return sc.read_h5ad(p)
    p = base / "08_annotated" / "all_samples.h5ad"
    if not p.is_file():
        sys.exit(f"ERROR: annotated h5ad not found: {p}\n  Run Phase 7 first.")
    print(f"  Loading: {p}")
    return sc.read_h5ad(p)


def build_comms_subtype(adata, cfg, tissue, min_node_cells=300):
    """Add obs['comms_subtype']: focal-subcluster substate where a cell was focally
    subclustered AND that substate has >= min_node_cells, else the parent broad
    label. Contamination_* / unresolved / unassigned* cells are dropped (returned
    adata is subset). Mapping is by obs_names (verified 100% match). Reads the focal
    subcluster objects + parent map from config/stress_pathways_8e.yaml."""
    import yaml
    spec_path = Path("config/stress_pathways_8e.yaml")
    if not spec_path.is_file():
        sys.exit(f"ERROR: {spec_path} not found (needed for --node-scheme subtype).")
    spec = yaml.safe_load(spec_path.read_text())
    focal = spec.get(tissue, {}).get("focal_subclusters", {})
    if not focal:
        sys.exit(f"ERROR: no focal_subclusters for '{tissue}' in {spec_path}.")

    broad_key = CANONICAL_KEY.get(tissue)
    if broad_key not in adata.obs.columns:
        sys.exit(f"ERROR: canonical key '{broad_key}' missing from main object.")
    # start from the broad label
    comms = adata.obs[broad_key].astype(str).copy()

    base = Path(cfg["results_dir"]) / "h5ad" / "08c_subclustered"
    for slug, parent in focal.items():
        p = base / f"{slug}.h5ad"
        if not p.is_file():
            print(f"  [warn] subcluster object missing: {p} — skipping {slug}")
            continue
        s = sc.read_h5ad(p)
        col = "subcluster_name" if "subcluster_name" in s.obs else "subcluster"
        lab = s.obs[col].astype(str)
        # drop contamination / unresolved within the substate labels
        keepmask = ~(lab.str.startswith("Contamination_") | (lab == "unresolved"))
        lab = lab[keepmask]
        # substates that clear the cell-count floor become nodes; others -> parent
        vc = lab.value_counts()
        node_states = set(vc[vc >= min_node_cells].index)
        lab = lab.map(lambda x: x if x in node_states else parent)
        # write onto the main label by barcode
        common = adata.obs_names.intersection(lab.index)
        comms.loc[common] = lab.loc[common].values
        print(f"  {slug}: {len(node_states)} substate nodes "
              f"({sorted(node_states)}) | {len(common):,} cells mapped; "
              f"<{min_node_cells} or contamination -> '{parent}'")

    adata.obs["comms_subtype"] = comms.values
    # drop unassigned / contamination that may remain at the broad level
    drop = (adata.obs["comms_subtype"].astype(str).str.startswith("unassigned") |
            adata.obs["comms_subtype"].astype(str).str.startswith("Contamination_") |
            (adata.obs["comms_subtype"].astype(str) == "unresolved"))
    if drop.any():
        print(f"  Dropping {int(drop.sum()):,} unassigned/contamination cells from comms_subtype.")
        adata = adata[~drop.values].copy()
    print(f"  comms_subtype nodes ({adata.obs['comms_subtype'].nunique()}): "
          f"{sorted(adata.obs['comms_subtype'].astype(str).unique())}")
    return adata


def load_de_results(cfg, subcluster_slug=None):
    """Load 08b DE master CSV (differential arm). Filtered downstream to
    level=='whole' so the join key matches the canonical whole-level vocabulary."""
    tdir = Path(cfg["results_dir"]) / "tables" / "08b_de"
    if subcluster_slug:
        p = tdir / f"08b_de_results_subcluster_{subcluster_slug}.csv"
    else:
        p = tdir / "08b_de_results.csv"
    if not p.is_file():
        print(f"  [warn] 08b DE results not found: {p.name} — differential arm skipped.")
        return pd.DataFrame()
    df = pd.read_csv(p, low_memory=False)
    print(f"  Loaded DE results: {len(df):,} rows from {p.name}")
    return df


# ============================================================================
# Arm 1: Baseline — rank_aggregate per group×age (pooled). DESCRIPTIVE.
# ============================================================================

def _baseline_one_slice(adata, cfg, args, celltype_key, level_label, region_key=None):
    """Run rank_aggregate per group×age within one level (whole, or one region).
    Returns list of per-(group,age) liana_res frames tagged with level."""
    import liana as li
    out = []
    if level_label != "whole" and region_key is not None:
        adata = adata[adata.obs[region_key].astype(str) == level_label].copy()
        if adata.n_obs < args.min_cells * 4:
            print(f"  SKIP level={level_label}: {adata.n_obs} cells")
            return out
    n_perms = args.n_perms if args.n_perms > 0 else None
    for age in sorted(adata.obs["age"].astype(str).unique()):
        for group in sorted(adata.obs["group"].astype(str).unique()):
            mask = ((adata.obs["age"].astype(str) == age) &
                    (adata.obs["group"].astype(str) == group)).values
            if mask.sum() < args.min_cells * 2:
                continue
            sub = adata[mask].copy()
            vc = sub.obs[celltype_key].value_counts()
            valid = vc[vc >= args.min_cells].index
            if len(valid) < 2:
                continue
            sub = sub[sub.obs[celltype_key].isin(valid)].copy()
            print(f"  rank_aggregate: level={level_label} {group}/{age} "
                  f"({sub.n_obs:,} cells, {len(valid)} types, n_perms={n_perms})")
            try:
                li.mt.rank_aggregate(
                    sub, groupby=celltype_key, resource_name="mouseconsensus",
                    expr_prop=args.expr_prop, min_cells=args.min_cells,
                    use_raw=False, layer="lognorm",
                    n_perms=n_perms, seed=cfg.get("random_seed", 42),
                    n_jobs=1, verbose=False, inplace=True,
                )
            except Exception as e:
                print(f"  [warn] level={level_label} {group}/{age} failed: {e}")
                continue
            res = sub.uns["liana_res"].copy()
            res["group"] = group; res["age"] = age; res["level"] = level_label
            out.append(res)
    return out


def run_baseline(adata, cfg, args, tdir, celltype_key):
    from statsmodels.stats.multitest import multipletests
    print("\n[baseline] rank_aggregate per group×age (DESCRIPTIVE landscape)")

    # Levels: 'whole' always; 'regional' adds one slice per celltypist_region
    # (brain only — placenta has no region column).
    levels = [s.strip() for s in args.levels.split(",") if s.strip()]
    region_key = "celltypist_region"
    do_regional = ("regional" in levels) and (region_key in adata.obs.columns)
    if "regional" in levels and not do_regional:
        print(f"  [info] regional requested but no '{region_key}' column — whole only")

    all_results = []
    if "whole" in levels:
        all_results += _baseline_one_slice(adata, cfg, args, celltype_key, "whole")
    if do_regional:
        regions = [r for r in sorted(adata.obs[region_key].astype(str).unique())
                   if r not in ("whole", "nan", "non-regional")]
        print(f"  regional levels: {regions}")
        for reg in regions:
            all_results += _baseline_one_slice(adata, cfg, args, celltype_key,
                                               reg, region_key)

    if not all_results:
        print("  [warn] No baseline results.")
        return pd.DataFrame()

    combined = pd.concat(all_results, ignore_index=True)

    # BH-FDR on the specificity p-value within group×age×level (descriptive
    # significance — cell-as-unit, never a stress claim). Only with perms on.
    # LIANA's permutation p-value column is 'cellphone_pvals' (rank_aggregate with
    # n_perms>0). Match it explicitly, plus a generic fallback.
    spec_p = "cellphone_pvals" if "cellphone_pvals" in combined.columns else None
    if spec_p is None:
        spec_p = next((c for c in combined.columns
                       if "pval" in c.lower() and combined[c].notna().any()), None)
    if spec_p:
        parts = []
        for _, grp in combined.groupby(["group", "age", "level"]):
            _, fdr, _, _ = multipletests(grp[spec_p].fillna(1.0), method="fdr_bh")
            parts.append(pd.Series(fdr, index=grp.index))
        combined["specificity_fdr"] = pd.concat(parts)
        print(f"  Added specificity_fdr (BH within group×age×level) from '{spec_p}'")
    else:
        print("  [info] no specificity p-value column (n_perms=0?) — specificity_fdr skipped")

    out = tdir / "08e_lr_baseline.csv"
    combined.to_csv(out, index=False)
    print(f"  Saved: {out.name}  ({len(combined):,} rows)")

    sr = compute_sender_receiver(combined)
    if not sr.empty:
        sr.to_csv(tdir / "08e_sender_receiver.csv", index=False)
        print(f"  Saved: 08e_sender_receiver.csv  ({len(sr):,} rows)")
    return combined


def compute_sender_receiver(baseline_df, magnitude_cutoff=0.05):
    if baseline_df.empty or "magnitude_rank" not in baseline_df.columns:
        return pd.DataFrame()
    records = []
    for (group, age), grp in baseline_df.groupby(["group", "age"]):
        active = grp[grp["magnitude_rank"] < magnitude_cutoff]
        for ct in pd.concat([active["source"], active["target"]]).unique():
            src = active[active["source"] == ct]
            tgt = active[active["target"] == ct]
            records.append({
                "group": group, "age": age, "cell_type": ct,
                "n_sent": len(src), "n_received": len(tgt),
                "sender_score": float((1 - src["magnitude_rank"]).mean()) if len(src) else 0.0,
                "receiver_score": float((1 - tgt["magnitude_rank"]).mean()) if len(tgt) else 0.0,
            })
    return pd.DataFrame(records)


# ============================================================================
# Arm 2: Differential — df_to_lr on 8b Wald stats (PRIMARY inferential).
# ============================================================================

def _contrast_family(name):
    """Resolve a contrast name (brain *_per_age / placenta *_E12.5) to a family
    prefix. Never match contrast names by exact string across tissues."""
    n = str(name).lower()
    if n.startswith("early_vs_relaxed"):
        return "EVR"
    if n.startswith("late_vs_relaxed"):
        return "LVR"
    if n.startswith("early_vs_late"):
        return "EVL"
    return None


def run_differential(adata, de_df, cfg, args, tdir, celltype_key):
    import liana as li
    from statsmodels.stats.multitest import multipletests
    print("\n[differential] df_to_lr on 8b Wald stats (PRIMARY inferential arm)")

    if de_df.empty:
        print("  Skipped — no DE results.")
        return pd.DataFrame()

    required = {"contrast", "celltype", "group_level", "level", "gene",
                "stat", "padj", "log2FC", "flag"}
    missing = required - set(de_df.columns)
    if missing:
        print(f"  [warn] DE table missing columns {missing} — skipping. "
              f"Available: {sorted(de_df.columns)}")
        return pd.DataFrame()

    # Filter to whole-level, gene-level, primary/secondary contrasts. The
    # whole-level celltype vocabulary == the canonical join key (celltypist_broad
    # / celltype_majority); region-level rows use a different vocabulary and must
    # NOT enter the join.
    use = de_df[de_df["level"].astype(str) == "whole"].copy()
    use = use.dropna(subset=["gene", "stat"])
    use = use[use["flag"].isin(["primary", "secondary"])]
    if use.empty:
        print("  [warn] No whole-level primary/secondary gene rows. Skipping.")
        return pd.DataFrame()

    # group_level holds the AGE directly (verified: 'P1'/'4W'/'3mo' or
    # 'E12.5'/'E18.5'); NOT an 'age-...' encoding.
    use["age"] = use["group_level"].astype(str)
    valid_ages = set(adata.obs["age"].astype(str).unique())
    use = use[use["age"].isin(valid_ages)]
    if use.empty:
        print(f"  [warn] No DE rows whose group_level matches data ages {sorted(valid_ages)}.")
        return pd.DataFrame()

    use["cfam"] = use["contrast"].map(_contrast_family)
    use = use[use["cfam"].notna()]
    print(f"  Contrasts: {sorted(use['contrast'].unique())}")
    print(f"  Ages:      {sorted(use['age'].unique())}")

    # Sanity: do the DE whole-level celltypes overlap the h5ad key categories?
    de_cts = set(use["celltype"].astype(str).unique())
    ad_cts = set(adata.obs[celltype_key].astype(str).unique())
    overlap = de_cts & ad_cts
    if not overlap:
        sys.exit(f"ERROR: no overlap between 8b whole-level celltype values and "
                 f"adata['{celltype_key}']. DE has {sorted(de_cts)[:5]}..., "
                 f"adata has {sorted(ad_cts)[:5]}.... Wrong join key?")
    print(f"  Celltype join overlap: {len(overlap)}/{len(de_cts)} DE types match adata.")

    all_lr = []
    for contrast_name in sorted(use["contrast"].unique()):
        for age in sorted(use["age"].unique()):
            sub_de = use[(use["contrast"] == contrast_name) & (use["age"] == age)]
            if sub_de.empty:
                continue
            mask = (adata.obs["age"].astype(str) == age).values
            if mask.sum() < args.min_cells * 2:
                continue
            sub_adata = adata[mask].copy()
            dea = sub_de[["celltype", "gene", "stat", "padj", "log2FC"]].rename(
                columns={"log2FC": "log2fc", "celltype": celltype_key})
            # df_to_lr requires the frame INDEXED BY GENE (matching var_names),
            # not a 'gene' column.
            dea = dea.set_index("gene")
            print(f"  df_to_lr: {contrast_name}/{age} "
                  f"({dea[celltype_key].nunique()} cell types)")
            try:
                lr = li.multi.df_to_lr(
                    sub_adata, dea_df=dea, resource_name="mouseconsensus",
                    expr_prop=args.expr_prop, groupby=celltype_key,
                    stat_keys=["stat", "padj", "log2fc"],
                    use_raw=False, layer="lognorm",
                    complex_col="stat", verbose=False, return_all_lrs=False,
                )
            except Exception as e:
                print(f"  [warn] df_to_lr failed {contrast_name}/{age}: {e}")
                continue
            lr["contrast_name"] = contrast_name
            lr["contrast_family"] = _contrast_family(contrast_name)
            lr["age"] = age
            all_lr.append(lr)

    if not all_lr:
        print("  [warn] No differential LR results.")
        return pd.DataFrame()

    combined = pd.concat(all_lr, ignore_index=True)
    pval_col = next((c for c in combined.columns if "pvalue" in c.lower()), None)
    if pval_col:
        parts = []
        for _, grp in combined.groupby(["contrast_name", "age", "source", "target"]):
            _, fdr, _, _ = multipletests(grp[pval_col].fillna(1.0), method="fdr_bh")
            parts.append(pd.Series(fdr, index=grp.index))
        combined["interaction_fdr"] = pd.concat(parts)

    out = tdir / "08e_lr_differential.csv"
    combined.to_csv(out, index=False)
    print(f"  Saved: {out.name}  ({len(combined):,} rows, "
          f"{combined['contrast_name'].nunique()} contrasts)")
    return combined


# ============================================================================
# Arm 3: Per-donor — rank_aggregate per donor → MW-U across donors.
# ============================================================================

def _per_donor_worker(payload):
    """Top-level (picklable) worker: run rank_aggregate on one donor's PRE-SLICED
    h5ad, optionally restricted to one region. Returns the liana_res df tagged with
    donor/group/age/level, or None / ('ERR', donor, msg)."""
    import liana as li
    (donor_h5ad, donor, group, age, level, region_key, celltype_key,
     expr_prop, min_cells, n_perms, seed) = payload
    sub = sc.read_h5ad(donor_h5ad)
    if level != "whole" and region_key is not None and region_key in sub.obs.columns:
        sub = sub[sub.obs[region_key].astype(str) == level].copy()
        if sub.n_obs < min_cells * 2:
            return None
    add_lognorm(sub)
    vc = sub.obs[celltype_key].value_counts()
    valid = vc[vc >= min_cells].index
    if len(valid) < 2:
        return None
    sub = sub[sub.obs[celltype_key].isin(valid)].copy()
    try:
        li.mt.rank_aggregate(
            sub, groupby=celltype_key, resource_name="mouseconsensus",
            expr_prop=expr_prop, min_cells=min_cells,
            use_raw=False, layer="lognorm",
            n_perms=(n_perms if n_perms > 0 else None), seed=seed,
            n_jobs=1, verbose=False, inplace=True,
        )
    except Exception as e:
        return ("ERR", f"{donor}/{level}", str(e))
    res = sub.uns["liana_res"].copy()
    res["donor_id"] = donor
    res["group"] = group
    res["age"] = age
    res["level"] = level
    return res


def run_per_donor(h5ad_path, adata, cfg, args, tdir, celltype_key):
    from scipy import stats as scipy_stats
    from statsmodels.stats.multitest import multipletests
    print("\n[perdonor] rank_aggregate per donor → MW-U across donors (corroboration)")

    donors = sorted(adata.obs["donor_id"].astype(str).unique())

    # Levels: 'whole' always; 'regional' adds per-region per-donor slices
    # (brain only — placenta has no region column).
    levels = [s.strip() for s in args.levels.split(",") if s.strip()]
    region_key = "celltypist_region"
    do_regional = ("regional" in levels) and (region_key in adata.obs.columns)
    region_list = ["whole"]
    if do_regional:
        region_list += [r for r in sorted(adata.obs[region_key].astype(str).unique())
                        if r not in ("whole", "nan", "non-regional")]
    print(f"  {len(donors)} donors × {len(region_list)} levels "
          f"({region_list}) | n_jobs={args.n_jobs}")

    # Pre-slice each donor to a small h5ad ONCE (cheap, serial); region subsetting
    # happens inside the worker so we don't write donor×region files.
    slice_dir = tdir / "_per_donor_slices"
    slice_dir.mkdir(parents=True, exist_ok=True)
    rk = region_key if do_regional else None
    payloads = []
    for d in donors:
        m = (adata.obs["donor_id"].astype(str) == d).values
        sub = adata[m].copy()
        meta = sub.obs[["group", "age"]].iloc[0]
        sp = slice_dir / f"{d}.h5ad"
        sub.write_h5ad(sp)
        for lvl in region_list:
            payloads.append((str(sp), d, str(meta["group"]), str(meta["age"]),
                             lvl, rk, celltype_key, args.expr_prop, args.min_cells,
                             args.n_perms, cfg.get("random_seed", 42)))
        del sub
    print(f"  Pre-sliced {len(donors)} donor h5ads → {len(payloads)} (donor×level) jobs")

    all_results = []
    for _, res, err in parallel_map(
            _per_donor_worker, payloads, n_jobs=args.n_jobs,
            use_threads=False, desc="per-donor LIANA"):
        if err:
            print(f"    [warn] worker error: {err}")
            continue
        if res is None:
            continue
        if isinstance(res, tuple) and res[0] == "ERR":
            print(f"    SKIP {res[1]}: {res[2]}")
            continue
        all_results.append(res)
        print(f"    {res['donor_id'].iloc[0]}: {res['group'].iloc[0]}/"
              f"{res['age'].iloc[0]}/{res['level'].iloc[0]}  {len(res):,} LR pairs")

    # Clean up the slice files (they're transient scratch).
    import shutil
    shutil.rmtree(slice_dir, ignore_errors=True)

    if not all_results:
        print("  [warn] No per-donor results.")
        return pd.DataFrame()

    pdf = pd.concat(all_results, ignore_index=True)
    pdf["activity"] = 1 - pdf["magnitude_rank"]
    pdf["lr_pair"] = pdf["ligand_complex"] + "→" + pdf["receptor_complex"]
    pdf["ct_pair"] = pdf["source"] + "→" + pdf["target"]
    pdf.to_csv(tdir / "08e_lr_per_donor.csv", index=False)
    print(f"  Saved: 08e_lr_per_donor.csv  ({len(pdf):,} rows)")

    # All pairwise group comparisons per LR×ct_pair×age (Relaxed = ctrl side).
    groups = sorted(pdf["group"].unique())
    ref = cfg.get("group_reference", "Relaxed")
    pairs = []
    for a, b in combinations(groups, 2):
        if b == ref:
            pairs.append((a, b))
        elif a == ref:
            pairs.append((b, a))
        else:
            pairs.append((a, b))

    records = []
    for (lr, ct, age, level), grp in pdf.groupby(["lr_pair", "ct_pair", "age", "level"]):
        for test_g, ctrl_g in pairs:
            a = grp.loc[grp["group"] == test_g, "activity"].values
            b = grp.loc[grp["group"] == ctrl_g, "activity"].values
            if len(a) < 2 or len(b) < 2:
                continue
            stat, p = scipy_stats.mannwhitneyu(a, b, alternative="two-sided")
            records.append({
                "lr_pair": lr, "ct_pair": ct, "age": age, "level": level,
                "test_group": test_g, "ctrl_group": ctrl_g,
                "contrast": f"{test_g}_vs_{ctrl_g}",
                "mean_test": float(a.mean()), "mean_ctrl": float(b.mean()),
                "delta_activity": float(a.mean() - b.mean()),
                "n_test": len(a), "n_ctrl": len(b),
                "mannwhitney_stat": float(stat), "pvalue": float(p),
                "reliability": "ok" if (len(a) >= 3 and len(b) >= 3) else "low_n",
            })

    if not records:
        print("  [warn] No group comparisons (too few donors per group).")
        return pdf

    quant = pd.DataFrame(records)
    parts = []
    for _, sl in quant.groupby(["age", "contrast", "level"]):
        _, fdr, _, _ = multipletests(sl["pvalue"].fillna(1.0), method="fdr_bh")
        parts.append(pd.Series(fdr, index=sl.index))
    quant["fdr"] = pd.concat(parts)
    quant["significant"] = quant["fdr"] < 0.05
    quant.to_csv(tdir / "08e_lr_quantified.csv", index=False)
    n_sig = int(quant["significant"].sum())
    print(f"  Saved: 08e_lr_quantified.csv  ({len(quant):,} tests, {n_sig} FDR<0.05)")
    return pdf


# ============================================================================
# Main
# ============================================================================

def main():
    args = parse_args()

    # Force 'spawn' so the per-donor ProcessPoolExecutor starts fresh interpreters
    # instead of fork()ing a process that has already initialised OpenMP/BLAS
    # (fork-after-OpenMP is unsafe and aborts the workers).
    if "perdonor" in args.arms:
        import multiprocessing as mp
        try:
            mp.set_start_method("spawn", force=True)
        except RuntimeError:
            pass

    cfg = load_config(args.config)
    tissue = cfg.get("tissue", "unknown")
    arms = {a.strip() for a in args.arms.split(",") if a.strip()}
    valid_arms = {"baseline", "differential", "perdonor"}
    bad = arms - valid_arms
    if bad:
        sys.exit(f"ERROR: unknown --arms {bad}. Valid: {sorted(valid_arms)}")

    print(f"\n{'='*60}\nPhase 8e compute: cell-cell communication  [{tissue}]")
    print(f"  arms={sorted(arms)}  n_perms={args.n_perms}  "
          f"subcluster={args.subcluster or '(main)'}\n{'='*60}")

    try:
        import liana as li
        print(f"  liana {li.__version__}")
    except ImportError:
        sys.exit("ERROR: liana not installed. Run: uv add liana")

    # Output dir (suffixed for subcluster runs, like 8b/8c; and for subtype scheme).
    label = PHASE if not args.subcluster else f"{PHASE}_subcluster_{args.subcluster}"
    if args.node_scheme == "subtype" and not args.subcluster:
        label = f"{PHASE}_subtype"
    tdir = phase_table_dir(cfg, label)
    tdir.mkdir(parents=True, exist_ok=True)

    h5ad_path = (Path(cfg["results_dir"]) / "h5ad" /
                 ("08c_subclustered/%s.h5ad" % args.subcluster if args.subcluster
                  else "08_annotated/all_samples.h5ad"))

    adata = load_annotated_h5ad(cfg, args.subcluster)
    check_var_names_are_symbols(adata)
    if args.head:
        args.n_perms = 0
        adata = adata[:args.head].copy()
        print(f"  *** HEAD MODE: kept first {adata.n_obs:,} cells, n_perms forced to 0. "
              f"NUMBERS MEANINGLESS — code-path check only. ***")
    if args.node_scheme == "subtype" and not args.subcluster:
        print("\n  Building comms_subtype label (focal substates + broad-for-rest)...")
        adata = build_comms_subtype(adata, cfg, tissue, args.min_node_cells)
        celltype_key = "comms_subtype"
        print(f"  Cell-type key (node-scheme=subtype): '{celltype_key}'")
    else:
        celltype_key = resolve_celltype_key(adata, tissue, args.celltype_key, args.subcluster)
    if "lognorm" not in adata.layers:
        print("  No lognorm layer — computing via add_lognorm().")
        add_lognorm(adata)

    if args.smoke:
        args.n_perms = 0
        rng = np.random.default_rng(cfg.get("random_seed", 42))
        keep = []
        for _, idx in adata.obs.groupby(["group", "age"], observed=True).indices.items():
            take = idx if len(idx) <= args.smoke else rng.choice(idx, args.smoke, replace=False)
            keep.extend(take.tolist())
        adata = adata[sorted(keep)].copy()
        print(f"  *** SMOKE MODE: subsampled to {adata.n_obs:,} cells "
              f"(≤{args.smoke}/group×age), n_perms forced to 0. NUMBERS MEANINGLESS. ***")

    print(f"  {adata.n_obs:,} cells | groups={sorted(adata.obs['group'].astype(str).unique())} "
          f"| ages={sorted(adata.obs['age'].astype(str).unique())}")

    if "baseline" in arms:
        run_baseline(adata, cfg, args, tdir, celltype_key)
    if "differential" in arms:
        de_df = load_de_results(cfg, args.subcluster)
        run_differential(adata, de_df, cfg, args, tdir, celltype_key)
    if "perdonor" in arms:
        run_per_donor(h5ad_path, adata, cfg, args, tdir, celltype_key)

    print(f"\n{'='*60}\nPhase 8e compute complete. CSVs in {tdir}:")
    for f in sorted(tdir.glob("08e_*.csv")):
        print(f"  {f.name}")
    print(f"\n  Next: uv run python scripts/08e_comms_summary.py --config {args.config}"
          + (f" --subcluster {args.subcluster}" if args.subcluster else ""))
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
