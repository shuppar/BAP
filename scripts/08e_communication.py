#!/usr/bin/env python
"""
08e_communication.py — Phase 8e: Cell-cell communication (liana-py)

Entry point. Orchestrates three analytical arms:
  1. Baseline   — rank_aggregate per group×age (pooled cells per group)
  2. Differential — df_to_lr on 8b Wald stats (all 3 contrasts: ES-v-Rel,
                    LS-v-Rel, ES-v-LS)
  3. Per-donor  — rank_aggregate per donor → proper group statistics

Plot modules (leading underscore = not entry points):
  _08e_plots_baseline.py     — chord, network, dotplot, count heatmap, pathway
  _08e_plots_differential.py — volcano, diff dotplot, stress signature, delta heatmaps
  _08e_plots_perdonor.py     — stripplot, delta bar (quantified with stats)

Output structure:
  plots/08e_communication/
    01_overview/           — chord comparison, pathway heatmap, LR trajectory
    02_baseline_per_group/ — per group×age: chord, network, dotplots, heatmaps
    03_differential/       — per contrast×age + cross-contrast heatmap
    03_differential/delta_heatmaps/ — Δ activity score heatmaps (all group pairs)
    04_sender_receiver/    — bubble plot, Δ sender/receiver heatmaps (all pairs)
    05_per_donor/          — per age: donor stripplots + Δ bar with FDR

  tables/08e_communication/
    08e_lr_baseline.csv           — pooled LR scores per group×age
    08e_lr_differential.csv       — df_to_lr Wald stats (all 3 contrasts)
    08e_sender_receiver.csv       — cell-type sender/receiver scores
    08e_lr_per_donor.csv          — raw per-donor scores (offline audit)
    08e_lr_quantified.csv         — Δ activity + Mann-Whitney FDR (all pairs)
    08e_lr_stress_signature_pivot.csv — top LR × contrasts×ages matrix

Usage:
  uv run python scripts/08e_communication.py --config config/dev_split.yaml
  uv run python scripts/08e_communication.py --config config/brain.yaml \\
      --n-perms 1000 --top-lr 50 --top-lr-large 200

Deps: uv add liana networkx mpl-chord-diagram
"""

import argparse
import sys
import warnings
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from _utils import load_config, add_lognorm, phase_table_dir

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

PHASE = "08e_communication"

LABEL_KEY_PRIORITY = [
    "subcluster_name",       # Phase 7d (named subclusters from CellTypist+markers)
    "subcluster",            # Phase 7b (integer subclusters; fallback if 7d not run)
    "manual_annotation", "scanvi_celltype",
    "celltypist_majority", "provisional_celltype",
]


# ============================================================================
# CLI
# ============================================================================

def parse_args():
    p = argparse.ArgumentParser(description="Phase 8e: cell-cell communication")
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--top-lr", type=int, default=30,
                   help="Top N LR pairs for standard plots (default 30)")
    p.add_argument("--top-lr-large", type=int, default=100,
                   help="Top N for large supplementary dotplot (default 100)")
    p.add_argument("--n-perms", type=int, default=0,
                   help="CellPhoneDB permutations (0=skip; use 500-1000 on workstation)")
    p.add_argument("--min-cells", type=int, default=5,
                   help="Min cells per cell type (default 5)")
    p.add_argument("--expr-prop", type=float, default=0.1,
                   help="Min expression proportion per LR gene (default 0.1)")
    p.add_argument("--magnitude-cutoff", type=float, default=0.05,
                   help="magnitude_rank cutoff for active interactions (default 0.05)")
    p.add_argument("--celltype-key", default=None,
                   help="obs column for cell types (default: auto-detect)")
    p.add_argument("--skip-per-donor", action="store_true",
                   help="Skip per-donor arm (faster; loses quantification)")
    p.add_argument("--focus-celltypes", default=None,
                   help="Comma-separated cell types for focused heatmaps. "
                        "Default: union of YAML stress_focused_cell_types and "
                        "top-10 cell types by max |Δ| between groups.")
    p.add_argument("--zscore-rows", action="store_true",
                   help="Also produce row-Z-scored versions of clustered heatmaps "
                        "(pattern view; complements absolute Δ).")
    p.add_argument("--subcluster", default=None,
                   help="Run on 7b subcluster output instead of Phase 7 (use slug "
                        "matching 7b folder, e.g. 'excitatory_neurons').")
    return p.parse_args()


# ============================================================================
# Shared helpers (also imported by plot modules via sys.path)
# ============================================================================

def plot_dir(cfg: dict, subdir: str = "", phase_override: str = None) -> Path:
    phase = phase_override or PHASE
    d = Path(cfg["results_dir"]) / "plots" / phase / subdir
    d.mkdir(parents=True, exist_ok=True)
    return d


def slug(s: str) -> str:
    return s.replace(" ", "_").replace("/", "-").replace(".", "")


def resolve_celltype_key(adata, explicit=None):
    if explicit:
        if explicit not in adata.obs.columns:
            sys.exit(f"ERROR: --celltype-key '{explicit}' not in adata.obs.")
        return explicit
    for key in LABEL_KEY_PRIORITY:
        if key in adata.obs.columns:
            if key == "manual_annotation" and adata.obs[key].astype(str).eq("").all():
                continue
            print(f"  Cell-type column: '{key}'")
            return key
    sys.exit(
        f"ERROR: no cell-type label column found. Looked for: {LABEL_KEY_PRIORITY}\n"
        "  Run Phase 7 first or pass --celltype-key."
    )


def check_var_names_are_symbols(adata):
    if any(v.startswith("ENSMUSG") for v in list(adata.var_names[:5])):
        sys.exit("ERROR: var_names look like Ensembl IDs — liana needs gene symbols.")


def load_annotated_h5ad(cfg, subcluster_slug=None):
    base = Path(cfg["results_dir"]) / "h5ad"
    if subcluster_slug:
        p = base / "08c_subclustered" / f"{subcluster_slug}.h5ad"
        if not p.is_file():
            sys.exit(f"ERROR: subcluster h5ad not found: {p}\n"
                     "  Run scripts/07b_subcluster.py first.")
        print(f"  Loading subcluster: {p}")
        return sc.read_h5ad(p)
    candidates = [
        base / "08b_label_transferred" / "all_samples.h5ad",
        base / "08_annotated" / "all_samples.h5ad",
        base / "08d_trajectory" / "all_samples.h5ad",
    ]
    for p in candidates:
        if p.is_file():
            print(f"  Loading: {p}")
            return sc.read_h5ad(p)
    sys.exit(
        "ERROR: annotated h5ad not found. Tried:\n"
        + "\n".join(f"  {p}" for p in candidates)
    )


def resolve_focus_celltypes(baseline_df, cfg, args, celltype_key, adata):
    """Build the focus celltype set.

    Default = union of:
      1. YAML stress_focused_cell_types (prior knowledge)
      2. Top 10 cell types by max |Δ activity| across all group pairs and ages

    User can override entirely via --focus-celltypes.
    Returns sorted list of strings matching actual labels in baseline_df.
    """
    # User explicit override → just split and validate
    if args.focus_celltypes:
        wanted = [c.strip() for c in args.focus_celltypes.split(",")]
        observed = set(adata.obs[celltype_key].astype(str).unique())
        kept = [c for c in wanted if c in observed]
        missing = [c for c in wanted if c not in observed]
        if missing:
            print(f"  [warn] --focus-celltypes: not in data, ignoring: {missing}")
        print(f"  Focus cell types (user-specified): {kept}")
        return sorted(set(kept))

    # Auto: top-10 by max |Δ| + YAML stress list
    from itertools import combinations
    observed = set(adata.obs[celltype_key].astype(str).unique())

    yaml_focus_raw = cfg.get("stress_focused_cell_types", []) or []
    # YAML uses snake_case; match against actual labels case-insensitively
    yaml_focus = set()
    for label in observed:
        norm_label = label.lower().replace(" ", "_").replace("-", "_")
        for y in yaml_focus_raw:
            y_norm = y.lower().replace(" ", "_").replace("-", "_")
            if y_norm in norm_label or norm_label in y_norm:
                yaml_focus.add(label)
                break

    # Data-driven: per cell type, find max |Δ activity| across all group pairs
    auto_focus = set()
    if not baseline_df.empty and "magnitude_rank" in baseline_df.columns:
        df = baseline_df.copy()
        df["activity"] = 1 - df["magnitude_rank"]
        df["ct_pair"] = df["source"] + "→" + df["target"]

        groups = sorted(df["group"].unique())
        max_abs_per_ct = {}
        for grp_a, grp_b in combinations(groups, 2):
            for age in df["age"].unique():
                mean_a = (df[(df["group"] == grp_a) & (df["age"] == age)]
                          .groupby(["source", "target"])["activity"].mean())
                mean_b = (df[(df["group"] == grp_b) & (df["age"] == age)]
                          .groupby(["source", "target"])["activity"].mean())
                delta = mean_a.subtract(mean_b, fill_value=0).abs()
                # Attribute |Δ| to BOTH source and target (cell type appears in interactions)
                for (src, tgt), v in delta.items():
                    max_abs_per_ct[src] = max(max_abs_per_ct.get(src, 0), v)
                    max_abs_per_ct[tgt] = max(max_abs_per_ct.get(tgt, 0), v)
        top10 = sorted(max_abs_per_ct.items(), key=lambda x: -x[1])[:10]
        auto_focus = {ct for ct, _ in top10}

    focus = sorted((yaml_focus | auto_focus) & observed)
    yaml_in = sorted(yaml_focus & observed)
    auto_in = sorted(auto_focus & observed)
    print(f"  Focus cell types ({len(focus)} total):")
    print(f"    from YAML stress_focused_cell_types ({len(yaml_in)}): {yaml_in}")
    print(f"    from top-10 by max |Δ| ({len(auto_in)}): {auto_in}")
    return focus


def load_de_results(cfg, subcluster_slug=None):
    """Load 08b DE results table. In subcluster mode, loads the matching
    08b_de_results_subcluster_{slug}.csv written by 08b's --subcluster mode.
    """
    tdir = Path(cfg["results_dir"]) / "tables" / "08b_de"
    if subcluster_slug:
        p = tdir / f"08b_de_results_subcluster_{subcluster_slug}.csv"
        if not p.is_file():
            print(f"  [warn] subcluster DE results not found: {p.name}")
            print(f"         Run: scripts/08b_de.py --subcluster {subcluster_slug}")
            print(f"         Differential CCC will be skipped.")
            return pd.DataFrame()
    else:
        p = tdir / "08b_de_results.csv"
        if not p.is_file():
            print("  [warn] 08b DE results not found — differential CCC skipped.")
            return pd.DataFrame()
    df = pd.read_csv(p)
    print(f"  Loaded DE results: {len(df):,} rows from {p.name}")
    return df


def compute_sender_receiver(baseline_df: pd.DataFrame) -> pd.DataFrame:
    if baseline_df.empty:
        return pd.DataFrame()
    records = []
    for (group, age), grp in baseline_df.groupby(["group", "age"]):
        active = grp[grp["magnitude_rank"] < 0.05] if "magnitude_rank" in grp.columns else grp
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
# Part 1: Baseline — rank_aggregate per group×age
# ============================================================================

def run_baseline(adata, cfg, args, tdir, celltype_key):
    import liana as li
    print("\n[Part 1] Baseline signaling (rank_aggregate per group×age)")
    all_results = []

    for age in sorted(adata.obs["age"].unique()):
        for group in sorted(adata.obs["group"].unique()):
            mask = (adata.obs["age"] == age) & (adata.obs["group"] == group)
            if mask.sum() < args.min_cells * 2:
                print(f"  SKIP {group}/{age}: {mask.sum()} cells")
                continue
            sub = adata[mask].copy()
            valid_cts = sub.obs[celltype_key].value_counts()
            valid_cts = valid_cts[valid_cts >= args.min_cells].index
            if len(valid_cts) < 2:
                print(f"  SKIP {group}/{age}: <2 valid cell types")
                continue
            sub = sub[sub.obs[celltype_key].isin(valid_cts)].copy()
            print(f"  rank_aggregate: {group}/{age} "
                  f"({sub.n_obs:,} cells, {len(valid_cts)} types)")
            n_perms = args.n_perms if args.n_perms > 0 else None
            try:
                li.mt.rank_aggregate(
                    sub, groupby=celltype_key,
                    resource_name="mouseconsensus",
                    expr_prop=args.expr_prop, min_cells=args.min_cells,
                    use_raw=False, layer="lognorm",
                    n_perms=n_perms, seed=cfg.get("random_seed", 42),
                    n_jobs=1, verbose=False, inplace=True,
                )
            except Exception as e:
                print(f"  [warn] {group}/{age} failed: {e} — skipping")
                continue
            res = sub.uns["liana_res"].copy()
            res["group"] = group
            res["age"] = age
            all_results.append(res)

    if not all_results:
        print("  [warn] No baseline results.")
        return pd.DataFrame()

    combined = pd.concat(all_results, ignore_index=True)
    out = tdir / "08e_lr_baseline.csv"
    combined.to_csv(out, index=False)
    print(f"  Saved: {out.name}  ({len(combined):,} rows)")
    return combined


# ============================================================================
# Part 2: Differential — df_to_lr on 8b Wald stats (all 3 contrasts)
# ============================================================================

def run_differential(adata, de_df, cfg, args, tdir, celltype_key):
    import liana as li
    from statsmodels.stats.multitest import multipletests

    print("\n[Part 2] Differential CCC (df_to_lr — ES-v-Rel, LS-v-Rel, ES-v-LS)")

    if de_df.empty:
        print("  Skipped — no DE results.")
        return pd.DataFrame()

    # 08b actual schema:
    #   contrast (not contrast_name), celltype (not cell_type),
    #   log2FC (not log2FoldChange), group_level (encodes age as "age-4W"),
    #   flag, gene, stat, padj, pvalue, direction, ...
    required = {"contrast", "celltype", "group_level", "gene", "stat", "padj",
                "log2FC", "flag"}
    missing = required - set(de_df.columns)
    if missing:
        print(f"  [warn] DE table missing columns: {missing}. Skipping.")
        print(f"         Available: {sorted(de_df.columns)}")
        return pd.DataFrame()

    # Filter to gene-level rows (some rows may be DESeq2-failure placeholders)
    de_df = de_df.dropna(subset=["gene", "stat"]).copy()
    if de_df.empty:
        print("  [warn] No gene-level rows in DE table.")
        return pd.DataFrame()

    # Extract age from group_level like "age-4W" (or "age-4W_sex-M" for stratified)
    def _extract_age(gl):
        if not isinstance(gl, str):
            return None
        for part in gl.split("_"):
            if part.startswith("age-"):
                return part[len("age-"):]
        return None
    de_df["age"] = de_df["group_level"].map(_extract_age)
    # Drop rows where age couldn't be extracted (e.g. across-age contrasts)
    de_df = de_df.dropna(subset=["age"])
    if de_df.empty:
        print("  [warn] No DE rows with parseable age. Skipping.")
        return pd.DataFrame()

    # primary = ES-v-Rel, LS-v-Rel; secondary = ES-v-LS
    use = de_df[de_df["flag"].isin(["primary", "secondary"])].copy()
    if use.empty:
        print("  [warn] No primary/secondary contrasts found.")
        return pd.DataFrame()

    print(f"  Contrasts found: {sorted(use['contrast'].unique())}")
    print(f"  Ages found:      {sorted(use['age'].unique())}")

    all_lr = []
    for contrast_name in sorted(use["contrast"].unique()):
        for age in sorted(use["age"].unique()):
            sub_de = use[(use["contrast"] == contrast_name) &
                         (use["age"] == age)]
            if sub_de.empty:
                continue
            mask = adata.obs["age"] == age
            if mask.sum() < args.min_cells * 2:
                continue
            sub_adata = adata[mask].copy()

            dea_df = sub_de[["celltype", "gene", "stat", "padj", "log2FC"]].rename(
                columns={"log2FC": "log2fc", "celltype": celltype_key}
            )
            print(f"  df_to_lr: {contrast_name}/{age} "
                  f"({dea_df[celltype_key].nunique()} cell types)")
            try:
                lr_res = li.multi.df_to_lr(
                    sub_adata, dea_df=dea_df,
                    resource_name="mouseconsensus",
                    expr_prop=args.expr_prop, groupby=celltype_key,
                    stat_keys=["stat", "padj", "log2fc"],
                    use_raw=False, layer="lognorm",
                    complex_col="stat", verbose=False, return_all_lrs=False,
                )
            except Exception as e:
                print(f"  [warn] df_to_lr failed {contrast_name}/{age}: {e}")
                continue
            lr_res["contrast_name"] = contrast_name
            lr_res["age"] = age
            all_lr.append(lr_res)

    if not all_lr:
        print("  [warn] No differential LR results.")
        return pd.DataFrame()

    combined = pd.concat(all_lr, ignore_index=True)

    pval_col = next((c for c in combined.columns
                     if "pvalue" in c.lower()), None)
    if pval_col:
        fdr_list = []
        for _, grp in combined.groupby(["contrast_name", "age", "source", "target"]):
            _, fdr, _, _ = multipletests(grp[pval_col].fillna(1.0), method="fdr_bh")
            fdr_list.append(pd.Series(fdr, index=grp.index))
        combined["interaction_fdr"] = pd.concat(fdr_list)
    else:
        combined["interaction_rank_within_slice"] = combined.groupby(
            ["contrast_name", "age", "source", "target"]
        )["interaction_stat"].rank(ascending=False)

    out = tdir / "08e_lr_differential.csv"
    combined.to_csv(out, index=False)
    n_contrasts = combined["contrast_name"].nunique()
    print(f"  Saved: {out.name}  ({len(combined):,} rows, {n_contrasts} contrasts)")
    return combined


# ============================================================================
# Part 3: Per-donor — rank_aggregate per donor → all pairwise group stats
# ============================================================================

def run_per_donor(adata, cfg, args, tdir, celltype_key):
    import liana as li
    from scipy import stats as scipy_stats
    from statsmodels.stats.multitest import multipletests

    print("\n[Part 3] Per-donor LR scores (quantified — all group pairs)")
    donors = sorted(adata.obs["donor_id"].unique())
    print(f"  {len(donors)} donors")
    all_results = []

    for donor in donors:
        mask = adata.obs["donor_id"] == donor
        sub = adata[mask].copy()
        meta = adata.obs.loc[mask, ["group", "age"]].iloc[0]
        group, age = meta["group"], meta["age"]

        valid_cts = sub.obs[celltype_key].value_counts()
        valid_cts = valid_cts[valid_cts >= args.min_cells].index
        if len(valid_cts) < 2:
            print(f"    SKIP {donor}: <2 cell types")
            continue
        sub = sub[sub.obs[celltype_key].isin(valid_cts)].copy()

        n_perms = args.n_perms if args.n_perms > 0 else None
        try:
            li.mt.rank_aggregate(
                sub, groupby=celltype_key,
                resource_name="mouseconsensus",
                expr_prop=args.expr_prop, min_cells=args.min_cells,
                use_raw=False, layer="lognorm",
                n_perms=n_perms, seed=cfg.get("random_seed", 42),
                n_jobs=1, verbose=False, inplace=True,
            )
        except Exception as e:
            print(f"    SKIP {donor}: {e}")
            continue

        res = sub.uns["liana_res"].copy()
        res["donor_id"] = donor
        res["group"] = group
        res["age"] = age
        all_results.append(res)
        print(f"    {donor}: {group}/{age}  {len(res):,} LR pairs")

    if not all_results:
        print("  [warn] No per-donor results.")
        return pd.DataFrame(), pd.DataFrame()

    per_donor_df = pd.concat(all_results, ignore_index=True)
    per_donor_df["activity"] = 1 - per_donor_df["magnitude_rank"]
    per_donor_df["lr_pair"] = (per_donor_df["ligand_complex"] + "→"
                                + per_donor_df["receptor_complex"])
    per_donor_df["ct_pair"] = per_donor_df["source"] + "→" + per_donor_df["target"]
    per_donor_df.to_csv(tdir / "08e_lr_per_donor.csv", index=False)
    print(f"  Saved: 08e_lr_per_donor.csv  ({len(per_donor_df):,} rows)")

    # All pairwise group comparisons per LR×ct_pair×age
    groups = sorted(per_donor_df["group"].unique())
    ref_group = cfg.get("group_reference", "Relaxed")
    group_pairs = []
    for a, b in combinations(groups, 2):
        if b == ref_group:
            group_pairs.append((a, b))
        elif a == ref_group:
            group_pairs.append((b, a))
        else:
            group_pairs.append((a, b))

    records = []
    for (lr, ct, age), grp in per_donor_df.groupby(["lr_pair", "ct_pair", "age"]):
        for test_grp, ctrl_grp in group_pairs:
            a_sc = grp.loc[grp["group"] == test_grp, "activity"].values
            b_sc = grp.loc[grp["group"] == ctrl_grp, "activity"].values
            if len(a_sc) < 2 or len(b_sc) < 2:
                continue
            stat, pval = scipy_stats.mannwhitneyu(a_sc, b_sc, alternative="two-sided")
            records.append({
                "lr_pair": lr, "ct_pair": ct, "age": age,
                "test_group": test_grp, "ctrl_group": ctrl_grp,
                "contrast": f"{test_grp}_vs_{ctrl_grp}",
                "mean_test": float(a_sc.mean()),
                "mean_ctrl": float(b_sc.mean()),
                "delta_activity": float(a_sc.mean() - b_sc.mean()),
                "n_test": len(a_sc), "n_ctrl": len(b_sc),
                "mannwhitney_stat": float(stat), "pvalue": float(pval),
            })

    if not records:
        print("  [warn] No group comparisons (too few donors per group).")
        return per_donor_df, pd.DataFrame()

    quant_df = pd.DataFrame(records)
    fdr_list = []
    for _, sl in quant_df.groupby(["age", "contrast"]):
        _, fdr, _, _ = multipletests(sl["pvalue"].fillna(1.0), method="fdr_bh")
        fdr_list.append(pd.Series(fdr, index=sl.index))
    quant_df["fdr"] = pd.concat(fdr_list)
    quant_df["significant"] = quant_df["fdr"] < 0.05
    quant_df.to_csv(tdir / "08e_lr_quantified.csv", index=False)
    n_sig = quant_df["significant"].sum()
    print(f"  Saved: 08e_lr_quantified.csv  ({len(quant_df):,} tests, {n_sig} FDR<0.05)")
    print(f"  Contrasts in quantification: {sorted(quant_df['contrast'].unique())}")
    return per_donor_df, quant_df


# ============================================================================
# Main
# ============================================================================

def main():
    args = parse_args()
    cfg = load_config(args.config)
    tissue = cfg.get("tissue", "unknown")
    ref_group = cfg.get("group_reference", "Relaxed")

    print(f"\n{'='*60}")
    print(f"Phase 8e: Cell-cell communication  [{tissue}]")
    print(f"{'='*60}")

    try:
        import liana as li
        print(f"  liana {li.__version__}")
    except ImportError:
        sys.exit("ERROR: liana not installed. Run: uv add liana")

    has_networkx = False
    try:
        import networkx  # noqa: F401
        has_networkx = True
    except ImportError:
        print("  [warn] networkx not installed — network plots skipped (uv add networkx)")

    has_chord = False
    try:
        from mpl_chord_diagram import chord_diagram  # noqa: F401
        has_chord = True
    except ImportError:
        print("  [warn] mpl-chord-diagram not installed — chord plots skipped "
              "(uv add mpl-chord-diagram)")

    from _08e_plots_baseline import (
        plot_chord_diagram, plot_delta_chord_diagram,
        plot_network_graph, plot_interaction_count_heatmap,
        plot_baseline_dotplot, plot_large_lr_dotplot,
        plot_top_lr_per_celltype_pair, plot_interaction_counts_barplot,
        plot_pathway_activity_heatmap, plot_lr_persistence_across_ages,
    )
    from _08e_plots_differential import (
        plot_differential_dotplot, plot_differential_volcano,
        plot_stress_signature_heatmap, plot_delta_lr_heatmap,
        plot_sender_receiver_bubble, plot_sender_receiver_heatmap,
        plot_delta_sender_receiver_heatmap,
        plot_rank_rank_scatter, get_pathway_map,
    )
    from _08e_plots_perdonor import plot_per_donor_quantification

    # Subcluster runs go to a separate output tree to avoid overwriting the
    # main run's plots and tables. Phase 7b/08b use the same convention.
    phase_name = PHASE + (f"_subcluster_{args.subcluster}" if args.subcluster else "")
    tdir = phase_table_dir(cfg, phase_name)

    # Local helper bound to phase_name (overrides PHASE default)
    def _pdir(subdir=""):
        return plot_dir(cfg, subdir, phase_override=phase_name)

    print("\n[0] Loading data...")
    adata = load_annotated_h5ad(cfg, subcluster_slug=args.subcluster)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    for col in ("group", "age", "donor_id"):
        if col not in adata.obs.columns:
            sys.exit(f"ERROR: adata.obs missing '{col}'")
    celltype_key = resolve_celltype_key(adata, args.celltype_key)
    check_var_names_are_symbols(adata)
    if "lognorm" not in adata.layers:
        print("  Computing lognorm layer...")
        add_lognorm(adata)
    from scipy.sparse import issparse
    if not issparse(adata.X):
        print("  [warn] adata.X is dense.")

    print("\n[0b] Loading 08b DE results...")
    de_df = load_de_results(cfg, subcluster_slug=args.subcluster)

    # --- Run ---
    baseline_df = run_baseline(adata, cfg, args, tdir, celltype_key)
    diff_df = run_differential(adata, de_df, cfg, args, tdir, celltype_key)

    sr_df = compute_sender_receiver(baseline_df)
    if not sr_df.empty:
        sr_df.to_csv(tdir / "08e_sender_receiver.csv", index=False)
        print("  Saved: 08e_sender_receiver.csv")

    if args.skip_per_donor:
        per_donor_df, quant_df = pd.DataFrame(), pd.DataFrame()
        print("\n[Part 3] Per-donor skipped (--skip-per-donor)")
    else:
        per_donor_df, quant_df = run_per_donor(adata, cfg, args, tdir, celltype_key)

    # --- Build group pair list (all three comparisons) ---
    groups = sorted(baseline_df["group"].unique()) if not baseline_df.empty else []
    ages   = sorted(baseline_df["age"].unique())   if not baseline_df.empty else []
    all_group_pairs = []
    for a, b in combinations(groups, 2):
        if b == ref_group:
            all_group_pairs.append((a, b))
        elif a == ref_group:
            all_group_pairs.append((b, a))
        else:
            all_group_pairs.append((a, b))

    # --- Plots ---
    print("\n[Plots] Generating figures...")

    # 01_overview
    pdir_ov = _pdir("01_overview")
    if not baseline_df.empty:
        plot_interaction_counts_barplot(baseline_df, ref_group, pdir_ov)
        plot_pathway_activity_heatmap(baseline_df, args.magnitude_cutoff, pdir_ov)
        plot_lr_persistence_across_ages(baseline_df, args.magnitude_cutoff,
                                        args.top_lr, pdir_ov)
        for age in ages:
            if has_chord:
                plot_delta_chord_diagram(baseline_df, age,
                                         args.magnitude_cutoff, pdir_ov)

    # 02_baseline_per_group
    if not baseline_df.empty:
        for group in groups:
            for age in ages:
                pdir_grp = _pdir(f"02_baseline_per_group/{slug(group)}_{slug(age)}")
                plot_baseline_dotplot(baseline_df, group, age,
                                      args.top_lr, pdir_grp)
                plot_large_lr_dotplot(baseline_df, group, age,
                                      args.top_lr_large, pdir_grp)
                plot_interaction_count_heatmap(baseline_df, group, age,
                                               args.magnitude_cutoff, pdir_grp)
                plot_top_lr_per_celltype_pair(baseline_df, group, age,
                                              args.top_lr, pdir_grp)
                if has_networkx:
                    plot_network_graph(baseline_df, group, age,
                                       args.magnitude_cutoff, pdir_grp)
                if has_chord:
                    plot_chord_diagram(baseline_df, group, age,
                                       args.magnitude_cutoff, pdir_grp)

    # 03_differential — per contrast×age
    if not diff_df.empty:
        for contrast_name in sorted(diff_df["contrast_name"].unique()):
            for age in sorted(diff_df["age"].unique()):
                sub = diff_df[(diff_df["contrast_name"] == contrast_name) &
                              (diff_df["age"] == age)]
                if sub.empty:
                    continue
                pdir_d = _pdir(f"03_differential/{slug(contrast_name)}_{slug(age)}")
                plot_differential_dotplot(diff_df, contrast_name, age,
                                          args.top_lr, pdir_d)
                plot_differential_volcano(diff_df, contrast_name, age,
                                          args.top_lr, pdir_d)
        plot_stress_signature_heatmap(diff_df, tdir,
                                      _pdir("03_differential"),
                                      top_n=args.top_lr_large)

    # 03_differential/delta_heatmaps — clustered Δ activity heatmaps + rank-rank
    if not baseline_df.empty:
        pdir_delta = _pdir("03_differential/delta_heatmaps")

        # Resolve focus celltype set (auto = YAML stress list ∪ top-10 by |Δ|)
        focus_cts = resolve_focus_celltypes(baseline_df, cfg, args,
                                            celltype_key, adata)

        # Fetch pathway annotations once (used as row colour bar)
        print("  Fetching pathway annotations for row colour bars...")
        pathway_map = get_pathway_map()
        if pathway_map.empty:
            print("    (no pathway map available — row colour bar will be skipped)")

        for age in ages:
            for grp_a, grp_b in all_group_pairs:
                # Full unfocused, absolute Δ
                plot_delta_lr_heatmap(
                    baseline_df, grp_a, grp_b, age,
                    args.magnitude_cutoff, pdir_delta,
                    focus_celltypes=None, zscore_rows=False,
                    pathway_map=pathway_map,
                )
                # Focused subset (always — focus_cts always non-empty by design)
                if focus_cts:
                    plot_delta_lr_heatmap(
                        baseline_df, grp_a, grp_b, age,
                        args.magnitude_cutoff, pdir_delta,
                        focus_celltypes=focus_cts, zscore_rows=False,
                        pathway_map=pathway_map,
                    )
                # Optional z-scored versions (pattern view)
                if args.zscore_rows:
                    plot_delta_lr_heatmap(
                        baseline_df, grp_a, grp_b, age,
                        args.magnitude_cutoff, pdir_delta,
                        focus_celltypes=None, zscore_rows=True,
                        pathway_map=pathway_map,
                    )
                    if focus_cts:
                        plot_delta_lr_heatmap(
                            baseline_df, grp_a, grp_b, age,
                            args.magnitude_cutoff, pdir_delta,
                            focus_celltypes=focus_cts, zscore_rows=True,
                            pathway_map=pathway_map,
                        )

        # Rank-rank scatters: pairwise concordance between Δ-signatures.
        # The central biological question: do ES and LS hit the same programs?
        # Pivotal comparisons (in this order):
        #   ES-v-Rel  vs  LS-v-Rel   — main signature concordance
        #   ES-v-Rel  vs  ES-v-LS    — does ES signal also differ from LS?
        #   LS-v-Rel  vs  ES-v-LS    — same for LS
        pdir_rr = _pdir("03_differential/rank_rank")
        # Group contrasts as ordered (test, ctrl) tuples
        contrast_tuples = list(all_group_pairs)
        if len(contrast_tuples) >= 2:
            for age in ages:
                for ca, cb in combinations(contrast_tuples, 2):
                    plot_rank_rank_scatter(
                        baseline_df, ca, cb, age,
                        args.magnitude_cutoff, pdir_rr,
                        focus_celltypes=None,
                    )
                    if focus_cts:
                        plot_rank_rank_scatter(
                            baseline_df, ca, cb, age,
                            args.magnitude_cutoff, pdir_rr,
                            focus_celltypes=focus_cts,
                        )

    # 04_sender_receiver
    pdir_sr = _pdir("04_sender_receiver")
    if not baseline_df.empty:
        for age in ages:
            plot_sender_receiver_bubble(baseline_df, age,
                                        args.magnitude_cutoff, pdir_sr)
        if not sr_df.empty:
            for age in ages:
                plot_sender_receiver_heatmap(sr_df, age, ref_group, pdir_sr)
            for age in ages:
                for grp_a, grp_b in all_group_pairs:
                    plot_delta_sender_receiver_heatmap(
                        sr_df, grp_a, grp_b, age, pdir_sr)

    # 05_per_donor
    if not per_donor_df.empty:
        for age in sorted(per_donor_df["age"].unique()):
            pdir_pd = _pdir(f"05_per_donor/{slug(age)}")
            plot_per_donor_quantification(
                per_donor_df[per_donor_df["age"] == age],
                quant_df[quant_df["age"] == age] if not quant_df.empty else pd.DataFrame(),
                ref_group, args.top_lr, pdir_pd,
            )

    # --- Summary ---
    print(f"\n{'='*60}")
    print("Phase 8e complete.")
    print(f"\n  Output folders (under plots/{phase_name}/):")
    print(f"    01_overview/")
    print(f"    02_baseline_per_group/{{group}}_{{age}}/")
    print(f"    03_differential/{{contrast}}_{{age}}/")
    print(f"    03_differential/delta_heatmaps/   ← clustered Δ heatmaps "
          f"(absolute + focused; pathway+celltype colour bars)")
    print(f"    03_differential/rank_rank/        ← LR signature concordance "
          f"scatters (Spearman ρ)")
    print(f"    04_sender_receiver/")
    print(f"    05_per_donor/{{age}}/")
    print(f"\n  Offline CSVs:")
    for f in sorted(tdir.glob("08e_*.csv")):
        print(f"    {f.name}")
    if not quant_df.empty:
        n_sig = quant_df["significant"].sum()
        contrasts = sorted(quant_df["contrast"].unique())
        print(f"\n  Per-donor stats: {n_sig} FDR<0.05 across {contrasts}")
    print()


if __name__ == "__main__":
    main()
