#!/usr/bin/env python
"""
07_annotation.py — Phase 7: cell type annotation.

Two-track approach:
  1. Reference-based: CellTypist label transfer
  2. Marker-based: rank_genes_groups per Leiden cluster + curated marker dotplot

All marker scoring and plotting runs on lognorm (recomputed at load time —
it was dropped after Phase 5 to save disk).

Composition plots (diagnostic, not the final biological result):
  - cluster_composition_by_sample.png : Leiden cluster × sample stacked bar
    → catches single-sample clusters (batch artifacts)
  - celltype_composition_by_sample.png: cell type × sample stacked bar
    → uses CellTypist majority label if available, Leiden number otherwise
    → first look at whether stress groups differ in composition
  - celltype_composition_by_group.png : same data, grouped by condition
    → Early/Late/Relaxed side by side per cell type

The quantitative composition analysis with statistics (scCODA + propeller)
is in Phase 8a.

Usage:
  uv run python scripts/07_annotation.py --config config/dev.yaml
  uv run python scripts/07_annotation.py --config config/brain.yaml

Inputs:
  {results_dir}/h5ad/07_clustered/all_samples.h5ad  (from Phase 6)

Outputs:
  {results_dir}/h5ad/08_annotated/all_samples.h5ad
  {results_dir}/plots/07_annotation/
    - umap_leiden_for_annotation.png
    - marker_dotplot.png
    - marker_heatmap_top10.png
    - umap_marker_scores.png
    - umap_celltypist.png                  (if model configured)
    - umap_celltypist_confidence.png       (if model configured)
    - cluster_composition_by_sample.png    : Leiden × sample
    - celltype_composition_by_sample.png   : cell type × sample
    - celltype_composition_by_group.png    : cell type × group
  {results_dir}/tables/
    - marker_genes_per_cluster.csv
    - annotation_summary.csv
    - celltype_composition.csv
    - celltypist_predictions.csv           (if model configured)
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from _utils import load_config, add_lognorm


# ---------------------------------------------------------------------------
# Built-in marker sets — fallback if not in YAML.
# PROVENANCE NOTE: consensus markers from the mouse snRNA-seq literature
# (Allen BCA, Di Bella 2021, Marsh & Blelloch 2020) assembled from training
# knowledge. For production: replace with gene lists from the supplementary
# tables of Di Bella et al. 2021 (brain) and Marsh & Blelloch 2020 (placenta).
# ---------------------------------------------------------------------------

BRAIN_MARKERS = {
    "Excitatory neurons":  ["Slc17a7", "Slc17a6", "Neurod2", "Neurod6", "Tbr1"],
    "Inhibitory neurons":  ["Gad1", "Gad2", "Slc32a1", "Lhx6", "Sst", "Pvalb", "Vip"],
    "Astrocytes":          ["Aqp4", "Gfap", "Aldh1l1", "Sox9", "Slc1a3"],
    "Oligodendrocytes":    ["Mbp", "Plp1", "Mog", "Cnp", "Mal"],
    "OPC":                 ["Pdgfra", "Cspg4", "Olig2", "Sox10"],
    "Microglia":           ["Cx3cr1", "P2ry12", "Tmem119", "Hexb"],
    "Endothelial":         ["Cldn5", "Pecam1", "Ly6c1"],
    "Pericytes":           ["Pdgfrb", "Rgs5", "Vtn"],
    "Radial glia / NPCs":  ["Nes", "Sox2", "Pax6", "Vim", "Fabp7"],
    "Choroid plexus":      ["Ttr", "Folr1", "Clic6"],
}

PLACENTA_MARKERS = {
    "Trophoblast (SpT)":     ["Prl3d1", "Prl8a8", "Tpbpa"],
    "Trophoblast (LaT)":     ["Gcm1", "Syna", "Synb"],
    "Trophoblast (TGC)":     ["Prl3b1", "Hand1"],
    "Decidual stromal":      ["Igfbp1", "Foxo1"],
    "Endothelial":           ["Cldn5", "Pecam1"],
    "Hofbauer / Macrophage": ["Cd68", "Adgre1", "Cx3cr1"],
    "NK cells":              ["Ncr1", "Klrb1c", "Gzma"],
    "Erythroblasts":         ["Hbb-bt", "Hba-a1"],
}


def get_markers(cfg: dict) -> dict:
    custom = cfg.get("annotation", {}).get("markers")
    if custom:
        return custom
    return BRAIN_MARKERS if cfg["tissue"] == "brain" else PLACENTA_MARKERS


# ---------------------------------------------------------------------------
# Track 1: CellTypist
# ---------------------------------------------------------------------------

def run_celltypist(adata, model_source: str):
    """Run CellTypist label transfer. Returns per-cell DataFrame or None."""
    try:
        import celltypist
        from celltypist import models
    except ImportError:
        print("  [skip] celltypist not installed. Run: uv add celltypist")
        return None

    tmp = adata.copy()
    tmp.X = tmp.layers["lognorm"].copy()

    try:
        if Path(model_source).is_file():
            model = models.Model.load(model_source)
            print(f"  Loaded local model: {model_source}")
        else:
            print(f"  Downloading/loading built-in model: {model_source}")
            models.download_models(model=model_source, force_update=False)
            model = models.Model.load(model_source)
    except Exception as e:
        print(f"  [warn] Could not load CellTypist model '{model_source}': {e}")
        return None

    predictions = celltypist.annotate(tmp, model=model, majority_voting=True)
    result = predictions.predicted_labels.copy()
    result["conf_score"] = predictions.probability_matrix.max(axis=1).values
    return result


# ---------------------------------------------------------------------------
# Track 2: marker-based
# ---------------------------------------------------------------------------

def run_marker_genes(adata, obs_key: str = "leiden") -> pd.DataFrame:
    """Wilcoxon rank_genes_groups on lognorm. Returns top-20 per cluster."""
    sc.tl.rank_genes_groups(
        adata, groupby=obs_key, method="wilcoxon",
        layer="lognorm", use_raw=False,
        key_added="rank_genes_groups",
    )
    result = sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_groups")
    return (result.sort_values("scores", ascending=False)
                  .groupby("group").head(20)
                  .reset_index(drop=True))


def score_marker_sets(adata, markers: dict) -> None:
    """score_genes on lognorm for each curated cell type marker set."""
    for ct, genes in markers.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            continue
        key = ("score_" + ct
               .replace(" ", "_").replace("/", "_")
               .replace("(", "").replace(")", ""))
        sc.tl.score_genes(adata, present, score_name=key, use_raw=False)


# ---------------------------------------------------------------------------
# Plots: marker / UMAP
# ---------------------------------------------------------------------------

def plot_leiden_for_reference(adata, obs_key: str, out: Path) -> None:
    n = adata.obs[obs_key].nunique()
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(adata, color=obs_key, ax=ax, show=False, frameon=False,
               legend_loc="on data", legend_fontsize=7, size=6,
               title=f"Leiden clusters ({n}) — for annotation reconciliation")
    fig.tight_layout()
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_umap_celltypist(adata, out_label: Path, out_conf: Path) -> None:
    if "celltypist_majority" not in adata.obs.columns:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color="celltypist_majority", ax=ax, show=False,
               frameon=False, legend_fontsize=7, size=6,
               title="CellTypist majority-vote labels")
    fig.tight_layout()
    fig.savefig(out_label, dpi=140, bbox_inches="tight")
    plt.close(fig)

    if "celltypist_conf_score" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(7, 5))
        sc.pl.umap(adata, color="celltypist_conf_score", ax=ax, show=False,
                   frameon=False, color_map="viridis", size=6,
                   title="CellTypist confidence score")
        fig.tight_layout()
        fig.savefig(out_conf, dpi=140, bbox_inches="tight")
        plt.close(fig)


def plot_marker_dotplot(adata, markers: dict, obs_key: str, out: Path) -> None:
    markers_present = {ct: [g for g in genes if g in adata.var_names]
                       for ct, genes in markers.items()}
    markers_present = {ct: g for ct, g in markers_present.items() if g}
    if not markers_present:
        print("  [skip] marker_dotplot: no curated marker genes in adata.var_names")
        return
    seen, gene_list = set(), []
    for genes in markers_present.values():
        for g in genes:
            if g not in seen:
                gene_list.append(g)
                seen.add(g)
    adata.obs[obs_key] = adata.obs[obs_key].astype("category")
    fig = sc.pl.dotplot(adata, var_names=gene_list, groupby=obs_key,
                        layer="lognorm", show=False, return_fig=True,
                        title=f"Curated markers × {obs_key}")
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close()


def plot_top_marker_heatmap(adata, top_markers: pd.DataFrame,
                             obs_key: str, out: Path, n: int = 10) -> None:
    top_genes = (top_markers.groupby("group")
                             .apply(lambda x: x.nlargest(n, "scores"))
                             .reset_index(drop=True)["names"]
                             .unique().tolist())
    top_genes = [g for g in top_genes if g in adata.var_names]
    if not top_genes:
        print("  [skip] marker_heatmap: no genes to plot")
        return
    adata.obs[obs_key] = adata.obs[obs_key].astype("category")
    fig = sc.pl.matrixplot(adata, var_names=top_genes, groupby=obs_key,
                           layer="lognorm", show=False, return_fig=True,
                           title=f"Top {n} markers per cluster")
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close()


def plot_marker_score_umaps(adata, markers: dict, out: Path) -> None:
    """One UMAP panel per cell type colored by its marker score."""
    score_keys = []
    for ct in markers:
        key = ("score_" + ct
               .replace(" ", "_").replace("/", "_")
               .replace("(", "").replace(")", ""))
        if key in adata.obs.columns:
            score_keys.append((ct, key))
    if not score_keys:
        return
    n = len(score_keys)
    ncols = min(4, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows))
    axes = np.array(axes).flatten()
    for ax, (ct, key) in zip(axes, score_keys):
        sc.pl.umap(adata, color=key, ax=ax, show=False, frameon=False,
                   color_map="Reds", size=6, title=ct)
    for ax in axes[len(score_keys):]:
        ax.set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)



def assign_provisional_celltype(adata, markers: dict) -> str:
    """Assign a provisional cell type label to each cell by taking whichever
    curated marker score is highest. Stored in adata.obs["provisional_celltype"].

    This gives real names (Microglia, Excitatory neurons, ...) for composition
    plots even when CellTypist is not configured. Labels are explicitly marked
    PRELIMINARY — override via annotation_summary.csv manual_annotation column.

    Returns the obs column name ("provisional_celltype").
    """
    score_cols = []
    ct_names = []
    for ct in markers:
        key = ("score_" + ct
               .replace(" ", "_").replace("/", "_")
               .replace("(", "").replace(")", ""))
        if key in adata.obs.columns:
            score_cols.append(key)
            ct_names.append(ct)

    if not score_cols:
        print("  [warn] No marker scores found — provisional labels unavailable.")
        return "leiden"

    scores = adata.obs[score_cols].values          # n_cells × n_types
    best_idx = scores.argmax(axis=1)               # index of highest score per cell
    adata.obs["provisional_celltype"] = [ct_names[i] for i in best_idx]
    n_types = adata.obs["provisional_celltype"].nunique()
    print(f"  Provisional cell type labels assigned: {n_types} types")
    print(f"  Distribution:")
    dist = adata.obs["provisional_celltype"].value_counts()
    for ct, n in dist.items():
        print(f"    {ct}: {n:,} cells ({100*n/len(adata):.1f}%)")
    print(f"  Edit annotation_summary.csv to correct any misassignments.")
    return "provisional_celltype"

# ---------------------------------------------------------------------------
# Plots: composition
# ---------------------------------------------------------------------------

def _stacked_bar(ct_table: pd.DataFrame, title: str, xlabel: str,
                 out: Path, figwidth_per_bar: float = 0.4,
                 min_width: float = 7.0) -> None:
    """Generic stacked bar helper. ct_table rows=groups, cols=categories (fractions)."""
    n = len(ct_table)
    width = max(min_width, figwidth_per_bar * n)
    fontsize = 9 if n <= 15 else (8 if n <= 30 else 6)
    cmap = plt.get_cmap("tab20")
    colors = [cmap(i % 20) for i in range(len(ct_table.columns))]

    fig, ax = plt.subplots(figsize=(width, 5))
    ct_table.plot(kind="bar", stacked=True, ax=ax, width=0.8,
                  color=colors, edgecolor="none", legend=True)
    ax.set_ylabel("fraction of cells")
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=6, ncol=2)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=fontsize)
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_cluster_composition_by_sample(adata, obs_key: str, out: Path) -> None:
    """Stacked bar: fraction of each Leiden cluster's cells from each sample.
    Single-sample clusters (one bar ~100% one color) = likely batch artifact."""
    ct = pd.crosstab(adata.obs[obs_key], adata.obs["sample_id"], normalize="index")
    _stacked_bar(
        ct,
        title=f"Sample composition per cluster — single-sample clusters = potential batch artifact",
        xlabel=obs_key,
        out=out,
    )


def plot_celltype_composition_by_sample(adata, celltype_key: str, out: Path) -> None:
    """Stacked bar: fraction of each cell type's cells from each sample.
    Tells you whether any cell type is dominated by one sample."""
    ct = pd.crosstab(adata.obs[celltype_key], adata.obs["sample_id"], normalize="index")
    _stacked_bar(
        ct,
        title=(f"Cell type composition per sample [PRELIMINARY: {celltype_key} labels]\n"
               f"Check annotation_summary.csv — correct manual_annotation if labels are wrong."),
        xlabel="cell type",
        out=out,
    )


def plot_celltype_composition_by_group(adata, celltype_key: str, out: Path) -> None:
    """Stacked bar: fraction of each group's cells in each cell type.
    One bar per sample, colored by cell type — first look at whether
    Early/Late/Relaxed differ in composition. This is diagnostic; the
    quantitative test is in Phase 8a (scCODA + propeller)."""
    if "group" not in adata.obs.columns:
        return

    # Fraction of each cell type within each sample
    ct = pd.crosstab(adata.obs["sample_id"], adata.obs[celltype_key], normalize="index")

    # Sort samples by group then sample_id for a clean layout
    sample_meta = (adata.obs[["sample_id", "group", "age"]]
                   .drop_duplicates()
                   .set_index("sample_id")
                   .reindex(ct.index))
    ct = ct.loc[sample_meta.sort_values(["group", "age"]).index]

    n = len(ct)
    width = max(8, 0.4 * n)
    fontsize = 9 if n <= 15 else (8 if n <= 30 else 6)
    cmap = plt.get_cmap("tab20")
    colors = [cmap(i % 20) for i in range(len(ct.columns))]

    fig, ax = plt.subplots(figsize=(width, 5))
    ct.plot(kind="bar", stacked=True, ax=ax, width=0.8,
            color=colors, edgecolor="none", legend=True)
    ax.set_ylabel("fraction of cells")
    ax.set_xlabel("sample (sorted by group × age)")
    ax.set_title(
        f"Cell type composition per sample [PRELIMINARY: {celltype_key} labels]\n"
        f"Check annotation_summary.csv if labels look wrong. Quantitative test in Phase 8a (scCODA)."
    )
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=6, ncol=2)

    # Add group-level x-axis separators
    group_order = sample_meta.sort_values(["group", "age"])["group"]
    boundaries = np.where(np.diff(group_order.values != group_order.values))[0] + 0.5
    for b in boundaries:
        ax.axvline(b, color="black", lw=1.5, ls="--")

    # Label groups below bars
    groups = group_order.unique()
    ticks = {g: [] for g in groups}
    for i, (sid, row) in enumerate(group_order.items()):
        ticks[row].append(i)
    for g, idxs in ticks.items():
        mid = np.mean(idxs)
        ax.text(mid, -0.18, g, ha="center", va="top",
                transform=ax.get_xaxis_transform(), fontsize=8, fontweight="bold")

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=fontsize)
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Annotation summary table
# ---------------------------------------------------------------------------

def build_annotation_summary(adata, top_markers: pd.DataFrame,
                              obs_key: str = "leiden") -> pd.DataFrame:
    clusters = sorted(adata.obs[obs_key].unique(), key=int)
    rows = []
    for c in clusters:
        mask = adata.obs[obs_key] == c
        row = {"cluster": c, "n_cells": int(mask.sum())}
        if "celltypist_majority" in adata.obs.columns:
            row["celltypist_majority"] = adata.obs.loc[mask, "celltypist_majority"].mode()[0]
            if "celltypist_conf_score" in adata.obs.columns:
                row["celltypist_conf_median"] = round(
                    float(adata.obs.loc[mask, "celltypist_conf_score"].median()), 3)
        top3 = top_markers[top_markers["group"] == c].head(3)["names"].tolist()
        row["top_markers"] = ", ".join(top3)
        row["manual_annotation"] = ""   # fill in notebook
        rows.append(row)
    return pd.DataFrame(rows)


def save_composition_table(adata, celltype_key: str, out: Path) -> None:
    """Save cell type × sample counts + fractions to CSV."""
    counts = pd.crosstab(adata.obs["sample_id"], adata.obs[celltype_key])
    # crosstab may produce a CategoricalIndex on columns — flatten to strings
    counts.columns = counts.columns.astype(str)
    fracs = counts.div(counts.sum(axis=1), axis=0)
    meta = (adata.obs[["sample_id", "group", "age", "sex", "pool"]]
            .drop_duplicates().set_index("sample_id"))
    out_df = meta.join(fracs)
    out_df.to_csv(out)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Phase 7: cell type annotation")
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()

    print(f"\n=== Phase 7: Annotation ===")
    print(f"Config: {args.config}")

    cfg = load_config(args.config)
    tissue = cfg["tissue"]

    in_path = Path(cfg["results_dir"]) / "h5ad" / "07_clustered" / "all_samples.h5ad"
    if not in_path.is_file():
        sys.exit(f"ERROR: missing {in_path}. Run 06_clustering.py first.")

    out_dir = Path(cfg["results_dir"]) / "h5ad" / "08_annotated"
    plot_dir = Path(cfg["results_dir"]) / "plots" / "07_annotation"
    table_dir = Path(cfg["results_dir"]) / "tables"
    for d in (out_dir, plot_dir, table_dir):
        d.mkdir(parents=True, exist_ok=True)

    print(f"\n[1/4] Loading {in_path}...")
    adata = sc.read_h5ad(in_path)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    if "leiden" not in adata.obs.columns:
        sys.exit("ERROR: 'leiden' not in adata.obs. Run 06_clustering.py first.")
    print(f"  Leiden clusters: {adata.obs['leiden'].nunique()}")

    # Recompute lognorm — dropped after Phase 5. All marker steps use this.
    print(f"  Recomputing lognorm layer...")
    add_lognorm(adata)

    # --- Track 1: CellTypist ---
    print(f"\n[2/4] Track 1: CellTypist reference annotation...")
    model_source = cfg.get("annotation", {}).get("celltypist_model")
    if not model_source:
        print(f"  No celltypist_model in config — skipping.")
        print(f"  To enable, add to YAML:")
        print(f"    annotation:")
        if tissue == "brain":
            print(f"      celltypist_model: Developing_Mouse_Brain")
        else:
            print(f"      celltypist_model: /path/to/placenta_model.pkl")
    else:
        ct_result = run_celltypist(adata, model_source)
        if ct_result is not None:
            adata.obs["celltypist_majority"] = ct_result["majority_voting"].values
            adata.obs["celltypist_predicted"] = ct_result["predicted_labels"].values
            adata.obs["celltypist_conf_score"] = ct_result["conf_score"].values
            ct_result.index = adata.obs_names
            ct_result.to_csv(table_dir / "celltypist_predictions.csv")
            print(f"  {adata.obs['celltypist_majority'].nunique()} cell types predicted")

    # --- Track 2: Marker-based ---
    print(f"\n[3/4] Track 2: marker genes + scoring (on lognorm)...")
    markers = get_markers(cfg)
    print(f"  Marker sets: {list(markers.keys())}")
    score_marker_sets(adata, markers)
    top_markers = run_marker_genes(adata, obs_key="leiden")
    top_markers.to_csv(table_dir / "marker_genes_per_cluster.csv", index=False)
    print(f"  Top markers written: {len(top_markers)} rows")

    summary = build_annotation_summary(adata, top_markers)
    summary.to_csv(table_dir / "annotation_summary.csv", index=False)
    print(f"\n  Annotation summary:")
    print(summary.to_string(index=False))

    # Decide which label column to use for composition plots.
    # Priority: CellTypist majority > provisional (highest marker score) > leiden
    if "celltypist_majority" in adata.obs.columns:
        celltype_key = "celltypist_majority"
        print(f"\n  Composition plots will use CellTypist majority labels.")
    else:
        print(f"\n  No CellTypist model configured — assigning provisional labels from marker scores...")
        celltype_key = assign_provisional_celltype(adata, markers)

    # --- Plots ---
    print(f"\n[4/4] Generating plots...")

    # UMAP + markers
    plot_leiden_for_reference(adata, "leiden",
                               plot_dir / "umap_leiden_for_annotation.png")
    plot_marker_dotplot(adata, markers, "leiden",
                        plot_dir / "marker_dotplot.png")
    plot_top_marker_heatmap(adata, top_markers, "leiden",
                             plot_dir / "marker_heatmap_top10.png")
    plot_marker_score_umaps(adata, markers,
                             plot_dir / "umap_marker_scores.png")
    plot_umap_celltypist(adata,
                         plot_dir / "umap_celltypist.png",
                         plot_dir / "umap_celltypist_confidence.png")

    # Composition diagnostics
    plot_cluster_composition_by_sample(
        adata, "leiden",
        plot_dir / "cluster_composition_by_sample.png")
    plot_celltype_composition_by_sample(
        adata, celltype_key,
        plot_dir / "celltype_composition_by_sample.png")
    plot_celltype_composition_by_group(
        adata, celltype_key,
        plot_dir / "celltype_composition_by_group.png")

    # Composition CSV
    save_composition_table(adata, celltype_key,
                           table_dir / "celltype_composition.csv")

    # Drop lognorm before saving (same policy as Phase 5)
    if "lognorm" in adata.layers:
        del adata.layers["lognorm"]
    adata.obs["manual_annotation"] = ""
    adata.write_h5ad(out_dir / "all_samples.h5ad")

    print(f"\n  Written: {out_dir / 'all_samples.h5ad'}")
    print(f"  Plots:   {plot_dir}")
    print(f"\n✓ Phase 7 complete.")
    print(f"\nKey plots to review:")
    print(f"  marker_dotplot.png                    — which cluster is which cell type?")
    print(f"  umap_marker_scores.png                — per-cell type score on UMAP")
    print(f"  cluster_composition_by_sample.png     — any single-sample clusters (batch artifacts)?")
    print(f"  celltype_composition_by_sample.png    — PRELIMINARY: is any cell type one-sample-dominated?")
    print(f"  celltype_composition_by_group.png     — PRELIMINARY: first look at stress vs composition")
    print(f"")
    print(f"  Composition plots use {celltype_key!r} labels.")
    print(f"  If labels look wrong, edit the manual_annotation column in:")
    print(f"    {table_dir / 'annotation_summary.csv'}")
    print(f"  Then re-run this script after transferring labels to adata.obs.")
    print(f"\nNext steps:")
    print(f"  1. Review annotation_summary.csv — correct manual_annotation where needed")
    print(f"  2. Transfer labels to adata.obs in a notebook")
    print(f"  3. Run 07b_subcluster.py for microglia + oligodendrocyte lineage")
    print(f"\nNext automated step: 07b_subcluster.py\n")


if __name__ == "__main__":
    main()
