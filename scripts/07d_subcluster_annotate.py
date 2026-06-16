#!/usr/bin/env python
"""
07d_subcluster_annotate.py — Phase 7d: annotate integer subcluster IDs with
cell type names, using CellTypist (where a model exists) and/or marker scoring
from a literature-curated YAML, then write the name back into the h5ad.

Why this exists: 07b_subcluster.py produces `obs['subcluster']` as Leiden
integers (0, 1, 2, …). Downstream scripts (08b DE, 08c pathway) label their
outputs by that column, so figures show "subcluster=2" instead of "DAM" unless
we map those integers to names here. This is the mapping step.

Two annotation tracks (same philosophy as 07_annotation.py):

  Track A — CellTypist reference transfer
    Runs when `annotation.celltypist_models` has an entry for this age AND
    celltypist is importable. Uses majority vote per integer subcluster.
    Produces `obs['celltypist_subcluster']`.

  Track B — Marker scoring (always runs as fallback / cross-check)
    Reads `subcluster_markers:` from the YAML (see config/subcluster_markers.yaml).
    Scores each gene set with `sc.tl.score_genes` on lognorm, aggregates per
    integer subcluster (mean score), assigns the top-scoring name.
    Produces `obs['marker_subcluster']`.

Final name:
  `obs['subcluster_name']` = CellTypist if Track A ran and confidence >= threshold;
  else marker-based top hit. Always filled — never NaN.

YAML marker format (see config/subcluster_markers.yaml for full example):
  subcluster_markers:
    microglia:                      # must match --celltype slug (lowercased)
      Homeostatic_Microglia:
        markers: [P2ry12, Tmem119, Cx3cr1, Hexb, Tgfbr1]
        refs: "Keren-Shaul 2017; Bennett 2016"
      DAM:
        markers: [Trem2, Apoe, Lpl, Cst7, Tyrobp, Axl, Cd9]
        refs: "Keren-Shaul 2017; Deczkowska 2018"
      ...

Usage:
  uv run python scripts/07d_subcluster_annotate.py \\
      --config config/brain.yaml \\
      --celltype Microglia \\
      --markers config/subcluster_markers.yaml

  # with explicit age (for CellTypist model selection):
  uv run python scripts/07d_subcluster_annotate.py \\
      --config config/brain.yaml \\
      --celltype Oligodendrocytes \\
      --markers config/subcluster_markers.yaml \\
      --age P1

  # force marker-only (skip CellTypist):
  uv run python scripts/07d_subcluster_annotate.py \\
      --config config/placenta.yaml \\
      --celltype trophoblast \\
      --markers config/subcluster_markers.yaml \\
      --no-celltypist

Inputs:
  {results_dir}/h5ad/08c_subclustered/{slug}.h5ad  (from 07b_subcluster.py)

Outputs (in place — overwrites the subcluster h5ad):
  obs['subcluster_name']       : final name (CellTypist or marker-based)
  obs['celltypist_subcluster'] : CellTypist majority label per cluster (if run)
  obs['marker_subcluster']     : top marker-score label per cluster (always)
  obs['subcluster_confidence'] : 'celltypist' | 'marker' | 'unresolved'

  plots/{slug}/
    subcluster_names_umap.png       : UMAP coloured by subcluster_name
    subcluster_marker_scores.png    : heatmap of marker scores per integer cluster
    subcluster_celltypist_umap.png  : CellTypist labels (if run)
  tables/
    subcluster_{slug}_annotation.csv : integer → name mapping + scores + source
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

from _utils import load_config, add_lognorm, phase_table_dir


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def slugify(name: str) -> str:
    return re.sub(r"[^0-9a-zA-Z]+", "_", name).strip("_").lower()


def load_marker_cfg(markers_yaml: Path, slug: str) -> dict:
    """Load the subcluster_markers block for this cell type slug.

    Returns {} if the file has no entry for slug (marker track will be skipped
    with a clear warning rather than crashing).
    """
    with markers_yaml.open() as f:
        cfg = yaml.safe_load(f)
    block = cfg.get("subcluster_markers", {})
    # Try exact slug match, then case-insensitive
    entry = block.get(slug)
    if entry is None:
        for k, v in block.items():
            if k.lower() == slug.lower():
                entry = v
                break
    if entry is None:
        print(f"  [warn] No subcluster_markers entry for '{slug}' in {markers_yaml}.")
        print(f"         Available: {sorted(block.keys())}")
        print(f"         Marker track will be skipped. Add an entry to enable it.")
    return entry or {}


# ---------------------------------------------------------------------------
# Track A: CellTypist
# ---------------------------------------------------------------------------

def run_celltypist(adata, model_name: str, subcluster_key: str = "subcluster",
                   confidence_threshold: float = 0.5) -> pd.Series:
    """Run CellTypist majority-vote per integer subcluster.

    Returns a Series indexed by obs_names with the majority label per cluster.
    Fills 'unresolved' where confidence < threshold.
    """
    try:
        import celltypist
        from celltypist import models
    except ImportError:
        print("  [warn] celltypist not importable. Skipping Track A.")
        return None

    print(f"  CellTypist: loading model '{model_name}'...")
    try:
        model = models.Model.load(model=model_name)
    except Exception as e:
        print(f"  [warn] Could not load CellTypist model '{model_name}': {e}")
        print(f"         Skipping Track A.")
        return None

    # CellTypist needs lognorm in X for prediction
    tmp = adata.copy()
    if "lognorm" in tmp.layers:
        tmp.X = tmp.layers["lognorm"].copy()
    else:
        add_lognorm(tmp)
        tmp.X = tmp.layers["lognorm"].copy()

    print(f"  CellTypist: predicting on {tmp.n_obs:,} nuclei...")
    pred = celltypist.annotate(tmp, model=model, majority_voting=True,
                               over_clustering=subcluster_key)
    # majority_voting produces pred.predicted_labels['majority_voting'] per cell
    labels = pred.predicted_labels["majority_voting"]
    conf = pred.probability_matrix.max(axis=1)

    # Map low-confidence cells to 'unresolved'
    labels = labels.copy()
    labels[conf < confidence_threshold] = "unresolved"
    return labels


# ---------------------------------------------------------------------------
# Track B: Marker scoring
# ---------------------------------------------------------------------------

def run_marker_scoring(adata, marker_cfg: dict,
                       subcluster_key: str = "subcluster") -> tuple[pd.Series, pd.DataFrame]:
    """Score each marker gene set per cell, then aggregate per integer cluster.

    Returns:
      - per_cell_labels : Series indexed by obs_names, top-scoring name per cell
      - score_table     : DataFrame (clusters × cell_types) of mean scores
    """
    if not marker_cfg:
        return None, None

    # Need lognorm
    if "lognorm" not in adata.layers:
        add_lognorm(adata)

    scores = {}
    for name, spec in marker_cfg.items():
        raw_genes = spec.get("markers", [])
        present = [g for g in raw_genes if g in adata.var_names]
        missing = [g for g in raw_genes if g not in adata.var_names]
        if missing:
            print(f"    [{name}] {len(missing)} marker(s) not in var_names (skipped): "
                  f"{missing[:5]}{'…' if len(missing) > 5 else ''}")
        if len(present) < 2:
            print(f"    [{name}] fewer than 2 markers found — skipping this set.")
            continue
        # score_genes modifies adata.obs in place; use a key per name
        key = f"_mscore_{slugify(name)}"
        # score on lognorm LAYER — adata.X is raw counts in 08c objects, and
        # use_raw=False alone would score on raw (library-size-dominated, wrong).
        sc.tl.score_genes(adata, gene_list=present, score_name=key,
                          use_raw=False, layer="lognorm")
        scores[name] = adata.obs[key].values
        # Clean up obs
        del adata.obs[key]

    if not scores:
        print("  [warn] No marker sets scored (all had <2 genes in var_names).")
        return None, None

    score_df = pd.DataFrame(scores, index=adata.obs_names)

    # Per-cell: top-scoring name (argmax across sets)
    per_cell = score_df.idxmax(axis=1)

    # Per-cluster mean scores → for the annotation table and heatmap
    clusters = adata.obs[subcluster_key].astype(str)
    cluster_scores = score_df.copy()
    cluster_scores["_cluster"] = clusters.values
    cluster_mean = cluster_scores.groupby("_cluster")[list(scores.keys())].mean()

    return per_cell, cluster_mean


# ---------------------------------------------------------------------------
# Merge tracks → subcluster_name
# ---------------------------------------------------------------------------

def merge_tracks(adata, ct_labels, marker_labels, marker_cluster_scores,
                 subcluster_key: str = "subcluster",
                 confidence_threshold: float = 0.5):
    """Build obs['subcluster_name'] from the two tracks.

    Priority: CellTypist per-cluster majority label (if it ran and isn't
    'unresolved') > marker top-hit > 'unresolved'.
    """
    clusters = adata.obs[subcluster_key].astype(str)
    unique_clusters = sorted(clusters.unique(), key=lambda x: int(x) if x.isdigit() else x)

    # Per-cluster decisions
    cluster_name = {}
    cluster_confidence = {}

    for cl in unique_clusters:
        mask = clusters == cl

        # CellTypist majority vote for this cluster
        ct_name = None
        if ct_labels is not None:
            ct_votes = ct_labels[mask]
            mode = ct_votes.value_counts()
            top = mode.index[0] if len(mode) else "unresolved"
            frac = mode.iloc[0] / mask.sum() if mask.sum() > 0 else 0
            if top != "unresolved" and frac >= confidence_threshold:
                ct_name = top

        # Marker top-hit for this cluster
        mk_name = None
        if marker_cluster_scores is not None and cl in marker_cluster_scores.index:
            row = marker_cluster_scores.loc[cl]
            mk_name = row.idxmax() if row.max() > 0 else None

        # Decision
        if ct_name:
            cluster_name[cl] = ct_name
            cluster_confidence[cl] = "celltypist"
        elif mk_name:
            cluster_name[cl] = mk_name
            cluster_confidence[cl] = "marker"
        else:
            cluster_name[cl] = "unresolved"
            cluster_confidence[cl] = "unresolved"

    # Map back to per-cell
    adata.obs["subcluster_name"] = clusters.map(cluster_name).fillna("unresolved")
    adata.obs["subcluster_confidence"] = clusters.map(cluster_confidence).fillna("unresolved")

    if ct_labels is not None:
        adata.obs["celltypist_subcluster"] = ct_labels.reindex(adata.obs_names).fillna("unresolved")
    if marker_labels is not None:
        adata.obs["marker_subcluster"] = marker_labels.reindex(adata.obs_names).fillna("unresolved")

    return cluster_name, cluster_confidence


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_name_umap(adata, plot_dir: Path, slug: str):
    if "X_umap" not in adata.obsm:
        print("  [skip] no X_umap in obsm — run 07b first to get UMAP coordinates.")
        return
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(adata, color="subcluster_name", ax=ax, show=False,
               frameon=False, legend_loc="on data", legend_fontsize=7,
               title=f"{slug}: subcluster names")
    fig.tight_layout()
    fig.savefig(plot_dir / "subcluster_names_umap.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_celltypist_umap(adata, plot_dir: Path, slug: str):
    if "celltypist_subcluster" not in adata.obs or "X_umap" not in adata.obsm:
        return
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(adata, color="celltypist_subcluster", ax=ax, show=False,
               frameon=False, legend_loc="on data", legend_fontsize=7,
               title=f"{slug}: CellTypist labels")
    fig.tight_layout()
    fig.savefig(plot_dir / "subcluster_celltypist_umap.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_marker_scores(cluster_scores: pd.DataFrame, plot_dir: Path, slug: str):
    if cluster_scores is None or cluster_scores.empty:
        return
    fig, ax = plt.subplots(figsize=(max(6, 0.7 * len(cluster_scores.columns)),
                                    max(3, 0.5 * len(cluster_scores))))
    import seaborn as sns
    sns.heatmap(cluster_scores, annot=True, fmt=".2f", cmap="YlOrRd",
                ax=ax, linewidths=0.5)
    ax.set_title(f"{slug}: mean marker scores per subcluster")
    ax.set_xlabel("cell type")
    ax.set_ylabel("subcluster")
    fig.tight_layout()
    fig.savefig(plot_dir / "subcluster_marker_scores.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Phase 7d: annotate subcluster integers")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--celltype", required=True,
                        help="Cell type slug from 07b (e.g. 'Microglia', 'trophoblast')")
    parser.add_argument("--markers", required=True, type=Path,
                        help="Path to subcluster_markers.yaml")
    parser.add_argument("--age", default=None,
                        help="Age for CellTypist model selection (e.g. P1, 4W, E12.5). "
                             "If omitted, uses the most common age in obs.")
    parser.add_argument("--celltypist-confidence", type=float, default=0.5,
                        help="Min CellTypist majority fraction to trust the label (default 0.5)")
    parser.add_argument("--no-celltypist", action="store_true",
                        help="Skip CellTypist entirely (marker track only)")
    args = parser.parse_args()

    print(f"\n=== Phase 7d: Subcluster annotation for '{args.celltype}' ===")
    cfg = load_config(args.config)
    tissue = cfg["tissue"]
    slug = slugify(args.celltype)

    base = Path(cfg["results_dir"]) / "h5ad" / "08c_subclustered"
    in_path = base / f"{slug}.h5ad"
    if not in_path.is_file():
        sys.exit(f"ERROR: subcluster h5ad not found: {in_path}\n"
                 f"  Run 07b_subcluster.py --celltype '{args.celltype}' first.")

    plot_dir = Path(cfg["results_dir"]) / "plots" / "07b_subcluster" / slug
    table_dir = phase_table_dir(cfg, "07d_subcluster_annotate")
    plot_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    if not args.markers.is_file():
        sys.exit(f"ERROR: --markers file not found: {args.markers}\n"
                 f"  Copy config/subcluster_markers.yaml from the repo and fill it in.")

    print(f"  Input:   {in_path}")
    print(f"  Markers: {args.markers}")
    print(f"  Tissue:  {tissue}")

    print(f"\n[1/4] Loading subcluster object...")
    adata = sc.read_h5ad(in_path)
    print(f"  {adata.n_obs:,} cells, {adata.obs['subcluster'].nunique()} subclusters")

    if "subcluster" not in adata.obs.columns:
        sys.exit("ERROR: 'subcluster' not in obs. Was 07b_subcluster.py run on this object?")

    # Determine age for model selection
    age = args.age
    if age is None and "age" in adata.obs.columns:
        age = adata.obs["age"].value_counts().index[0]
        print(f"  Age (inferred from obs): {age}")
    elif age:
        print(f"  Age (from --age): {age}")

    # -----------------------------------------------------------------------
    # Track A: CellTypist
    # -----------------------------------------------------------------------
    ct_labels = None
    if not args.no_celltypist:
        print(f"\n[2/4] Track A: CellTypist...")
        models_cfg = cfg.get("annotation", {}).get("celltypist_models", {})
        model_name = models_cfg.get(age) if age else None
        if model_name:
            print(f"  Model for age {age}: {model_name}")
            ct_labels = run_celltypist(adata, model_name,
                                       confidence_threshold=args.celltypist_confidence)
        else:
            ages_with_models = list(models_cfg.keys())
            print(f"  No CellTypist model for age='{age}'. "
                  f"Ages with models: {ages_with_models or 'none'}.")
            print(f"  Skipping Track A. Add a model path to annotation.celltypist_models "
                  f"in the YAML to enable it.")
    else:
        print(f"\n[2/4] Track A: skipped (--no-celltypist).")

    # -----------------------------------------------------------------------
    # Track B: Marker scoring
    # -----------------------------------------------------------------------
    print(f"\n[3/4] Track B: Marker scoring...")
    marker_cfg = load_marker_cfg(args.markers, slug)
    if marker_cfg:
        print(f"  {len(marker_cfg)} cell type(s) defined for '{slug}':")
        for name, spec in marker_cfg.items():
            n = len(spec.get("markers", []))
            refs = spec.get("refs", "")
            print(f"    {name}: {n} markers  [{refs}]")
        add_lognorm(adata)
        marker_labels, cluster_scores = run_marker_scoring(adata, marker_cfg)
    else:
        marker_labels, cluster_scores = None, None

    if ct_labels is None and marker_labels is None:
        sys.exit(
            "ERROR: both CellTypist and marker scoring produced no labels.\n"
            "  Check that:\n"
            "    (a) a CellTypist model is set for this age in the YAML, OR\n"
            "    (b) subcluster_markers.yaml has an entry for this cell type slug.\n"
            f"  slug = '{slug}', age = '{age}'"
        )

    # -----------------------------------------------------------------------
    # Merge → subcluster_name
    # -----------------------------------------------------------------------
    print(f"\n[4/4] Merging tracks → subcluster_name...")
    cluster_name, cluster_confidence = merge_tracks(
        adata, ct_labels, marker_labels, cluster_scores,
        confidence_threshold=args.celltypist_confidence,
    )

    print(f"\n  Integer → name mapping:")
    for cl in sorted(cluster_name.keys(), key=lambda x: int(x) if x.isdigit() else x):
        print(f"    {cl:>3} → {cluster_name[cl]}  [{cluster_confidence[cl]}]")

    n_unresolved = sum(1 for v in cluster_name.values() if v == "unresolved")
    if n_unresolved:
        print(f"\n  [warn] {n_unresolved} subcluster(s) are 'unresolved' — no confident "
              f"label from either track. Consider adding more markers to the YAML or "
              f"inspecting the marker dotplot from 07b.")

    # -----------------------------------------------------------------------
    # Plots
    # -----------------------------------------------------------------------
    plot_name_umap(adata, plot_dir, slug)
    plot_celltypist_umap(adata, plot_dir, slug)
    plot_marker_scores(cluster_scores, plot_dir, slug)

    # -----------------------------------------------------------------------
    # Annotation table
    # -----------------------------------------------------------------------
    rows = []
    for cl in sorted(cluster_name.keys(), key=lambda x: int(x) if x.isdigit() else x):
        n_cells = int((adata.obs["subcluster"].astype(str) == cl).sum())
        row = {
            "subcluster_int": cl,
            "subcluster_name": cluster_name[cl],
            "confidence_source": cluster_confidence[cl],
            "n_cells": n_cells,
        }
        # Add per-cluster marker scores if available
        if cluster_scores is not None and cl in cluster_scores.index:
            for ct, score in cluster_scores.loc[cl].items():
                row[f"score_{slugify(ct)}"] = round(score, 4)
        rows.append(row)

    ann_table = pd.DataFrame(rows)
    ann_path = table_dir / f"07d_subcluster_{slug}_annotation.csv"
    ann_table.to_csv(ann_path, index=False)
    print(f"\n  Annotation table: {ann_path}")

    # -----------------------------------------------------------------------
    # Write back to h5ad (in place)
    # -----------------------------------------------------------------------
    # Drop lognorm layer before saving (keep disk lean — recomputed on demand)
    if "lognorm" in adata.layers:
        del adata.layers["lognorm"]

    # Cast new obs columns to category
    for col in ("subcluster_name", "subcluster_confidence",
                "celltypist_subcluster", "marker_subcluster"):
        if col in adata.obs.columns:
            adata.obs[col] = adata.obs[col].astype("category")

    adata.write_h5ad(in_path)
    print(f"  Updated h5ad (in place): {in_path}")
    print(f"  New obs columns: subcluster_name, subcluster_confidence"
          + (", celltypist_subcluster" if ct_labels is not None else "")
          + (", marker_subcluster" if marker_labels is not None else ""))

    print(f"\n  Plots: {plot_dir}")
    print(f"\n✓ Phase 7d complete for '{args.celltype}'.")
    print(f"  Downstream 08b/08c --subcluster runs will now show named subtypes.")
    print(f"  Inspect 'subcluster_names_umap.png' and 'subcluster_marker_scores.png'.")
    print(f"  Edit subcluster_markers.yaml and re-run to refine unresolved clusters.\n")


if __name__ == "__main__":
    main()
