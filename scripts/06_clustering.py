#!/usr/bin/env python
"""
06_clustering.py — Phase 6: multi-resolution Leiden clustering.

Builds a neighbor graph on the scVI latent space and runs Leiden at
resolutions [0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0]. Auto-selects resolution
via the geometric knee on the n_clusters vs resolution curve, then produces
diagnostic plots so you can confirm or override.

Resolution selection (in order of priority):
  1. --resolution CLI flag — use this after inspecting resolution_selection.png
  2. clustering.resolution in YAML
  3. Auto: geometric knee on the n_clusters vs resolution curve

Geometric knee: find the point of maximum perpendicular distance from the line
connecting the first and last points of the curve. Same idea as the PCA elbow —
where does adding more resolution stop meaningfully splitting clusters?

Workflow:
  - First run: let auto-selection pick a candidate.
  - Open resolution_selection.png — inspect the knee curve and confirm the pick.
  - If you want a different value: re-run with --resolution 0.6 (or whatever).

Usage:
  uv run python scripts/06_clustering.py --config config/dev.yaml
  uv run python scripts/06_clustering.py --config config/brain.yaml
  uv run python scripts/06_clustering.py --config config/brain.yaml --resolution 0.6

Inputs:
  {results_dir}/h5ad/06_integrated/all_samples.h5ad  (from Phase 5)

Outputs:
  {results_dir}/h5ad/07_clustered/all_samples.h5ad
  {results_dir}/plots/06_clustering/
    - resolution_selection.png         : knee curve — inspect this first
    - clustree.png                     : cluster splitting diagram
    - umap_leiden_res{X}.png           : UMAP per resolution
    - umap_leiden_chosen.png           : UMAP at chosen resolution
    - cluster_qc_metrics.png           : per-cluster QC (catches junk clusters)
    - cluster_composition_by_sample.png: per-cluster sample composition
  {results_dir}/tables/
    - summary_clustering.csv           : n_clusters + knee_distance per resolution
    - cluster_qc_per_resolution.csv    : per-cluster cell counts + QC medians
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

from _utils import load_config, phase_paths, phase_table_dir


RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0]


# ---------------------------------------------------------------------------
# Clustering + resolution selection
# ---------------------------------------------------------------------------

def run_leiden(adata, resolutions: list[float], seed: int) -> pd.DataFrame:
    """Run Leiden at each resolution, store results in adata.obs.
    Returns DataFrame with [resolution, n_clusters, obs_key]."""
    rows = []
    for res in resolutions:
        key = f"leiden_r{res:.1f}".replace(".", "_")
        sc.tl.leiden(adata, resolution=res, random_state=seed, key_added=key, flavor="igraph", n_iterations=2, directed=False)
        n = adata.obs[key].nunique()
        rows.append({"resolution": res, "n_clusters": n, "obs_key": key})
        print(f"  res={res:.1f} → {n} clusters")
    return pd.DataFrame(rows)


def knee_distances(resolutions: np.ndarray, n_clusters: np.ndarray) -> np.ndarray:
    """Perpendicular distance of each point from the line connecting endpoints.

    Normalize both axes to [0,1] first so axis scales don't bias the result.
    The point with maximum distance is the knee.
    """
    x = np.array(resolutions, dtype=float)
    y = np.array(n_clusters, dtype=float)
    x_n = (x - x.min()) / (x.max() - x.min() + 1e-9)
    y_n = (y - y.min()) / (y.max() - y.min() + 1e-9)
    dx, dy = x_n[-1] - x_n[0], y_n[-1] - y_n[0]
    line_len = np.sqrt(dx**2 + dy**2) + 1e-9
    dists = np.abs(dy * x_n - dx * y_n + x_n[-1] * y_n[0] - y_n[-1] * x_n[0]) / line_len
    return dists


def choose_resolution(resolution_df: pd.DataFrame,
                      forced: float | None = None) -> tuple[float, str]:
    """Return (resolution, method_label)."""
    if forced is not None:
        return forced, "manual override"
    res = resolution_df["resolution"].values
    nc = resolution_df["n_clusters"].values
    if len(res) >= 3:
        dists = knee_distances(res, nc)
        return float(res[np.argmax(dists)]), "geometric knee"
    return 0.6, "fallback default"


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_resolution_selection(resolution_df: pd.DataFrame, chosen: float,
                               method: str, out: Path) -> None:
    """Knee curve with chosen resolution marked and a distance inset.

    The main panel shows n_clusters vs resolution — look for where the curve
    flattens. The inset shows the perpendicular distances the algorithm uses,
    so you can see why it picked what it picked.
    """
    res = resolution_df["resolution"].values
    nc = resolution_df["n_clusters"].values
    dists = resolution_df["knee_distance"].values

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(res, nc, marker="o", color="steelblue", lw=1.5, label="n clusters")
    ax.axvline(chosen, color="red", ls="--", lw=1.5,
               label=f"chosen: {chosen}  ({method})")

    chosen_idx = np.argmin(np.abs(res - chosen))
    ax.scatter([res[chosen_idx]], [nc[chosen_idx]], color="red", s=100, zorder=5)

    for x, y in zip(res, nc):
        ax.annotate(str(y), (x, y), xytext=(0, 7), textcoords="offset points",
                    ha="center", fontsize=8, color="dimgray")

    ax.set_xlabel("Leiden resolution")
    ax.set_ylabel("n clusters")
    ax.set_title(
        f"Resolution selection  (auto: {method})\n"
        f"To override: re-run with --resolution <value>"
    )
    ax.legend(fontsize=9)

    # Inset: the raw knee distances
    ax_in = ax.inset_axes([0.58, 0.08, 0.38, 0.32])
    ax_in.bar(res, dists, width=np.diff(res).min() * 0.7,
              color="lightcoral", edgecolor="gray", lw=0.5)
    ax_in.axvline(chosen, color="red", ls="--", lw=1)
    ax_in.set_xticks(res)
    ax_in.set_xticklabels([f"{r:.1f}" for r in res], fontsize=6, rotation=45)
    ax_in.set_yticks([])
    ax_in.set_title("knee distance", fontsize=7)

    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  → Open {out.name} to confirm or choose a different resolution")


def plot_umap_clusters(adata, obs_key: str, title: str, out: Path) -> None:
    n = adata.obs[obs_key].nunique()
    fig, ax = plt.subplots(figsize=(6, 5))
    sc.pl.umap(adata, color=obs_key, ax=ax, show=False, frameon=False,
               title=f"{title} ({n} clusters)", legend_loc="on data",
               legend_fontsize=6, size=6, alpha=0.7)
    # Rasterize the scatter points (keeps text/axes vector) — standard for
    # large-cell-count figures; keeps PDF size sane on 600K+ points.
    for coll in ax.collections:
        coll.set_rasterized(True)
    fig.tight_layout()
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_clustree(adata, resolution_df: pd.DataFrame, out: Path) -> None:
    """Lightweight clustree: nodes = clusters per resolution, edges = majority
    parent at the previous resolution."""
    keys = resolution_df["obs_key"].tolist()
    resolutions = resolution_df["resolution"].tolist()
    n_levels = len(keys)

    edges = []
    for i in range(n_levels - 1):
        cross = pd.crosstab(adata.obs[keys[i]], adata.obs[keys[i + 1]])
        for child in cross.columns:
            col = cross[child]
            edges.append((i, str(col.idxmax()), i + 1, str(child), int(col.max())))

    max_clusters = resolution_df["n_clusters"].max()
    fig, ax = plt.subplots(figsize=(max(8, max_clusters * 0.5), n_levels * 1.5))
    ax.set_xlim(-1, max_clusters + 1)
    ax.set_ylim(-0.5, n_levels - 0.5)
    ax.set_yticks(range(n_levels))
    ax.set_yticklabels([f"res={r:.1f}" for r in resolutions])
    ax.set_xticks([])
    ax.set_title("Clustree — cluster splitting across Leiden resolutions")

    positions = {}
    for i, key in enumerate(keys):
        for j, c in enumerate(sorted(adata.obs[key].unique(), key=int)):
            positions[(i, c)] = j

    for (lvl_from, c_from, lvl_to, c_to, n) in edges:
        x0 = positions.get((lvl_from, c_from), 0)
        x1 = positions.get((lvl_to, c_to), 0)
        ax.plot([x0, x1], [lvl_from, lvl_to], color="gray", alpha=0.5,
                lw=max(0.3, np.log1p(n) / 4), zorder=1)

    cmap = plt.get_cmap("tab20")
    scale = max(1, adata.n_obs / 2000)
    for (lvl, c), x in positions.items():
        n_cells = int((adata.obs[keys[lvl]] == c).sum())
        ax.scatter(x, lvl, s=max(30, min(400, n_cells / scale)),
                   color=cmap(int(c) % 20), zorder=2, edgecolors="k", lw=0.4)
        ax.text(x, lvl + 0.08, str(c), ha="center", va="bottom", fontsize=5)

    fig.tight_layout()
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_cluster_qc(adata, obs_key: str, out: Path) -> None:
    """Per-cluster violin of QC metrics. Junk clusters = low UMI / high %mt."""
    metrics = [m for m in ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
               if m in adata.obs.columns]
    if not metrics:
        print("  [skip] cluster_qc_metrics.png: QC metrics not in adata.obs")
        return
    adata.obs[obs_key] = adata.obs[obs_key].astype("category")
    fig, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 4))
    if len(metrics) == 1:
        axes = [axes]
    for ax, m in zip(axes, metrics):
        sc.pl.violin(adata, keys=m, groupby=obs_key, ax=ax, show=False,
                     rotation=90, stripplot=False)
        ax.set_title(m)
    fig.suptitle(f"Per-cluster QC ({obs_key})")
    fig.tight_layout()
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)


def plot_cluster_composition(adata, obs_key: str, out: Path) -> None:
    """Stacked bar: sample fraction per cluster. ~100% one sample = batch artifact."""
    ct = pd.crosstab(adata.obs[obs_key], adata.obs["sample_id"], normalize="index")
    fig, ax = plt.subplots(figsize=(max(6, 0.3 * len(ct)), 4))
    ct.plot(kind="bar", stacked=True, ax=ax, width=0.8,
            colormap="tab20", edgecolor="none", legend=True)
    ax.set_ylabel("fraction of cells")
    ax.set_title("Sample composition per cluster — single-sample clusters = potential batch artifact")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=6, ncol=2)
    plt.setp(ax.get_xticklabels(), rotation=0, fontsize=7)
    fig.tight_layout()
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def cluster_qc_table(adata, obs_key: str) -> pd.DataFrame:
    obs = adata.obs
    metrics = [m for m in ["n_genes_by_counts", "total_counts", "pct_counts_mt",
                            "pct_counts_hemo"] if m in obs.columns]
    has_cycling = "cycling" in obs.columns
    rows = []
    for c in sorted(obs[obs_key].unique(), key=int):
        mask = obs[obs_key] == c
        row = {"cluster": c, "n_cells": int(mask.sum()),
               "n_samples": obs.loc[mask, "sample_id"].nunique()}
        for m in metrics:
            row[f"median_{m}"] = round(float(obs.loc[mask, m].median()), 3)
        if has_cycling:
            # Fraction of this cluster that is cycling — a cluster near 1.0 here
            # may be a proliferating subset of an otherwise-known cell type
            # rather than a distinct identity. Flag for review at annotation.
            row["frac_cycling"] = round(
                float((obs.loc[mask, "cycling"] == "cycling").mean()), 3)
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Phase 6: Leiden clustering")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--resolution", type=float, default=None,
                        help="Override auto-selection. Inspect resolution_selection.png first.")
    args = parser.parse_args()

    print(f"\n=== Phase 6: Clustering ===")
    print(f"Config: {args.config}")

    cfg = load_config(args.config)
    seed = int(cfg.get("random_seed", 42))

    in_path = Path(cfg["results_dir"]) / "h5ad" / "06_integrated" / "all_samples.h5ad"
    if not in_path.is_file():
        sys.exit(f"ERROR: missing {in_path}. Run 05_integration.py first.")

    out_dir = Path(cfg["results_dir"]) / "h5ad" / "07_clustered"
    plot_dir = Path(cfg["results_dir"]) / "plots" / "06_clustering"
    table_dir = phase_table_dir(cfg, "06_clustering")
    for d in (out_dir, plot_dir, table_dir):
        d.mkdir(parents=True, exist_ok=True)

    forced_res = args.resolution or \
                 cfg.get("clustering", {}).get("resolution") or \
                 cfg.get("integration", {}).get("leiden_resolution")

    print(f"\n[1/4] Loading {in_path}...")
    adata = sc.read_h5ad(in_path)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    if "X_scVI" not in adata.obsm:
        sys.exit("ERROR: X_scVI not in adata.obsm. Re-run Phase 5.")

    # Binary cycling label derived from the per-cell phase (G1/S/G2M from Phase 4).
    # This is a LABEL only — it does not feed scVI or clustering; it just lets us
    # color UMAPs and compute a per-cluster cycling fraction (catches clusters
    # that are merely proliferating cells of an otherwise-known type). The
    # per-cell S_score/G2M_score/phase/cc_difference remain in obs untouched.
    if "phase" in adata.obs.columns:
        adata.obs["cycling"] = (
            adata.obs["phase"].astype(str).ne("G1")
            .map({True: "cycling", False: "non_cycling"})
            .astype("category")
        )
        frac = (adata.obs["cycling"] == "cycling").mean()
        print(f"  Added 'cycling' label: {frac:.1%} of cells cycling (S or G2M)")
    else:
        print("  [warn] 'phase' not in obs — skipping 'cycling' label. Re-run Phase 4 "
              "if you want per-cell cycle annotation.")

    # Reuse the neighbor graph computed in Phase 5 (Option B): clustering runs
    # on the SAME graph the UMAP was built from, so cluster labels and the
    # embedding are consistent. Phase 5 persists the graph in adata.obsp.
    if "neighbors" not in adata.uns or "connectivities" not in adata.obsp:
        sys.exit(
            "ERROR: no neighbor graph found in the integrated h5ad.\n"
            "  Phase 6 reuses the graph built in Phase 5 (Option B). Re-run "
            "Phase 5 (05_integration.py) — it computes and saves the graph."
        )
    n_neighbors = adata.uns["neighbors"]["params"].get("n_neighbors", "?")
    print(f"\n[2/4] Reusing Phase 5 neighbor graph (n_neighbors={n_neighbors})")

    print(f"\n[3/4] Multi-resolution Leiden sweep: {RESOLUTIONS}...")
    resolution_df = run_leiden(adata, RESOLUTIONS, seed)
    dists = knee_distances(resolution_df["resolution"].values,
                           resolution_df["n_clusters"].values)
    resolution_df["knee_distance"] = dists

    chosen_res, method = choose_resolution(resolution_df, forced=forced_res)
    print(f"  Chosen: res={chosen_res}  ({method})")

    # If chosen_res isn't in the sweep (e.g. manual --resolution 0.5), run it now
    swept = resolution_df["resolution"].values
    nearest = float(swept[np.argmin(np.abs(swept - chosen_res))])
    chosen_key = f"leiden_r{nearest:.1f}".replace(".", "_")
    if chosen_key not in adata.obs.columns:
        print(f"  Running Leiden at resolution {chosen_res} (not in sweep)...")
        chosen_key = f"leiden_r{chosen_res:.1f}".replace(".", "_")
        sc.tl.leiden(adata, resolution=chosen_res, random_state=seed, key_added=chosen_key, flavor="igraph", n_iterations=2, directed=False)

    adata.obs["leiden"] = adata.obs[chosen_key].astype("category").copy()
    n_chosen = adata.obs["leiden"].nunique()
    print(f"  Clusters at chosen resolution: {n_chosen}")
    # NOTE: UMAP is NOT recomputed here. The embedding (X_umap) comes from
    # Phase 5, built on the same neighbor graph clustering just used (Option B).

    resolution_df.to_csv(table_dir / "06_clustering_summary.csv", index=False)
    cluster_qc_table(adata, "leiden").to_csv(
        table_dir / "06_clustering_cluster_qc_per_resolution.csv", index=False)

    print(f"\n[4/4] Generating plots...")
    plot_resolution_selection(resolution_df, chosen_res, method,
                               plot_dir / "resolution_selection.png")

    for _, row in resolution_df.iterrows():
        res_str = f"{row['resolution']:.1f}".replace(".", "p")
        plot_umap_clusters(adata, row["obs_key"],
                           title=f"res={row['resolution']:.1f}",
                           out=plot_dir / f"umap_leiden_res{res_str}.png")

    plot_umap_clusters(adata, "leiden",
                       title=f"Leiden res={chosen_res} (chosen)",
                       out=plot_dir / "umap_leiden_chosen.png")
    # Cycling label UMAP — visual companion to frac_cycling in the QC table.
    if "cycling" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(6, 5))
        sc.pl.umap(adata, color="cycling", ax=ax, show=False, frameon=False,
                   title="Cycling (S/G2M) vs non-cycling (G1)", size=6, alpha=0.7,
                   palette={"cycling": "#d62728", "non_cycling": "#cccccc"})
        for coll in ax.collections:
            coll.set_rasterized(True)
        fig.tight_layout()
        fig.savefig(plot_dir / "umap_cycling.png", dpi=130, bbox_inches="tight")
        plt.close(fig)
    plot_clustree(adata, resolution_df, plot_dir / "clustree.png")
    plot_cluster_qc(adata, "leiden", plot_dir / "cluster_qc_metrics.png")
    plot_cluster_composition(adata, "leiden",
                              plot_dir / "cluster_composition_by_sample.png")

    adata.write_h5ad(out_dir / "all_samples.h5ad")

    print(f"\n  Written: {out_dir / 'all_samples.h5ad'}")
    print(f"  Clusters: {n_chosen} (res={chosen_res}, {method})")
    print(f"  Plots: {plot_dir}")
    print(f"\n✓ Phase 6 complete.")
    print(f"\nWorkflow:")
    print(f"  1. Open resolution_selection.png — confirm auto-pick or choose manually")
    print(f"  2. If overriding: re-run with --resolution <value>")
    print(f"  3. Check clustree.png, cluster_qc_metrics.png, cluster_composition_by_sample.png")
    print(f"\nNext step: Phase 7 annotation (07_annotation.py)\n")


if __name__ == "__main__":
    main()
