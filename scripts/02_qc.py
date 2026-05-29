#!/usr/bin/env python
"""
02_qc.py — Phase 2 per-sample QC.

Runs AFTER Phase 0 validation. On the laptop, this consumes raw Cell Ranger h5
files directly (no CellBender available locally). On the workstation, swap the
input path in the config to point at CellBender-corrected files.

What it does, per sample:
  1. Load 10x h5
  2. (Dev mode only) subsample cells to subset.max_cells_per_sample
  3. Compute QC metrics: n_genes, total_counts, pct_counts_mt/ribo/hemo, top-20
  4. Determine thresholds:
       - Hard caps from config: pct_mt_max, pct_hemo_max
       - MAD-based per-sample: median ± n_mads * MAD on n_genes and total_counts
  5. Filter cells, save filtered .h5ad
  6. Write violin (pre/post), scatter, threshold-histogram plots
  7. Append a row to summary_qc.csv

Usage:
  uv run python scripts/02_qc.py --config config/dev.yaml
  uv run python scripts/02_qc.py --config config/brain.yaml

Outputs (in {results_dir}/h5ad/03_qc_filtered/ and {results_dir}/plots/02_qc/):
  - {sample_id}.h5ad                 : filtered AnnData per sample
  - {sample_id}_violin_prepost.png   : QC metric violins, before vs. after
  - {sample_id}_scatter.png          : counts vs n_genes, threshold lines
  - {sample_id}_thresholds.png       : histograms with cutoffs marked
  - summary_qc.csv                   : pre/post counts + applied thresholds
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml


# ----------------------------------------------------------------------------
# Config loader (matches the one in 01_validate.py — kept in sync deliberately)
# ----------------------------------------------------------------------------

def load_config(path: Path) -> dict:
    """Same load logic as 01_validate.py: handles samples_from + subset.sample_ids."""
    with path.open() as f:
        cfg = yaml.safe_load(f)
    if "samples_from" in cfg:
        with Path(cfg["samples_from"]).open() as f:
            src = yaml.safe_load(f)
        cfg["samples"] = src["samples"]
    subset = cfg.get("subset", {})
    if subset.get("enabled", False):
        ids = set(subset.get("sample_ids", []))
        before = len(cfg["samples"])
        cfg["samples"] = [s for s in cfg["samples"] if s["id"] in ids]
        missing = ids - {s["id"] for s in cfg["samples"]}
        if missing:
            sys.exit(f"ERROR: subset.sample_ids not in manifest: {sorted(missing)}")
        print(f"  Subset: {len(cfg['samples'])}/{before} samples")
    cwd = Path.cwd()
    for s in cfg["samples"]:
        h5 = Path(s["h5"])
        if not h5.is_absolute():
            s["h5"] = str((cwd / h5).resolve())
    return cfg


# ----------------------------------------------------------------------------
# QC core
# ----------------------------------------------------------------------------

def annotate_gene_categories(adata) -> None:
    """Flag mt/ribo/hemo genes in adata.var. Mouse gene naming.
    Bias toward lowercase-prefixed (mouse) but accept uppercase as a fallback."""
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "MT-"))
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl", "RPS", "RPL"))
    adata.var["hemo"] = adata.var_names.str.startswith(("Hbb", "Hba", "HBB", "HBA"))


def compute_thresholds(adata, qc_cfg: dict) -> dict:
    """Determine per-sample QC thresholds. Returns a dict the filter step consumes.

    Hard caps (from config, snRNA-specific):
      pct_mt_max     — typically 1.0 (nuclei should have near-zero %mt)
      pct_hemo_max   — typically 5.0 (catches hemoglobin contamination)

    Hard floors (optional; AND'd with MAD bounds — cell must pass BOTH):
      min_counts     — absolute UMI floor (e.g. 500). Catches debris even if a
                       sample's distribution is shifted so its MAD lower bound
                       lands below this. Omit from config to skip.
      min_genes      — absolute gene floor (e.g. 200). Same logic.

    MAD-based (per-sample, adaptive):
      n_genes:       [median - n_mads*MAD, median + n_mads*MAD]
      total_counts:  same, in log1p space (counts are heavily right-skewed)
    """
    n_mads = qc_cfg.get("n_mads", 5)

    def mad_bounds(x, log_space=False):
        x = np.log1p(x) if log_space else np.asarray(x)
        med = np.median(x)
        mad = np.median(np.abs(x - med))
        lo, hi = med - n_mads * mad, med + n_mads * mad
        if log_space:
            lo, hi = np.expm1(lo), np.expm1(hi)
        return float(max(lo, 0)), float(hi)

    n_genes_lo, n_genes_hi = mad_bounds(adata.obs["n_genes_by_counts"])
    counts_lo, counts_hi = mad_bounds(adata.obs["total_counts"], log_space=True)

    return {
        "pct_mt_max": float(qc_cfg.get("pct_mt_max", 1.0)),
        "pct_hemo_max": float(qc_cfg.get("pct_hemo_max", 5.0)),
        "n_genes_lo": n_genes_lo,
        "n_genes_hi": n_genes_hi,
        "total_counts_lo": counts_lo,
        "total_counts_hi": counts_hi,
        # Hard floors — None means "don't apply". The filter step uses max() of
        # these and the MAD bounds, so a cell must clear BOTH.
        "min_counts_floor": qc_cfg.get("min_counts"),
        "min_genes_floor": qc_cfg.get("min_genes"),
        "n_mads": n_mads,
    }


def apply_filters(adata, thr: dict):
    """Return (filtered_adata, mask) where mask is a per-cell bool array (True = keep).

    Hard floors (min_counts_floor, min_genes_floor) AND'd with MAD bounds:
    a cell must clear BOTH the MAD lower bound AND the absolute floor.
    """
    o = adata.obs
    # Effective lower bounds: max of MAD bound and absolute floor (if set)
    n_genes_lo = max(thr["n_genes_lo"], thr["min_genes_floor"] or 0)
    counts_lo  = max(thr["total_counts_lo"], thr["min_counts_floor"] or 0)
    mask = (
        (o["pct_counts_mt"] <= thr["pct_mt_max"])
        & (o["pct_counts_hemo"] <= thr["pct_hemo_max"])
        & (o["n_genes_by_counts"] >= n_genes_lo)
        & (o["n_genes_by_counts"] <= thr["n_genes_hi"])
        & (o["total_counts"] >= counts_lo)
        & (o["total_counts"] <= thr["total_counts_hi"])
    )
    return adata[mask].copy(), mask


# ----------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------

def plot_violin_prepost(adata_pre, adata_post, sample_id: str, out: Path) -> None:
    """Side-by-side violins of the 4 key QC metrics, pre vs post filter."""
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_hemo"]
    fig, axes = plt.subplots(2, 4, figsize=(14, 6))
    for j, m in enumerate(metrics):
        for i, (a, label) in enumerate([(adata_pre, "pre"), (adata_post, "post")]):
            ax = axes[i, j]
            ax.violinplot(a.obs[m].values, showmedians=True)
            ax.set_title(f"{label}: {m}")
            ax.set_xticks([])
    fig.suptitle(f"{sample_id} — QC pre/post  (pre n={adata_pre.n_obs}, post n={adata_post.n_obs})")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


def plot_scatter(adata_pre, thr: dict, sample_id: str, out: Path) -> None:
    """total_counts vs n_genes, colored by %mt, with threshold lines."""
    fig, ax = plt.subplots(figsize=(7, 6))
    sc_obj = ax.scatter(
        adata_pre.obs["total_counts"], adata_pre.obs["n_genes_by_counts"],
        c=adata_pre.obs["pct_counts_mt"], s=2, alpha=0.5, cmap="viridis",
    )
    plt.colorbar(sc_obj, ax=ax, label="pct_counts_mt")
    ax.axhline(thr["n_genes_lo"], color="red", ls="--", lw=0.8, label=f"n_genes ∈ [{thr['n_genes_lo']:.0f}, {thr['n_genes_hi']:.0f}]")
    ax.axhline(thr["n_genes_hi"], color="red", ls="--", lw=0.8)
    ax.axvline(thr["total_counts_lo"], color="orange", ls="--", lw=0.8, label=f"counts ∈ [{thr['total_counts_lo']:.0f}, {thr['total_counts_hi']:.0f}]")
    ax.axvline(thr["total_counts_hi"], color="orange", ls="--", lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("total_counts (UMI per cell)")
    ax.set_ylabel("n_genes_by_counts")
    ax.set_title(f"{sample_id} — counts vs genes (pre-filter)")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


def plot_thresholds(adata_pre, thr: dict, sample_id: str, out: Path) -> None:
    """Four histograms with threshold lines marked."""
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    specs = [
        ("n_genes_by_counts", [thr["n_genes_lo"], thr["n_genes_hi"]], "red"),
        ("total_counts",      [thr["total_counts_lo"], thr["total_counts_hi"]], "orange"),
        ("pct_counts_mt",     [thr["pct_mt_max"]], "purple"),
        ("pct_counts_hemo",   [thr["pct_hemo_max"]], "brown"),
    ]
    for ax, (metric, cuts, color) in zip(axes.flat, specs):
        vals = adata_pre.obs[metric].values
        # Log-scale x for counts/genes; linear for percentages
        if metric in ("total_counts", "n_genes_by_counts"):
            bins = np.logspace(np.log10(max(vals.min(), 1)), np.log10(vals.max()), 60)
            ax.set_xscale("log")
        else:
            bins = 60
        ax.hist(vals, bins=bins, color="lightgray", edgecolor="k", lw=0.3)
        for c in cuts:
            ax.axvline(c, color=color, ls="--", lw=1)
        ax.set_title(metric)
    fig.suptitle(f"{sample_id} — QC distributions + thresholds")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Per-sample driver
# ----------------------------------------------------------------------------

def process_sample(sample: dict, qc_cfg: dict, max_cells: int | None,
                   seed: int, out_h5ad: Path, out_plot: Path) -> dict:
    """Run QC for one sample. Returns a summary row dict."""
    sid = sample["id"]
    print(f"  [{sid}] loading {Path(sample['h5']).name}")
    adata = sc.read_10x_h5(sample["h5"])
    adata.var_names_make_unique()

    # Attach metadata up front — downstream phases expect these on .obs
    for k in ["donor_id", "age", "group", "sex", "pool", "library"]:
        adata.obs[k] = sample[k]
    adata.obs["sample_id"] = sid

    # Dev mode: subsample cells before any computation (saves RAM on laptop)
    n_pre_subset = adata.n_obs
    if max_cells is not None and adata.n_obs > max_cells:
        rng = np.random.default_rng(seed)
        keep = rng.choice(adata.n_obs, size=max_cells, replace=False)
        adata = adata[np.sort(keep)].copy()
        print(f"  [{sid}] subsampled {n_pre_subset} -> {adata.n_obs} cells (dev mode)")

    # QC metrics
    annotate_gene_categories(adata)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hemo"],
        percent_top=[20], log1p=False, inplace=True,
    )
    adata_pre = adata.copy()  # for plots

    # Thresholds + filter
    thr = compute_thresholds(adata, qc_cfg)
    adata_post, mask = apply_filters(adata, thr)
    n_pre, n_post = adata_pre.n_obs, adata_post.n_obs
    print(f"  [{sid}] kept {n_post}/{n_pre} cells ({100*n_post/n_pre:.1f}%)")

    # Save filtered AnnData
    out_h5ad.mkdir(parents=True, exist_ok=True)
    adata_post.write_h5ad(out_h5ad / f"{sid}.h5ad")

    # Plots
    out_plot.mkdir(parents=True, exist_ok=True)
    plot_violin_prepost(adata_pre, adata_post, sid, out_plot / f"{sid}_violin_prepost.png")
    plot_scatter(adata_pre, thr, sid, out_plot / f"{sid}_scatter.png")
    plot_thresholds(adata_pre, thr, sid, out_plot / f"{sid}_thresholds.png")

    return {
        "sample_id": sid, "n_pre": n_pre, "n_post": n_post,
        "pct_kept": round(100 * n_post / n_pre, 2),
        "median_umi_post": float(np.median(adata_post.obs["total_counts"])),
        "median_genes_post": float(np.median(adata_post.obs["n_genes_by_counts"])),
        "pct_mt_median_post": float(np.median(adata_post.obs["pct_counts_mt"])),
        "pct_hemo_median_post": float(np.median(adata_post.obs["pct_counts_hemo"])),
        **thr,
    }


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Phase 2: per-sample QC")
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()

    print(f"\n=== Phase 2: Per-sample QC ===")
    print(f"Config: {args.config}")

    cfg = load_config(args.config)
    samples = cfg["samples"]
    print(f"Tissue: {cfg['tissue']}")
    print(f"Samples: {len(samples)}")

    results_dir = Path(cfg["results_dir"])
    out_h5ad = results_dir / "h5ad" / "03_qc_filtered"
    out_plot = results_dir / "plots" / "02_qc"

    # Dev cell-cap (None for full runs)
    max_cells = cfg.get("subset", {}).get("max_cells_per_sample")
    seed = cfg.get("random_seed", 42)

    print(f"\nProcessing {len(samples)} samples...")
    rows = []
    for s in samples:
        row = process_sample(s, cfg["qc"], max_cells, seed, out_h5ad, out_plot)
        rows.append(row)

    summary = pd.DataFrame(rows)
    summary_path = results_dir / "tables" / "summary_qc.csv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, index=False)

    print(f"\n  Summary:")
    print(summary[["sample_id", "n_pre", "n_post", "pct_kept",
                   "median_umi_post", "pct_mt_median_post"]].to_string(index=False))
    print(f"\n  Filtered h5ads:  {out_h5ad}")
    print(f"  Plots:           {out_plot}")
    print(f"  Summary CSV:     {summary_path}")
    print(f"\n✓ Phase 2 complete.\n")


if __name__ == "__main__":
    main()
