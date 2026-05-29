#!/usr/bin/env python
"""
03_doublets.py — Phase 3 doublet detection.

Groups samples by `pool` (multiplexed library), combines their post-QC matrices,
runs scDblFinder once per pool (via R subprocess), splits results back to
samples, removes called doublets, writes a new .h5ad per sample.

Why per pool: doublets form within a single physical capture. With 10x Flex,
samples are multiplexed into the same pool, so the relevant "library" for
scDblFinder is the pool, not the sample. Passing `samples=` to scDblFinder
ensures simulated doublets respect within-sample boundaries (avoids fake
cross-sample doublets).

Usage:
  uv run python scripts/03_doublets.py --config config/dev.yaml
  uv run python scripts/03_doublets.py --config config/brain.yaml

Inputs:
  Per-sample h5ad in {results_dir}/h5ad/03_qc_filtered/{sample_id}.h5ad
  (output of 02_qc.py)

Outputs:
  Per-sample h5ad in {results_dir}/h5ad/04_doublets_removed/{sample_id}.h5ad
  Plots in           {results_dir}/plots/03_doublets/
  Summary CSV at     {results_dir}/tables/summary_doublets.csv
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp
import yaml


# ----------------------------------------------------------------------------
# Config loader — same as other phases
# ----------------------------------------------------------------------------

def load_config(path: Path) -> dict:
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
    return cfg


# ----------------------------------------------------------------------------
# R subprocess: prepare inputs, call Rscript, read TSV back
# ----------------------------------------------------------------------------

def run_scdblfinder_for_pool(pool_id: str, adatas: dict, rscript_path: Path,
                              seed: int) -> pd.DataFrame:
    """Combine same-pool samples, call scDblFinder via R, return per-cell results.

    adatas: dict of {sample_id: AnnData} for samples in this pool.
    Returns DataFrame with columns [barcode, sample_id, doublet_score, doublet_class].
    Barcodes are kept as the original .obs_names — we don't add a sample suffix here,
    so the join back on the Python side uses (sample_id, barcode).
    """
    # Concatenate along cells. inner join on genes — they should match (same Flex panel)
    # but inner is defensive in case any sample has been gene-filtered differently.
    sids = sorted(adatas.keys())
    combined = ad.concat(
        [adatas[s] for s in sids], axis=0, join="inner",
        keys=sids, label="_pool_sid", index_unique=None,
    )
    n_cells = combined.n_obs
    n_genes = combined.n_vars
    print(f"  [{pool_id}] combined {len(sids)} samples → {n_cells} cells × {n_genes} genes")

    # Write matrix (genes x cells, MM format) + barcode/feature/sample TSVs
    with tempfile.TemporaryDirectory(prefix=f"scdbl_{pool_id}_") as td:
        td = Path(td)
        mtx_path = td / "counts.mtx"
        bc_path = td / "barcodes.tsv"
        feat_path = td / "features.tsv"
        samp_path = td / "samples.tsv"
        out_path = td / "doublets.tsv"

        # Transpose to genes x cells; ensure sparse + integer-friendly
        X = combined.X
        if not sp.issparse(X):
            X = sp.csr_matrix(X)
        X_T = X.T.tocoo()  # genes x cells, COO for MM writer
        sio.mmwrite(str(mtx_path), X_T, field="integer")

        bc_path.write_text("\n".join(combined.obs_names) + "\n")
        feat_path.write_text("\n".join(combined.var_names) + "\n")
        samp_path.write_text("\n".join(combined.obs["_pool_sid"].astype(str)) + "\n")

        cmd = [
            "Rscript", str(rscript_path),
            "--matrix",   str(mtx_path),
            "--barcodes", str(bc_path),
            "--features", str(feat_path),
            "--samples",  str(samp_path),
            "--output",   str(out_path),
            "--seed",     str(seed),
        ]
        print(f"  [{pool_id}] calling Rscript...")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            print("---- R stdout ----\n" + proc.stdout)
            print("---- R stderr ----\n" + proc.stderr)
            sys.exit(f"ERROR: scDblFinder failed for pool {pool_id}")
        # Echo R's progress lines so the user can see what happened
        for line in proc.stdout.splitlines():
            print(f"    {line}")

        return pd.read_csv(out_path, sep="\t")


# ----------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------

def plot_doublet_score_hist(scores: pd.Series, classes: pd.Series,
                             pool_id: str, out: Path) -> None:
    """Histogram of doublet scores colored by class."""
    fig, ax = plt.subplots(figsize=(7, 4))
    bins = np.linspace(0, 1, 50)
    for cls, color in [("singlet", "steelblue"), ("doublet", "salmon")]:
        ax.hist(scores[classes == cls], bins=bins, color=color, alpha=0.7,
                label=f"{cls} (n={int((classes==cls).sum())})", edgecolor="k", lw=0.3)
    ax.set_xlabel("scDblFinder.score")
    ax.set_ylabel("n cells")
    ax.set_title(f"Pool {pool_id} — doublet score distribution")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


def plot_rate_per_sample(summary: pd.DataFrame, out: Path) -> None:
    """Bar chart of doublet rate per sample."""
    fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(summary)), 4))
    ax.bar(summary["sample_id"], summary["pct_doublet"], color="salmon", edgecolor="k")
    ax.set_ylabel("% doublets")
    ax.set_title("Doublet rate per sample")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Phase 3: doublet detection (scDblFinder)")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--rscript", type=Path, default=Path("scripts/run_scdblfinder.R"),
                       help="Path to the scDblFinder R script")
    args = parser.parse_args()

    print(f"\n=== Phase 3: Doublet detection ===")
    print(f"Config: {args.config}")

    if not args.rscript.is_file():
        sys.exit(f"ERROR: Rscript not found: {args.rscript}")
    if shutil.which("Rscript") is None:
        sys.exit("ERROR: Rscript not in PATH. Install R (brew install r) and retry.")

    cfg = load_config(args.config)
    results_dir = Path(cfg["results_dir"])
    in_dir = results_dir / "h5ad" / "03_qc_filtered"
    out_dir = results_dir / "h5ad" / "04_doublets_removed"
    plot_dir = results_dir / "plots" / "03_doublets"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    seed = cfg.get("random_seed", 42)

    # Group samples by pool, load post-QC AnnDatas
    by_pool: dict[str, dict[str, ad.AnnData]] = defaultdict(dict)
    for s in cfg["samples"]:
        sid, pool = s["id"], s["pool"]
        h5ad_path = in_dir / f"{sid}.h5ad"
        if not h5ad_path.is_file():
            sys.exit(f"ERROR: missing input {h5ad_path}. Run 02_qc.py first.")
        adata = sc.read_h5ad(h5ad_path)
        by_pool[pool][sid] = adata
    print(f"\nGrouping: {len(by_pool)} pool(s)")
    for p, d in by_pool.items():
        print(f"  {p}: {sorted(d.keys())}")

    # Run scDblFinder per pool, then split + write per-sample outputs
    rows = []
    for pool_id, adatas in by_pool.items():
        results = run_scdblfinder_for_pool(pool_id, adatas, args.rscript, seed)
        plot_doublet_score_hist(
            results["doublet_score"], results["doublet_class"],
            pool_id, plot_dir / f"pool_{pool_id}_score_hist.png",
        )

        # Split results back to samples and write filtered h5ads
        for sid, adata in adatas.items():
            sub = results[results["sample_id"] == sid].set_index("barcode")
            # Reindex to .obs_names order (defensive: scDblFinder doesn't reorder, but be safe)
            sub = sub.reindex(adata.obs_names)
            assert sub["doublet_class"].notna().all(), \
                f"barcode mismatch for {sid}: scDblFinder result missing cells"

            adata.obs["doublet_score"] = sub["doublet_score"].values
            adata.obs["doublet_class"] = sub["doublet_class"].values

            n_pre = adata.n_obs
            keep = adata.obs["doublet_class"] == "singlet"
            adata_filt = adata[keep].copy()
            n_post = adata_filt.n_obs
            n_dbl = n_pre - n_post

            adata_filt.write_h5ad(out_dir / f"{sid}.h5ad")
            print(f"  [{sid}] removed {n_dbl} doublets, kept {n_post}/{n_pre} ({100*n_post/n_pre:.1f}%)")
            rows.append({
                "sample_id": sid, "pool": pool_id, "n_pre": n_pre,
                "n_doublets": int(n_dbl), "n_post": n_post,
                "pct_doublet": round(100 * n_dbl / n_pre, 2),
            })

    summary = pd.DataFrame(rows)
    summary_path = results_dir / "tables" / "summary_doublets.csv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, index=False)
    plot_rate_per_sample(summary, plot_dir / "doublet_rate_per_sample.png")

    print(f"\n  Summary:")
    print(summary.to_string(index=False))
    print(f"\n  Filtered h5ads: {out_dir}")
    print(f"  Plots:          {plot_dir}")
    print(f"  Summary CSV:    {summary_path}")
    print(f"\n✓ Phase 3 complete.\n")


if __name__ == "__main__":
    main()
