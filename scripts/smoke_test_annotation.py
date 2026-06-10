"""
smoke_test_annotation.py

Smoke test for reference-based correlation labelling of placenta snRNA-seq.

For each Leiden cluster in the integrated h5ad:
  1. Compute per-cell Spearman correlation against each STAMP reference type
     (over the union of significant DEGs, ~6800 genes).
  2. Per cell: pick top match, compute gap to runner-up.
  3. Per cluster: aggregate to fractions table + majority label + mean gap +
     low-confidence flag (gap<0.05 or purity<0.5).

Outputs:
  - stdout: per-cluster summary table
  - tables/07_annotation_smoke/cluster_summary.csv
  - tables/07_annotation_smoke/cluster_fractions.csv  (cluster x ref_type, %)

This script does NOT write to obs or modify the h5ad. It's diagnostic only.

Run on workstation. Takes ~minutes for ~400k cells with the vectorized
Spearman implementation.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
import scanpy as sc
import h5py
from scipy import sparse


# ---------------------------- helpers --------------------------------------- #

def load_reference(path: Path) -> pd.DataFrame:
    """Load STAMP reference matrix (ref_types x genes of log2FC)."""
    with h5py.File(path, "r") as f:
        mat = f["log2fc"][:]
        types = [s.decode() if isinstance(s, bytes) else s for s in f["cell_types"][:]]
        genes = [s.decode() if isinstance(s, bytes) else s for s in f["genes"][:]]
    return pd.DataFrame(mat, index=pd.Index(types, name="cell_type"), columns=genes)


def rank_along_axis_fast(X: np.ndarray, axis: int = 1) -> np.ndarray:
    """
    Vectorized average-rank transform using scipy.stats.rankdata.
    Ties get mean rank (correct for Spearman).
    For snRNA-seq, most ties are zeros — rankdata handles this correctly.
    """
    from scipy.stats import rankdata
    return rankdata(X, method="average", axis=axis).astype(np.float32)


def spearman_matrix(query_ranked: np.ndarray, ref_ranked: np.ndarray) -> np.ndarray:
    """
    Spearman correlation between every query row and every reference row.

    query_ranked: (n_cells, n_genes) average ranks
    ref_ranked:   (n_types, n_genes) average ranks

    Returns (n_cells, n_types) correlations.

    Pearson on ranks = Spearman. We center and normalize each row, then matmul.
    """
    def standardize(R: np.ndarray) -> np.ndarray:
        R = R - R.mean(axis=1, keepdims=True)
        denom = np.linalg.norm(R, axis=1, keepdims=True)
        denom[denom == 0] = 1.0
        return (R / denom).astype(np.float32)

    Q = standardize(query_ranked)
    T = standardize(ref_ranked)
    return Q @ T.T  # (n_cells, n_types)


# ---------------------------- main ------------------------------------------ #

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", type=Path, required=True,
                    help="Integrated placenta h5ad with .obs[leiden_key] set")
    ap.add_argument("--reference", type=Path, required=True,
                    help="STAMP reference .h5 from build_placenta_reference.py")
    ap.add_argument("--leiden-key", type=str, default="leiden",
                    help="obs column with cluster labels (default: leiden)")
    ap.add_argument("--out-dir", type=Path, default=Path("tables/07_annotation_smoke"))
    ap.add_argument("--gap-threshold", type=float, default=0.05,
                    help="Per-cell gap below this = low-confidence (default: 0.05)")
    ap.add_argument("--purity-threshold", type=float, default=0.5,
                    help="Per-cluster majority fraction below this = low-confidence (default: 0.5)")
    ap.add_argument("--max-cells", type=int, default=None,
                    help="Subsample to this many cells for quick test (default: all)")
    args = ap.parse_args()

    t0 = time.time()

    print(f"Loading h5ad: {args.h5ad}", flush=True)
    adata = sc.read_h5ad(args.h5ad)
    print(f"  shape: {adata.shape}", flush=True)

    if args.leiden_key not in adata.obs:
        print(f"ERROR: obs['{args.leiden_key}'] not found. Available: "
              f"{list(adata.obs.columns)}", file=sys.stderr)
        return 1

    if args.max_cells and adata.n_obs > args.max_cells:
        print(f"Subsampling to {args.max_cells} cells", flush=True)
        sc.pp.subsample(adata, n_obs=args.max_cells, random_state=42)

    print(f"Loading reference: {args.reference}", flush=True)
    ref = load_reference(args.reference)
    print(f"  reference: {ref.shape[0]} cell types x {ref.shape[1]:,} genes",
          flush=True)

    # Resolve gene names — STAMP is symbols; check what adata.var_names looks like
    if adata.var_names[0].startswith("ENSMUS"):
        # need to map ensembl -> symbol via var['symbol']
        if "symbol" not in adata.var.columns:
            print("ERROR: adata uses Ensembl IDs but no 'symbol' column in .var",
                  file=sys.stderr)
            return 1
        symbol_map = adata.var["symbol"].astype(str)
        # Build a position lookup from symbol -> var_index
        sym_to_idx = pd.Series(np.arange(adata.n_vars), index=symbol_map.values)
        # Drop duplicate symbols (keep first)
        sym_to_idx = sym_to_idx[~sym_to_idx.index.duplicated(keep="first")]
        gene_in_query = sym_to_idx.index
    else:
        # assume symbols already
        gene_in_query = adata.var_names
        sym_to_idx = pd.Series(np.arange(adata.n_vars), index=gene_in_query)

    # Intersect reference genes with query genes
    common = ref.columns.intersection(gene_in_query)
    print(f"Common genes (ref ∩ query): {len(common):,} / {ref.shape[1]:,} ref genes",
          flush=True)
    if len(common) < 500:
        print(f"WARN: only {len(common)} common genes — annotation will be unreliable",
              file=sys.stderr)

    ref_common = ref[common].values.astype(np.float32)  # (n_types, n_genes)
    query_idx = sym_to_idx.loc[common].values
    # Extract query expression on the common genes
    # Use lognorm if available, else .X (assumed lognorm in Phase 5 output)
    from importlib import import_module
    utils_path = Path(__file__).parent / "_utils.py"
    if utils_path.exists():
        sys.path.insert(0, str(utils_path.parent))
        try:
            _u = import_module("_utils")
            _u.add_lognorm(adata)
            X = adata.layers["lognorm"]
        except Exception:
            X = adata.X
    else:
        X = adata.X

    # Slice columns
    if sparse.issparse(X):
        Xs = X[:, query_idx].toarray().astype(np.float32)
    else:
        Xs = X[:, query_idx].astype(np.float32)

    print(f"Query expression matrix: {Xs.shape}  (densified for ranking)",
          flush=True)
    print(f"  elapsed so far: {time.time()-t0:.1f}s", flush=True)

    # Rank-transform — chunk over cells to avoid blowing memory
    print("Rank-transforming reference and query (Spearman setup)...", flush=True)
    ref_ranked = rank_along_axis_fast(ref_common, axis=1)
    # query: chunked
    chunk = 5000
    n_cells = Xs.shape[0]
    corrs = np.empty((n_cells, ref.shape[0]), dtype=np.float32)
    for start in range(0, n_cells, chunk):
        end = min(start + chunk, n_cells)
        q_ranked = rank_along_axis_fast(Xs[start:end], axis=1)
        corrs[start:end] = spearman_matrix(q_ranked, ref_ranked)
        if (start // chunk) % 10 == 0:
            print(f"  ranked {end:>7,} / {n_cells:,}  "
                  f"({(time.time()-t0):.0f}s elapsed)", flush=True)

    print(f"Correlation matrix done: {corrs.shape} in {time.time()-t0:.1f}s",
          flush=True)

    # Per-cell top match + gap
    type_names = list(ref.index)
    top_idx = corrs.argmax(axis=1)
    top_corr = corrs.max(axis=1)
    # Runner-up: mask out top, max again
    corrs_for_gap = corrs.copy()
    corrs_for_gap[np.arange(n_cells), top_idx] = -np.inf
    runner_corr = corrs_for_gap.max(axis=1)
    gap = top_corr - runner_corr
    top_label = np.array(type_names, dtype=object)[top_idx]

    # Per-cluster aggregation
    clusters = adata.obs[args.leiden_key].astype(str).values

    # Fractions per cluster
    print("\nAggregating per cluster...", flush=True)
    df_cell = pd.DataFrame({
        "cluster": clusters,
        "top_match": top_label,
        "top_corr": top_corr,
        "gap": gap,
        "low_conf_cell": gap < args.gap_threshold,
    })

    # Cluster summary
    summary_rows = []
    fractions_rows = []
    for cl, sub in df_cell.groupby("cluster"):
        n = len(sub)
        counts = sub["top_match"].value_counts()
        fractions = counts / n
        majority = counts.index[0]
        purity = fractions.iloc[0]
        mean_gap = sub["gap"].mean()
        median_top_corr = sub["top_corr"].median()
        n_low_conf = sub["low_conf_cell"].sum()
        low_conf = (purity < args.purity_threshold) or (mean_gap < args.gap_threshold)
        flag = "low_confidence" if low_conf else ""
        # add runner-up to summary
        runner_up = counts.index[1] if len(counts) > 1 else ""
        runner_frac = fractions.iloc[1] if len(counts) > 1 else 0.0
        summary_rows.append({
            "cluster": cl,
            "n_cells": n,
            "majority_label": majority,
            "purity": round(purity, 3),
            "runner_up": runner_up,
            "runner_up_frac": round(runner_frac, 3),
            "mean_gap": round(mean_gap, 4),
            "median_top_corr": round(median_top_corr, 3),
            "n_low_conf_cells": int(n_low_conf),
            "pct_low_conf": round(100.0 * n_low_conf / n, 1),
            "flag": flag,
        })
        # Fractions for top 5 types
        for ct, frac in fractions.head(5).items():
            fractions_rows.append({
                "cluster": cl,
                "cell_type": ct,
                "fraction": round(float(frac), 4),
                "n_cells": int(counts[ct]),
            })

    summary = pd.DataFrame(summary_rows).sort_values(
        "n_cells", ascending=False
    ).reset_index(drop=True)
    fractions = pd.DataFrame(fractions_rows)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.out_dir / "cluster_summary.csv", index=False)
    fractions.to_csv(args.out_dir / "cluster_fractions.csv", index=False)

    # Print to stdout — the diagnostic the user actually needs to see
    print(f"\n=== Per-cluster summary ({len(summary)} clusters) ===")
    with pd.option_context("display.max_rows", None,
                           "display.max_columns", None,
                           "display.width", 200):
        print(summary.to_string(index=False))

    print(f"\nWrote:")
    print(f"  {args.out_dir/'cluster_summary.csv'}")
    print(f"  {args.out_dir/'cluster_fractions.csv'}")
    print(f"\nTotal time: {time.time()-t0:.1f}s")

    # Top-line health checks
    n_flagged = (summary["flag"] == "low_confidence").sum()
    print(f"\nHealth check:")
    print(f"  clusters flagged low-confidence: {n_flagged}/{len(summary)}")
    print(f"  unique labels assigned: {summary['majority_label'].nunique()}")
    if n_flagged / max(len(summary), 1) > 0.3:
        print("  WARN: >30% of clusters low-confidence — review gap/purity thresholds")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
