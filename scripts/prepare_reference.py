#!/usr/bin/env python
"""
prepare_reference.py — validate & prepare a labeled reference .h5ad for Phase 7c.

Tissue-agnostic. Takes ANY reference AnnData you've obtained (ABC Atlas subset,
a published atlas, your own labeled data), checks it has the columns the tissue
YAML declares (labels_key, and region_key if regional claims are wanted),
optionally subsets it to keep it small, and writes it to the path the YAML's
reference.ref_h5ad expects.

This does NOT auto-download. Reference atlases are large (multi-GB) and their
metadata joins are source-specific; a silent auto-downloader would be fragile
and untestable. Instead: you obtain the reference, this script makes sure it's
actually usable by Phase 7c BEFORE you spend GPU hours — failing loudly if not.

How to obtain references (guidance, not run here):
  Brain — Allen Brain Cell Atlas via abc_atlas_access (workstation only):
    uv pip install "git+https://github.com/alleninstitute/abc_atlas_access.git"
    Then in Python:
      from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
      cache = AbcProjectCache.from_cache_dir("data/abc_atlas")
      meta  = cache.get_metadata_dataframe(
                  directory="WMB-10X",
                  file_name="cell_metadata_with_cluster_annotation")
      expr_path = cache.get_file_path(directory="WMB-10Xv3-<region>",
                                      file_name="WMB-10Xv3-<region>/raw")
      # read expr_path with anndata, join `meta` onto .obs so the cell-type and
      # region columns are present, write out a single labeled .h5ad, then point
      # this script at it with --source.
  Placenta — no standard spatial atlas; use any labeled placenta reference
    (e.g. Marsh & Blelloch) with a cell-type column; leave region_key null.

Usage:
  uv run python scripts/prepare_reference.py --config config/brain.yaml \\
      --source /path/to/raw_reference.h5ad
  # optional subsampling to keep the reference small/fast:
  uv run python scripts/prepare_reference.py --config config/brain.yaml \\
      --source /path/to/raw_reference.h5ad --max-cells-per-label 2000

Reads from the tissue YAML `reference:` block:
  - ref_h5ad   : destination path (where Phase 7c will look)
  - labels_key : required cell-type column that must exist in --source .obs
  - region_key : if non-null, must also exist in --source .obs
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

from _utils import load_config


def main():
    parser = argparse.ArgumentParser(description="Validate & prepare a reference h5ad for Phase 7c")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path,
                        help="Path to the raw reference .h5ad you obtained")
    parser.add_argument("--max-cells-per-label", type=int, default=None,
                        help="Optional: downsample to at most N cells per cell type")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    print(f"\n=== Prepare reference for Phase 7c ===")
    cfg = load_config(args.config)
    tissue = cfg["tissue"]

    ref_block = cfg.get("reference")
    if ref_block is None:
        sys.exit("ERROR: no 'reference:' block in config. Nothing to prepare for.")

    dest = ref_block.get("ref_h5ad")
    if not dest:
        sys.exit(
            "ERROR: reference.ref_h5ad is null in the YAML. Set it to the path where\n"
            "  you want the prepared reference written (e.g. data/refs/brain_ref.h5ad),\n"
            "  then re-run."
        )
    labels_key = ref_block.get("labels_key")
    region_key = ref_block.get("region_key")   # may be None
    if not labels_key:
        sys.exit("ERROR: reference.labels_key is required.")

    if not args.source.is_file():
        sys.exit(f"ERROR: --source not found: {args.source}")

    print(f"  Tissue:     {tissue}")
    print(f"  Source:     {args.source}")
    print(f"  Destination:{dest}")
    print(f"  labels_key: {labels_key}")
    print(f"  region_key: {region_key if region_key else '(none — cell-type only)'}")

    print(f"\n[1/4] Loading source...")
    adata = sc.read_h5ad(args.source)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # --- Validate required columns (hard fail, loud) ---
    print(f"\n[2/4] Validating columns...")
    missing = []
    if labels_key not in adata.obs.columns:
        missing.append(labels_key)
    if region_key and region_key not in adata.obs.columns:
        missing.append(region_key)
    if missing:
        sys.exit(
            f"ERROR: source reference is missing required .obs column(s): {missing}\n"
            f"  Available columns: {list(adata.obs.columns)}\n"
            f"  Either fix the source (join the cell-type/region metadata onto .obs)\n"
            f"  or update labels_key/region_key in {args.config}."
        )
    n_labels = adata.obs[labels_key].nunique()
    print(f"  ✓ labels_key '{labels_key}': {n_labels} cell types")
    if region_key:
        n_regions = adata.obs[region_key].nunique()
        print(f"  ✓ region_key '{region_key}': {n_regions} regions")

    # --- Sanity: raw counts in .X (scVI/scANVI need them) ---
    print(f"\n[3/4] Checking expression matrix...")
    X = adata.X
    # Heuristic: integer-valued max and no negatives ⇒ raw counts. Loud warning
    # (not silent) if it looks log-normalized, since scANVI wants raw counts.
    sample = X[:1000] if adata.n_obs > 1000 else X
    sample = sample.toarray() if hasattr(sample, "toarray") else np.asarray(sample)
    looks_integer = np.allclose(sample, np.round(sample)) and sample.min() >= 0
    if not looks_integer:
        print(f"  [warn] .X does not look like raw integer counts (max={sample.max():.2f}).")
        print(f"         scANVI in Phase 7c expects RAW COUNTS in .X. If this reference")
        print(f"         is log-normalized, point .X at the raw layer before proceeding.")
        print(f"         Continuing, but verify — this is the one thing that silently")
        print(f"         degrades label transfer quality.")
    else:
        print(f"  ✓ .X looks like raw counts")

    # --- Optional subsample per label ---
    if args.max_cells_per_label:
        print(f"\n[4/4] Subsampling to ≤{args.max_cells_per_label} cells per '{labels_key}'...")
        rng = np.random.default_rng(args.seed)
        keep_idx = []
        for lab, idx in adata.obs.groupby(labels_key, observed=True).indices.items():
            if len(idx) > args.max_cells_per_label:
                idx = rng.choice(idx, size=args.max_cells_per_label, replace=False)
            keep_idx.append(idx)
        keep_idx = np.sort(np.concatenate(keep_idx))
        adata = adata[keep_idx].copy()
        print(f"  Subsampled to {adata.n_obs:,} cells")
    else:
        print(f"\n[4/4] No subsampling (pass --max-cells-per-label to enable)")

    # --- Write ---
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(dest_path)
    print(f"\n✓ Reference prepared and written to: {dest_path}")
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes, {n_labels} cell types")
    print(f"  Phase 7c will now find it at the path in reference.ref_h5ad.")
    print(f"\nNext: uv run python scripts/07c_label_transfer.py --config {args.config}\n")


if __name__ == "__main__":
    main()
