#!/usr/bin/env python
"""
prepare_brain_reference.py — build the brain reference for Phase 7 / 7c.

Downloads the Allen Brain Cell Atlas (WMB-10Xv3, ~2.3M adult mouse brain
nuclei, 10x v3 chemistry — closest match to the user's 10x Flex), subsamples
per `subclass` label for tractable size, concatenates regional h5ads into one
labeled reference, and trains a CellTypist .pkl on it.

Outputs (per default paths; override with --ref-out / --pkl-out):
  refs/abc_brain_ref.h5ad   — labeled reference used by 07c_label_transfer.py
  refs/celltypist_brain_adult.pkl — adult-brain model used by 07_annotation.py
                                    for both 4W and 3mo (one model fits both;
                                    ABC atlas is adult P56, biologically close).

Notes:
  - Skips download for regions whose h5ad is already in the cache.
  - Skips the whole step if both outputs already exist (idempotent).
  - WMB-10Xv3 is split by anatomical dissection region; we discover the
    region list from the cache rather than hardcoding it.
  - Labels: uses the `subclass` level of the WMB taxonomy (~330 cell types,
    sweet spot between resolution and CellTypist trainability for a 600K-cell
    query). Override with --label-key.
  - Genes are written as gene SYMBOLS in var_names (CellTypist trains on
    symbols; the query data also uses symbols).
  - Trained CellTypist model uses log2(CPM+1)-style normalization that
    matches ABC atlas's `log2` files — but we train from raw counts → lognorm
    inside this script (CellTypist's `train()` expects lognorm input).

Usage:
  uv run python scripts/prepare_brain_reference.py
  # or with options:
  uv run python scripts/prepare_brain_reference.py \\
      --cache-dir refs/abc_atlas --cells-per-label 300 \\
      --label-key subclass

API verified against:
  https://alleninstitute.github.io/abc_atlas_access/
  (AbcProjectCache, get_metadata_dataframe, get_file_path, list_data_files)
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cache-dir", type=Path, default=Path("refs/abc_atlas"),
                   help="Where AbcProjectCache stores downloaded files.")
    p.add_argument("--ref-out", type=Path, default=Path("refs/abc_brain_ref.h5ad"),
                   help="Destination for the labeled reference h5ad.")
    p.add_argument("--pkl-out", type=Path, default=Path("refs/celltypist_brain_adult.pkl"),
                   help="Destination for the trained CellTypist .pkl.")
    p.add_argument("--label-key", default="subclass",
                   choices=["class", "subclass", "supertype", "cluster"],
                   help="Which taxonomy level to use as labels (default: subclass).")
    p.add_argument("--cells-per-label", type=int, default=300,
                   help="Max cells per label after subsampling (default 300).")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--force", action="store_true",
                   help="Re-run even if outputs already exist.")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)

    print(f"\n{'='*60}")
    print("prepare_brain_reference.py — ABC atlas → reference + CellTypist .pkl")
    print(f"{'='*60}")

    # -- Idempotency ---------------------------------------------------------
    if not args.force and args.ref_out.exists() and args.pkl_out.exists():
        print(f"\n  Both outputs already exist (use --force to rebuild):")
        print(f"    {args.ref_out}")
        print(f"    {args.pkl_out}")
        print(f"  Skipping.")
        return

    # -- Imports (heavy; defer until needed) ---------------------------------
    try:
        from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
    except ImportError:
        sys.exit(
            "ERROR: abc_atlas_access not installed. From the main venv:\n"
            "  uv pip install 'git+https://github.com/AllenInstitute/abc_atlas_access.git'"
        )
    try:
        import celltypist
    except ImportError:
        sys.exit(
            "ERROR: celltypist not installed. Should be in pyproject.toml main deps.\n"
            "  uv add celltypist"
        )

    # -- Cache setup ---------------------------------------------------------
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n[1/6] Initialising AbcProjectCache at {args.cache_dir}")
    abc_cache = AbcProjectCache.from_cache_dir(args.cache_dir)
    print(f"  Manifest: {abc_cache.current_manifest}")

    # -- Cell metadata with cluster annotation (small, ~few hundred MB) ------
    print(f"\n[2/6] Downloading cell metadata + cluster annotations...")
    cell_meta = abc_cache.get_metadata_dataframe(
        directory="WMB-10X",
        file_name="cell_metadata_with_cluster_annotation",
        dtype={"cell_label": str},
    )
    cell_meta.set_index("cell_label", inplace=True)
    print(f"  Total cells in WMB-10X metadata: {len(cell_meta):,}")

    # Restrict to 10Xv3 cells (matches user's 10x Flex chemistry more closely
    # than 10Xv2; also reduces download volume).
    if "library_method" in cell_meta.columns:
        v3 = cell_meta[cell_meta["library_method"].str.contains("v3", na=False)]
    elif "feature_matrix_label" in cell_meta.columns:
        v3 = cell_meta[cell_meta["feature_matrix_label"].str.contains("10Xv3", na=False)]
    else:
        # UNVERIFIED — column name may differ; we keep all 10X cells if not present.
        print(f"  [warn] No library_method or feature_matrix_label column found.")
        print(f"         Keeping ALL cells (will mix 10Xv2 and 10Xv3).")
        v3 = cell_meta
    print(f"  After 10Xv3 filter: {len(v3):,} cells")

    if args.label_key not in v3.columns:
        sys.exit(
            f"ERROR: label key '{args.label_key}' not in cell metadata columns.\n"
            f"  Available: {sorted(v3.columns.tolist())}"
        )

    # -- Subsample at metadata level (before any expression download) --------
    print(f"\n[3/6] Subsampling to {args.cells_per_label} cells per "
          f"'{args.label_key}' (before download)...")
    keep_idx = []
    for lab, sub in v3.groupby(args.label_key, observed=True):
        if len(sub) > args.cells_per_label:
            picked = rng.choice(sub.index.values, size=args.cells_per_label, replace=False)
        else:
            picked = sub.index.values
        keep_idx.extend(picked)
    selected = v3.loc[keep_idx]
    n_labels = selected[args.label_key].nunique()
    print(f"  Selected {len(selected):,} cells across {n_labels} {args.label_key}s")

    # -- Group by region for targeted download -------------------------------
    if "feature_matrix_label" not in selected.columns:
        sys.exit(
            "ERROR: 'feature_matrix_label' column not in cell metadata.\n"
            "  Needed to map cells back to regional h5ad files."
        )
    by_region = selected.groupby("feature_matrix_label", observed=True)
    regions = sorted(by_region.groups.keys())
    print(f"\n[4/6] Need to download {len(regions)} regional expression matrices:")
    for r in regions:
        print(f"    {r}: {len(by_region.get_group(r)):,} cells")

    # -- Download + extract per region, concatenate --------------------------
    parts = []
    for i, region in enumerate(regions, 1):
        cells_in_region = by_region.get_group(region).index.values
        print(f"\n  [{i}/{len(regions)}] {region} ({len(cells_in_region):,} cells)...")

        # Cache directory is the chemistry (e.g. 'WMB-10Xv3'); the region is
        # encoded only in the file_name (e.g. 'WMB-10Xv3-PAL/raw').
        # feature_matrix_label like 'WMB-10Xv3-Isocortex-1' → chemistry
        # 'WMB-10Xv3' (first 2 hyphen-separated tokens).
        chemistry = "-".join(region.split("-")[:2])  # 'WMB-10Xv3'
        try:
            h5_path = abc_cache.get_file_path(
                directory=chemistry,
                file_name=f"{region}/raw",
            )
        except Exception as e:
            print(f"    [skip] could not download {region}: {type(e).__name__}: {e}")
            continue

        a = sc.read_h5ad(h5_path, backed="r")
        present = a.obs_names.intersection(cells_in_region)
        if len(present) == 0:
            print(f"    [skip] no cells from {region} found in h5ad index.")
            continue
        sub = a[present].to_memory().copy()  # bring just the selected cells in-memory
        parts.append(sub)
        print(f"    {sub.n_obs:,} cells extracted, {sub.n_vars:,} genes")

    if not parts:
        sys.exit("ERROR: no regions returned data. Aborting.")

    print(f"\n[5/6] Concatenating {len(parts)} regional pieces...")
    adata = ad.concat(parts, axis=0, join="outer", merge="first",
                      label="_region_source", keys=regions[:len(parts)])
    print(f"  Combined: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # Attach the rich label metadata we already have for these cells.
    # Dedupe columns (args.label_key may be 'subclass', already in the list)
    # and drop any colliding columns from existing obs — the metadata we're
    # joining is the authoritative source (some ABC regional h5ads ship with
    # 'class'/'subclass' already in .obs, which would cause duplicate columns).
    join_cols = list(dict.fromkeys([args.label_key, "class", "subclass"]))
    extra_cols = [c for c in ["region_of_interest_acronym",
                              "anatomical_division_label"]
                  if c in selected.columns]
    all_join_cols = list(dict.fromkeys(join_cols + extra_cols))

    collisions = [c for c in all_join_cols if c in adata.obs.columns]
    if collisions:
        print(f"  Dropping {len(collisions)} colliding obs columns "
              f"(will be replaced by metadata join): {collisions}")
        adata.obs = adata.obs.drop(columns=collisions)

    adata.obs = adata.obs.join(selected[all_join_cols], how="left")
    # Standardise the cell-type column name to what brain.yaml expects.
    adata.obs["cell_type"] = adata.obs[args.label_key].astype(str)
    if "anatomical_division_label" in adata.obs.columns:
        adata.obs["region"] = adata.obs["anatomical_division_label"].astype(str)
    elif "region_of_interest_acronym" in adata.obs.columns:
        adata.obs["region"] = adata.obs["region_of_interest_acronym"].astype(str)

    # Drop rows with missing labels (rare but possible after the join).
    n_pre = adata.n_obs
    adata = adata[adata.obs["cell_type"].notna() & (adata.obs["cell_type"] != "nan")]
    if adata.n_obs < n_pre:
        print(f"  Dropped {n_pre - adata.n_obs} cells with missing labels.")

    # Map gene IDs → symbols if needed. ABC atlas h5ads index .var by Ensembl ID
    # and carry gene_symbol in .var. CellTypist trains on the index of .var,
    # so we set the index to symbols.
    if "gene_symbol" in adata.var.columns:
        # Drop genes whose symbol is missing/dup; CellTypist needs unique names.
        sym = adata.var["gene_symbol"].astype(str)
        mask = sym.notna() & (sym != "nan") & ~sym.duplicated(keep="first")
        adata = adata[:, mask].copy()
        adata.var_names = adata.var["gene_symbol"].astype(str).values
        adata.var_names_make_unique()
        print(f"  Set var_names to gene symbols ({adata.n_vars:,} unique).")

    # Save reference h5ad (raw counts in .X, labeled).
    args.ref_out.parent.mkdir(parents=True, exist_ok=True)
    print(f"  Writing reference → {args.ref_out}")
    adata.write_h5ad(args.ref_out)

    # -- Train CellTypist .pkl -----------------------------------------------
    print(f"\n[6/6] Training CellTypist on reference (label='cell_type')...")
    # CellTypist expects lognorm in .X. Add it now.
    lognorm = adata.copy()
    sc.pp.normalize_total(lognorm, target_sum=1e4)
    sc.pp.log1p(lognorm)

    # Train. n_jobs=-1 uses all CPUs; mini_batch=True for memory safety on
    # large reference. SGD solver is faster for >100K cells.
    model = celltypist.train(
        lognorm, labels="cell_type",
        n_jobs=-1, use_SGD=True, mini_batch=True,
        check_expression=False,  # we already normalized
    )
    args.pkl_out.parent.mkdir(parents=True, exist_ok=True)
    model.write(args.pkl_out)
    print(f"  Saved CellTypist model → {args.pkl_out}")

    print(f"\n{'='*60}")
    print("✓ Brain reference build complete.")
    print(f"{'='*60}")
    print(f"  Reference h5ad: {args.ref_out}  ({adata.n_obs:,} cells)")
    print(f"  CellTypist:     {args.pkl_out}")
    print(f"\nNext: update config/brain.yaml reference + annotation blocks to point")
    print(f"here (run scripts/build_yaml.py to regenerate from sample_metadata.csv).")


if __name__ == "__main__":
    main()
