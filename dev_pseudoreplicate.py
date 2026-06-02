#!/usr/bin/env python
"""
dev_pseudoreplicate.py — LOCAL DEV SMOKE-TEST HELPER. NOT part of the pipeline.

>>> DO NOT run this on the server / real data. DO NOT commit results from it. <<<

Purpose: the dev subset has 1 sample per group, so composition (8a) and
pseudobulk DE (8b) can't run — they need >=2-3 donors per group. This splits
each real sample's cells into N pseudo-donors so the DOWNSTREAM CODE PATHS can
be exercised locally before transferring the pipeline to the server.

It SPLITS (partitions) cells, it does NOT duplicate them: each pseudo-donor is a
distinct random subset of one real sample. That avoids the zero-variance
degeneracy that literal copies would cause in scCODA / DESeq2. Total cell count
is unchanged; cells are just relabeled with new donor_id / sample_id.

WHAT THIS DOES NOT DO: produce meaningful statistics. Pseudo-donors are not real
biological replicates — there is no true between-animal variance. Any "credible
effect" or "DEG" from pseudo-replicated data is an artifact. This only confirms
the scripts run and emit well-formed tables/plots.

Usage (from Analysis/):
  python dev_pseudoreplicate.py                       # 3 pseudo-donors, default paths
  python dev_pseudoreplicate.py --n 3 \\
      --in  results/dev/h5ad/08_annotated/all_samples.h5ad \\
      --out results/dev/h5ad/08_annotated/all_samples.h5ad   # overwrites (backs up first)

By default it backs up the input to <input>.orig_backup (once) and overwrites
the input path, so 08a/08b read it with no path changes. Restore with:
  cp results/dev/h5ad/08_annotated/all_samples.h5ad.orig_backup \\
     results/dev/h5ad/08_annotated/all_samples.h5ad
"""

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np
import scanpy as sc


def main():
    ap = argparse.ArgumentParser(description="DEV ONLY: split samples into pseudo-donors")
    ap.add_argument("--in", dest="inp", type=Path,
                    default=Path("results/dev/h5ad/08_annotated/all_samples.h5ad"))
    ap.add_argument("--out", dest="out", type=Path, default=None,
                    help="default: overwrite --in (after backing up to .orig_backup)")
    ap.add_argument("--n", type=int, default=3, help="pseudo-donors per real sample")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    print("=" * 70)
    print("DEV SMOKE-TEST: splitting samples into pseudo-donors.")
    print("Results are NOT statistically meaningful. Local testing only.")
    print("=" * 70)

    if not args.inp.is_file():
        sys.exit(f"ERROR: input not found: {args.inp}")
    out = args.out or args.inp

    adata = sc.read_h5ad(args.inp)
    print(f"\nLoaded {adata.n_obs:,} cells from {args.inp}")
    if "donor_id" not in adata.obs or "sample_id" not in adata.obs:
        sys.exit("ERROR: obs needs donor_id and sample_id.")

    rng = np.random.default_rng(args.seed)
    new_donor = np.empty(adata.n_obs, dtype=object)
    new_sample = np.empty(adata.n_obs, dtype=object)

    # Partition each real sample's cells into args.n pseudo-donors
    for sid in adata.obs["sample_id"].unique():
        idx = np.where(adata.obs["sample_id"].values == sid)[0]
        rng.shuffle(idx)
        parts = np.array_split(idx, args.n)
        for r, part in enumerate(parts, start=1):
            new_donor[part] = f"{sid}_r{r}"
            new_sample[part] = f"{sid}_r{r}"
        sizes = [len(p) for p in parts]
        print(f"  {sid}: {len(idx)} cells -> {args.n} pseudo-donors {sizes}")

    adata.obs["donor_id"] = new_donor
    adata.obs["sample_id"] = new_sample
    # Cast to category to match the rest of the pipeline's obs dtypes
    for c in ("donor_id", "sample_id"):
        adata.obs[c] = adata.obs[c].astype("category")

    n_donors = adata.obs["donor_id"].nunique()
    print(f"\nNow {n_donors} pseudo-donors "
          f"(was {len(adata.obs['sample_id'].unique()) // args.n} real samples).")
    # Sanity: per-group donor counts (what 8a/8b will see)
    if "group" in adata.obs:
        print("  donors per group:")
        gd = adata.obs.groupby("group", observed=True)["donor_id"].nunique()
        for g, n in gd.items():
            print(f"    {g}: {n}")

    # Back up the input once, then write
    if out == args.inp:
        backup = args.inp.with_suffix(args.inp.suffix + ".orig_backup")
        if not backup.exists():
            shutil.copy2(args.inp, backup)
            print(f"\nBacked up original -> {backup}")
        else:
            print(f"\nBackup already exists ({backup}); not overwriting it.")

    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    print(f"Wrote pseudo-replicated object -> {out}")
    print(f"\nNow run 8a/8b on dev with min_donors at its default. Remember: numbers")
    print(f"are meaningless; you're checking the code runs and output is well-formed.")
    print(f"Restore real data with:")
    print(f"  cp {out}.orig_backup {out}\n")


if __name__ == "__main__":
    main()
