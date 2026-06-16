#!/usr/bin/env python
"""
patch_broad_regionfree.py — ONE-TIME patch: recompute region-free celltypist_broad
on already-annotated h5ad objects, without re-running Phase 7.

Why: the first Phase 7 run produced region-TAGGED broad for adults ("Excitatory
neurons (CB)") but plain broad for P1 ("Excitatory neurons"), so broad didn't
align across ages. The permanent fix is in 07_annotation.py (derive_brain_broad
now strips the region suffix); this patches the existing object so we don't have
to re-run the ~1h Phase 7. Future runs are correct automatically.

Idempotent: re-running just rewrites the same region-free broad.

Parallelism: the broad map is a vectorized .map() (sub-second per object), so
the per-object work isn't worth threading. What IS parallel: patching multiple
tissues' objects concurrently (independent file read+write). With one tissue
(brain) it's effectively serial; pass --tissues brain placenta to do both at once.

Usage:
  uv run python scripts/patch_broad_regionfree.py --tissues brain
  uv run python scripts/patch_broad_regionfree.py --tissues brain placenta
"""

import argparse
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import scanpy as sc


def _coarse(b: str) -> str:
    return re.sub(r"\s*\([^)]*\)\s*$", "", str(b)).strip()


def _load_map_csv(path: Path, key_col: str, val_col: str) -> dict:
    if not path.is_file():
        sys.exit(f"ERROR: mapping CSV not found: {path}")
    df = pd.read_csv(path)
    return dict(zip(df[key_col].astype(str), df[val_col].astype(str)))


def patch_one(tissue: str, results_root: Path, config_dir: Path,
              abc_csv: Path, ros_csv: Path) -> str:
    h5 = results_root / tissue / "h5ad" / "08_annotated" / "all_samples.h5ad"
    if not h5.is_file():
        return f"[{tissue}] SKIP — no annotated object at {h5}"

    a = sc.read_h5ad(h5)
    if "celltypist_class" not in a.obs.columns:
        return f"[{tissue}] SKIP — no celltypist_class (not a CellTypist/scANVI tissue)"

    abc_map = _load_map_csv(abc_csv, "class", "broad") if abc_csv.is_file() else {}
    ros_map = _load_map_csv(ros_csv, "class", "broad") if ros_csv.is_file() else {}
    merged = {**abc_map, **ros_map}
    merged = {k: _coarse(v) for k, v in merged.items()}
    if not merged:
        return f"[{tissue}] SKIP — no class->broad maps found"

    before = (a.obs["celltypist_broad"].nunique()
              if "celltypist_broad" in a.obs.columns else None)
    cls = a.obs["celltypist_class"].astype(str)
    broad = cls.map(merged).fillna("unassigned")
    a.obs["celltypist_broad"] = pd.Categorical(broad)
    n_unmapped = int((broad == "unassigned").sum())

    a.write_h5ad(h5)
    cats = sorted(a.obs["celltypist_broad"].cat.categories.tolist())
    return (f"[{tissue}] patched {h5}\n"
            f"    broad classes: {before} -> {len(cats)} (region-free)\n"
            f"    {n_unmapped} unmapped->unassigned\n"
            f"    classes: {cats}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissues", nargs="+", default=["brain"])
    ap.add_argument("--results-root", type=Path, default=Path("results"))
    ap.add_argument("--config-dir", type=Path, default=Path("config"))
    ap.add_argument("--abc-csv", type=Path,
                    default=Path("refs/abc_class_to_broad.csv"))
    args = ap.parse_args()

    ros_csv = args.config_dir / "rosenberg_class_to_broad.csv"

    jobs = {}
    with ProcessPoolExecutor(max_workers=len(args.tissues)) as ex:
        for t in args.tissues:
            jobs[ex.submit(patch_one, t, args.results_root, args.config_dir,
                           args.abc_csv, ros_csv)] = t
        for fut in as_completed(jobs):
            print(fut.result())


if __name__ == "__main__":
    main()
