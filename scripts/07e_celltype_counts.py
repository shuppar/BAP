#!/usr/bin/env python
"""
07e_celltype_counts.py — per-donor × cell-type cell counts (diagnostic CSV).

Reads the Phase 7 integrated + annotated h5ad and writes a long-form CSV with
the absolute number of cells per (donor × annotation tier × level × cell type).
Companion to 07_annotation_class_per_cluster_age.csv (gate audit) and
08a_dropped_cells_per_donor.csv (dropped contaminants/unassigned). This file
is the raw count substrate that 8a's propeller test, 8f cross-tissue, and 8g
cross-age operate on — useful for QA, ad-hoc pivots, and supplementary tables.

Long form (one row per donor × granularity × level × region × celltype):

    tissue, donor_id, sample_id, age, group, sex, pool,
    granularity, level, region, celltype, n_cells, category

  - granularity:
        brain    -> celltypist_broad | celltypist_class | celltypist_subclass
        placenta -> celltype_majority
  - level:
        brain    -> whole (sums regions) + region (per celltypist_region value)
        placenta -> whole only
  - region:
        brain    -> 'whole' for level=whole, else the region label
        placenta -> 'whole' always
  - category:
        'assigned'   for real cell types
        'unassigned' for gate demotions (unassigned_immune / _glia / _vascular
                      / _erythroid) and sentinels (no_region_model,
                      no_subclass_model)

Notes
-----
- Includes ALL cells in the integrated h5ad (unassigned rows are preserved so
  the table can serve as a sanity check; filter by `category=='assigned'` if
  you want just the real types).
- Does NOT include `subcluster_name` (those labels live in the 7b focal-type
  h5ads, not in the main integrated object). Add a separate dump if needed.
- Uses `assigned_sex` if present (project source-of-truth), else falls back
  to `sex`.

Usage:
    uv run python scripts/07e_celltype_counts.py --config config/brain.yaml
    uv run python scripts/07e_celltype_counts.py --config config/placenta.yaml
    uv run python scripts/07e_celltype_counts.py --config config/brain.yaml \\
        --input results/brain/h5ad/07_annotated.h5ad
"""

import argparse
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

from _utils import load_config, phase_table_dir


BRAIN_GRANULARITIES    = ["celltypist_broad", "celltypist_class",
                          "celltypist_subclass"]
PLACENTA_GRANULARITIES = ["celltype_majority"]

UNASSIGNED_PREFIXES = ("unassigned",)
UNASSIGNED_LITERALS = {"unassigned", "no_region_model", "no_subclass_model"}


def is_unassigned(label) -> bool:
    if pd.isna(label):
        return True
    s = str(label)
    return s.startswith(UNASSIGNED_PREFIXES) or s in UNASSIGNED_LITERALS


def _donor_cols_present(obs_cols) -> list[str]:
    """Donor / metadata columns to carry through to every row. Uses
    assigned_sex (project source-of-truth) when available."""
    candidates = ["donor_id", "sample_id", "age", "group", "pool"]
    cols = [c for c in candidates if c in obs_cols]
    # sex: prefer assigned_sex, fall back to sex
    if "assigned_sex" in obs_cols:
        cols.append("assigned_sex")
    elif "sex" in obs_cols:
        cols.append("sex")
    return cols


def dump_counts(obs: pd.DataFrame, tissue: str,
                granularities: list[str],
                region_col: str | None) -> pd.DataFrame:
    """Build the long-form count table.

    region_col=None  -> whole-only (placenta).
    region_col=str   -> whole + per-region rows (brain).
    """
    donor_cols = _donor_cols_present(obs.columns)
    if "donor_id" not in donor_cols:
        sys.exit("ERROR: obs has no 'donor_id' column — can't build "
                 "per-donor counts.")
    print(f"  donor metadata columns: {donor_cols}")

    rows = []
    for gran in granularities:
        if gran not in obs.columns:
            print(f"  [skip] obs has no '{gran}' column")
            continue
        # whole-tissue counts: group by donor metadata + celltype
        g_whole = (obs.groupby(donor_cols + [gran], observed=True)
                      .size().reset_index(name="n_cells"))
        g_whole["granularity"] = gran
        g_whole["level"]       = "whole"
        g_whole["region"]      = "whole"
        g_whole = g_whole.rename(columns={gran: "celltype"})
        rows.append(g_whole)
        print(f"  {gran:<22s} whole: {len(g_whole):>7,} rows")

        if region_col and region_col in obs.columns:
            g_reg = (obs.groupby(donor_cols + [region_col, gran], observed=True)
                        .size().reset_index(name="n_cells"))
            g_reg["granularity"] = gran
            g_reg["level"]       = "region"
            g_reg = g_reg.rename(columns={region_col: "region", gran: "celltype"})
            rows.append(g_reg)
            print(f"  {gran:<22s} region: {len(g_reg):>7,} rows")
        elif region_col:
            print(f"  [note] region_col '{region_col}' not in obs; "
                  f"skipping per-region rows for {gran}")

    if not rows:
        sys.exit("ERROR: no granularities found in obs — nothing to dump.")

    out = pd.concat(rows, ignore_index=True)
    out["tissue"]   = tissue
    out["category"] = out["celltype"].apply(
        lambda x: "unassigned" if is_unassigned(x) else "assigned")

    # Normalize the sex column name
    if "assigned_sex" in out.columns and "sex" not in out.columns:
        out = out.rename(columns={"assigned_sex": "sex"})

    col_order = ["tissue", "donor_id", "sample_id", "age", "group", "sex",
                 "pool", "granularity", "level", "region", "celltype",
                 "n_cells", "category"]
    return out[[c for c in col_order if c in out.columns]]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--input", default=None, type=Path,
                    help="Override path to the integrated h5ad "
                         "(default: results/{tissue}/h5ad/07_annotated.h5ad)")
    args = ap.parse_args()

    print("\n=== 07e: per-donor cell-type counts ===")
    cfg = load_config(args.config)
    tissue = cfg.get("tissue")
    if tissue not in ("brain", "placenta"):
        sys.exit(f"ERROR: unknown tissue '{tissue}' in config")

    if args.input:
        h5ad_path = args.input
    else:
        # Phase 7 annotation has historically written to several different
        # subfolders. Try them in order; first hit wins.
        candidates = [
            Path(cfg["results_dir"]) / "h5ad" / "08_annotated" / "all_samples.h5ad",
            Path(cfg["results_dir"]) / "h5ad" / "07_annotated.h5ad",
            Path(cfg["results_dir"]) / "h5ad" / "07_annotated" / "all_samples.h5ad",
            Path(cfg["results_dir"]) / "h5ad" / "07_clustered" / "all_samples.h5ad",
        ]
        h5ad_path = next((p for p in candidates if p.is_file()), candidates[0])
    if not h5ad_path.is_file():
        sys.exit(f"ERROR: integrated h5ad not found.\n"
                 f"  Looked at:\n    "
                 + "\n    ".join(str(c) for c in [h5ad_path] if True)
                 + f"\n  Pass --input <path> to override.")

    print(f"Reading {h5ad_path} ...")
    adata = ad.read_h5ad(h5ad_path, backed="r")
    print(f"  {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    obs = adata.obs.copy()        # only .obs is needed; release the rest
    del adata

    if tissue == "brain":
        grans = BRAIN_GRANULARITIES
        region_col = "celltypist_region"
    else:
        grans = PLACENTA_GRANULARITIES
        region_col = None

    print(f"\nBuilding counts for granularities: {grans}")
    out = dump_counts(obs, tissue, grans, region_col)

    table_dir = phase_table_dir(cfg, "07_annotation")
    table_dir.mkdir(parents=True, exist_ok=True)
    out_path = table_dir / "07e_celltype_counts.csv"
    out.to_csv(out_path, index=False)
    print(f"\n  Wrote {len(out):,} rows -> {out_path}")
    print(f"  ({out_path.stat().st_size / 1e6:.2f} MB)")

    # Quick sanity summaries
    print("\nRows per (granularity × level):")
    print(out.groupby(["granularity", "level"], observed=True)
             .size().to_string())

    print("\nCells per (granularity × category):")
    cat_sum = (out.groupby(["granularity", "level", "category"], observed=True)
                  ["n_cells"].sum().unstack("category", fill_value=0))
    print(cat_sum.to_string())

    # Whole-tissue rows should sum to n_obs per granularity. Per-region rows
    # should also sum to n_obs per granularity (each cell has exactly one
    # region label). Confirm.
    print("\nSanity: per-granularity cell totals (each row should equal n_obs):")
    tot = (out.groupby(["granularity", "level"], observed=True)
              ["n_cells"].sum().unstack("level", fill_value=0))
    print(tot.to_string())

    print(f"\n  Tip: filter to assigned cells only with "
          f"`df[df['category']=='assigned']`\n"
          f"        pivot per-donor: "
          f"`df.pivot_table(index='donor_id', columns='celltype', "
          f"values='n_cells', aggfunc='sum', fill_value=0)`\n")


if __name__ == "__main__":
    main()
