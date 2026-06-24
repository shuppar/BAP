#!/usr/bin/env python
"""One-off patch: add specificity_fdr to existing 08e_lr_baseline.csv files.

The perms ran (cellphone_pvals present) but specificity_fdr was skipped due to a
column-name mismatch. This recomputes it by BH-correcting cellphone_pvals within
each (group, age, level) — no LIANA re-run. Idempotent (re-runnable).

Usage:
  uv run python scripts/fix_08e_specificity_fdr.py
"""
import sys
from pathlib import Path
import pandas as pd
from statsmodels.stats.multitest import multipletests

DIRS = [
    "results/placenta/tables/08e_communication/08e_lr_baseline.csv",
    "results/placenta/tables/08e_communication_subtype/08e_lr_baseline.csv",
    "results/brain/tables/08e_communication/08e_lr_baseline.csv",
    "results/brain/tables/08e_communication_subtype/08e_lr_baseline.csv",
]


def patch(path):
    p = Path(path)
    if not p.is_file():
        print(f"  SKIP (missing): {path}")
        return
    df = pd.read_csv(p, low_memory=False)
    if "cellphone_pvals" not in df.columns:
        print(f"  SKIP (no cellphone_pvals — perms off?): {path}")
        return
    grp_cols = [c for c in ("group", "age", "level") if c in df.columns]
    parts = []
    for _, sl in df.groupby(grp_cols):
        _, fdr, _, _ = multipletests(sl["cellphone_pvals"].fillna(1.0), method="fdr_bh")
        parts.append(pd.Series(fdr, index=sl.index))
    df["specificity_fdr"] = pd.concat(parts)
    df.to_csv(p, index=False)
    n_sig = int((df["specificity_fdr"] < 0.05).sum())
    print(f"  OK: {path}")
    print(f"      rows={len(df):,}  specificity_fdr<0.05={n_sig:,}  "
          f"(BH within {'+'.join(grp_cols)})")


if __name__ == "__main__":
    print("Adding specificity_fdr from cellphone_pvals (BH within group×age×level):")
    for d in DIRS:
        patch(d)
    print("Done.")
