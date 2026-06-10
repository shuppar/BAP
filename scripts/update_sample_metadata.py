"""
update_sample_metadata.py — one-pass updater.

What it does (idempotent):
  1. Fix placenta path columns: 'per_sample_outs/XXX.N/' → 'per_sample_outs/XXX_N/'
     (dots in sample folder names → underscores, matching on-disk reality).
  2. Add/refresh `assigned_sex` column by merging from both
     results/{brain,placenta}/validation/sex_check.csv.
  3. Write updated sample_metadata.csv in place.

Run from project root after Phase 0 has produced both sex_check.csv files.

After running, re-run build_yaml.py to regenerate brain.yaml + placenta.yaml.
"""
import re
import sys
import pandas as pd
from pathlib import Path

META_PATH = Path("sample_metadata.csv")
BRAIN_SX = Path("results/brain/validation/sex_check.csv")
PLAC_SX = Path("results/placenta/validation/sex_check.csv")

if not META_PATH.exists():
    sys.exit(f"ERROR: {META_PATH} not found. Run from project root.")
if not BRAIN_SX.exists() or not PLAC_SX.exists():
    sys.exit(
        f"ERROR: sex_check.csv missing for one or both tissues:\n"
        f"  {BRAIN_SX}: {'present' if BRAIN_SX.exists() else 'MISSING'}\n"
        f"  {PLAC_SX}: {'present' if PLAC_SX.exists() else 'MISSING'}\n"
        f"Run Phase 0 first."
    )

meta = pd.read_csv(META_PATH)
print(f"[update] loaded {META_PATH} — {len(meta)} samples")

# --- 1. Path normalization (placenta dots → underscores) ---------------------
# Matches both h5_path and raw_h5_path columns, and the cellranger_folder
# column (just the folder name, e.g. "MES2.3" → "MES2_3").

def fix_subpath(s):
    """Replace dots between letters and digits in per_sample_outs sample folder."""
    if pd.isna(s):
        return s
    return re.sub(r"(per_sample_outs/[A-Za-z0-9]+)\.([0-9]+)/", r"\1_\2/", str(s))

def fix_folder(s):
    """Convert standalone 'XXX.N' → 'XXX_N' (cellranger_folder column)."""
    if pd.isna(s):
        return s
    return re.sub(r"^([A-Za-z0-9]+)\.([0-9]+)$", r"\1_\2", str(s))

n_before = (meta["h5_path"].astype(str).str.contains(r"per_sample_outs/[A-Za-z0-9]+\.[0-9]", regex=True)).sum()
meta["h5_path"] = meta["h5_path"].apply(fix_subpath)
meta["raw_h5_path"] = meta["raw_h5_path"].apply(fix_subpath)
meta["cellranger_folder"] = meta["cellranger_folder"].apply(fix_folder)
n_after = (meta["h5_path"].astype(str).str.contains(r"per_sample_outs/[A-Za-z0-9]+\.[0-9]", regex=True)).sum()
print(f"[update] paths normalized: {n_before} rows had dots in path → {n_after} remaining")

# --- 2. Merge assigned_sex from both sex_check.csv files ---------------------
brain_sx = pd.read_csv(BRAIN_SX)[["sample_id", "assigned_sex"]]
plac_sx = pd.read_csv(PLAC_SX)[["sample_id", "assigned_sex"]]
sx_all = pd.concat([brain_sx, plac_sx], ignore_index=True)

# Sanity: every metadata sample should appear in exactly one sex_check
n_missing = (~meta["sample_id"].isin(sx_all["sample_id"])).sum()
if n_missing:
    missing = meta.loc[~meta["sample_id"].isin(sx_all["sample_id"]), "sample_id"].tolist()
    print(f"[update] WARN: {n_missing} samples in metadata have no sex_check row: {missing}")

# Drop any pre-existing assigned_sex column (idempotency); then merge fresh.
if "assigned_sex" in meta.columns:
    meta = meta.drop(columns=["assigned_sex"])
meta = meta.merge(sx_all, on="sample_id", how="left")

n_assigned = meta["assigned_sex"].notna().sum()
n_swaps = (
    meta["sex_declared"].isin(["M", "F"]) &
    (meta["sex_declared"] != meta["assigned_sex"])
).sum()
print(f"[update] assigned_sex filled for {n_assigned}/{len(meta)} samples")
print(f"[update] declared/assigned mismatches: {n_swaps}")
if n_swaps:
    swap_rows = meta[
        meta["sex_declared"].isin(["M", "F"]) &
        (meta["sex_declared"] != meta["assigned_sex"])
    ][["sample_id", "tissue", "sex_declared", "assigned_sex"]]
    print(swap_rows.to_string(index=False))

# --- 3. Save in place --------------------------------------------------------
meta.to_csv(META_PATH, index=False)
print(f"[update] wrote {META_PATH} ({len(meta)} rows, {len(meta.columns)} columns)")
print(f"[update] new column order:")
print(f"  {list(meta.columns)}")
print(f"\nNext: uv run python scripts/build_yaml.py")
