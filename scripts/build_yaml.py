#!/usr/bin/env python
"""
build_yaml.py — generate config/brain.yaml and config/placenta.yaml from sample_metadata.csv

Run from the Analysis/ directory:
    uv run python scripts/build_yaml.py

This script is idempotent — re-run it any time the metadata CSV changes.
Writes two YAML files (one per tissue) with:
  - sample manifest (id, paths, group, age, sex, pool, donor_id, library)
  - sex marker genes for the Y-chromosome inferred sex check (Phase 0)
  - QC thresholds
  - Random seed

Dev-mode subsetting lives in config/dev.yaml (hand-maintained), which pulls
sample records from brain.yaml via `samples_from:` and selects an explicit
subset by `subset.sample_ids`.
"""

from pathlib import Path
import sys

import pandas as pd
import yaml


# -----------------------------------------------------------------------------
# Paths and config — adjust if you move things around
# -----------------------------------------------------------------------------
METADATA_CSV = Path("sample_metadata.csv")
CONFIG_DIR = Path("config")


# -----------------------------------------------------------------------------
# Shared config blocks (same for brain and placenta)
# -----------------------------------------------------------------------------
SHARED_CONFIG = {
    "qc": {
        "pct_mt_max": 1.0,          # snRNA: should be near zero
        "pct_hemo_max": 5.0,        # critical for placenta
        "n_mads": 5,                # MAD-based gene/UMI thresholds
    },
    "sex_markers": {
        "y_linked": ["Ddx3y", "Uty", "Eif2s3y", "Kdm5d"],
        "x_linked": ["Xist"],
    },
    "random_seed": 42,
}


def build_sample_entry(row: pd.Series) -> dict:
    """Convert one CSV row into one sample dict for the YAML.

    The sex field is normalized: 'TBD' -> 'unknown' (will be inferred in Phase 0).
    """
    sex = row["sex_declared"]
    if sex == "TBD":
        sex = "unknown"
    return {
        "id": row["sample_id"],
        "donor_id": row["donor_id"],
        "h5": row["h5_path"],
        "raw_h5": row["raw_h5_path"],
        "age": row["age"],
        "group": row["group"],
        "sex": sex,
        "pool": row["pool"],
        "library": row["library"],
    }


def build_yaml_for_tissue(df: pd.DataFrame, tissue: str) -> dict:
    """Build the config dict for one tissue."""
    sub = df[df["tissue"] == tissue].copy()

    cfg = {
        "tissue": tissue,
        "group_reference": "Relaxed",       # Relaxed is baseline; +logFC = upregulated in stress
        "results_dir": f"results/{tissue}",
        "samples": [build_sample_entry(r) for _, r in sub.iterrows()],
        **SHARED_CONFIG,
    }
    return cfg


def main():
    if not METADATA_CSV.is_file():
        sys.exit(f"ERROR: {METADATA_CSV} not found. Run from the Analysis/ directory.")

    df = pd.read_csv(METADATA_CSV)
    CONFIG_DIR.mkdir(exist_ok=True)

    for tissue in ["brain", "placenta"]:
        cfg = build_yaml_for_tissue(df, tissue)
        out = CONFIG_DIR / f"{tissue}.yaml"
        with out.open("w") as f:
            # default_flow_style=False -> block style (one key per line)
            # sort_keys=False -> preserve our intentional ordering
            yaml.safe_dump(cfg, f, default_flow_style=False, sort_keys=False, width=120)
        n = len(cfg["samples"])
        print(f"  Wrote {out}  ({n} samples)")

    # Summary
    print()
    print("Summary:")
    for tissue in ["brain", "placenta"]:
        sub = df[df["tissue"] == tissue]
        print(f"  {tissue} ({len(sub)} samples):")
        print(sub.groupby(["age", "group"]).size().to_string().replace("\n", "\n    "))
        print()


if __name__ == "__main__":
    main()
