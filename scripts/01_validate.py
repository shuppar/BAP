#!/usr/bin/env python
"""
01_validate.py — Phase 0 validation. Run this BEFORE any compute-heavy phase.

What it does:
  1. Loads the YAML config and checks every sample's h5 file exists and is readable
  2. Checks required metadata fields are present, sample IDs unique
  3. Prints a balance table: samples per (age × condition × sex)
  4. Prints a library × condition × sex contingency table (confounding check)
  5. For each sample: loads counts, computes per-sample fingerprint
     (n_cells, median UMI, %mt, %hemo, %ribo, top 20 genes)
  6. Y-chromosome inferred sex per sample, compared to declared sex
  7. Writes everything to results/{tissue}/validation/ as CSVs + plots
  8. Exits 0 if all checks pass, non-zero if any fail (so it can gate downstream)

Usage:
  python scripts/01_validate.py --config config/brain.yaml
  python scripts/01_validate.py --config config/brain.yaml --no-fingerprints  # skip slow per-sample loading

Output files (in results/{tissue}/validation/):
  - manifest_balance.csv          : samples per (age, group, sex)
  - library_confound.csv          : library × group × sex table
  - sample_fingerprints.csv       : per-sample QC summary
  - sex_check.csv                 : declared vs. inferred sex per sample
  - balance_heatmap.png
  - library_confound_heatmap.png
  - sex_check_scatter.png
  - validation_report.txt         : human-readable summary
"""

import argparse
import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import yaml

# ----------------------------------------------------------------------------
# Helpers (kept short and inline — easier to read than importing from elsewhere)
# ----------------------------------------------------------------------------

REQUIRED_SAMPLE_FIELDS = ["id", "donor_id", "h5", "age", "group", "sex", "pool", "library"]
VALID_GROUPS = {"Early_Stress", "Late_Stress", "Relaxed"}
VALID_SEXES = {"M", "F", "unknown"}


def load_config(path: Path) -> dict:
    """Load YAML config. Returns a plain dict — no schema validation, fail-on-use.

    Supports the dev.yaml indirection pattern:
      - If the config has `samples_from: <other.yaml>`, samples are loaded from that file.
      - If `subset.enabled: true` with `subset.sample_ids: [...]`, samples are filtered
        to that ID list (order preserved from the source manifest).

    Relative h5 paths are resolved relative to the current working directory.
    """
    with path.open() as f:
        cfg = yaml.safe_load(f)

    # Indirection: pull sample records from another YAML if requested
    if "samples_from" in cfg:
        src = Path(cfg["samples_from"])
        with src.open() as f:
            src_cfg = yaml.safe_load(f)
        cfg["samples"] = src_cfg["samples"]

    # Subset by explicit ID list (dev mode)
    subset = cfg.get("subset", {})
    if subset.get("enabled", False):
        ids = set(subset.get("sample_ids", []))
        if not ids:
            sys.exit("ERROR: subset.enabled=true but subset.sample_ids is empty")
        before = len(cfg["samples"])
        cfg["samples"] = [s for s in cfg["samples"] if s["id"] in ids]
        missing = ids - {s["id"] for s in cfg["samples"]}
        if missing:
            sys.exit(f"ERROR: subset.sample_ids not found in manifest: {sorted(missing)}")
        print(f"  Subset: {len(cfg['samples'])}/{before} samples ({sorted(ids)})")

    # Resolve h5 paths
    cwd = Path.cwd()
    for s in cfg["samples"]:
        h5 = Path(s["h5"])
        if not h5.is_absolute():
            s["h5"] = str((cwd / h5).resolve())
    return cfg


def validate_manifest(samples: list[dict]) -> list[str]:
    """Return list of errors. Empty list = manifest is valid."""
    errors = []
    ids_seen = set()
    for i, s in enumerate(samples):
        # Required fields
        missing = [f for f in REQUIRED_SAMPLE_FIELDS if f not in s]
        if missing:
            errors.append(f"sample[{i}]: missing fields {missing}")
            continue
        # Duplicate IDs
        if s["id"] in ids_seen:
            errors.append(f"duplicate sample id: {s['id']}")
        ids_seen.add(s["id"])
        # Enum-like fields
        if s["group"] not in VALID_GROUPS:
            errors.append(f"{s['id']}: group={s['group']!r}, expected one of {VALID_GROUPS}")
        if s["sex"] not in VALID_SEXES:
            errors.append(f"{s['id']}: sex={s['sex']!r}, expected one of {VALID_SEXES}")
        # File existence
        if not Path(s["h5"]).is_file():
            errors.append(f"{s['id']}: h5 file not found at {s['h5']}")
    return errors


def balance_table(samples: list[dict]) -> pd.DataFrame:
    """Sample counts per (age, group, sex). Useful sanity-check on the design."""
    df = pd.DataFrame(samples)
    counts = df.groupby(["age", "group", "sex"]).size().reset_index(name="n_samples")
    return counts


def library_confound_table(samples: list[dict]) -> pd.DataFrame:
    """library × group × sex counts. If group perfectly co-varies with library,
    integration can't separate the two and the experiment is fundamentally confounded."""
    df = pd.DataFrame(samples)
    return pd.crosstab(df["library"], [df["group"], df["sex"]])


def compute_fingerprint(h5_path: Path, sample_id: str) -> dict:
    """Per-sample QC summary. Loads only what's needed; doesn't keep the AnnData."""
    adata = sc.read_10x_h5(str(h5_path))
    adata.var_names_make_unique()

    # Compute %mt, %ribo, %hemo using gene-name patterns
    adata.var["mt"] = adata.var_names.str.startswith("mt-") | adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl", "RPS", "RPL"))
    adata.var["hemo"] = adata.var_names.str.startswith(("Hbb", "Hba", "HBB", "HBA"))

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hemo"],
                                percent_top=[20], log1p=False, inplace=True)

    # Top 20 genes by total counts
    gene_sums = np.asarray(adata.X.sum(axis=0)).flatten()
    top_idx = np.argsort(gene_sums)[-20:][::-1]
    top_genes = adata.var_names[top_idx].tolist()

    return {
        "sample_id": sample_id,
        "n_cells": adata.n_obs,
        "n_genes_detected": int((gene_sums > 0).sum()),
        "median_umi": float(np.median(adata.obs["total_counts"])),
        "median_genes": float(np.median(adata.obs["n_genes_by_counts"])),
        "pct_mt_median": float(np.median(adata.obs["pct_counts_mt"])),
        "pct_ribo_median": float(np.median(adata.obs["pct_counts_ribo"])),
        "pct_hemo_median": float(np.median(adata.obs["pct_counts_hemo"])),
        "top_genes": ", ".join(top_genes[:5]),  # just top 5 for the CSV
    }


def infer_sex(h5_path: Path, y_markers: list[str], x_markers: list[str]) -> dict:
    """Infer sex from Y-linked + Xist expression per sample.
    Returns {y_score, x_score, inferred_sex}."""
    adata = sc.read_10x_h5(str(h5_path))
    adata.var_names_make_unique()

    # Score = log1p mean expression of marker set across cells
    def gene_score(genes):
        present = [g for g in genes if g in adata.var_names]
        if not present:
            return 0.0
        x = adata[:, present].X
        x = x.toarray() if hasattr(x, "toarray") else x
        return float(np.log1p(x.mean()))

    y_score = gene_score(y_markers)
    x_score = gene_score(x_markers)

    # Simple heuristic: high Y, low Xist -> M; low Y, high Xist -> F.
    # Boundary cases get flagged as "ambiguous" — usually means contamination.
    if y_score > 0.05 and x_score < 0.05:
        inferred = "M"
    elif y_score < 0.05 and x_score > 0.05:
        inferred = "F"
    else:
        inferred = "ambiguous"

    return {"y_score": y_score, "x_score": x_score, "inferred_sex": inferred}


def plot_balance(df: pd.DataFrame, out: Path) -> None:
    """Heatmap: rows=group×sex, cols=age, cells=n_samples."""
    pivot = df.pivot_table(index=["group", "sex"], columns="age",
                          values="n_samples", fill_value=0).astype(int)
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.heatmap(pivot, annot=True, fmt="d", cmap="Blues", ax=ax, cbar_kws={"label": "n samples"})
    ax.set_title("Sample balance: group × sex × age")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


def plot_library_confound(df: pd.DataFrame, out: Path) -> None:
    """Heatmap of library × (group, sex). Warn if any row has zeros."""
    fig, ax = plt.subplots(figsize=(8, max(3, 0.5 * len(df))))
    sns.heatmap(df, annot=True, fmt="d", cmap="Reds", ax=ax)
    ax.set_title("Library confound check (samples per library × group × sex)")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


def plot_sex_check(df: pd.DataFrame, out: Path) -> None:
    """Y-score vs Xist-score per sample, colored by declared sex.
    'unknown' declared (E12.5 placenta) plotted in gray — no mismatch check.
    Mismatches stand out: e.g., a 'declared M' point in the upper-left (high Xist, low Y)."""
    fig, ax = plt.subplots(figsize=(6, 5))
    for sex, color in zip(["M", "F", "unknown"], ["steelblue", "salmon", "lightgray"]):
        sub = df[df["declared_sex"] == sex]
        if len(sub):
            ax.scatter(sub["y_score"], sub["x_score"], c=color, label=f"declared={sex}", s=80, alpha=0.7)
    # Flag mismatches only where sex was declared (not 'unknown')
    declared = df[df["declared_sex"].isin(["M", "F"])]
    mismatch = declared[declared["declared_sex"] != declared["inferred_sex"]]
    if len(mismatch):
        ax.scatter(mismatch["y_score"], mismatch["x_score"],
                  marker="x", c="black", s=200, linewidths=3, label="mismatch")
        for _, row in mismatch.iterrows():
            ax.annotate(row["sample_id"], (row["y_score"], row["x_score"]),
                       xytext=(5, 5), textcoords="offset points", fontsize=8)
    ax.set_xlabel("Y-linked score (Ddx3y, Uty, ...)")
    ax.set_ylabel("X-inactivation score (Xist)")
    ax.set_title("Sex check: inferred vs declared")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Phase 0: validate manifest, sex, fingerprints")
    parser.add_argument("--config", required=True, type=Path, help="Path to YAML config")
    parser.add_argument("--no-fingerprints", action="store_true",
                       help="Skip per-sample fingerprint loading (just validate manifest)")
    args = parser.parse_args()

    print(f"\n=== Phase 0: Validation ===")
    print(f"Config: {args.config}")

    cfg = load_config(args.config)
    samples = cfg["samples"]
    print(f"Tissue: {cfg['tissue']}")
    print(f"Samples: {len(samples)}")

    # --- 1. Manifest validation ---
    print("\n[1/4] Validating manifest...")
    errors = validate_manifest(samples)
    if errors:
        print(f"  ✗ {len(errors)} error(s):")
        for e in errors:
            print(f"     - {e}")
        sys.exit(1)
    print(f"  ✓ All {len(samples)} samples have required fields and h5 files exist")

    # --- 2. Output dir ---
    out_dir = Path(cfg["results_dir"]) / "validation"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"  Output dir: {out_dir}")

    # --- 3. Balance + confound tables ---
    print("\n[2/4] Computing balance + library-confound tables...")
    bal = balance_table(samples)
    bal.to_csv(out_dir / "manifest_balance.csv", index=False)
    print(f"  Balance table ({len(bal)} groups):")
    print(bal.to_string(index=False))
    plot_balance(bal, out_dir / "balance_heatmap.png")

    confound = library_confound_table(samples)
    confound.to_csv(out_dir / "library_confound.csv")
    plot_library_confound(confound, out_dir / "library_confound_heatmap.png")
    # Flag perfect confounding: any library with only one (condition, sex) combo
    per_lib_combos = Counter()
    for s in samples:
        per_lib_combos[s["library"]] += 1
    single_combo_libs = [lib for lib, n in per_lib_combos.items()
                        if confound.loc[lib].astype(bool).sum() == 1]
    if single_combo_libs:
        print(f"  ⚠ Libraries with only ONE (group, sex) combination: {single_combo_libs}")
        print(f"     (this would confound library effects with biology — check carefully)")

    # --- 4. Per-sample fingerprints + sex check ---
    if args.no_fingerprints:
        print("\n[3/4] Skipping fingerprints (--no-fingerprints)")
    else:
        print(f"\n[3/4] Computing per-sample fingerprints ({len(samples)} samples)...")
        fingerprints = []
        sex_checks = []
        for s in samples:
            print(f"  Loading {s['id']}...")
            fp = compute_fingerprint(Path(s["h5"]), s["id"])
            fingerprints.append(fp)
            sx = infer_sex(Path(s["h5"]),
                          cfg["sex_markers"]["y_linked"],
                          cfg["sex_markers"]["x_linked"])
            sex_checks.append({"sample_id": s["id"], "declared_sex": s["sex"], **sx})

        fp_df = pd.DataFrame(fingerprints)
        fp_df.to_csv(out_dir / "sample_fingerprints.csv", index=False)
        print(f"\n  Fingerprints:")
        print(fp_df[["sample_id", "n_cells", "median_umi", "pct_mt_median",
                     "pct_hemo_median"]].to_string(index=False))

        # --- 5. Sex check ---
        print(f"\n[4/4] Sex check (Y-linked vs Xist)...")
        sx_df = pd.DataFrame(sex_checks)
        sx_df.to_csv(out_dir / "sex_check.csv", index=False)
        # Only flag mismatches where sex was actually declared (not 'unknown' E12.5 placenta)
        declared = sx_df[sx_df["declared_sex"].isin(["M", "F"])]
        mismatches = declared[declared["declared_sex"] != declared["inferred_sex"]]
        n_unknown = (sx_df["declared_sex"] == "unknown").sum()
        if len(mismatches):
            print(f"  ⚠ {len(mismatches)} sex MISMATCH(es):")
            print(mismatches.to_string(index=False))
        else:
            print(f"  ✓ All {len(declared)} declared-sex samples: declared matches inferred")
        if n_unknown:
            print(f"  ℹ {n_unknown} sample(s) with declared sex='unknown' (E12.5 placenta) — sex inferred from Y/Xist")
        plot_sex_check(sx_df, out_dir / "sex_check_scatter.png")

    # --- Summary ---
    with open(out_dir / "validation_report.txt", "w") as f:
        f.write(f"Phase 0 validation: PASSED\n")
        f.write(f"Tissue: {cfg['tissue']}\n")
        f.write(f"Samples: {len(samples)}\n")
        f.write(f"Outputs: {out_dir}\n")

    print(f"\n✓ Validation complete. Outputs in: {out_dir}")
    print(f"  Review the plots and CSVs before launching Phase 1.\n")


if __name__ == "__main__":
    main()
