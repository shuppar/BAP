#!/usr/bin/env python
"""
04_integration_prep.py — Phase 4: concatenate per tissue + normalize + HVG selection.

Combines all post-doublet sample .h5ads into one AnnData per tissue, log-normalizes
(in a layer — raw counts stay in .X for scVI), selects highly variable genes,
and flags HVG exclusions (hemoglobin, mito, ribo, sex-linked, plus pregnancy
genes for placenta).

Layer/slot conventions for the output AnnData:
  - .X                       : raw counts (sparse, int) — what scVI consumes
  - .layers["lognorm"]       : log1p(normalize_total(X, 1e4)) — for plotting/QC
  - .var["highly_variable"]  : bool, set by sc.pp.highly_variable_genes
  - .var["hvg_excluded"]     : bool, True if gene is in exclusion list
  - .var["use_for_scvi"]     : highly_variable AND NOT hvg_excluded

Usage:
  uv run python scripts/04_integration_prep.py --config config/dev.yaml
  uv run python scripts/04_integration_prep.py --config config/brain.yaml
  uv run python scripts/04_integration_prep.py --config config/placenta.yaml

Inputs:
  Per-sample h5ads in {results_dir}/h5ad/04_doublets_removed/{sample_id}.h5ad

Outputs:
  Concatenated h5ad at {results_dir}/h5ad/05_integration_ready/all_samples.h5ad
  Plots in             {results_dir}/plots/04_integration_prep/
  Summary table at     {results_dir}/tables/summary_integration_prep.csv
"""

import argparse
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml


# ----------------------------------------------------------------------------
# Config loader (same as other phases)
# ----------------------------------------------------------------------------

def load_config(path: Path) -> dict:
    with path.open() as f:
        cfg = yaml.safe_load(f)
    if "samples_from" in cfg:
        with Path(cfg["samples_from"]).open() as f:
            src = yaml.safe_load(f)
        cfg["samples"] = src["samples"]
    subset = cfg.get("subset", {})
    if subset.get("enabled", False):
        ids = set(subset.get("sample_ids", []))
        before = len(cfg["samples"])
        cfg["samples"] = [s for s in cfg["samples"] if s["id"] in ids]
        missing = ids - {s["id"] for s in cfg["samples"]}
        if missing:
            sys.exit(f"ERROR: subset.sample_ids not in manifest: {sorted(missing)}")
        print(f"  Subset: {len(cfg['samples'])}/{before} samples")
    return cfg


# ----------------------------------------------------------------------------
# HVG exclusion lists
# ----------------------------------------------------------------------------

# Gene-name patterns excluded from HVGs regardless of tissue. These genes are
# either technical artifacts (mt, ribo) or dominate variance for non-biological
# reasons we don't want driving integration (hemoglobin contamination, sex).
GENERIC_EXCLUDE_PREFIXES = ("mt-", "MT-", "Rps", "Rpl", "RPS", "RPL",
                            "Hbb", "Hba", "HBB", "HBA")
GENERIC_EXCLUDE_EXACT = {"Xist", "Ddx3y", "Uty", "Eif2s3y", "Kdm5d", "Tsix"}

# Placenta-specific: prolactin family + pregnancy-specific glycoproteins dominate
# variance in mouse placenta snRNA-seq and would otherwise eat HVG slots.
# (Prl3d1, Prl8a8, Prl7a1, ... — entire Prl* family. Psg* same story.)
PLACENTA_EXCLUDE_PREFIXES = ("Prl", "Psg", "Cgb")
PLACENTA_EXCLUDE_EXACT = {"Cga"}


def build_exclusion_mask(var_names: pd.Index, tissue: str) -> pd.Series:
    """Return bool Series, True = exclude from HVG. Logged for the user."""
    excl_prefixes = list(GENERIC_EXCLUDE_PREFIXES)
    excl_exact = set(GENERIC_EXCLUDE_EXACT)
    if tissue == "placenta":
        excl_prefixes += list(PLACENTA_EXCLUDE_PREFIXES)
        excl_exact |= PLACENTA_EXCLUDE_EXACT

    mask = var_names.str.startswith(tuple(excl_prefixes)) | var_names.isin(excl_exact)
    return pd.Series(mask, index=var_names, name="hvg_excluded")


# ----------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------

def plot_hvg_dispersion(adata, out: Path) -> None:
    """scanpy HVG plot — variance vs mean, HVGs highlighted."""
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig(out, dpi=120, bbox_inches="tight")
    plt.close()


def plot_cells_per_sample(adata, out: Path) -> None:
    """Bar of cells per sample in the concatenated object — sanity check."""
    counts = adata.obs["sample_id"].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(counts)), 4))
    ax.bar(counts.index, counts.values, color="steelblue", edgecolor="k")
    ax.set_ylabel("n cells")
    ax.set_title(f"Cells per sample after concat (total: {adata.n_obs:,})")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


def plot_excluded_summary(adata, out: Path) -> None:
    """How many HVGs were knocked out by the exclusion list."""
    n_hv = int(adata.var["highly_variable"].sum())
    n_use = int(adata.var["use_for_scvi"].sum())
    n_excl = n_hv - n_use
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(["raw HVGs", "after exclusion"], [n_hv, n_use],
           color=["lightcoral", "steelblue"], edgecolor="k")
    ax.set_ylabel("n genes")
    ax.set_title(f"HVG exclusion ({n_excl} removed)")
    for i, v in enumerate([n_hv, n_use]):
        ax.text(i, v, str(v), ha="center", va="bottom")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Phase 4: concat + normalize + HVG")
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()

    print(f"\n=== Phase 4: Integration prep (concat + normalize + HVG) ===")
    print(f"Config: {args.config}")

    cfg = load_config(args.config)
    tissue = cfg["tissue"]
    results_dir = Path(cfg["results_dir"])
    in_dir = results_dir / "h5ad" / "04_doublets_removed"
    out_dir = results_dir / "h5ad" / "05_integration_ready"
    plot_dir = results_dir / "plots" / "04_integration_prep"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # HVG count: tissue-specific default, YAML-overridable
    integ_cfg = cfg.get("integration", {})
    default_n_hvg = 2000 if tissue == "placenta" else 3000
    n_hvg = int(integ_cfg.get("n_hvg", default_n_hvg))
    print(f"Tissue: {tissue}  |  target n_hvg: {n_hvg}")

    # --- 1. Load + concat ---
    print(f"\n[1/4] Loading + concatenating {len(cfg['samples'])} samples...")
    adatas = {}
    for s in cfg["samples"]:
        path = in_dir / f"{s['id']}.h5ad"
        if not path.is_file():
            sys.exit(f"ERROR: missing input {path}. Run 03_doublets.py first.")
        adatas[s["id"]] = sc.read_h5ad(path)

    # Concat: outer join keeps all genes across samples. Should be identical for
    # Flex (same probe panel) but outer is defensive. index_unique="-" appends
    # the sample_id to barcodes so they're globally unique.
    combined = ad.concat(
        adatas, axis=0, join="outer", merge="same",
        label="sample_id_concat", index_unique="-",
    )
    # sample_id was already on .obs from 02_qc.py; the concat-added label is redundant
    combined.obs.drop(columns=["sample_id_concat"], errors="ignore", inplace=True)
    print(f"  Combined: {combined.n_obs:,} cells × {combined.n_vars:,} genes")

    # Ensure .X is sparse + integer counts (scVI requirement)
    if not sp.issparse(combined.X):
        combined.X = sp.csr_matrix(combined.X)
    # Cast to int32 if it's float (concat may upcast); scVI tolerates float but
    # int is cleaner and saves memory.
    if combined.X.dtype.kind == "f":
        combined.X = combined.X.astype(np.int32)

    # --- 2. Normalize → layer (keep raw in .X for scVI) ---
    print(f"\n[2/4] Log-normalizing → .layers['lognorm'] (raw counts stay in .X)")
    combined.layers["counts"] = combined.X.copy()  # explicit alias of raw counts
    # Do normalize/log1p on a temp copy so .X stays raw
    tmp = combined.copy()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    combined.layers["lognorm"] = tmp.X
    del tmp

    # --- 3. HVG selection (seurat_v3 on raw counts, batch-aware by pool) ---
    print(f"\n[3/4] Selecting {n_hvg} HVGs (seurat_v3, batch_key=pool)...")
    sc.pp.highly_variable_genes(
        combined, n_top_genes=n_hvg, flavor="seurat_v3",
        batch_key="pool", layer="counts",
    )
    n_hv_raw = int(combined.var["highly_variable"].sum())
    print(f"  Raw HVGs: {n_hv_raw}")

    # Exclusion list — flag in .var, compute the use_for_scvi mask
    combined.var["hvg_excluded"] = build_exclusion_mask(combined.var_names, tissue).values
    combined.var["use_for_scvi"] = combined.var["highly_variable"] & ~combined.var["hvg_excluded"]
    n_excl = int((combined.var["highly_variable"] & combined.var["hvg_excluded"]).sum())
    n_use = int(combined.var["use_for_scvi"].sum())
    print(f"  Excluded from HVGs: {n_excl} (mito/ribo/hemo/sex" +
          (" + Prl/Psg/Cga" if tissue == "placenta" else "") + ")")
    print(f"  Final HVG set for scVI: {n_use}")

    # --- 4. Write + plots + summary ---
    print(f"\n[4/4] Writing outputs...")
    out_path = out_dir / "all_samples.h5ad"
    combined.write_h5ad(out_path)
    print(f"  Wrote {out_path}  ({combined.n_obs:,} cells × {combined.n_vars:,} genes)")

    plot_cells_per_sample(combined, plot_dir / "cells_per_sample.png")
    plot_hvg_dispersion(combined, plot_dir / "hvg_dispersion.png")
    plot_excluded_summary(combined, plot_dir / "hvg_exclusion_summary.png")

    # Summary table: per-sample cell counts in concatenated object
    summary = combined.obs.groupby("sample_id").size().reset_index(name="n_cells")
    summary_path = results_dir / "tables" / "summary_integration_prep.csv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, index=False)

    print(f"\n  Summary:")
    print(summary.to_string(index=False))
    print(f"\n  Output h5ad:  {out_path}")
    print(f"  Plots:        {plot_dir}")
    print(f"  Summary CSV:  {summary_path}")
    print(f"\n✓ Phase 4 complete. Ready for scVI (Phase 5) on the workstation.\n")


if __name__ == "__main__":
    main()
