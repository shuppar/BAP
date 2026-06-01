#!/usr/bin/env python
"""
04_integration_prep.py — Phase 4: concatenate + normalize + HVG + cell cycle scoring.

Changes from original:
  - Cell cycle scoring added (Step 3): computes S_score, G2M_score, phase, and
    cc_difference per cell using the Tirosh/Regev mouse cell cycle gene list.
    Scores are stored in adata.obs and carried through to all downstream phases.
    They are NOT used to regress anything here — see Phase 5 for the integration
    decision, and Phase 8 for proliferation analysis.

Layer/slot conventions for the output AnnData:
  - .X                         : raw counts (sparse, int) — what scVI consumes
  - .layers["lognorm"]         : log1p(normalize_total(X, 1e4)) — for scoring + plots
  - .layers["counts"]          : explicit alias of raw counts
  - .var["highly_variable"]    : bool, from sc.pp.highly_variable_genes
  - .var["hvg_excluded"]       : bool, True if gene is in exclusion list
  - .var["use_for_scvi"]       : highly_variable AND NOT hvg_excluded
  - .obs["S_score"]            : cell cycle S-phase score (Tirosh et al. converted to mouse)
  - .obs["G2M_score"]          : cell cycle G2/M-phase score
  - .obs["phase"]              : inferred phase: "G1", "S", or "G2M"
  - .obs["cc_difference"]      : S_score - G2M_score (Di Bella CC.Difference analog)

Usage:
  uv run python scripts/04_integration_prep.py --config config/dev.yaml
  uv run python scripts/04_integration_prep.py --config config/brain.yaml
  uv run python scripts/04_integration_prep.py --config config/placenta.yaml

Inputs:
  Per-sample h5ads in {results_dir}/h5ad/04_doublets_removed/{sample_id}.h5ad

Outputs:
  {results_dir}/h5ad/05_integration_ready/all_samples.h5ad
  {results_dir}/plots/04_integration_prep/
    - cells_per_sample.png
    - hvg_dispersion.png
    - hvg_exclusion_summary.png
    - cell_cycle_scores.png        : S vs G2M score per cell, colored by phase
    - cell_cycle_phase_bar.png     : fraction of cells per phase per sample
  {results_dir}/tables/summary_integration_prep.csv
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

from _utils import load_config, add_lognorm, phase_paths


# ----------------------------------------------------------------------------
# HVG exclusion lists
# ----------------------------------------------------------------------------

GENERIC_EXCLUDE_PREFIXES = ("mt-", "MT-", "Rps", "Rpl", "RPS", "RPL",
                             "Hbb", "Hba", "HBB", "HBA")
GENERIC_EXCLUDE_EXACT = {"Xist", "Ddx3y", "Uty", "Eif2s3y", "Kdm5d", "Tsix"}
PLACENTA_EXCLUDE_PREFIXES = ("Prl", "Psg", "Cgb")
PLACENTA_EXCLUDE_EXACT = {"Cga"}


def build_exclusion_mask(var_names: pd.Index, tissue: str) -> pd.Series:
    excl_prefixes = list(GENERIC_EXCLUDE_PREFIXES)
    excl_exact = set(GENERIC_EXCLUDE_EXACT)
    if tissue == "placenta":
        excl_prefixes += list(PLACENTA_EXCLUDE_PREFIXES)
        excl_exact |= PLACENTA_EXCLUDE_EXACT
    mask = var_names.str.startswith(tuple(excl_prefixes)) | var_names.isin(excl_exact)
    return pd.Series(mask, index=var_names, name="hvg_excluded")


# ----------------------------------------------------------------------------
# Cell cycle gene lists (mouse)
#
# Source: Tirosh et al. 2016 (Science) — originally human; converted to mouse
# title-case with manual curation. This is the de-facto standard for mouse
# single-cell cell cycle scoring (used in Seurat, scanpy, and Di Bella 2021).
#
# Production note: for maximum rigor, replace these with the gene lists from
# Di Bella et al. 2021 Supplementary Table or Kowalczyk et al. 2015 (mouse-
# specific cell cycle genes, PMID 26430063). The lists below are a well-
# validated approximation but not formally published as a standalone resource.
# ----------------------------------------------------------------------------

S_GENES_MOUSE = [
    "Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung", "Gins2",
    "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1", "Mlf1ip", "Hells", "Rfc2",
    "Rpa2", "Nasp", "Rad51ap1", "Gmnn", "Wdc", "Slbp", "Ccne2", "Ubr7",
    "Pold3", "Msh2", "Atad2", "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1",
    "Tipin", "Dscc1", "Blm", "Casp8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b",
    "Brip1", "E2f8",
]

G2M_GENES_MOUSE = [
    "Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", "Ndc80",
    "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo", "Cenpf", "Tacc3", "Fam64a",
    "Smc4", "Ccnb2", "Ckap2l", "Ckap2", "Aurkb", "Bub1", "Kif11", "Anp32e",
    "Tubb4b", "Gtse1", "Kif20b", "Hjurp", "Cdca3", "Hn1", "Cdc20", "Ttk",
    "Cdc25c", "Kif2c", "Rangap1", "Ncapd2", "Dlgap5", "Cdca2", "Cdca8",
    "Ect2", "Kif23", "Hmmr", "Aurka", "Psrc1", "Anln", "Lbr", "Ckap5",
    "Cenpe", "Ctcf", "Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa",
]


def score_cell_cycle(adata) -> None:
    """Score cell cycle phase per cell using scanpy's score_genes_cell_cycle.

    Operates on the lognorm layer (must be present in adata.layers).
    Adds to adata.obs: S_score, G2M_score, phase, cc_difference.

    cc_difference = S_score - G2M_score, following Di Bella et al. 2021.
    Cells with high cc_difference are in S phase; low (negative) are in G2M;
    near-zero are G1/quiescent.
    """
    # score_genes_cell_cycle needs lognorm in .X — use a temp copy
    tmp = adata.copy()
    tmp.X = tmp.layers["lognorm"].copy()

    s_present = [g for g in S_GENES_MOUSE if g in tmp.var_names]
    g2m_present = [g for g in G2M_GENES_MOUSE if g in tmp.var_names]
    print(f"  S-phase genes found: {len(s_present)}/{len(S_GENES_MOUSE)}")
    print(f"  G2M genes found:     {len(g2m_present)}/{len(G2M_GENES_MOUSE)}")

    if len(s_present) < 5 or len(g2m_present) < 5:
        print("  [warn] Too few cell cycle genes found — scores will be unreliable.")
        print("         Check that gene names match your annotation (mouse title-case).")

    sc.tl.score_genes_cell_cycle(tmp, s_genes=s_present, g2m_genes=g2m_present)

    adata.obs["S_score"] = tmp.obs["S_score"].values
    adata.obs["G2M_score"] = tmp.obs["G2M_score"].values
    adata.obs["phase"] = tmp.obs["phase"].values
    adata.obs["cc_difference"] = adata.obs["S_score"] - adata.obs["G2M_score"]
    del tmp


# ----------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------

def plot_hvg_dispersion(adata, out: Path) -> None:
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig(out, dpi=120, bbox_inches="tight")
    plt.close()


def plot_cells_per_sample(adata, out: Path) -> None:
    counts = adata.obs["sample_id"].value_counts().sort_index()
    n = len(counts)
    width = max(6, 0.35 * n)
    fontsize = 9 if n <= 12 else (8 if n <= 24 else 7)
    fig, ax = plt.subplots(figsize=(width, 4.5))
    ax.bar(counts.index, counts.values, color="steelblue", edgecolor="k")
    ax.set_ylabel("n cells")
    ax.set_title(f"Cells per sample after concat (n={n} samples, total: {adata.n_obs:,})")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=fontsize)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def plot_excluded_summary(adata, out: Path) -> None:
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


def plot_cell_cycle_scores(adata, out_scatter: Path, out_bar: Path) -> None:
    """Two diagnostic plots for cell cycle scoring.

    Scatter: S_score vs G2M_score colored by phase — well-separated clusters
    indicate the scoring worked. Overlapping blobs = low signal (expected for
    mostly post-mitotic tissues like 4W/3mo brain).

    Bar: fraction of cells per phase per sample — shows whether any sample is
    unusually enriched for cycling cells (could indicate a biology or QC issue).
    """
    if "phase" not in adata.obs.columns:
        return

    # Scatter
    fig, ax = plt.subplots(figsize=(6, 5))
    colors = {"S": "steelblue", "G2M": "salmon", "G1": "lightgray"}
    for phase, color in colors.items():
        sub = adata.obs[adata.obs["phase"] == phase]
        ax.scatter(sub["S_score"], sub["G2M_score"], c=color, s=3, alpha=0.5,
                   label=f"{phase} (n={len(sub):,})", rasterized=True)
    ax.set_xlabel("S score")
    ax.set_ylabel("G2M score")
    ax.set_title("Cell cycle scores (Tirosh mouse gene list)")
    ax.legend(markerscale=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_scatter, dpi=130)
    plt.close(fig)

    # Stacked bar: phase fraction per sample
    phase_counts = pd.crosstab(adata.obs["sample_id"], adata.obs["phase"],
                                normalize="index")
    # Ensure consistent column order
    for col in ["G1", "S", "G2M"]:
        if col not in phase_counts.columns:
            phase_counts[col] = 0.0
    phase_counts = phase_counts[["G1", "S", "G2M"]]

    n = len(phase_counts)
    width = max(6, 0.35 * n)
    fontsize = 9 if n <= 12 else (8 if n <= 24 else 7)
    fig, ax = plt.subplots(figsize=(width, 4.5))
    phase_counts.plot(kind="bar", stacked=True, ax=ax, width=0.8,
                      color=["lightgray", "steelblue", "salmon"], edgecolor="none")
    ax.set_ylabel("fraction of cells")
    ax.set_title("Cell cycle phase per sample")
    ax.legend(title="phase", bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=fontsize)
    fig.tight_layout()
    fig.savefig(out_bar, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Phase 4: concat + normalize + HVG + cell cycle")
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()

    print(f"\n=== Phase 4: Integration prep ===")
    print(f"Config: {args.config}")

    cfg = load_config(args.config)
    tissue = cfg["tissue"]
    in_dir = phase_paths(cfg, "doublets")["h5ad"]
    paths = phase_paths(cfg, "integration_prep")
    out_dir = paths["h5ad"]
    plot_dir = paths["plots"]

    integ_cfg = cfg.get("integration", {})
    default_n_hvg = 2000 if tissue == "placenta" else 3000
    n_hvg = int(integ_cfg.get("n_hvg", default_n_hvg))
    print(f"Tissue: {tissue}  |  target n_hvg: {n_hvg}")

    # --- 1. Load + concat ---
    print(f"\n[1/5] Loading + concatenating {len(cfg['samples'])} samples...")
    adatas = {}
    for s in cfg["samples"]:
        path = in_dir / f"{s['id']}.h5ad"
        if not path.is_file():
            sys.exit(f"ERROR: missing input {path}. Run 03_doublets.py first.")
        adatas[s["id"]] = sc.read_h5ad(path)

    combined = ad.concat(
        adatas, axis=0, join="outer", merge="same",
        label="sample_id_concat", index_unique="-",
    )
    combined.obs.drop(columns=["sample_id_concat"], errors="ignore", inplace=True)
    print(f"  Combined: {combined.n_obs:,} cells × {combined.n_vars:,} genes")

    if not sp.issparse(combined.X):
        combined.X = sp.csr_matrix(combined.X)
    if combined.X.dtype.kind == "f":
        combined.X = combined.X.astype(np.int32)

    # --- 2. Normalize ---
    print(f"\n[2/5] Log-normalizing → .layers['lognorm'] (raw counts stay in .X)")
    combined.layers["counts"] = combined.X.copy()
    add_lognorm(combined)

    # --- 3. Cell cycle scoring ---
    # Runs on lognorm layer. Scores stored in .obs; not used to modify .X or
    # the latent here. Phase 5 does NOT condition on these by default — see
    # Phase 5 docstring. Phase 8 uses them for proliferation analysis.
    print(f"\n[3/5] Cell cycle scoring (Tirosh mouse gene list)...")
    score_cell_cycle(combined)
    phase_counts = combined.obs["phase"].value_counts()
    print(f"  Phase distribution: {phase_counts.to_dict()}")

    # --- 4. HVG selection ---
    print(f"\n[4/5] Selecting {n_hvg} HVGs (seurat_v3, batch_key=pool)...")
    sc.pp.highly_variable_genes(
        combined, n_top_genes=n_hvg, flavor="seurat_v3",
        batch_key="pool", layer="counts",
    )
    n_hv_raw = int(combined.var["highly_variable"].sum())
    print(f"  Raw HVGs: {n_hv_raw}")

    combined.var["hvg_excluded"] = build_exclusion_mask(combined.var_names, tissue).values
    combined.var["use_for_scvi"] = combined.var["highly_variable"] & ~combined.var["hvg_excluded"]
    n_excl = int((combined.var["highly_variable"] & combined.var["hvg_excluded"]).sum())
    n_use = int(combined.var["use_for_scvi"].sum())
    print(f"  Excluded from HVGs: {n_excl} (mito/ribo/hemo/sex" +
          (" + Prl/Psg/Cga" if tissue == "placenta" else "") + ")")
    print(f"  Final HVG set for scVI: {n_use}")

    # --- 5. Write + plots ---
    print(f"\n[5/5] Writing outputs...")
    out_path = out_dir / "all_samples.h5ad"
    combined.write_h5ad(out_path)
    print(f"  Wrote {out_path}  ({combined.n_obs:,} cells × {combined.n_vars:,} genes)")

    plot_cells_per_sample(combined, plot_dir / "cells_per_sample.png")
    plot_hvg_dispersion(combined, plot_dir / "hvg_dispersion.png")
    plot_excluded_summary(combined, plot_dir / "hvg_exclusion_summary.png")
    plot_cell_cycle_scores(combined,
                           plot_dir / "cell_cycle_scores.png",
                           plot_dir / "cell_cycle_phase_bar.png")

    summary = combined.obs.groupby("sample_id").size().reset_index(name="n_cells")
    summary_path = paths["tables"] / "summary_integration_prep.csv"
    summary.to_csv(summary_path, index=False)

    print(f"\n  Summary:")
    print(summary.to_string(index=False))
    print(f"\n  Output h5ad : {out_path}")
    print(f"  Plots       : {plot_dir}")
    print(f"\n  Cell cycle obs columns added: S_score, G2M_score, phase, cc_difference")
    print(f"  These carry through to all downstream phases.")
    print(f"  To condition scVI on cell cycle: set scvi.condition_cell_cycle: true in YAML.")
    print(f"\n✓ Phase 4 complete. Ready for Phase 5.\n")


if __name__ == "__main__":
    main()
