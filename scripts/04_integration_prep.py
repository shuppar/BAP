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
# Ensembl -> symbol var_names swap (post-concat)
#
# SoupX (Phase 1) writes Ensembl IDs as var_names. The rest of the pipeline
# (cell cycle scoring, HVG exclusion, marker dotplots, DEG/volcano labels)
# was written assuming gene SYMBOLS as var_names, matching scanpy's default
# 10x reader and field convention. We swap once, post-concat, on the combined
# object so concat itself still aligns genes on the stable Ensembl key.
#
# Guards:
#   - blank/NA symbols fall back to their Ensembl ID (never become "")
#   - var_names_make_unique() suffixes the remaining duplicate symbols
#   - focal gene lists (cell cycle + brain marker-gate genes) are asserted
#     present and UNSUFFIXED -- if a gene we rely on got a "-1" suffix from a
#     collision, fail loud rather than silently miss it downstream.
# ----------------------------------------------------------------------------

def swap_var_names_to_symbols(adata) -> None:
    """Swap var_names from Ensembl IDs to gene symbols, IN PLACE.

    Keeps Ensembl IDs in var['gene_ids']. No-op if var_names already symbols.
    """
    if not str(adata.var_names[0]).startswith("ENSMUS"):
        print("  var_names already symbols (no swap needed)")
        return

    # Find the symbol column (SoupX writes 'symbol'; 10x reader uses other names)
    sym_col = None
    for cand in ("symbol", "Symbol", "gene_symbols", "feature_name"):
        if cand in adata.var.columns:
            sym_col = cand
            break
    if sym_col is None:
        raise ValueError(
            "var_names are Ensembl IDs but no symbol column found in var. "
            f"Columns: {list(adata.var.columns)}"
        )

    ensembl = adata.var_names.tolist()
    symbols = adata.var[sym_col].astype(str).values.copy()

    # Blank/NA symbols -> fall back to Ensembl ID
    n_blank = 0
    for i, s in enumerate(symbols):
        if s in ("", "nan", "None") or pd.isna(s):
            symbols[i] = ensembl[i]
            n_blank += 1

    adata.var["gene_ids"] = ensembl
    adata.var_names = pd.Index(symbols)
    n_dup_before = int(adata.var_names.duplicated().sum())
    adata.var_names_make_unique()

    print(f"  swapped var_names: Ensembl → symbol "
          f"({n_blank} blank→Ensembl, {n_dup_before} duplicate symbols suffixed)")

    # Guard: focal genes must be present and unsuffixed
    present = set(adata.var_names)
    missing = sorted(g for g in _FOCAL_SYMBOLS if g not in present)
    suffixed = sorted(
        g for g in _FOCAL_SYMBOLS
        if g not in present and any(v.startswith(g + "-") for v in adata.var_names)
    )
    if suffixed:
        raise ValueError(
            f"Focal gene(s) got make_unique suffixes (symbol collision): {suffixed}.\n"
            f"  These are matched by exact symbol downstream and would be silently "
            f"missed. Investigate the collision before proceeding."
        )
    if missing:
        # missing-but-not-suffixed = simply absent from the panel; that's fine
        # for cell cycle (score_genes ignores misses) but worth a note.
        print(f"  note: {len(missing)} focal gene(s) absent from panel "
              f"(not in Flex probe set): {missing[:8]}{'...' if len(missing) > 8 else ''}")


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
    "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1", "Mlf1ip", "Cenpu", "Hells", "Rfc2",
    "Rpa2", "Nasp", "Rad51ap1", "Gmnn", "Wdc", "Slbp", "Ccne2", "Ubr7",
    "Pold3", "Msh2", "Atad2", "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1",
    "Tipin", "Dscc1", "Blm", "Casp8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b",
    "Brip1", "E2f8",
]   # Mlf1ip/Cenpu: same gene, old/new MGI symbol — both kept so the list
    # survives either annotation version (score_genes_cell_cycle ignores misses).

G2M_GENES_MOUSE = [
    "Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", "Ndc80",
    "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo", "Cenpf", "Tacc3",
    "Fam64a", "Pimreg", "Smc4", "Ccnb2", "Ckap2l", "Ckap2", "Aurkb", "Bub1",
    "Kif11", "Anp32e", "Tubb4b", "Gtse1", "Kif20b", "Hjurp", "Cdca3",
    "Hn1", "Jpt1", "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1", "Ncapd2",
    "Dlgap5", "Cdca2", "Cdca8", "Ect2", "Kif23", "Hmmr", "Aurka", "Psrc1",
    "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf", "Nek2", "G2e3", "Gas2l3",
    "Cbx5", "Cenpa",
]   # Fam64a/Pimreg and Hn1/Jpt1: old/new MGI symbols for the same genes —
    # both kept so the list survives either annotation version.


# Genes the pipeline matches by exact symbol downstream. If any of these got
# suffixed by make_unique (i.e. collided with another symbol), we want to know.
_FOCAL_SYMBOLS = set(S_GENES_MOUSE) | set(G2M_GENES_MOUSE) | {
    # brain marker-gate genes (Phase 7 BRAIN_GATE_CONFIG)
    "Cx3cr1", "P2ry12", "Tmem119", "Csf1r", "Aif1",
    "Aqp4", "Gja1", "Slc1a3", "Aldh1l1",
    "Mbp", "Mog", "Plp1", "Mag",
    "Cldn5", "Pecam1", "Cdh5",
}


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

    # HARD FAIL, not a warning: if too few genes match, the scores would be
    # garbage but would still look like valid numbers downstream. The most
    # likely cause is a gene-naming mismatch (e.g. var_names are Ensembl IDs,
    # not mouse symbols). Stop here rather than produce plausible-but-wrong
    # cell cycle assignments that silently corrupt Phase 8 proliferation analysis.
    MIN_GENES = 10
    if len(s_present) < MIN_GENES or len(g2m_present) < MIN_GENES:
        example_vars = list(tmp.var_names[:5])
        raise ValueError(
            f"Too few cell cycle genes matched (S={len(s_present)}, "
            f"G2M={len(g2m_present)}; need >={MIN_GENES} each).\n"
            f"  adata.var_names look like: {example_vars}\n"
            f"  Expected mouse gene symbols (e.g. 'Mcm5', 'Top2a').\n"
            f"  If var_names are Ensembl IDs, map them to symbols before this phase,\n"
            f"  or pass symbol-based gene lists. Refusing to compute unreliable\n"
            f"  cell cycle scores that would silently corrupt downstream analysis."
        )

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

    # Capture Ensembl->symbol mapping from the union of all samples BEFORE concat.
    # merge="same" can drop the var['symbol'] column if per-sample var frames
    # differ at all, so we rebuild it on the combined object from this mapping.
    ens2sym: dict[str, str] = {}
    for a in adatas.values():
        for cand in ("symbol", "Symbol", "gene_symbols", "feature_name"):
            if cand in a.var.columns:
                for ens, sym in zip(a.var_names, a.var[cand].astype(str).values):
                    if ens not in ens2sym and sym not in ("", "nan", "None"):
                        ens2sym[ens] = sym
                break

    combined = ad.concat(
        adatas, axis=0, join="outer", merge="same",
        label="sample_id_concat", index_unique="-",
    )
    combined.obs.drop(columns=["sample_id_concat"], errors="ignore", inplace=True)
    print(f"  Combined: {combined.n_obs:,} cells × {combined.n_vars:,} genes")

    # Ensure a 'symbol' column exists on the combined object (rebuild from the
    # pre-concat mapping if merge="same" dropped it).
    if "symbol" not in combined.var.columns and ens2sym:
        combined.var["symbol"] = [
            ens2sym.get(ens, ens) for ens in combined.var_names
        ]

    # Swap Ensembl IDs -> symbols (concat aligned on stable Ensembl key above;
    # everything downstream wants symbols). Ensembl kept in var['gene_ids'].
    swap_var_names_to_symbols(combined)

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
    summary_path = paths["tables"] / "04_integration_prep_summary.csv"
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
