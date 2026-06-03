#!/usr/bin/env python
"""
08d_trajectory.py — Phase 8d: trajectory analysis.

Two components:

  1. PAGA (always runs)
     Partition-based graph abstraction on the scVI latent. Shows which cell
     types are transcriptionally adjacent and how that connectivity differs
     between Early Stress, Late Stress, and Relaxed groups per age.
     No directionality assumptions — robust and reviewer-friendly.

  2. Diffusion pseudotime (DPT)
     Orders cells along a transcriptional continuum anchored at a progenitor
     root. The study question is whether prenatal stress shifts cells along a
     lineage axis — valid at every age, so all ages are treated identically.
     The per-lineage group comparison (Kruskal-Wallis + η²) runs both pooled
     across ages and split per age, so no age is dropped. Age-split rows carry
     a pool_age_confound caveat (project doc §2). Small-n rows carry a note to
     read η² over the p-value.

What is NOT here and why:
  - RNA velocity: fundamentally incompatible with 10x Flex. Probes target
    exons only; intronic/unspliced reads are not captured. 10x Genomics
    explicitly states velocity is not recommended for Flex data
    (kb.10xgenomics.com article 25938615598477).
  - CellRank: its primary value is combining velocity directionality with
    graph connectivity. Without velocity, ConnectivityKernel-only CellRank
    duplicates what PAGA already provides (transcriptional adjacency +
    terminal state identification) with more complexity and a sidecar venv.
    Not worth it for this dataset.

Design notes:
  - Pool-age confounding (project doc §2): P1/4W/3mo are dominated by
    different pools. scVI corrects for pool in the latent, but any cross-age
    developmental trend must note this caveat in methods.
  - Focal lineages: oligodendrocyte (OPC→MOL), microglia, astrocyte
    maturation. Auto-detected from cell type labels.

Usage:
  uv run python scripts/08d_trajectory.py --config config/dev.yaml
  uv run python scripts/08d_trajectory.py --config config/brain.yaml
  uv run python scripts/08d_trajectory.py --config config/brain.yaml \\
      --root-celltype "Radial glia / NPCs"
  uv run python scripts/08d_trajectory.py --config config/placenta.yaml \\
      --root-celltype "LaTP"

Inputs (first that exists):
  {results_dir}/h5ad/08b_label_transferred/all_samples.h5ad  (Phase 7c)
  {results_dir}/h5ad/08_annotated/all_samples.h5ad           (Phase 7)

Outputs:
  {results_dir}/plots/08d_trajectory/
    paga/
      paga_by_celltype.png
      paga_by_group_{age}.png              : one per age
      paga_transitions_heatmap.png
      umap_paga_init.png
    pseudotime/
      dpt_umap.png                         : pseudotime on UMAP, all cells
      dpt_violin_by_celltype.png           : all ages
      dpt_violin_by_group_{lineage}.png    : per focal lineage, ages pooled
      dpt_group_comparison.png             : Kruskal-Wallis + η², all lineage×age rows
  {results_dir}/tables/
    trajectory_paga_connectivities.csv
    trajectory_dpt_group_comparison.csv    : group_level col = all_ages | age=<X>;
                                             note col carries pool_age + small_n caveats
"""

import argparse
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scipy.stats as stats

from _utils import load_config, add_lognorm, phase_table_dir

# RNA velocity is NOT used in this script.
# 10x Flex probes target exons only — intronic/unspliced reads are not captured
# and velocity is explicitly not recommended for Flex data by 10x Genomics.
# Reference: kb.10xgenomics.com article 25938615598477


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

LABEL_KEY_PRIORITY = [
    "manual_annotation", "scanvi_celltype", "celltypist_majority", "provisional_celltype",
]

# Default progenitor cell types (root anchoring for DPT).
# Ordered by preference — first match in adata wins.
DEFAULT_ROOT_CELLTYPES = {
    "brain": [
        "Radial glia / NPCs", "Radial_Glia", "NPC", "OPC",
        "Radial glia", "NPCs", "Progenitor",
    ],
    "placenta": [
        "LaTP", "Labyrinth_Trophoblast_Progenitor", "Trophoblast progenitor",
        "Stem trophoblast",
    ],
}

FOCAL_LINEAGES = {
    "brain": {
        "Oligodendrocyte lineage": ["OPC", "COP", "NFOL", "MFOL", "MOL",
                                    "Oligodendrocyte", "oligodendro"],
        "Microglia":               ["Microglia", "microglia", "DAM", "Homeostatic"],
        "Astrocyte maturation":    ["Astrocyte", "astrocyte", "Radial glia", "NPC"],
    },
    "placenta": {
        "Trophoblast differentiation": ["LaTP", "SynT", "S-TGC", "Trophoblast",
                                        "trophoblast", "SpT", "GlyT"],
    },
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def resolve_celltype_key(adata, explicit=None):
    if explicit:
        if explicit not in adata.obs.columns:
            sys.exit(f"ERROR: --celltype-key '{explicit}' not in adata.obs.")
        return explicit
    for key in LABEL_KEY_PRIORITY:
        if key in adata.obs.columns:
            if key == "manual_annotation" and adata.obs[key].astype(str).eq("").all():
                continue
            return key
    sys.exit("ERROR: no usable cell-type label column. Run Phase 7 first.")


def find_root_cell(adata, celltype_key, root_candidates):
    """Return the obs_name of one cell to use as DPT root.

    Picks the cell closest to the centroid of the first matching root cell type
    in the scVI latent (X_scVI). Falls back to the first cell of the first
    cluster if nothing matches.
    """
    labels = adata.obs[celltype_key].astype(str)
    for cand in root_candidates:
        mask = labels.str.contains(cand, case=False, regex=False)
        if mask.sum() > 0:
            sub = adata[mask]
            centroid = sub.obsm["X_scVI"].mean(axis=0)
            dists = np.linalg.norm(sub.obsm["X_scVI"] - centroid, axis=1)
            root_idx = np.argmin(dists)
            root_name = sub.obs_names[root_idx]
            print(f"  Root cell type matched: '{cand}' → cell {root_name}")
            return root_name, cand
    # Fallback
    root_name = adata.obs_names[0]
    print(f"  [warn] No root cell type matched. Falling back to first cell: {root_name}")
    print(f"         Pass --root-celltype to set a biologically meaningful root.")
    return root_name, "fallback"


def cells_in_lineage(adata, celltype_key, fragments):
    """Boolean mask for cells whose label contains any fragment (case-insensitive)."""
    labels = adata.obs[celltype_key].astype(str)
    mask = pd.Series(False, index=adata.obs_names)
    for frag in fragments:
        mask = mask | labels.str.contains(frag, case=False, regex=False)
    return mask.values


def safe_savefig(fig, path, dpi=140):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 1. PAGA
# ---------------------------------------------------------------------------

def run_paga(adata, celltype_key, plot_dir, table_dir, seed):
    """PAGA on scVI neighbors, grouped by cell type.

    Produces: paga_by_celltype.png, paga_by_group_{age}.png,
              paga_transitions_heatmap.png, trajectory_paga_connectivities.csv
    """
    print("\n  [PAGA] Computing...")
    paga_dir = plot_dir / "paga"
    paga_dir.mkdir(parents=True, exist_ok=True)

    adata.obs[celltype_key] = adata.obs[celltype_key].astype("category")
    sc.tl.paga(adata, groups=celltype_key)

    # Save connectivities
    conn = pd.DataFrame(
        adata.uns["paga"]["connectivities"].toarray()
        if sp.issparse(adata.uns["paga"]["connectivities"])
        else adata.uns["paga"]["connectivities"],
        index=adata.obs[celltype_key].cat.categories,
        columns=adata.obs[celltype_key].cat.categories,
    )
    conn.to_csv(table_dir / "08d_trajectory_paga_connectivities.csv")

    # PAGA coloured by cell type
    fig, ax = plt.subplots(figsize=(8, 7))
    sc.pl.paga(adata, color=celltype_key, ax=ax, show=False,
               title="PAGA — cell type connectivity",
               node_size_scale=1.5, edge_width_scale=1.0,
               fontsize=7, frameon=False)
    fig.tight_layout()
    safe_savefig(fig, paga_dir / "paga_by_celltype.png")

    # PAGA coloured by group, per age
    ages = adata.obs["age"].unique() if "age" in adata.obs.columns else ["all"]
    for age in ages:
        sub = adata[adata.obs["age"] == age] if age != "all" else adata
        if sub.n_obs < 20:
            continue
        sub_copy = sub.copy()
        sub_copy.obs[celltype_key] = sub_copy.obs[celltype_key].astype("category")
        try:
            sc.tl.paga(sub_copy, groups=celltype_key)
            for color_key in ("group", "pool"):
                if color_key not in sub_copy.obs.columns:
                    continue
                fig, ax = plt.subplots(figsize=(7, 6))
                sc.pl.paga(sub_copy, color=color_key, ax=ax, show=False,
                           title=f"PAGA age={age}, coloured by {color_key}",
                           node_size_scale=1.2, frameon=False, fontsize=7)
                fig.tight_layout()
                safe_savefig(fig, paga_dir / f"paga_by_{color_key}_age{age}.png")
        except Exception as e:
            print(f"    [warn] PAGA for age={age} failed: {e}")

    # Connectivity heatmap
    fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(conn)), max(5, 0.4 * len(conn))))
    import seaborn as sns
    sns.heatmap(conn, ax=ax, cmap="Blues", vmin=0, vmax=1,
                linewidths=0.3, annot=len(conn) <= 15, fmt=".2f",
                annot_kws={"size": 6})
    ax.set_title("PAGA connectivities (cell type)")
    fig.tight_layout()
    safe_savefig(fig, paga_dir / "paga_transitions_heatmap.png")

    # UMAP coloured by PAGA pos (initialise UMAP from PAGA)
    sc.pl.paga(adata, show=False)  # needed to init adata.uns['paga']['pos']
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.umap(adata, init_pos="paga", random_state=seed)

    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(adata, color=celltype_key, ax=ax, show=False, frameon=False,
               legend_fontsize=6, size=8, title="UMAP (PAGA-initialised)")
    fig.tight_layout()
    safe_savefig(fig, paga_dir / "umap_paga_init.png")

    print(f"    PAGA done. Connectivities → {table_dir}/trajectory_paga_connectivities.csv")

    # Edge diagnostics — audit suspicious edges offline (no workstation needed)
    write_paga_edge_diagnostics(adata, celltype_key, conn, table_dir)

    return adata   # UMAP coords updated


# Ambient / contamination markers whose presence on an edge suggests the
# connection is driven by soup, not biology (project doc: snRNA ambient is severe).
AMBIENT_MARKERS = {
    "brain":    ["Malat1", "Meg3", "mt-Co1", "mt-Co3", "mt-Atp6", "mt-Nd1",
                 "Hbb-bs", "Hba-a1", "Hbb-bt", "Hba-a2"],
    "placenta": ["Hbb-bs", "Hba-a1", "Hbb-bt", "Hba-a2", "Prl3b1", "Prl8a8",
                 "Psg17", "Psg18", "mt-Co1", "Malat1"],
}


def write_paga_edge_diagnostics(adata, celltype_key, conn, table_dir,
                                 min_connectivity=0.1, n_shared=15):
    """One row per cell-type pair edge above min_connectivity, with the info
    needed to judge offline whether the edge is biology or artifact:

      - connectivity                  : PAGA edge weight
      - n_cells_A / n_cells_B         : cluster sizes (tiny clusters = noisy edges)
      - doublet_rate_A / _B           : high doublet rate → edge may be doublets
      - top_shared_genes              : genes highly expressed in BOTH endpoints
      - ambient_driven                : True if shared genes are mostly ambient markers
      - n_ambient_in_shared           : how many shared genes are ambient
      - mean_pct_mt_A / _B            : high %mt → contamination

    Audit recipe (offline): for a surprising edge, check ambient_driven first,
    then doublet_rate, then whether top_shared_genes are real lineage markers.
    """
    tissue = adata.uns.get("tissue", "brain")
    ambient = set(AMBIENT_MARKERS.get(tissue, []))

    # Map var_names → symbols if available (edges report human-readable genes)
    symbol_map = None
    for c in ("symbol", "gene_symbol", "Symbol"):
        if c in adata.var.columns:
            symbol_map = dict(zip(adata.var_names.astype(str), adata.var[c].astype(str)))
            break

    cats = list(conn.index)
    labels = adata.obs[celltype_key].astype(str)

    # Per-cluster mean expression (lognorm) for shared-gene detection
    if "lognorm" not in adata.layers:
        add_lognorm(adata)
    import scipy.sparse as _sp
    L = adata.layers["lognorm"]
    L = L.toarray() if _sp.issparse(L) else np.asarray(L)

    # Precompute per-cluster mean expression + QC summaries
    cl_mean, cl_meta = {}, {}
    for c in cats:
        m = (labels == c).values
        if m.sum() == 0:
            continue
        cl_mean[c] = L[m].mean(axis=0)
        meta = {"n_cells": int(m.sum())}
        if "doublet_class" in adata.obs.columns:
            meta["doublet_rate"] = float(
                (adata.obs.loc[m, "doublet_class"].astype(str) == "doublet").mean())
        else:
            meta["doublet_rate"] = np.nan
        meta["mean_pct_mt"] = (float(adata.obs.loc[m, "pct_counts_mt"].mean())
                               if "pct_counts_mt" in adata.obs.columns else np.nan)
        cl_meta[c] = meta

    var_names = adata.var_names.astype(str).values
    rows = []
    for i, a in enumerate(cats):
        for j, b in enumerate(cats):
            if j <= i:
                continue
            w = float(conn.iloc[i, j])
            if w < min_connectivity or a not in cl_mean or b not in cl_mean:
                continue
            # Shared signal: genes high in both endpoints (min of the two means)
            shared_score = np.minimum(cl_mean[a], cl_mean[b])
            top_idx = np.argsort(shared_score)[::-1][:n_shared]
            top_genes_raw = var_names[top_idx]
            top_genes = [symbol_map.get(g, g) if symbol_map else g for g in top_genes_raw]
            n_amb = sum(1 for g in top_genes if g in ambient)
            rows.append({
                "celltype_A": a, "celltype_B": b,
                "connectivity": round(w, 4),
                "n_cells_A": cl_meta[a]["n_cells"], "n_cells_B": cl_meta[b]["n_cells"],
                "doublet_rate_A": round(cl_meta[a]["doublet_rate"], 4)
                                  if not np.isnan(cl_meta[a]["doublet_rate"]) else np.nan,
                "doublet_rate_B": round(cl_meta[b]["doublet_rate"], 4)
                                  if not np.isnan(cl_meta[b]["doublet_rate"]) else np.nan,
                "mean_pct_mt_A": round(cl_meta[a]["mean_pct_mt"], 3)
                                 if not np.isnan(cl_meta[a]["mean_pct_mt"]) else np.nan,
                "mean_pct_mt_B": round(cl_meta[b]["mean_pct_mt"], 3)
                                 if not np.isnan(cl_meta[b]["mean_pct_mt"]) else np.nan,
                "n_ambient_in_shared": n_amb,
                "ambient_driven": n_amb >= (n_shared // 3),   # ≥1/3 shared are ambient
                "top_shared_genes": ", ".join(top_genes),
            })

    if rows:
        df = (pd.DataFrame(rows).sort_values("connectivity", ascending=False)
              .reset_index(drop=True))
        out = table_dir / "08d_trajectory_paga_edge_diagnostics.csv"
        df.to_csv(out, index=False)
        n_amb = int(df["ambient_driven"].sum())
        print(f"    Edge diagnostics → {out} ({len(df)} edges, {n_amb} ambient-driven)")
    else:
        print(f"    [info] No edges above connectivity {min_connectivity} to diagnose.")


# ---------------------------------------------------------------------------
# 2. Diffusion pseudotime
# ---------------------------------------------------------------------------

def run_dpt(adata, celltype_key, root_name, plot_dir, table_dir):
    """Diffusion map + DPT anchored at root_name.

    DPT orders cells along a transcriptional continuum from the progenitor root.
    The question for this study is whether prenatal stress alters that ordering
    within a cell-type lineage — a valid question at every age (whether the
    tissue is actively differentiating or in steady state, the lineage axis
    still exists and stress can shift cells along it).

    All ages are treated identically. The per-lineage group comparison
    (Kruskal-Wallis + η²) is run twice:
      - pooled across all ages (group_level = "all_ages")
      - split by age (group_level = each age value)
    so every age that has the lineage gets its own row. No age is excluded.

    Caveat tags carried in the output (informative, not gating):
      - small_n note when total cells < 50
      - pool_age_confound note (project doc §2): age co-varies with pool;
        any age-split result inherits that confound.

    Produces: dpt_umap.png, dpt_violin_by_celltype.png,
              dpt_violin_by_group_{lineage}.png (pooled across ages),
              dpt_group_comparison.png, trajectory_dpt_group_comparison.csv
    """
    print("\n  [DPT] Computing diffusion map + pseudotime...")
    dpt_dir = plot_dir / "pseudotime"
    dpt_dir.mkdir(parents=True, exist_ok=True)

    tissue = adata.uns.get("tissue", "brain")
    ages_in_data = (sorted(adata.obs["age"].unique())
                    if "age" in adata.obs.columns else [])
    print(f"    Ages present (all treated identically): {ages_in_data or '(no age column)'}")

    # Diffusion map on scVI neighbors
    sc.tl.diffmap(adata)

    # Set root
    root_idx = adata.obs_names.get_loc(root_name)
    adata.uns["iroot"] = root_idx
    sc.tl.dpt(adata)

    # UMAP coloured by DPT
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    sc.pl.umap(adata, color="dpt_pseudotime", ax=axes[0], show=False,
               frameon=False, color_map="viridis", size=8,
               title="Diffusion pseudotime")
    sc.pl.umap(adata, color=celltype_key, ax=axes[1], show=False,
               frameon=False, legend_fontsize=6, size=8, title="Cell types")
    fig.tight_layout()
    safe_savefig(fig, dpt_dir / "dpt_umap.png")

    # Violin: DPT per cell type (all cells)
    order = (adata.obs.groupby(celltype_key, observed=True)["dpt_pseudotime"]
             .median().sort_values().index.tolist())
    adata.obs[celltype_key] = pd.Categorical(adata.obs[celltype_key], categories=order)
    fig, ax = plt.subplots(figsize=(max(8, 0.6 * adata.obs[celltype_key].nunique()), 5))
    sc.pl.violin(adata, keys="dpt_pseudotime", groupby=celltype_key,
                 ax=ax, show=False, rotation=45)
    ax.set_title("Pseudotime distribution per cell type (all ages)")
    fig.tight_layout()
    safe_savefig(fig, dpt_dir / "dpt_violin_by_celltype.png")

    # Per focal lineage: Kruskal-Wallis across groups, both pooled and per-age.
    lineages = FOCAL_LINEAGES.get(tissue, {})
    rows = []

    def kw_for(sub, lineage_name, group_level, age_confound):
        """Run KW + η² on one (lineage, age-slice) subset, append a row."""
        if "group" not in sub.obs.columns or sub.obs["group"].nunique() < 2:
            return
        groups = sorted(sub.obs["group"].unique())
        data = [sub.obs.loc[sub.obs["group"] == g, "dpt_pseudotime"].dropna().values
                for g in groups]
        if len(data) < 2 or not all(len(d) > 0 for d in data):
            return
        try:
            stat, pval = stats.kruskal(*data)
            k = len(groups)
            n = sum(len(d) for d in data)
            eta2 = max(0.0, (stat - k + 1) / (n - k)) if n > k else np.nan
        except Exception as e:
            print(f"    [warn] KW {lineage_name} / {group_level}: {e}")
            return
        notes = []
        if n < 50:
            notes.append("small_n — report η², not just p-value")
        if age_confound:
            notes.append("pool_age_confound (§2): age co-varies with pool")
        rows.append({
            "lineage": lineage_name,
            "group_level": group_level,
            "n_cells": int(sum(len(d) for d in data)),
            "groups": str(groups),
            "kruskal_stat": round(stat, 4),
            "kruskal_pval": round(pval, 6),
            "eta2_effect_size": round(eta2, 4),
            "note": "; ".join(notes),
        })

    for lineage_name, fragments in lineages.items():
        mask = cells_in_lineage(adata, celltype_key, fragments)
        if mask.sum() < 20:
            continue
        lin = adata[mask]

        # Pooled across all ages
        kw_for(lin, lineage_name, "all_ages", age_confound=False)

        # Per age (no age dropped) — carries pool-age confound caveat
        if "age" in lin.obs.columns:
            for age in sorted(lin.obs["age"].unique()):
                age_sub = lin[lin.obs["age"] == age]
                if age_sub.n_obs < 20:
                    continue
                kw_for(age_sub, lineage_name, f"age={age}", age_confound=True)

        # Violin: pooled across ages, by group
        if "group" in lin.obs.columns and lin.obs["group"].nunique() >= 2:
            groups = sorted(lin.obs["group"].unique())
            data = [lin.obs.loc[lin.obs["group"] == g, "dpt_pseudotime"].dropna().values
                    for g in groups]
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.violinplot(data, showmedians=True)
            ax.set_xticks(range(1, len(groups) + 1))
            ax.set_xticklabels(groups, rotation=30, ha="right")
            ax.set_ylabel("DPT pseudotime")
            ax.set_title(f"{lineage_name} — pseudotime by group (all ages pooled)")
            fig.tight_layout()
            slug = lineage_name.replace(" ", "_").replace("/", "_")
            safe_savefig(fig, dpt_dir / f"dpt_violin_by_group_{slug}.png")

    if rows:
        kw_df = pd.DataFrame(rows)
        kw_df.to_csv(table_dir / "08d_trajectory_dpt_group_comparison.csv", index=False)
        print(f"    DPT group comparison: {len(rows)} rows "
              f"({kw_df['lineage'].nunique()} lineages × age-levels)")
        print(kw_df[["lineage", "group_level", "kruskal_pval",
                      "eta2_effect_size"]].to_string(index=False))

        # Summary bar: all rows, labelled lineage + group_level
        fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(kw_df)), max(4, 0.4 * len(kw_df))))
        colors = ["salmon" if p < 0.05 else "lightgray" for p in kw_df["kruskal_pval"]]
        ax.barh(kw_df["lineage"] + " (" + kw_df["group_level"] + ")",
                -np.log10(kw_df["kruskal_pval"].clip(lower=1e-10)),
                color=colors)
        ax.axvline(-np.log10(0.05), color="k", ls="--", lw=0.8)
        ax.set_xlabel("-log10(p-value, Kruskal-Wallis)")
        ax.set_title("DPT group comparison per focal lineage × age\n"
                     "salmon = p<0.05 | check η² + note column for caveats")
        fig.tight_layout()
        safe_savefig(fig, dpt_dir / "dpt_group_comparison.png")

    print(f"    DPT done. obs['dpt_pseudotime'] written for all cells.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Phase 8d: trajectory analysis")
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--celltype-key", default=None,
                    help="obs column for cell type labels (auto-detected if omitted)")
    ap.add_argument("--root-celltype", default=None,
                    help="Cell type label fragment for DPT root (e.g. 'Radial glia / NPCs')")
    args = ap.parse_args()

    print(f"\n=== Phase 8d: Trajectory analysis (PAGA + DPT) ===")
    print(f"  RNA velocity: NOT run (10x Flex incompatible — kb.10xgenomics.com/25938615598477)")
    print(f"  CellRank: NOT run (no velocity → ConnectivityKernel duplicates PAGA)")
    cfg = load_config(args.config)
    tissue = cfg["tissue"]
    seed   = int(cfg.get("random_seed", 42))

    base = Path(cfg["results_dir"]) / "h5ad"
    candidates = [base / "08b_label_transferred" / "all_samples.h5ad",
                  base / "08_annotated" / "all_samples.h5ad"]
    in_path = next((p for p in candidates if p.is_file()), None)
    if in_path is None:
        sys.exit("ERROR: no annotated input. Checked:\n  " +
                 "\n  ".join(str(p) for p in candidates))
    print(f"\n  Input:  {in_path}")
    print(f"  Tissue: {tissue}")

    plot_dir  = Path(cfg["results_dir"]) / "plots" / "08d_trajectory"
    table_dir = phase_table_dir(cfg, "08d_trajectory")
    plot_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n[1/2] Loading data...")
    adata = sc.read_h5ad(in_path)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    if "X_scVI" not in adata.obsm:
        sys.exit("ERROR: 'X_scVI' not in obsm. Run Phase 5 first.")
    if "X_umap" not in adata.obsm:
        sys.exit("ERROR: 'X_umap' not in obsm. Run Phase 5 first.")

    celltype_key = resolve_celltype_key(adata, args.celltype_key)
    print(f"  Cell type column: '{celltype_key}' ({adata.obs[celltype_key].nunique()} types)")

    adata.uns["tissue"] = tissue
    add_lognorm(adata)

    if "neighbors" not in adata.uns:
        print("  Recomputing neighbors on X_scVI...")
        sc.pp.neighbors(adata, use_rep="X_scVI", random_state=seed)

    # -----------------------------------------------------------------------
    # 1. PAGA
    # -----------------------------------------------------------------------
    print(f"\n[1/2] PAGA...")
    adata = run_paga(adata, celltype_key, plot_dir, table_dir, seed)

    # -----------------------------------------------------------------------
    # 2. Diffusion pseudotime (age-scoped)
    # -----------------------------------------------------------------------
    print(f"\n[2/2] Diffusion pseudotime...")
    root_candidates = ([args.root_celltype] if args.root_celltype
                       else DEFAULT_ROOT_CELLTYPES.get(tissue, []))
    root_name, _ = find_root_cell(adata, celltype_key, root_candidates)
    run_dpt(adata, celltype_key, root_name, plot_dir, table_dir)

    # -----------------------------------------------------------------------
    # Save
    # -----------------------------------------------------------------------
    if "lognorm" in adata.layers:
        del adata.layers["lognorm"]
    out_path = Path(cfg["results_dir"]) / "h5ad" / "08d_trajectory" / "all_samples.h5ad"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path)

    print(f"\n  Written: {out_path}")
    print(f"  Plots:   {plot_dir}")
    print(f"\n✓ Phase 8d complete.")
    print(f"\nKey outputs:")
    print(f"  paga/paga_by_celltype.png          — cell type connectivity graph")
    print(f"  paga/paga_transitions_heatmap.png  — connectivity matrix")
    print(f"  pseudotime/dpt_umap.png            — pseudotime on UMAP")
    print(f"  pseudotime/dpt_group_comparison.png — Kruskal-Wallis + η², per lineage×age")
    print(f"  tables/trajectory_dpt_group_comparison.csv")
    print(f"    → 'group_level' col: 'all_ages' (pooled) or 'age=<X>' (per age)")
    print(f"    → 'note' col carries small_n + pool_age_confound caveats")
    print(f"\nPool-age caveat (project doc §2): P1/4W/3mo dominated by different")
    print(f"  pools. Age-split DPT rows inherit this confound (flagged in 'note').\n")


if __name__ == "__main__":
    main()
