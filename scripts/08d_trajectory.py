#!/usr/bin/env python
"""
08d_trajectory.py — Phase 8d: trajectory analysis (PAGA + diffusion pseudotime).

Supplementary / mechanistic, NOT a headline. Complements the 8a composition
shifts and the 8b disruption finding by asking whether maturation *progression
itself* changes under prenatal stress.

TWO MODES
---------
A. Whole-tissue PAGA (no --subcluster):
     PAGA on the cross-age-aligned broad tier (brain: celltypist_broad;
     placenta: celltype_majority). Connectivity graph + heatmap + edge
     diagnostics. NO DPT (pseudotime across unrelated broad types is
     uninterpretable). Does NOT touch the saved X_umap.

B. Focal lineage (--subcluster <lineage>):
     Loads the 08c subcluster object (brain) or subsets the main object to a
     trophoblast cell-type list (placenta). Drops contaminants + unassigned,
     recomputes neighbors + diffmap on the lineage's own X_scVI, then per root:
       - DPT pseudotime, UMAP, marker-vs-pseudotime trends, per-celltype violin
       - per-DONOR summary (median pseudotime + mature fraction)
       - per-(sex stratum) group comparison: animal is the statistical unit
         (NOT pooled cells — that is pseudoreplication)

PER-DONOR GROUP COMPARISON (the part that earns 8d its place)
  Each donor -> a scalar (median lineage pseudotime; "mature fraction" past a
  quantile threshold). Then, on those donor scalars:
    - Mann-Whitney U for the primary pairwise contrasts
      (brain: Early-vs-Relaxed + Late-vs-Relaxed [+ Early-vs-Late secondary];
       placenta within-age: E12.5=Early-vs-Relaxed, E18.5=Late-vs-Relaxed)
    - Kruskal-Wallis omnibus where >=3 groups (brain only)
  Effect sizes reported (rank-biserial for MW-U, eta^2 for KW). With ~4
  donors/group this is underpowered: low_n flagged, read effect size + the
  per-donor scatter, not the p-value. Sex strata {combined, M, F} iterate on the
  comparison only (the pseudotime axis is computed once on all cells); M/F low_n.

AGE HANDLING
  brain: one DPT axis across all ages; comparison run per-age (groups within an
         age, pool ~constant) AND all-ages-pooled (flagged pool_age_confound,
         project doc §2 — age co-varies with pool).
  placenta: within-age by design (age = stress window, pool-confounded); DPT
         recomputed separately for E12.5 / E18.5, comparison within each age.

NOT HERE: RNA velocity (10x Flex is probe/exon-only — no spliced/unspliced;
  10x explicitly does not recommend velocity for Flex). CellRank without
  velocity only duplicates PAGA. See project doc / INSTRUCTIONS "Trajectory (8d)".

Usage:
  # whole-tissue PAGA
  uv run python scripts/08d_trajectory.py --config config/brain.yaml
  uv run python scripts/08d_trajectory.py --config config/placenta.yaml
  # focal lineage (brain subcluster object)
  uv run python scripts/08d_trajectory.py --config config/brain.yaml \\
      --subcluster opc_oligodendrocytes
  # placenta trophoblast (subset of main object, within-age)
  uv run python scripts/08d_trajectory.py --config config/placenta.yaml \\
      --subcluster trophoblast

Config block (per tissue YAML):
  trajectory:
    whole_paga_key: celltypist_broad     # placenta: celltype_majority
    mature_quantile: 0.66
    min_donors: 2                        # per group to run; <3 -> low_n
    reliable_donors: 3
    within_age: true                     # placenta only
    lineages:
      opc_oligodendrocytes:
        root: OPC                        # str or list of str (immune has two)
        order:   [OPC, COP, MFOL, MOL]
        markers: [Pdgfra, Cspg4, ...]
      ...
      trophoblast:                       # placenta: needs `celltypes`
        celltypes: [LaTP, SynTI, ...]
        root: LaTP
        ...
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

from _utils import load_config, add_lognorm, phase_table_dir, iter_strata


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Ambient / contamination markers whose dominance on a PAGA edge suggests the
# connection is soup-driven, not biology (snRNA ambient is severe — project doc).
AMBIENT_MARKERS = {
    "brain":    ["Malat1", "Meg3", "mt-Co1", "mt-Co3", "mt-Atp6", "mt-Nd1",
                 "Hbb-bs", "Hba-a1", "Hbb-bt", "Hba-a2"],
    "placenta": ["Hbb-bs", "Hba-a1", "Hbb-bt", "Hba-a2", "Prl3b1", "Prl8a8",
                 "Psg17", "Psg18", "mt-Co1", "Malat1"],
}

# Immune microglia "states" are not a clean differentiation lineage — DPT there
# is a state/activation axis, flagged so it is read accordingly.
STATE_AXIS_LINEAGES = {"immune"}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def safe_savefig(fig, path, dpi=150):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def slugify(s):
    return str(s).replace(" ", "_").replace("/", "_").replace(".", "p")


def classify_group(g):
    """Map a raw group label to canonical {Relaxed, Early, Late} or None."""
    gl = str(g).lower()
    if "relax" in gl:
        return "Relaxed"
    if "early" in gl:
        return "Early"
    if "late" in gl:
        return "Late"
    return None


def symbol_to_var(adata):
    """Return dict mapping gene symbol -> var_name, for marker lookup.

    var_names may be Ensembl IDs; var['symbol'] (if present) holds the symbol.
    Falls back to identity when var_names are already symbols.
    """
    for c in ("symbol", "gene_symbol", "Symbol"):
        if c in adata.var.columns:
            d = {}
            for vn, sym in zip(adata.var_names.astype(str), adata.var[c].astype(str)):
                if sym and sym != "nan" and sym not in d:
                    d[sym] = vn
            # also allow direct var_name hits
            for vn in adata.var_names.astype(str):
                d.setdefault(vn, vn)
            return d
    return {vn: vn for vn in adata.var_names.astype(str)}


def drop_noncelltypes(adata, label_col):
    """Drop Contamination_* and unassigned* from a subcluster/lineage object."""
    labels = adata.obs[label_col].astype(str)
    keep = ~(labels.str.startswith("Contamination_")
             | labels.str.lower().eq("unresolved")
             | labels.str.lower().str.startswith("unassigned"))
    n_drop = int((~keep).sum())
    if n_drop:
        print(f"    Dropped {n_drop:,} contaminant/unassigned cells "
              f"({sorted(labels[~keep].unique())})")
    return adata[keep.values].copy()


# ---------------------------------------------------------------------------
# PAGA
# ---------------------------------------------------------------------------

def run_paga(adata, label_col, plot_dir, table_dir, prefix, seed,
             draw_group_age=True):
    """PAGA on the current neighbors graph, grouped by `label_col`.

    Writes connectivities CSV, cell-type PAGA, optional per-age group/pool PAGA,
    connectivity heatmap, edge diagnostics. Returns nothing (no X_umap clobber).
    """
    paga_dir = plot_dir / "paga"
    paga_dir.mkdir(parents=True, exist_ok=True)

    adata.obs[label_col] = adata.obs[label_col].astype("category")
    sc.tl.paga(adata, groups=label_col)

    conn_raw = adata.uns["paga"]["connectivities"]
    conn = pd.DataFrame(
        conn_raw.toarray() if sp.issparse(conn_raw) else np.asarray(conn_raw),
        index=adata.obs[label_col].cat.categories,
        columns=adata.obs[label_col].cat.categories,
    )
    conn.to_csv(table_dir / f"{prefix}_paga_connectivities.csv")

    # PAGA coloured by cell type. sc.pl.paga needs uns['paga']['pos']; compute it.
    sc.pl.paga(adata, show=False)
    fig, ax = plt.subplots(figsize=(8, 7))
    sc.pl.paga(adata, color=label_col, ax=ax, show=False,
               title=f"PAGA — {label_col} connectivity",
               node_size_scale=1.5, edge_width_scale=1.0, fontsize=7, frameon=False)
    safe_savefig(fig, paga_dir / "paga_by_celltype.png")

    # Per-age PAGA coloured by group / pool (recompute on each age slice).
    if draw_group_age and "age" in adata.obs.columns:
        for age in sorted(adata.obs["age"].unique()):
            sub = adata[adata.obs["age"] == age].copy()
            if sub.n_obs < 50:
                continue
            sub.obs[label_col] = sub.obs[label_col].astype("category")
            try:
                sc.tl.paga(sub, groups=label_col)
                sc.pl.paga(sub, show=False)
                for ck in ("group", "pool"):
                    if ck not in sub.obs.columns:
                        continue
                    fig, ax = plt.subplots(figsize=(7, 6))
                    sc.pl.paga(sub, color=ck, ax=ax, show=False, frameon=False,
                               title=f"PAGA age={age}, by {ck}",
                               node_size_scale=1.2, fontsize=7)
                    safe_savefig(fig, paga_dir / f"paga_by_{ck}_age{slugify(age)}.png")
            except Exception as e:
                print(f"    [warn] per-age PAGA age={age}: {e}")
        # restore full-object paga for any downstream use
        adata.obs[label_col] = adata.obs[label_col].astype("category")
        sc.tl.paga(adata, groups=label_col)

    # Connectivity heatmap
    import seaborn as sns
    fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(conn)), max(5, 0.4 * len(conn))))
    sns.heatmap(conn, ax=ax, cmap="Blues", vmin=0, vmax=1, linewidths=0.3,
                annot=len(conn) <= 15, fmt=".2f", annot_kws={"size": 6})
    ax.set_title(f"PAGA connectivities ({label_col})")
    safe_savefig(fig, paga_dir / "paga_transitions_heatmap.png")

    write_paga_edge_diagnostics(adata, label_col, conn, table_dir, prefix)
    print(f"    PAGA done -> {table_dir}/{prefix}_paga_connectivities.csv")


def write_paga_edge_diagnostics(adata, label_col, conn, table_dir, prefix,
                                min_connectivity=0.1, n_shared=15):
    """One row per cell-type edge above min_connectivity with offline-audit info:
    connectivity, cluster sizes, doublet rate, %mt, top shared genes, and an
    ambient_driven flag (>=1/3 of shared top genes are ambient markers)."""
    tissue = adata.uns.get("tissue", "brain")
    ambient = set(AMBIENT_MARKERS.get(tissue, []))
    sym = symbol_to_var(adata)
    var_to_sym = {v: k for k, v in sym.items()}  # var_name -> symbol (best effort)

    if "lognorm" not in adata.layers:
        add_lognorm(adata)
    L = adata.layers["lognorm"]
    L = L.toarray() if sp.issparse(L) else np.asarray(L)
    labels = adata.obs[label_col].astype(str)
    var_names = adata.var_names.astype(str).values

    cl_mean, cl_meta = {}, {}
    for c in conn.index:
        m = (labels == c).values
        if m.sum() == 0:
            continue
        cl_mean[c] = L[m].mean(axis=0)
        meta = {"n_cells": int(m.sum())}
        meta["doublet_rate"] = (
            float((adata.obs.loc[m, "doublet_class"].astype(str) == "doublet").mean())
            if "doublet_class" in adata.obs.columns else np.nan)
        meta["mean_pct_mt"] = (float(adata.obs.loc[m, "pct_counts_mt"].mean())
                               if "pct_counts_mt" in adata.obs.columns else np.nan)
        cl_meta[c] = meta

    cats = list(conn.index)
    rows = []
    for i, a in enumerate(cats):
        for j, b in enumerate(cats):
            if j <= i:
                continue
            w = float(conn.iloc[i, j])
            if w < min_connectivity or a not in cl_mean or b not in cl_mean:
                continue
            shared = np.minimum(cl_mean[a], cl_mean[b])
            top_idx = np.argsort(shared)[::-1][:n_shared]
            top_genes = [var_to_sym.get(g, g) for g in var_names[top_idx]]
            n_amb = sum(1 for g in top_genes if g in ambient)
            rows.append({
                "celltype_A": a, "celltype_B": b, "connectivity": round(w, 4),
                "n_cells_A": cl_meta[a]["n_cells"], "n_cells_B": cl_meta[b]["n_cells"],
                "doublet_rate_A": _r(cl_meta[a]["doublet_rate"]),
                "doublet_rate_B": _r(cl_meta[b]["doublet_rate"]),
                "mean_pct_mt_A": _r(cl_meta[a]["mean_pct_mt"], 3),
                "mean_pct_mt_B": _r(cl_meta[b]["mean_pct_mt"], 3),
                "n_ambient_in_shared": n_amb,
                "ambient_driven": n_amb >= (n_shared // 3),
                "top_shared_genes": ", ".join(top_genes),
            })
    if rows:
        df = (pd.DataFrame(rows).sort_values("connectivity", ascending=False)
              .reset_index(drop=True))
        df.to_csv(table_dir / f"{prefix}_paga_edge_diagnostics.csv", index=False)
        print(f"    Edge diagnostics -> {prefix}_paga_edge_diagnostics.csv "
              f"({len(df)} edges, {int(df['ambient_driven'].sum())} ambient-driven)")
    else:
        print(f"    [info] no edges above connectivity {min_connectivity}")


def _r(x, nd=4):
    return round(x, nd) if (x is not None and not (isinstance(x, float) and np.isnan(x))) else np.nan


# ---------------------------------------------------------------------------
# DPT on one lineage subset (axis computed once on all cells in the subset)
# ---------------------------------------------------------------------------

def compute_dpt(lin, label_col, root_label, seed, order=None):
    """Recompute neighbors+diffmap on lin's own X_scVI, root DPT at the cell
    closest to root_label's centroid. Writes lin.obs['dpt_pseudotime'] in place.

    Orientation guard: if `order` is given and the root cell type's mean
    pseudotime exceeds the terminal type's, the axis is inverted -> flip it so
    pseudotime increases root->terminal. Returns (root_obs_name, flipped) or
    (None, False) if the root label is absent."""
    sc.pp.neighbors(lin, use_rep="X_scVI", random_state=seed)
    sc.tl.diffmap(lin)
    labels = lin.obs[label_col].astype(str)
    mask = labels.str.contains(root_label, case=False, regex=False).values
    if mask.sum() == 0:
        print(f"    [warn] root label '{root_label}' not found in lineage; "
              f"DPT skipped for this root")
        return None, False
    sub = lin[mask]
    centroid = sub.obsm["X_scVI"].mean(axis=0)
    dists = np.linalg.norm(sub.obsm["X_scVI"] - centroid, axis=1)
    root_name = sub.obs_names[int(np.argmin(dists))]
    lin.uns["iroot"] = lin.obs_names.get_loc(root_name)
    sc.tl.dpt(lin)

    flipped = False
    if order and len(order) >= 2:
        means = lin.obs.groupby(label_col, observed=True)["dpt_pseudotime"].mean()
        root_ct, term_ct = order[0], order[-1]
        if root_ct in means.index and term_ct in means.index:
            if means[root_ct] > means[term_ct]:
                pt = lin.obs["dpt_pseudotime"].values
                lin.obs["dpt_pseudotime"] = float(np.nanmax(pt)) - pt
                flipped = True
                print(f"    [orient] flipped pseudotime: '{root_ct}' had higher "
                      f"DPT than '{term_ct}' (axis was inverted)")
    return root_name, flipped


def plot_dpt_axis(lin, label_col, markers, sym, plot_dir, root_label, tissue, age_tag):
    """UMAP coloured by pseudotime + cell type, per-celltype violin, and
    marker-vs-pseudotime trends. Uses the lineage's existing X_umap.
    Path includes age_tag so within-age slices don't overwrite each other."""
    pt_dir = plot_dir / "pseudotime" / slugify(age_tag) / slugify(root_label)
    pt_dir.mkdir(parents=True, exist_ok=True)

    if "X_umap" in lin.obsm:
        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        sc.pl.umap(lin, color="dpt_pseudotime", ax=axes[0], show=False, frameon=False,
                   color_map="viridis", size=10, title=f"DPT (root={root_label})")
        sc.pl.umap(lin, color=label_col, ax=axes[1], show=False, frameon=False,
                   legend_fontsize=6, size=10, title="Subtypes")
        safe_savefig(fig, pt_dir / "dpt_umap.png")

    # Per-celltype violin, ordered by median pseudotime
    order = (lin.obs.groupby(label_col, observed=True)["dpt_pseudotime"]
             .median().sort_values().index.tolist())
    lin.obs["_ct_ord"] = pd.Categorical(lin.obs[label_col].astype(str), categories=order)
    fig, ax = plt.subplots(figsize=(max(6, 0.8 * len(order)), 4.5))
    sc.pl.violin(lin, keys="dpt_pseudotime", groupby="_ct_ord", ax=ax, show=False,
                 rotation=45)
    ax.set_title(f"Pseudotime by subtype (root={root_label})")
    safe_savefig(fig, pt_dir / "dpt_violin_by_subtype.png")
    del lin.obs["_ct_ord"]

    # Marker trends: lognorm vs pseudotime, binned mean. Warn-skip missing markers.
    if "lognorm" not in lin.layers:
        add_lognorm(lin)
    present, missing = [], []
    for mk in (markers or []):
        vn = sym.get(mk)
        (present.append((mk, vn)) if vn in set(lin.var_names.astype(str)) else missing.append(mk))
    if missing:
        print(f"    [warn] markers not in panel (skipped): {missing}")
    if present:
        L = lin.layers["lognorm"]
        pt = lin.obs["dpt_pseudotime"].values
        order_idx = np.argsort(pt)
        nbin = 20
        bins = np.linspace(pt.min(), pt.max(), nbin + 1)
        centers = 0.5 * (bins[:-1] + bins[1:])
        which = np.clip(np.digitize(pt, bins) - 1, 0, nbin - 1)
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        for mk, vn in present:
            col = list(lin.var_names.astype(str)).index(vn)
            vals = L[:, col]
            vals = vals.toarray().ravel() if sp.issparse(vals) else np.asarray(vals).ravel()
            means = [vals[which == b].mean() if (which == b).any() else np.nan
                     for b in range(nbin)]
            ax.plot(centers, means, marker="o", ms=3, lw=1.4, label=mk)
        ax.set_xlabel("DPT pseudotime"); ax.set_ylabel("mean lognorm expr")
        ax.set_title(f"Marker trends along pseudotime (root={root_label})")
        ax.legend(fontsize=7, ncol=2)
        safe_savefig(fig, pt_dir / "dpt_marker_trends.png")


# ---------------------------------------------------------------------------
# Per-donor summary + group comparison
# ---------------------------------------------------------------------------

def donor_summary(obs, mature_threshold):
    """Per-donor median pseudotime + mature fraction over the cells in `obs`."""
    rows = []
    for donor, d in obs.groupby("donor_id", observed=True):
        pt = d["dpt_pseudotime"].dropna().values
        if len(pt) == 0:
            continue
        g = classify_group(d["group"].iloc[0])
        rows.append({
            "donor_id": donor,
            "group": g,
            "raw_group": str(d["group"].iloc[0]),
            "age": str(d["age"].iloc[0]) if "age" in d else "NA",
            "sex": str(d["sex"].iloc[0]) if "sex" in d else "NA",
            "n_cells": int(len(pt)),
            "median_pseudotime": float(np.median(pt)),
            "mature_fraction": float((pt > mature_threshold).mean()),
        })
    return pd.DataFrame(rows)


def compare_groups(summ, contrasts, sex_label, group_level, min_donors,
                   reliable_donors, base_note):
    """Run MW-U pairwise (per contrast) + KW omnibus on donor scalars.

    `summ` is a donor-summary frame already filtered to the relevant slice
    (sex stratum, age). Returns a list of result rows."""
    out = []
    metrics = ["median_pseudotime", "mature_fraction"]
    present_groups = [g for g in ("Relaxed", "Early", "Late")
                      if (summ["group"] == g).sum() > 0]

    # Pairwise MW-U
    for cname, (a, b) in contrasts.items():   # a=test, b=reference
        ga = summ.loc[summ["group"] == a]
        gb = summ.loc[summ["group"] == b]
        na, nb = len(ga), len(gb)
        if na < min_donors or nb < min_donors:
            continue
        rel = "ok" if (na >= reliable_donors and nb >= reliable_donors) else "low_n"
        for metric in metrics:
            x, y = ga[metric].values, gb[metric].values
            try:
                U, p = stats.mannwhitneyu(x, y, alternative="two-sided")
                rbc = 1.0 - 2.0 * U / (na * nb)   # rank-biserial
            except Exception as e:
                U, p, rbc = np.nan, np.nan, np.nan
                print(f"    [warn] MW-U {cname}/{metric}: {e}")
            out.append({
                "sex_stratum": sex_label, "group_level": group_level,
                "contrast": cname, "test": "mann_whitney_u", "metric": metric,
                "n_test": na, "n_ref": nb,
                "median_test": _r(float(np.median(x)), 4),
                "median_ref": _r(float(np.median(y)), 4),
                "statistic": _r(float(U), 3), "pvalue": _r(float(p), 6),
                "effect_size": _r(float(rbc), 4), "effect_type": "rank_biserial",
                "reliability": rel, "note": base_note,
            })

    # KW omnibus (>=3 groups)
    if len(present_groups) >= 3:
        sizes = {g: (summ["group"] == g).sum() for g in present_groups}
        if all(n >= min_donors for n in sizes.values()):
            rel = "ok" if all(n >= reliable_donors for n in sizes.values()) else "low_n"
            for metric in metrics:
                data = [summ.loc[summ["group"] == g, metric].values
                        for g in present_groups]
                try:
                    H, p = stats.kruskal(*data)
                    k = len(present_groups); n = sum(len(d) for d in data)
                    eta2 = max(0.0, (H - k + 1) / (n - k)) if n > k else np.nan
                except Exception as e:
                    H, p, eta2 = np.nan, np.nan, np.nan
                    print(f"    [warn] KW {metric}: {e}")
                out.append({
                    "sex_stratum": sex_label, "group_level": group_level,
                    "contrast": "omnibus_3group", "test": "kruskal_wallis",
                    "metric": metric, "n_test": int(sum(sizes.values())), "n_ref": np.nan,
                    "median_test": np.nan, "median_ref": np.nan,
                    "statistic": _r(float(H), 3), "pvalue": _r(float(p), 6),
                    "effect_size": _r(float(eta2), 4), "effect_type": "eta2",
                    "reliability": rel, "note": base_note,
                })
    return out


def plot_per_donor(summ, contrasts, sex_label, plot_dir, root_label, group_level):
    """Donor-level scatter/box: each donor a point, x=group, two panels
    (median pseudotime, mature fraction)."""
    pd_dir = plot_dir / "per_donor" / sex_label
    pd_dir.mkdir(parents=True, exist_ok=True)
    order = [g for g in ("Relaxed", "Early", "Late") if (summ["group"] == g).any()]
    if not order:
        return
    color = {"Relaxed": "gray", "Early": "tab:red", "Late": "tab:blue"}
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, metric, title in zip(
            axes, ["median_pseudotime", "mature_fraction"],
            ["Median pseudotime / donor", "Mature fraction / donor"]):
        for xi, g in enumerate(order):
            vals = summ.loc[summ["group"] == g, metric].values
            if len(vals):
                ax.boxplot(vals, positions=[xi], widths=0.5,
                           showfliers=False, patch_artist=False)
                jit = (np.random.RandomState(0).rand(len(vals)) - 0.5) * 0.18
                ax.scatter(np.full(len(vals), xi) + jit, vals, s=36,
                           color=color[g], edgecolor="k", lw=0.4, zorder=3)
        ax.set_xticks(range(len(order))); ax.set_xticklabels(order)
        ax.set_title(title); ax.set_ylabel(metric)
    fig.suptitle(f"{root_label} — per-donor ({group_level}, sex={sex_label})")
    fig.tight_layout()
    safe_savefig(fig, pd_dir / f"per_donor_{slugify(group_level)}_root_{slugify(root_label)}.png")


# ---------------------------------------------------------------------------
# Lineage driver
# ---------------------------------------------------------------------------

def contrasts_for_groups(present):
    """Available pairwise contrasts given the groups present (canonical labels)."""
    c = {}
    if "Early" in present and "Relaxed" in present:
        c["early_vs_relaxed"] = ("Early", "Relaxed")
    if "Late" in present and "Relaxed" in present:
        c["late_vs_relaxed"] = ("Late", "Relaxed")
    if "Early" in present and "Late" in present:
        c["early_vs_late"] = ("Early", "Late")
    return c


def run_lineage(lin, label_col, lin_name, lin_cfg, plot_dir, table_dir, prefix,
                cfg, seed, within_age):
    """Full Mode-B analysis for one lineage object (already cleaned)."""
    tissue = lin.uns.get("tissue", cfg["tissue"])
    add_lognorm(lin)
    sym = symbol_to_var(lin)
    roots = lin_cfg["root"]
    roots = roots if isinstance(roots, list) else [roots]
    markers = lin_cfg.get("markers", [])
    order = lin_cfg.get("order", [])
    mature_q = float(cfg["trajectory"].get("mature_quantile", 0.66))
    min_donors = int(cfg["trajectory"].get("min_donors", 2))
    reliable = int(cfg["trajectory"].get("reliable_donors", 3))
    strata = iter_strata(cfg, "sex")
    state_axis = lin_name in STATE_AXIS_LINEAGES

    # Lineage PAGA (subtype connectivity) — recompute neighbors first.
    sc.pp.neighbors(lin, use_rep="X_scVI", random_state=seed)
    run_paga(lin, label_col, plot_dir, table_dir, prefix, seed,
             draw_group_age=not within_age)

    # Determine the age-slices to compute DPT on.
    if within_age:
        age_slices = [(a, lin[lin.obs["age"] == a].copy())
                      for a in sorted(lin.obs["age"].unique())
                      if (lin.obs["age"] == a).sum() >= 20]
    else:
        age_slices = [("all_ages", lin)]   # single cross-age axis

    donor_rows, cmp_rows = [], []
    for age_tag, sub in age_slices:
        if within_age:
            print(f"  [age={age_tag}] {sub.n_obs:,} cells")
        for root_label in roots:
            root_name, flipped = compute_dpt(sub, label_col, root_label, seed, order=order)
            if root_name is None:
                continue
            print(f"    DPT root='{root_label}' -> {root_name}"
                  f"{' [flipped]' if flipped else ''}")
            thr = float(np.quantile(sub.obs["dpt_pseudotime"].dropna(), mature_q))

            plot_dpt_axis(sub, label_col, markers, sym, plot_dir, root_label, tissue, age_tag)

            summ = donor_summary(sub.obs, thr)
            summ["lineage"] = lin_name
            summ["root"] = root_label
            summ["age_slice"] = age_tag
            summ["mature_threshold"] = round(thr, 4)
            donor_rows.append(summ)

            # base note assembly
            note_bits = []
            if flipped:
                note_bits.append("pseudotime_axis_flipped_to_root_terminal")
            if state_axis:
                note_bits.append("state_axis_not_lineage")
            if lin_name == "trophoblast":
                note_bits.append("trophoblast_branched_DPT_labyrinth_arm")

            for sex_label, sex_val in strata:
                s = summ if sex_val is None else summ[summ["sex"] == sex_val]
                if len(s) == 0:
                    continue
                sex_note = note_bits + ([] if sex_label == "combined" else ["sex_stratum_low_n"])

                if within_age:
                    present = [g for g in ("Relaxed", "Early", "Late")
                               if (s["group"] == g).any()]
                    cons = contrasts_for_groups(present)
                    cmp_rows += compare_groups(
                        s, cons, sex_label, f"age={age_tag}", min_donors, reliable,
                        "; ".join(sex_note))
                else:
                    # brain: per-age (clean) + all-ages pooled (pool_age_confound)
                    for age in sorted(s["age"].unique()):
                        sa = s[s["age"] == age]
                        present = [g for g in ("Relaxed", "Early", "Late")
                                   if (sa["group"] == g).any()]
                        cons = contrasts_for_groups(present)
                        cmp_rows += compare_groups(
                            sa, cons, sex_label, f"age={age}", min_donors, reliable,
                            "; ".join(sex_note))
                    present = [g for g in ("Relaxed", "Early", "Late")
                               if (s["group"] == g).any()]
                    cons = contrasts_for_groups(present)
                    cmp_rows += compare_groups(
                        s, cons, sex_label, "all_ages", min_donors, reliable,
                        "; ".join(sex_note + ["pool_age_confound (§2): age co-varies with pool"]))

                # per-donor figure (combined per slice is enough; draw for each stratum)
                plot_per_donor(s, None, sex_label, plot_dir, root_label,
                               age_tag if within_age else "all_ages")

    if donor_rows:
        dsum = pd.concat(donor_rows, ignore_index=True)
        dsum.to_csv(table_dir / f"{prefix}_dpt_per_donor_summary.csv", index=False)
        print(f"  Per-donor summary: {len(dsum)} donor-rows -> "
              f"{prefix}_dpt_per_donor_summary.csv")
    if cmp_rows:
        cdf = pd.DataFrame(cmp_rows)
        cdf.to_csv(table_dir / f"{prefix}_dpt_group_comparison.csv", index=False)
        print(f"  Group comparison: {len(cdf)} rows -> {prefix}_dpt_group_comparison.csv")
        show = cdf[cdf["sex_stratum"] == "combined"]
        if len(show):
            print(show[["group_level", "contrast", "metric", "pvalue",
                        "effect_size", "reliability"]].to_string(index=False))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Phase 8d: trajectory (PAGA + DPT)")
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--subcluster", default=None,
                    help="focal lineage name (Mode B). Omit for whole-tissue PAGA (Mode A).")
    args = ap.parse_args()

    print("\n=== Phase 8d: Trajectory (PAGA + DPT) ===")
    print("  RNA velocity / CellRank: NOT run (10x Flex exon-only; no velocity).")
    cfg = load_config(args.config)
    tissue = cfg["tissue"]
    seed = int(cfg.get("random_seed", 42))
    if "trajectory" not in cfg:
        sys.exit("ERROR: no 'trajectory:' block in config. Add it (see script docstring).")
    tcfg = cfg["trajectory"]
    within_age = bool(tcfg.get("within_age", False))

    results_dir = Path(cfg["results_dir"])
    h5dir = results_dir / "h5ad"

    # -----------------------------------------------------------------------
    # Mode A — whole-tissue PAGA
    # -----------------------------------------------------------------------
    if args.subcluster is None:
        in_path = h5dir / "08_annotated" / "all_samples.h5ad"
        if not in_path.is_file():
            sys.exit(f"ERROR: {in_path} not found (run Phase 7).")
        print(f"\n[Mode A: whole-tissue PAGA]\n  Input: {in_path}")
        adata = sc.read_h5ad(in_path)
        adata.uns["tissue"] = tissue
        key = tcfg.get("whole_paga_key",
                       "celltypist_broad" if tissue == "brain" else "celltype_majority")
        if key not in adata.obs.columns:
            sys.exit(f"ERROR: whole_paga_key '{key}' not in obs.")
        adata = drop_noncelltypes(adata, key)
        print(f"  {adata.n_obs:,} cells | key='{key}' "
              f"({adata.obs[key].nunique()} types)")
        if "X_scVI" not in adata.obsm:
            sys.exit("ERROR: X_scVI missing (run Phase 5).")
        sc.pp.neighbors(adata, use_rep="X_scVI", random_state=seed)
        plot_dir = results_dir / "plots" / "08d_trajectory"
        table_dir = phase_table_dir(cfg, "08d_trajectory")
        run_paga(adata, key, plot_dir, table_dir, "08d_trajectory", seed,
                 draw_group_age=True)
        print(f"\n✓ Mode A done. Plots: {plot_dir}/paga/  (X_umap NOT modified)\n")
        return

    # -----------------------------------------------------------------------
    # Mode B — focal lineage
    # -----------------------------------------------------------------------
    lin_name = args.subcluster
    lineages = tcfg.get("lineages", {})
    if lin_name not in lineages:
        sys.exit(f"ERROR: lineage '{lin_name}' not in trajectory.lineages "
                 f"({sorted(lineages)}).")
    lin_cfg = lineages[lin_name]
    prefix = f"08d_trajectory_subcluster_{lin_name}"
    plot_dir = results_dir / "plots" / prefix
    table_dir = phase_table_dir(cfg, prefix)
    print(f"\n[Mode B: lineage '{lin_name}']")

    if tissue == "brain":
        in_path = h5dir / "08c_subclustered" / f"{lin_name}.h5ad"
        if not in_path.is_file():
            sys.exit(f"ERROR: {in_path} not found.")
        lin = sc.read_h5ad(in_path)
        label_col = "subcluster_name" if "subcluster_name" in lin.obs.columns else "subcluster"
    else:
        # placenta: subset the main object to the trophoblast cell-type list
        in_path = h5dir / "08_annotated" / "all_samples.h5ad"
        if not in_path.is_file():
            sys.exit(f"ERROR: {in_path} not found (run Phase 7).")
        cts = lin_cfg.get("celltypes")
        if not cts:
            sys.exit(f"ERROR: placenta lineage '{lin_name}' needs a 'celltypes' list.")
        full = sc.read_h5ad(in_path)
        label_col = "celltype_majority"
        keep = full.obs[label_col].astype(str).isin(cts)
        lin = full[keep.values].copy()
        del full
        print(f"  Subset to {sorted(cts)}: {lin.n_obs:,} cells")

    lin.uns["tissue"] = tissue
    if "X_scVI" not in lin.obsm:
        sys.exit("ERROR: X_scVI missing in lineage object (run Phase 5).")
    lin = drop_noncelltypes(lin, label_col)
    print(f"  {lin.n_obs:,} cells | label='{label_col}' "
          f"({lin.obs[label_col].nunique()} subtypes) | within_age={within_age}")

    run_lineage(lin, label_col, lin_name, lin_cfg, plot_dir, table_dir, prefix,
                cfg, seed, within_age)

    print(f"\n✓ Mode B done ('{lin_name}'). Plots: {plot_dir}")
    print(f"  Tables: {table_dir}/{prefix}_dpt_*.csv\n")


if __name__ == "__main__":
    main()
