#!/usr/bin/env python
"""
08c_pathways_summary.py — Phase 8c summary plots (CSV-only inputs, no recompute).

Mirrors the 8b / 8b_summary split: 08c_pathways.py writes CSVs; this file plots.

ITERATION: per (sex × level × celltype). The ages (group_level) are ROWS inside
each multi-panel figure (developmental story visible in one figure).

VOLCANOS (headline; pathways AND TFs):
  Per (sex × level × celltype):
    1. solo early_vs_relaxed_per_age   (rows=age × cols=collection)
    2. solo late_vs_relaxed_per_age
    3. solo early_vs_late_per_age
    4. DUAL overlay E-v-R + L-v-R (rows=age × cols=collection): Early on
       bottom-x/left-y (circle), Late on top-x/right-y (square), SHARED NES
       range, Y dynamic (per-panel max, floor 3, cap 20; off-scale = ^ at cap).
       Dots colored by sig+direction (gray ns / red sig&NES>0 / blue sig&NES<0).
       Labels ONLY concordant-sig (sig in BOTH + same NES sign).
  TF versions: same, single collection column.

  M5 panels for placenta-main cell types: labels/rows/cols restricted to
  cell-type-relevant GO terms (config/celltype_go_relevance.yaml); other cell
  types unrestricted (top-N by FDR).

OTHER PLOTS:
  dotplot panels (rows=age × cols=collection), celltype×pathway NES heatmap
  (rows=age × cols=collection), TF×celltype heatmap (rows=age),
  pathway + TF concordance (sex×age 3×3 matrix), per-cell ridges,
  per-cell pathway-localization grid (one panel per cell type, coloured by that
  cell type's top pathway — validates that scores localize sensibly),
  pathway×age trajectory (sex columns; within_group_across_age, brain-only),
  leading-edge heatmaps. (Bubble dropped — duplicated the heatmap.)

Default --sex-strata combined for per-slice plots; the sex-matrix plots
(concordance, per-cell UMAP, trajectory) always use all sexes.

Output: plots/08c_pathways{suffix}/...  (PNG @ 300 DPI; constrained_layout)

Usage:
  uv run python scripts/08c_pathways_summary.py --config config/brain.yaml
  uv run python scripts/08c_pathways_summary.py --config config/brain.yaml --subcluster immune
  uv run python scripts/08c_pathways_summary.py --config config/placenta.yaml
"""

import argparse
import ast
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
import seaborn as sns
import anndata as ad

try:
    from adjustText import adjust_text
    HAVE_ADJUSTTEXT = True
except ImportError:
    HAVE_ADJUSTTEXT = False
    print("[note] adjustText not installed — labels will use fallback placement")

from _utils import load_config


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SIG_FDR = 0.05
FDR_FLOOR = 1e-20
YCAP = -np.log10(FDR_FLOOR)            # 20.0
YFLOOR = 3.0                           # minimum y-top so faint slices aren't flat
MAX_LABEL_CHARS = 42
COLLECTIONS = ["MH", "M2", "M5"]   # M8 (Tabula Muris Senis cell-identity sets) dropped — noise, not pathways

CONTRAST_EVR = "EVR"
CONTRAST_LVR = "LVR"
CONTRAST_EVL = "EVL"
WAA = "WAA"
SOLO_CONTRASTS = [CONTRAST_EVR, CONTRAST_LVR, CONTRAST_EVL]
PRIMARY_CONTRASTS = {CONTRAST_EVR, CONTRAST_LVR}

SLICE_ID = ["sex", "level", "celltype"]     # age is now a row dimension, not a key

AGE_ORDER = ["P1", "4W", "3mo", "E12.5", "E18.5"]
SEX_ORDER = ["combined", "M", "F"]
GROUP_ORDER = ["Relaxed", "Early_Stress", "Late_Stress"]

GROUP_COLORS = {"Relaxed": "#7f7f7f", "Early_Stress": "#d62728", "Late_Stress": "#1f77b4"}
EARLY_C = "#d62728"
LATE_C = "#1f77b4"
NS_C = "lightgray"
DIRECTION_COLORS = {"up": "#d62728", "down": "#1f77b4"}
CONTRAST_LABEL = {CONTRAST_EVR: "Early vs Relaxed",
                  CONTRAST_LVR: "Late vs Relaxed",
                  CONTRAST_EVL: "Early vs Late"}
CONTRAST_SLUG = {CONTRAST_EVR: "early_vs_relaxed",
                 CONTRAST_LVR: "late_vs_relaxed",
                 CONTRAST_EVL: "early_vs_late"}


def contrast_family(c):
    """Map a raw contrast name to a tissue-agnostic family token.

    Brain uses generic names (early_vs_relaxed_per_age); placenta bakes the age
    into the name (early_vs_relaxed_E12.5). Match by prefix so both work.
    """
    c = str(c)
    if c.startswith("early_vs_relaxed"):
        return CONTRAST_EVR
    if c.startswith("late_vs_relaxed"):
        return CONTRAST_LVR
    if c.startswith("early_vs_late"):
        return CONTRAST_EVL
    if c.startswith("within_group_across_age"):
        return WAA
    return "OTHER"

# Collection-name prefixes stripped from pathway labels (longest first).
STRIP_PREFIXES = [
    "TABULA_MURIS_SENIS_", "HALLMARK_", "REACTOME_", "BIOCARTA_",
    "GOBP_", "GOMF_", "GOCC_", "KEGGLEGACY_", "KEGGMEDICUS_", "KEGG_",
    "WP_", "PID_", "NABA_", "DESCARTES_",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe(name) -> str:
    return re.sub(r"[^0-9A-Za-z._-]+", "_", str(name)).strip("_")


def parse_pair(pair_repr):
    if pair_repr is None or (isinstance(pair_repr, float) and pd.isna(pair_repr)):
        return []
    s = str(pair_repr).strip()
    if not s or s.lower() == "nan":
        return []
    try:
        v = ast.literal_eval(s)
        return [str(x) for x in v] if isinstance(v, (list, tuple)) else []
    except (ValueError, SyntaxError):
        return []


def strip_prefix(s):
    s = str(s)
    for p in STRIP_PREFIXES:
        if s.startswith(p):
            return s[len(p):]
    return s


def disp(s, n=MAX_LABEL_CHARS):
    """Display form for a pathway name: strip collection prefix + truncate."""
    s = strip_prefix(s)
    return s if len(s) <= n else s[: n - 1] + "…"


def truncate(s, n=MAX_LABEL_CHARS):
    s = str(s)
    return s if len(s) <= n else s[: n - 1] + "…"


def order_ages(ages):
    ages = [a for a in ages if a is not None and str(a) != "nan"]
    known = [a for a in AGE_ORDER if a in ages]
    rest = sorted([a for a in ages if a not in AGE_ORDER])
    return known + rest


def order_groups(groups):
    """group_level for within_group_across_age = stress group (Relaxed/Early/Late)."""
    g = set(groups)
    return [x for x in GROUP_ORDER if x in g] + [x for x in groups if x not in GROUP_ORDER]


def neglog10_fdr(fdr_series):
    f = pd.to_numeric(fdr_series, errors="coerce").to_numpy(dtype=float)
    offscale = f < FDR_FLOOR
    f = np.clip(f, FDR_FLOOR, 1.0)
    y = np.nan_to_num(-np.log10(f), nan=0.0)
    return y, offscale


def ylim_top(*yarrays):
    arrs = [np.asarray(y, dtype=float).ravel() for y in yarrays
            if y is not None and len(y)]
    if not arrs:
        return YFLOOR
    vals = np.concatenate(arrs)
    if vals.size == 0:
        return YFLOOR
    m = float(np.nanmax(vals))
    return max(YFLOOR, min(YCAP, m) * 1.05)


def sigdir_colors(x, fdr):
    """Per-point color by significance + direction. x = effect (NES/activity)."""
    x = np.asarray(x, dtype=float)
    fdr = pd.to_numeric(pd.Series(fdr), errors="coerce").to_numpy(dtype=float)
    sig = fdr < SIG_FDR
    cols = np.full(len(x), NS_C, dtype=object)
    cols[sig & (x > 0)] = DIRECTION_COLORS["up"]
    cols[sig & (x < 0)] = DIRECTION_COLORS["down"]
    return cols, sig


def safe_vmax(arr):
    """Finite, strictly-positive vmax for a diverging colour norm."""
    a = np.abs(np.asarray(arr, dtype=float))
    m = np.nanmax(a) if np.isfinite(a).any() else 1.0
    if not np.isfinite(m) or m <= 0:
        m = 1.0
    return float(m)


def safe_savefig(fig, path: Path, dpi=300):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def add_labels(ax, xs, ys, names, fontsize=5.8, bold=False, prefix_strip=True):
    fn = disp if prefix_strip else (lambda s: truncate(str(s)))
    texts = []
    for x, y, nm in zip(xs, ys, names):
        texts.append(ax.text(x, y, fn(nm), fontsize=fontsize, ha="left",
                             va="bottom", fontweight="bold" if bold else "normal"))
    if HAVE_ADJUSTTEXT and len(texts) > 1:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.4, alpha=0.7),
                    expand=(1.3, 1.4), force_text=(0.4, 0.5))
    return texts


def slice_dir(plot_root, sex, level, celltype=None):
    d = plot_root / _safe(sex) / _safe(level)
    if celltype is not None:
        d = d / _safe(celltype)
    return d


def slice_title(sex, level, celltype=None, extra=""):
    parts = [f"sex={sex}", f"level={level}"]
    if celltype is not None:
        parts.append(f"ct={celltype}")
    s = " | ".join(parts)
    return f"{s}{('  ·  ' + extra) if extra else ''}"


def _tf_fdr_col(df):
    return "FDR_ctx_celltype" if "FDR_ctx_celltype" in df.columns else "FDR"


# ---------------------------------------------------------------------------
# GO relevance (placenta-main M5 restriction)
# ---------------------------------------------------------------------------

def load_go_relevance(tissue, cfg_dir):
    """Return ({celltype: set(go_terms)}, tissue_union_set or None).

    Empty dict / None means no restriction (all cell types -> top-N by FDR).
    """
    path = cfg_dir / "celltype_go_relevance.yaml"
    if not path.is_file():
        print(f"  [GO relevance] {path} absent — M5 unrestricted")
        return {}, None
    with open(path) as fh:
        rel = yaml.safe_load(fh) or {}
    block = rel.get(tissue, {})
    go_map = {}
    for ct, spec in block.items():
        terms = set(spec.get("go_terms", [])) if isinstance(spec, dict) else set()
        if terms:
            go_map[ct] = terms
    if not go_map:
        print(f"  [GO relevance] no entries for tissue '{tissue}' — M5 unrestricted")
        return {}, None
    union = set().union(*go_map.values())
    print(f"  [GO relevance] {len(go_map)} cell types, union {len(union)} GO terms")
    return go_map, union


def _coll_label_policy(coll, celltype, go_map):
    """Return (do_label, label_allow_set_or_None) for a volcano collection panel."""
    if coll == "M5" and celltype in go_map:
        return True, go_map[celltype]           # placenta-main: relevant GO only
    return True, None                           # all sig


# ---------------------------------------------------------------------------
# Volcano scatter (single axis)
# ---------------------------------------------------------------------------

def _scatter_volcano(ax, df, x_col, fdr_col, name_col, max_labels,
                     label=True, label_allow=None, label_bold=False):
    if df is None or df.empty:
        ax.text(0.5, 0.5, "—", ha="center", va="center",
                transform=ax.transAxes, color="gray", fontsize=8)
        return 0
    y, off = neglog10_fdr(df[fdr_col])
    x = pd.to_numeric(df[x_col], errors="coerce").to_numpy(dtype=float)
    cols, sig = sigdir_colors(x, df[fdr_col])
    ytop = ylim_top(y)
    ax.scatter(x[~sig], y[~sig], s=7, c=NS_C, alpha=0.4, linewidth=0)
    ax.scatter(x[sig], y[sig], s=22, c=cols[sig], alpha=0.85,
               edgecolors="black", linewidth=0.25)
    osig = sig & off
    if osig.any():
        ax.scatter(x[osig], np.full(osig.sum(), ytop), marker="^", s=34,
                   c=cols[osig], edgecolors="black", linewidth=0.3, zorder=5)
    ax.axhline(-np.log10(SIG_FDR), color="gray", linestyle="--", linewidth=0.5)
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.5)
    ax.set_ylim(0, ytop)
    n_sig = int(sig.sum())
    if label and n_sig:
        names = df[name_col].to_numpy()
        sig_idx = np.where(sig)[0]
        if label_allow is not None:
            sig_idx = np.array([i for i in sig_idx if names[i] in label_allow], dtype=int)
        if sig_idx.size:
            fdr_vals = pd.to_numeric(df[fdr_col], errors="coerce").to_numpy()
            order = sig_idx[np.argsort(fdr_vals[sig_idx])][:max_labels]
            add_labels(ax, x[order], y[order], names[order], fontsize=5.6,
                       bold=label_bold, prefix_strip=(name_col == "source"))
    return n_sig


# ---------------------------------------------------------------------------
# 1–3. Solo pathway volcano (rows=age × cols=collection)
# ---------------------------------------------------------------------------

def plot_solo_volcano_pathways(slice_block, contrast, outpath, max_labels, go_map):
    df = slice_block[slice_block["cfam"] == contrast]
    if df.empty or not (df["FDR"] < SIG_FDR).any():
        return
    ct = df.iloc[0]["celltype"]
    ages = order_ages(df["group_level"].unique())
    nrow, ncol = len(ages), len(COLLECTIONS)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.4 * ncol, 3.3 * nrow),
                             constrained_layout=True, squeeze=False)
    for i, age in enumerate(ages):
        for j, coll in enumerate(COLLECTIONS):
            ax = axes[i][j]
            d = df[(df["group_level"] == age) & (df["collection"] == coll)]
            do_label, allow = _coll_label_policy(coll, ct, go_map)
            n_sig = _scatter_volcano(ax, d, "NES", "FDR", "source", max_labels,
                                     label=do_label, label_allow=allow)
            if i == 0:
                ax.set_title(coll, fontsize=10)
            if i == nrow - 1:
                ax.set_xlabel("NES", fontsize=8)
            if j == 0:
                ax.set_ylabel(f"{age}\n−log₁₀(FDR)", fontsize=8)
            ax.tick_params(labelsize=6)
            ax.text(0.97, 0.95, f"sig={n_sig}", transform=ax.transAxes,
                    ha="right", va="top", fontsize=6)
    r = df.iloc[0]
    fig.suptitle(slice_title(r["sex"], r["level"], ct,
                 extra=CONTRAST_LABEL.get(contrast, contrast)), fontsize=11)
    safe_savefig(fig, outpath)


# ---------------------------------------------------------------------------
# 4. Dual pathway volcano (rows=age × cols=collection, twin axes)
# ---------------------------------------------------------------------------

def _dual_cell(ax, ec, lc, max_labels, allow_m5=None):
    """One dual-overlay cell. Early=circle (primary), Late=square (twin); both
    colored by sig+direction. Returns n_concordant."""
    have_e = ec is not None and len(ec)
    have_l = lc is not None and len(lc)
    if not have_e and not have_l:
        ax.text(0.5, 0.5, "—", ha="center", va="center",
                transform=ax.transAxes, color="gray", fontsize=8)
        return 0
    nes_parts = []
    if have_e:
        nes_parts.append(ec["NES"].astype(float))
    if have_l:
        nes_parts.append(lc["NES"].astype(float))
    nes_all = pd.concat(nes_parts)
    nmax = max(abs(nes_all.min()), abs(nes_all.max()), 1.0) * 1.1

    ye, offe = neglog10_fdr(ec["FDR"]) if have_e else (np.array([]), np.array([], bool))
    yl, offl = neglog10_fdr(lc["FDR"]) if have_l else (np.array([]), np.array([], bool))
    ytop = ylim_top(ye, yl)

    if have_e:
        xe = ec["NES"].astype(float).to_numpy()
        ce, se = sigdir_colors(xe, ec["FDR"])
        ax.scatter(xe, ye, s=20, c=ce, marker="o", alpha=0.8,
                   edgecolors="black", linewidth=0.2)
        if (se & offe).any():
            o = se & offe
            ax.scatter(xe[o], np.full(o.sum(), ytop), marker="^", s=30,
                       c=ce[o], edgecolors="black", linewidth=0.25, zorder=5)
    ax.set_xlim(-nmax, nmax)
    ax.set_ylim(0, ytop)
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.4)
    ax.axhline(-np.log10(SIG_FDR), color="gray", linestyle=":", linewidth=0.4)
    ax.tick_params(labelsize=6)

    ax_r = ax.twinx()
    ax_t = ax_r.twiny()
    if have_l:
        xl = lc["NES"].astype(float).to_numpy()
        cl, sl = sigdir_colors(xl, lc["FDR"])
        ax_t.scatter(xl, yl, s=20, c=cl, marker="s", alpha=0.8,
                     edgecolors="black", linewidth=0.2)
        if (sl & offl).any():
            o = sl & offl
            ax_t.scatter(xl[o], np.full(o.sum(), ytop), marker="^", s=30,
                         c=cl[o], edgecolors="black", linewidth=0.25, zorder=5)
    ax_r.set_ylim(0, ytop)
    ax_t.set_xlim(-nmax, nmax)
    ax_r.tick_params(labelsize=6)
    ax_t.tick_params(labelsize=6)

    n_conc = 0
    if have_e and have_l:
        m = pd.merge(
            ec[["source", "NES", "FDR"]].rename(columns={"NES": "NES_e", "FDR": "FDR_e"}),
            lc[["source", "NES", "FDR"]].rename(columns={"NES": "NES_l", "FDR": "FDR_l"}),
            on="source", how="inner")
        m = m[(m["FDR_e"] < SIG_FDR) & (m["FDR_l"] < SIG_FDR)
              & (np.sign(m["NES_e"]) == np.sign(m["NES_l"]))]
        if allow_m5 is not None:
            m = m[m["source"].isin(allow_m5)]
        n_conc = len(m)
        if n_conc:
            m = m.reindex(m[["FDR_e", "FDR_l"]].max(axis=1).sort_values().index).head(max_labels)
            yy, _ = neglog10_fdr(m["FDR_e"])
            add_labels(ax, m["NES_e"].astype(float).to_numpy(), yy,
                       m["source"].to_numpy(), fontsize=5.6, bold=True)
    return n_conc


def plot_dual_volcano_pathways(slice_block, outpath, max_labels, go_map):
    e = slice_block[slice_block["cfam"] == CONTRAST_EVR]
    l = slice_block[slice_block["cfam"] == CONTRAST_LVR]
    if e.empty and l.empty:
        return
    e_sig = (e["FDR"] < SIG_FDR).any() if not e.empty else False
    l_sig = (l["FDR"] < SIG_FDR).any() if not l.empty else False
    if not (e_sig or l_sig):
        return
    ct = (e.iloc[0] if not e.empty else l.iloc[0])["celltype"]
    ages = order_ages(pd.concat([e["group_level"], l["group_level"]]).unique())
    nrow, ncol = len(ages), len(COLLECTIONS)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.6 * ncol, 3.4 * nrow),
                             constrained_layout=True, squeeze=False)
    for i, age in enumerate(ages):
        for j, coll in enumerate(COLLECTIONS):
            ax = axes[i][j]
            ec = e[(e["group_level"] == age) & (e["collection"] == coll)]
            lc = l[(l["group_level"] == age) & (l["collection"] == coll)]
            if coll == "M5":
                allow_m5 = go_map.get(ct)        # None -> all concordant
            else:
                allow_m5 = None
            n_conc = _dual_cell(ax, ec, lc, max_labels, allow_m5=allow_m5)
            if i == 0:
                ax.set_title(coll, fontsize=10)
            if j == 0:
                ax.set_ylabel(f"{age}\n−log₁₀FDR (E)", fontsize=7, color=EARLY_C)
            ax.text(0.97, 0.96, f"conc={n_conc}", transform=ax.transAxes,
                    ha="right", va="top", fontsize=6)
    legend_el = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
               markeredgecolor="black", markersize=8, label="Early vs Rel (●)"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="gray",
               markeredgecolor="black", markersize=8, label="Late vs Rel (■)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=DIRECTION_COLORS["up"],
               markersize=8, label="sig up"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=DIRECTION_COLORS["down"],
               markersize=8, label="sig down"),
    ]
    fig.legend(handles=legend_el, loc="lower center", ncol=4, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, -0.01))
    r = (e.iloc[0] if not e.empty else l.iloc[0])
    fig.suptitle(slice_title(r["sex"], r["level"], ct,
                 extra="DUAL · ●Early/■Late · color=sig+dir · labels=concordant-sig"),
                 fontsize=11)
    safe_savefig(fig, outpath)


# ---------------------------------------------------------------------------
# TF volcanos (rows=age, single column)
# ---------------------------------------------------------------------------

def plot_solo_volcano_tf(tf_block, contrast, outpath, max_labels):
    df = tf_block[tf_block["cfam"] == contrast]
    if df.empty:
        return
    fdr_col = _tf_fdr_col(df)
    if not (pd.to_numeric(df[fdr_col], errors="coerce") < SIG_FDR).any():
        return
    ages = order_ages(df["group_level"].unique())
    fig, axes = plt.subplots(len(ages), 1, figsize=(7.5, 4.6 * len(ages)),
                             constrained_layout=True, squeeze=False)
    for i, age in enumerate(ages):
        ax = axes[i][0]
        d = df[df["group_level"] == age]
        n_sig = _scatter_volcano(ax, d, "activity_score", fdr_col, "TF", max_labels,
                                 label=True, label_bold=True)
        ax.set_ylabel(f"{age}\n−log₁₀({fdr_col})", fontsize=8)
        if i == len(ages) - 1:
            ax.set_xlabel("TF activity (ULM t-stat)", fontsize=9)
        ax.set_title(f"sig={n_sig}", fontsize=8, loc="right")
    r = df.iloc[0]
    fig.suptitle(slice_title(r["sex"], r["level"], r["celltype"],
                 extra=f"TF · {CONTRAST_LABEL.get(contrast, contrast)}"), fontsize=11)
    safe_savefig(fig, outpath)


def plot_dual_volcano_tf(tf_block, outpath, max_labels):
    e = tf_block[tf_block["cfam"] == CONTRAST_EVR]
    l = tf_block[tf_block["cfam"] == CONTRAST_LVR]
    if e.empty and l.empty:
        return
    fdr_col = _tf_fdr_col(tf_block)
    e_sig = (pd.to_numeric(e[fdr_col], errors="coerce") < SIG_FDR).any() if not e.empty else False
    l_sig = (pd.to_numeric(l[fdr_col], errors="coerce") < SIG_FDR).any() if not l.empty else False
    if not (e_sig or l_sig):
        return
    ct = (e.iloc[0] if not e.empty else l.iloc[0])["celltype"]
    ages = order_ages(pd.concat([e["group_level"], l["group_level"]]).unique())
    fig, axes = plt.subplots(len(ages), 1, figsize=(8.0, 5.2 * len(ages)),
                             constrained_layout=True, squeeze=False)
    for i, age in enumerate(ages):
        ax = axes[i][0]
        ec = e[e["group_level"] == age]
        lc = l[l["group_level"] == age]
        ec = (ec[["TF", "activity_score", fdr_col]]
              .rename(columns={"TF": "source", "activity_score": "NES", fdr_col: "FDR"})
              if len(ec) else ec)
        lc = (lc[["TF", "activity_score", fdr_col]]
              .rename(columns={"TF": "source", "activity_score": "NES", fdr_col: "FDR"})
              if len(lc) else lc)
        _dual_cell(ax, ec, lc, max_labels, allow_m5=None)
        ax.set_ylabel(f"{age}\n−log₁₀FDR (E)", fontsize=8, color=EARLY_C)
    legend_el = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
               markeredgecolor="black", markersize=8, label="Early (●)"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="gray",
               markeredgecolor="black", markersize=8, label="Late (■)"),
    ]
    fig.legend(handles=legend_el, loc="lower center", ncol=2, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, -0.01))
    r = (e.iloc[0] if not e.empty else l.iloc[0])
    fig.suptitle(slice_title(r["sex"], r["level"], ct,
                 extra="TF DUAL · ●Early/■Late · labels=concordant-sig"), fontsize=11)
    safe_savefig(fig, outpath)


# ---------------------------------------------------------------------------
# Dotplot panels (rows=age × cols=collection)
# ---------------------------------------------------------------------------

def plot_dotplot_panels(slice_block, contrast, outpath, top_n, go_map):
    df = slice_block[slice_block["cfam"] == contrast]
    if df.empty or not (df["FDR"] < SIG_FDR).any():
        return
    ct = df.iloc[0]["celltype"]
    ages = order_ages(df["group_level"].unique())
    nrow, ncol = len(ages), len(COLLECTIONS)
    fig, axes = plt.subplots(nrow, ncol,
                             figsize=(5.6 * ncol, max(3.0, 0.22 * top_n) * nrow),
                             constrained_layout=True, squeeze=False)
    used = False
    for i, age in enumerate(ages):
        for j, coll in enumerate(COLLECTIONS):
            ax = axes[i][j]
            d = df[(df["group_level"] == age) & (df["collection"] == coll)
                   & (df["FDR"] < SIG_FDR)].copy()
            if coll == "M5" and ct in go_map:
                d = d[d["source"].isin(go_map[ct])]
            if d.empty:
                ax.text(0.5, 0.5, "—", ha="center", va="center",
                        transform=ax.transAxes, color="gray", fontsize=8)
                ax.set_xticks([]); ax.set_yticks([])
                if i == 0:
                    ax.set_title(coll, fontsize=10)
                if j == 0:
                    ax.set_ylabel(age, fontsize=9)
                continue
            used = True
            top = d.sort_values("FDR").head(top_n).iloc[::-1]
            nlogq, _ = neglog10_fdr(top["FDR"])
            y = np.arange(len(top))
            cols, _ = sigdir_colors(top["NES"].to_numpy(), top["FDR"])
            sizes = 30 + 22 * np.clip(nlogq, 0, YCAP)
            ax.scatter(top["NES"], y, s=sizes, c=cols, edgecolors="black",
                       linewidth=0.3, alpha=0.9)
            ax.set_yticks(y)
            ax.set_yticklabels([disp(s, 38) for s in top["source"]], fontsize=6)
            ax.axvline(0, color="gray", linestyle="--", linewidth=0.5)
            ax.tick_params(axis="x", labelsize=6)
            ax.grid(axis="x", linestyle=":", linewidth=0.4, alpha=0.6)
            if i == 0:
                ax.set_title(coll, fontsize=10)
            if j == 0:
                ax.set_ylabel(age, fontsize=9)
            if i == nrow - 1:
                ax.set_xlabel("NES", fontsize=8)
    if not used:
        plt.close(fig); return
    r = df.iloc[0]
    fig.suptitle(slice_title(r["sex"], r["level"], ct,
                 extra=CONTRAST_LABEL.get(contrast, contrast)), fontsize=11)
    safe_savefig(fig, outpath)


# ---------------------------------------------------------------------------
# Cell-type × pathway NES heatmap (rows=age × cols=collection)
# ---------------------------------------------------------------------------

def plot_celltype_pathway_heatmap(block, outpath, top_n, go_union):
    if not (block["FDR"] < SIG_FDR).any():
        return
    ages = order_ages(block["group_level"].unique())
    nct = max(block["celltype"].nunique(), 1)
    nrow, ncol = len(ages), len(COLLECTIONS)
    fig, axes = plt.subplots(nrow, ncol,
                             figsize=(6.5 * ncol, max(3.0, 0.32 * nct + 1.5) * nrow),
                             constrained_layout=True, squeeze=False)
    used = False
    for i, age in enumerate(ages):
        for j, coll in enumerate(COLLECTIONS):
            ax = axes[i][j]
            d = block[(block["group_level"] == age) & (block["collection"] == coll)]
            sig = d[d["FDR"] < SIG_FDR]
            if coll == "M5" and go_union is not None:
                sig = sig[sig["source"].isin(go_union)]
            if sig.empty:
                ax.text(0.5, 0.5, "—", ha="center", va="center",
                        transform=ax.transAxes, color="gray", fontsize=8)
                ax.set_xticks([]); ax.set_yticks([])
                if i == 0:
                    ax.set_title(coll, fontsize=10)
                if j == 0:
                    ax.set_ylabel(age, fontsize=9)
                continue
            used = True
            chosen = (sig.groupby("source")["FDR"].min()
                         .sort_values().head(top_n).index.tolist())
            sub = d[d["source"].isin(chosen)]
            nes = sub.pivot_table(index="celltype", columns="source", values="NES", aggfunc="first")
            fdr = sub.pivot_table(index="celltype", columns="source", values="FDR", aggfunc="first")
            rows = nes.abs().mean(axis=1).sort_values(ascending=False).index
            nes = nes.reindex(index=rows, columns=chosen)
            fdr = fdr.reindex(index=rows, columns=chosen)
            vmax = safe_vmax(nes.values)
            im = ax.imshow(nes.values, cmap="RdBu_r",
                           norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), aspect="auto")
            for (yi, xj), f in np.ndenumerate(fdr.values):
                if not np.isnan(f) and f < SIG_FDR:
                    ax.add_patch(mpatches.Rectangle((xj - 0.5, yi - 0.5), 1, 1,
                                 fill=False, edgecolor="black", linewidth=1.1))
            ax.set_xticks(range(len(chosen)))
            ax.set_xticklabels([disp(c, 32) for c in chosen],
                               rotation=70, ha="right", fontsize=5.8)
            ax.set_yticks(range(len(rows)))
            ax.set_yticklabels(rows, fontsize=6.5)
            if i == 0:
                ax.set_title(coll, fontsize=10)
            if j == 0:
                ax.set_ylabel(age, fontsize=9)
            fig.colorbar(im, ax=ax, shrink=0.5, label="NES")
    if not used:
        plt.close(fig); return
    r = block.iloc[0]
    fig.suptitle(slice_title(r["sex"], r["level"],
                 extra=f"{CONTRAST_LABEL.get(r['cfam'], r['cfam'])} · "
                       f"black=FDR<{SIG_FDR}"), fontsize=11)
    safe_savefig(fig, outpath)


def plot_tf_celltype_heatmap(tf_block, outpath, top_n):
    if tf_block.empty:
        return
    fdr_col = _tf_fdr_col(tf_block)
    if not (pd.to_numeric(tf_block[fdr_col], errors="coerce") < SIG_FDR).any():
        return
    ages = order_ages(tf_block["group_level"].unique())
    ntf = max(tf_block["TF"].nunique(), 1)
    nct = max(tf_block["celltype"].nunique(), 1)
    fig, axes = plt.subplots(len(ages), 1,
                             figsize=(max(8, 0.34 * ntf + 4),
                                      max(3.0, 0.32 * nct + 1.5) * len(ages)),
                             constrained_layout=True, squeeze=False)
    used = False
    for i, age in enumerate(ages):
        ax = axes[i][0]
        d = tf_block[tf_block["group_level"] == age]
        sig = d[d[fdr_col] < SIG_FDR]
        if sig.empty:
            ax.text(0.5, 0.5, "—", ha="center", va="center",
                    transform=ax.transAxes, color="gray", fontsize=8)
            ax.set_xticks([]); ax.set_yticks([]); ax.set_ylabel(age, fontsize=9)
            continue
        used = True
        chosen = sig.groupby("TF")[fdr_col].min().sort_values().head(top_n).index.tolist()
        sub = d[d["TF"].isin(chosen)]
        act = sub.pivot_table(index="celltype", columns="TF", values="activity_score", aggfunc="first")
        fdr = sub.pivot_table(index="celltype", columns="TF", values=fdr_col, aggfunc="first")
        rows = act.abs().mean(axis=1).sort_values(ascending=False).index
        act = act.reindex(index=rows, columns=chosen)
        fdr = fdr.reindex(index=rows, columns=chosen)
        vmax = safe_vmax(act.values)
        im = ax.imshow(act.values, cmap="RdBu_r",
                       norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), aspect="auto")
        for (yi, xj), f in np.ndenumerate(fdr.values):
            if not np.isnan(f) and f < SIG_FDR:
                ax.add_patch(mpatches.Rectangle((xj - 0.5, yi - 0.5), 1, 1,
                             fill=False, edgecolor="black", linewidth=1.1))
        ax.set_xticks(range(len(chosen)))
        ax.set_xticklabels(chosen, rotation=70, ha="right", fontsize=6.5)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(rows, fontsize=7)
        ax.set_ylabel(age, fontsize=9)
        fig.colorbar(im, ax=ax, shrink=0.6, label="activity")
    if not used:
        plt.close(fig); return
    r = tf_block.iloc[0]
    fig.suptitle(slice_title(r["sex"], r["level"],
                 extra=f"TF · {CONTRAST_LABEL.get(r['cfam'], r['cfam'])} · "
                       f"black={fdr_col}<{SIG_FDR}"), fontsize=11)
    safe_savefig(fig, outpath)


# ---------------------------------------------------------------------------
# Concordance: sex × age 3×3 matrix (per level × celltype)
# ---------------------------------------------------------------------------

def plot_concordance_matrix(df_evr, df_lvr, key, val, fdr, outpath, max_labels,
                            title, prefix_strip):
    sexes = [s for s in SEX_ORDER if s in set(df_evr["sex"]).union(df_lvr["sex"])]
    ages = order_ages(pd.concat([df_evr["group_level"], df_lvr["group_level"]]).unique())
    if not sexes or not ages:
        return
    fig, axes = plt.subplots(len(sexes), len(ages),
                             figsize=(4.0 * len(ages), 4.0 * len(sexes)),
                             constrained_layout=True, squeeze=False)
    any_cell = False
    for i, sx in enumerate(sexes):
        for j, age in enumerate(ages):
            ax = axes[i][j]
            a = df_evr[(df_evr["sex"] == sx) & (df_evr["group_level"] == age)]
            b = df_lvr[(df_lvr["sex"] == sx) & (df_lvr["group_level"] == age)]
            if i == 0:
                ax.set_title(age, fontsize=10)
            if j == 0:
                ax.set_ylabel(f"sex={sx}\n{val} · Late", fontsize=8)
            if a.empty or b.empty:
                ax.text(0.5, 0.5, "—", ha="center", va="center",
                        transform=ax.transAxes, color="gray", fontsize=8)
                ax.set_xticks([]); ax.set_yticks([])
                continue
            m = pd.merge(a[[key, val, fdr]].rename(columns={val: "va", fdr: "fa"}),
                         b[[key, val, fdr]].rename(columns={val: "vb", fdr: "fb"}),
                         on=key, how="inner")
            if m.empty:
                ax.text(0.5, 0.5, "—", ha="center", va="center",
                        transform=ax.transAxes, color="gray", fontsize=8)
                continue
            any_cell = True
            m["sa"] = m["fa"] < SIG_FDR
            m["sb"] = m["fb"] < SIG_FDR
            m["state"] = "ns"
            m.loc[m["sa"] & m["sb"], "state"] = "both"
            m.loc[m["sa"] & ~m["sb"], "state"] = "early_only"
            m.loc[~m["sa"] & m["sb"], "state"] = "late_only"
            palette = {"ns": NS_C, "both": "#9467bd",
                       "early_only": EARLY_C, "late_only": LATE_C}
            for st, col in palette.items():
                s = m[m["state"] == st]
                ax.scatter(s["va"], s["vb"], c=col, s=10 if st == "ns" else 22,
                           alpha=0.6 if st == "ns" else 0.85,
                           edgecolors="none" if st == "ns" else "black", linewidth=0.25)
            lo = float(min(m["va"].min(), m["vb"].min()))
            hi = float(max(m["va"].max(), m["vb"].max()))
            pad = 0.05 * (hi - lo) if hi > lo else 0.5
            lo, hi = lo - pad, hi + pad
            ax.plot([lo, hi], [lo, hi], color="gray", linestyle="--", linewidth=0.6)
            ax.axhline(0, color="gray", linestyle=":", linewidth=0.4)
            ax.axvline(0, color="gray", linestyle=":", linewidth=0.4)
            ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
            pool = m[m["state"] == "both"].copy()
            if not pool.empty:
                pool["ext"] = pool["va"].abs() + pool["vb"].abs()
                top = pool.sort_values("ext", ascending=False).head(max_labels)
                add_labels(ax, top["va"], top["vb"], top[key], fontsize=5.6,
                           bold=True, prefix_strip=prefix_strip)
            ax.tick_params(labelsize=6)
            ax.text(0.04, 0.96, f"both={int((m['state']=='both').sum())}",
                    transform=ax.transAxes, ha="left", va="top", fontsize=6)
            if i == len(sexes) - 1:
                ax.set_xlabel(f"{val} · Early", fontsize=8)
    if not any_cell:
        plt.close(fig); return
    fig.suptitle(title, fontsize=11)
    safe_savefig(fig, outpath)


# ---------------------------------------------------------------------------
# Per-cell ridges (per sex × celltype; not GO-filtered)
# ---------------------------------------------------------------------------

def _pc_celltype_col(pc):
    for c in ("celltype", "celltypist_broad", "celltype_majority", "subcluster_name"):
        if c in pc.obs.columns:
            return c
    return None


def plot_per_cell_ridges(pc, gsea_df, plot_root, top_n, group_col="group"):
    ct_col = _pc_celltype_col(pc)
    if ct_col is None or group_col not in pc.obs.columns:
        print("  [skip ridges] missing celltype/group col"); return
    chosen = [p for p in gsea_df[gsea_df["FDR"] < SIG_FDR]["source"].value_counts().index[:top_n]
              if p in pc.var_names]
    if not chosen:
        print("  [skip ridges] no chosen pathway present"); return
    X = pc[:, chosen].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    scores = pd.DataFrame(X, index=pc.obs_names, columns=chosen)
    meta = pc.obs[[c for c in ["sex", group_col, ct_col] if c in pc.obs.columns]].copy()
    meta = meta.rename(columns={ct_col: "celltype"})
    sexes = sorted(meta["sex"].dropna().unique()) if "sex" in meta else ["combined"]
    for sx in sexes:
        for ct in sorted(meta["celltype"].dropna().unique()):
            mask = (meta["celltype"] == ct)
            if "sex" in meta:
                mask &= (meta["sex"] == sx)
            if int(mask.sum()) < 50:
                continue
            long = scores.loc[meta.index[mask]].stack().reset_index()
            long.columns = ["cell", "pathway", "score"]
            long["group"] = long["cell"].map(meta[group_col])
            long["pathway"] = long["pathway"].map(lambda p: disp(p, 38))
            order = [disp(p, 38) for p in chosen]
            hue_order = [g for g in GROUP_ORDER if g in long["group"].unique()]
            fig, ax = plt.subplots(figsize=(10, max(4, 0.4 * len(chosen) + 1)),
                                   constrained_layout=True)
            sns.violinplot(data=long, x="score", y="pathway", hue="group", order=order,
                           hue_order=hue_order,
                           palette={k: v for k, v in GROUP_COLORS.items() if k in hue_order},
                           inner="quartile", linewidth=0.5, cut=0, density_norm="width", ax=ax)
            ax.set_xlabel("AUCell score"); ax.set_ylabel("")
            ax.set_title(f"sex={sx} | celltype={ct} | n={int(mask.sum()):,} cells", fontsize=10)
            ax.legend(loc="best", fontsize=8, frameon=False)
            safe_savefig(fig, plot_root / "per_cell" / "ridges" / _safe(sx)
                         / f"{_safe(ct)}_top_pathways.png")


# ---------------------------------------------------------------------------
# Per-cell pathway-localization grid (one figure/tissue; validation panel)
# ---------------------------------------------------------------------------

def plot_per_cell_localization_grid(pc, gsea_df, plot_root, max_cells=100000):
    """One figure per tissue: a grid of mini-UMAPs, one panel per major cell type,
    coloured by THAT cell type's top relevant pathway (lowest FDR) AUCell score.
    Sanity check that pathway scores localize to sensible cell types — every
    cell type's top pathway should light up that cell type's island. Combined
    sex, groups pooled, cells subsampled for speed."""
    if "X_umap" not in pc.obsm:
        print("  [skip localization grid] missing X_umap"); return
    ct_col = _pc_celltype_col(pc)
    if ct_col is None:
        print("  [skip localization grid] missing celltype col"); return

    # subsample once (deterministic) for all panels
    n = pc.n_obs
    if n > max_cells:
        rng = np.random.default_rng(0)
        idx = np.sort(rng.choice(n, max_cells, replace=False))
    else:
        idx = np.arange(n)
    umap = pc.obsm["X_umap"][idx]
    ct_sub = pc.obs[ct_col].to_numpy()[idx]

    # per cell type: top relevant pathway present in the scored h5ad
    picks = []  # (celltype, pathway, fell_back)
    for ct in [c for c in pd.unique(ct_sub) if pd.notna(c)]:
        sub = gsea_df[(gsea_df["celltype"] == ct) & (gsea_df["FDR"] < SIG_FDR)]
        if sub.empty:
            continue
        ranked = sub.sort_values("FDR")["source"].tolist()
        present = [p for p in ranked if p in pc.var_names]
        if present:
            picks.append((ct, present[0], present[0] != ranked[0]))
    if not picks:
        print("  [skip localization grid] no cell type has a sig pathway present in h5ad")
        return

    ncols = min(5, len(picks))
    nrows = int(np.ceil(len(picks) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.4 * ncols, 3.4 * nrows),
                             constrained_layout=True, squeeze=False)
    for k, (ct, path, fell_back) in enumerate(picks):
        r, c = divmod(k, ncols)
        ax = axes[r][c]
        arr = pc[:, path].X
        arr = arr.toarray().ravel() if hasattr(arr, "toarray") else np.asarray(arr).ravel()
        arr = arr[idx]
        vmin = float(np.percentile(arr, 1)); vmax = float(np.percentile(arr, 99))
        if vmax <= vmin:
            vmax = vmin + 1e-3
        sc = ax.scatter(umap[:, 0], umap[:, 1], c=arr, cmap="viridis",
                        vmin=vmin, vmax=vmax, s=2, alpha=0.55, linewidth=0)
        # outline this cell type's cells so you can see the bright region coincide
        m = (ct_sub == ct)
        ax.scatter(umap[m, 0], umap[m, 1], facecolors="none", edgecolors="black",
                   s=5, linewidth=0.18, alpha=0.6)
        ax.set_xticks([]); ax.set_yticks([])
        tag = " (fallback)" if fell_back else ""
        ax.set_title(f"{truncate(str(ct), 24)}{tag}\n{disp(path, 36)}", fontsize=7)
        fig.colorbar(sc, ax=ax, shrink=0.6)
    for k in range(len(picks), nrows * ncols):
        r, c = divmod(k, ncols)
        axes[r][c].set_axis_off()
    fig.suptitle("Pathway-score localization — each panel = that cell type's top pathway "
                 f"(black outline = the cell type; n≤{max_cells:,} cells)", fontsize=11)
    safe_savefig(fig, plot_root / "per_cell" / "pathway_localization_grid.png")


# ---------------------------------------------------------------------------
# Pathway × age trajectory (sex columns; within_group_across_age, brain-only)
# ---------------------------------------------------------------------------

def plot_pathway_age_trajectory(gsea_df, plot_root, top_n):
    wa = gsea_df[gsea_df["cfam"] == WAA].copy()
    if wa.empty:
        return
    wa["pair_s"] = wa["pair"].map(lambda p: "_".join(parse_pair(p)) or "na")
    sexes_all = [s for s in SEX_ORDER if s in wa["sex"].unique()]
    for (lv, ct), block in wa.groupby(["level", "celltype"], observed=True):
        if not (block["FDR"] < SIG_FDR).any():
            continue
        chosen = (block[block["FDR"] < SIG_FDR].groupby("source")["FDR"].min()
                  .sort_values().head(top_n).index.tolist())
        if not chosen:
            continue
        nrow, ncol = len(chosen), len(sexes_all)
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 2.2 * nrow),
                                 constrained_layout=True, squeeze=False)
        for i, path in enumerate(chosen):
            for j, sx in enumerate(sexes_all):
                ax = axes[i][j]
                sp = block[(block["source"] == path) & (block["sex"] == sx)]
                pair_order = sorted(sp["pair_s"].unique())
                if not pair_order:
                    ax.set_axis_off()
                else:
                    for g in order_groups(sp["group_level"].unique()):
                        gg = sp[sp["group_level"] == g].set_index("pair_s").reindex(pair_order)
                        ax.plot(range(len(pair_order)), gg["NES"].values, marker="o",
                                linewidth=1.3, color=GROUP_COLORS.get(g, "black"),
                                label=g.replace("_Stress", ""))
                        sm = (gg["FDR"] < SIG_FDR).fillna(False).to_numpy()
                        if sm.any():
                            xi = np.where(sm)[0]
                            ax.scatter(xi, gg["NES"].values[xi], s=50, facecolors="none",
                                       edgecolors=GROUP_COLORS.get(g, "black"), linewidth=1.1)
                    ax.axhline(0, color="gray", linestyle="--", linewidth=0.5)
                    ax.set_xticks(range(len(pair_order)))
                    ax.set_xticklabels(pair_order, rotation=30, ha="right", fontsize=6)
                    ax.tick_params(axis="y", labelsize=6)
                if i == 0:
                    ax.set_title(f"sex={sx}", fontsize=9)
                if j == 0:
                    ax.set_ylabel(disp(path, 28), fontsize=6.5)
                if i == 0 and j == ncol - 1 and pair_order:
                    ax.legend(fontsize=6, frameon=False, loc="best")
        fig.suptitle(f"Pathway × age trajectory · level={lv} | ct={ct} "
                     f"(open circle = FDR<{SIG_FDR})", fontsize=11)
        safe_savefig(fig, plot_root / "trajectory" / _safe(lv) / _safe(ct)
                     / "pathway_x_age_by_sex.png")


# ---------------------------------------------------------------------------
# Leading-edge heatmaps
# ---------------------------------------------------------------------------

def plot_leading_edge(le_df, gsea_df, plot_root, top_n_paths, top_n_genes, scope_sex):
    head = gsea_df[(gsea_df["sex"] == scope_sex)
                   & (gsea_df["cfam"].isin(PRIMARY_CONTRASTS))
                   & (gsea_df["level"] == "whole") & (gsea_df["FDR"] < SIG_FDR)]
    if head.empty:
        return
    chosen = head["source"].value_counts().head(top_n_paths).index.tolist()
    le = le_df[(le_df["sex"] == scope_sex)
               & (le_df["cfam"].isin(PRIMARY_CONTRASTS))
               & (le_df["level"] == "whole") & (le_df["pathway"].isin(chosen))]
    if le.empty:
        return
    for path in chosen:
        sub = le[le["pathway"] == path]
        if sub.empty:
            continue
        genes = (sub.groupby("gene")["log2FC"].agg(lambda s: s.abs().max())
                    .sort_values(ascending=False).head(top_n_genes).index.tolist())
        sub_g = sub[sub["gene"].isin(genes)]
        contrasts = [c for c in [CONTRAST_EVR, CONTRAST_LVR] if c in sub_g["cfam"].unique()]
        if not contrasts:
            continue
        fig, axes = plt.subplots(1, len(contrasts),
                                 figsize=(6 + 3 * (len(contrasts) - 1), max(5, 0.3 * len(genes) + 2)),
                                 constrained_layout=True, sharey=True, squeeze=False)
        axes = axes[0]
        vmax = safe_vmax(sub_g["log2FC"].to_numpy())
        for ax, c in zip(axes, contrasts):
            piv = (sub_g[sub_g["cfam"] == c]
                   .pivot_table(index="gene", columns="celltype", values="log2FC", aggfunc="mean")
                   .reindex(index=genes))
            im = ax.imshow(piv.values, cmap="RdBu_r",
                           norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), aspect="auto")
            ax.set_xticks(range(piv.shape[1]))
            ax.set_xticklabels(piv.columns, rotation=60, ha="right", fontsize=7)
            ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=7)
            ax.set_title(CONTRAST_LABEL.get(c, c), fontsize=9)
            fig.colorbar(im, ax=ax, shrink=0.6, label="log₂FC")
        fig.suptitle(f"{disp(path, 70)} · leading-edge × celltype · sex={scope_sex}",
                     fontsize=11)
        safe_savefig(fig, plot_root / "leading_edge" / _safe(scope_sex)
                     / f"{_safe(path)}_genes_x_celltype.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--subcluster", default=None)
    ap.add_argument("--sex-strata", default="combined",
                    help="Sexes for per-slice plots (default combined). "
                         "Concordance/UMAP/trajectory matrices always use all sexes.")
    ap.add_argument("--top-n", type=int, default=25)
    ap.add_argument("--max-volcano-labels", type=int, default=25)
    ap.add_argument("--skip-volcano", action="store_true")
    ap.add_argument("--skip-dotplot", action="store_true")
    ap.add_argument("--skip-cross-celltype", action="store_true")
    ap.add_argument("--skip-concordance", action="store_true")
    ap.add_argument("--skip-per-cell", action="store_true")
    ap.add_argument("--skip-trajectory", action="store_true")
    ap.add_argument("--skip-leading-edge", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config)
    tissue = cfg.get("tissue")
    rdir = Path(cfg["results_dir"])
    cfg_dir = args.config.parent
    suffix = f"_subcluster_{args.subcluster}" if args.subcluster else ""

    tbl = rdir / "tables" / f"08c_pathways{suffix}"
    gsea_path = tbl / f"08c_pathway_results{suffix}.csv"
    tf_path = tbl / f"08c_tf_activity{suffix}.csv"
    le_path = tbl / f"08c_pathway_leading_edge{suffix}.csv"
    pc_path = rdir / "h5ad" / "08c_pathway_scores" / f"{tissue}{suffix}_per_cell_scores.h5ad"
    plot_root = rdir / "plots" / f"08c_pathways{suffix}"
    plot_root.mkdir(parents=True, exist_ok=True)

    if not gsea_path.is_file():
        sys.exit(f"ERROR: {gsea_path} not found. Run 08c_pathways.py first.")
    print(f"\n=== 08c summary plots ===  tissue={tissue} subcluster={args.subcluster or '-'}")
    print(f"  loading {gsea_path}  ({gsea_path.stat().st_size/1e6:.0f} MB)")

    gsea = pd.read_csv(gsea_path, low_memory=False)   # ALL sexes
    gsea = gsea[gsea["collection"].isin(COLLECTIONS)].copy()   # drop M8 etc. from EVERY plot
    gsea["cfam"] = gsea["contrast"].map(contrast_family)
    print(f"  {len(gsea):,} GSEA rows (collections {COLLECTIONS})")
    print(f"  contrast families: {gsea['cfam'].value_counts().to_dict()}")
    go_map, go_union = load_go_relevance(tissue, cfg_dir)

    sex_strata = [s.strip() for s in args.sex_strata.split(",")]
    gsea_ps = gsea[gsea["sex"].isin(sex_strata)].copy()
    print(f"  per-slice sexes {sex_strata}: {len(gsea_ps):,} rows")

    tf = pd.DataFrame()
    if tf_path.is_file():
        tf = pd.read_csv(tf_path, low_memory=False)
        tf["cfam"] = tf["contrast"].map(contrast_family)
        print(f"  {len(tf):,} TF rows (all sexes)")
    tf_ps = tf[tf["sex"].isin(sex_strata)].copy() if not tf.empty else tf

    # ---- Volcanos (per-slice sexes) ----
    if not args.skip_volcano:
        n = 0
        for keys, block in gsea_ps.groupby(SLICE_ID, observed=True, dropna=False):
            sx, lv, ct = keys
            if not (block["FDR"] < SIG_FDR).any():
                continue
            d = slice_dir(plot_root, sx, lv, ct)
            for contrast in SOLO_CONTRASTS:
                plot_solo_volcano_pathways(block, contrast,
                                           d / f"volcano_{CONTRAST_SLUG[contrast]}.png",
                                           args.max_volcano_labels, go_map)
            plot_dual_volcano_pathways(block, d / "volcano_DUAL_early_late.png",
                                       args.max_volcano_labels, go_map)
            if not tf_ps.empty:
                tfb = tf_ps[(tf_ps["sex"] == sx) & (tf_ps["level"] == lv)
                            & (tf_ps["celltype"] == ct)]
                if not tfb.empty:
                    for contrast in SOLO_CONTRASTS:
                        plot_solo_volcano_tf(tfb, contrast,
                                             d / f"tf_volcano_{CONTRAST_SLUG[contrast]}.png",
                                             args.max_volcano_labels)
                    plot_dual_volcano_tf(tfb, d / "tf_volcano_DUAL_early_late.png",
                                         args.max_volcano_labels)
            n += 1
            if n % 25 == 0:
                print(f"    volcano slices: {n} ...")
        print(f"  volcano slices plotted: {n}")

    # ---- Dotplots (per-slice sexes) ----
    if not args.skip_dotplot:
        n = 0
        for keys, block in gsea_ps.groupby(SLICE_ID, observed=True, dropna=False):
            sx, lv, ct = keys
            d = slice_dir(plot_root, sx, lv, ct)
            for contrast in SOLO_CONTRASTS:
                plot_dotplot_panels(block, contrast,
                                    d / f"dotplot_{CONTRAST_SLUG[contrast]}.png",
                                    args.top_n, go_map)
            n += 1
        print(f"  dotplot slices attempted: {n}")

    # ---- Cross-celltype heatmaps (per sex×level×contrast) ----
    if not args.skip_cross_celltype:
        n = 0
        for keys, block in gsea_ps.groupby(["sex", "level", "cfam"],
                                           observed=True, dropna=False):
            sx, lv, contrast = keys
            if contrast not in SOLO_CONTRASTS:
                continue
            if not (block["FDR"] < SIG_FDR).any():
                continue
            d = slice_dir(plot_root, sx, lv) / "_celltype_heatmaps"
            plot_celltype_pathway_heatmap(
                block, d / f"celltype_x_pathway_{CONTRAST_SLUG[contrast]}.png",
                args.top_n, go_union)
            if not tf_ps.empty:
                tfb = tf_ps[(tf_ps["sex"] == sx) & (tf_ps["level"] == lv)
                            & (tf_ps["cfam"] == contrast)]
                if not tfb.empty:
                    plot_tf_celltype_heatmap(
                        tfb, d / f"tf_x_celltype_{CONTRAST_SLUG[contrast]}.png", args.top_n)
            n += 1
        print(f"  cross-celltype heatmap files: {n}")

    # ---- Concordance: sex×age matrices (ALL sexes) ----
    if not args.skip_concordance:
        nc = nct = 0
        prim = gsea[gsea["cfam"].isin(PRIMARY_CONTRASTS)]
        for (lv, ct), block in prim.groupby(["level", "celltype"], observed=True):
            e = block[block["cfam"] == CONTRAST_EVR]
            l = block[block["cfam"] == CONTRAST_LVR]
            if e.empty or l.empty:
                continue
            cc = plot_root / "concordance" / _safe(lv) / _safe(ct)
            plot_concordance_matrix(
                e, l, "source", "NES", "FDR",
                cc / "concordance_pathway_sex_x_age.png", args.top_n,
                f"Pathway NES concordance (Early vs Late) · level={lv} | ct={ct}",
                prefix_strip=True)
            nc += 1
            if not tf.empty:
                te = tf[(tf["cfam"] == CONTRAST_EVR) & (tf["level"] == lv)
                        & (tf["celltype"] == ct)]
                tl = tf[(tf["cfam"] == CONTRAST_LVR) & (tf["level"] == lv)
                        & (tf["celltype"] == ct)]
                if not te.empty and not tl.empty:
                    fc = _tf_fdr_col(te)
                    plot_concordance_matrix(
                        te, tl, "TF", "activity_score", fc,
                        cc / "concordance_tf_sex_x_age.png", args.top_n,
                        f"TF activity concordance (Early vs Late) · level={lv} | ct={ct}",
                        prefix_strip=False)
                    nct += 1
        print(f"  concordance matrices: {nc} pathway, {nct} TF")

    # ---- Per-cell ----
    if not args.skip_per_cell:
        if pc_path.is_file():
            print(f"  loading per-cell {pc_path}")
            pc = ad.read_h5ad(pc_path)
            plot_per_cell_ridges(pc, gsea, plot_root, top_n=20)
            plot_per_cell_localization_grid(pc, gsea, plot_root, max_cells=100000)
            print("  per-cell ridges + localization grid done")
        else:
            print(f"  [skip per-cell] {pc_path} not found")

    # ---- Trajectory (brain-only) ----
    if not args.skip_trajectory:
        plot_pathway_age_trajectory(gsea, plot_root, top_n=12)
        print("  pathway×age trajectories done (within_group_across_age only)")

    # ---- Leading-edge ----
    if not args.skip_leading_edge:
        if le_path.is_file():
            print(f"  loading leading-edge {le_path}  ({le_path.stat().st_size/1e6:.0f} MB)")
            le = pd.read_csv(le_path, low_memory=False)
            le["cfam"] = le["contrast"].map(contrast_family)
            scope = sex_strata[0]
            plot_leading_edge(le, gsea, plot_root, top_n_paths=50, top_n_genes=30,
                              scope_sex=scope)
            print(f"  leading-edge heatmaps done (sex={scope})")
        else:
            print(f"  [skip leading-edge] {le_path} not found")

    print(f"\n✓ Phase 8c summary complete. Plots under {plot_root}")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    main()
