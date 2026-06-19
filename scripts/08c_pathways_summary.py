#!/usr/bin/env python
"""
08c_pathways_summary.py — Phase 8c summary plots (CSV-only inputs, no recompute).

Mirrors the 8b / 8b_summary split: 08c_pathways.py writes CSVs; this file plots.

VOLCANOS (headline; pathways AND TFs):
  Per (sex × group_level × level × celltype), four volcanos each:
    1. solo early_vs_relaxed_per_age
    2. solo late_vs_relaxed_per_age
    3. solo early_vs_late_per_age
    4. DUAL overlay E-v-R + L-v-R: left Y=-log10(FDR_early), right Y=-log10(FDR_late);
       bottom X=NES_early, top X=NES_late (NES axes SHARED range). Labels ONLY
       pathways/TFs sig in BOTH (FDR<0.05 each) AND same NES/activity sign.

  Y-axis: -log10(FDR) capped at 20; off-scale (FDR<1e-20) drawn as ▲ at y=20.

OTHER PLOTS (iterate all sex strata):
  dotplot panels, celltype×pathway NES heatmap, TF×celltype heatmap,
  bubble across contrasts, pathway + TF concordance scatters,
  per-cell ridges, per-cell UMAPs, leading-edge heatmaps,
  pathway×age trajectory (within_group_across_age, brain-only).

Output: plots/08c_pathways{suffix}/...  (PNG @ 300 DPI; constrained_layout)

Usage:
  uv run python scripts/08c_pathways_summary.py --config config/brain.yaml
  uv run python scripts/08c_pathways_summary.py --config config/brain.yaml --subcluster immune
"""

import argparse
import ast
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
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


SIG_FDR = 0.05
FDR_FLOOR = 1e-20
YCAP = -np.log10(FDR_FLOOR)
MAX_LABEL_CHARS = 42
COLLECTIONS = ["MH", "M2", "M5", "M8"]

CONTRAST_EVR = "early_vs_relaxed_per_age"
CONTRAST_LVR = "late_vs_relaxed_per_age"
CONTRAST_EVL = "early_vs_late_per_age"
SOLO_CONTRASTS = [CONTRAST_EVR, CONTRAST_LVR, CONTRAST_EVL]
PRIMARY_CONTRASTS = {CONTRAST_EVR, CONTRAST_LVR}

SLICE_ID = ["sex", "group_level", "level", "celltype"]

GROUP_COLORS = {"Relaxed": "#7f7f7f", "Early_Stress": "#d62728", "Late_Stress": "#1f77b4"}
EARLY_C = "#d62728"
LATE_C = "#1f77b4"
DIRECTION_COLORS = {"up": "#d62728", "down": "#1f77b4"}
CONTRAST_LABEL = {CONTRAST_EVR: "Early vs Relaxed",
                  CONTRAST_LVR: "Late vs Relaxed",
                  CONTRAST_EVL: "Early vs Late"}


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


def truncate(s, n=MAX_LABEL_CHARS):
    s = str(s)
    return s if len(s) <= n else s[: n - 1] + "…"


def safe_vmax(arr, floor=1.0):
    """Max |value| with NaN-safety; never returns nan/0 (TwoSlopeNorm needs vmax>0)."""
    a = np.abs(np.asarray(arr, dtype=float))
    if not np.isfinite(a).any():
        return floor
    v = np.nanmax(a)
    return max(float(v), floor) if np.isfinite(v) else floor


def neglog10_fdr(fdr_series):
    f = pd.to_numeric(fdr_series, errors="coerce").to_numpy(dtype=float)
    offscale = f < FDR_FLOOR
    f = np.clip(f, FDR_FLOOR, 1.0)
    y = -np.log10(f)
    y = np.nan_to_num(y, nan=0.0)
    return y, offscale


def safe_savefig(fig, path: Path, dpi=300):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def add_labels(ax, xs, ys, names, fontsize=6.0, bold=False):
    texts = []
    for x, y, nm in zip(xs, ys, names):
        texts.append(ax.text(x, y, truncate(nm), fontsize=fontsize, ha="left",
                             va="bottom", fontweight="bold" if bold else "normal"))
    if HAVE_ADJUSTTEXT and len(texts) > 1:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.4, alpha=0.7),
                    expand=(1.3, 1.4), force_text=(0.4, 0.5))
    return texts


def slice_dir(plot_root, sex, group_level, level, celltype=None):
    d = plot_root / _safe(sex) / _safe(group_level) / _safe(level)
    if celltype is not None:
        d = d / _safe(celltype)
    return d


def slice_title(sex, group_level, level, celltype=None, extra=""):
    parts = [f"sex={sex}", f"age={group_level}", f"level={level}"]
    if celltype is not None:
        parts.append(f"ct={celltype}")
    s = " | ".join(parts)
    return f"{s}{('  ·  ' + extra) if extra else ''}"


def _scatter_volcano(ax, df, x_col, fdr_col, name_col, color_up, color_down,
                     label_all_sig, max_labels, label_bold=False):
    if df.empty:
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes, color="gray")
        return 0
    y, off = neglog10_fdr(df[fdr_col])
    x = pd.to_numeric(df[x_col], errors="coerce").to_numpy(dtype=float)
    sig = (pd.to_numeric(df[fdr_col], errors="coerce") < SIG_FDR).to_numpy()
    ax.scatter(x[~sig], y[~sig], s=7, c="lightgray", alpha=0.4, linewidth=0)
    up = sig & (x > 0)
    dn = sig & (x < 0)
    ax.scatter(x[up], y[up], s=24, c=color_up, alpha=0.85, edgecolors="black", linewidth=0.25)
    ax.scatter(x[dn], y[dn], s=24, c=color_down, alpha=0.85, edgecolors="black", linewidth=0.25)
    osig = sig & off
    if osig.any():
        ax.scatter(x[osig], np.full(osig.sum(), YCAP), marker="^", s=40,
                   c=[color_up if xx > 0 else color_down for xx in x[osig]],
                   edgecolors="black", linewidth=0.3, zorder=5)
    ax.axhline(-np.log10(SIG_FDR), color="gray", linestyle="--", linewidth=0.5)
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.5)
    n_sig = int(sig.sum())
    if label_all_sig and n_sig:
        names = df[name_col].to_numpy()
        sig_idx = np.where(sig)[0]
        fdr_vals = pd.to_numeric(df[fdr_col], errors="coerce").to_numpy()
        order = sig_idx[np.argsort(fdr_vals[sig_idx])][:max_labels]
        add_labels(ax, x[order], y[order], names[order], fontsize=5.8, bold=label_bold)
    return n_sig


def plot_solo_volcano_pathways(slice_df, contrast, outpath, max_labels):
    df = slice_df[slice_df["contrast"] == contrast]
    if df.empty:
        return
    fig, axes = plt.subplots(1, 4, figsize=(22, 5.6), constrained_layout=True)
    for ax, coll in zip(axes, COLLECTIONS):
        d = df[df["collection"] == coll]
        n_sig = _scatter_volcano(ax, d, "NES", "FDR", "source",
                                 DIRECTION_COLORS["up"], DIRECTION_COLORS["down"],
                                 True, max_labels)
        cap = "" if n_sig <= max_labels else f" (top {max_labels} labeled)"
        ax.set_xlabel("NES")
        ax.set_ylabel(r"-log$_{10}$(FDR)")
        ax.set_title(f"{coll}  n={len(d):,}  sig={n_sig}{cap}", fontsize=9)
    r = df.iloc[0]
    fig.suptitle(slice_title(r["sex"], r["group_level"], r["level"], r["celltype"],
                             extra=CONTRAST_LABEL.get(contrast, contrast)), fontsize=10)
    safe_savefig(fig, outpath)


def plot_dual_volcano_pathways(slice_df, outpath, max_labels):
    e = slice_df[slice_df["contrast"] == CONTRAST_EVR]
    l = slice_df[slice_df["contrast"] == CONTRAST_LVR]
    if e.empty or l.empty:
        return
    fig, axes = plt.subplots(1, 4, figsize=(22, 5.8), constrained_layout=True)
    any_panel = False
    for ax, coll in zip(axes, COLLECTIONS):
        ec = e[e["collection"] == coll]
        lc = l[l["collection"] == coll]
        if ec.empty and lc.empty:
            ax.set_axis_off(); ax.set_title(f"no {coll}", fontsize=9); continue
        any_panel = True
        nes_all = pd.concat([ec["NES"], lc["NES"]]).astype(float)
        nmax = safe_vmax(nes_all.values) * 1.1
        ye, offe = neglog10_fdr(ec["FDR"])
        yl, offl = neglog10_fdr(lc["FDR"])
        xe = ec["NES"].astype(float).to_numpy()
        xl = lc["NES"].astype(float).to_numpy()
        ax.scatter(xe, ye, s=16, c=EARLY_C, alpha=0.55, linewidth=0)
        if offe.any():
            ax.scatter(xe[offe], np.full(offe.sum(), YCAP), marker="^", s=30,
                       c=EARLY_C, edgecolors="black", linewidth=0.25, zorder=5)
        ax.set_xlim(-nmax, nmax); ax.set_ylim(0, YCAP * 1.05)
        ax.set_xlabel("NES (Early vs Relaxed)", color=EARLY_C)
        ax.set_ylabel(r"-log$_{10}$(FDR) Early", color=EARLY_C)
        ax.tick_params(axis="x", colors=EARLY_C); ax.tick_params(axis="y", colors=EARLY_C)
        ax_r = ax.twinx(); ax_t = ax_r.twiny()
        ax_t.scatter(xl, yl, s=16, c=LATE_C, alpha=0.55, linewidth=0)
        if offl.any():
            ax_t.scatter(xl[offl], np.full(offl.sum(), YCAP), marker="^", s=30,
                         c=LATE_C, edgecolors="black", linewidth=0.25, zorder=5)
        ax_r.set_ylim(0, YCAP * 1.05); ax_t.set_xlim(-nmax, nmax)
        ax_r.set_ylabel(r"-log$_{10}$(FDR) Late", color=LATE_C)
        ax_t.set_xlabel("NES (Late vs Relaxed)", color=LATE_C)
        ax_r.tick_params(axis="y", colors=LATE_C); ax_t.tick_params(axis="x", colors=LATE_C)
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.4)
        ax.axhline(-np.log10(SIG_FDR), color="gray", linestyle=":", linewidth=0.4)
        m = pd.merge(
            ec[["source", "NES", "FDR"]].rename(columns={"NES": "NES_e", "FDR": "FDR_e"}),
            lc[["source", "NES", "FDR"]].rename(columns={"NES": "NES_l", "FDR": "FDR_l"}),
            on="source", how="inner")
        m = m[(m["FDR_e"] < SIG_FDR) & (m["FDR_l"] < SIG_FDR)
              & (np.sign(m["NES_e"]) == np.sign(m["NES_l"]))]
        n_conc = len(m)
        if n_conc:
            m = m.reindex(m[["FDR_e", "FDR_l"]].max(axis=1).sort_values().index).head(max_labels)
            ye_l, _ = neglog10_fdr(m["FDR_e"])
            texts = []
            for xx, yy, nm in zip(m["NES_e"].astype(float), ye_l, m["source"]):
                texts.append(ax.text(xx, yy, truncate(nm), fontsize=5.8, ha="left",
                                     va="bottom", fontweight="bold", color="black"))
            if HAVE_ADJUSTTEXT and len(texts) > 1:
                adjust_text(texts, ax=ax,
                            arrowprops=dict(arrowstyle="-", color="gray", lw=0.4),
                            expand=(1.3, 1.4))
        ax.set_title(f"{coll}  concordant-sig={n_conc}", fontsize=9)
    if not any_panel:
        plt.close(fig); return
    r = e.iloc[0]
    fig.suptitle(slice_title(r["sex"], r["group_level"], r["level"], r["celltype"],
                             extra="DUAL · shared NES range · labels = sig in BOTH + same NES sign"),
                 fontsize=10)
    safe_savefig(fig, outpath)


def _tf_fdr_col(df):
    return "FDR_ctx_celltype" if "FDR_ctx_celltype" in df.columns else "FDR"


def plot_solo_volcano_tf(tf_slice, contrast, outpath, max_labels):
    df = tf_slice[tf_slice["contrast"] == contrast]
    if df.empty:
        return
    fdr_col = _tf_fdr_col(df)
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    n_sig = _scatter_volcano(ax, df, "activity_score", fdr_col, "TF",
                             DIRECTION_COLORS["up"], DIRECTION_COLORS["down"],
                             True, max_labels, label_bold=True)
    ax.set_xlabel("TF activity score (ULM t-stat)")
    ax.set_ylabel(rf"-log$_{{10}}$({fdr_col})")
    cap = "" if n_sig <= max_labels else f" (top {max_labels} labeled)"
    r = df.iloc[0]
    ax.set_title(f"TF activity · {CONTRAST_LABEL.get(contrast, contrast)} · "
                 f"sig({fdr_col}<{SIG_FDR})={n_sig}{cap}", fontsize=10)
    fig.suptitle(slice_title(r["sex"], r["group_level"], r["level"], r["celltype"]), fontsize=10)
    safe_savefig(fig, outpath)


def plot_dual_volcano_tf(tf_slice, outpath, max_labels):
    e = tf_slice[tf_slice["contrast"] == CONTRAST_EVR]
    l = tf_slice[tf_slice["contrast"] == CONTRAST_LVR]
    if e.empty or l.empty:
        return
    fdr_col = _tf_fdr_col(e)
    fig, ax = plt.subplots(figsize=(8.5, 6.5), constrained_layout=True)
    act_all = pd.concat([e["activity_score"], l["activity_score"]]).astype(float)
    amax = safe_vmax(act_all.values) * 1.1
    ye, offe = neglog10_fdr(e[fdr_col])
    yl, offl = neglog10_fdr(l[fdr_col])
    xe = e["activity_score"].astype(float).to_numpy()
    xl = l["activity_score"].astype(float).to_numpy()
    ax.scatter(xe, ye, s=22, c=EARLY_C, alpha=0.55, linewidth=0)
    if offe.any():
        ax.scatter(xe[offe], np.full(offe.sum(), YCAP), marker="^", s=36,
                   c=EARLY_C, edgecolors="black", linewidth=0.25, zorder=5)
    ax.set_xlim(-amax, amax); ax.set_ylim(0, YCAP * 1.05)
    ax.set_xlabel("activity (Early vs Relaxed)", color=EARLY_C)
    ax.set_ylabel(rf"-log$_{{10}}$({fdr_col}) Early", color=EARLY_C)
    ax.tick_params(axis="x", colors=EARLY_C); ax.tick_params(axis="y", colors=EARLY_C)
    ax_r = ax.twinx(); ax_t = ax_r.twiny()
    ax_t.scatter(xl, yl, s=22, c=LATE_C, alpha=0.55, linewidth=0)
    if offl.any():
        ax_t.scatter(xl[offl], np.full(offl.sum(), YCAP), marker="^", s=36,
                     c=LATE_C, edgecolors="black", linewidth=0.25, zorder=5)
    ax_r.set_ylim(0, YCAP * 1.05); ax_t.set_xlim(-amax, amax)
    ax_r.set_ylabel(rf"-log$_{{10}}$({fdr_col}) Late", color=LATE_C)
    ax_t.set_xlabel("activity (Late vs Relaxed)", color=LATE_C)
    ax_r.tick_params(axis="y", colors=LATE_C); ax_t.tick_params(axis="x", colors=LATE_C)
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.4)
    ax.axhline(-np.log10(SIG_FDR), color="gray", linestyle=":", linewidth=0.4)
    m = pd.merge(
        e[["TF", "activity_score", fdr_col]].rename(columns={"activity_score": "a_e", fdr_col: "f_e"}),
        l[["TF", "activity_score", fdr_col]].rename(columns={"activity_score": "a_l", fdr_col: "f_l"}),
        on="TF", how="inner")
    m = m[(m["f_e"] < SIG_FDR) & (m["f_l"] < SIG_FDR)
          & (np.sign(m["a_e"]) == np.sign(m["a_l"]))]
    n_conc = len(m)
    if n_conc:
        m = m.reindex(m[["f_e", "f_l"]].max(axis=1).sort_values().index).head(max_labels)
        ye_l, _ = neglog10_fdr(m["f_e"])
        texts = []
        for xx, yy, nm in zip(m["a_e"].astype(float), ye_l, m["TF"]):
            texts.append(ax.text(xx, yy, str(nm), fontsize=6.2, ha="left",
                                 va="bottom", fontweight="bold"))
        if HAVE_ADJUSTTEXT and len(texts) > 1:
            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle="-", color="gray", lw=0.4),
                        expand=(1.3, 1.4))
    legend_el = [Line2D([0], [0], marker="o", color="w", markerfacecolor=EARLY_C,
                        markersize=8, label="Early vs Rel"),
                 Line2D([0], [0], marker="o", color="w", markerfacecolor=LATE_C,
                        markersize=8, label="Late vs Rel")]
    ax.legend(handles=legend_el, loc="upper center", fontsize=8, frameon=False, ncol=2)
    r = e.iloc[0]
    ax.set_title(f"TF DUAL · concordant-sig={n_conc} · labels = sig in BOTH + same sign",
                 fontsize=10)
    fig.suptitle(slice_title(r["sex"], r["group_level"], r["level"], r["celltype"]), fontsize=10)
    safe_savefig(fig, outpath)


def plot_dotplot_panels(slice_df, contrast, outpath, top_n):
    df = slice_df[slice_df["contrast"] == contrast]
    if df.empty:
        return
    fig, axes = plt.subplots(1, 4, figsize=(22, max(4.5, 0.25 * top_n)), constrained_layout=True)
    used = False
    for ax, coll in zip(axes, COLLECTIONS):
        d = df[(df["collection"] == coll) & (df["FDR"] < SIG_FDR)].copy()
        if d.empty:
            ax.text(0.5, 0.5, f"no sig {coll}", ha="center", va="center",
                    transform=ax.transAxes, color="gray"); ax.set_axis_off(); continue
        used = True
        top = d.sort_values("FDR").head(top_n).iloc[::-1]
        top["nlogq"], _ = neglog10_fdr(top["FDR"])
        y = np.arange(len(top))
        colors = [DIRECTION_COLORS["up"] if n > 0 else DIRECTION_COLORS["down"] for n in top["NES"]]
        sizes = 40 + 26 * top["nlogq"].clip(upper=YCAP)
        ax.scatter(top["NES"], y, s=sizes, c=colors, edgecolors="black", linewidth=0.3, alpha=0.9)
        ax.set_yticks(y); ax.set_yticklabels([truncate(s) for s in top["source"]], fontsize=6.5)
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.5)
        ax.set_xlabel("NES"); ax.set_title(f"{coll}  top {len(top)} sig", fontsize=9)
        ax.grid(axis="x", linestyle=":", linewidth=0.4, alpha=0.6)
    if not used:
        plt.close(fig); return
    r = df.iloc[0]
    fig.suptitle(slice_title(r["sex"], r["group_level"], r["level"], r["celltype"],
                             extra=CONTRAST_LABEL.get(contrast, contrast)), fontsize=10)
    safe_savefig(fig, outpath)


def plot_celltype_pathway_heatmap(block, outpath, top_n):
    fig, axes = plt.subplots(1, 4, figsize=(24, 6), constrained_layout=True)
    used = False
    for ax, coll in zip(axes, COLLECTIONS):
        d = block[block["collection"] == coll]
        sig = d[d["FDR"] < SIG_FDR]
        if sig.empty:
            ax.text(0.5, 0.5, f"no sig {coll}", ha="center", va="center",
                    transform=ax.transAxes, color="gray"); ax.set_axis_off(); continue
        used = True
        chosen = sig.groupby("source")["FDR"].min().sort_values().head(top_n).index.tolist()
        sub = d[d["source"].isin(chosen)]
        nes = sub.pivot_table(index="celltype", columns="source", values="NES", aggfunc="first")
        fdr = sub.pivot_table(index="celltype", columns="source", values="FDR", aggfunc="first")
        rows = nes.abs().mean(axis=1).sort_values(ascending=False).index
        nes = nes.reindex(index=rows, columns=chosen); fdr = fdr.reindex(index=rows, columns=chosen)
        vmax = safe_vmax(nes.values)
        im = ax.imshow(nes.values, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), aspect="auto")
        for (i, j), f in np.ndenumerate(fdr.values):
            if not np.isnan(f) and f < SIG_FDR:
                ax.add_patch(mpatches.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                                edgecolor="black", linewidth=1.2))
        ax.set_xticks(range(len(chosen)))
        ax.set_xticklabels([truncate(c, 35) for c in chosen], rotation=70, ha="right", fontsize=6.5)
        ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=7.5)
        ax.set_title(f"{coll}  ({len(chosen)}p × {len(rows)}ct; black=FDR<{SIG_FDR})", fontsize=9)
        fig.colorbar(im, ax=ax, shrink=0.6, label="NES")
    if not used:
        plt.close(fig); return
    r = block.iloc[0]
    fig.suptitle(slice_title(r["sex"], r["group_level"], r["level"],
                             extra=CONTRAST_LABEL.get(r["contrast"], r["contrast"])), fontsize=10)
    safe_savefig(fig, outpath)


def plot_tf_celltype_heatmap(tf_block, outpath, top_n):
    df = tf_block
    if df.empty:
        return
    fdr_col = _tf_fdr_col(df)
    sig = df[df[fdr_col] < SIG_FDR]
    if sig.empty:
        return
    chosen = sig.groupby("TF")[fdr_col].min().sort_values().head(top_n).index.tolist()
    sub = df[df["TF"].isin(chosen)]
    act = sub.pivot_table(index="celltype", columns="TF", values="activity_score", aggfunc="first")
    fdr = sub.pivot_table(index="celltype", columns="TF", values=fdr_col, aggfunc="first")
    rows = act.abs().mean(axis=1).sort_values(ascending=False).index
    act = act.reindex(index=rows, columns=chosen); fdr = fdr.reindex(index=rows, columns=chosen)
    vmax = safe_vmax(act.values)
    fig, ax = plt.subplots(figsize=(max(8, 0.35 * len(chosen) + 4),
                                     max(4, 0.35 * len(rows) + 2)), constrained_layout=True)
    im = ax.imshow(act.values, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), aspect="auto")
    for (i, j), f in np.ndenumerate(fdr.values):
        if not np.isnan(f) and f < SIG_FDR:
            ax.add_patch(mpatches.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                            edgecolor="black", linewidth=1.2))
    ax.set_xticks(range(len(chosen))); ax.set_xticklabels(chosen, rotation=70, ha="right", fontsize=7)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=8)
    ax.set_title(f"TF×celltype top {len(chosen)} (black={fdr_col}<{SIG_FDR})", fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.7, label="activity")
    r = df.iloc[0]
    fig.suptitle(slice_title(r["sex"], r["group_level"], r["level"],
                             extra=CONTRAST_LABEL.get(r["contrast"], r["contrast"])), fontsize=10)
    safe_savefig(fig, outpath)


def plot_bubble(gsea_sl, outpath, top_n):
    contrasts = [c for c in [CONTRAST_EVR, CONTRAST_LVR] if c in gsea_sl["contrast"].unique()]
    if not contrasts:
        return
    chosen = set()
    for c in contrasts:
        sub = gsea_sl[(gsea_sl["contrast"] == c) & (gsea_sl["FDR"] < SIG_FDR)]
        if sub.empty:
            continue
        chosen |= set(sub.groupby("source")["FDR"].min().sort_values().head(top_n).index)
    if not chosen:
        return
    chosen = list(chosen)
    celltypes = sorted(gsea_sl["celltype"].unique())
    fig, axes = plt.subplots(1, len(contrasts),
                              figsize=(8 + 3 * (len(contrasts) - 1), max(5, 0.3 * len(chosen) + 2)),
                              constrained_layout=True, sharey=True)
    if len(contrasts) == 1:
        axes = [axes]
    for ax, c in zip(axes, contrasts):
        sub = gsea_sl[(gsea_sl["contrast"] == c) & gsea_sl["source"].isin(chosen)].copy()
        if sub.empty:
            ax.set_axis_off(); ax.set_title(f"{CONTRAST_LABEL[c]}\n(none)", fontsize=9); continue
        sub["nlogq"], _ = neglog10_fdr(sub["FDR"])
        ct_idx = {ct: i for i, ct in enumerate(celltypes)}
        p_idx = {p: i for i, p in enumerate(chosen)}
        xs = sub["celltype"].map(ct_idx); ys = sub["source"].map(p_idx)
        sizes = 20 + 24 * sub["nlogq"].clip(upper=8)
        sc = ax.scatter(xs, ys, c=sub["NES"], cmap="RdBu_r", vmin=-2.5, vmax=2.5,
                        s=sizes, edgecolors="black", linewidth=0.3, alpha=0.85)
        sig = sub["FDR"] < SIG_FDR
        if sig.any():
            ax.scatter(xs[sig], ys[sig], facecolors="none", edgecolors="black",
                       linewidth=1.0, s=sizes[sig] + 30)
        ax.set_xticks(range(len(celltypes)))
        ax.set_xticklabels(celltypes, rotation=60, ha="right", fontsize=8)
        ax.set_yticks(range(len(chosen)))
        ax.set_yticklabels([truncate(p) for p in chosen], fontsize=7)
        ax.set_title(CONTRAST_LABEL[c], fontsize=9)
        ax.grid(linestyle=":", linewidth=0.4, alpha=0.5)
        fig.colorbar(sc, ax=ax, shrink=0.6, label="NES")
    s, lv = gsea_sl.iloc[0]["sex"], gsea_sl.iloc[0]["level"]
    fig.suptitle(f"sex={s} | level={lv}  ·  top {len(chosen)} pathways "
                 f"(black outline = FDR<{SIG_FDR})", fontsize=10)
    safe_savefig(fig, outpath)


def plot_concordance(df_a, df_b, key, val, fdr, outpath, max_labels, title_prefix):
    a = df_a[[key, val, fdr]].rename(columns={val: "va", fdr: "fa"})
    b = df_b[[key, val, fdr]].rename(columns={val: "vb", fdr: "fb"})
    m = pd.merge(a, b, on=key, how="inner")
    if m.empty:
        return
    m["sa"] = m["fa"] < SIG_FDR; m["sb"] = m["fb"] < SIG_FDR
    m["state"] = "ns"
    m.loc[m["sa"] & m["sb"], "state"] = "both"
    m.loc[m["sa"] & ~m["sb"], "state"] = "early_only"
    m.loc[~m["sa"] & m["sb"], "state"] = "late_only"
    fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)
    palette = {"ns": "lightgray", "both": "#9467bd", "early_only": EARLY_C, "late_only": LATE_C}
    for st, col in palette.items():
        s = m[m["state"] == st]
        ax.scatter(s["va"], s["vb"], c=col, s=15 if st == "ns" else 28,
                   alpha=0.65 if st == "ns" else 0.85,
                   edgecolors="none" if st == "ns" else "black", linewidth=0.3,
                   label=f"{st} (n={len(s)})")
    lo = float(min(m["va"].min(), m["vb"].min())); hi = float(max(m["va"].max(), m["vb"].max()))
    pad = 0.05 * (hi - lo) if hi > lo else 0.5
    lo, hi = lo - pad, hi + pad
    ax.plot([lo, hi], [lo, hi], color="gray", linestyle="--", linewidth=0.7)
    ax.axhline(0, color="gray", linestyle=":", linewidth=0.4)
    ax.axvline(0, color="gray", linestyle=":", linewidth=0.4)
    pool = m[m["state"] == "both"].copy()
    if not pool.empty:
        pool["ext"] = pool["va"].abs() + pool["vb"].abs()
        top = pool.sort_values("ext", ascending=False).head(max_labels)
        add_labels(ax, top["va"], top["vb"], top[key], fontsize=6.5, bold=True)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel(f"{val}  ·  Early vs Relaxed")
    ax.set_ylabel(f"{val}  ·  Late vs Relaxed")
    ax.legend(loc="best", fontsize=8, frameon=False)
    ax.set_title(f"{title_prefix}  n_shared={len(m):,}", fontsize=9)
    safe_savefig(fig, outpath)


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
            long["pathway"] = long["pathway"].map(lambda p: truncate(p, 38))
            order = [truncate(p, 38) for p in chosen]
            hue_order = [g for g in ["Relaxed", "Early_Stress", "Late_Stress"]
                         if g in long["group"].unique()]
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


def plot_per_cell_umaps(pc, gsea_df, plot_root, top_n, group_col="group"):
    if "X_umap" not in pc.obsm or group_col not in pc.obs.columns:
        print("  [skip UMAPs] missing X_umap/group"); return
    chosen = [p for p in gsea_df[gsea_df["FDR"] < SIG_FDR]["source"].value_counts().index[:top_n]
              if p in pc.var_names]
    if not chosen:
        return
    obs = pc.obs; umap = pc.obsm["X_umap"]
    groups = [g for g in ["Relaxed", "Early_Stress", "Late_Stress"] if g in obs[group_col].unique()]
    sexes = sorted(obs["sex"].dropna().unique()) if "sex" in obs.columns else ["combined"]
    for sx in sexes:
        smask = (obs["sex"] == sx).to_numpy() if "sex" in obs.columns else np.ones(len(obs), bool)
        if int(smask.sum()) < 100:
            continue
        for path in chosen:
            sc_arr = pc[:, path].X
            sc_arr = sc_arr.toarray().ravel() if hasattr(sc_arr, "toarray") else np.asarray(sc_arr).ravel()
            vmin = float(np.percentile(sc_arr[smask], 1)); vmax = float(np.percentile(sc_arr[smask], 99))
            if vmax <= vmin:
                vmax = vmin + 1e-3
            fig, axes = plt.subplots(1, len(groups), figsize=(4 * len(groups), 4.3),
                                      constrained_layout=True, sharex=True, sharey=True)
            if len(groups) == 1:
                axes = [axes]
            sc = None
            for ax, g in zip(axes, groups):
                m = smask & (obs[group_col] == g).to_numpy()
                if int(m.sum()) == 0:
                    ax.set_axis_off(); ax.set_title(f"{g} (n=0)", fontsize=9); continue
                sc = ax.scatter(umap[m, 0], umap[m, 1], c=sc_arr[m], cmap="viridis",
                                vmin=vmin, vmax=vmax, s=2, alpha=0.6, linewidth=0)
                ax.set_title(f"{g} (n={int(m.sum()):,})", fontsize=9)
                ax.set_xticks([]); ax.set_yticks([])
            if sc is not None:
                fig.colorbar(sc, ax=axes, shrink=0.7, label="AUCell")
            fig.suptitle(f"{truncate(path, 60)}  |  sex={sx}", fontsize=10)
            safe_savefig(fig, plot_root / "per_cell" / "umaps" / _safe(sx)
                         / f"{_safe(path)}_by_group.png")


def plot_pathway_age_trajectory(gsea_df, plot_root, top_n):
    wa = gsea_df[gsea_df["contrast"] == "within_group_across_age"].copy()
    if wa.empty:
        return
    wa["pair_s"] = wa["pair"].map(lambda p: "_".join(parse_pair(p)) or "na")
    for keys, block in wa.groupby(["sex", "level", "celltype"], observed=True):
        sx, lv, ct = keys
        sig = block[block["FDR"] < SIG_FDR]
        if sig.empty:
            continue
        chosen = sig.groupby("source")["FDR"].min().sort_values().head(top_n).index.tolist()
        sub = block[block["source"].isin(chosen)]
        pair_order = sorted(sub["pair_s"].unique())
        groups = sorted(sub["group_level"].unique())
        n_cols = 4; n_rows = (len(chosen) + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 2.6 * n_rows),
                                  constrained_layout=True, squeeze=False)
        for k, path in enumerate(chosen):
            r, c = divmod(k, n_cols); ax = axes[r][c]
            sp = sub[sub["source"] == path]
            for g in groups:
                gg = sp[sp["group_level"] == g].set_index("pair_s").reindex(pair_order)
                ax.plot(range(len(pair_order)), gg["NES"].values, marker="o", linewidth=1.4,
                        color=GROUP_COLORS.get(g, "black"), label=g)
                sm = (gg["FDR"] < SIG_FDR).fillna(False).to_numpy()
                if sm.any():
                    xi = np.where(sm)[0]
                    ax.scatter(xi, gg["NES"].values[xi], s=60, facecolors="none",
                               edgecolors=GROUP_COLORS.get(g, "black"), linewidth=1.2)
            ax.axhline(0, color="gray", linestyle="--", linewidth=0.5)
            ax.set_xticks(range(len(pair_order)))
            ax.set_xticklabels(pair_order, rotation=30, ha="right", fontsize=7)
            ax.set_ylabel("NES", fontsize=8); ax.set_title(truncate(path, 38), fontsize=8)
            if k == 0:
                ax.legend(fontsize=7, frameon=False)
        for k in range(len(chosen), n_rows * n_cols):
            r, c = divmod(k, n_cols); axes[r][c].set_axis_off()
        fig.suptitle(f"Pathway × age  ·  sex={sx} | level={lv} | ct={ct}  "
                     f"(open circle = FDR<{SIG_FDR})", fontsize=10)
        safe_savefig(fig, plot_root / "trajectory" / _safe(sx) / _safe(lv) / _safe(ct)
                     / "pathway_x_age_top.png")


def plot_leading_edge(le_df, gsea_df, plot_root, top_n_paths, top_n_genes, scope_sex):
    head = gsea_df[(gsea_df["sex"] == scope_sex)
                   & (gsea_df["contrast"].isin(PRIMARY_CONTRASTS))
                   & (gsea_df["level"] == "whole") & (gsea_df["FDR"] < SIG_FDR)]
    if head.empty:
        return
    chosen = head["source"].value_counts().head(top_n_paths).index.tolist()
    le = le_df[(le_df["sex"] == scope_sex)
               & (le_df["contrast"].isin(PRIMARY_CONTRASTS))
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
        contrasts = [c for c in [CONTRAST_EVR, CONTRAST_LVR] if c in sub_g["contrast"].unique()]
        if not contrasts:
            continue
        fig, axes = plt.subplots(1, len(contrasts),
                                  figsize=(6 + 3 * (len(contrasts) - 1), max(5, 0.3 * len(genes) + 2)),
                                  constrained_layout=True, sharey=True)
        if len(contrasts) == 1:
            axes = [axes]
        vmax = safe_vmax(sub_g["log2FC"].values)
        for ax, c in zip(axes, contrasts):
            piv = (sub_g[sub_g["contrast"] == c]
                   .pivot_table(index="gene", columns="celltype", values="log2FC", aggfunc="mean")
                   .reindex(index=genes))
            im = ax.imshow(piv.values, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), aspect="auto")
            ax.set_xticks(range(piv.shape[1]))
            ax.set_xticklabels(piv.columns, rotation=60, ha="right", fontsize=7)
            ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=7)
            ax.set_title(CONTRAST_LABEL.get(c, c), fontsize=9)
            fig.colorbar(im, ax=ax, shrink=0.6, label=r"log$_2$FC")
        fig.suptitle(f"{truncate(path, 70)}  ·  leading-edge × celltype  ·  sex={scope_sex}",
                     fontsize=10)
        safe_savefig(fig, plot_root / "leading_edge" / _safe(scope_sex)
                     / f"{_safe(path)}_genes_x_celltype.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--subcluster", default=None)
    ap.add_argument("--sex-strata", default="combined,M,F")
    ap.add_argument("--top-n", type=int, default=30)
    ap.add_argument("--max-volcano-labels", type=int, default=60)
    ap.add_argument("--skip-volcano", action="store_true")
    ap.add_argument("--skip-dotplot", action="store_true")
    ap.add_argument("--skip-cross-celltype", action="store_true")
    ap.add_argument("--skip-cross-contrast", action="store_true")
    ap.add_argument("--skip-per-cell", action="store_true")
    ap.add_argument("--skip-trajectory", action="store_true")
    ap.add_argument("--skip-leading-edge", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config)
    tissue = cfg.get("tissue")
    rdir = Path(cfg["results_dir"])
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
    gsea = pd.read_csv(gsea_path, low_memory=False)
    sex_strata = [s.strip() for s in args.sex_strata.split(",")]
    gsea = gsea[gsea["sex"].isin(sex_strata)].copy()
    print(f"  {len(gsea):,} GSEA rows (sex {sex_strata})")

    tf = pd.DataFrame()
    if tf_path.is_file():
        tf = pd.read_csv(tf_path, low_memory=False)
        tf = tf[tf["sex"].isin(sex_strata)].copy()
        print(f"  {len(tf):,} TF rows")

    if not args.skip_volcano:
        n = 0
        for keys, block in gsea.groupby(SLICE_ID, observed=True, dropna=False):
            sx, gl, lv, ct = keys
            if not (block["FDR"] < SIG_FDR).any():
                continue
            d = slice_dir(plot_root, sx, gl, lv, ct)
            for contrast in SOLO_CONTRASTS:
                if (block["contrast"] == contrast).any():
                    plot_solo_volcano_pathways(block, contrast,
                                               d / f"volcano_{_safe(contrast)}.png",
                                               args.max_volcano_labels)
            plot_dual_volcano_pathways(block, d / "volcano_DUAL_early_late.png",
                                       args.max_volcano_labels)
            if not tf.empty:
                tfb = tf[(tf["sex"] == sx) & (tf["group_level"] == gl)
                         & (tf["level"] == lv) & (tf["celltype"] == ct)]
                if not tfb.empty:
                    for contrast in SOLO_CONTRASTS:
                        if (tfb["contrast"] == contrast).any():
                            plot_solo_volcano_tf(tfb, contrast,
                                                 d / f"tf_volcano_{_safe(contrast)}.png",
                                                 args.max_volcano_labels)
                    plot_dual_volcano_tf(tfb, d / "tf_volcano_DUAL_early_late.png",
                                         args.max_volcano_labels)
            n += 1
            if n % 50 == 0:
                print(f"    volcano slices: {n} ...")
        print(f"  volcano slices plotted: {n}")

    if not args.skip_dotplot:
        n = 0
        for keys, block in gsea.groupby(SLICE_ID, observed=True, dropna=False):
            sx, gl, lv, ct = keys
            d = slice_dir(plot_root, sx, gl, lv, ct)
            for contrast in SOLO_CONTRASTS:
                sub = block[block["contrast"] == contrast]
                if not sub.empty and (sub["FDR"] < SIG_FDR).any():
                    plot_dotplot_panels(block, contrast,
                                        d / f"dotplot_{_safe(contrast)}.png", args.top_n)
                    n += 1
        print(f"  dotplot panels: {n}")

    if not args.skip_cross_celltype:
        n = 0
        for keys, block in gsea.groupby(["sex", "group_level", "level", "contrast"],
                                         observed=True, dropna=False):
            sx, gl, lv, contrast = keys
            if not (block["FDR"] < SIG_FDR).any():
                continue
            d = slice_dir(plot_root, sx, gl, lv) / "_celltype_heatmaps" / _safe(contrast)
            plot_celltype_pathway_heatmap(block, d / "celltype_x_pathway_heatmap.png", args.top_n)
            if not tf.empty:
                tfb = tf[(tf["sex"] == sx) & (tf["group_level"] == gl)
                         & (tf["level"] == lv) & (tf["contrast"] == contrast)]
                if not tfb.empty:
                    plot_tf_celltype_heatmap(tfb, d / "tf_x_celltype_heatmap.png", args.top_n)
            n += 1
        print(f"  cross-celltype heatmaps: {n}")

    if not args.skip_cross_contrast:
        primary = gsea[gsea["contrast"].isin(PRIMARY_CONTRASTS)]
        nb = nc = nct = 0
        for keys, block in primary.groupby(["sex", "group_level", "level"],
                                            observed=True, dropna=False):
            sx, gl, lv = keys
            cc = plot_root / "cross_contrast" / _safe(sx) / _safe(gl) / _safe(lv)
            plot_bubble(block, cc / "bubble_pathways_x_celltype.png", args.top_n)
            nb += 1
            for ct, sub in block.groupby("celltype", observed=True):
                e = sub[sub["contrast"] == CONTRAST_EVR]
                l = sub[sub["contrast"] == CONTRAST_LVR]
                if e.empty or l.empty:
                    continue
                plot_concordance(e, l, "source", "NES", "FDR",
                                 cc / f"concordance_pathway_{_safe(ct)}.png", args.top_n,
                                 f"Pathway NES concordance · sex={sx}|age={gl}|level={lv}|ct={ct}")
                nc += 1
                if not tf.empty:
                    te = tf[(tf["sex"] == sx) & (tf["group_level"] == gl) & (tf["level"] == lv)
                            & (tf["celltype"] == ct) & (tf["contrast"] == CONTRAST_EVR)]
                    tl = tf[(tf["sex"] == sx) & (tf["group_level"] == gl) & (tf["level"] == lv)
                            & (tf["celltype"] == ct) & (tf["contrast"] == CONTRAST_LVR)]
                    if not te.empty and not tl.empty:
                        fc = _tf_fdr_col(te)
                        plot_concordance(te, tl, "TF", "activity_score", fc,
                                         cc / f"concordance_tf_{_safe(ct)}.png", args.top_n,
                                         f"TF activity concordance · sex={sx}|age={gl}|level={lv}|ct={ct}")
                        nct += 1
        print(f"  cross-contrast: {nb} bubble, {nc} pathway-conc, {nct} TF-conc")

    if not args.skip_per_cell:
        if pc_path.is_file():
            print(f"  loading per-cell {pc_path}")
            pc = ad.read_h5ad(pc_path)
            plot_per_cell_ridges(pc, gsea, plot_root, top_n=20)
            plot_per_cell_umaps(pc, gsea, plot_root, top_n=15)
            print(f"  per-cell ridges + UMAPs done")
        else:
            print(f"  [skip per-cell] {pc_path} not found")

    if not args.skip_trajectory:
        plot_pathway_age_trajectory(gsea, plot_root, top_n=20)
        print(f"  pathway×age trajectories done (within_group_across_age only)")

    if not args.skip_leading_edge:
        if le_path.is_file():
            print(f"  loading leading-edge {le_path}  ({le_path.stat().st_size/1e6:.0f} MB)")
            le = pd.read_csv(le_path, low_memory=False)
            le = le[le["sex"].isin(sex_strata)].copy()
            plot_leading_edge(le, gsea, plot_root, top_n_paths=50, top_n_genes=30,
                              scope_sex=sex_strata[0])
            print(f"  leading-edge heatmaps done")
        else:
            print(f"  [skip leading-edge] {le_path} not found")

    print(f"\n✓ Phase 8c summary complete. Plots under {plot_root}")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    main()
