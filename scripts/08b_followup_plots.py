#!/usr/bin/env python
"""
08b_followup_plots.py — post-hoc plots from the 8b master DE CSV.

Two visualizations, both reading the master CSV directly (no recomputation):

  (1) DISRUPTION — developmental-disruption story from within_group_across_age.
      One figure per (sex × level). Two panels side-by-side:
        Panel A. Mirror bar chart: per cell type, # genes whose age trajectory
                 is LOST under stress (left, red) vs INDUCED by stress (right,
                 blue). The asymmetry is the headline.
        Panel B. Paired |log2FC| boxplots for the LOST-trajectory genes —
                 three boxes per cell type (Relaxed / Early / Late). Shows
                 effect-size collapse from ~1.1 in Relaxed to ~0.6 in stress.

  (2) CONSISTENCY — stress-consistency story from per-age contrasts.
      One figure per (sex × level × age). Stacked horizontal bar per cell
      type:  Early-only sig | Sig in both (convergent) | Late-only sig.
      The gray middle segment IS the convergent stress signature.

Both default to brain-only (within_group_across_age is brain-only by design,
and per-age stress contrasts are only complete in brain — placenta has only
one stress contrast per age so the consistency plot reduces to "all early"
or "all late" and is uninformative).

Sig thresholds (LOCKED, match 08b_de.py): padj<0.05 AND |log2FC|>1.

Output:
  results/{tissue}/plots/08b_de{_subcluster_X}/summary/disruption/{sex}/{level}.png
  results/{tissue}/plots/08b_de{_subcluster_X}/summary/consistency/{sex}/{level}/{age}.png

Usage:
  uv run python scripts/08b_followup_plots.py --config config/brain.yaml
  uv run python scripts/08b_followup_plots.py --config config/brain.yaml --subcluster immune
  uv run python scripts/08b_followup_plots.py --config config/brain.yaml --plots disruption
  uv run python scripts/08b_followup_plots.py --config config/brain.yaml --plots consistency
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _utils import load_config, phase_table_dir


# ---------------------------------------------------------------------------
# Constants (must match 08b_de.py / 08b_de_summary.py)
# ---------------------------------------------------------------------------

PADJ_THR = 0.05
LFC_THR = 1.0

GROUPS = ["Relaxed", "Early_Stress", "Late_Stress"]
SEX_ORDER = ["combined", "M", "F"]
AGE_ORDER = {
    "brain":    ["P1", "4W", "3mo"],
    "placenta": ["E12.5", "E18.5"],
}

# Slimmest columns needed for both plots
NEEDED_COLS = ["contrast", "test_method", "sex", "group_level", "pair",
               "level", "celltype", "gene", "log2FC", "padj"]

# Color palette (consistent across panels)
COL_LOST   = "#c0392b"  # trajectory lost (red)
COL_GAINED = "#2980b9"  # trajectory induced (blue)
COL_BOTH   = "#34495e"  # convergent (dark gray)
COL_RELAX  = "#7f8c8d"  # Relaxed baseline (mid gray)
COL_EARLY  = "#c0392b"
COL_LATE   = "#2980b9"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def slugify(s):
    return re.sub(r"[^0-9a-zA-Z]+", "_", str(s)).strip("_").lower()


def safe_fig(fig, out: Path, dpi=140):
    out.parent.mkdir(parents=True, exist_ok=True)
    # constrained_layout figures don't play with bbox_inches='tight'. Detect
    # and skip if so. (See matplotlib docs — using both is contradictory.)
    using_constrained = getattr(fig, "get_constrained_layout",
                                lambda: False)()
    if using_constrained:
        fig.savefig(out, dpi=dpi)
    else:
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def is_sig(df):
    """Boolean Series — padj<0.05 & |log2FC|>1, NaN-safe."""
    return (df["padj"].notna() & (df["padj"] < PADJ_THR)
            & df["log2FC"].notna() & (df["log2FC"].abs() > LFC_THR))


def collapse_age_pairs(w):
    """within_group_across_age has multiple pairwise rows per gene per
    (slice × group). Collapse to one row per (slice × group × gene) by
    keeping the most-significant pair (smallest padj). Without this,
    `set_index('gene')` would produce a non-unique index."""
    return (w.sort_values("padj")
             .drop_duplicates(
                 ["sex", "level", "celltype", "group_level", "gene"],
                 keep="first"))


def classify_disruption(rel, early, late):
    """Three gene-indexed frames -> dict direction-class -> set(gene)."""
    R = set(rel.index[is_sig(rel)])
    E = set(early.index[is_sig(early)])
    L = set(late.index[is_sig(late)])
    return {
        "universal":     R & E & L,
        "relaxed_only":  R - E - L,
        "stress_shared": (E & L) - R,
        "early_only":    E - R - L,
        "late_only":     L - R - E,
    }


# ---------------------------------------------------------------------------
# Plot 1 — Developmental disruption (mirror bar + effect-size boxplots)
# ---------------------------------------------------------------------------

def plot_disruption(w, out_dir, tissue, sex_label, level, min_n_lost=5):
    """Two-panel paper-quality figure for one (sex × level) slice.

    Layout fixes vs initial version:
      - constrained_layout=True + dedicated suptitle space
      - filter cell types with n_lost < min_n_lost (drops degenerate boxes)
      - in-bar count labels rendered at fontsize 10 with white text on red/blue
      - Panel B has NO redundant y-tick labels (shares cell-type order with A);
        each cell type's box-group has a tiny "n=..." annotation at the left
      - legends moved OUT of the plotting area via bbox_to_anchor
      - pool-confound caveat demoted to a small figure-bottom footnote
    """
    sub = w[(w["sex"] == sex_label) & (w["level"] == level)]
    if sub.empty:
        return False

    celltypes, n_lost, n_gained = [], [], []
    lfc_by_ct_grp = {}

    for ct, g in sub.groupby("celltype", observed=True):
        per_grp = {grp: gg.set_index("gene")[["log2FC", "padj"]]
                   for grp, gg in g.groupby("group_level", observed=True)}
        if not all(k in per_grp for k in GROUPS):
            continue
        rel, early, late = (per_grp["Relaxed"], per_grp["Early_Stress"],
                            per_grp["Late_Stress"])
        classes = classify_disruption(rel, early, late)
        nl, ng = len(classes["relaxed_only"]), len(classes["stress_shared"])
        if nl < min_n_lost:
            continue  # noise — skip
        celltypes.append(ct)
        n_lost.append(nl)
        n_gained.append(ng)
        ro_genes = list(classes["relaxed_only"])
        for grp_lbl, frame in [("Relaxed", rel), ("Early_Stress", early),
                                ("Late_Stress", late)]:
            present = frame.index.intersection(ro_genes)
            arr = (frame.loc[present, "log2FC"].abs().values
                   if len(present) else np.array([]))
            lfc_by_ct_grp[(ct, grp_lbl)] = arr

    if not celltypes:
        return False

    # Sort by descending n_lost
    order = np.argsort([-n for n in n_lost])
    celltypes = [celltypes[i] for i in order]
    n_lost = [n_lost[i] for i in order]
    n_gained = [n_gained[i] for i in order]

    nrow = len(celltypes)
    fig_w = 15
    fig_h = max(5.0, 0.6 * nrow + 2.4)
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(fig_w, fig_h),
        gridspec_kw=dict(width_ratios=[1, 1.6], wspace=0.25),
        constrained_layout=True,
    )

    # ---------- Panel A: mirror bar ----------
    y_pos = np.arange(nrow)
    bar_h = 0.7
    axA.barh(y_pos, [-n for n in n_lost], bar_h, color=COL_LOST,
             edgecolor="black", lw=0.4)
    axA.barh(y_pos, n_gained, bar_h, color=COL_GAINED,
             edgecolor="black", lw=0.4)
    max_n = max(max(n_lost), max(n_gained) if n_gained else 0, 1)
    label_offset = max_n * 0.015        # tiny inset from bar tip
    for i, (nl, ng) in enumerate(zip(n_lost, n_gained)):
        # LOST: bar extends to -nl. Label sits inside the bar near the tip.
        if nl > 0:
            inside = nl > max_n * 0.10
            axA.text(-nl + label_offset if inside else -nl - label_offset,
                     i, f"{nl}",
                     ha="left" if inside else "right",
                     va="center", fontsize=10,
                     color="white" if inside else "black",
                     fontweight="bold")
        if ng > 0:
            inside = ng > max_n * 0.10
            axA.text(ng - label_offset if inside else ng + label_offset,
                     i, f"{ng}",
                     ha="right" if inside else "left",
                     va="center", fontsize=10,
                     color="white" if inside else "black",
                     fontweight="bold")
    axA.axvline(0, color="k", lw=0.7)
    axA.set_yticks(y_pos)
    axA.set_yticklabels(celltypes, fontsize=10)
    axA.invert_yaxis()
    # Symmetric x-labels (absolute values)
    xt = axA.get_xticks()
    axA.set_xticklabels([f"{abs(int(x))}" for x in xt])
    axA.set_xlabel("# age-DE genes", fontsize=10)
    axA.set_title(
        "A. Asymmetry: developmental trajectories\n"
        "← LOST under stress     GAINED under stress →",
        fontsize=11)
    # Legend OUTSIDE the plot (top-right of the panel)
    axA.legend(handles=[
        mpatches.Patch(color=COL_LOST, label="LOST (sig in Relaxed only)"),
        mpatches.Patch(color=COL_GAINED, label="GAINED (sig in Early ∩ Late)"),
    ], fontsize=8, loc="upper left", bbox_to_anchor=(0.0, -0.10),
       frameon=False, ncol=2)
    axA.spines[["top", "right"]].set_visible(False)

    # ---------- Panel B: paired |LFC| boxplots ----------
    box_positions, box_data, box_colors = [], [], []
    grp_colors = {"Relaxed": COL_RELAX,
                  "Early_Stress": COL_LOST,
                  "Late_Stress": COL_LATE}
    for i, ct in enumerate(celltypes):
        for j, grp in enumerate(GROUPS):
            arr = lfc_by_ct_grp.get((ct, grp), np.array([]))
            box_positions.append(i * 4 + j)
            box_data.append(arr if len(arr) > 0 else np.array([np.nan]))
            box_colors.append(grp_colors[grp])
    bp = axB.boxplot(box_data, positions=box_positions, widths=0.75,
                     patch_artist=True, showfliers=False,
                     medianprops=dict(color="black", lw=1.2),
                     boxprops=dict(lw=0.5),
                     whiskerprops=dict(lw=0.5),
                     capprops=dict(lw=0.5),
                     vert=False)
    for patch, col in zip(bp["boxes"], box_colors):
        patch.set_facecolor(col)
        patch.set_alpha(0.78)
    # Light horizontal separators between cell-type groups for visual grouping
    for i in range(1, nrow):
        axB.axhline(i * 4 - 0.5, color="0.85", lw=0.4, zorder=0)
    # NO duplicate y-tick labels — Panel A already names the cell types.
    # Hide tick labels but keep the ticks for alignment.
    axB.set_yticks([i * 4 + 1 for i in range(nrow)])
    axB.set_yticklabels([])
    axB.tick_params(axis="y", length=0)
    # Add "n=NNN" annotation at the left edge of each cell-type group
    xmin = -0.05  # we'll set xlim after
    axB.invert_yaxis()
    axB.axvline(LFC_THR, color="k", lw=0.5, ls="--", alpha=0.5)
    # Right-side small "n=..." labels for each cell type's box group
    for i, ct in enumerate(celltypes):
        axB.text(0.01, i * 4 + 1, f"n={n_lost[i]}",
                 transform=axB.get_yaxis_transform(),
                 ha="left", va="center", fontsize=8, color="0.35",
                 fontstyle="italic")
    axB.set_xlim(left=0)
    axB.set_xlabel("|log2 fold change|  (LOST-trajectory genes only)",
                   fontsize=10)
    axB.set_title(
        "B. Effect-size collapse on lost-trajectory genes\n"
        "(per cell type, three boxes top→bottom = Relaxed / Early / Late)",
        fontsize=11)
    axB.legend(handles=[
        mpatches.Patch(color=COL_RELAX, label="Relaxed (baseline)"),
        mpatches.Patch(color=COL_LOST, label="Early stress"),
        mpatches.Patch(color=COL_LATE, label="Late stress"),
    ], fontsize=8, loc="upper left", bbox_to_anchor=(0.0, -0.10),
       frameon=False, ncol=3)
    axB.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        f"{tissue} | sex={sex_label} | level={level} | "
        f"within_group_across_age "
        f"(padj<{PADJ_THR} & |log2FC|>{LFC_THR})",
        fontsize=12)
    # Pool-confound caveat as a small italic footnote INSIDE the figure
    fig.text(0.5, 0.005,
             "Pool-confounded contrast — interpret with care. "
             f"Cell types with n_lost < {min_n_lost} are hidden.",
             ha="center", fontsize=7, style="italic", color="0.4")

    out = out_dir / sex_label / f"{slugify(level)}.png"
    safe_fig(fig, out)
    return True


# ---------------------------------------------------------------------------
# Plot 2 — Stress-consistency (Early-only / Both / Late-only stacked bar)
# ---------------------------------------------------------------------------

def plot_consistency(df_per_age, out_dir, tissue):
    """Per (sex × level × age), one stacked-bar figure across cell types."""
    drawn = 0
    for (sex_label, level, age), g in df_per_age.groupby(
            ["sex", "level", "group_level"], observed=True):
        early = g[g["contrast"].str.startswith("early_vs_relaxed")]
        late  = g[g["contrast"].str.startswith("late_vs_relaxed")]
        if early.empty or late.empty:
            # Need BOTH contrasts at this age to compute "convergent" overlap
            continue
        e_sig = early[is_sig(early)]
        l_sig = late[is_sig(late)]

        rows = []
        for ct in sorted(set(g["celltype"].astype(str).unique())):
            E = set(e_sig.loc[e_sig.celltype == ct, "gene"].dropna().astype(str))
            L = set(l_sig.loc[l_sig.celltype == ct, "gene"].dropna().astype(str))
            both = E & L
            early_only = E - L
            late_only  = L - E
            if (len(both) + len(early_only) + len(late_only)) == 0:
                continue
            rows.append(dict(celltype=ct,
                            early_only=len(early_only),
                            both=len(both),
                            late_only=len(late_only)))
        if not rows:
            continue
        d = pd.DataFrame(rows)
        d["total"] = d["early_only"] + d["both"] + d["late_only"]
        d = d.sort_values("total", ascending=True)  # smallest at top, biggest at bottom

        n = len(d)
        fig_h = max(3.5, 0.5 * n + 1.8)
        fig, ax = plt.subplots(figsize=(10, fig_h))
        y = np.arange(n)
        bar_h = 0.7

        ax.barh(y, d["early_only"], bar_h, color=COL_EARLY,
                edgecolor="black", lw=0.4)
        ax.barh(y, d["both"], bar_h, left=d["early_only"], color=COL_BOTH,
                edgecolor="black", lw=0.4)
        ax.barh(y, d["late_only"], bar_h,
                left=d["early_only"] + d["both"], color=COL_LATE,
                edgecolor="black", lw=0.4)

        # In-bar segment labels (only if segment is big enough to be readable)
        for i, row in enumerate(d.itertuples()):
            x_offset = 0
            for seg, val, color in [
                ("E", row.early_only, "white"),
                ("∩", row.both,       "white"),
                ("L", row.late_only,  "white"),
            ]:
                if val > 0 and val > row.total * 0.06:
                    ax.text(x_offset + val / 2, i, f"{seg}:{val}",
                            ha="center", va="center", fontsize=7,
                            color=color, fontweight="bold")
                x_offset += val
            # Total label at end of bar
            ax.text(row.total, i, f"  {row.total}",
                    ha="left", va="center", fontsize=8, color="black")

        ax.set_yticks(y)
        ax.set_yticklabels(d["celltype"], fontsize=9)
        ax.set_xlabel("# significant DEGs (vs. Relaxed)", fontsize=9)
        ax.set_title(
            f"{tissue} | sex={sex_label} | level={level} | age={age}\n"
            f"Stress-consistency: Early-only | Sig in BOTH (convergent) | Late-only\n"
            f"(padj<{PADJ_THR} & |log2FC|>{LFC_THR}; central dark segment = "
            f"the convergent stress signature)",
            fontsize=10)
        ax.legend(handles=[
            mpatches.Patch(color=COL_EARLY, label="Early-only sig"),
            mpatches.Patch(color=COL_BOTH,  label="Sig in BOTH (convergent)"),
            mpatches.Patch(color=COL_LATE,  label="Late-only sig"),
        ], fontsize=8, loc="lower right", frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
        out = out_dir / sex_label / slugify(level) / f"{slugify(age)}.png"
        safe_fig(fig, out)
        drawn += 1
    return drawn


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--subcluster", default=None,
                    help="Read the subcluster master CSV instead of the main one.")
    ap.add_argument("--plots", default="all",
                    choices=["all", "disruption", "consistency"],
                    help="Which plot type(s) to draw (default: both).")
    args = ap.parse_args()

    print("\n=== 08b followup plots ===")
    cfg = load_config(args.config)
    tissue = cfg.get("tissue")
    print(f"  Tissue: {tissue}")
    if args.subcluster:
        print(f"  Subcluster mode: {args.subcluster}")

    table_dir = phase_table_dir(cfg, "08b_de")
    suffix = f"_subcluster_{args.subcluster}" if args.subcluster else ""
    csv_path = table_dir / f"08b_de_results{suffix}.csv"
    if not csv_path.is_file():
        sys.exit(f"ERROR: master CSV not found: {csv_path}")
    print(f"Reading {csv_path} ({csv_path.stat().st_size / 1e6:.1f} MB)...")
    df = pd.read_csv(csv_path, usecols=lambda c: c in NEEDED_COLS,
                     low_memory=False)
    print(f"  {len(df):,} rows loaded.")

    plot_root = (Path(cfg["results_dir"]) / "plots"
                 / ("08b_de" + suffix) / "summary")

    # ---------- Disruption ----------
    if args.plots in ("all", "disruption"):
        w = df[(df["test_method"] == "Wald")
               & (df["contrast"] == "within_group_across_age")]
        if w.empty:
            print("\n  [skip] DISRUPTION: no within_group_across_age Wald rows "
                  "(expected for placenta).")
        else:
            print(f"\n  DISRUPTION: {len(w):,} within_group_across_age rows.")
            w = collapse_age_pairs(w)
            print(f"             after collapsing age-pairs: {len(w):,} rows.")
            out_dir = plot_root / "disruption"
            print(f"  -> {out_dir}")
            n_drawn = 0
            for sex_label in SEX_ORDER:
                if sex_label not in w["sex"].unique():
                    continue
                for level in sorted(w["level"].astype(str).unique()):
                    if plot_disruption(w, out_dir, tissue, sex_label, level):
                        n_drawn += 1
            print(f"  disruption: drew {n_drawn} figure(s).")

    # ---------- Consistency ----------
    if args.plots in ("all", "consistency"):
        per_age_contrasts = ("early_vs_relaxed_per_age",
                             "late_vs_relaxed_per_age",
                             "early_vs_relaxed_E12.5",
                             "late_vs_relaxed_E18.5")
        d = df[(df["test_method"] == "Wald")
               & df["contrast"].isin(per_age_contrasts)]
        if d.empty:
            print("\n  [skip] CONSISTENCY: no per-age Wald rows.")
        else:
            print(f"\n  CONSISTENCY: {len(d):,} per-age Wald rows.")
            out_dir = plot_root / "consistency"
            print(f"  -> {out_dir}")
            n_drawn = plot_consistency(d, out_dir, tissue)
            print(f"  consistency: drew {n_drawn} figure(s) "
                  f"(skipped (sex,level,age) groups missing one of the two "
                  f"stress contrasts).")

    print("\n  ✓ Done.\n")


if __name__ == "__main__":
    main()
