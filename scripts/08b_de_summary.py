#!/usr/bin/env python
"""
08b_de_summary.py — cross-slice summary plots for Phase 8b DE results.

Reads the master CSV(s) written by 08b_de.py and produces the cross-slice
visualisations agreed in the project doc:

  1. Top-DEG multi-panel heatmap (rows = top union genes, cols = cell type × age)
  2. UpSet plot (overlap of significant gene sets across slices)
  3. n_DEGs facet bar chart (height = #sig DEGs, faceted by contrast, stacked
     up vs down)
  4. Bubble / dot plot (variant of heatmap; size = -log10(padj), color = log2FC)
  5. RRHO-lite cross-contrast scatter (Early vs Late log2FC concordance per slice)
  6. Top-gene dotplot (Scanpy-style: gene × group × age, color = mean lognorm)
  7. Per-celltype volcano grid (small multiples: panel per age × contrast)
  8. Venn diagrams (directional: all-sig + up-only + down-only) for
        - Per-age Early ∩ Late at each age (2-way Venns)
        - Across-age P1 ∩ 4W ∩ 3mo within Early-vs-Relaxed and Late-vs-Relaxed
          (3-way Venns)

Independent from 08b_de.py — re-run this to tweak figures without re-fitting DE.

Significance thresholds (LOCKED, match 08b_de.py):
    padj < 0.05 AND |log2FC| > 1   (i.e., >= 2x fold change)

Main + subcluster modes mirror 08b_de.py:
    main:        reads results/{tissue}/tables/08b_de/08b_de_results.csv
                 writes results/{tissue}/plots/08b_de/summary/...
    subcluster:  reads ...08b_de_results_subcluster_{slug}.csv
                 writes results/{tissue}/plots/08b_de_subcluster_{slug}/summary/...

Master CSV is large (multi-GB for brain main). Read only the columns each plot
needs via `usecols=` and reuse across plot functions; never load full rows.

Usage:
    uv run python scripts/08b_de_summary.py --config config/brain.yaml --n-jobs 16
    uv run python scripts/08b_de_summary.py --config config/placenta.yaml --n-jobs 16
    uv run python scripts/08b_de_summary.py --config config/brain.yaml --subcluster immune
    uv run python scripts/08b_de_summary.py --config config/brain.yaml --plots heatmap,venn,bar
"""

import argparse
import re
import sys
import warnings
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _utils import load_config, phase_table_dir, parallel_map


# ---------------------------------------------------------------------------
# Constants (must match 08b_de.py)
# ---------------------------------------------------------------------------

PADJ_THR = 0.05
LFC_THR = 1.0

AGE_ORDER = {
    "brain":    ["P1", "4W", "3mo"],
    "placenta": ["E12.5", "E18.5"],
}
GROUP_ORDER = ["Relaxed", "Early_Stress", "Late_Stress"]
SEX_ORDER = ["combined", "M", "F"]

# Per-age pairwise contrasts the script knows how to summarise. The two
# "primary stress" contrasts feed Venns / heatmaps / RRHO-lite; the
# omnibus/interaction/cross-age contrasts feed the bar chart (n-DEGs view) and
# the master heatmap.
STRESS_CONTRASTS = ("early_vs_relaxed_per_age", "late_vs_relaxed_per_age",
                    "early_vs_relaxed_E12.5", "late_vs_relaxed_E18.5")

# -----------------------------------------------------------------------------
# Visualization blocklist — genes that survive Phase 7 gating but are known to
# dominate low-signal slices for technical / developmental reasons. They stay
# in the master CSV (data preserved), still appear as gray dots in volcanoes /
# scatters / Venn sets, and still count toward "n sig DEGs" totals. They are
# only excluded from CURATED TOP-N selections used for labels and heatmap rows
# so that real biology surfaces. Override at CLI with --no-blocklist.
#
# Categories:
#   - Hemoglobin / heme synthesis: P1 brain carries developmental erythroid
#     ambient burden (nucleated erythroblasts in vasculature); SoupX strips
#     ambient but not residual cellular Hb. Diagnosed 2026-06: 679 sig Hb
#     rows, ~650 in within_group_across_age (P1 vs older within group),
#     direction consistent across all groups including Relaxed -> ambient.
#   - Sex-linked: pop in combined-sex analyses when sex stratification is
#     imperfect. Already on the HVG exclusion list but they survive into DE.
#   - Mitochondrial: snRNA-seq should have ~0% mito but ambient contributes
#     a baseline.
# -----------------------------------------------------------------------------
BLOCKLIST_FOR_VIZ = {
    # Hemoglobin / heme
    "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Hbb-b1", "Hbb-b2",
    "Hbb-y", "Hbb-bh1", "Hbb-bh2", "Alas2",
    # Sex-linked
    "Xist", "Tsix", "Ddx3y", "Uty", "Eif2s3y", "Kdm5d", "Eif2s3x",
}
# mt-* genes match by prefix (catch all mitochondrial transcripts at once)
BLOCKLIST_PREFIXES = ("mt-",)


def is_blocklisted(gene: str) -> bool:
    if not isinstance(gene, str):
        return False
    if gene in BLOCKLIST_FOR_VIZ:
        return True
    return any(gene.startswith(p) for p in BLOCKLIST_PREFIXES)


def filter_blocklist(genes_or_df, gene_col: str = "gene",
                     enabled: bool = True):
    """Drop blocklisted genes. Accepts either an iterable of gene symbols
    (returns filtered list) or a DataFrame (returns filtered DataFrame)."""
    if not enabled:
        return genes_or_df
    if isinstance(genes_or_df, pd.DataFrame):
        return genes_or_df[~genes_or_df[gene_col].apply(is_blocklisted)]
    return [g for g in genes_or_df if not is_blocklisted(g)]

# Subset of columns we need from the master CSV. log2FC + padj are central;
# stat / lfcSE / n_donors_* / flag / note are not used by these plots.
NEEDED_COLS = ["contrast", "test_method", "sex", "group_level", "pair", "level",
               "celltype", "gene", "log2FC", "padj", "direction", "reliability"]

ALL_PLOTS = ["heatmap", "upset", "bar", "bubble", "rrho", "dotplot", "grid", "venn"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def slugify(name: str) -> str:
    return re.sub(r"[^0-9a-zA-Z]+", "_", str(name)).strip("_").lower()


def ordered(values, order):
    vals = list(dict.fromkeys(values))
    return [v for v in order if v in vals] + sorted(v for v in vals if v not in order)


def is_sig_mask(df, lfc_thr=LFC_THR, padj_thr=PADJ_THR):
    """Project-wide significance: Wald rows with padj < padj_thr AND
    |log2FC| > lfc_thr. LRT rows (NaN log2FC by design) are excluded — they
    can't be put on a volcano axis."""
    return ((df["test_method"] == "Wald")
            & df["padj"].notna()
            & (df["padj"] < padj_thr)
            & df["log2FC"].notna()
            & (df["log2FC"].abs() > lfc_thr))


def safe_fig(fig, out: Path, dpi=140):
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def read_master(cfg: dict, subcluster: str | None) -> tuple[pd.DataFrame, Path, Path]:
    """Read the master CSV (only NEEDED_COLS) and return (df, csv_path,
    plot_root). plot_root mirrors 08b_de.py's plot_root + a `summary/`
    subdirectory."""
    table_dir = phase_table_dir(cfg, "08b_de")
    suffix = f"_subcluster_{subcluster}" if subcluster else ""
    csv_path = table_dir / f"08b_de_results{suffix}.csv"
    if not csv_path.is_file():
        sys.exit(f"ERROR: master CSV not found: {csv_path}\n"
                 f"  Run 08b_de.py {'--subcluster ' + subcluster if subcluster else ''} first.")
    print(f"  Reading {csv_path} (size {csv_path.stat().st_size / 1e6:.1f} MB)...")
    # Use only the columns we need (saves several minutes on multi-GB CSVs)
    df = pd.read_csv(csv_path, usecols=lambda c: c in NEEDED_COLS,
                     dtype={"reliability": "string", "direction": "string"})
    print(f"    {len(df):,} rows loaded.")
    plot_root = (Path(cfg["results_dir"]) / "plots"
                 / ("08b_de" + suffix) / "summary")
    plot_root.mkdir(parents=True, exist_ok=True)
    return df, csv_path, plot_root


def pair_to_test_ref(pair_str):
    """The `pair` column is stored as the string repr of a list, e.g.
    "['Early_Stress', 'Relaxed']". Recover (test, ref). Returns (None, None)
    if not parseable."""
    if pair_str is None or pd.isna(pair_str):
        return None, None
    s = str(pair_str).strip()
    # eval-free parse — strip [] '' and split by comma
    inner = s.strip("[]")
    parts = [p.strip().strip("'").strip('"') for p in inner.split(",")]
    if len(parts) == 2:
        return parts[0], parts[1]
    return None, None


# ---------------------------------------------------------------------------
# Plot 1 — Top-DEG multi-panel heatmap
# ---------------------------------------------------------------------------

def plot_top_deg_heatmap(df_sig_subset, df_all_subset, title, out, max_genes=50,
                        vlim=3.0, ct_order=None, use_blocklist=True):
    """Top-N-union-gene heatmap.

    df_sig_subset : already filtered to sig Wald rows for one (sex × contrast).
                    Used to pick the top-N union gene list (by min padj).
    df_all_subset : the same (sex × contrast) slice but unfiltered. Used to
                    fill log2FC values for the picked genes — even in slices
                    where they're not significant — so the heatmap shows
                    direction-of-effect rather than just presence.
    use_blocklist : when True, blocklisted genes (Hb / sex-linked / mito)
                    are excluded from the top-N selection but the unfiltered
                    `df_all_subset` is left alone, so they could still appear
                    if a blocklisted gene leaks in via another route.
    Rows = top max_genes. Cols = (celltype, group_level).
    Cell text = log2FC (rounded). Bold border = padj < PADJ_THR & |LFC| > LFC_THR.
    """
    if df_sig_subset.empty:
        return
    candidates = (df_sig_subset
                  .sort_values("padj")
                  .drop_duplicates("gene"))
    if use_blocklist:
        candidates = filter_blocklist(candidates)
    top_genes = candidates.head(max_genes)["gene"].tolist()
    if not top_genes:
        return

    # Restrict the unfiltered subset to those genes for the heatmap fill
    df_plot = df_all_subset[df_all_subset["gene"].isin(top_genes)].copy()
    df_plot["group_level"] = df_plot["group_level"].astype(str)
    df_plot["celltype"] = df_plot["celltype"].astype(str)
    # Build composite column label
    df_plot["_col"] = df_plot["celltype"] + " | " + df_plot["group_level"]

    lfc = (df_plot.pivot_table(index="gene", columns="_col",
                               values="log2FC", aggfunc="first")
                 .reindex(index=top_genes))
    padj = (df_plot.pivot_table(index="gene", columns="_col",
                                values="padj", aggfunc="first")
                  .reindex(index=top_genes))

    # Order cols: by celltype (optional supplied order), then age
    cols = list(lfc.columns)
    if ct_order:
        cols = sorted(cols, key=lambda c: (
            ct_order.index(c.split(" | ")[0]) if c.split(" | ")[0] in ct_order
            else len(ct_order), c))
    lfc = lfc.reindex(columns=cols)
    padj = padj.reindex(columns=cols)

    nrow, ncol = lfc.shape
    if nrow == 0 or ncol == 0:
        return
    fig_w = max(6.0, 0.45 * ncol + 3.5)
    fig_h = max(4.0, 0.22 * nrow + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    data = lfc.values.astype(float)
    im = ax.imshow(np.ma.masked_invalid(data), cmap="RdBu_r",
                   vmin=-vlim, vmax=vlim, aspect="auto")
    ax.set_xticks(range(ncol))
    ax.set_xticklabels(cols, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(nrow))
    ax.set_yticklabels(lfc.index, fontsize=7)

    for i in range(nrow):
        for j in range(ncol):
            v = data[i, j]
            if np.isnan(v):
                ax.text(j, i, "·", ha="center", va="center", fontsize=6, color="0.7")
                continue
            p = padj.values[i, j]
            sig = bool(pd.notna(p) and p < PADJ_THR and abs(v) > LFC_THR)
            ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                    fontsize=5.5,
                    fontweight="bold" if sig else "normal",
                    color="white" if abs(v) > vlim * 0.6 else "black")
            if sig:
                ax.add_patch(plt.Rectangle((j - 0.46, i - 0.46), 0.92, 0.92,
                                           fill=False, edgecolor="black",
                                           lw=1.6, zorder=5))
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("log2 fold change", fontsize=7)
    ax.set_title(title, fontsize=10)
    fig.text(0.5, -0.02,
             f"Top {len(top_genes)} union genes (by min padj). "
             f"Bold-outlined cells: padj<{PADJ_THR} & |log2FC|>{LFC_THR}.",
             ha="center", fontsize=6, style="italic")
    safe_fig(fig, out)


# ---------------------------------------------------------------------------
# Plot 2 — UpSet plot of significant gene sets
# ---------------------------------------------------------------------------

def plot_upset(df_sig_subset, title, out, max_intersections=15, min_set_size=5):
    """Top intersections of sig gene sets across (cell type × age) slices.
    Falls back gracefully if upsetplot isn't installed.

    NOTE: UpSet builds its panels with a custom gridspec. Saving with
    bbox_inches='tight' (our default in safe_fig) re-tightens the layout in
    a way that conflicts with that gridspec and either errors or produces
    a clipped figure. We bypass safe_fig here and save with default bbox.
    """
    try:
        from upsetplot import UpSet, from_contents
    except ImportError:
        return  # quietly skip; main() already warned
    if df_sig_subset.empty:
        return
    # Build {slice_label -> {gene set}}
    sets = {}
    df_sig_subset = df_sig_subset.assign(
        _slice=df_sig_subset["celltype"].astype(str)
               + " | " + df_sig_subset["group_level"].astype(str))
    for label, g in df_sig_subset.groupby("_slice"):
        genes = set(g["gene"].dropna().astype(str))
        if len(genes) >= min_set_size:
            sets[label] = genes
    if len(sets) < 2:
        return
    # Cap the number of sets to avoid unreadable UpSet (most useful at <12)
    if len(sets) > 10:
        sets = dict(sorted(sets.items(), key=lambda kv: -len(kv[1]))[:10])

    data = from_contents(sets)
    fig = plt.figure(figsize=(max(8, 0.6 * len(sets) + 4), 6))
    # show_counts=True triggers a matplotlib >=3.8 incompatibility in
    # upsetplot 0.9.0 (passes array-shape positions to ax.text -> TypeError
    # "only 0-dimensional arrays can be converted"). Set False and add the
    # counts ourselves below.
    upset = UpSet(data, subset_size="count",
                  show_counts=False, sort_by="cardinality",
                  max_subset_rank=max_intersections)
    axes_dict = upset.plot(fig=fig)

    # Re-add counts on both bar charts. Use ax.bar_label on the actual
    # BarContainers (skips background-shading patches that my earlier
    # ax.patches loop was picking up by mistake).
    for ax_name in ("intersections", "totals"):
        ax = axes_dict.get(ax_name)
        if ax is None:
            continue
        for container in ax.containers:
            try:
                ax.bar_label(container, fontsize=7, padding=2, fmt="%d")
            except Exception as e:
                # bar_label can fail on degenerate containers (e.g. all NaN);
                # log and move on rather than aborting the whole figure.
                print(f"    [upset bar_label] {ax_name}: "
                      f"{type(e).__name__}: {e}")

    fig.suptitle(title, fontsize=10, y=1.02)
    out.parent.mkdir(parents=True, exist_ok=True)
    # Do NOT use bbox_inches='tight' here — it fights UpSet's internal
    # gridspec layout. pad_inches sets a small margin instead.
    fig.savefig(out, dpi=140, bbox_inches=None, pad_inches=0.3)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Plot 3 — n_DEGs facet bar chart (stacked up vs down)
# ---------------------------------------------------------------------------

def plot_n_degs_bar(df_sig_master, out_dir, tissue):
    """One figure per sex stratum, faceted by contrast.

    Bars: one per slice (celltype × group_level × level). Height = #sig DEGs.
    Stacked by direction (up = warm color, down = cool). One panel per contrast.
    """
    if df_sig_master.empty:
        return
    df = df_sig_master.assign(
        _slice=df_sig_master["celltype"].astype(str)
               + " | " + df_sig_master["group_level"].astype(str)
               + " | " + df_sig_master["level"].astype(str))
    counts = (df.groupby(["sex", "contrast", "_slice", "direction"],
                         observed=True).size()
               .unstack("direction", fill_value=0))
    if "up" not in counts.columns:
        counts["up"] = 0
    if "down" not in counts.columns:
        counts["down"] = 0
    counts = counts[["up", "down"]]

    for sex_label in ordered(counts.index.get_level_values("sex").unique(), SEX_ORDER):
        sub = counts.xs(sex_label, level="sex")
        contrasts = sorted(sub.index.get_level_values("contrast").unique())
        if not contrasts:
            continue
        n = len(contrasts)
        fig, axes = plt.subplots(n, 1, figsize=(max(8, 0.16 * 50), 2.6 * n + 1),
                                 squeeze=False)
        for ax, cname in zip(axes[:, 0], contrasts):
            d = sub.xs(cname, level="contrast")
            # Keep only slices with ANY hits, sort by total
            d = d[d.sum(axis=1) > 0].sort_values(by=["up", "down"],
                                                 ascending=False).head(60)
            if d.empty:
                ax.text(0.5, 0.5, "no significant slices", ha="center",
                        va="center", transform=ax.transAxes, fontsize=8)
                ax.set_title(cname, fontsize=9, loc="left")
                ax.set_axis_off()
                continue
            x = np.arange(len(d))
            ax.bar(x, d["up"], color="#c0392b", label="up", edgecolor="white", lw=0.3)
            ax.bar(x, d["down"], bottom=d["up"], color="#2980b9", label="down",
                   edgecolor="white", lw=0.3)
            ax.set_xticks(x)
            ax.set_xticklabels(d.index, rotation=70, ha="right", fontsize=6)
            ax.set_ylabel("# sig DEGs", fontsize=8)
            ax.set_title(f"{cname}  (top {len(d)} slices)", fontsize=9, loc="left")
            ax.spines[["top", "right"]].set_visible(False)
            ax.legend(fontsize=6, frameon=False)
        fig.suptitle(f"{tissue} | sex={sex_label} | sig DEGs per slice "
                     f"(padj<{PADJ_THR} & |log2FC|>{LFC_THR})", fontsize=11)
        safe_fig(fig, out_dir / f"n_DEGs_bar_{sex_label}.png")


# ---------------------------------------------------------------------------
# Plot 4 — Bubble (variant of heatmap; smaller, easier to scan)
# ---------------------------------------------------------------------------

def plot_bubble(df_all_subset, df_sig_subset, title, out, max_genes=30,
                use_blocklist=True):
    if df_sig_subset.empty:
        return
    candidates = df_sig_subset.sort_values("padj").drop_duplicates("gene")
    if use_blocklist:
        candidates = filter_blocklist(candidates)
    top_genes = candidates.head(max_genes)["gene"].tolist()
    if not top_genes:
        return
    df = df_all_subset[df_all_subset["gene"].isin(top_genes)].copy()
    df["_col"] = (df["celltype"].astype(str) + " | "
                  + df["group_level"].astype(str))
    df = df.dropna(subset=["log2FC", "padj"])
    if df.empty:
        return
    cats_y = top_genes
    cats_x = sorted(df["_col"].unique())
    x_idx = {c: i for i, c in enumerate(cats_x)}
    y_idx = {g: i for i, g in enumerate(cats_y)}
    df["_x"] = df["_col"].map(x_idx)
    df["_y"] = df["gene"].map(y_idx)
    df["_size"] = -np.log10(df["padj"].clip(lower=1e-300))
    df["_size"] = (df["_size"] / df["_size"].max() * 200).clip(lower=8)

    fig_w = max(6, 0.5 * len(cats_x) + 4)
    fig_h = max(4, 0.22 * len(cats_y) + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sc = ax.scatter(df["_x"], df["_y"], c=df["log2FC"], cmap="RdBu_r",
                    s=df["_size"], vmin=-3, vmax=3, edgecolor="0.3", lw=0.3)
    ax.set_xticks(range(len(cats_x)))
    ax.set_xticklabels(cats_x, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(cats_y)))
    ax.set_yticklabels(cats_y, fontsize=7)
    ax.set_xlim(-0.5, len(cats_x) - 0.5)
    ax.set_ylim(-0.5, len(cats_y) - 0.5)
    cb = fig.colorbar(sc, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("log2 fold change", fontsize=7)
    ax.set_title(title, fontsize=10)
    fig.text(0.5, -0.02,
             f"Size = -log10(padj). Top {len(cats_y)} union genes.",
             ha="center", fontsize=6, style="italic")
    safe_fig(fig, out)


# ---------------------------------------------------------------------------
# Plot 5 — RRHO-lite cross-contrast scatter (Early vs Late log2FC)
# ---------------------------------------------------------------------------

def plot_rrho_scatter(df_master, out_dir, tissue, n_labels_both=8,
                     n_labels_single=3, use_blocklist=True):
    """For each matched (sex × cell type × age × level), scatter log2FC of
    Early-vs-Relaxed vs Late-vs-Relaxed. Quadrants reveal whether the two
    stress timings perturb genes concordantly (top-right + bottom-left) or
    divergently (top-left + bottom-right). Labels:
      - up to n_labels_both genes that are sig in BOTH contrasts (the
        publication-grade story — convergent stress signature)
      - up to n_labels_single genes per single-contrast set, picked by
        largest |log2FC|"""
    df = df_master[df_master["test_method"] == "Wald"]
    df = df[df["contrast"].isin(("early_vs_relaxed_per_age",
                                 "late_vs_relaxed_per_age"))]
    if df.empty:
        return
    try:
        from adjustText import adjust_text
        have_adjust_text = True
    except ImportError:
        have_adjust_text = False

    keys = ["sex", "celltype", "group_level", "level"]
    grouped = df.groupby(keys + ["contrast"], observed=True)

    slices = {}
    for (sex_label, ct, age, level, cname), g in grouped:
        k = (sex_label, ct, age, level)
        slices.setdefault(k, {})[cname] = g.set_index("gene")[["log2FC", "padj"]]
    drawn = 0
    for (sex_label, ct, age, level), d in slices.items():
        if ("early_vs_relaxed_per_age" not in d
                or "late_vs_relaxed_per_age" not in d):
            continue
        e = d["early_vs_relaxed_per_age"]
        l = d["late_vs_relaxed_per_age"]
        merged = e.join(l, lsuffix="_E", rsuffix="_L", how="inner").dropna(
                  subset=["log2FC_E", "log2FC_L"])
        if len(merged) < 30:
            continue
        from scipy.stats import spearmanr
        rho, _ = spearmanr(merged["log2FC_E"], merged["log2FC_L"])
        sig_E = (merged["padj_E"] < PADJ_THR) & (merged["log2FC_E"].abs() > LFC_THR)
        sig_L = (merged["padj_L"] < PADJ_THR) & (merged["log2FC_L"].abs() > LFC_THR)
        sig_both = sig_E & sig_L

        fig, ax = plt.subplots(figsize=(6.5, 6.0))
        ax.scatter(merged["log2FC_E"], merged["log2FC_L"],
                   s=3, color="lightgray", rasterized=True, alpha=0.5)
        ax.scatter(merged.loc[sig_E & ~sig_L, "log2FC_E"],
                   merged.loc[sig_E & ~sig_L, "log2FC_L"],
                   s=10, color="#c0392b", edgecolor="0.3", lw=0.2,
                   label=f"Early sig ({int((sig_E & ~sig_L).sum())})", alpha=0.75)
        ax.scatter(merged.loc[sig_L & ~sig_E, "log2FC_E"],
                   merged.loc[sig_L & ~sig_E, "log2FC_L"],
                   s=10, color="#2980b9", edgecolor="0.3", lw=0.2,
                   label=f"Late sig ({int((sig_L & ~sig_E).sum())})", alpha=0.75)
        ax.scatter(merged.loc[sig_both, "log2FC_E"],
                   merged.loc[sig_both, "log2FC_L"],
                   s=22, color="black", edgecolor="white", lw=0.4,
                   label=f"Both sig ({int(sig_both.sum())})", alpha=0.95)
        ax.axhline(0, color="k", lw=0.5)
        ax.axvline(0, color="k", lw=0.5)

        # Label the both-sig genes (the convergent stress signature) first
        texts = []
        if int(sig_both.sum()) > 0:
            both_df = merged.loc[sig_both].copy()
            both_df["_score"] = (both_df["log2FC_E"].abs()
                                 + both_df["log2FC_L"].abs())
            if use_blocklist:
                both_df = both_df[~both_df.index.to_series().apply(is_blocklisted)]
            for gene, row in both_df.nlargest(n_labels_both, "_score").iterrows():
                texts.append(ax.text(row["log2FC_E"], row["log2FC_L"], str(gene),
                                     fontsize=7, ha="left", va="bottom",
                                     color="black", fontweight="bold",
                                     bbox=dict(boxstyle="round,pad=0.15",
                                               fc="white", ec="0.5", lw=0.3,
                                               alpha=0.85)))
        # Then a few standout single-contrast genes (top |log2FC|)
        for mask, col in [(sig_E & ~sig_L, "#c0392b"),
                          (sig_L & ~sig_E, "#2980b9")]:
            if int(mask.sum()) == 0:
                continue
            sub_df = merged.loc[mask].copy()
            sub_df["_score"] = (sub_df["log2FC_E"].abs()
                                + sub_df["log2FC_L"].abs())
            if use_blocklist:
                sub_df = sub_df[~sub_df.index.to_series().apply(is_blocklisted)]
            for gene, row in sub_df.nlargest(n_labels_single, "_score").iterrows():
                texts.append(ax.text(row["log2FC_E"], row["log2FC_L"], str(gene),
                                     fontsize=6, ha="left", va="bottom",
                                     color=col, alpha=0.85,
                                     bbox=dict(boxstyle="round,pad=0.1",
                                               fc="white", ec="none", alpha=0.6)))
        if have_adjust_text and len(texts) > 1:
            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle="-", color="0.6", lw=0.3),
                        expand_points=(1.3, 1.3),
                        expand_text=(1.2, 1.3))

        ax.set_xlabel("log2FC Early vs Relaxed")
        ax.set_ylabel("log2FC Late vs Relaxed")
        ax.set_title(f"{tissue} | {ct} | {age} | {level} | sex={sex_label}\n"
                     f"Spearman ρ={rho:.2f}  (padj<{PADJ_THR} & "
                     f"|log2FC|>{LFC_THR})", fontsize=9)
        ax.legend(fontsize=7, frameon=False, loc="lower right")
        ax.spines[["top", "right"]].set_visible(False)
        out = (out_dir / sex_label / slugify(ct) / age / slugify(level)
               / "rrho_scatter.png")
        safe_fig(fig, out)
        drawn += 1
    print(f"    rrho: drew {drawn} scatter(s) -> {out_dir}"
          + (" (with adjustText)" if have_adjust_text
             else " (adjustText not installed — labels may overlap)"))


# ---------------------------------------------------------------------------
# Plot 6 — Top-gene dotplot (Scanpy-style; uses expression-matrix CSV)
# ---------------------------------------------------------------------------

def plot_top_gene_dotplot(df_sig_master, cfg, out_dir, tissue, subcluster,
                         n_top=20, use_blocklist=True):
    """For each (sex × contrast × cell type), pick top n_top sig genes by
    padj. Read the per-sample expression matrix written by 08b_de.py and
    draw a Scanpy-style dotplot: rows = group × age, cols = gene, color =
    mean lognorm, size = n_cells.

    Falls back to a heatmap if all genes have the same n_cells (no size
    variability) — common when one slice dominates.
    """
    if df_sig_master.empty:
        return
    table_dir = phase_table_dir(cfg, "08b_de")
    suffix = f"_subcluster_{subcluster}" if subcluster else ""
    expr_path = table_dir / f"08b_de_gene_expression_per_sample{suffix}.csv"
    meta_path = table_dir / f"08b_sample_metadata{suffix}.csv"
    if not expr_path.is_file() or not meta_path.is_file():
        print(f"    [dotplot skip] expression matrix not found at "
              f"{expr_path.name} / {meta_path.name}.")
        return
    expr = pd.read_csv(expr_path)
    meta = pd.read_csv(meta_path)
    if "sample_id" not in meta.columns and meta.columns[0] != "sample_id":
        meta = meta.rename(columns={meta.columns[0]: "sample_id"})
    meta = meta.set_index("sample_id")

    for (sex_label, cname), g in df_sig_master.groupby(["sex", "contrast"],
                                                       observed=True):
        if not str(cname).startswith(("early_vs_relaxed", "late_vs_relaxed",
                                       "early_vs_late")):
            continue
        # Pick the top sig genes for THIS contrast, per cell type
        g_sorted = g.sort_values("padj")
        if use_blocklist:
            g_sorted = filter_blocklist(g_sorted)
        per_ct = g_sorted.groupby("celltype", observed=True).head(n_top)
        if per_ct.empty:
            continue
        for ct, gg in per_ct.groupby("celltype", observed=True):
            genes = gg["gene"].tolist()
            sub = expr[(expr["celltype"] == ct) & (expr["gene"].isin(genes))]
            if sub.empty:
                continue
            # Restrict samples to this sex stratum (if relevant)
            if sex_label != "combined" and "sex" in meta.columns:
                samples_in_sex = set(meta.index[meta["sex"] == sex_label])
                sub = sub[sub["sample_id"].astype(str).isin(samples_in_sex)]
                if sub.empty:
                    continue
            # Build (group × age) -> gene mean
            sub = sub.merge(meta, left_on="sample_id", right_index=True,
                            how="left")
            if "group" not in sub.columns or "age" not in sub.columns:
                continue
            sub["_row"] = sub["group"].astype(str) + " | " + sub["age"].astype(str)
            mean_p = sub.pivot_table(index="_row", columns="gene",
                                     values="mean_lognorm", aggfunc="mean")
            ncell_p = sub.pivot_table(index="_row", columns="gene",
                                      values="n_cells", aggfunc="mean")
            # Order rows: GROUP_ORDER × AGE_ORDER
            ages = AGE_ORDER.get(tissue, sorted(sub["age"].astype(str).unique()))
            rows = [f"{gp} | {ag}" for gp in GROUP_ORDER for ag in ages
                    if f"{gp} | {ag}" in mean_p.index]
            mean_p = mean_p.reindex(index=rows, columns=genes)
            ncell_p = ncell_p.reindex(index=rows, columns=genes)
            if mean_p.empty:
                continue
            # Draw dotplot
            nr, nc = mean_p.shape
            fig, ax = plt.subplots(figsize=(max(6, 0.42 * nc + 2),
                                            max(3, 0.35 * nr + 1.5)))
            sizes = ncell_p.values.astype(float)
            # Normalise sizes to a sensible range
            smax = np.nanmax(sizes) if np.isfinite(np.nanmax(sizes)) else 1.0
            s_norm = (np.clip(sizes / max(smax, 1), 0, 1) * 220 + 12)
            xs, ys = np.meshgrid(np.arange(nc), np.arange(nr))
            valid = ~np.isnan(mean_p.values)
            sc = ax.scatter(xs[valid], ys[valid],
                            c=mean_p.values[valid], cmap="viridis",
                            s=s_norm[valid], edgecolor="0.3", lw=0.25)
            ax.set_xticks(range(nc))
            ax.set_xticklabels(mean_p.columns, rotation=45, ha="right", fontsize=7)
            ax.set_yticks(range(nr))
            ax.set_yticklabels(mean_p.index, fontsize=7)
            ax.set_xlim(-0.5, nc - 0.5)
            ax.set_ylim(-0.5, nr - 0.5)
            cb = fig.colorbar(sc, ax=ax, fraction=0.025, pad=0.02)
            cb.set_label("mean lognorm", fontsize=7)
            ax.set_title(f"{tissue} | {cname} | {ct} | sex={sex_label}\n"
                         f"top {len(genes)} DEGs (size = n_cells)",
                         fontsize=9)
            out = (out_dir / sex_label / slugify(cname) / slugify(ct)
                   / "top_gene_dotplot.png")
            safe_fig(fig, out)


# ---------------------------------------------------------------------------
# Plot 7 — Per-celltype volcano grid (small multiples)
# ---------------------------------------------------------------------------

def plot_volcano_grid(df_master_wald, out_dir, tissue, n_labels=5,
                     use_blocklist=True):
    """For each (sex × cell type), draw a grid of volcanoes:
       rows = age, cols = the stress contrasts at that age.
    Lets you compare signal magnitude across ages / contrasts at a glance.
    Labels the top `n_labels` significant genes per panel (by smallest padj,
    with ties broken by larger |log2FC|)."""
    df = df_master_wald[df_master_wald["contrast"].isin(STRESS_CONTRASTS)]
    if df.empty:
        return
    ages = ordered(df["group_level"].astype(str).unique(), AGE_ORDER.get(tissue, []))
    contrasts = ordered(df["contrast"].astype(str).unique(), list(STRESS_CONTRASTS))

    # Optional non-overlapping label placement
    try:
        from adjustText import adjust_text
        have_adjust_text = True
    except ImportError:
        have_adjust_text = False

    drawn = 0
    for (sex_label, ct), g in df.groupby(["sex", "celltype"], observed=True):
        nrow, ncol = len(ages), len(contrasts)
        if nrow == 0 or ncol == 0:
            continue
        fig, axes = plt.subplots(nrow, ncol,
                                 figsize=(max(6, 3.0 * ncol),
                                          max(4, 2.8 * nrow)),
                                 squeeze=False)
        for i, age in enumerate(ages):
            for j, cname in enumerate(contrasts):
                ax = axes[i, j]
                sub = g[(g["group_level"].astype(str) == age)
                        & (g["contrast"] == cname)]
                if sub.empty:
                    ax.text(0.5, 0.5, "—", ha="center", va="center",
                            transform=ax.transAxes, color="0.7")
                    ax.set_xticks([]); ax.set_yticks([])
                    if i == 0:
                        ax.set_title(cname.replace("_per_age", ""), fontsize=8)
                    if j == 0:
                        ax.set_ylabel(age, fontsize=8)
                    continue
                sub = sub.dropna(subset=["log2FC", "padj"])
                if sub.empty:
                    continue
                nlp = -np.log10(sub["padj"].clip(lower=1e-300))
                sig_mask = (sub["padj"] < PADJ_THR) & (sub["log2FC"].abs() > LFC_THR)
                ax.scatter(sub.loc[~sig_mask, "log2FC"], nlp[~sig_mask],
                           s=3, color="lightgray", rasterized=True)
                ax.scatter(sub.loc[sig_mask, "log2FC"], nlp[sig_mask],
                           s=8, color="salmon", edgecolor="0.3", lw=0.2,
                           rasterized=True)
                ax.axhline(-np.log10(PADJ_THR), color="k", lw=0.4, ls="--")
                ax.axvline(LFC_THR, color="k", lw=0.4, ls="--")
                ax.axvline(-LFC_THR, color="k", lw=0.4, ls="--")

                # Label top sig genes
                sig_rows = sub[sig_mask].copy()
                if use_blocklist:
                    sig_rows = filter_blocklist(sig_rows)
                if len(sig_rows) > 0:
                    sig_rows["_score"] = (-np.log10(sig_rows["padj"].clip(lower=1e-300))
                                          + sig_rows["log2FC"].abs())
                    top = sig_rows.nlargest(n_labels, "_score")
                    texts = []
                    for _, row in top.iterrows():
                        nlp_pt = -np.log10(max(row["padj"], 1e-300))
                        t = ax.text(row["log2FC"], nlp_pt, str(row["gene"]),
                                    fontsize=6, ha="left", va="bottom",
                                    color="black",
                                    bbox=dict(boxstyle="round,pad=0.15",
                                              fc="white", ec="none", alpha=0.7))
                        texts.append(t)
                    if have_adjust_text and len(texts) > 1:
                        adjust_text(texts, ax=ax,
                                    arrowprops=dict(arrowstyle="-",
                                                    color="0.5", lw=0.3),
                                    expand_points=(1.2, 1.2),
                                    expand_text=(1.1, 1.2),
                                    only_move={"points": "y", "text": "xy"})

                ax.text(0.04, 0.96, f"{int(sig_mask.sum())} sig",
                        transform=ax.transAxes, fontsize=7, va="top",
                        color="0.4", fontweight="bold")
                if i == 0:
                    ax.set_title(cname.replace("_per_age", ""), fontsize=8)
                if j == 0:
                    ax.set_ylabel(age, fontsize=8)
                ax.tick_params(labelsize=6)
        fig.suptitle(f"{tissue} | {ct} | sex={sex_label}\n"
                     f"(thresholds padj<{PADJ_THR} & |log2FC|>{LFC_THR}; "
                     f"top {n_labels} sig genes labeled per panel)",
                     fontsize=10)
        out = out_dir / sex_label / f"{slugify(ct)}_volcano_grid.png"
        safe_fig(fig, out)
        drawn += 1
    print(f"    volcano-grid: drew {drawn} grid(s) -> {out_dir}"
          + (" (with adjustText)" if have_adjust_text
             else " (adjustText not installed — labels may overlap)"))


# ---------------------------------------------------------------------------
# Plot 8 — Venn diagrams (directional)
# ---------------------------------------------------------------------------

def _build_sig_index(df_wald):
    """Pre-index significant Wald genes by (contrast, sex, ct, level,
    group_level) -> {'up': set, 'down': set, 'all': set}. Replaces a
    full-DataFrame scan per lookup with O(1) dict access. With 600+ Venn
    jobs × 18 lookups each, this is the difference between minutes and
    hours."""
    sig = df_wald[(df_wald["test_method"] == "Wald")
                  & df_wald["padj"].notna()
                  & (df_wald["padj"] < PADJ_THR)
                  & df_wald["log2FC"].notna()
                  & (df_wald["log2FC"].abs() > LFC_THR)]
    if sig.empty:
        return {}
    grp_keys = ["contrast", "sex", "celltype", "level"]
    idx = {}
    # Use itertuples for speed; group by full key then split up/down inside
    sig = sig.assign(_grp=sig["group_level"].astype(str))
    for (cn, sx, ct, lv, age), g in sig.groupby(grp_keys + ["_grp"],
                                                observed=True):
        genes = g["gene"].dropna().astype(str)
        lfc   = g["log2FC"].values
        all_genes  = set(genes)
        up_genes   = set(genes[lfc > 0])
        down_genes = set(genes[lfc < 0])
        idx[(cn, sx, str(ct), str(lv), str(age))] = {
            "all": all_genes, "up": up_genes, "down": down_genes,
        }
    return idx


def _gene_sets_for_slice(sig_idx, contrast, sex_label, ct, level, age,
                        direction=None):
    """O(1) lookup. `sig_idx` is built by _build_sig_index. `direction` in
    {'up','down',None (= all sig)}."""
    key = (contrast, sex_label, str(ct), str(level), str(age))
    bucket = sig_idx.get(key)
    if not bucket:
        return set()
    return bucket["all"] if direction is None else bucket[direction]


def plot_venns(df_wald, out_dir, tissue, ct_key=None, n_jobs=1):
    """Directional Venns per (cell type × region × sex).

       (a) Per-age 2-way: Early ∩ Late at each age in AGE_ORDER[tissue]
           figure = one row of 3 sub-Venns (P1, 4W, 3mo) × 3 directions
                    (all-sig, up, down) = up to 9 sub-panels
       (b) Across-age 3-way: P1 ∩ 4W ∩ 3mo within Early-vs-Relaxed and
           within Late-vs-Relaxed (brain only — placenta has no cross-age)
           figure = one row of 2 sub-Venns × 3 directions = 6 sub-panels
    """
    try:
        from matplotlib_venn import venn2, venn3
    except ImportError:
        return  # quietly skip; main() already warned

    ages = AGE_ORDER.get(tissue, [])
    if not ages:
        return

    df_wald = df_wald[df_wald["test_method"] == "Wald"]
    # Determine contrasts present (handles brain & placenta with their
    # tissue-specific contrast names).
    cnames = set(df_wald["contrast"].unique())
    early_c = next((c for c in cnames if "early_vs_relaxed" in c), None)
    late_c  = next((c for c in cnames if "late_vs_relaxed"  in c), None)
    if not (early_c and late_c):
        return

    # Pre-build the gene-set index ONCE before fanning out. Each Venn job runs
    # ~18 set lookups; without this index they'd each do a full-DataFrame
    # boolean scan, turning a minutes-long phase into hours.
    print("    venn: pre-indexing significant Wald genes...")
    sig_idx = _build_sig_index(df_wald)
    print(f"    venn: indexed {len(sig_idx):,} (contrast, sex, ct, level, age) keys.")

    sexes  = ordered(df_wald["sex"].unique(), SEX_ORDER)
    levels = ordered(df_wald["level"].astype(str).unique(), ["whole"])
    cts    = sorted(df_wald["celltype"].astype(str).unique())
    directions = [("all", None), ("up", "up"), ("down", "down")]

    jobs = []
    for sex_label in sexes:
        for ct in cts:
            for level in levels:
                jobs.append(("per_age", sex_label, ct, level, early_c, late_c, ages))
                if len(ages) >= 3:
                    jobs.append(("cross_age", sex_label, ct, level,
                                 early_c, late_c, ages))

    def _draw(job):
        kind, sex_label, ct, level, early_c, late_c, ages = job
        if kind == "per_age":
            # 3 ages × 3 directions = up to 9 sub-Venns
            fig, axes = plt.subplots(len(directions), len(ages),
                                     figsize=(3.5 * len(ages),
                                              3.0 * len(directions)),
                                     squeeze=False)
            any_drawn = False
            for di, (dlab, dval) in enumerate(directions):
                for aj, age in enumerate(ages):
                    ax = axes[di, aj]
                    s_e = _gene_sets_for_slice(sig_idx, early_c, sex_label, ct,
                                              level, age, dval)
                    s_l = _gene_sets_for_slice(sig_idx, late_c, sex_label, ct,
                                              level, age, dval)
                    if not (s_e or s_l):
                        ax.text(0.5, 0.5, "—", ha="center", va="center",
                                transform=ax.transAxes, color="0.7")
                        ax.set_axis_off()
                    else:
                        any_drawn = True
                        venn2([s_e, s_l],
                              set_labels=(f"Early\n(n={len(s_e)})",
                                          f"Late\n(n={len(s_l)})"),
                              ax=ax,
                              set_colors=("#c0392b", "#2980b9"), alpha=0.55)
                    if di == 0:
                        ax.set_title(age, fontsize=10)
                    if aj == 0:
                        ax.set_ylabel(dlab, fontsize=9, rotation=0,
                                      ha="right", va="center", labelpad=22)
            if any_drawn:
                fig.suptitle(f"{tissue} | {ct} | level={level} | "
                             f"sex={sex_label}\nEarly ∩ Late per age "
                             f"(padj<{PADJ_THR} & |log2FC|>{LFC_THR})",
                             fontsize=10)
                out = (out_dir / sex_label / slugify(level) / slugify(ct)
                       / "venn_early_vs_late_per_age.png")
                safe_fig(fig, out)
            else:
                plt.close(fig)
        elif kind == "cross_age":
            cnames = [("Early-vs-Relaxed", early_c),
                      ("Late-vs-Relaxed",  late_c)]
            fig, axes = plt.subplots(len(directions), len(cnames),
                                     figsize=(4.5 * len(cnames),
                                              3.5 * len(directions)),
                                     squeeze=False)
            any_drawn = False
            for di, (dlab, dval) in enumerate(directions):
                for cj, (clabel, cname) in enumerate(cnames):
                    ax = axes[di, cj]
                    sets = [
                        _gene_sets_for_slice(sig_idx, cname, sex_label, ct,
                                            level, age, dval)
                        for age in ages
                    ]
                    if not any(sets):
                        ax.text(0.5, 0.5, "—", ha="center", va="center",
                                transform=ax.transAxes, color="0.7")
                        ax.set_axis_off()
                    else:
                        any_drawn = True
                        venn3(sets,
                              set_labels=tuple(f"{a}\n(n={len(s)})"
                                               for a, s in zip(ages, sets)),
                              ax=ax,
                              set_colors=("#27ae60", "#e67e22", "#8e44ad"),
                              alpha=0.55)
                    if di == 0:
                        ax.set_title(clabel, fontsize=10)
                    if cj == 0:
                        ax.set_ylabel(dlab, fontsize=9, rotation=0,
                                      ha="right", va="center", labelpad=22)
            if any_drawn:
                fig.suptitle(f"{tissue} | {ct} | level={level} | "
                             f"sex={sex_label}\nAcross-age overlap per "
                             f"contrast (padj<{PADJ_THR} & |log2FC|>{LFC_THR})",
                             fontsize=10)
                out = (out_dir / sex_label / slugify(level) / slugify(ct)
                       / "venn_across_age_per_contrast.png")
                safe_fig(fig, out)
            else:
                plt.close(fig)

    print(f"    venn: queued {len(jobs)} jobs across {n_jobs} workers...")
    # plt is not thread-safe -> serialize (n_jobs=1) unless caller asked otherwise
    n_jobs_eff = max(1, n_jobs) if n_jobs <= 4 else 4   # plot draw is the bottleneck
    completed = 0
    for _job, _res, err in parallel_map(_draw, jobs, n_jobs=n_jobs_eff,
                                        use_threads=True, desc="venn"):
        if err:
            print(f"    [venn warn] {err[:120]}")
        else:
            completed += 1
    print(f"    venn: completed {completed}/{len(jobs)} -> {out_dir}")


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def parse_plots(s):
    if not s or s.lower() == "all":
        return ALL_PLOTS
    requested = [p.strip().lower() for p in s.split(",") if p.strip()]
    unknown = [p for p in requested if p not in ALL_PLOTS]
    if unknown:
        sys.exit(f"ERROR: unknown plot(s) {unknown}. Valid: {ALL_PLOTS}")
    return requested


def main():
    ap = argparse.ArgumentParser(description="Phase 8b summary plots.")
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--subcluster", default=None,
                    help="Mirror 08b_de.py: read the master CSV for this "
                         "subcluster slug (e.g. 'immune') and write summary "
                         "plots under plots/08b_de_subcluster_<slug>/summary/.")
    ap.add_argument("--plots", default="all",
                    help="Comma-separated subset of plot types to draw "
                         f"(default: all). Valid: {ALL_PLOTS}")
    ap.add_argument("--n-jobs", type=int, default=8,
                    help="Workers for the Venn parallel loop (capped at 4 "
                         "because matplotlib draw isn't fully thread-safe).")
    ap.add_argument("--max-genes-heatmap", type=int, default=50)
    ap.add_argument("--max-genes-bubble", type=int, default=30)
    ap.add_argument("--top-n-dotplot", type=int, default=20)
    ap.add_argument("--no-blocklist", action="store_true",
                    help="Disable the visualization gene blocklist (Hb / "
                         "sex-linked / mito). They'll be allowed back into "
                         "top-N selections for heatmaps, bubbles, dotplots, "
                         "volcano labels, and RRHO labels. The master CSV, "
                         "n_DEGs counts, UpSet sets, and Venn sets are never "
                         "filtered (data preserved). Use for QA.")
    args = ap.parse_args()
    use_blocklist = not args.no_blocklist

    print("\n=== Phase 8b summary plots ===")
    cfg = load_config(args.config)
    tissue = cfg.get("tissue")
    if tissue not in AGE_ORDER:
        print(f"  [note] unknown tissue '{tissue}'; falling back to sorted ages.")

    plots = parse_plots(args.plots)
    print(f"  Plots requested: {plots}")
    if use_blocklist:
        print(f"  Viz blocklist ACTIVE ({len(BLOCKLIST_FOR_VIZ)} genes + "
              f"prefixes {BLOCKLIST_PREFIXES}): "
              f"{sorted(BLOCKLIST_FOR_VIZ)[:5]}...  "
              f"(use --no-blocklist to disable for QA)")
    else:
        print(f"  Viz blocklist DISABLED — all genes eligible for top-N "
              f"selections.")

    # Soft availability check for optional libs
    try:
        import matplotlib_venn  # noqa: F401
    except ImportError:
        if "venn" in plots:
            print("  [warn] matplotlib_venn not installed; venn plots will be skipped.")
    try:
        import upsetplot  # noqa: F401
    except ImportError:
        if "upset" in plots:
            print("  [warn] upsetplot not installed; upset plots will be skipped.")

    df, csv_path, plot_root = read_master(cfg, args.subcluster)
    print(f"  Plot root: {plot_root}")

    # Filter to columns / rows we need most. Wald-only view used by most plots.
    df_wald = df[df["test_method"] == "Wald"].copy()
    sig_mask = is_sig_mask(df_wald)
    df_sig = df_wald[sig_mask].copy()
    print(f"  Wald rows: {len(df_wald):,}; sig: {len(df_sig):,} "
          f"(padj<{PADJ_THR} & |log2FC|>{LFC_THR})")

    if df_sig.empty and set(plots) - {"bar"}:
        print("  [note] no significant Wald rows — summary plots will mostly be empty.")

    # ---- 1. heatmap ------------------------------------------------------
    if "heatmap" in plots and not df_sig.empty:
        out_dir = plot_root / "heatmap"
        print("\n  [1/8] heatmaps...")
        for (sex_label, cname), g_sig in df_sig.groupby(["sex", "contrast"],
                                                       observed=True):
            g_all = df_wald[(df_wald["sex"] == sex_label)
                            & (df_wald["contrast"] == cname)]
            for level, g_lvl_sig in g_sig.groupby("level", observed=True):
                g_lvl_all = g_all[g_all["level"] == level]
                title = (f"{tissue} | {cname} | level={level} | "
                         f"sex={sex_label}")
                out = (out_dir / sex_label / slugify(cname)
                       / f"{slugify(level)}.png")
                plot_top_deg_heatmap(g_lvl_sig, g_lvl_all, title, out,
                                     max_genes=args.max_genes_heatmap,
                                     use_blocklist=use_blocklist)
        print(f"    -> {out_dir}")

    # ---- 2. upset --------------------------------------------------------
    if "upset" in plots and not df_sig.empty:
        out_dir = plot_root / "upset"
        print("\n  [2/8] upset plots...")
        for (sex_label, cname), g in df_sig.groupby(["sex", "contrast"],
                                                    observed=True):
            for level, g_lvl in g.groupby("level", observed=True):
                title = (f"{tissue} | {cname} | level={level} | sex={sex_label}")
                out = (out_dir / sex_label / slugify(cname)
                       / f"{slugify(level)}.png")
                plot_upset(g_lvl, title, out)
        print(f"    -> {out_dir}")

    # ---- 3. n_DEGs bar ---------------------------------------------------
    if "bar" in plots:
        out_dir = plot_root / "n_DEGs_bar"
        print("\n  [3/8] n_DEGs bar charts...")
        plot_n_degs_bar(df_sig, out_dir, tissue)
        print(f"    -> {out_dir}")

    # ---- 4. bubble -------------------------------------------------------
    if "bubble" in plots and not df_sig.empty:
        out_dir = plot_root / "bubble"
        print("\n  [4/8] bubble plots...")
        for (sex_label, cname), g_sig in df_sig.groupby(["sex", "contrast"],
                                                       observed=True):
            g_all = df_wald[(df_wald["sex"] == sex_label)
                            & (df_wald["contrast"] == cname)]
            for level, g_lvl_sig in g_sig.groupby("level", observed=True):
                g_lvl_all = g_all[g_all["level"] == level]
                title = (f"{tissue} | {cname} | level={level} | sex={sex_label}")
                out = (out_dir / sex_label / slugify(cname)
                       / f"{slugify(level)}.png")
                plot_bubble(g_lvl_all, g_lvl_sig, title, out,
                            max_genes=args.max_genes_bubble,
                            use_blocklist=use_blocklist)
        print(f"    -> {out_dir}")

    # ---- 5. rrho-lite scatter -------------------------------------------
    if "rrho" in plots:
        out_dir = plot_root / "rrho_scatter"
        print("\n  [5/8] cross-contrast (Early vs Late) scatter...")
        plot_rrho_scatter(df_wald, out_dir, tissue,
                          use_blocklist=use_blocklist)
        print(f"    -> {out_dir}")

    # ---- 6. dotplot ------------------------------------------------------
    if "dotplot" in plots and not df_sig.empty:
        out_dir = plot_root / "top_gene_dotplot"
        print("\n  [6/8] top-gene dotplots (reads expression matrix)...")
        plot_top_gene_dotplot(df_sig, cfg, out_dir, tissue, args.subcluster,
                              n_top=args.top_n_dotplot,
                              use_blocklist=use_blocklist)
        print(f"    -> {out_dir}")

    # ---- 7. volcano grid ------------------------------------------------
    if "grid" in plots:
        out_dir = plot_root / "volcano_grid"
        print("\n  [7/8] per-celltype volcano grids...")
        plot_volcano_grid(df_wald, out_dir, tissue,
                          use_blocklist=use_blocklist)
        print(f"    -> {out_dir}")

    # ---- 8. venn --------------------------------------------------------
    if "venn" in plots:
        out_dir = plot_root / "venn"
        print("\n  [8/8] Venn diagrams (directional)...")
        plot_venns(df_wald, out_dir, tissue, n_jobs=args.n_jobs)
        print(f"    -> {out_dir}")

    print(f"\n✓ Phase 8b summary plots written under {plot_root}")
    print(f"  (thresholds project-wide: padj<{PADJ_THR}, |log2FC|>{LFC_THR}; "
          f"Wald only — LRT rows are excluded because log2FC is NaN by design.)\n")


if __name__ == "__main__":
    main()
