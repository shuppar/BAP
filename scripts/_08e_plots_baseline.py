"""
_08e_plots_baseline.py — baseline CCC plot functions for Phase 8e.

Imported by 08e_communication.py. Not a standalone entry point.

All functions take baseline_df (08e_lr_baseline.csv) as input and are
fully reproducible offline from that CSV.

Functions:
  plot_chord_diagram              — single group×age chord
  plot_delta_chord_diagram        — all groups side-by-side for one age
  plot_network_graph              — per-sender network (reference image style)
  plot_interaction_count_heatmap  — source×target count matrix
  plot_baseline_dotplot           — top LR × ct_pair scatter
  plot_large_lr_dotplot           — top N (50/100/200) dotplot
  plot_top_lr_per_celltype_pair   — bar chart for most-active pairs
  plot_interaction_counts_barplot — total counts per group×age
  plot_pathway_activity_heatmap   — LR pairs collapsed to pathways
  plot_lr_persistence_across_ages — LR activity trajectory P1→4W→3mo
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path


def _slug(s: str) -> str:
    return s.replace(" ", "_").replace("/", "-").replace(".", "")


def _lr_label(row) -> str:
    return f"{row['ligand_complex']}→{row['receptor_complex']}"


# ============================================================================
# Chord diagrams
# ============================================================================

def plot_chord_diagram(baseline_df, group, age, magnitude_cutoff, pdir):
    """Single chord diagram for one group×age."""
    try:
        from mpl_chord_diagram import chord_diagram
    except ImportError:
        return

    sub = baseline_df[
        (baseline_df["group"] == group) & (baseline_df["age"] == age) &
        (baseline_df["magnitude_rank"] < magnitude_cutoff)
    ]
    if sub.empty:
        return

    ct_pairs = sub.groupby(["source", "target"]).size().reset_index(name="n")
    cell_types = sorted(set(ct_pairs["source"]) | set(ct_pairs["target"]))
    if len(cell_types) < 2:
        return

    ct_idx = {c: i for i, c in enumerate(cell_types)}
    mat = np.zeros((len(cell_types), len(cell_types)))
    for _, row in ct_pairs.iterrows():
        mat[ct_idx[row["source"]], ct_idx[row["target"]]] += row["n"]

    fig, ax = plt.subplots(figsize=(8, 8))
    colors = plt.cm.tab20(np.linspace(0, 1, len(cell_types)))
    chord_diagram(mat, cell_types, ax=ax, colors=colors,
                  fontsize=9, rotate_names=True, gap=0.03,
                  use_gradient=True, chord_width=0.7)
    ax.set_title(f"{int(mat.sum())} L-R interactions\n{group} | {age}", fontsize=11, pad=20)
    out = pdir / f"chord_{_slug(group)}_{_slug(age)}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}")


def plot_delta_chord_diagram(baseline_df, age, magnitude_cutoff, pdir):
    """Side-by-side chord for all groups at one age — direct visual comparison."""
    try:
        from mpl_chord_diagram import chord_diagram
    except ImportError:
        return

    active = baseline_df[
        (baseline_df["age"] == age) &
        (baseline_df["magnitude_rank"] < magnitude_cutoff)
    ]
    groups = sorted(active["group"].unique())
    if len(groups) < 2:
        return

    cell_types = sorted(set(active["source"]) | set(active["target"]))
    if len(cell_types) < 2:
        return
    ct_idx = {c: i for i, c in enumerate(cell_types)}
    colors = plt.cm.tab20(np.linspace(0, 1, len(cell_types)))

    fig, axes = plt.subplots(1, len(groups), figsize=(8 * len(groups), 8))
    if len(groups) == 1:
        axes = [axes]

    for ax, grp in zip(axes, groups):
        grp_active = active[active["group"] == grp]
        ct_pairs = grp_active.groupby(["source", "target"]).size().reset_index(name="n")
        mat = np.zeros((len(cell_types), len(cell_types)))
        for _, row in ct_pairs.iterrows():
            mat[ct_idx[row["source"]], ct_idx[row["target"]]] += row["n"]
        chord_diagram(mat, cell_types, ax=ax, colors=colors,
                      fontsize=8, rotate_names=True, gap=0.03,
                      use_gradient=True, chord_width=0.7)
        ax.set_title(f"{int(mat.sum())} interactions\n{grp} | {age}",
                     fontsize=10, pad=20)

    fig.suptitle(f"CCC comparison across groups: {age}", fontsize=12)
    out = pdir / f"chord_comparison_{_slug(age)}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}")


# ============================================================================
# Network graph (per-sender style, like reference image)
# ============================================================================

def plot_network_graph(baseline_df, group, age, magnitude_cutoff, pdir, level=None):
    """Network graph: one panel per source cell type, edges to targets.
    Edge width = n active LR pairs. Matches reference image style.
    If `level` is given, the slice is restricted to that level (e.g. a brain region
    or 'whole'); the level is appended to the title + filename."""
    try:
        import networkx as nx
    except ImportError:
        return

    active = baseline_df[
        (baseline_df["group"] == group) & (baseline_df["age"] == age) &
        (baseline_df["magnitude_rank"] < magnitude_cutoff)
    ]
    if level is not None and "level" in active.columns:
        active = active[active["level"].astype(str) == str(level)]
    if active.empty:
        return

    edges = (active.groupby(["source", "target"])
             .size().reset_index(name="n_lr"))
    cell_types = sorted(set(edges["source"]) | set(edges["target"]))
    if len(cell_types) < 2:
        return

    G_full = nx.DiGraph()
    G_full.add_nodes_from(cell_types)
    pos = nx.circular_layout(G_full)
    sources = sorted(edges["source"].unique())
    max_n = edges["n_lr"].max()

    ncols = min(3, len(sources))
    nrows = int(np.ceil(len(sources) / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(5 * ncols, 4.5 * nrows),
                             constrained_layout=True)
    axes_flat = np.array(axes).flatten() if len(sources) > 1 else [axes]
    cmap = plt.cm.tab10
    source_colors = {s: cmap(i / max(len(sources) - 1, 1))
                     for i, s in enumerate(sources)}

    for ax, src in zip(axes_flat, sources):
        src_edges = edges[edges["source"] == src]
        G = nx.DiGraph()
        G.add_nodes_from(cell_types)
        for _, row in src_edges.iterrows():
            G.add_edge(row["source"], row["target"], weight=row["n_lr"])
        color = source_colors[src]
        edge_widths = [G[u][v]["weight"] / max_n * 4 for u, v in G.edges()]
        nx.draw_networkx_nodes(G, pos, ax=ax, node_size=300,
                               node_color=["white"] * len(cell_types),
                               edgecolors="gray", linewidths=0.8)
        nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=[src],
                               node_size=500, node_color=[color], alpha=0.9)
        nx.draw_networkx_labels(G, pos, ax=ax, font_size=7)
        nx.draw_networkx_edges(G, pos, ax=ax, width=edge_widths,
                               edge_color=[color] * len(G.edges()),
                               arrows=True, arrowsize=12,
                               connectionstyle="arc3,rad=0.1", alpha=0.8)
        ax.set_title(src, fontsize=9, color=color, fontweight="bold")
        ax.axis("off")

    for ax in axes_flat[len(sources):]:
        ax.set_visible(False)

    lvl_tag = f" | level={level}" if level is not None else ""
    fig.suptitle(f"LR network: {group} | {age}{lvl_tag}\n"
                 f"(magnitude_rank<{magnitude_cutoff}; edge width=n active LR pairs)",
                 fontsize=10)
    lvl_slug = f"_{_slug(level)}" if level is not None else ""
    out = pdir / f"network_graph_{_slug(group)}_{_slug(age)}{lvl_slug}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Plot: {out.name}")


def _drop_unassigned(df):
    """Remove unassigned_* / contamination pseudo-types from source+target."""
    if df.empty:
        return df
    bad = df["source"].astype(str).str.startswith(("unassigned", "Contamination", "unresolved")) | \
          df["target"].astype(str).str.startswith(("unassigned", "Contamination", "unresolved"))
    return df[~bad]


def _draw_focal_fan_grid(edges, cell_types, sources, value_col, sign_col,
                         title, out_png, sqrt_scale=True):
    """Shared focal-fan grid renderer. edges has columns source,target,<value_col>,
    <sign_col>. width = sqrt(value/max) (field-standard for network edges);
    colour = red(+)/blue(−) from sign_col."""
    try:
        import networkx as nx
    except ImportError:
        return False
    if edges.empty or len(cell_types) < 2:
        return False
    G_full = nx.DiGraph(); G_full.add_nodes_from(cell_types)
    pos = nx.circular_layout(G_full)
    vmax = float(edges[value_col].abs().max()) or 1.0

    ncols = min(3, len(sources))
    nrows = int(np.ceil(len(sources) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows),
                             constrained_layout=True)
    axes_flat = np.array(axes).flatten() if len(sources) > 1 else [axes]

    for ax, src in zip(axes_flat, sources):
        se = edges[edges["source"] == src]
        G = nx.DiGraph(); G.add_nodes_from(cell_types)
        for _, r in se.iterrows():
            G.add_edge(r["source"], r["target"], v=r[value_col], s=r[sign_col])
        def _w(val):
            frac = abs(val) / vmax
            frac = np.sqrt(frac) if sqrt_scale else frac
            return frac * 5.5 + 0.4
        widths = [_w(G[u][v]["v"]) for u, v in G.edges()]
        ecolors = ["#d73027" if G[u][v]["s"] > 0 else "#4575b4" for u, v in G.edges()]
        nx.draw_networkx_nodes(G, pos, ax=ax, node_size=300,
                               node_color=["white"] * len(cell_types),
                               edgecolors="gray", linewidths=0.8)
        nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=[src],
                               node_size=520, node_color="0.3", alpha=0.9)
        nx.draw_networkx_labels(G, pos, ax=ax, font_size=7)
        nx.draw_networkx_edges(G, pos, ax=ax, width=widths, edge_color=ecolors,
                               arrows=True, arrowsize=12,
                               connectionstyle="arc3,rad=0.1", alpha=0.85)
        ax.set_title(src, fontsize=9, fontweight="bold")
        ax.axis("off")
    for ax in axes_flat[len(sources):]:
        ax.set_visible(False)

    fig.suptitle(title, fontsize=10)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return True


def _baseline_edge_metrics(b, test_group, ref_group, spec_fdr):
    """Per (source,target): n_changed (count of specificity-sig pairs with |Δ|>0)
    and sum_abs_delta (Σ|Δ score|), plus net sign (Σ Δ). Δ per pair =
    mean(test score) − mean(ref score), score = 1 − magnitude_rank."""
    b = b.copy()
    if "specificity_fdr" in b.columns:
        key = ["source", "target", "ligand_complex", "receptor_complex"]
        b["_ok"] = b["specificity_fdr"] <= spec_fdr
        b = b[b.groupby(key)["_ok"].transform("any")]
    if b.empty:
        return pd.DataFrame()
    b["score"] = 1.0 - b["magnitude_rank"].astype(float)
    pp = b.pivot_table(index=["source", "target", "ligand_complex", "receptor_complex"],
                       columns="group", values="score", aggfunc="mean")
    if test_group not in pp.columns or ref_group not in pp.columns:
        return pd.DataFrame()
    pp["d"] = pp[test_group].fillna(0) - pp[ref_group].fillna(0)
    pp = pp[pp["d"].abs() > 1e-9]
    if pp.empty:
        return pd.DataFrame()
    g = pp.reset_index().groupby(["source", "target"])
    out = pd.DataFrame({
        "n_changed": g.size(),
        "sum_abs_delta": g["d"].apply(lambda x: np.abs(x).sum()),
        "net_delta": g["d"].sum(),
    }).reset_index()
    return out


def _differential_edge_metrics(d, contrast_name, age, fdr):
    """Per (source,target): n_sig (FDR<fdr LR pairs), sum_abs_stat (Σ|interaction_stat|),
    net_stat (Σ interaction_stat for sign)."""
    d = d[(d["contrast_name"] == contrast_name) & (d["age"].astype(str) == str(age))].copy()
    d = d.dropna(subset=["interaction_stat", "interaction_padj"])
    d = d[d["interaction_padj"] < fdr]
    if d.empty:
        return pd.DataFrame()
    g = d.groupby(["source", "target"])
    out = pd.DataFrame({
        "n_changed": g.size(),
        "sum_abs_stat": g["interaction_stat"].apply(lambda x: np.abs(x).sum()),
        "net_stat": g["interaction_stat"].sum(),
    }).reset_index()
    return out


def plot_delta_network_grid(baseline_df, test_group, ref_group, age, magnitude_cutoff,
                            pdir, spec_fdr=0.05, level=None,
                            metric="count", arm="baseline", differential_df=None,
                            contrast_name=None, fdr=0.05):
    """Δ focal-fan grid (one panel per source cell type), sqrt-scaled edge widths,
    unassigned_* dropped. Four variants:

      arm='baseline'     metric='count'      → width=√(# changed specificity-sig pairs)
      arm='baseline'     metric='magnitude'  → width=√(Σ|Δ score|)
      arm='differential' metric='count'      → width=√(# FDR<fdr LR pairs)   [placenta]
      arm='differential' metric='magnitude'  → width=√(Σ|interaction_stat|)  [placenta]

    Colour = net direction (red=up / blue=down in stress). Baseline = DESCRIPTIVE
    (pooled cells); differential = FDR-backed."""
    if arm == "baseline":
        b = baseline_df[(baseline_df["age"] == age) &
                        (baseline_df["group"].isin([test_group, ref_group]))].copy()
        if level is not None and "level" in b.columns:
            b = b[b["level"].astype(str) == str(level)]
        b = _drop_unassigned(b)
        if b.empty or "magnitude_rank" not in b.columns:
            return
        b = b[b["magnitude_rank"] < magnitude_cutoff]
        m = _baseline_edge_metrics(b, test_group, ref_group, spec_fdr)
        value_col = "n_changed" if metric == "count" else "sum_abs_delta"
        sign_col = "net_delta"
        arm_tag = "DESCRIPTIVE (baseline, pooled cells, specificity-filtered)"
        wlabel = "# changed sig pairs" if metric == "count" else "Σ|Δ score|"
    else:  # differential
        if differential_df is None or differential_df.empty or contrast_name is None:
            return
        m = _differential_edge_metrics(_drop_unassigned(differential_df),
                                       contrast_name, age, fdr)
        value_col = "n_changed" if metric == "count" else "sum_abs_stat"
        sign_col = "net_stat"
        arm_tag = f"FDR<{fdr} (differential arm)"
        wlabel = f"# sig pairs (FDR<{fdr})" if metric == "count" else "Σ|interaction_stat|"

    if m.empty:
        return
    cell_types = sorted(set(m["source"]) | set(m["target"]))
    sources = sorted(m["source"].unique())
    lvl_tag = f" | level={level}" if level is not None else ""
    title = (f"Δ LR network ({metric}): {test_group} − {ref_group} | {age}{lvl_tag}\n"
             f"edge width=√({wlabel}); red=up / blue=down in stress; {arm_tag}")
    lvl_slug = f"_{_slug(level)}" if level is not None else ""
    out = (pdir / f"delta_grid_{arm}_{metric}_"
                  f"{_slug(test_group)}_vs_{_slug(ref_group)}_{_slug(age)}{lvl_slug}.png")
    if _draw_focal_fan_grid(m, cell_types, sources, value_col, sign_col,
                            title, out, sqrt_scale=True):
        print(f"  Plot: {out.name}")


# ============================================================================
# Count heatmap
# ============================================================================

def plot_interaction_count_heatmap(baseline_df, group, age, magnitude_cutoff, pdir):
    """Source×target heatmap of n active LR pairs."""
    active = baseline_df[
        (baseline_df["group"] == group) & (baseline_df["age"] == age) &
        (baseline_df["magnitude_rank"] < magnitude_cutoff)
    ]
    if active.empty:
        return

    mat = active.groupby(["source", "target"]).size().unstack(fill_value=0)
    if mat.shape[0] < 2:
        return

    fig, ax = plt.subplots(figsize=(max(5, mat.shape[1] * 0.7 + 1),
                                    max(4, mat.shape[0] * 0.6 + 1)))
    im = ax.imshow(mat.values, cmap="YlOrRd", aspect="auto")
    plt.colorbar(im, ax=ax, label="n active LR pairs")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=8)
    ax.set_xlabel("Target"); ax.set_ylabel("Source")
    ax.set_title(f"Interaction count heatmap\n{group} | {age}")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            if v > 0:
                ax.text(j, i, str(v), ha="center", va="center", fontsize=7,
                        color="white" if v > mat.values.max() * 0.6 else "black")
    fig.tight_layout()
    out = pdir / f"interaction_count_heatmap_{_slug(group)}_{_slug(age)}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}")


# ============================================================================
# Dotplots
# ============================================================================

def plot_baseline_dotplot(baseline_df, group, age, top_n, pdir):
    """Standard dotplot: top LR pairs, size=magnitude, color=specificity."""
    sub = baseline_df[
        (baseline_df["group"] == group) & (baseline_df["age"] == age) &
        (baseline_df["magnitude_rank"] < 0.2)
    ].copy()
    if sub.empty:
        return

    sub["lr_pair"] = sub.apply(_lr_label, axis=1)
    sub["ct_pair"] = sub["source"] + " → " + sub["target"]
    top = sub.nsmallest(top_n, "magnitude_rank")
    if top.empty:
        return

    ct_pairs = sorted(top["ct_pair"].unique())
    lr_pairs = list(top.groupby("lr_pair")["magnitude_rank"].min().sort_values().index)
    ct_idx = {c: i for i, c in enumerate(ct_pairs)}
    lr_idx = {l: i for i, l in enumerate(lr_pairs)}

    fig, ax = plt.subplots(
        figsize=(max(8, len(ct_pairs) * 0.6 + 2),
                 max(6, len(lr_pairs) * 0.3 + 2)))
    sc = ax.scatter(
        top["ct_pair"].map(ct_idx), top["lr_pair"].map(lr_idx),
        s=(1 - top["magnitude_rank"]) * 300,
        c=top["specificity_rank"] if "specificity_rank" in top.columns else 0.5,
        cmap="RdYlBu_r", vmin=0, vmax=1, alpha=0.8)
    plt.colorbar(sc, ax=ax, label="specificity_rank")
    ax.set_xticks(range(len(ct_pairs)))
    ax.set_xticklabels(ct_pairs, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(lr_pairs)))
    ax.set_yticklabels(lr_pairs, fontsize=7)
    ax.set_title(f"Top {top_n} LR pairs\n{group} | {age}")
    fig.tight_layout()
    out = pdir / f"baseline_dotplot_{_slug(group)}_{_slug(age)}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}")


def plot_large_lr_dotplot(baseline_df, group, age, top_n, pdir):
    """Large dotplot (top 50/100/200) — supplementary sensemaking figure."""
    sub = baseline_df[
        (baseline_df["group"] == group) & (baseline_df["age"] == age)
    ].copy()
    if sub.empty or "magnitude_rank" not in sub.columns:
        return

    sub["lr_pair"] = sub.apply(_lr_label, axis=1)
    sub["ct_pair"] = sub["source"] + " → " + sub["target"]
    top = sub.nsmallest(top_n, "magnitude_rank").copy()
    if top.empty:
        return

    lr_order = list(top.groupby("lr_pair")["magnitude_rank"].min().sort_values().index)
    ct_order = sorted(top["ct_pair"].unique())
    lr_idx = {l: i for i, l in enumerate(lr_order)}
    ct_idx = {c: i for i, c in enumerate(ct_order)}

    fig, ax = plt.subplots(
        figsize=(max(6, len(ct_order) * 0.55 + 2),
                 max(8, len(lr_order) * 0.22 + 2)))
    sc = ax.scatter(
        top["ct_pair"].map(ct_idx), top["lr_pair"].map(lr_idx),
        s=(1 - top["magnitude_rank"]) * 250 + 10,
        c=top["specificity_rank"] if "specificity_rank" in top.columns else 0.5,
        cmap="RdYlBu_r", vmin=0, vmax=1, alpha=0.85, linewidths=0)
    plt.colorbar(sc, ax=ax, label="specificity_rank")
    ax.set_xticks(range(len(ct_order)))
    ax.set_xticklabels(ct_order, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(lr_order)))
    ax.set_yticklabels(lr_order, fontsize=6)
    ax.set_title(f"Top {top_n} LR pairs — {group} | {age}\n"
                 f"(size=1−magnitude_rank, color=specificity_rank)", fontsize=9)
    fig.tight_layout()
    out = pdir / f"lr_dotplot_top{top_n}_{_slug(group)}_{_slug(age)}.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}  ({top_n} LR pairs, {len(ct_order)} cell-type pairs)")


def plot_top_lr_per_celltype_pair(baseline_df, group, age, top_n, pdir):
    """Bar chart of top LR pairs for the 4 most-active source→target pairs."""
    active = baseline_df[
        (baseline_df["group"] == group) & (baseline_df["age"] == age) &
        (baseline_df["magnitude_rank"] < 0.2)
    ].copy()
    if active.empty:
        return

    active["lr_pair"] = active.apply(_lr_label, axis=1)
    top_pairs = (active.groupby(["source", "target"])
                 .size().nlargest(4).index.tolist())
    if not top_pairs:
        return

    ncols = min(2, len(top_pairs))
    nrows = int(np.ceil(len(top_pairs) / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(7 * ncols, max(4, top_n * 0.3 + 1) * nrows),
                             constrained_layout=True)
    axes_flat = np.array(axes).flatten() if len(top_pairs) > 1 else [axes]

    for ax, (src, tgt) in zip(axes_flat, top_pairs):
        pair_df = active[(active["source"] == src) & (active["target"] == tgt)]
        pair_df = pair_df.nsmallest(top_n, "magnitude_rank").copy()
        if pair_df.empty:
            ax.set_visible(False)
            continue
        pair_df = pair_df.sort_values("magnitude_rank")
        scores = 1 - pair_df["magnitude_rank"]
        colors = plt.cm.RdYlBu_r(
            pair_df["specificity_rank"].values
            if "specificity_rank" in pair_df.columns
            else np.linspace(0.2, 0.8, len(pair_df))
        )
        ax.barh(range(len(pair_df)), scores, color=colors)
        ax.set_yticks(range(len(pair_df)))
        ax.set_yticklabels(pair_df["lr_pair"].values, fontsize=7)
        ax.set_xlabel("1 − magnitude_rank")
        ax.set_title(f"{src} → {tgt}", fontsize=9, fontweight="bold")
        ax.invert_yaxis()

    for ax in axes_flat[len(top_pairs):]:
        ax.set_visible(False)

    fig.suptitle(f"Top {top_n} LR pairs per cell-type pair\n{group} | {age}", fontsize=10)
    out = pdir / f"top_lr_per_pair_{_slug(group)}_{_slug(age)}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Plot: {out.name}")


# ============================================================================
# Summary / global plots
# ============================================================================

def plot_interaction_counts_barplot(baseline_df, ref_group, pdir):
    """Bar chart: total active LR interactions per group×age."""
    if baseline_df.empty or "magnitude_rank" not in baseline_df.columns:
        return

    active = baseline_df[baseline_df["magnitude_rank"] < 0.05]
    if active.empty:
        return

    counts = active.groupby(["age", "group"]).size().reset_index(name="n")
    ages = sorted(counts["age"].unique())
    groups = sorted(counts["group"].unique())
    group_colors = {g: plt.cm.Set2(i / max(len(groups) - 1, 1))
                    for i, g in enumerate(groups)}

    import matplotlib.patches as mpatches
    fig, axes = plt.subplots(1, len(ages), figsize=(4 * len(ages), 4),
                             sharey=False, constrained_layout=True)
    if len(ages) == 1:
        axes = [axes]

    for ax, age in zip(axes, ages):
        sub = counts[counts["age"] == age]
        for grp in groups:
            row = sub[sub["group"] == grp]
            n = int(row["n"].values[0]) if not row.empty else 0
            ax.bar(grp, n, color=group_colors[grp])
        ax.set_title(age)
        ax.set_ylabel("# active LR pairs (magnitude_rank<0.05)")
        ax.tick_params(axis="x", rotation=30)

    handles = [mpatches.Patch(color=group_colors[g], label=g) for g in groups]
    fig.legend(handles=handles, loc="upper right", fontsize=8)
    fig.suptitle("Active LR interactions per group × age", fontsize=11)
    out = pdir / "interaction_counts_by_group_age.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Plot: {out.name}")


def plot_pathway_activity_heatmap(baseline_df, magnitude_cutoff, pdir):
    """Collapse LR pairs to named pathways — pathway × group×age heatmap."""
    if baseline_df.empty or "magnitude_rank" not in baseline_df.columns:
        return
    try:
        import liana as li
        res = li.rs.select_resource("mouseconsensus")
        pathway_cols = [c for c in res.columns
                        if any(k in c.lower() for k in ("pathway", "category", "family"))]
        if not pathway_cols:
            print("  [skip] pathway heatmap: no pathway column in liana resource")
            return
        pc = pathway_cols[0]
        l_col = next((c for c in res.columns if "ligand" in c.lower()), None)
        r_col = next((c for c in res.columns if "receptor" in c.lower()), None)
        if not l_col or not r_col:
            return
        pw_map = res[[l_col, r_col, pc]].copy()
        pw_map.columns = ["ligand_complex", "receptor_complex", "pathway"]
        pw_map = pw_map.dropna(subset=["pathway"]).drop_duplicates()
    except Exception:
        print("  [skip] pathway heatmap: could not fetch liana resource")
        return

    merged = baseline_df.merge(pw_map, on=["ligand_complex", "receptor_complex"], how="left")
    merged = merged.dropna(subset=["pathway"])
    active = merged[merged["magnitude_rank"] < magnitude_cutoff].copy()
    if active.empty:
        return

    active["activity"] = 1 - active["magnitude_rank"]
    active["group_age"] = active["group"] + "\n" + active["age"]
    pivot = (active.groupby(["pathway", "group_age"])["activity"]
             .mean().unstack(fill_value=0))
    pivot = pivot.loc[pivot.mean(axis=1).nlargest(30).index]
    if pivot.empty:
        return

    fig, ax = plt.subplots(figsize=(max(8, pivot.shape[1] * 0.9 + 2),
                                    max(6, pivot.shape[0] * 0.35 + 2)))
    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
    plt.colorbar(im, ax=ax, label="Mean LR activity (1−magnitude_rank)")
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, fontsize=8, rotation=30, ha="right")
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_title(f"Pathway-level LR activity (top 30 pathways)", fontsize=10)
    fig.tight_layout()
    out = pdir / "pathway_activity_heatmap.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}")


def plot_lr_persistence_across_ages(baseline_df, magnitude_cutoff, top_n, pdir):
    """LR activity trajectory across ages — the developmental arc of signaling."""
    if baseline_df.empty or "magnitude_rank" not in baseline_df.columns:
        return
    ages = sorted(baseline_df["age"].unique())
    if len(ages) < 2:
        return  # single-age dev data; silently skip

    active = baseline_df[baseline_df["magnitude_rank"] < magnitude_cutoff].copy()
    active["activity"] = 1 - active["magnitude_rank"]
    active["lr_pair"] = active.apply(_lr_label, axis=1)
    active["ct_pair"] = active["source"] + " → " + active["target"]

    top_lrs = (active.groupby("lr_pair")["activity"]
               .mean().nlargest(top_n).index.tolist())
    df = active[active["lr_pair"].isin(top_lrs)]
    ct_pairs = sorted(df["ct_pair"].unique())[:12]
    groups = sorted(df["group"].unique())
    group_colors = {g: plt.cm.Set1(i / max(len(groups) - 1, 1))
                    for i, g in enumerate(groups)}

    ncols = min(3, len(ct_pairs))
    nrows = int(np.ceil(len(ct_pairs) / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(6 * ncols, 4 * nrows),
                             constrained_layout=True)
    axes_flat = np.array(axes).flatten() if len(ct_pairs) > 1 else [axes]

    for ax, ct in zip(axes_flat, ct_pairs):
        sub = df[df["ct_pair"] == ct]
        for grp, color in group_colors.items():
            age_means = (sub[sub["group"] == grp]
                         .groupby("age")["activity"].mean().reindex(ages))
            ax.plot(ages, age_means.values, color=color, marker="o",
                    label=grp, linewidth=1.8)
        ax.set_title(ct, fontsize=8)
        ax.set_ylim(0, 1)
        ax.set_ylabel("Activity", fontsize=7)
        ax.tick_params(axis="x", labelsize=7, rotation=20)

    for ax in axes_flat[len(ct_pairs):]:
        ax.set_visible(False)

    handles = [plt.Line2D([0], [0], color=group_colors[g], marker="o", label=g)
               for g in groups]
    fig.legend(handles=handles, loc="lower right", fontsize=8)
    fig.suptitle(f"LR activity trajectory across ages (top {top_n} LR pairs)", fontsize=10)
    out = pdir / "lr_persistence_across_ages.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Plot: {out.name}")
