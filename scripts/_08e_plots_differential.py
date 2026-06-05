"""
_08e_plots_differential.py — differential CCC plot functions for Phase 8e.

Imported by 08e_communication.py. Not a standalone entry point.

Key addition: ALL plots cover all three group comparisons:
  ES-v-Relaxed, LS-v-Relaxed, ES-v-LS (not just stress-vs-reference).

Functions:
  plot_differential_dotplot        — top deregulated LR pairs per contrast×age
  plot_differential_volcano        — Wald stat vs expression level
  plot_stress_signature_heatmap    — top LR pairs × contrasts×ages (persistence)
  plot_delta_lr_heatmap            — Δ mean activity score: grp_a vs grp_b
                                     (the "delta color" heatmap requested)
  plot_sender_receiver_bubble      — x=sender, y=receiver, size=n LR pairs
  plot_sender_receiver_heatmap     — Δ sender/receiver vs ref per cell type
  plot_delta_sender_receiver_heatmap — any two groups comparison
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
# Differential LR plots (from df_to_lr / 08e_lr_differential.csv)
# ============================================================================

def plot_differential_dotplot(diff_df, contrast_name, age, top_n, pdir):
    """Dotplot of top differential LR pairs by |interaction_stat|.
    Covers all three contrasts (ES-v-Rel, LS-v-Rel, ES-v-LS) as called."""
    sub = diff_df[
        (diff_df["contrast_name"] == contrast_name) &
        (diff_df["age"] == age)
    ].copy()
    if sub.empty or "interaction_stat" not in sub.columns:
        return

    sub["lr_pair"] = sub.apply(_lr_label, axis=1)
    sub["ct_pair"] = sub["source"] + " → " + sub["target"]
    sub["abs_stat"] = sub["interaction_stat"].abs()
    top = sub.nlargest(top_n, "abs_stat")
    if top.empty:
        return

    ct_pairs = sorted(top["ct_pair"].unique())
    lr_pairs = list(top.groupby("lr_pair")["abs_stat"].max()
                    .sort_values(ascending=False).index)
    ct_idx = {c: i for i, c in enumerate(ct_pairs)}
    lr_idx = {l: i for i, l in enumerate(lr_pairs)}

    fig, ax = plt.subplots(
        figsize=(max(8, len(ct_pairs) * 0.6 + 2),
                 max(6, len(lr_pairs) * 0.3 + 2)))
    colors = top["interaction_stat"]
    max_abs = max(colors.abs().max(), 1e-6)
    sc = ax.scatter(
        top["ct_pair"].map(ct_idx), top["lr_pair"].map(lr_idx),
        s=150, c=colors, cmap="RdBu_r", vmin=-max_abs, vmax=max_abs, alpha=0.85)
    plt.colorbar(sc, ax=ax, label="interaction_stat (red=up in test group)")
    ax.set_xticks(range(len(ct_pairs)))
    ax.set_xticklabels(ct_pairs, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(lr_pairs)))
    ax.set_yticklabels(lr_pairs, fontsize=7)
    ax.set_title(f"Top {top_n} differential LR pairs\n{contrast_name} | {age}")
    fig.tight_layout()
    out = pdir / f"differential_dotplot_{_slug(contrast_name)}_{_slug(age)}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}")


def plot_differential_volcano(diff_df, contrast_name, age, top_n, pdir):
    """Volcano: x=interaction_stat, y=expression level (lr_means).
    High stat + high expression = biologically plausible signal."""
    sub = diff_df[
        (diff_df["contrast_name"] == contrast_name) &
        (diff_df["age"] == age)
    ].copy()
    if sub.empty or "interaction_stat" not in sub.columns:
        return

    has_expr = "lr_means" in sub.columns
    sub["lr_pair"] = sub.apply(_lr_label, axis=1)
    sub["ct_pair"] = sub["source"] + " → " + sub["target"]
    sub["label"] = sub["lr_pair"] + "\n[" + sub["ct_pair"] + "]"
    sub["abs_stat"] = sub["interaction_stat"].abs()

    fig, ax = plt.subplots(figsize=(9, 6))
    y_vals = np.log1p(sub["lr_means"]) if has_expr else sub["abs_stat"]
    y_label = "log1p(lr_means) — expression" if has_expr else "|interaction_stat|"
    colors = ["#d73027" if s > 0 else "#4575b4" for s in sub["interaction_stat"]]
    ax.scatter(sub["interaction_stat"], y_vals, c=colors, alpha=0.6, s=30, linewidths=0)

    for _, r in sub.nlargest(top_n, "abs_stat").iterrows():
        y = np.log1p(r["lr_means"]) if has_expr else r["abs_stat"]
        ax.annotate(r["label"], (r["interaction_stat"], y),
                    fontsize=5.5, ha="center", va="bottom",
                    xytext=(0, 3), textcoords="offset points",
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))

    ax.axvline(0, color="k", lw=0.8, ls="--")
    ax.set_xlabel("interaction_stat (positive = up in test group)")
    ax.set_ylabel(y_label)
    ax.set_title(f"Differential LR volcano\n{contrast_name} | {age}\n"
                 f"red=up, blue=down (top {top_n} labeled)")
    fig.tight_layout()
    out = pdir / f"differential_volcano_{_slug(contrast_name)}_{_slug(age)}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}")


def plot_stress_signature_heatmap(diff_df, tdir, pdir, top_n=100):
    """Top N differential LR pairs × contrast×age — persistence heatmap.

    Shows which LR changes are consistent across all 3 contrasts and ages
    (persistent signal) vs specific to one contrast or age.
    Saved as plot + offline CSV.
    """
    if diff_df.empty or "interaction_stat" not in diff_df.columns:
        return

    df = diff_df.copy()
    df["lr_pair"] = df["ligand_complex"] + "→" + df["receptor_complex"]
    df["ct_pair"] = df["source"] + "→" + df["target"]
    df["lr_ct"] = df["lr_pair"] + "  [" + df["ct_pair"] + "]"
    df["contrast_age"] = df["contrast_name"] + "\n" + df["age"]

    top_lrct = (df.groupby("lr_ct")["interaction_stat"]
                .apply(lambda x: x.abs().max())
                .nlargest(top_n).index)
    pivot = (df[df["lr_ct"].isin(top_lrct)]
             .groupby(["lr_ct", "contrast_age"])["interaction_stat"]
             .mean().unstack(fill_value=0)
             .loc[top_lrct])

    if pivot.empty:
        return

    max_abs = max(pivot.values.__abs__().max(), 1e-6)
    fig, ax = plt.subplots(figsize=(max(8, pivot.shape[1] * 1.2 + 2),
                                    max(8, pivot.shape[0] * 0.28 + 2)))
    im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r",
                   vmin=-max_abs, vmax=max_abs)
    plt.colorbar(im, ax=ax, label="interaction_stat (red=up in test)")
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, fontsize=8, rotation=30, ha="right")
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=7)
    ax.set_title(f"Top {top_n} differential LR pairs × contrast × age\n"
                 f"(all 3 comparisons: ES-v-Rel, LS-v-Rel, ES-v-LS; "
                 f"persistent = consistently red/blue across columns)", fontsize=9)
    fig.tight_layout()
    out = pdir / "stress_signature_heatmap.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}")

    pivot.to_csv(tdir / "08e_lr_stress_signature_pivot.csv")
    print(f"  Table: 08e_lr_stress_signature_pivot.csv")


# ============================================================================
# Delta LR activity heatmap (the "color = delta" heatmap requested)
# Reproducible from 08e_lr_baseline.csv
# ============================================================================

def plot_delta_lr_heatmap(baseline_df, grp_a, grp_b, age, magnitude_cutoff, pdir,
                          focus_celltypes=None, zscore_rows=False,
                          pathway_map=None, max_rows=80):
    """Clustered Δ LR activity heatmap: grp_a − grp_b at one age.

    Hierarchical clustering on both axes:
      - rows (LR pairs): correlation distance, average linkage (patterns)
      - cols (cell-type pairs): euclidean distance, average linkage (magnitudes)

    Annotations:
      - row colors: signaling pathway from liana mouseconsensus (if pathway_map given)
      - col colors: source cell type / target cell type (two bars)

    Filters:
      - focus_celltypes: list of cell types — keeps interactions where source
        OR target is in the list. If None, no filtering.
      - zscore_rows: if True, z-score each row (pattern view); else absolute Δ.

    Reproducible from 08e_lr_baseline.csv (+ optional liana resource for pathway).
    """
    try:
        import seaborn as sns
    except ImportError:
        print("  [skip] clustered heatmap: seaborn not installed")
        return

    if baseline_df.empty or "magnitude_rank" not in baseline_df.columns:
        return

    def _slice(grp):
        return baseline_df[
            (baseline_df["group"] == grp) & (baseline_df["age"] == age)
        ].copy()

    df_a = _slice(grp_a)
    df_b = _slice(grp_b)
    if df_a.empty or df_b.empty:
        return

    for df in (df_a, df_b):
        df["lr_pair"] = df.apply(_lr_label, axis=1)
        df["ct_pair"] = df["source"] + "→" + df["target"]
        df["activity"] = 1 - df["magnitude_rank"]

    # Apply focus filter (source OR target in list)
    if focus_celltypes:
        focus_set = set(focus_celltypes)
        df_a = df_a[df_a["source"].isin(focus_set) | df_a["target"].isin(focus_set)]
        df_b = df_b[df_b["source"].isin(focus_set) | df_b["target"].isin(focus_set)]
        if df_a.empty or df_b.empty:
            print(f"  [skip] delta_lr_heatmap {grp_a}_vs_{grp_b}_{age}_focused: "
                  f"no interactions after focus filter")
            return

    active_lrs = set(
        df_a.loc[df_a["magnitude_rank"] < magnitude_cutoff, "lr_pair"]
    ) | set(
        df_b.loc[df_b["magnitude_rank"] < magnitude_cutoff, "lr_pair"]
    )
    if not active_lrs:
        return

    def _mean_mat(df):
        return (df[df["lr_pair"].isin(active_lrs)]
                .groupby(["lr_pair", "ct_pair"])["activity"]
                .mean().unstack(fill_value=0))

    mat_a = _mean_mat(df_a)
    mat_b = _mean_mat(df_b)
    all_lrs = sorted(active_lrs)
    all_cts = sorted(set(mat_a.columns) | set(mat_b.columns))
    mat_a = mat_a.reindex(index=all_lrs, columns=all_cts, fill_value=0)
    mat_b = mat_b.reindex(index=all_lrs, columns=all_cts, fill_value=0)
    delta = mat_a - mat_b
    delta = delta.loc[(delta != 0).any(axis=1), (delta != 0).any(axis=0)]
    if delta.empty:
        return

    # Cap rows by |Δ| for readability while preserving clustering coherence
    top_rows = delta.abs().max(axis=1).nlargest(max_rows).index
    delta = delta.loc[top_rows]

    title_suffix = " (focused)" if focus_celltypes else ""
    out_suffix = "_focused" if focus_celltypes else ""

    # Z-score rows if requested (pattern view)
    plot_mat = delta.copy()
    if zscore_rows:
        row_means = plot_mat.mean(axis=1).values[:, None]
        row_sds = plot_mat.std(axis=1).values[:, None]
        row_sds[row_sds == 0] = 1
        plot_mat[:] = (plot_mat.values - row_means) / row_sds
        cbar_label = "Row-Z-scored Δ activity"
    else:
        cbar_label = f"Δ activity ({grp_a} − {grp_b})"

    # NaN-safe for clustering
    plot_mat = plot_mat.fillna(0)

    # Build row annotation: signaling pathway
    row_colors = None
    if pathway_map is not None and not pathway_map.empty:
        # pathway_map has columns: lr_pair, pathway (built upstream)
        pw_lookup = dict(zip(pathway_map["lr_pair"], pathway_map["pathway"]))
        row_pw = plot_mat.index.map(lambda lr: pw_lookup.get(lr.split("  [")[0], "unknown"))
        unique_pw = sorted(set(row_pw) - {"unknown"})
        # Cap at 15 top pathways; rest as "other"
        if len(unique_pw) > 15:
            pw_counts = pd.Series(row_pw).value_counts()
            top_pw = set(pw_counts.head(15).index) - {"unknown"}
            row_pw = [p if p in top_pw else "other" for p in row_pw]
            unique_pw = sorted(set(row_pw) - {"unknown", "other"})
        cmap_pw = plt.cm.tab20
        pw_color = {p: cmap_pw(i / max(len(unique_pw), 1))
                    for i, p in enumerate(unique_pw)}
        pw_color["unknown"] = "lightgray"
        pw_color["other"] = "darkgray"
        row_colors = pd.Series([pw_color[p] for p in row_pw],
                                index=plot_mat.index, name="pathway")

    # Build column annotation: source + target cell type
    col_colors = None
    sources = [c.split("→")[0] for c in plot_mat.columns]
    targets = [c.split("→")[1] for c in plot_mat.columns]
    all_ct_names = sorted(set(sources) | set(targets))
    cmap_ct = plt.cm.tab20
    ct_color = {ct: cmap_ct(i / max(len(all_ct_names) - 1, 1))
                for i, ct in enumerate(all_ct_names)}
    col_colors = pd.DataFrame({
        "source": [ct_color[s] for s in sources],
        "target": [ct_color[t] for t in targets],
    }, index=plot_mat.columns)

    max_abs = max(float(np.abs(plot_mat.values).max()), 1e-6)

    try:
        g = sns.clustermap(
            plot_mat,
            cmap="RdBu_r",
            center=0,
            vmin=-max_abs, vmax=max_abs,
            row_cluster=True, col_cluster=True,
            metric="correlation",  # rows; sns uses same for cols by default
            method="average",
            row_colors=row_colors,
            col_colors=col_colors,
            figsize=(max(7, plot_mat.shape[1] * 0.45 + 4),
                     max(8, plot_mat.shape[0] * 0.22 + 4)),
            xticklabels=True, yticklabels=True,
            cbar_kws={"label": cbar_label},
            dendrogram_ratio=(0.12, 0.08),
        )
    except Exception as e:
        print(f"  [warn] clustermap failed ({e}); plotting unclustered")
        fig, ax = plt.subplots(figsize=(max(6, plot_mat.shape[1] * 0.45 + 2),
                                        max(8, plot_mat.shape[0] * 0.22 + 2)))
        im = ax.imshow(plot_mat.values, aspect="auto", cmap="RdBu_r",
                       vmin=-max_abs, vmax=max_abs)
        plt.colorbar(im, ax=ax, label=cbar_label)
        ax.set_xticks(range(plot_mat.shape[1]))
        ax.set_xticklabels(plot_mat.columns, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(plot_mat.shape[0]))
        ax.set_yticklabels(plot_mat.index, fontsize=6)
        ax.set_title(f"Δ LR activity: {grp_a} − {grp_b} | {age}{title_suffix}")
        out = pdir / f"delta_lr_heatmap_{_slug(grp_a)}_vs_{_slug(grp_b)}_{_slug(age)}{out_suffix}.png"
        fig.tight_layout()
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return

    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(),
                                  rotation=45, ha="right", fontsize=7)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=6)
    g.ax_heatmap.set_xlabel("Source→Target cell-type pair")
    g.ax_heatmap.set_ylabel("LR pair")
    cluster_note = " (clustered: corr/avg)"
    z_note = " [row-Z]" if zscore_rows else ""
    g.fig.suptitle(
        f"Δ LR activity: {grp_a} − {grp_b}  |  {age}{title_suffix}{z_note}{cluster_note}\n"
        f"(top {len(plot_mat)} LR pairs by |Δ|; "
        f"red={grp_a} stronger, blue={grp_b} stronger)",
        fontsize=10, y=1.02)

    z_suffix = "_zscore" if zscore_rows else ""
    out = (pdir / f"delta_lr_heatmap_{_slug(grp_a)}_vs_{_slug(grp_b)}_{_slug(age)}"
                  f"{out_suffix}{z_suffix}.png")
    g.fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(g.fig)
    print(f"  Plot: {out.name}")


# ============================================================================
# Pathway lookup helper (shared with overview pathway plot)
# ============================================================================

def get_pathway_map():
    """Return DataFrame[lr_pair, pathway] from liana mouseconsensus, or empty."""
    try:
        import liana as li
        res = li.rs.select_resource("mouseconsensus")
        pathway_cols = [c for c in res.columns
                        if any(k in c.lower() for k in ("pathway", "category", "family"))]
        if not pathway_cols:
            return pd.DataFrame()
        pc = pathway_cols[0]
        l_col = next((c for c in res.columns if "ligand" in c.lower()), None)
        r_col = next((c for c in res.columns if "receptor" in c.lower()), None)
        if not l_col or not r_col:
            return pd.DataFrame()
        out = res[[l_col, r_col, pc]].dropna(subset=[pc]).drop_duplicates()
        out.columns = ["ligand", "receptor", "pathway"]
        out["lr_pair"] = out["ligand"] + "→" + out["receptor"]
        return out[["lr_pair", "pathway"]]
    except Exception:
        return pd.DataFrame()


# ============================================================================
# Rank-rank scatter — central question for this dataset:
# "Do early and late stress hit the same LR programs?"
# ============================================================================

def plot_rank_rank_scatter(baseline_df, contrast_a, contrast_b, age,
                            magnitude_cutoff, pdir, top_n_label=20,
                            focus_celltypes=None):
    """Rank-rank scatter comparing two group-pair Δ-activity signatures.

    contrast_a, contrast_b: each a tuple (grp_test, grp_ctrl).
      e.g. (("Early_Stress","Relaxed"), ("Late_Stress","Relaxed")) → ES vs LS overlap
      e.g. (("Early_Stress","Relaxed"), ("Early_Stress","Late_Stress")) — etc.

    Plots Δ activity (test − ctrl) per LR×ct_pair for both contrasts. Each point
    is one LR×ct_pair. Points coloured by directional concordance:
      red    = both Δ > 0 (gained in both stress types)
      blue   = both Δ < 0 (lost in both)
      orange = discordant
      gray   = small (|Δ| below cutoff in both)
    Spearman rho displayed on plot.

    Reproducible from 08e_lr_baseline.csv.
    """
    from scipy.stats import spearmanr

    if baseline_df.empty or "magnitude_rank" not in baseline_df.columns:
        return

    def _activity_mean(grp):
        sub = baseline_df[
            (baseline_df["group"] == grp) & (baseline_df["age"] == age)
        ].copy()
        if sub.empty:
            return None
        sub["lr_pair"] = sub.apply(_lr_label, axis=1)
        sub["ct_pair"] = sub["source"] + "→" + sub["target"]
        sub["activity"] = 1 - sub["magnitude_rank"]
        if focus_celltypes:
            focus_set = set(focus_celltypes)
            sub = sub[sub["source"].isin(focus_set) | sub["target"].isin(focus_set)]
            if sub.empty:
                return None
        return (sub.groupby(["lr_pair", "ct_pair"])["activity"]
                .mean().reset_index())

    def _delta(test_grp, ctrl_grp):
        t = _activity_mean(test_grp)
        c = _activity_mean(ctrl_grp)
        if t is None or c is None:
            return None
        m = t.merge(c, on=["lr_pair", "ct_pair"], how="outer",
                    suffixes=("_t", "_c")).fillna(0)
        m["delta"] = m["activity_t"] - m["activity_c"]
        return m[["lr_pair", "ct_pair", "delta"]]

    d_a = _delta(*contrast_a)
    d_b = _delta(*contrast_b)
    if d_a is None or d_b is None:
        return

    merged = d_a.merge(d_b, on=["lr_pair", "ct_pair"], how="outer",
                       suffixes=("_a", "_b")).fillna(0)
    if merged.empty:
        return

    # Concordance color
    cutoff = 0.02  # tiny threshold to call something a "change"
    def _color(da, db):
        if abs(da) < cutoff and abs(db) < cutoff:
            return "lightgray"
        if da > 0 and db > 0:
            return "#d73027"  # both up
        if da < 0 and db < 0:
            return "#4575b4"  # both down
        return "#fdae61"      # discordant

    merged["color"] = [
        _color(a, b) for a, b in zip(merged["delta_a"], merged["delta_b"])
    ]

    label_a = f"{contrast_a[0]} − {contrast_a[1]}"
    label_b = f"{contrast_b[0]} − {contrast_b[1]}"

    rho, pval = spearmanr(merged["delta_a"], merged["delta_b"])

    fig, ax = plt.subplots(figsize=(7.5, 7))
    ax.axhline(0, color="k", lw=0.6, alpha=0.4)
    ax.axvline(0, color="k", lw=0.6, alpha=0.4)
    # Diagonal reference (concordant signature line)
    lim = float(max(merged[["delta_a", "delta_b"]].abs().max().max(), 0.05))
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.6, alpha=0.3)
    ax.scatter(merged["delta_a"], merged["delta_b"],
               c=merged["color"], s=22, alpha=0.75, linewidths=0)

    # Label top points by combined |Δ|
    merged["abs_combined"] = merged["delta_a"].abs() + merged["delta_b"].abs()
    top = merged.nlargest(top_n_label, "abs_combined")
    for _, r in top.iterrows():
        ax.annotate(f"{r['lr_pair']}\n[{r['ct_pair']}]",
                    (r["delta_a"], r["delta_b"]),
                    fontsize=5.5, ha="center", va="bottom",
                    xytext=(0, 4), textcoords="offset points",
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.4))

    ax.set_xlabel(f"Δ activity: {label_a}")
    ax.set_ylabel(f"Δ activity: {label_b}")
    sig_marker = "***" if pval < 1e-3 else "**" if pval < 1e-2 else "*" if pval < 0.05 else ""
    ax.set_title(
        f"LR signature concordance: {label_a}  vs  {label_b}\n"
        f"{age}   Spearman ρ = {rho:.2f}{sig_marker}   "
        f"(p={pval:.1e}; n={len(merged):,} LR×ct_pair)",
        fontsize=9)

    # Legend
    from matplotlib.patches import Patch
    handles = [
        Patch(color="#d73027", label="Concordant ↑ (both gained)"),
        Patch(color="#4575b4", label="Concordant ↓ (both lost)"),
        Patch(color="#fdae61", label="Discordant"),
        Patch(color="lightgray", label="Negligible Δ"),
    ]
    ax.legend(handles=handles, fontsize=7, loc="best")

    focus_suffix = "_focused" if focus_celltypes else ""
    out = (pdir / f"rank_rank_{_slug(contrast_a[0])}-{_slug(contrast_a[1])}"
                  f"_vs_{_slug(contrast_b[0])}-{_slug(contrast_b[1])}"
                  f"_{_slug(age)}{focus_suffix}.png")
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}  (ρ={rho:.2f})")


# ============================================================================
# Sender / receiver plots
# ============================================================================

def _compute_sender_receiver(baseline_df, magnitude_cutoff):
    """Inline copy — avoids circular import with main module."""
    records = []
    for (group, age), grp in baseline_df.groupby(["group", "age"]):
        active = grp[grp["magnitude_rank"] < magnitude_cutoff]
        for ct in pd.concat([active["source"], active["target"]]).unique():
            src = active[active["source"] == ct]
            tgt = active[active["target"] == ct]
            records.append({
                "group": group, "age": age, "cell_type": ct,
                "n_sent": len(src), "n_received": len(tgt),
                "sender_score": float((1 - src["magnitude_rank"]).mean()) if len(src) else 0.0,
                "receiver_score": float((1 - tgt["magnitude_rank"]).mean()) if len(tgt) else 0.0,
            })
    return pd.DataFrame(records)


def plot_sender_receiver_bubble(baseline_df, age, magnitude_cutoff, pdir):
    """Bubble plot: x=sender score, y=receiver score, one bubble per cell type per group.
    Shows whether a cell type is primarily sender, receiver, or hub,
    and how stress shifts that role."""
    active = baseline_df[
        (baseline_df["age"] == age) &
        (baseline_df["magnitude_rank"] < magnitude_cutoff)
    ]
    if active.empty:
        return

    sr = _compute_sender_receiver(active, magnitude_cutoff)
    if sr.empty:
        return

    groups = sorted(sr["group"].unique())
    fig, axes = plt.subplots(1, len(groups), figsize=(5 * len(groups), 5),
                             sharey=True, sharex=True, constrained_layout=True)
    if len(groups) == 1:
        axes = [axes]

    all_cts = sorted(sr["cell_type"].unique())
    ct_colors = {ct: plt.cm.tab20(i / max(len(all_cts) - 1, 1))
                 for i, ct in enumerate(all_cts)}

    for ax, grp in zip(axes, groups):
        grp_sr = sr[sr["group"] == grp]
        for _, row in grp_sr.iterrows():
            ct = row["cell_type"]
            size = (row["n_sent"] + row["n_received"]) * 20 + 50
            ax.scatter(row["sender_score"], row["receiver_score"],
                       s=size, color=ct_colors.get(ct, "gray"),
                       alpha=0.8, edgecolors="k", linewidths=0.5, zorder=3)
            ax.annotate(ct, (row["sender_score"], row["receiver_score"]),
                        fontsize=7, ha="center", va="bottom",
                        xytext=(0, 5), textcoords="offset points")
        lim = max(grp_sr[["sender_score", "receiver_score"]].max().max(), 0.01)
        ax.plot([0, lim], [0, lim], "k--", lw=0.8, alpha=0.4)
        ax.set_xlabel("Sender score")
        ax.set_ylabel("Receiver score")
        ax.set_title(grp, fontsize=9)
        ax.set_xlim(left=0); ax.set_ylim(bottom=0)

    fig.suptitle(f"Sender vs receiver roles: {age}\n(size = total active LR pairs)",
                 fontsize=10)
    out = pdir / f"sender_receiver_bubble_{_slug(age)}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Plot: {out.name}")


def plot_sender_receiver_heatmap(sr_df, age, ref_group, pdir):
    """Δ sender/receiver score (all groups − ref) per cell type, for one age."""
    sub = sr_df[sr_df["age"] == age].copy()
    if sub.empty:
        return

    ref = sub[sub["group"] == ref_group].set_index("cell_type")[
        ["sender_score", "receiver_score"]]
    if ref.empty:
        return

    groups = [g for g in sub["group"].unique() if g != ref_group]
    if not groups:
        return

    cts = sorted(ref.index)
    fig, axes = plt.subplots(1, len(groups),
                             figsize=(6 * len(groups), max(4, len(cts) * 0.4 + 2)),
                             sharey=True, constrained_layout=True)
    if len(groups) == 1:
        axes = [axes]

    for ax, grp in zip(axes, sorted(groups)):
        test = sub[sub["group"] == grp].set_index("cell_type")[
            ["sender_score", "receiver_score"]]
        delta = test.reindex(cts).fillna(0) - ref.reindex(cts).fillna(0)
        data = delta[["sender_score", "receiver_score"]].values
        max_abs = max(np.abs(data).max(), 1e-6)
        im = ax.imshow(data, aspect="auto", cmap="RdBu_r",
                       vmin=-max_abs, vmax=max_abs)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                     label=f"Δ score ({grp} − {ref_group})")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Δ sender", "Δ receiver"], fontsize=9)
        ax.set_yticks(range(len(cts)))
        ax.set_yticklabels(cts, fontsize=7)
        ax.set_title(f"{grp} − {ref_group}", fontsize=9)

    fig.suptitle(f"Δ sender/receiver vs {ref_group}: {age}", fontsize=11)
    out = pdir / f"sender_receiver_heatmap_{_slug(age)}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Plot: {out.name}")


def plot_delta_sender_receiver_heatmap(sr_df, grp_a, grp_b, age, pdir):
    """Δ sender/receiver score: grp_a − grp_b (any pair, not just vs ref).

    Covers ES-v-LS comparison that the ref-only heatmap misses.
    """
    sub = sr_df[sr_df["age"] == age].copy()
    if sub.empty:
        return

    df_a = sub[sub["group"] == grp_a].set_index("cell_type")[
        ["sender_score", "receiver_score"]]
    df_b = sub[sub["group"] == grp_b].set_index("cell_type")[
        ["sender_score", "receiver_score"]]
    if df_a.empty or df_b.empty:
        return

    cts = sorted(set(df_a.index) | set(df_b.index))
    delta = (df_a.reindex(cts).fillna(0) - df_b.reindex(cts).fillna(0))
    data = delta[["sender_score", "receiver_score"]].values
    max_abs = max(np.abs(data).max(), 1e-6)

    fig, ax = plt.subplots(figsize=(4, max(4, len(cts) * 0.4 + 2)))
    im = ax.imshow(data, aspect="auto", cmap="RdBu_r",
                   vmin=-max_abs, vmax=max_abs)
    plt.colorbar(im, ax=ax, label=f"Δ ({grp_a} − {grp_b})")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Δ sender", "Δ receiver"], fontsize=9)
    ax.set_yticks(range(len(cts)))
    ax.set_yticklabels(cts, fontsize=7)
    ax.set_title(f"Δ sender/receiver\n{grp_a} − {grp_b}  |  {age}", fontsize=9)
    fig.tight_layout()
    out = pdir / f"delta_sr_{_slug(grp_a)}_vs_{_slug(grp_b)}_{_slug(age)}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}")
