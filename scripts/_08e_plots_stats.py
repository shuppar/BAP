"""
_08e_plots_stats.py — statistically-grounded Phase 8e plots.

Imported by 08e_comms_summary.py. Not a standalone entry point.

These REPLACE the weak/degenerate plots (the old volcano put |stat| on y → a
mathematical 'V' with no significance axis) and ADD the Δ signalling network /
Δ cell-type heatmap. All stress figures here read REAL significance from the CSVs:
  - differential arm: `interaction_padj` (LIANA's combined interaction FDR)
  - per-donor arm:     handled in _08e_plots_perdonor.py (MW-U `fdr`)
The baseline arm stays DESCRIPTIVE (pooled cells, no stress test) — its plots
live in the unchanged _08e_plots_baseline.py.

Functions:
  plot_differential_volcano        — real volcano: x=interaction_stat,
                                      y=-log10(interaction_padj); FDR<0.05 coloured;
                                      leader-line labels (no overlap).
  plot_differential_dotplot        — top sig LR pairs; colour=direction,
                                      size=-log10(FDR), ★ if FDR<0.05.
  plot_delta_network               — cell-type→cell-type signalling graph;
                                      edges aggregate signed interaction_stat over
                                      FDR<thresh pairs (red=up/blue=down in test).
  plot_delta_celltype_heatmap      — source×target matrix of net signed stat over
                                      sig pairs (chord-equivalent, very readable).
  plot_sender_receiver_bubble      — fixed: area-scaled bubbles + leader-line labels.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FDR_THRESH = 0.05
_NEG_LOG_CAP = 50.0   # cap -log10(padj) so padj==0 doesn't blow up the axis


def _slug(s):
    return str(s).replace(" ", "_").replace("/", "-").replace(".", "")


def _lr_label(row):
    return f"{row['ligand_complex']}→{row['receptor_complex']}"


def _apply_stat_floor(sub, quantile):
    """Slice-specific effect floor: drop rows whose |interaction_stat| is below the
    `quantile`-th percentile WITHIN this already-significance-filtered slice. Returns
    the filtered frame (unchanged if quantile<=0 or slice too small)."""
    if quantile is None or quantile <= 0 or sub.empty:
        return sub
    floor = sub["interaction_stat"].abs().quantile(quantile)
    out = sub[sub["interaction_stat"].abs() >= floor]
    return out if not out.empty else sub


def _neg_log10(padj):
    p = np.asarray(padj, dtype=float)
    p = np.where(~np.isfinite(p) | (p <= 0), np.nan, p)
    out = -np.log10(p)
    floor = np.nanmin(p[p > 0]) if np.any(p > 0) else 1e-300
    out = np.where(np.isnan(out), -np.log10(floor), out)
    return np.clip(out, 0, _NEG_LOG_CAP)


def _stack_labels(ax, items, side, fontsize=6.5):
    """Place labels along one vertical margin with leader lines (no overlap).

    items: list of (x_data, y_data, text). side: 'left' or 'right'.
    Labels are stacked evenly in the upper portion of the axis on the chosen
    side; a thin line connects each label to its data point.
    """
    if not items:
        return
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    items = sorted(items, key=lambda t: t[1], reverse=True)  # by y desc
    n = len(items)
    # label x anchor just inside the spine; stack y evenly over the top 92%..40%
    if side == "right":
        lx = xlim[0] + 0.62 * (xlim[1] - xlim[0])
        ha = "left"
    else:
        lx = xlim[0] + 0.38 * (xlim[1] - xlim[0])
        ha = "right"
    y_hi = ylim[0] + 0.97 * (ylim[1] - ylim[0])
    y_lo = ylim[0] + 0.42 * (ylim[1] - ylim[0])
    ys = np.linspace(y_hi, y_lo, n) if n > 1 else [y_hi]
    for (xd, yd, txt), ly in zip(items, ys):
        ax.annotate(
            txt, xy=(xd, yd), xytext=(lx, ly), fontsize=fontsize, ha=ha, va="center",
            arrowprops=dict(arrowstyle="-", lw=0.4, color="0.6",
                            connectionstyle="arc3,rad=0.0"),
        )


# ---------------------------------------------------------------------------
# Differential volcano (REAL — uses interaction_padj)
# ---------------------------------------------------------------------------

def plot_differential_volcano(diff_df, contrast_name, age, top_n, pdir):
    sub = diff_df[(diff_df["contrast_name"] == contrast_name) &
                  (diff_df["age"].astype(str) == str(age))].copy()
    sub = sub.dropna(subset=["interaction_stat", "interaction_padj"])
    if sub.empty:
        return
    sub["nlp"] = _neg_log10(sub["interaction_padj"])
    sub["sig"] = sub["interaction_padj"] < FDR_THRESH
    sub["lr"] = sub.apply(_lr_label, axis=1)

    up = sub[sub["sig"] & (sub["interaction_stat"] > 0)]
    dn = sub[sub["sig"] & (sub["interaction_stat"] < 0)]
    ns = sub[~sub["sig"]]

    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.scatter(ns["interaction_stat"], ns["nlp"], s=6, c="0.8", alpha=0.5,
               linewidths=0, label="ns")
    ax.scatter(dn["interaction_stat"], dn["nlp"], s=10, c="#4575b4", alpha=0.7,
               linewidths=0, label=f"down (FDR<{FDR_THRESH})")
    ax.scatter(up["interaction_stat"], up["nlp"], s=10, c="#d73027", alpha=0.7,
               linewidths=0, label=f"up (FDR<{FDR_THRESH})")
    ax.axhline(-np.log10(FDR_THRESH), color="k", lw=0.7, ls="--")
    ax.axvline(0, color="k", lw=0.6)

    # leader-line labels: top_n significant pairs by |stat|, split by direction
    sig = sub[sub["sig"]].copy()
    sig["abs"] = sig["interaction_stat"].abs()
    top = sig.nlargest(top_n, "abs")
    left = [(r["interaction_stat"], r["nlp"], f"{r['lr']}\n[{r['source']}→{r['target']}]")
            for _, r in top[top["interaction_stat"] < 0].iterrows()]
    right = [(r["interaction_stat"], r["nlp"], f"{r['lr']}\n[{r['source']}→{r['target']}]")
             for _, r in top[top["interaction_stat"] >= 0].iterrows()]
    _stack_labels(ax, left, "left")
    _stack_labels(ax, right, "right")

    ax.set_xlabel("interaction_stat  (positive = up in test group)")
    ax.set_ylabel("−log10(interaction FDR)")
    ax.set_title(f"Differential LR volcano — {contrast_name} | {age}\n"
                 f"{len(up)} up / {len(dn)} down at FDR<{FDR_THRESH} "
                 f"(top {top_n} labelled)")
    ax.legend(fontsize=7, loc="lower center", ncol=3, frameon=False)
    fig.tight_layout()
    out = pdir / f"volcano_{_slug(contrast_name)}_{_slug(age)}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}  ({len(up)}↑/{len(dn)}↓ sig)")


# ---------------------------------------------------------------------------
# Differential dotplot (top sig pairs, ★ FDR<0.05)
# ---------------------------------------------------------------------------

def plot_differential_dotplot(diff_df, contrast_name, age, top_n, pdir, quantile=0.0):
    sub = diff_df[(diff_df["contrast_name"] == contrast_name) &
                  (diff_df["age"].astype(str) == str(age))].copy()
    sub = sub.dropna(subset=["interaction_stat", "interaction_padj"])
    if sub.empty:
        return
    sub["nlp"] = _neg_log10(sub["interaction_padj"])
    sub["abs"] = sub["interaction_stat"].abs()
    sub["sig"] = sub["interaction_padj"] < FDR_THRESH
    sig = _apply_stat_floor(sub[sub["sig"]], quantile)
    top = sig.nlargest(top_n, "abs").sort_values("interaction_stat")
    if top.empty:
        return
    labels = [f"{_lr_label(r)}  [{r['source']}→{r['target']}]" for _, r in top.iterrows()]
    colors = ["#d73027" if s > 0 else "#4575b4" for s in top["interaction_stat"]]
    sizes = 20 + 18 * top["nlp"].clip(0, 12)

    fig, ax = plt.subplots(figsize=(8.5, max(4, len(top) * 0.34 + 1)))
    ax.scatter(top["interaction_stat"], range(len(top)), s=sizes, c=colors,
               alpha=0.85, edgecolors="k", linewidths=0.3)
    for i, (_, r) in enumerate(top.iterrows()):
        if r["sig"]:
            x = r["interaction_stat"]
            ax.text(x + (0.05 if x >= 0 else -0.05), i, "★", va="center",
                    ha="left" if x >= 0 else "right", fontsize=8, color="k")
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("interaction_stat  (positive = up in test group)")
    ax.set_title(f"Top {top_n} differential LR pairs — {contrast_name} | {age}\n"
                 f"size = −log10(FDR); ★ = FDR<{FDR_THRESH}; red=up / blue=down")
    fig.tight_layout()
    out = pdir / f"dotplot_{_slug(contrast_name)}_{_slug(age)}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}")


# ---------------------------------------------------------------------------
# Δ signalling network (cell-type graph from differential, FDR-filtered)
# ---------------------------------------------------------------------------

def _aggregate_celltype_edges(diff_df, contrast_name, age, fdr_thresh, focus=None,
                              quantile=0.0):
    sub = diff_df[(diff_df["contrast_name"] == contrast_name) &
                  (diff_df["age"].astype(str) == str(age))].copy()
    sub = sub.dropna(subset=["interaction_stat", "interaction_padj"])
    sub = sub[sub["interaction_padj"] < fdr_thresh]
    sub = _apply_stat_floor(sub, quantile)
    if focus:
        sub = sub[sub["source"].isin(focus) | sub["target"].isin(focus)]
    if sub.empty:
        return pd.DataFrame()
    agg = (sub.groupby(["source", "target"])
           .agg(net_stat=("interaction_stat", "mean"),
                n_sig=("interaction_stat", "size"),
                n_up=("interaction_stat", lambda s: int((s > 0).sum())),
                n_down=("interaction_stat", lambda s: int((s < 0).sum())))
           .reset_index())
    return agg


def plot_delta_network(diff_df, contrast_name, age, pdir,
                       fdr_thresh=FDR_THRESH, focus=None, min_edges=1, quantile=0.0):
    agg = _aggregate_celltype_edges(diff_df, contrast_name, age, fdr_thresh, focus, quantile)
    if agg.empty or len(agg) < min_edges:
        print(f"  [info] delta_network {contrast_name}/{age}: no FDR<{fdr_thresh} edges")
        return
    nodes = sorted(set(agg["source"]) | set(agg["target"]))
    ang = {n: 2 * np.pi * i / len(nodes) for i, n in enumerate(nodes)}
    pos = {n: (np.cos(a), np.sin(a)) for n, a in ang.items()}

    wmax = agg["n_sig"].max() or 1            # breadth → alpha
    smax = agg["net_stat"].abs().max() or 1.0  # magnitude → width

    fig, ax = plt.subplots(figsize=(9, 9))
    for _, e in agg.iterrows():
        x0, y0 = pos[e["source"]]
        x1, y1 = pos[e["target"]]
        col = "#d73027" if e["net_stat"] > 0 else "#4575b4"
        lw = 0.6 + 5.0 * abs(e["net_stat"]) / smax     # WIDTH = effect magnitude
        alpha = 0.25 + 0.6 * e["n_sig"] / wmax          # ALPHA = breadth (#sig pairs)
        # curved edge to disambiguate direction (source→target)
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=lw, alpha=alpha,
                                    connectionstyle="arc3,rad=0.12", shrinkA=14, shrinkB=14))
    for n, (x, y) in pos.items():
        ax.scatter([x], [y], s=420, c="white", edgecolors="k", linewidths=1.2, zorder=3)
        ax.text(x, y, n, ha="center", va="center", fontsize=7, zorder=4, wrap=True)
    ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.35, 1.35)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"Δ signalling network — {contrast_name} | {age}\n"
                 f"edges = cell-type pairs with FDR<{fdr_thresh} LR pairs "
                 f"(width = |mean interaction_stat| effect magnitude; "
                 f"alpha = #sig pairs breadth;\n"
                 f"red = up / blue = down in test; arrow = source→target)", fontsize=10)
    fig.tight_layout()
    out = pdir / f"delta_network_{_slug(contrast_name)}_{_slug(age)}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}  ({len(agg)} sig cell-type edges)")


def plot_delta_celltype_heatmap(diff_df, contrast_name, age, pdir, fdr_thresh=FDR_THRESH, quantile=0.0):
    """source×target matrix of net signed interaction_stat over FDR<thresh pairs.
    Chord-equivalent, fully readable, carries direction + significance."""
    agg = _aggregate_celltype_edges(diff_df, contrast_name, age, fdr_thresh, quantile=quantile)
    if agg.empty:
        return
    cts = sorted(set(agg["source"]) | set(agg["target"]))
    mat = pd.DataFrame(0.0, index=cts, columns=cts)
    cnt = pd.DataFrame(0, index=cts, columns=cts)
    for _, e in agg.iterrows():
        mat.loc[e["source"], e["target"]] = e["net_stat"]
        cnt.loc[e["source"], e["target"]] = e["n_sig"]
    vmax = np.abs(mat.values).max() or 1.0

    fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(cts) + 2), max(5, 0.55 * len(cts) + 1)))
    im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(cts))); ax.set_xticklabels(cts, rotation=90, fontsize=7)
    ax.set_yticks(range(len(cts))); ax.set_yticklabels(cts, fontsize=7)
    ax.set_xlabel("target (receiver)"); ax.set_ylabel("source (sender)")
    for i in range(len(cts)):
        for j in range(len(cts)):
            c = cnt.values[i, j]
            if c > 0:
                ax.text(j, i, int(c), ha="center", va="center", fontsize=6, color="k")
    fig.colorbar(im, ax=ax, fraction=0.046, label="mean signed interaction_stat")
    ax.set_title(f"Δ cell-type signalling — {contrast_name} | {age}\n"
                 f"colour = mean signed stat (red=up/blue=down in test); "
                 f"number = #FDR<{fdr_thresh} LR pairs", fontsize=9)
    fig.tight_layout()
    out = pdir / f"delta_celltype_heatmap_{_slug(contrast_name)}_{_slug(age)}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}")


# ---------------------------------------------------------------------------
# Δ chord (directional ribbons; width = |Δ|; sig edges outlined)
# ---------------------------------------------------------------------------

def _bezier(p0, p1, p2, n=40):
    t = np.linspace(0, 1, n)[:, None]
    return (1 - t) ** 2 * p0 + 2 * (1 - t) * t * p1 + t ** 2 * p2


def _plot_delta_chord_unused(diff_df, contrast_name, age, pdir,
                     fdr_thresh=FDR_THRESH, focus=None, max_nodes_warn=12,
                     all_edges=True):
    """DEAD CODE — superseded by the plot_delta_chord defined later in this module.
    Kept (renamed) only because the matplotlib-Path helpers below it are shared by
    the active version; do not call."""
    sub = diff_df[(diff_df["contrast_name"] == contrast_name) &
                  (diff_df["age"].astype(str) == str(age))].copy()
    sub = sub.dropna(subset=["interaction_stat", "interaction_padj"])
    if focus:
        sub = sub[sub["source"].isin(focus) | sub["target"].isin(focus)]
    if sub.empty:
        return
    sub["sig_pair"] = sub["interaction_padj"] < fdr_thresh
    agg = (sub.groupby(["source", "target"])
           .agg(net_stat=("interaction_stat", "mean"),
                n_sig=("sig_pair", "sum"),
                n_tot=("sig_pair", "size")).reset_index())
    agg = agg[agg["net_stat"].abs() > 0]
    if not all_edges:
        agg = agg[agg["n_sig"] > 0]
    if agg.empty:
        print(f"  [info] delta_chord {contrast_name}/{age}: nothing to draw")
        return

    nodes = sorted(set(agg["source"]) | set(agg["target"]))
    if len(nodes) > max_nodes_warn:
        print(f"  [warn] delta_chord {contrast_name}/{age}: {len(nodes)} cell types "
              f"— chord will be dense (network/heatmap are clearer at this size)")

    N = len(nodes)
    gap = 0.012 * 2 * np.pi
    seg = (2 * np.pi - N * gap) / N
    arc = {}                      # node -> (start_angle, end_angle)
    a = 0.0
    for n in nodes:
        arc[n] = (a, a + seg)
        a += seg + gap
    R = 1.0
    smax = agg["net_stat"].abs().max() or 1.0

    fig, ax = plt.subplots(figsize=(10, 10))
    # node arcs + labels
    for n, (a0, a1) in arc.items():
        aa = np.linspace(a0, a1, 30)
        ax.plot(R * np.cos(aa), R * np.sin(aa), lw=6, color="0.7", solid_capstyle="butt")
        am = (a0 + a1) / 2
        rot = np.degrees(am)
        if 90 < rot < 270:
            rot += 180
        ax.text(1.12 * np.cos(am), 1.12 * np.sin(am), n, rotation=rot,
                ha="center", va="center", fontsize=7)

    # sort so significant (outlined) ribbons draw on top
    agg = agg.sort_values("n_sig")
    # allocate sub-arc slots per node so ribbons don't all stack at the arc centre
    src_count = agg.groupby("source").cumcount()
    src_tot = agg.groupby("source")["source"].transform("size")
    tgt_count = agg.groupby("target").cumcount()
    tgt_tot = agg.groupby("target")["target"].transform("size")
    agg = agg.assign(_si=src_count, _sn=src_tot, _ti=tgt_count, _tn=tgt_tot)

    for _, e in agg.iterrows():
        col = "#d73027" if e["net_stat"] > 0 else "#4575b4"
        sig = e["n_sig"] > 0
        w = 1.0 + 9.0 * abs(e["net_stat"]) / smax       # ribbon width ∝ |Δ|
        alpha = 0.85 if sig else 0.28
        a0s, a1s = arc[e["source"]]
        a0t, a1t = arc[e["target"]]
        # place this ribbon's source/target anchor within the node's arc
        sa = a0s + (e["_si"] + 0.5) / e["_sn"] * (a1s - a0s)
        ta = a0t + (e["_ti"] + 0.5) / e["_tn"] * (a1t - a0t)
        p_src = np.array([R * np.cos(sa), R * np.sin(sa)])
        p_tgt = np.array([R * np.cos(ta), R * np.sin(ta)])
        ctrl = np.array([0.0, 0.0])                      # curve through centre
        curve = _bezier(p_src, ctrl, p_tgt)
        # taper: thin at source, thick at target (encodes direction)
        lws = np.linspace(0.4 * w, w, len(curve))
        for i in range(len(curve) - 1):
            ax.plot(curve[i:i+2, 0], curve[i:i+2, 1], color=col, alpha=alpha,
                    lw=lws[i], solid_capstyle="round", zorder=2 if sig else 1)
        if sig:                                          # dark outline for sig edges
            ax.plot(curve[:, 0], curve[:, 1], color="k", alpha=0.5,
                    lw=w + 1.1, zorder=1.5)
        # arrowhead at target end (direction)
        d = curve[-1] - curve[-3]
        ax.annotate("", xy=p_tgt, xytext=p_tgt - 0.06 * d / (np.linalg.norm(d) + 1e-9),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.0,
                                    alpha=alpha), zorder=3)

    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3)
    ax.set_aspect("equal"); ax.axis("off")
    n_sig_edges = int((agg["n_sig"] > 0).sum())
    ax.set_title(f"Δ signalling chord — {contrast_name} | {age}\n"
                 f"ribbon width = |mean Δ interaction_stat|; thin→thick = source→target; "
                 f"red=up/blue=down in test;\n"
                 f"outlined = ≥1 FDR<{fdr_thresh} LR pair ({n_sig_edges} sig edges; "
                 f"faded = non-sig context)", fontsize=9)
    fig.tight_layout()
    out = pdir / f"delta_chord_{_slug(contrast_name)}_{_slug(age)}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}  ({len(agg)} edges, {n_sig_edges} sig)")


# ---------------------------------------------------------------------------
# Δ chord (directed tapered ribbons; sign-on-outline; FDR<thresh only)
# ---------------------------------------------------------------------------

def _ribbon_path(S, T, w_s, w_t, bend=0.45):
    """Filled directed ribbon: narrow (half-width w_s) at source S, wide (w_t) at
    target T. Two cubic Béziers curved toward centre. Returns a matplotlib Path."""
    from matplotlib.path import Path as MplPath
    S = np.asarray(S, float); T = np.asarray(T, float)
    ts = np.array([-S[1], S[0]]); ts /= (np.linalg.norm(ts) or 1)   # tangent @ S
    tt = np.array([-T[1], T[0]]); tt /= (np.linalg.norm(tt) or 1)   # tangent @ T
    Sl, Sr = S + w_s * ts, S - w_s * ts
    Tl, Tr = T + w_t * tt, T - w_t * tt
    cS, cT = S * bend, T * bend     # control points pulled toward origin
    verts = [tuple(Sl),
             tuple(cS), tuple(cT), tuple(Tr),     # Sl -> Tr (curve)
             tuple(Tl),                            # across target base
             tuple(cT), tuple(cS), tuple(Sr),     # Tl -> Sr (curve back)
             tuple(Sl)]                            # close
    codes = [MplPath.MOVETO,
             MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.LINETO,
             MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.CLOSEPOLY]
    return MplPath(verts, codes)


def plot_delta_chord(diff_df, contrast_name, age, pdir, fdr_thresh=FDR_THRESH, focus=None, quantile=0.0):
    """Directed chord of significant signalling changes.
    Encoding: taper = direction (narrow=source → wide=target);
              width  = #FDR<thresh LR pairs (strength of change);
              fill   = source cell-type colour (trace the sender);
              outline= red (net up) / blue (net down) in test group;
              only source→target pairs with ≥1 FDR<thresh LR pair are drawn.
    No cell-type-count cap (drawn even for many types)."""
    from matplotlib.patches import PathPatch, Patch
    agg = _aggregate_celltype_edges(diff_df, contrast_name, age, fdr_thresh, focus, quantile)
    if agg.empty:
        print(f"  [info] delta_chord {contrast_name}/{age}: no FDR<{fdr_thresh} edges")
        return
    nodes = sorted(set(agg["source"]) | set(agg["target"]))
    n = len(nodes)
    ang = {nd: 2 * np.pi * i / n for i, nd in enumerate(nodes)}
    pos = {nd: np.array([np.cos(a), np.sin(a)]) for nd, a in ang.items()}
    colmap = {nd: plt.cm.tab20(i / max(n - 1, 1)) for i, nd in enumerate(nodes)}
    nmax = agg["n_sig"].max() or 1
    smax = agg["net_stat"].abs().max() or 1.0

    if n > 12:
        print(f"  [warn] delta_chord {contrast_name}/{age}: {n} cell types — dense "
              f"(drawn anyway; network/heatmap are more legible at this size).")

    fig, ax = plt.subplots(figsize=(10, 10))
    # draw widest ribbons first so thin ones stay visible on top
    for _, e in agg.assign(_mag=agg["net_stat"].abs()).sort_values(
            "_mag", ascending=False).iterrows():
        s, t = e["source"], e["target"]
        w_t = 0.02 + 0.13 * abs(e["net_stat"]) / smax      # WIDTH = effect magnitude
        w_s = 0.012
        alpha = 0.2 + 0.6 * e["n_sig"] / nmax              # ALPHA = breadth (#sig pairs)
        sign_col = "#d73027" if e["net_stat"] > 0 else "#4575b4"
        if s == t:                       # autocrine: small loop just outside node
            P = pos[s]; r = 0.05 + 0.08 * abs(e["net_stat"]) / smax
            loop = plt.Circle(P * 1.12, r, facecolor=colmap[s], edgecolor=sign_col,
                              lw=1.4, alpha=alpha, zorder=2)
            ax.add_patch(loop)
            continue
        path = _ribbon_path(pos[s], pos[t], w_s, w_t)
        ax.add_patch(PathPatch(path, facecolor=colmap[s], edgecolor=sign_col,
                               lw=1.3, alpha=alpha, zorder=2))
    # nodes + labels
    for nd, P in pos.items():
        ax.scatter([P[0]], [P[1]], s=300, c=[colmap[nd]], edgecolors="k",
                   linewidths=1.0, zorder=4)
        lab = P * 1.18
        ax.text(lab[0], lab[1], nd, ha="center", va="center", fontsize=7, zorder=5)
    ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.4, 1.4)
    ax.set_aspect("equal"); ax.axis("off")
    leg = [Patch(facecolor="0.7", edgecolor="#d73027", label="net up in test (outline)"),
           Patch(facecolor="0.7", edgecolor="#4575b4", label="net down in test (outline)")]
    ax.legend(handles=leg, fontsize=8, loc="lower right", frameon=False)
    ax.set_title(f"Δ signalling chord — {contrast_name} | {age}\n"
                 f"taper = source→target; width = |mean interaction_stat| (effect "
                 f"magnitude); alpha = #FDR<{fdr_thresh} LR pairs (breadth);\n"
                 f"fill = source; outline = up(red)/down(blue) in test", fontsize=9)
    fig.tight_layout()
    out = pdir / f"delta_chord_{_slug(contrast_name)}_{_slug(age)}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}  ({len(agg)} directed sig ribbons)")


# ---------------------------------------------------------------------------
# Sender / receiver up-down bars (differential polarization, FDR-backed)
# ---------------------------------------------------------------------------

def plot_sender_receiver_updown_bars(diff_df, contrast_name, age, pdir, fdr=0.05, quantile=0.0):
    """Per cell type: # significantly UP vs DOWN differential LR pairs, split by role
    (as sender / as receiver). The honest summary of differential signalling
    polarization (matches the 8b per-celltype up/down structure). FDR-backed."""
    d = diff_df[(diff_df["contrast_name"] == contrast_name) &
                (diff_df["age"].astype(str) == str(age))].copy()
    d = d.dropna(subset=["interaction_stat", "interaction_padj"])
    d = d[d["interaction_padj"] < fdr]
    d = _apply_stat_floor(d, quantile)
    if d.empty:
        return
    d["dir"] = np.where(d["interaction_stat"] > 0, "up", "down")
    snd = d.groupby(["source", "dir"]).size().unstack("dir").fillna(0)
    rcv = d.groupby(["target", "dir"]).size().unstack("dir").fillna(0)
    for df_ in (snd, rcv):
        for col in ("up", "down"):
            if col not in df_.columns:
                df_[col] = 0
    cts = sorted(set(snd.index) | set(rcv.index),
                 key=lambda c: -(snd.reindex([c]).fillna(0).values.sum()
                                 + rcv.reindex([c]).fillna(0).values.sum()))
    if not cts:
        return
    y = np.arange(len(cts))
    fig, axes = plt.subplots(1, 2, figsize=(11, max(4, len(cts) * 0.34 + 1)), sharey=True)
    for ax, mat, role in ((axes[0], snd, "as sender"), (axes[1], rcv, "as receiver")):
        up = mat.reindex(cts)["up"].fillna(0).values
        dn = -mat.reindex(cts)["down"].fillna(0).values
        ax.barh(y, up, color="#d73027", label="up in stress")
        ax.barh(y, dn, color="#4575b4", label="down in stress")
        ax.axvline(0, color="k", lw=0.6)
        ax.set_yticks(y); ax.set_yticklabels(cts, fontsize=7)
        ax.set_title(role, fontsize=9)
        ax.set_xlabel("# sig LR pairs (← down | up →)")
    axes[0].invert_yaxis()
    axes[0].legend(fontsize=7, loc="lower right", frameon=False)
    fig.suptitle(f"Differential signalling polarization — {contrast_name} | {age}\n"
                 f"sig LR pairs at FDR<{fdr}, by cell-type role", fontsize=10)
    fig.tight_layout()
    out = pdir / f"updown_bars_{_slug(contrast_name)}_{_slug(age)}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}")


# ---------------------------------------------------------------------------
# Sender / receiver bubble (FIXED: area-scaled, leader-line labels)
# ---------------------------------------------------------------------------

def plot_sender_receiver_bubble(baseline_df, age, magnitude_cutoff, pdir):
    sub = baseline_df[baseline_df["age"].astype(str) == str(age)].copy()
    sub = baseline_df[baseline_df["age"].astype(str) == str(age)].copy()
    if sub.empty or "magnitude_rank" not in sub.columns:
        return
    active = sub[sub["magnitude_rank"] < magnitude_cutoff]
    if active.empty:
        return
    groups = sorted(active["group"].astype(str).unique())
    fig, axes = plt.subplots(1, len(groups), figsize=(6.0 * len(groups), 5.6),
                             squeeze=False)
    axes = axes[0]
    for ax, grp in zip(axes, groups):
        g = active[active["group"].astype(str) == grp]
        cts = sorted(set(g["source"]) | set(g["target"]))
        rows = []
        for ct in cts:
            sent = g[g["source"] == ct]
            recd = g[g["target"] == ct]
            rows.append({
                "ct": ct,
                "sender": float((1 - sent["magnitude_rank"]).mean()) if len(sent) else 0.0,
                "receiver": float((1 - recd["magnitude_rank"]).mean()) if len(recd) else 0.0,
                "n": len(sent) + len(recd),
            })
        d = pd.DataFrame(rows)
        nmax = d["n"].max() or 1
        sizes = 30 + 600 * d["n"] / nmax          # area-scaled, capped
        colors = plt.cm.tab20(np.linspace(0, 1, len(d)))
        ax.scatter(d["sender"], d["receiver"], s=sizes, c=colors,
                   alpha=0.7, edgecolors="k", linewidths=0.5, zorder=2)
        lim = max(d["sender"].max(), d["receiver"].max()) * 1.15 + 0.05
        ax.plot([0, lim], [0, lim], ls="--", c="0.6", lw=0.8, zorder=1)
        # leader-line labels along the right margin, ordered by total activity
        items = sorted([(r["sender"], r["receiver"], r["ct"]) for _, r in d.iterrows()],
                       key=lambda t: t[0] + t[1], reverse=True)
        ys = np.linspace(lim * 0.97, lim * 0.05, len(items)) if len(items) > 1 else [lim * 0.5]
        for (sx, sy, txt), ly in zip(items, ys):
            ax.annotate(txt, xy=(sx, sy), xytext=(lim * 1.02, ly), fontsize=6.5,
                        ha="left", va="center",
                        arrowprops=dict(arrowstyle="-", lw=0.4, color="0.6"))
        ax.set_xlim(0, lim * 1.35); ax.set_ylim(0, lim)
        ax.set_xlabel("sender score (1−magnitude_rank)")
        ax.set_ylabel("receiver score")
        ax.set_title(grp, fontsize=9)
    fig.suptitle(f"Sender vs receiver roles — {age}  "
                 f"(size = #active LR pairs; DESCRIPTIVE)", fontsize=10)
    fig.tight_layout()
    out = pdir / f"sender_receiver_bubble_{_slug(age)}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {out.name}")
