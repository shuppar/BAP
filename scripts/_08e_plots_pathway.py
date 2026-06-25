"""
_08e_plots_pathway.py — per-pathway cell-cell-communication graphs.

One chord + one network PER stress pathway (one file each), restricting LR pairs
to that pathway's 8c leading-edge genes. This is the noise-reduction figure:
instead of one hairball, ~10 readable per-pathway graphs.

Two arms feed it (both where available):
  baseline      — Δ(stress − relaxed) interaction score per cell-type edge.
                  DESCRIPTIVE (pooled cells). Works brain (whole+regional) + placenta.
                  Significance is attributed to the PATHWAY (8c FDR<0.05), not the edge.
  differential  — aggregated signed interaction_stat over FDR<0.05 LR pairs.
                  FDR-backed. Placenta only (brain differential is null).

Pathway → gene set: from config/stress_pathways_8e.yaml `graph_pathways` × the 8c
leading-edge genes (FDR<0.05, level=='whole') of that pathway, per tissue. An LR
pair is in-pathway if its ligand OR receptor is in that gene set.

Encoding (shared with _08e_plots_stats so the visual grammar is consistent):
  taper/arrow = source→target; width = |Δ score| (baseline) or |mean stat| (diff);
  alpha       = breadth (#contributing LR pairs);
  colour      = up(red)/down(blue) in the stress group.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Patch
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd

UP, DN = "#d73027", "#4575b4"


def _slug(s):
    return str(s).replace(" ", "_").replace("/", "-").replace(".", "").replace("HALLMARK_", "")


# ---------------------------------------------------------------------------
# Pathway gene-set loader (8c leading edge × config whitelist)
# ---------------------------------------------------------------------------

def load_pathway_genesets(cfg, tissue, spec, le_path):
    """Return {pathway_name: set(genes)} for the tissue's graph_pathways, from the
    8c leading-edge CSV (FDR<0.05, level=='whole'). Reads only needed columns."""
    want = set(spec.get(tissue, {}).get("graph_pathways", []))
    if not want:
        return {}
    cols = ["level", "pathway", "pathway_FDR", "gene"]
    le = pd.read_csv(le_path, usecols=cols, low_memory=False)
    le = le[(le.level == "whole") & (le.pathway_FDR < 0.05) & (le.pathway.isin(want))]
    return {pw: set(g.gene.unique()) for pw, g in le.groupby("pathway")}


def _in_pathway_mask(df, genes, lig_col="ligand_complex", rec_col="receptor_complex"):
    gs = set(genes)
    lig = df[lig_col].astype(str).str.split("_")
    rec = df[rec_col].astype(str).str.split("_")
    return np.array([any(x in gs for x in l) or any(x in gs for x in r)
                     for l, r in zip(lig, rec)])


# ---------------------------------------------------------------------------
# Edge aggregation
# ---------------------------------------------------------------------------

def _baseline_delta_edges(baseline_df, genes, age, level, test_group, ref_group,
                          spec_fdr=0.05, quantile=0.0):
    """Per (source,target): Δ = mean(test score) − mean(ref score), over in-pathway
    LR pairs. Score = (1 − magnitude_rank) so higher = stronger signalling.

    Edge-filter (field-standard, CellPhoneDB-style): keep only LR pairs whose
    specificity is significant (specificity_fdr <= spec_fdr) in EITHER the test or
    ref group. `quantile` adds a slice-specific effect floor: drop per-pair Δ below
    that percentile of |Δ| within this slice before aggregating to edges."""
    b = baseline_df[(baseline_df["age"].astype(str) == str(age)) &
                    (baseline_df.get("level", "whole").astype(str) == str(level))].copy()
    if b.empty or "magnitude_rank" not in b.columns:
        return pd.DataFrame()
    b = b[_in_pathway_mask(b, genes)]
    if b.empty:
        return pd.DataFrame()
    if "specificity_fdr" in b.columns:
        b = b[b["group"].isin([test_group, ref_group])].copy()
        b["_spec_ok"] = b["specificity_fdr"] <= spec_fdr
        key = ["source", "target", "ligand_complex", "receptor_complex"]
        ok = b.groupby(key)["_spec_ok"].transform("any")  # significant in EITHER group
        b = b[ok]
        if b.empty:
            return pd.DataFrame()
    b["score"] = 1.0 - b["magnitude_rank"].astype(float)
    # per-pair Δ across the two groups, for the slice-specific floor
    pp = (b.pivot_table(index=["source", "target", "ligand_complex", "receptor_complex"],
                        columns="group", values="score", aggfunc="mean"))
    if test_group not in pp.columns or ref_group not in pp.columns:
        return pd.DataFrame()
    pp["d"] = pp[test_group].fillna(0) - pp[ref_group].fillna(0)
    pp = pp[pp["d"].abs() > 1e-9]
    if pp.empty:
        return pd.DataFrame()
    if quantile and quantile > 0:
        floor = pp["d"].abs().quantile(quantile)
        pp = pp[pp["d"].abs() >= floor]
        if pp.empty:
            return pd.DataFrame()
    g = pp.reset_index().groupby(["source", "target"])
    out = pd.DataFrame({"delta": g["d"].mean(), "n_pairs": g["d"].size()}).reset_index()
    return out[out["delta"].abs() > 1e-9]


def _differential_edges(diff_df, genes, contrast_name, age, fdr=0.05, quantile=0.0):
    d = diff_df[(diff_df["contrast_name"] == contrast_name) &
                (diff_df["age"].astype(str) == str(age))].copy()
    d = d.dropna(subset=["interaction_stat", "interaction_padj"])
    d = d[d["interaction_padj"] < fdr]
    if d.empty:
        return pd.DataFrame()
    d = d[_in_pathway_mask(d, genes)]
    if d.empty:
        return pd.DataFrame()
    if quantile and quantile > 0:
        floor = d["interaction_stat"].abs().quantile(quantile)
        d = d[d["interaction_stat"].abs() >= floor]
        if d.empty:
            return pd.DataFrame()
    agg = (d.groupby(["source", "target"])
           .agg(delta=("interaction_stat", "sum"),       # Σ signed stat (magnitude moved)
                n_pairs=("interaction_stat", "size"))
           .reset_index())
    return agg


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------

def _circle_pos(nodes):
    ang = {n: 2 * np.pi * i / len(nodes) for i, n in enumerate(nodes)}
    return {n: np.array([np.cos(a), np.sin(a)]) for n, a in ang.items()}


def _ribbon(S, T, w_s, w_t, bend=0.45):
    S = np.asarray(S, float); T = np.asarray(T, float)
    ts = np.array([-S[1], S[0]]); ts /= (np.linalg.norm(ts) or 1)
    tt = np.array([-T[1], T[0]]); tt /= (np.linalg.norm(tt) or 1)
    Sl, Sr = S + w_s * ts, S - w_s * ts
    Tl, Tr = T + w_t * tt, T - w_t * tt
    cS, cT = S * bend, T * bend
    verts = [Sl, cS, cT, Tr, Tl, cT, cS, Sr, Sl]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.CLOSEPOLY]
    return MplPath([tuple(v) for v in verts], codes)


def _draw_graph(edges, title, subtitle, out_png, kind="chord"):
    """edges: DataFrame[source,target,delta,n_pairs]. kind in {chord,network}."""
    if edges.empty:
        return False
    nodes = sorted(set(edges["source"]) | set(edges["target"]))
    pos = _circle_pos(nodes)
    n = len(nodes)
    cmap = {nd: plt.cm.tab20(i / max(n - 1, 1)) for i, nd in enumerate(nodes)}
    dmax = edges["delta"].abs().max() or 1.0
    nmax = edges["n_pairs"].max() or 1

    fig, ax = plt.subplots(figsize=(9, 9))
    for _, e in edges.sort_values("delta", key=lambda s: s.abs(), ascending=False).iterrows():
        s, t = e["source"], e["target"]
        col = UP if e["delta"] > 0 else DN
        alpha = 0.25 + 0.6 * e["n_pairs"] / nmax
        if s == t:
            r = 0.04 + 0.08 * abs(e["delta"]) / dmax
            ax.add_patch(plt.Circle(pos[s] * 1.12, r, facecolor=cmap[s],
                                    edgecolor=col, lw=1.4, alpha=alpha, zorder=2))
            continue
        if kind == "chord":
            w_t = 0.02 + 0.13 * abs(e["delta"]) / dmax
            ax.add_patch(PathPatch(_ribbon(pos[s], pos[t], 0.012, w_t),
                                   facecolor=cmap[s], edgecolor=col, lw=1.2,
                                   alpha=alpha, zorder=2))
        else:
            lw = 0.6 + 5.0 * abs(e["delta"]) / dmax
            ax.annotate("", xy=pos[t], xytext=pos[s],
                        arrowprops=dict(arrowstyle="-|>", color=col, lw=lw, alpha=alpha,
                                        connectionstyle="arc3,rad=0.12",
                                        shrinkA=14, shrinkB=14))
    for nd, P in pos.items():
        ax.scatter([P[0]], [P[1]], s=320, c=[cmap[nd]], edgecolors="k",
                   linewidths=1.0, zorder=4)
        ax.text(P[0] * 1.18, P[1] * 1.18, nd, ha="center", va="center",
                fontsize=7, zorder=5)
    ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.4, 1.4)
    ax.set_aspect("equal"); ax.axis("off")
    ax.legend(handles=[Patch(facecolor="0.7", edgecolor=UP, label="up in stress"),
                       Patch(facecolor="0.7", edgecolor=DN, label="down in stress")],
              fontsize=8, loc="lower right", frameon=False)
    ax.set_title(f"{title}\n{subtitle}", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def plot_by_pathway(genesets, baseline_df, differential_df, contrasts, pdir,
                    q_delta=0.0, q_stat=0.0):
    """For each pathway × contrast × level, draw baseline Δ chord+network and (if
    available) differential chord+network. contrasts: list of dicts with keys
    name, test_group, ref_group, age, label. q_delta/q_stat = slice-specific effect
    floors for the baseline / differential edge builders."""
    made = 0
    for pw, genes in genesets.items():
        pslug = _slug(pw)
        for c in contrasts:
            # ---- baseline (descriptive) ----
            levels = sorted(baseline_df.get("level", pd.Series(["whole"])).astype(str).unique()) \
                if not baseline_df.empty else []
            for lvl in levels:
                e = _baseline_delta_edges(baseline_df, genes, c["age"], lvl,
                                          c["test_group"], c["ref_group"],
                                          quantile=q_delta)
                if e.empty:
                    continue
                sub = (f"baseline Δ(score) {c['label']} | {c['age']} | level={lvl} | "
                       f"{len(e)} edges  —  pathway 8c-FDR<0.05; edges DESCRIPTIVE")
                d = pdir / "baseline" / lvl
                d.mkdir(parents=True, exist_ok=True)
                for kind in ("chord", "network"):
                    if _draw_graph(e, f"{pw} — {kind}", sub,
                                   d / f"{pslug}_{kind}_{c['label']}_{c['age']}.png", kind):
                        made += 1
            # ---- differential (FDR-backed; placenta) ----
            if differential_df is not None and not differential_df.empty:
                e = _differential_edges(differential_df, genes, c["name"], c["age"],
                                        quantile=q_stat)
                if not e.empty:
                    sub = (f"differential Σ(stat) {c['label']} | {c['age']} | "
                           f"{len(e)} edges  —  FDR<0.05 LR pairs in pathway")
                    d = pdir / "differential"
                    d.mkdir(parents=True, exist_ok=True)
                    for kind in ("chord", "network"):
                        if _draw_graph(e, f"{pw} — {kind}", sub,
                                       d / f"{pslug}_{kind}_{c['label']}_{c['age']}.png", kind):
                            made += 1
    print(f"  by-pathway graphs written: {made}")
    return made


def _draw_graph_on_ax(ax, edges, title):
    """Network graph onto a provided ax (for side-by-side companion panels)."""
    if edges.empty:
        ax.text(0.5, 0.5, "no edges", ha="center", va="center", transform=ax.transAxes,
                fontsize=10, color="0.5")
        ax.set_title(title, fontsize=10); ax.axis("off")
        return
    nodes = sorted(set(edges["source"]) | set(edges["target"]))
    pos = _circle_pos(nodes)
    n = len(nodes)
    cmap = {nd: plt.cm.tab20(i / max(n - 1, 1)) for i, nd in enumerate(nodes)}
    dmax = edges["delta"].abs().max() or 1.0
    nmax = edges["n_pairs"].max() or 1
    for _, e in edges.iterrows():
        if e["source"] == e["target"]:
            continue
        col = UP if e["delta"] > 0 else DN
        lw = 0.6 + 5.0 * abs(e["delta"]) / dmax
        alpha = 0.25 + 0.6 * e["n_pairs"] / nmax
        ax.annotate("", xy=pos[e["target"]], xytext=pos[e["source"]],
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=lw, alpha=alpha,
                                    connectionstyle="arc3,rad=0.12", shrinkA=12, shrinkB=12))
    for nd, P in pos.items():
        ax.scatter([P[0]], [P[1]], s=240, c=[cmap[nd]], edgecolors="k",
                   linewidths=1.0, zorder=4)
        ax.text(P[0] * 1.2, P[1] * 1.2, nd, ha="center", va="center", fontsize=6, zorder=5)
    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.45, 1.45)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=10)


def plot_cross_scheme_companion(genesets, baseline_broad, baseline_subtype,
                                contrasts, pdir, spec_fdr=0.05, q_delta=0.0):
    """Per pathway × contrast (whole level): broad scheme (left) vs subtype scheme
    (right), side by side, so you can see whether a pathway's signalling localizes
    to a focal substate (e.g. Hofbauer vs whole Myeloid). Baseline Δ, descriptive."""
    made = 0
    for pw, genes in genesets.items():
        pslug = _slug(pw)
        for c in contrasts:
            eb = _baseline_delta_edges(baseline_broad, genes, c["age"], "whole",
                                       c["test_group"], c["ref_group"], spec_fdr,
                                       quantile=q_delta)
            es = _baseline_delta_edges(baseline_subtype, genes, c["age"], "whole",
                                       c["test_group"], c["ref_group"], spec_fdr,
                                       quantile=q_delta)
            if eb.empty and es.empty:
                continue
            fig, axes = plt.subplots(1, 2, figsize=(18, 9))
            _draw_graph_on_ax(axes[0], eb, f"broad  ({len(eb)} edges)")
            _draw_graph_on_ax(axes[1], es, f"subtype  ({len(es)} edges)")
            fig.suptitle(f"{pw} — broad vs subtype  |  {c['label']} | {c['age']}\n"
                         f"baseline Δ(score), specificity-filtered (either group); "
                         f"red=up / blue=down in stress; DESCRIPTIVE (pathway 8c-FDR<0.05)",
                         fontsize=11)
            d = pdir / "cross_scheme"
            d.mkdir(parents=True, exist_ok=True)
            fig.tight_layout()
            fig.savefig(d / f"{pslug}_{c['label']}_{c['age']}.png", dpi=150, bbox_inches="tight")
            plt.close(fig)
            made += 1
    print(f"  cross-scheme companion figures written: {made}")
    return made



# ---------------------------------------------------------------------------
# Per-pathway LR-pair detail: DOTPLOT + ranked LOLLIPOP
#   arm='differential' (whole only)  → value=interaction_stat, sig=interaction_padj
#   arm='baseline'      (whole+regional) → value=Δ score, sig=specificity_fdr (either)
# ---------------------------------------------------------------------------

def _lr_label(row):
    l = str(row.get("ligand", row.get("ligand_complex", "?")))
    r = str(row.get("receptor", row.get("receptor_complex", "?")))
    return f"{l}→{r}"


def _pathway_lr_table(df, genes, arm, contrast, fdr, spec_fdr, level=None,
                      quantile=0.0):
    """Return a tidy per-(LR pair × cell-type pair) table for one pathway × contrast,
    columns: lr, ctp, source, target, value, size.
    arm='differential': value=interaction_stat, size=-log10(interaction_padj), FDR<fdr.
    arm='baseline':     value=Δ score (test−ref), size=|Δ|, specificity-sig (either grp).

    SLICE-SPECIFIC FLOOR: after building this exact slice's significant pairs, drop
    those whose |value| is below the `quantile`-th percentile OF THIS SLICE. The
    cutoff is computed from the pairs this plot actually draws — not a global value."""
    d = df[_in_pathway_mask(df, genes)].copy()
    if d.empty:
        return pd.DataFrame()
    if arm == "differential":
        d = d[(d["contrast_name"] == contrast["name"]) &
              (d["age"].astype(str) == str(contrast["age"]))]
        d = d.dropna(subset=["interaction_stat", "interaction_padj"])
        d = d[d["interaction_padj"] < fdr]
        if d.empty:
            return pd.DataFrame()
        d["value"] = d["interaction_stat"].astype(float)
        if quantile and quantile > 0:
            floor = d["value"].abs().quantile(quantile)
            d = d[d["value"].abs() >= floor]
        if d.empty:
            return pd.DataFrame()
        d["lr"] = d.apply(_lr_label, axis=1)
        d["ctp"] = d["source"].astype(str) + "→" + d["target"].astype(str)
        d["size"] = -np.log10(d["interaction_padj"].clip(lower=1e-300))
        return d[["lr", "ctp", "source", "target", "value", "size"]]
    # baseline arm
    tg, rg = contrast["test_group"], contrast["ref_group"]
    b = d[(d["age"].astype(str) == str(contrast["age"])) &
          (d["group"].isin([tg, rg]))].copy()
    if level is not None and "level" in b.columns:
        b = b[b["level"].astype(str) == str(level)]
    if b.empty or "magnitude_rank" not in b.columns:
        return pd.DataFrame()
    key = ["source", "target", "ligand_complex", "receptor_complex"]
    if "specificity_fdr" in b.columns:
        b["_ok"] = b["specificity_fdr"] <= spec_fdr
        b = b[b.groupby(key)["_ok"].transform("any")]
        if b.empty:
            return pd.DataFrame()
    b["score"] = 1.0 - b["magnitude_rank"].astype(float)
    pp = b.pivot_table(index=key, columns="group", values="score", aggfunc="mean")
    if tg not in pp.columns or rg not in pp.columns:
        return pd.DataFrame()
    pp["value"] = pp[tg].fillna(0) - pp[rg].fillna(0)
    pp = pp[pp["value"].abs() > 1e-9].reset_index()
    if pp.empty:
        return pd.DataFrame()
    if quantile and quantile > 0:
        floor = pp["value"].abs().quantile(quantile)
        pp = pp[pp["value"].abs() >= floor]
        if pp.empty:
            return pd.DataFrame()
    pp["lr"] = pp.apply(_lr_label, axis=1)
    pp["ctp"] = pp["source"].astype(str) + "→" + pp["target"].astype(str)
    pp["size"] = pp["value"].abs()
    return pp[["lr", "ctp", "source", "target", "value", "size"]]


def plot_pathway_lr_dotplot(df, genesets, contrasts, pdir, arm="differential",
                            level=None, fdr=0.05, spec_fdr=0.05, top_n=30,
                            quantile=0.0):
    """Per pathway × contrast dot plot: rows=LR pairs, cols=cell-type pairs,
    colour=sign of effect (red up/blue down in stress), size=significance/|effect|."""
    lvl_slug = f"_{_slug(level)}" if level is not None else ""
    arm_tag = ("FDR<%g; size=-log10(FDR)" % fdr if arm == "differential"
               else "DESCRIPTIVE Δ score; specificity-filtered; size=|Δ|")
    made = 0
    for pw, genes in genesets.items():
        pslug = _slug(pw)
        for c in contrasts:
            d = _pathway_lr_table(df, genes, arm, c, fdr, spec_fdr, level, quantile=quantile)
            if d.empty:
                continue
            top_lr = (d.groupby("lr")["size"].max().sort_values(ascending=False)
                      .head(top_n).index.tolist())
            d = d[d["lr"].isin(top_lr)]
            lrs = d.groupby("lr")["size"].max().sort_values().index.tolist()
            ctps = sorted(d["ctp"].unique())
            xi = {x: i for i, x in enumerate(ctps)}
            yi = {y: i for i, y in enumerate(lrs)}
            smax = d["size"].max() or 1.0
            fig, ax = plt.subplots(figsize=(max(6, len(ctps) * 0.4 + 4),
                                            max(4, len(lrs) * 0.34 + 2)))
            for _, r in d.iterrows():
                col = UP if r["value"] > 0 else DN
                ax.scatter(xi[r["ctp"]], yi[r["lr"]], s=30 + 220 * r["size"] / smax,
                           c=col, edgecolors="k", linewidths=0.3, alpha=0.85)
            ax.set_xticks(range(len(ctps))); ax.set_xticklabels(ctps, rotation=90, fontsize=6)
            ax.set_yticks(range(len(lrs))); ax.set_yticklabels(lrs, fontsize=7)
            lvl_t = f" | level={level}" if level is not None else ""
            ax.set_title(f"{pw} — LR pairs changed | {c['label']} | {c['age']}{lvl_t}\n"
                         f"{arm_tag}; red=up/blue=down in stress", fontsize=9)
            ax.margins(0.04)
            d_out = pdir / "lr_dotplot"; d_out.mkdir(parents=True, exist_ok=True)
            fig.tight_layout()
            fig.savefig(d_out / f"{pslug}_{arm}_{c['label']}_{c['age']}{lvl_slug}.png",
                        dpi=150, bbox_inches="tight")
            plt.close(fig); made += 1
    print(f"  per-pathway LR dotplots ({arm}{lvl_slug}) written: {made}")
    return made


def plot_pathway_lr_ranked(df, genesets, contrasts, pdir, arm="differential",
                           level=None, fdr=0.05, spec_fdr=0.05, top_n=25,
                           quantile=0.0):
    """Per pathway × contrast horizontal lollipop: top-N LR pairs by |effect|,
    labelled 'Ligand→Receptor [source→target]', stem red(up)/blue(down) in stress."""
    lvl_slug = f"_{_slug(level)}" if level is not None else ""
    xlab = ("interaction_stat" if arm == "differential" else "Δ score (test−ref)")
    made = 0
    for pw, genes in genesets.items():
        pslug = _slug(pw)
        for c in contrasts:
            d = _pathway_lr_table(df, genes, arm, c, fdr, spec_fdr, level, quantile=quantile)
            if d.empty:
                continue
            d["lab"] = d.apply(lambda r: f"{r['lr']}  [{r['source']}→{r['target']}]", axis=1)
            d["absval"] = d["value"].abs()
            d = d.sort_values("absval", ascending=False).head(top_n).sort_values("value")
            y = np.arange(len(d))
            cols = [UP if v > 0 else DN for v in d["value"]]
            fig, ax = plt.subplots(figsize=(9, max(3, len(d) * 0.32 + 1.5)))
            ax.hlines(y, 0, d["value"], color=cols, lw=2, alpha=0.8)
            ax.scatter(d["value"], y, c=cols, s=45, edgecolors="k", linewidths=0.4, zorder=3)
            ax.axvline(0, color="k", lw=0.6)
            ax.set_yticks(y); ax.set_yticklabels(d["lab"], fontsize=7)
            lvl_t = f" | level={level}" if level is not None else ""
            ax.set_xlabel(f"{xlab}  (← down | up → in stress)")
            ax.set_title(f"{pw} — top {len(d)} changed LR pairs | {c['label']} | {c['age']}{lvl_t}\n"
                         f"ranked by |{xlab}|", fontsize=9)
            d_out = pdir / "lr_ranked"; d_out.mkdir(parents=True, exist_ok=True)
            fig.tight_layout()
            fig.savefig(d_out / f"{pslug}_{arm}_{c['label']}_{c['age']}{lvl_slug}.png",
                        dpi=150, bbox_inches="tight")
            plt.close(fig); made += 1
    print(f"  per-pathway ranked LR lollipops ({arm}{lvl_slug}) written: {made}")
    return made
