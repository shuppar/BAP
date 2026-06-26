#!/usr/bin/env python
"""h09_summary_plots.py -- figures for the placenta cross-species arm (CSV-only, no recompute).

Reads the saved h09e-h09h tables and produces Fig 4's placenta panels:
  1. rrho_grid.png        -- 4-compartment RRHO heatmap grid (recomputed from saved ranks),
                             annotated with permutation empirical p + Spearman.
  2. nes_concordance.png  -- mouse NES vs human NES for concordant pathways, HYPOXIA highlighted.
  3. pathway_dotplot_<comp>.png -- top concordant pathways per compartment (dot = -log10 FDR).
  4. hypoxia_leadingedge.png    -- shared HALLMARK_HYPOXIA leading-edge genes: mouse vs human stat.

Mirrors the 8c/8e compute-plot split. rrho_matrix imported verbatim from h09e.

Usage (from project root):
  uv run python scripts/h09_summary_plots.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from h09e_cross_species_rrho import rrho_matrix, OUT_DIR  # noqa: E402

TAB = OUT_DIR / "tables"
PLOT = OUT_DIR / "plots" / "h09_summary"


def _load_stat(path, value_col):
    """Saved per-compartment stat CSV (MultiIndex compartment/gene). Return {comp: Series}."""
    df = pd.read_csv(path)
    df.columns = ["compartment", "gene"] + list(df.columns[2:])
    df = df.rename(columns={df.columns[-1]: value_col})
    return {c: g.set_index("gene")[value_col] for c, g in df.groupby("compartment")}


def plot_rrho_grid(mouse_stat, human_stat, rrho_summ, null_summ):
    comps = [c for c in ["trophoblast", "decidua_stromal", "vascular", "immune"]
             if c in mouse_stat and c in human_stat]
    n = len(comps)
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    axes = axes.ravel()
    emp = dict(zip(null_summ["compartment"], null_summ["empirical_p"])) if not null_summ.empty else {}
    spr = dict(zip(rrho_summ["compartment"], rrho_summ["spearman_r"])) if not rrho_summ.empty else {}
    for ax, comp in zip(axes, comps):
        mat, cut = rrho_matrix(mouse_stat[comp], human_stat[comp])
        if mat is None:
            ax.axis("off"); continue
        im = ax.imshow(mat, origin="lower", aspect="auto", cmap="inferno")
        p = emp.get(comp, np.nan); r = spr.get(comp, np.nan)
        ptxt = "<1e-4" if p == 1e-4 else (f"{p:.1e}" if pd.notna(p) else "NA")
        ax.set_title(f"{comp}\nperm p={ptxt}  Spearman r={r:.2f}", fontsize=10)
        ax.set_xlabel("human obese-vs-lean rank", fontsize=8)
        ax.set_ylabel("mouse Late-vs-Relaxed rank", fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.046, label="-log10 p")
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle("Cross-species RRHO: mouse prenatal-stress vs human obese-stress placenta", fontsize=12)
    fig.tight_layout()
    fig.savefig(PLOT / "rrho_grid.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_nes_concordance(pw):
    if pw.empty:
        return
    fig, ax = plt.subplots(figsize=(7, 6.5))
    comps = pw["compartment"].unique()
    colors = dict(zip(comps, plt.cm.tab10(np.linspace(0, 1, len(comps)))))
    for comp in comps:
        s = pw[pw["compartment"] == comp]
        ax.scatter(s["NES_mouse"], s["NES_human"], s=28, alpha=0.7,
                   color=colors[comp], label=comp, edgecolor="none")
    # highlight + label hypoxia (adjustText repels the colliding labels)
    hyp = pw[pw["source"] == "HALLMARK_HYPOXIA"].reset_index(drop=True)
    for _, h in hyp.iterrows():
        ax.scatter(h["NES_mouse"], h["NES_human"], s=130, facecolor="none",
                   edgecolor="red", linewidth=1.8, zorder=5)
    _label_points(ax, hyp["NES_mouse"], hyp["NES_human"],
                  ["HYPOXIA: " + c for c in hyp["compartment"]], color="red")
    lim = np.nanmax(np.abs([pw["NES_mouse"], pw["NES_human"]])) * 1.1
    ax.axhline(0, color="gray", lw=0.5); ax.axvline(0, color="gray", lw=0.5)
    ax.plot([-lim, lim], [-lim, lim], "--", color="gray", lw=0.6)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("mouse NES (Late-vs-Relaxed)")
    ax.set_ylabel("human NES (obese-vs-lean)")
    ax.set_title("Concordant pathways (FDR<0.05 both species, same sign)")
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(PLOT / "nes_concordance.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_pathway_dotplots(pw, top_n=15):
    for comp in pw["compartment"].unique():
        s = pw[pw["compartment"] == comp].copy()
        s["min_fdr"] = s[["FDR_mouse", "FDR_human"]].max(axis=1)   # the weaker (conservative)
        s["score"] = -np.log10(s["min_fdr"].clip(lower=1e-300))
        s = s.nlargest(top_n, "score").iloc[::-1]
        fig, ax = plt.subplots(figsize=(7, 0.32 * len(s) + 1.2))
        col = np.where(s["direction"] == "up_both", "#c0392b", "#2471a3")
        ax.scatter(s["score"], range(len(s)), s=np.clip(s["score"] * 18, 30, 320),
                   c=col, alpha=0.8, edgecolor="k", linewidth=0.4)
        ax.set_yticks(range(len(s)))
        ax.set_yticklabels([p.replace("HALLMARK_", "H:").replace("REACTOME_", "R:")
                            .replace("GOBP_", "GO:")[:48] for p in s["source"]], fontsize=7)
        ax.set_xlabel("-log10 FDR (conservative, max of both species)")
        ax.set_title(f"{comp}: concordant pathways (red=up_both, blue=down_both)", fontsize=10)
        fig.tight_layout()
        fig.savefig(PLOT / f"pathway_dotplot_{comp}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)


def _label_points(ax, xs, ys, texts, color="black"):
    """Non-overlapping labels with leader lines via adjustText (fallback: plain annotate)."""
    try:
        from adjustText import adjust_text
        objs = [ax.text(x, y, t, fontsize=6, color=color) for x, y, t in zip(xs, ys, texts)]
        adjust_text(objs, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.4),
                    expand=(1.3, 1.6), force_text=(0.4, 0.6))
    except Exception:
        for x, y, t in zip(xs, ys, texts):
            ax.annotate(t, (x, y), fontsize=6, xytext=(3, 3), textcoords="offset points", color=color)


def plot_le_scatter(le, pw, pathway, fname, title):
    """Shared leading-edge genes (mouse vs human stat) for one pathway, across compartments."""
    h = le[(le["pathway"] == pathway) & (le["shared"])].copy()
    if h.empty:
        print(f"  [{pathway}] no shared leading-edge genes")
        return
    comps = h["compartment"].unique()
    fig, axes = plt.subplots(1, len(comps), figsize=(5.8 * len(comps), 5.4), squeeze=False)
    for ax, comp in zip(axes[0], comps):
        s = h[h["compartment"] == comp].assign(
            mag=lambda d: d["mouse_stat"].abs() + d["human_stat"].abs())
        ax.scatter(s["mouse_stat"], s["human_stat"], s=22, alpha=0.7, color="#8e44ad")
        top = s.nlargest(12, "mag")
        _label_points(ax, top["mouse_stat"], top["human_stat"], top["gene"])
        ax.axhline(0, color="gray", lw=0.5); ax.axvline(0, color="gray", lw=0.5)
        ax.set_xlabel("mouse Wald stat"); ax.set_ylabel("human Wald stat")
        ax.set_title(f"{comp}: shared\nleading-edge (n={len(s)})", fontsize=9)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(PLOT / fname, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_top_pathway_le(le, pw):
    """Hypoxia (headline) + the top concordant pathway per remaining compartment."""
    plot_le_scatter(le, pw, "HALLMARK_HYPOXIA", "le_HALLMARK_HYPOXIA.png",
                    "Conserved hypoxia program: shared leading-edge genes (mouse vs human stat)")
    # top concordant pathway per compartment by conservative FDR
    seen = {"HALLMARK_HYPOXIA"}
    for comp in pw["compartment"].unique():
        s = pw[pw["compartment"] == comp].copy()
        s["min_fdr"] = s[["FDR_mouse", "FDR_human"]].max(axis=1)
        for cand in s.sort_values("min_fdr")["source"]:
            if cand in seen:
                continue
            seen.add(cand)
            slug = cand.replace("/", "_")[:40]
            plot_le_scatter(le, pw, cand, f"le_{slug}.png",
                            f"{cand}: shared leading-edge genes ({comp})")
            break


def main():
    PLOT.mkdir(parents=True, exist_ok=True)
    mouse_stat = _load_stat(TAB / "h09e_mouse_stat_humansym.csv", "mouse_stat")
    human_stat = _load_stat(TAB / "h09e_human_stat.csv", "human_stat")
    rrho_summ = pd.read_csv(TAB / "h09e_rrho_summary.csv")
    null_summ = pd.read_csv(TAB / "h09f_permutation_null.csv") if (TAB / "h09f_permutation_null.csv").is_file() else pd.DataFrame()
    pw = pd.read_csv(TAB / "h09g_concordant_pathways.csv")
    le = pd.read_csv(TAB / "h09h_leading_edge.csv")

    print("[plots] RRHO grid")
    plot_rrho_grid(mouse_stat, human_stat, rrho_summ, null_summ)
    print("[plots] NES concordance scatter")
    plot_nes_concordance(pw)
    print("[plots] per-compartment pathway dotplots")
    plot_pathway_dotplots(pw)
    print("[plots] leading-edge scatters (hypoxia + top pathway per compartment)")
    plot_top_pathway_le(le, pw)
    print(f"\n[plots] wrote -> {PLOT}")
    for f in sorted(PLOT.glob("*.png")):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
