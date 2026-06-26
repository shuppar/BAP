#!/usr/bin/env python
"""h09k_plots.py -- figures for the Admati 2x2 PE arm (CSV-only, no recompute).

Resolved biology (two conserved axes):
  eoPE -> HYPOXIA / oxygen-response / glycolysis (UP), broad across compartments;
          shares the Gunter-Rahman conserved hypoxia genes (NDRG1, BNIP3L, ERRFI1, PLIN2...).
  loPE -> OXIDATIVE PHOSPHORYLATION / electron-transport (DOWN) + insulin/peptide-hormone
          response; coordinated mito-ETC shutdown (NDUFB4, CHCHD2, COX subunits...).
GA-matched diagonal NOT stronger (structure is by PE subtype x biology, not stage).

Reads h09k_concordant_pathways_2x2.csv (now with FDR) + h09k_leading_edge_2x2.csv.

Panels:
  1. peak_grid.png          -- 2x2 concordance-peak heatmap, all compartments (diagonal starred).
  2. subtype_diverging.png  -- diverging dotplot: headline pathways, eoPE(+) vs loPE(-) by NES,
                               size = -log10 FDR  (REPLACES the split scatters).
  3. pathway_dotplot_{onset}.png -- top named concordant pathways per subtype, size=-log10 FDR.
  4. genes_{onset}_{pathway}.png -- named shared leading-edge genes (mouse vs human stat).

Usage:  uv run python scripts/h09k_plots.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from h09e_cross_species_rrho import OUT_DIR  # noqa: E402
try:
    from h09_summary_plots import _label_points
except Exception:
    def _label_points(ax, xs, ys, texts, color="black"):
        for x, y, t in zip(xs, ys, texts):
            ax.annotate(t, (x, y), fontsize=6, xytext=(3, 3), textcoords="offset points", color=color)

TAB = OUT_DIR / "tables"
PLOT = OUT_DIR / "plots" / "h09k_admati_2x2"
GA_MATCHED = {("E12.5", "eoPE"), ("E18.5", "loPE")}
COMPS = ["trophoblast", "decidua_stromal", "vascular", "immune"]


def _short(p):
    return (p.replace("HALLMARK_", "H:").replace("REACTOME_", "R:")
            .replace("GOBP_", "GO:").replace("_", " "))[:42]


def _fdr_size(fdr):
    return np.clip(-np.log10(np.asarray(fdr, float).clip(1e-300)) * 14, 20, 320)


def plot_peak_grid(rrho):
    fig, axes = plt.subplots(1, len(COMPS), figsize=(3.5 * len(COMPS), 3.3), squeeze=False)
    for ax, comp in zip(axes[0], COMPS):
        s = rrho[rrho["compartment"] == comp]
        if s.empty:
            ax.axis("off"); continue
        piv = s.pivot(index="mouse_age", columns="human_onset", values="concordance_peak")
        im = ax.imshow(piv.values, cmap="viridis", aspect="auto")
        ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns)
        ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index)
        for i, a in enumerate(piv.index):
            for j, o in enumerate(piv.columns):
                star = "*" if (a, o) in GA_MATCHED else ""
                ax.text(j, i, f"{piv.values[i,j]:.0f}{star}", ha="center", va="center",
                        color="w", fontsize=11, fontweight="bold")
        ax.set_title(comp, fontsize=9)
        ax.set_xlabel("human"); ax.set_ylabel("mouse" if comp == COMPS[0] else "")
    fig.suptitle("Cross-species RRHO concordance peak (-log10 p). *=GA-matched (NOT stronger). "
                 "Peaks are n-sensitive \u2014 compare within column, not across.", fontsize=9)
    fig.tight_layout()
    fig.savefig(PLOT / "peak_grid.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_subtype_diverging(pw):
    """Headline figure: a few signature pathways per axis, eoPE vs loPE, NES on x (size=FDR)."""
    # representative headline pathways for each conserved axis
    eo_paths = ["HALLMARK_HYPOXIA", "GOBP_RESPONSE_TO_OXYGEN_LEVELS",
                "HALLMARK_GLYCOLYSIS", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"]
    lo_paths = ["HALLMARK_OXIDATIVE_PHOSPHORYLATION", "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
                "GOBP_RESPONSE_TO_INSULIN", "GOBP_RESPONSE_TO_PEPTIDE_HORMONE", "REACTOME_TRANSLATION"]
    order = eo_paths + lo_paths
    rows = []
    for p in order:
        for onset, col in [("eoPE", "#c0392b"), ("loPE", "#2471a3")]:
            sub = pw[(pw["pathway"] == p) & (pw["human_onset"] == onset)]
            if not sub.empty:
                # average NES_human + min FDR across the cells where it's concordant
                rows.append(dict(pathway=p, onset=onset, color=col,
                                 nes=sub["NES_human"].mean(),
                                 fdr=sub[["FDR_mouse", "FDR_human"]].max(axis=1).min()))
    d = pd.DataFrame(rows)
    if d.empty:
        print("  [diverging] no rows"); return
    ylab = [p for p in order if p in d["pathway"].values]
    ypos = {p: i for i, p in enumerate(ylab[::-1])}
    fig, ax = plt.subplots(figsize=(8.5, 0.5 * len(ylab) + 1.5))
    for _, r in d.iterrows():
        ax.scatter(r["nes"], ypos[r["pathway"]], s=_fdr_size(r["fdr"]),
                   color=r["color"], alpha=0.85, edgecolor="k", linewidth=0.4,
                   label=r["onset"])
    ax.axvline(0, color="gray", lw=0.8)
    ax.set_yticks(range(len(ylab))); ax.set_yticklabels([_short(p) for p in ylab[::-1]], fontsize=8)
    ax.set_xlabel("human NES  (left = suppressed, right = induced)")
    ax.set_title("Two conserved PE axes vs mouse prenatal stress\n"
                 "eoPE \u2192 hypoxia/glycolysis UP (red);  loPE \u2192 OXPHOS/translation DOWN (blue)\n"
                 "dot size = -log10 FDR", fontsize=10)
    h, l = ax.get_legend_handles_labels()
    seen = dict(zip(l, h)); ax.legend(seen.values(), seen.keys(), loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(PLOT / "subtype_diverging.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_pathway_dotplot(pw, onset, top_n=18):
    s = pw[pw["human_onset"] == onset].copy()
    if s.empty:
        return
    s["fdr"] = s[["FDR_mouse", "FDR_human"]].max(axis=1)
    # collapse to one row per pathway (best (min) FDR, mean NES across cells)
    agg = s.groupby("pathway").agg(nes=("NES_human", "mean"), fdr=("fdr", "min"),
                                   direction=("direction", "first")).reset_index()
    agg["score"] = -np.log10(agg["fdr"].clip(lower=1e-300))
    agg = agg.nlargest(top_n, "score").sort_values("nes")
    fig, ax = plt.subplots(figsize=(7.5, 0.34 * len(agg) + 1.2))
    col = np.where(agg["direction"] == "up_both", "#c0392b", "#2471a3")
    ax.scatter(agg["nes"], range(len(agg)), s=_fdr_size(agg["fdr"]),
               c=col, alpha=0.85, edgecolor="k", linewidth=0.4)
    ax.axvline(0, color="gray", lw=0.7)
    ax.set_yticks(range(len(agg))); ax.set_yticklabels([_short(p) for p in agg["pathway"]], fontsize=7)
    ax.set_xlabel("human NES")
    ax.set_title(f"{onset}: top concordant pathways (red=up_both, blue=down_both; size=-log10 FDR)",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(PLOT / f"pathway_dotplot_{onset}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_gene_panel(le, onset, pathway_filter, fname, title, comp="trophoblast"):
    h = le[(le["human_onset"] == onset) & (le["compartment"] == comp)
           & (le["pathway"].str.contains(pathway_filter))].drop_duplicates("gene")
    if h.empty:
        print(f"  [{onset}/{pathway_filter}] no genes"); return
    h = h.assign(mag=h["mouse_stat"].abs() + h["human_stat"].abs())
    fig, ax = plt.subplots(figsize=(6.2, 5.6))
    ax.scatter(h["mouse_stat"], h["human_stat"], s=24, alpha=0.7, color="#8e44ad")
    top = h.nlargest(18, "mag")
    _label_points(ax, top["mouse_stat"], top["human_stat"], top["gene"])
    ax.axhline(0, color="gray", lw=0.5); ax.axvline(0, color="gray", lw=0.5)
    ax.set_xlabel("mouse Wald stat"); ax.set_ylabel("human Wald stat")
    ax.set_title(f"{title}\n{comp}, n={len(h)} shared leading-edge genes", fontsize=9)
    fig.tight_layout()
    fig.savefig(PLOT / fname, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    PLOT.mkdir(parents=True, exist_ok=True)
    rrho = pd.read_csv(TAB / "h09k_rrho_2x2_summary.csv")
    pw = pd.read_csv(TAB / "h09k_concordant_pathways_2x2.csv")
    le = pd.read_csv(TAB / "h09k_leading_edge_2x2.csv")
    if "FDR_human" not in pw.columns:
        sys.exit("h09k_concordant_pathways_2x2.csv has no FDR columns -- re-run patched h09k first.")

    print("[plots] peak grid"); plot_peak_grid(rrho)
    print("[plots] subtype diverging (headline)"); plot_subtype_diverging(pw)
    print("[plots] pathway dotplots"); 
    plot_pathway_dotplot(pw, "eoPE"); plot_pathway_dotplot(pw, "loPE")
    print("[plots] gene panels")
    plot_gene_panel(le, "eoPE", "HALLMARK_HYPOXIA", "genes_eoPE_HYPOXIA.png",
                    "eoPE shared HYPOXIA genes (vs mouse prenatal stress)")
    plot_gene_panel(le, "loPE", "OXIDATIVE|ELECTRON|RESPIRATORY", "genes_loPE_OXPHOS.png",
                    "loPE shared OXPHOS/ETC genes (vs mouse prenatal stress)")
    print(f"\n[plots] wrote -> {PLOT}")
    for f in sorted(PLOT.glob("*.png")):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
