#!/usr/bin/env python
"""h10_summary_plots.py -- comprehensive brain cross-species figures (all 4 datasets).

Reads ONLY saved h10b outputs (rrho_summary, concordant_pathways, leading_edge,
concordant_tfs, rankings parquet) -- never recomputes DE. CSV/parquet-only, Mac-runnable,
serial (matplotlib not thread-safe).

Filter (locked): a cell is PLOTTED if empirical_p < 0.05 (above chance). robust_class is a
VISUAL encoding, not a gate -- directional (concordant_up/down) drawn solid/boxed, ambiguous
& discordant drawn lighter/flagged. Empty (no above-chance cell) panels are skipped.

Plots:
  1. compact overview      -- max peak per dataset x human-celltype (one small heatmap)
  2. combined master       -- mouse(contrast,age,ct) x dataset.human_ct, peak; sig only
  3. rrho maps             -- rank-rank heatmaps, p<0.05 cells only, per dataset
  4. concordant pathways   -- pathway x sig-cell heatmap, signed NES, per dataset
  5. leading-edge panels   -- shared genes (mouse vs human stat) for top pathways, headline cells
  6. concordant TFs        -- unique TF x sig-cell heatmap, signed activity, per dataset

Usage (Mac or WS, from project root):
  uv run python scripts/h10_summary_plots.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from h09e_cross_species_rrho import rrho_matrix  # noqa: E402
try:
    from h09_summary_plots import _label_points  # reuse per INSTRUCTIONS
except Exception:
    def _label_points(ax, xs, ys, texts, color="black"):
        for x, y, t in zip(xs, ys, texts):
            ax.annotate(t, (x, y), fontsize=6, color=color)

BRAIN = Path("data/human_validation/brain")
SYNTH = BRAIN / "_synthesis" / "plots"
DATASETS = {  # ds -> (subdir, pretty label, disorder)
    "velmeshev": ("velmeshev_2019_autism", "Velmeshev ASD", "ASD"),
    "maitra": ("maitra_2023_GSE213982", "Maitra MDD-F", "MDD (female)"),
    "nagy": ("nagy_2020_GSE144136", "Nagy MDD-M", "MDD (male)"),
    "macnair": ("macnair_2025_MS", "Macnair MS", "MS (stressed-glia ref)"),
}
P_CUT = 0.05
CT_ORDER = ["ExN", "InN", "Ast", "Oli", "OPC", "Mic", "Endo"]
AGE_ORDER = ["P1", "4W", "3mo"]

# Biological-thread regexes (verified against real pathway names in the concordant CSVs).
# These drive the PATHWAY-KEYED views, which are gated on GSEA FDR (the concordant-pathway
# definition) and NOT on RRHO peak -- because peak-based views are dominated by neurons and
# structurally hide small-but-coherent programs like microglial IFN. Two complementary lenses.
THREADS = {
    "IFN/immune": r"INTERFERON|INFLAMMAT|IMMUNE|CYTOKINE|ISG|INNATE|TNF|NFKB|COMPLEMENT|IL[0-9]|JAK_STAT",
    "ECM/mesench": r"MATRIX|COLLAGEN|MESENCHYM|EMT|EPITHELIAL_MESEN|FIBR|INTEGRIN|ADHESION",
    "gliogenesis": r"GLIOGEN|GLIAL|GLIA_DIFF|ASTROCYTE_DIFF|OLIGODENDROCYTE_DIFF|MYELIN",
    "synaptic": r"SYNAP|NEUROTRANSMIT|AXON|DENDRIT",
}
IMMUNE_PAT = THREADS["IFN/immune"]


def tabdir(ds):
    return BRAIN / DATASETS[ds][0] / "tables"


def plotdir(ds):
    d = BRAIN / DATASETS[ds][0] / "plots" / "h10_summary"
    d.mkdir(parents=True, exist_ok=True)
    return d


def load_summary(ds):
    f = tabdir(ds) / f"h10b_{ds}_rrho_summary.csv"
    if not f.exists():
        return None
    df = pd.read_csv(f)
    df["dataset"] = ds
    df["sig"] = df["empirical_p"] < P_CUT
    df["directional"] = df["robust_class"].astype(str).str.startswith("concordant")
    # signed peak: + for concordant_up, - for concordant_down (down-down programs)
    sign = np.where(df["robust_class"] == "concordant_down", -1.0, 1.0)
    df["signed_peak"] = df["concordance_peak"] * sign
    df["cell"] = (df["contrast"].str.replace("_vs_relaxed", "", regex=False)
                  + "/" + df["mouse_age"] + "/" + df["level"]
                  + "/" + df["mouse_ct"] + "\u2192" + df["human_ct"])
    return df


# ---------------------------------------------------------------------------
# 1. compact overview: max peak per dataset x human-celltype
# ---------------------------------------------------------------------------
def plot_overview(all_summ):
    rows = []
    for ds, df in all_summ.items():
        sig = df[df["sig"]]
        for ct in CT_ORDER:
            sub = sig[sig["human_ct"] == ct]
            rows.append({"dataset": DATASETS[ds][1], "human_ct": ct,
                         "max_peak": sub["concordance_peak"].max() if len(sub) else np.nan})
    piv = (pd.DataFrame(rows).pivot(index="human_ct", columns="dataset", values="max_peak")
           .reindex(CT_ORDER))
    fig, ax = plt.subplots(figsize=(1.4 * piv.shape[1] + 2, 0.6 * len(CT_ORDER) + 1.5))
    im = ax.imshow(piv.values, cmap="magma", aspect="auto")
    ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, rotation=30, ha="right")
    ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                        color="white" if v < np.nanmax(piv.values) * 0.6 else "black", fontsize=8)
    fig.colorbar(im, ax=ax, label="max RRHO peak (-log10 p), sig cells only")
    ax.set_title("Brain cross-species: strongest concordance per cell type\n"
                 f"(max over mouse grid; empirical p < {P_CUT}; blank = none above chance)",
                 fontsize=10)
    fig.tight_layout()
    SYNTH.mkdir(parents=True, exist_ok=True)
    fig.savefig(SYNTH / "01_overview_maxpeak.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] overview -> {SYNTH / '01_overview_maxpeak.png'}")


# ---------------------------------------------------------------------------
# 2. combined master: mouse(contrast,age,ct) x dataset.human_ct, signed peak
# ---------------------------------------------------------------------------
def plot_master(all_summ):
    big = pd.concat([df[df["sig"]] for df in all_summ.values()], ignore_index=True)
    if big.empty:
        print("[plot] master: no sig cells"); return
    big["mouse_key"] = (big["contrast"].str.replace("_vs_relaxed", "", regex=False)
                        + " " + big["mouse_age"] + " " + big["level"] + " " + big["mouse_ct"])
    big["col_key"] = big["dataset"].map(lambda d: DATASETS[d][1]) + " : " + big["human_ct"]
    piv = big.pivot_table(index="mouse_key", columns="col_key", values="signed_peak",
                          aggfunc="first")
    # order rows by max |signed peak|
    piv = piv.reindex(piv.abs().max(axis=1).sort_values(ascending=False).index)
    dirpiv = big.pivot_table(index="mouse_key", columns="col_key", values="directional",
                             aggfunc="first").reindex(piv.index)[piv.columns]

    vmax = np.nanmax(np.abs(piv.values))
    fig, ax = plt.subplots(figsize=(0.5 * piv.shape[1] + 4, 0.32 * piv.shape[0] + 2))
    im = ax.imshow(piv.values, cmap="RdBu_r", aspect="auto", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index, fontsize=6)
    # box directional cells
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            if bool(dirpiv.values[i, j]) and np.isfinite(piv.values[i, j]):
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                           edgecolor="black", lw=1.2))
    fig.colorbar(im, ax=ax, label="signed RRHO peak (+up/-down)")
    ax.set_title("Brain cross-species master: mouse prenatal-stress \u00d7 human disorder\n"
                 f"(sig cells, p<{P_CUT}; black box = directionally robust; red=up-up, blue=down-down)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(SYNTH / "02_master_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] master -> {SYNTH / '02_master_heatmap.png'}  ({piv.shape[0]}x{piv.shape[1]})")


# ---------------------------------------------------------------------------
# 3. rrho maps: p<0.05 cells, per dataset
# ---------------------------------------------------------------------------
def load_rankings(ds):
    rk = pd.read_parquet(tabdir(ds) / f"h10b_{ds}_rankings.parquet")
    mouse, human = {}, {}
    for (c, a, l, ct), g in rk[rk.side == "mouse"].groupby(["contrast", "age", "level", "celltype"]):
        mouse[(c, a, l, ct)] = g.set_index("gene")["stat"]
    for ct, g in rk[rk.side == "human"].groupby("celltype"):
        human[ct] = g.set_index("gene")["stat"]
    return mouse, human


def plot_maps(ds, summ):
    sig = summ[summ["sig"]].sort_values("concordance_peak", ascending=False)
    if sig.empty:
        print(f"[plot] {ds} maps: no sig cells -- skipped"); return
    mouse, human = load_rankings(ds)
    n = len(sig)
    nc = min(5, n); nr = int(np.ceil(n / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(2.4 * nc, 2.6 * nr), squeeze=False)
    for idx, (_, r) in enumerate(sig.iterrows()):
        ax = axes[idx // nc][idx % nc]
        mk = (r["contrast"], r["mouse_age"], r["level"], r["mouse_ct"])
        if mk not in mouse or r["human_ct"] not in human:
            ax.axis("off"); continue
        mat, _ = rrho_matrix(mouse[mk], human[r["human_ct"]])
        if mat is None:
            ax.axis("off"); continue
        # PER-PANEL scaling: each map to its own max so secondary cells show their structure
        # (the peak value in the title carries the cross-panel magnitude comparison; the
        # absolute scale lives in the summary/overview heatmaps).
        ax.imshow(mat, cmap="viridis", origin="upper", aspect="auto",
                  vmin=0, vmax=float(mat.max()))
        clr = "black" if r["directional"] else "darkorange"
        flag = "" if r["directional"] else " [amb]"
        ax.set_title(f"{r['mouse_ct']}\u2192{r['human_ct']} {r['mouse_age']}/{r['level']}\n"
                     f"peak {r['concordance_peak']:.0f} p={r['empirical_p']:.0e}{flag}\n"
                     f"{r['robust_class']}", fontsize=6.5, color=clr)
        ax.set_xticks([]); ax.set_yticks([])
    for k in range(n, nr * nc):
        axes[k // nc][k % nc].axis("off")
    fig.suptitle(f"{DATASETS[ds][1]} -- RRHO maps (p<{P_CUT}; PER-PANEL scaled; top-left=shared UP, "
                 f"bottom-right=shared DOWN; orange title=ambiguous direction)", fontsize=10, y=1.0)
    fig.tight_layout()
    out = plotdir(ds) / f"h10_{ds}_rrho_maps_sig.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] {ds} maps ({n} sig cells) -> {out}")


# ---------------------------------------------------------------------------
# 4. concordant pathways: pathway x sig-cell heatmap, signed NES
# ---------------------------------------------------------------------------
def plot_pathways(ds, summ):
    f = tabdir(ds) / f"h10b_{ds}_concordant_pathways.csv"
    if not f.exists():
        return
    pw = pd.read_csv(f)
    sig_cells = summ[summ["sig"]][["contrast", "mouse_age", "level", "mouse_ct", "human_ct"]]
    pw = pw.merge(sig_cells, on=["contrast", "mouse_age", "level", "mouse_ct", "human_ct"])
    if pw.empty:
        print(f"[plot] {ds} pathways: none in sig cells -- skipped"); return
    pw["cell"] = (pw["mouse_age"] + "/" + pw["level"].str[:4] + "/"
                  + pw["mouse_ct"] + "\u2192" + pw["human_ct"])
    pw["signed_nes"] = pw[["NES_mouse", "NES_human"]].mean(axis=1)
    # top pathways by recurrence x effect
    rank = (pw.groupby("pathway")
            .agg(n=("cell", "nunique"), eff=("signed_nes", lambda s: s.abs().mean()))
            .assign(score=lambda d: d["n"] * d["eff"]).sort_values("score", ascending=False))
    top = rank.head(25).index
    sub = pw[pw["pathway"].isin(top)]
    piv = sub.pivot_table(index="pathway", columns="cell", values="signed_nes", aggfunc="mean") \
        .reindex(top)
    if piv.empty or piv.shape[1] == 0:
        print(f"[plot] {ds} pathways: empty pivot -- skipped"); return
    vmax = np.nanmax(np.abs(piv.values))
    fig, ax = plt.subplots(figsize=(0.5 * piv.shape[1] + 5, 0.34 * piv.shape[0] + 2))
    im = ax.imshow(piv.values, cmap="RdBu_r", aspect="auto", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, rotation=90, fontsize=6.5)
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels([p.replace("HALLMARK_", "H:").replace("REACTOME_", "R:")
                        .replace("GOBP_", "GO:")[:48] for p in piv.index], fontsize=6.5)
    fig.colorbar(im, ax=ax, label="mean NES (+up / -down both species)")
    ax.set_title(f"{DATASETS[ds][1]} -- concordant pathways \u00d7 sig cells "
                 f"(FDR<0.05 both, same sign; top 25 by recurrence\u00d7effect)", fontsize=9)
    fig.tight_layout()
    out = plotdir(ds) / f"h10_{ds}_concordant_pathways.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] {ds} pathways ({len(top)} pw x {piv.shape[1]} cells) -> {out}")


# ---------------------------------------------------------------------------
# 5. leading-edge gene panels: shared genes (mouse vs human stat) for top pathways
#    in the strongest sig cells
# ---------------------------------------------------------------------------
def plot_leading_edge(ds, summ):
    f = tabdir(ds) / f"h10b_{ds}_leading_edge.csv"
    if not f.exists():
        return
    sig = summ[summ["sig"] & summ["directional"]].sort_values("concordance_peak", ascending=False)
    if sig.empty:
        print(f"[plot] {ds} leading-edge: no directional sig cells -- skipped"); return
    le = pd.read_csv(f, low_memory=False)
    # restrict to top-4 directional cells; one panel each (most-conserved pathway's genes)
    cells = sig.head(4)
    n = len(cells)
    fig, axes = plt.subplots(1, n, figsize=(3.4 * n, 3.4), squeeze=False)
    plotted = 0
    for j, (_, r) in enumerate(cells.iterrows()):
        ax = axes[0][j]
        cl = le[(le.contrast == r["contrast"]) & (le.mouse_age == r["mouse_age"])
                & (le.level == r["level"]) & (le.mouse_ct == r["mouse_ct"])
                & (le.human_ct == r["human_ct"])]
        if cl.empty:
            ax.axis("off"); continue
        # pick the pathway with the most shared genes in this cell
        top_pw = cl["pathway"].value_counts().idxmax()
        g = cl[cl.pathway == top_pw].copy()
        ax.scatter(g["mouse_stat"], g["human_stat"], s=14, alpha=0.6,
                   color="#2c7fb8", edgecolor="none")
        ax.axhline(0, color="grey", lw=0.5); ax.axvline(0, color="grey", lw=0.5)
        # label the genes farthest from origin (top ~12)
        g["mag"] = np.hypot(g["mouse_stat"], g["human_stat"])
        lab = g.nlargest(min(12, len(g)), "mag")
        _label_points(ax, lab["mouse_stat"].values, lab["human_stat"].values,
                      lab["gene"].values, color="black")
        ax.set_xlabel("mouse Wald stat", fontsize=8)
        ax.set_ylabel("human Wald stat", fontsize=8)
        ax.set_title(f"{r['mouse_ct']}\u2192{r['human_ct']} {r['mouse_age']}\n"
                     f"{top_pw.replace('HALLMARK_','H:').replace('REACTOME_','R:').replace('GOBP_','GO:')[:34]}\n"
                     f"({len(g)} shared LE genes)", fontsize=7.5)
        plotted += 1
    if plotted == 0:
        plt.close(fig); print(f"[plot] {ds} leading-edge: nothing plotted -- skipped"); return
    fig.suptitle(f"{DATASETS[ds][1]} -- leading-edge gene conservation "
                 f"(shared genes in top directional cells; both axes = signed Wald stat)",
                 fontsize=10, y=1.04)
    fig.tight_layout()
    out = plotdir(ds) / f"h10_{ds}_leading_edge.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] {ds} leading-edge ({plotted} panels) -> {out}")


# ---------------------------------------------------------------------------
# 6. concordant TFs: unique TF x sig-cell heatmap, signed activity
# ---------------------------------------------------------------------------
def plot_tfs(ds, summ):
    f = tabdir(ds) / f"h10b_{ds}_concordant_tfs.csv"
    if not f.exists():
        print(f"[plot] {ds} TFs: file missing -- skipped"); return
    tf = pd.read_csv(f)
    sig_cells = summ[summ["sig"]][["contrast", "mouse_age", "level", "mouse_ct", "human_ct"]]
    tf = tf.merge(sig_cells, on=["contrast", "mouse_age", "level", "mouse_ct", "human_ct"])
    if tf.empty:
        print(f"[plot] {ds} TFs: none in sig cells -- skipped"); return
    tf["cell"] = (tf["mouse_age"] + "/" + tf["level"].str[:4] + "/"
                  + tf["mouse_ct"] + "\u2192" + tf["human_ct"])
    tf["signed_act"] = tf[["act_mouse", "act_human"]].mean(axis=1)
    rank = (tf.groupby("TF").agg(n=("cell", "nunique"),
                                 eff=("signed_act", lambda s: s.abs().mean()))
            .assign(score=lambda d: d["n"] * d["eff"]).sort_values("score", ascending=False))
    top = rank.head(25).index
    piv = (tf[tf.TF.isin(top)].pivot_table(index="TF", columns="cell",
                                           values="signed_act", aggfunc="mean").reindex(top))
    if piv.empty or piv.shape[1] == 0:
        print(f"[plot] {ds} TFs: empty pivot -- skipped"); return
    vmax = np.nanmax(np.abs(piv.values))
    fig, ax = plt.subplots(figsize=(0.5 * piv.shape[1] + 4, 0.32 * piv.shape[0] + 2))
    im = ax.imshow(piv.values, cmap="PRGn", aspect="auto", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, rotation=90, fontsize=6.5)
    ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index, fontsize=7)
    fig.colorbar(im, ax=ax, label="mean TF activity (+/- both species)")
    ax.set_title(f"{DATASETS[ds][1]} -- concordant TFs \u00d7 sig cells "
                 f"(CollecTRI human ULM, FDR<0.05 both, same sign; top 25)", fontsize=9)
    fig.tight_layout()
    out = plotdir(ds) / f"h10_{ds}_concordant_tfs.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] {ds} TFs ({len(top)} TF x {piv.shape[1]} cells) -> {out}")


def _load_all_pathways():
    """Concat every dataset's concordant_pathways.csv (already GSEA-FDR<0.05 both, same sign)."""
    frames = []
    for ds in DATASETS:
        f = tabdir(ds) / f"h10b_{ds}_concordant_pathways.csv"
        if f.exists():
            d = pd.read_csv(f); d["dataset"] = ds
            frames.append(d)
    if not frames:
        return pd.DataFrame()
    pw = pd.concat(frames, ignore_index=True)
    pw["signed_nes"] = pw[["NES_mouse", "NES_human"]].mean(axis=1)
    pw["disorder"] = pw["dataset"].map(lambda d: DATASETS[d][1])
    return pw


# ---------------------------------------------------------------------------
# 7. THREAD SCANNER -- 8f/8g spine (IFN / ECM / gliogenesis / synaptic) across all
#    datasets, GSEA-FDR-gated, NO peak filter (the peak-independent lens).
# ---------------------------------------------------------------------------
def plot_thread_scanner(pw):
    if pw.empty:
        print("[plot] thread scanner: no pathways"); return
    rows = []
    for thread, pat in THREADS.items():
        sub = pw[pw["pathway"].str.contains(pat, case=False, regex=True)]
        for (disorder, ct), g in sub.groupby(["disorder", "human_ct"]):
            rows.append({"thread": thread, "disorder": disorder, "human_ct": ct,
                         "mean_nes": g["signed_nes"].mean(), "n": len(g)})
    if not rows:
        print("[plot] thread scanner: no thread hits"); return
    sc = pd.DataFrame(rows)
    threads = list(THREADS)
    fig, axes = plt.subplots(1, len(threads), figsize=(3.2 * len(threads), 3.6),
                             squeeze=False)
    vmax = np.nanmax(np.abs(sc["mean_nes"]))
    for j, thread in enumerate(threads):
        ax = axes[0][j]
        s = sc[sc["thread"] == thread]
        piv = s.pivot_table(index="human_ct", columns="disorder", values="mean_nes") \
            .reindex(CT_ORDER)
        npiv = s.pivot_table(index="human_ct", columns="disorder", values="n") \
            .reindex(CT_ORDER)
        im = ax.imshow(piv.values, cmap="RdBu_r", aspect="auto", vmin=-vmax, vmax=vmax)
        ax.set_xticks(range(piv.shape[1]))
        ax.set_xticklabels(piv.columns, rotation=90, fontsize=7)
        ax.set_yticks(range(len(CT_ORDER))); ax.set_yticklabels(CT_ORDER, fontsize=8)
        # annotate n (number of concordant pathways behind each cell)
        for a in range(piv.shape[0]):
            for b in range(piv.shape[1]):
                v, nn = piv.values[a, b], npiv.values[a, b]
                if np.isfinite(v):
                    ax.text(b, a, f"{int(nn)}", ha="center", va="center", fontsize=6,
                            color="white" if abs(v) > vmax * 0.5 else "black")
        ax.set_title(thread, fontsize=9)
    fig.colorbar(axes[0][-1].images[0], ax=axes, fraction=0.02, pad=0.02,
                 label="mean NES (+up/-down both species)")
    fig.suptitle("8f/8g thread scanner in human cortex -- concordant pathways by celltype "
                 "(GSEA FDR<0.05 both, same sign; NOT peak-gated; n = #pathways)", fontsize=11, y=1.04)
    SYNTH.mkdir(parents=True, exist_ok=True)
    fig.savefig(SYNTH / "03_thread_scanner.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] thread scanner -> {SYNTH / '03_thread_scanner.png'}")


# ---------------------------------------------------------------------------
# 8. Mic/IFN -- microglia-focused immune pathway heatmap across datasets x age
# ---------------------------------------------------------------------------
def plot_ifn_microglia(pw):
    imm = pw[pw["pathway"].str.contains(IMMUNE_PAT, case=False, regex=True)
             & (pw["human_ct"] == "Mic")].copy()
    if imm.empty:
        print("[plot] Mic/IFN: no microglial immune pathways"); return
    imm["col"] = imm["disorder"] + "\n" + imm["mouse_age"]
    # top immune pathways by recurrence across datasets
    top = imm["pathway"].value_counts().head(25).index
    sub = imm[imm["pathway"].isin(top)]
    piv = sub.pivot_table(index="pathway", columns="col", values="signed_nes", aggfunc="mean") \
        .reindex(top)
    if piv.shape[1] == 0:
        print("[plot] Mic/IFN: empty"); return
    vmax = np.nanmax(np.abs(piv.values))
    fig, ax = plt.subplots(figsize=(0.55 * piv.shape[1] + 6, 0.36 * piv.shape[0] + 2))
    im = ax.imshow(piv.values, cmap="RdBu_r", aspect="auto", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels([p.replace("HALLMARK_", "H:").replace("REACTOME_", "R:")
                        .replace("GOBP_", "GO:")[:50] for p in piv.index], fontsize=7)
    fig.colorbar(im, ax=ax, label="mean NES (+up / -down both species)")
    ax.set_title("MICROGLIA -- conserved immune/IFN programs (mouse prenatal-stress \u00d7 human)\n"
                 "down in MDD/ASD (perinatal IFN co-suppression, the 8f/8g thread); "
                 "UP in MS (neuroinflammation) -- the directional inversion", fontsize=9)
    fig.tight_layout()
    fig.savefig(SYNTH / "04_microglia_ifn.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] Mic/IFN microglia -> {SYNTH / '04_microglia_ifn.png'}")


# ---------------------------------------------------------------------------
# 9. Mic/IFN -- IFN/immune across ALL celltypes (Mic/Endo/ExN carry it)
# ---------------------------------------------------------------------------
def plot_ifn_allcelltypes(pw):
    imm = pw[pw["pathway"].str.contains(IMMUNE_PAT, case=False, regex=True)].copy()
    if imm.empty:
        print("[plot] IFN all-ct: none"); return
    # mean NES per (celltype x disorder), collapsing pathways+ages -> the cross-celltype view
    g = (imm.groupby(["human_ct", "disorder"])
         .agg(mean_nes=("signed_nes", "mean"), n=("pathway", "size")).reset_index())
    piv = g.pivot(index="human_ct", columns="disorder", values="mean_nes").reindex(CT_ORDER)
    npiv = g.pivot(index="human_ct", columns="disorder", values="n").reindex(CT_ORDER)
    vmax = np.nanmax(np.abs(piv.values))
    fig, ax = plt.subplots(figsize=(1.5 * piv.shape[1] + 2, 0.6 * len(CT_ORDER) + 1.5))
    im = ax.imshow(piv.values, cmap="RdBu_r", aspect="auto", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, rotation=20, ha="right")
    ax.set_yticks(range(len(CT_ORDER))); ax.set_yticklabels(CT_ORDER)
    for a in range(piv.shape[0]):
        for b in range(piv.shape[1]):
            v, nn = piv.values[a, b], npiv.values[a, b]
            if np.isfinite(v):
                ax.text(b, a, f"{v:+.2f}\n(n={int(nn)})", ha="center", va="center",
                        fontsize=7, color="white" if abs(v) > vmax * 0.5 else "black")
    fig.colorbar(im, ax=ax, label="mean NES, immune/IFN pathways")
    ax.set_title("IFN/immune concordance across ALL celltypes\n"
                 "(mean NES of immune pathways; n = #concordant pathways; "
                 "Mic/Endo/ExN carry the thread)", fontsize=9)
    fig.tight_layout()
    fig.savefig(SYNTH / "05_ifn_all_celltypes.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] IFN all-celltypes -> {SYNTH / '05_ifn_all_celltypes.png'}")


def main():
    all_summ = {}
    for ds in DATASETS:
        s = load_summary(ds)
        if s is None:
            print(f"[plot] {ds}: no summary -- skipped"); continue
        all_summ[ds] = s
        nsig = int(s["sig"].sum()); ndir = int((s["sig"] & s["directional"]).sum())
        print(f"[plot] {ds}: {nsig} sig cells ({ndir} directional)")

    if not all_summ:
        sys.exit("no datasets with summaries found")

    # cross-dataset (peak-keyed)
    plot_overview(all_summ)
    plot_master(all_summ)
    # cross-dataset (pathway-keyed, peak-independent -- the lens that surfaces 8f/8g threads)
    pw_all = _load_all_pathways()
    plot_thread_scanner(pw_all)
    plot_ifn_microglia(pw_all)
    plot_ifn_allcelltypes(pw_all)
    # per-dataset
    for ds, summ in all_summ.items():
        plot_maps(ds, summ)
        plot_pathways(ds, summ)
        plot_leading_edge(ds, summ)
        plot_tfs(ds, summ)
    print(f"\n[plot] cross-dataset figures -> {SYNTH}")
    print("[plot] per-dataset figures -> data/human_validation/brain/<ds>/plots/h10_summary/")


if __name__ == "__main__":
    main()
