#!/usr/bin/env python
"""h09k_admati_2x2.py -- 2x2 cross-species PE validation (Admati sc, all compartments).

Grid: mouse {E12.5 Early-vs-Relaxed, E18.5 Late-vs-Relaxed} x human {eoPE, loPE} x 4
compartments (trophoblast, decidua_stromal, vascular, immune). GA-matched DIAGONAL
(eoPE<->E12.5, loPE<->E18.5) is the headline; off-diagonal is the cross-check. The
stressor-timing confound on the mouse axis (E12.5=Early, E18.5=Late) is ACCEPTED + stated.

Per 2x2-cell x compartment: RRHO (+ permutation null) + concordant GSEA + leading-edge.
Everything lifted from the existing placenta arm (h09e/h09g/h09h) -- same methods.

Human DE from the h09j pseudobulk parquet (author-annotated sc, no SoupX -- flagged).
Mouse DE recomputed per age via h09e.pseudobulk_de (~ sex + group).

Usage (from project root):
  uv run python scripts/h09k_admati_2x2.py --n-perm 5000 --n-jobs 24
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from _utils import unassigned_mask, parallel_map  # noqa: E402
from h09e_cross_species_rrho import (  # noqa: E402
    pseudobulk_de, rrho_matrix, classify_rrho_concordance,
    MOUSE_H5AD, ORTHO, MAP_YAML, OUT_DIR,
)
from h09g_pathways_tf import (  # noqa: E402
    run_gsea_on_ranks, add_fdr, load_genesets_tsv, HUMAN_MSIGDB, COLLECTIONS,
    MIN_GENES, FDR_CUT,
)
from h09h_leading_edge import compute_leading_edge  # noqa: E402

PB_PATH = OUT_DIR / "tables" / "h09j_admati_pseudobulk.parquet"
META_PATH = OUT_DIR / "tables" / "h09j_admati_group_meta.csv"
TAB = OUT_DIR / "tables"
PLOT = OUT_DIR / "plots" / "h09k_admati_2x2"

MOUSE_ARMS = {"E12.5": dict(test="Early_Stress", ref="Relaxed"),
              "E18.5": dict(test="Late_Stress", ref="Relaxed")}
HUMAN_ARMS = {"eoPE": dict(test="early_PE", ref="early_control"),
              "loPE": dict(test="late_PE", ref="late_control")}
GA_MATCHED = {("E12.5", "eoPE"), ("E18.5", "loPE")}   # diagonal
MIN_DONORS = 2


# ---- RRHO concordance peak + permutation null (concordant quadrant, either sign) ----
def concordance_peak(mat):
    k = mat.shape[0]; h = k // 2
    return float(max(mat[:h, :h].max(), mat[h:, h:].max()))


def _shuffle_chunk(args):
    m_vals, h_vals, n_sh, seed = args
    rng = np.random.default_rng(seed)
    idx = pd.RangeIndex(len(m_vals))
    m_ser = pd.Series(m_vals, index=idx)
    peaks = np.empty(n_sh)
    for s in range(n_sh):
        mat, _ = rrho_matrix(m_ser, pd.Series(rng.permutation(h_vals), index=idx))
        peaks[s] = concordance_peak(mat) if mat is not None else 0.0
    return peaks


def permutation_null(m_stat, h_stat, obs_peak, n_perm, n_jobs):
    common = m_stat.index.intersection(h_stat.index)
    m_vals, h_vals = m_stat.loc[common].to_numpy(), h_stat.loc[common].to_numpy()
    n_chunks = max(n_jobs, 1) * 2
    base = n_perm // n_chunks
    jobs = [(m_vals, h_vals, base + (1 if k < n_perm - base * n_chunks else 0), 7000 + k)
            for k in range(n_chunks)]
    jobs = [j for j in jobs if j[2] > 0]
    peaks = []
    for _j, res, err in parallel_map(_shuffle_chunk, jobs, n_jobs=n_jobs,
                                     use_threads=False, desc="null"):
        if err:
            print(f"    [warn] null chunk: {err.splitlines()[-1]}"); continue
        peaks.append(res)
    null = np.concatenate(peaks) if peaks else np.array([])
    emp_p = (1 + int((null >= obs_peak).sum())) / (1 + len(null)) if len(null) else np.nan
    return emp_p, len(null)


# ---- DE ----
def mouse_rankings(m2h, comps):
    print("[h09k] mouse rankings (E12.5 Early-vs-Relaxed, E18.5 Late-vs-Relaxed)")
    mo = sc.read_h5ad(MOUSE_H5AD)
    mo = mo[~unassigned_mask(mo.obs, ["celltype_majority"])].copy()
    cmap = yaml.safe_load(MAP_YAML.read_text())["placenta_compartments"]["mouse"]
    mo.obs["compartment"] = mo.obs["celltype_majority"].astype(str).map(cmap)
    if "counts" not in mo.layers:
        mo.layers["counts"] = mo.X.copy()
    out = {}
    for age, spec in MOUSE_ARMS.items():
        sub = mo[(mo.obs["age"] == age) & mo.obs["compartment"].isin(comps)].copy()
        de = pseudobulk_de(sub, "compartment", "donor_id", "group",
                           spec["ref"], spec["test"], "sex", f"mouse/{age}")
        out[age] = {c: s.rename(index=m2h)[~s.rename(index=m2h).index.duplicated()].dropna()
                    for c, s in de.items()}
    del mo
    return out


def human_rankings(comps):
    print("[h09k] human rankings (eoPE vs early_control, loPE vs late_control)")
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    pb = pd.read_parquet(PB_PATH)
    meta = pd.read_csv(META_PATH, index_col=0)
    out = {}
    for onset, spec in HUMAN_ARMS.items():
        out[onset] = {}
        for comp in comps:
            sel = meta[(meta["compartment"] == comp)
                       & (meta["condition"].isin([spec["ref"], spec["test"]]))]
            gc = sel["condition"].value_counts()
            if gc.get(spec["ref"], 0) < MIN_DONORS or gc.get(spec["test"], 0) < MIN_DONORS:
                print(f"    [human/{onset}/{comp}] too few donors {dict(gc)} -- skipped")
                continue
            cmat = pb.loc[sel.index].astype(int)
            cmat = cmat.loc[:, cmat.sum(axis=0) > 0]
            md = sel[["condition"]].copy()
            md["condition"] = pd.Categorical(md["condition"],
                                             categories=[spec["ref"], spec["test"]])
            dds = DeseqDataSet(counts=cmat, metadata=md, design_factors=["condition"], quiet=True)
            dds.deseq2()
            st = DeseqStats(dds, contrast=["condition", spec["test"], spec["ref"]], quiet=True)
            st.summary()
            res = st.results_df.dropna(subset=["stat"])
            out[onset][comp] = res["stat"]
            print(f"    [human/{onset}/{comp}] {dict(gc)}, {len(res)} genes")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=5000)
    ap.add_argument("--n-jobs", type=int, default=24)
    args = ap.parse_args()
    PLOT.mkdir(parents=True, exist_ok=True)

    cmapy = yaml.safe_load(MAP_YAML.read_text())["placenta_compartments"]
    comps = cmapy["rrho_compartments"]
    m2h = dict(pd.read_csv(ORTHO, sep="\t")[["mouse_symbol", "human_symbol"]].values)

    mouse = mouse_rankings(m2h, comps)
    human = human_rankings(comps)

    # GSEA per side x arm x compartment (computed once, reused across the 2x2 cells)
    net = load_genesets_tsv(HUMAN_MSIGDB, COLLECTIONS, MIN_GENES)
    coll_map = dict(net[["source", "collection"]].drop_duplicates().values)
    print("[h09k] GSEA per ranking")
    mgsea, hgsea = {}, {}
    for age in MOUSE_ARMS:
        for comp in comps:
            if comp in mouse[age]:
                mgsea[(age, comp)] = add_fdr(run_gsea_on_ranks(mouse[age][comp], net, MIN_GENES), coll_map)
    for onset in HUMAN_ARMS:
        for comp in comps:
            if comp in human[onset]:
                hgsea[(onset, comp)] = add_fdr(run_gsea_on_ranks(human[onset][comp], net, MIN_GENES), coll_map)

    members = pd.read_csv(HUMAN_MSIGDB, sep="\t").groupby("gs_name")["gene_symbol"].apply(set).to_dict()

    # persist the rankings (so h09k_rrho_maps / replots never recompute DE)
    rank_rows = []
    for age in MOUSE_ARMS:
        for comp in mouse[age]:
            for g, v in mouse[age][comp].items():
                rank_rows.append({"side": "mouse", "arm": age, "compartment": comp, "gene": g, "stat": v})
    for onset in HUMAN_ARMS:
        for comp in human[onset]:
            for g, v in human[onset][comp].items():
                rank_rows.append({"side": "human", "arm": onset, "compartment": comp, "gene": g, "stat": v})
    pd.DataFrame(rank_rows).to_parquet(TAB / "h09k_rankings.parquet")
    print(f"[h09k] rankings saved -> {TAB/'h09k_rankings.parquet'}")

    rrho_rows, pw_rows, le_rows = [], [], []
    for age in MOUSE_ARMS:
        for onset in HUMAN_ARMS:
            ga = (age, onset) in GA_MATCHED
            for comp in comps:
                if comp not in mouse[age] or comp not in human[onset]:
                    continue
                m, h = mouse[age][comp], human[onset][comp]
                mat, cut = rrho_matrix(m, h)
                if mat is None:
                    continue
                klass, _ = classify_rrho_concordance(mat, cut)
                peak = concordance_peak(mat)
                emp_p, n_null = permutation_null(m, h, peak, args.n_perm, args.n_jobs)
                common = m.index.intersection(h.index)
                rho, _ = spearmanr(m.loc[common], h.loc[common])
                rrho_rows.append({"mouse_age": age, "human_onset": onset, "compartment": comp,
                                  "ga_matched": ga, "rrho_class": klass,
                                  "concordance_peak": round(peak, 2), "empirical_p": emp_p,
                                  "n_perm": n_null, "spearman_r": round(rho, 3),
                                  "n_shared_genes": len(common)})
                print(f"  [{age} x {onset} x {comp}] {'DIAG' if ga else 'off'} "
                      f"{klass} peak={peak:.1f} p={emp_p:.1e} r={rho:.2f}")

                # concordant pathways (intersect precomputed GSEA, same sign, FDR<cut)
                if (age, comp) in mgsea and (onset, comp) in hgsea:
                    mg, hg = mgsea[(age, comp)], hgsea[(onset, comp)]
                    mrg = mg.merge(hg, on="source", suffixes=("_mouse", "_human"))
                    conc = mrg[(mrg["FDR_mouse"] < FDR_CUT) & (mrg["FDR_human"] < FDR_CUT)
                               & (np.sign(mrg["NES_mouse"]) == np.sign(mrg["NES_human"]))]
                    for _, g in conc.iterrows():
                        direction = "up_both" if g["NES_mouse"] > 0 else "down_both"
                        pw_rows.append({"mouse_age": age, "human_onset": onset,
                                        "compartment": comp, "ga_matched": ga,
                                        "pathway": g["source"], "direction": direction,
                                        "NES_mouse": round(g["NES_mouse"], 3),
                                        "NES_human": round(g["NES_human"], 3),
                                        "FDR_mouse": g["FDR_mouse"], "FDR_human": g["FDR_human"]})
                        # leading edge (shared) for this pathway
                        mem = members.get(g["source"], set())
                        m_le = dict(compute_leading_edge(m, mem, g["NES_mouse"]))
                        h_le = dict(compute_leading_edge(h, mem, g["NES_human"]))
                        shared = set(m_le) & set(h_le)
                        for gene in sorted(shared):
                            le_rows.append({"mouse_age": age, "human_onset": onset,
                                            "compartment": comp, "pathway": g["source"],
                                            "gene": gene, "mouse_stat": round(m_le[gene], 3),
                                            "human_stat": round(h_le[gene], 3)})

    rrho = pd.DataFrame(rrho_rows)
    pw = pd.DataFrame(pw_rows)
    le = pd.DataFrame(le_rows)
    rrho.to_csv(TAB / "h09k_rrho_2x2_summary.csv", index=False)
    pw.to_csv(TAB / "h09k_concordant_pathways_2x2.csv", index=False)
    le.to_csv(TAB / "h09k_leading_edge_2x2.csv", index=False)

    # 2x2 peak heatmap per compartment (rows=mouse age, cols=human onset)
    for comp in comps:
        s = rrho[rrho["compartment"] == comp]
        if s.empty:
            continue
        piv = s.pivot(index="mouse_age", columns="human_onset", values="concordance_peak")
        pv = s.pivot(index="mouse_age", columns="human_onset", values="empirical_p")
        fig, ax = plt.subplots(figsize=(4.2, 3.6))
        im = ax.imshow(piv.values, cmap="viridis", aspect="auto")
        ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns)
        ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index)
        for i, a in enumerate(piv.index):
            for j, o in enumerate(piv.columns):
                star = "*" if (a, o) in GA_MATCHED else ""
                ax.text(j, i, f"{piv.values[i,j]:.1f}{star}\np={pv.values[i,j]:.0e}",
                        ha="center", va="center", color="w", fontsize=8)
        ax.set_title(f"{comp}: RRHO concordance peak\n(*=GA-matched diagonal)", fontsize=9)
        ax.set_xlabel("human onset"); ax.set_ylabel("mouse age")
        fig.colorbar(im, ax=ax, label="-log10 p peak")
        fig.tight_layout()
        fig.savefig(PLOT / f"peak_2x2_{comp}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    print(f"\n[h09k] RRHO 2x2 -> {TAB/'h09k_rrho_2x2_summary.csv'}")
    print(rrho.to_string(index=False))
    print(f"\n[h09k] diagonal vs off-diagonal mean concordance peak:")
    print(rrho.groupby("ga_matched")["concordance_peak"].mean().to_string())
    print(f"[h09k] concordant pathways -> {TAB/'h09k_concordant_pathways_2x2.csv'} ({len(pw)})")
    print(f"[h09k] leading edge -> {TAB/'h09k_leading_edge_2x2.csv'} ({len(le)})")
    print(f"[h09k] plots -> {PLOT}")


if __name__ == "__main__":
    main()
