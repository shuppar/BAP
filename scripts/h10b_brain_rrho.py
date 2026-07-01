#!/usr/bin/env python
"""h10b_brain_rrho.py -- brain cross-species RRHO engine (dataset-agnostic).

Grid: mouse {early_vs_relaxed, late_vs_relaxed} x {P1, 4W, 3mo} x {whole, Isocortex}
      x 7 broad celltypes   <->   human {disorder vs control} x 7 broad.
Per cell: RRHO (+ permutation null) + concordant GSEA + leading-edge. All cross-species
machinery lifted VERBATIM from the placenta arm (h09e/h09g/h09h); only the I/O wrappers
are new. PATHFINDER dataset = Velmeshev (--dataset velmeshev).

MOUSE side: reuse the existing 08b Wald `stat` (results/brain/tables/08b_de/08b_de_results.csv)
  -- NO recompute, so rankings are identical to the paper DE (Figs 2/3). Filter to the
  pairwise per-age contrasts, bridge mouse->human symbols via refs/mouse_human_orthologs.tsv.
  Isocortex level carries ONLY neurons (ExN/InN) by construction -> T2 = neurons; glia T1-only.

HUMAN side: pseudobulk parquet from h10a (groups x genes) -> PyDESeq2 (~ sex + diagnosis)
  per broad celltype, exactly as h09k's human_rankings (read parquet, DESeq2 directly).

Oli/OPC: mouse broad 'OPC/Oligodendrocytes' is one class; compared against human Oli AND
  human OPC separately (merged-primary). Subcluster MOL/OPC split = a later secondary pass.

Outputs (data/human_validation/brain/<dataset>/tables/):
  h10b_<ds>_rankings.parquet         -- mouse+human ranking vectors (maps/replots reuse)
  h10b_<ds>_rrho_summary.csv         -- per (contrast,age,level,mouse_ct,human_ct): peak,p,class
  h10b_<ds>_concordant_pathways.csv  -- FDR<cut same-sign intersect, with NES+FDR
  h10b_<ds>_leading_edge.csv         -- shared leading-edge genes per concordant pathway

Usage (WS, from project root):
  uv run python scripts/h10b_brain_rrho.py --dataset velmeshev --variant primary \\
      --n-perm 5000 --n-jobs 24
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from _utils import parallel_map  # noqa: E402
from h09e_cross_species_rrho import (  # noqa: E402
    rrho_matrix, classify_rrho_concordance, ORTHO,
)
from h09g_pathways_tf import (  # noqa: E402
    run_gsea_on_ranks, add_fdr, load_genesets_tsv, HUMAN_MSIGDB, COLLECTIONS,
    MIN_GENES, FDR_CUT, run_tf_ulm, load_collectri_human, bh,
)
from h09h_leading_edge import compute_leading_edge  # noqa: E402
# reuse h09k's null verbatim (concordance_peak + permutation_null are dataset-agnostic)
from h09k_admati_2x2 import concordance_peak, permutation_null  # noqa: E402

MARGIN_FRAC = 0.25  # winning quadrant must beat runner-up by >=25% of peak to trust the label


def robust_label(mat, spearman_r):
    """Return (robust_class, concordance_margin, spearman_agrees).

    The argmax `rrho_class` flips on noise when quadrants are near-tied, and can contradict
    the global Spearman sign. We trust the directional label only when BOTH hold:
      (a) winning quadrant beats runner-up by >= MARGIN_FRAC * peak, AND
      (b) Spearman sign agrees: concordant_* expects r>0, discordant expects r<0.
    Otherwise -> 'ambiguous'. (The peak magnitude is still reported + shuffle-tested;
    only the DIRECTION claim is withheld. The RRHO map is the ground truth.)
    """
    if mat is None:
        return "none", 0.0, False
    k = mat.shape[0]; h = k // 2
    quads = {"concordant_up": float(mat[:h, :h].max()),
             "concordant_down": float(mat[h:, h:].max()),
             "discordant": float(max(mat[:h, h:].max(), mat[h:, :h].max()))}
    ordered = sorted(quads.values(), reverse=True)
    peak = ordered[0]
    margin = peak - ordered[1]
    best = max(quads, key=quads.get)
    if peak < 2:
        return "none", round(margin, 2), False
    expect_pos = best.startswith("concordant")          # concordant -> r>0, discordant -> r<0
    spearman_agrees = (spearman_r > 0) if expect_pos else (spearman_r < 0)
    robust = best if (margin >= MARGIN_FRAC * peak and spearman_agrees) else "ambiguous"
    return robust, round(margin, 2), bool(spearman_agrees)

# ----------------------------------------------------------------------------
# Paths + constants
# ----------------------------------------------------------------------------
BRAIN = Path("data/human_validation/brain")
MOUSE_DE = Path("results/brain/tables/08b_de/08b_de_results.csv")
MIN_DONORS = 3  # human pseudobulk inclusion per arm (matches h10a powering floor)

# mouse 08b celltype label -> broad bridge token
MOUSE_CT2BROAD = {
    "Astrocytes/Ependymal": "Ast",
    "Excitatory neurons": "ExN",
    "Inhibitory neurons": "InN",
    "Immune": "Mic",
    "OPC/Oligodendrocytes": "Oli_OPC",   # merged; bridged to human Oli AND OPC
    "Vascular": "Endo",
    # dropped (no homolog in the human sets): 'Dopaminergic neurons',
    # 'Olfactory ensheathing cells'
}
# which human broad types a mouse broad maps to (Oli_OPC -> both Oli and OPC)
MOUSE2HUMAN_BROAD = {
    "Ast": ["Ast"], "ExN": ["ExN"], "InN": ["InN"], "Mic": ["Mic"],
    "Endo": ["Endo"], "Oli_OPC": ["Oli", "OPC"],
}
MOUSE_CONTRASTS = {
    "early_vs_relaxed": "early_vs_relaxed_per_age",
    "late_vs_relaxed": "late_vs_relaxed_per_age",
}
MOUSE_AGES = ["P1", "4W", "3mo"]
MOUSE_LEVELS = ["whole", "Isocortex"]   # T1 whole (all types), T2 Isocortex (neurons only)

# age-matched highlight (biologically appropriate human<->mouse age); a label, not a filter
AGE_MATCH = {
    "velmeshev": ["4W", "3mo"],   # peds/adolescent ASD
    "maitra": ["3mo"], "nagy": ["3mo"],   # adult MDD
    "macnair": ["3mo"],                   # adult MS
    "hwang": ["3mo"],                     # adult PTSD/MDD
}

# per-dataset human config: pseudobulk parquet + meta + condition contrast
DATASETS = {
    "velmeshev": dict(
        sub="velmeshev_2019_autism",
        pb="h10a_velmeshev_pseudobulk_{variant}.parquet",
        meta="h10a_velmeshev_group_meta_{variant}.csv",
        cond_col="diagnosis", test="ASD", ref="Control", covar="sex",
    ),
    "maitra": dict(   # MDD female, dlPFC BA9 (sex constant=F -> covar auto-drops)
        sub="maitra_2023_GSE213982",
        pb="h10c_maitra_pseudobulk_{variant}.parquet",
        meta="h10c_maitra_group_meta_{variant}.csv",
        cond_col="diagnosis", test="MDD", ref="Control", covar="sex",
    ),
    "nagy": dict(     # MDD male, dlPFC BA9 (sex constant=M -> covar auto-drops)
        sub="nagy_2020_GSE144136",
        pb="h10d_nagy_pseudobulk_{variant}.parquet",
        meta="h10d_nagy_group_meta_{variant}.csv",
        cond_col="diagnosis", test="MDD", ref="Control", covar="sex",
    ),
    "macnair": dict(  # MS (stressed-glia reference, NOT etiology); both sexes -> ~ sex + diagnosis
        sub="macnair_2025_MS",
        pb="h10e_macnair_pseudobulk_{variant}.parquet",
        meta="h10e_macnair_group_meta_{variant}.csv",
        cond_col="diagnosis", test="MS", ref="Control", covar="sex",
    ),
    "hwang": dict(    # PTSD + MDD dlPFC (Hwang/Girgenti 2025); 3-level diagnosis {CON,MDD,PTSD}.
        # test is set at runtime by --dx-contrast (PTSD or MDD); reader's isin([ref,test]) drops
        # the third arm -> a PTSD-vs-CON / MDD-vs-CON internal directional control in one dataset.
        sub="hwang_ptsd",
        pb="h10f_hwang_pseudobulk_{variant}.parquet",
        meta="h10f_hwang_group_meta_{variant}.csv",
        cond_col="diagnosis", test=None, ref="CON", covar="sex",
    ),
}


def mouse_rankings(m2h):
    """Read 08b stat, filter to the per-age pairwise grid, bridge to human symbols.

    Returns nested dict: out[contrast][age][level][broad] = Series(stat, index=human symbol).
    Only sex='combined' here (primary). Isocortex yields only ExN/InN (data-enforced).
    """
    print(f"[h10b] reading mouse 08b DE {MOUSE_DE}")
    usecols = ["contrast", "group_level", "sex", "level", "celltype", "gene", "stat"]
    df = pd.read_csv(MOUSE_DE, usecols=usecols, low_memory=False)
    df = df[(df["sex"] == "combined")
            & (df["contrast"].isin(MOUSE_CONTRASTS.values()))
            & (df["level"].isin(MOUSE_LEVELS))
            & (df["celltype"].isin(MOUSE_CT2BROAD))].copy()
    df["broad"] = df["celltype"].map(MOUSE_CT2BROAD)

    out = {}
    for cname, cval in MOUSE_CONTRASTS.items():
        out[cname] = {}
        for age in MOUSE_AGES:
            out[cname][age] = {}
            for level in MOUSE_LEVELS:
                out[cname][age][level] = {}
                sl = df[(df["contrast"] == cval) & (df["group_level"] == age)
                        & (df["level"] == level)]
                for broad, g in sl.groupby("broad"):
                    s = g.dropna(subset=["stat"]).set_index("gene")["stat"]
                    s = s[~s.index.duplicated()]
                    s = s.rename(index=m2h)
                    s = s[~s.index.duplicated()].dropna()
                    if len(s) >= 200:
                        out[cname][age][level][broad] = s
    return out


def human_rankings(ds, variant):
    """Pseudobulk parquet -> PyDESeq2 (~ covar + cond) per broad. Mirrors h09k.human_rankings.

    Returns out[broad] = Series(stat indexed by human symbol).
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    tab = BRAIN / ds["sub"] / "tables"
    pb = pd.read_parquet(tab / ds["pb"].format(variant=variant))
    meta = pd.read_csv(tab / ds["meta"].format(variant=variant), index_col=0)
    cc, test, ref, covar = ds["cond_col"], ds["test"], ds["ref"], ds["covar"]
    print(f"[h10b] human DE ({cc}: {test} vs {ref}, covar={covar})")

    out = {}
    for broad, sel in meta.groupby("broad"):
        sel = sel[sel[cc].isin([ref, test])]
        gc = sel[cc].value_counts()
        if gc.get(ref, 0) < MIN_DONORS or gc.get(test, 0) < MIN_DONORS:
            print(f"    [human/{broad}] too few donors {dict(gc)} -- skipped")
            continue
        cmat = pb.loc[sel.index].astype(int)
        cmat = cmat.loc[:, cmat.sum(axis=0) > 0]
        md = sel[[cc]].copy()
        md[cc] = pd.Categorical(md[cc], categories=[ref, test])
        design = [cc]
        if covar in sel.columns and sel[covar].nunique() > 1:
            md[covar] = sel[covar].values
            design = [covar, cc]
        dds = DeseqDataSet(counts=cmat, metadata=md, design_factors=design, quiet=True)
        dds.deseq2()
        st = DeseqStats(dds, contrast=[cc, test, ref], quiet=True)
        st.summary()
        res = st.results_df.dropna(subset=["stat"])
        out[broad] = res["stat"]
        print(f"    [human/{broad}] {dict(gc)}, {len(res)} genes")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=sorted(DATASETS))
    ap.add_argument("--variant", default="primary", choices=["primary", "sensitivity"])
    ap.add_argument("--n-perm", type=int, default=5000)
    ap.add_argument("--n-jobs", type=int, default=24)
    ap.add_argument("--dx-contrast", choices=["PTSD", "MDD"], default=None,
                    help="hwang only: which arm vs CON (ref). Required for --dataset hwang; "
                         "run once per arm (outputs suffixed _PTSD / _MDD so they don't clobber).")
    ap.add_argument("--tf", action="store_true",
                    help="also compute concordant TF activity (CollecTRI human ULM); "
                         "off by default since placenta showed single-vector ULM nulls")
    args = ap.parse_args()

    ds = dict(DATASETS[args.dataset])   # copy: dx-contrast override must not mutate the registry
    if args.dx_contrast:
        if args.dataset != "hwang":
            sys.exit("--dx-contrast only applies to --dataset hwang")
        ds["test"] = args.dx_contrast
    elif ds.get("test") is None:
        sys.exit(f"--dataset {args.dataset} requires --dx-contrast {{PTSD,MDD}} (3-level diagnosis)")
    out_tag = args.dataset + (f"_{args.dx_contrast}" if args.dx_contrast else "")
    tab = BRAIN / ds["sub"] / "tables"
    tab.mkdir(parents=True, exist_ok=True)
    age_match = set(AGE_MATCH.get(args.dataset, []))

    m2h = dict(pd.read_csv(ORTHO, sep="\t")[["mouse_symbol", "human_symbol"]].values)
    print(f"[h10b] {len(m2h)} mouse->human 1:1 orthologs")

    mouse = mouse_rankings(m2h)
    human = human_rankings(ds, args.variant)
    if not human:
        sys.exit("ERROR: no human rankings produced (all arms too thin?)")

    # GSEA per ranking, computed once and reused across grid cells
    net = load_genesets_tsv(HUMAN_MSIGDB, COLLECTIONS, MIN_GENES)
    coll_map = dict(net[["source", "collection"]].drop_duplicates().values)
    members = (pd.read_csv(HUMAN_MSIGDB, sep="\t")
               .groupby("gs_name")["gene_symbol"].apply(set).to_dict())
    print("[h10b] GSEA per ranking")
    hgsea = {b: add_fdr(run_gsea_on_ranks(human[b], net, MIN_GENES), coll_map) for b in human}
    mgsea = {}
    for cname in mouse:
        for age in mouse[cname]:
            for level in mouse[cname][age]:
                for broad, s in mouse[cname][age][level].items():
                    mgsea[(cname, age, level, broad)] = add_fdr(
                        run_gsea_on_ranks(s, net, MIN_GENES), coll_map)

    # TF activity per ranking (CollecTRI human ULM), once + reused; concordant = FDR<cut both
    # + same activity sign. Off by default (placenta single-vector ULM nulled); --tf to enable.
    htf, mtf = {}, {}
    if args.tf:
        print("[h10b] TF ULM per ranking (CollecTRI human)")
        collectri = load_collectri_human()
        for b in human:
            t = run_tf_ulm(human[b], collectri); t["FDR"] = bh(t["pvalue"]); htf[b] = t
        for cname in mouse:
            for age in mouse[cname]:
                for level in mouse[cname][age]:
                    for broad, s in mouse[cname][age][level].items():
                        t = run_tf_ulm(s, collectri); t["FDR"] = bh(t["pvalue"])
                        mtf[(cname, age, level, broad)] = t

    # persist rankings (maps/replots never recompute DE)
    rank_rows = []
    for cname in mouse:
        for age in mouse[cname]:
            for level in mouse[cname][age]:
                for broad, s in mouse[cname][age][level].items():
                    for g, v in s.items():
                        rank_rows.append({"side": "mouse", "contrast": cname, "age": age,
                                          "level": level, "celltype": broad, "gene": g, "stat": v})
    for broad, s in human.items():
        for g, v in s.items():
            rank_rows.append({"side": "human", "contrast": ds["test"], "age": "human",
                              "level": "human", "celltype": broad, "gene": g, "stat": v})
    pd.DataFrame(rank_rows).to_parquet(tab / f"h10b_{out_tag}_rankings.parquet")
    print(f"[h10b] rankings -> {tab / f'h10b_{out_tag}_rankings.parquet'}")

    rrho_rows, pw_rows, le_rows, tf_rows = [], [], [], []
    for cname in mouse:
        for age in mouse[cname]:
            for level in mouse[cname][age]:
                for m_broad, m_stat in mouse[cname][age][level].items():
                    for h_broad in MOUSE2HUMAN_BROAD[m_broad]:
                        if h_broad not in human:
                            continue
                        h_stat = human[h_broad]
                        mat, cut = rrho_matrix(m_stat, h_stat)
                        if mat is None:
                            continue
                        klass, _ = classify_rrho_concordance(mat, cut)
                        peak = concordance_peak(mat)
                        emp_p, n_null = permutation_null(m_stat, h_stat, peak,
                                                         args.n_perm, args.n_jobs)
                        common = m_stat.index.intersection(h_stat.index)
                        rho, _ = spearmanr(m_stat.loc[common], h_stat.loc[common])
                        robust, margin, sp_ok = robust_label(mat, rho)
                        rrho_rows.append({
                            "contrast": cname, "mouse_age": age, "level": level,
                            "mouse_ct": m_broad, "human_ct": h_broad,
                            "age_matched": age in age_match,
                            "rrho_class": klass, "robust_class": robust,
                            "concordance_margin": margin, "spearman_agrees": sp_ok,
                            "concordance_peak": round(peak, 2),
                            "empirical_p": emp_p, "n_perm": n_null,
                            "spearman_r": round(rho, 3), "n_shared_genes": len(common)})
                        tag = "AGEMATCH" if age in age_match else ""
                        print(f"  [{cname} {age} {level} {m_broad}->{h_broad}] {tag} "
                              f"{robust}({klass}) peak={peak:.1f} margin={margin:.1f} "
                              f"p={emp_p:.1e} r={rho:.2f}")

                        # concordant pathways (intersect precomputed GSEA, same sign)
                        mk = (cname, age, level, m_broad)
                        if mk in mgsea and h_broad in hgsea:
                            mg, hg = mgsea[mk], hgsea[h_broad]
                            mrg = mg.merge(hg, on="source", suffixes=("_mouse", "_human"))
                            conc = mrg[(mrg["FDR_mouse"] < FDR_CUT) & (mrg["FDR_human"] < FDR_CUT)
                                       & (np.sign(mrg["NES_mouse"]) == np.sign(mrg["NES_human"]))]
                            for _, gp in conc.iterrows():
                                direction = "up_both" if gp["NES_mouse"] > 0 else "down_both"
                                pw_rows.append({
                                    "contrast": cname, "mouse_age": age, "level": level,
                                    "mouse_ct": m_broad, "human_ct": h_broad,
                                    "age_matched": age in age_match,
                                    "pathway": gp["source"], "direction": direction,
                                    "NES_mouse": round(gp["NES_mouse"], 3),
                                    "NES_human": round(gp["NES_human"], 3),
                                    "FDR_mouse": gp["FDR_mouse"], "FDR_human": gp["FDR_human"]})
                                mem = members.get(gp["source"], set())
                                m_le = dict(compute_leading_edge(m_stat, mem, gp["NES_mouse"]))
                                h_le = dict(compute_leading_edge(h_stat, mem, gp["NES_human"]))
                                for gene in sorted(set(m_le) & set(h_le)):
                                    le_rows.append({
                                        "contrast": cname, "mouse_age": age, "level": level,
                                        "mouse_ct": m_broad, "human_ct": h_broad,
                                        "pathway": gp["source"], "gene": gene,
                                        "mouse_stat": round(m_le[gene], 3),
                                        "human_stat": round(h_le[gene], 3)})

                        # concordant TFs (FDR<cut both species, same activity sign)
                        if args.tf and mk in mtf and h_broad in htf:
                            tm, th = mtf[mk], htf[h_broad]
                            tmh = tm.merge(th, on="source", suffixes=("_mouse", "_human"))
                            tconc = tmh[(tmh["FDR_mouse"] < FDR_CUT) & (tmh["FDR_human"] < FDR_CUT)
                                        & (np.sign(tmh["activity_score_mouse"])
                                           == np.sign(tmh["activity_score_human"]))]
                            for _, tr in tconc.iterrows():
                                tf_rows.append({
                                    "contrast": cname, "mouse_age": age, "level": level,
                                    "mouse_ct": m_broad, "human_ct": h_broad,
                                    "age_matched": age in age_match,
                                    "TF": tr["source"],
                                    "direction": "up_both" if tr["activity_score_mouse"] > 0
                                                 else "down_both",
                                    "act_mouse": round(tr["activity_score_mouse"], 3),
                                    "act_human": round(tr["activity_score_human"], 3),
                                    "FDR_mouse": tr["FDR_mouse"], "FDR_human": tr["FDR_human"]})

    rrho = pd.DataFrame(rrho_rows)
    pd.DataFrame(rrho_rows).to_csv(tab / f"h10b_{out_tag}_rrho_summary.csv", index=False)
    pd.DataFrame(pw_rows).to_csv(tab / f"h10b_{out_tag}_concordant_pathways.csv", index=False)
    pd.DataFrame(le_rows).to_csv(tab / f"h10b_{out_tag}_leading_edge.csv", index=False)
    if args.tf:
        pd.DataFrame(tf_rows).to_csv(tab / f"h10b_{out_tag}_concordant_tfs.csv", index=False)
        print(f"[h10b] concordant TFs: {len(tf_rows)} -> "
              f"{tab / f'h10b_{out_tag}_concordant_tfs.csv'}")

    print(f"\n[h10b] RRHO summary ({len(rrho)} cells) -> "
          f"{tab / f'h10b_{out_tag}_rrho_summary.csv'}")
    if not rrho.empty:
        print(rrho.sort_values("concordance_peak", ascending=False)
              .head(20)[["contrast", "mouse_age", "level", "mouse_ct", "human_ct",
                         "age_matched", "robust_class", "rrho_class", "concordance_peak",
                         "concordance_margin", "empirical_p", "spearman_r"]]
              .to_string(index=False))
        print("\n[h10b] robust_class census:")
        print(rrho["robust_class"].value_counts().to_string())
        print("\n[h10b] mean concordance peak, age-matched vs not:")
        print(rrho.groupby("age_matched")["concordance_peak"].mean().to_string())
    print(f"[h10b] concordant pathways: {len(pw_rows)} | leading-edge rows: {len(le_rows)}")


if __name__ == "__main__":
    main()
