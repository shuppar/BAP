#!/usr/bin/env python
"""h09g_pathways_tf.py -- cross-species concordant GSEA + TF activity for the placenta RRHO.

Two single-species runs per compartment, then intersect (mirrors 8c machinery exactly;
two-GSEA-intersect design locked with the user):
  GSEA: fgsea-multilevel (run_fgsea.R) on each species' compartment Wald-stat ranking,
        per-collection BH FDR (add_fdr). Concordant pathway = FDR<0.05 in BOTH species,
        same NES sign.
  TF:   decoupler ULM (CollecTRI HUMAN network) on each species' ranking, analytical
        t-stat p (run_tf_ulm), BH within species. Concordant TF = FDR<0.05 in BOTH,
        same activity-score sign.

GSEA/TF helper functions lifted VERBATIM from 08c_pathways.py (mouse->human swaps only:
human msigdb TSV + CollecTRI organism='human'). DE re-run via h09e's pseudobulk_de.

Usage (from project root):
  uv run python scripts/h09g_pathways_tf.py
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy.stats import false_discovery_control

sys.path.insert(0, str(Path(__file__).parent))
from _utils import unassigned_mask  # noqa: E402
from h09e_cross_species_rrho import (  # noqa: E402
    pseudobulk_de, MOUSE_H5AD, HUMAN_H5AD, ORTHO, MAP_YAML, OUT_DIR,
)

HUMAN_MSIGDB = Path("refs/msigdb_human.tsv")
COLLECTIONS = ["H", "C2", "C5"]
MIN_GENES, MAX_GENES = 15, 500
FDR_CUT = 0.05


# ============================================================================
# Lifted VERBATIM from 08c_pathways.py (human swaps noted). Do not diverge.
# ============================================================================
def load_genesets_tsv(path, collections, min_genes):
    df = pd.read_csv(path, sep="\t")
    needed = {"collection", "gs_name", "gene_symbol"}
    if not needed.issubset(df.columns):
        sys.exit(f"ERROR: {path} missing columns {needed - set(df.columns)}.")
    if collections:
        df = df[df["collection"].isin(collections)]
    net = (df.rename(columns={"gs_name": "source", "gene_symbol": "target"})
             [["source", "target", "collection"]].drop_duplicates())
    sizes = net.groupby("source").size()
    net = net[net["source"].isin(sizes[sizes >= min_genes].index)]
    return net.reset_index(drop=True)


def load_collectri_human():  # 08c load_collectri, organism swapped mouse->human
    import decoupler as dc
    if hasattr(dc, "op") and hasattr(dc.op, "collectri"):
        return dc.op.collectri(organism="human")
    return dc.get_collectri(organism="human", split_complexes=False)


def run_gsea_on_ranks(rank_series, net, min_genes, seed=42, max_genes=500):
    if rank_series is None or rank_series.empty:
        return pd.DataFrame(columns=["source", "NES", "pvalue"])
    worker = Path(__file__).resolve().parent / "run_fgsea.R"
    if not worker.is_file():
        raise FileNotFoundError(f"fgsea worker not found at {worker}")
    with tempfile.TemporaryDirectory(prefix="fgsea_") as td:
        td = Path(td)
        ranks_path, pathways_path, out_path = td / "ranks.tsv", td / "pw.tsv", td / "out.tsv"
        pd.DataFrame({"gene": rank_series.index.astype(str),
                      "stat": pd.to_numeric(rank_series.values, errors="coerce")}
                     ).to_csv(ranks_path, sep="\t", index=False)
        (net[["source", "target"]].rename(columns={"source": "pathway_name", "target": "gene"})
            .to_csv(pathways_path, sep="\t", index=False))
        cmd = ["Rscript", str(worker), str(ranks_path), str(pathways_path), str(out_path),
               str(int(min_genes)), str(int(max_genes)), str(int(seed))]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
        if proc.returncode != 0:
            raise RuntimeError(f"fgsea worker failed:\n{proc.stderr[-1500:]}")
        if (not out_path.exists()) or out_path.stat().st_size == 0:
            return pd.DataFrame(columns=["source", "NES", "pvalue"])
        out = pd.read_csv(out_path, sep="\t")
    if out.empty:
        return pd.DataFrame(columns=["source", "NES", "pvalue"])
    return pd.DataFrame({"source": out["pathway"].astype(str).values,
                         "NES": pd.to_numeric(out["NES"], errors="coerce").values,
                         "pvalue": pd.to_numeric(out["pval"], errors="coerce").values})


def run_tf_ulm(rank_series, collectri, min_targets=5):
    import decoupler as dc
    from scipy.stats import t as student_t
    mat = rank_series.to_frame().T
    mat.index = ["contrast"]
    df_t = max(len(rank_series) - 2, 1)
    if hasattr(dc, "mt") and hasattr(dc.mt, "ulm"):
        out = dc.mt.ulm(data=mat, net=collectri, tmin=min_targets)
        if isinstance(out, tuple):
            est, pval = out[0], out[1]
            res = pd.DataFrame({"source": est.columns, "activity_score": est.iloc[0].values,
                                "pvalue_dc": pval.iloc[0].values})
        else:
            try:
                est = dc.pp.get_obsm(out, key="score_ulm")
                pv = dc.pp.get_obsm(out, key="padj_ulm")
                res = pd.DataFrame({"source": est.var_names, "activity_score": est.X[0],
                                    "pvalue_dc": pv.X[0] if pv is not None else np.nan})
            except Exception:
                res = pd.DataFrame(columns=["source", "activity_score", "pvalue_dc"])
    else:
        acts, pvals = dc.run_ulm(mat=mat, net=collectri, min_n=min_targets)
        res = pd.DataFrame({"source": acts.columns, "activity_score": acts.iloc[0].values,
                            "pvalue_dc": pvals.iloc[0].values})
    if not res.empty:
        t_abs = np.abs(res["activity_score"].astype(float).values)
        res["pvalue"] = 2.0 * student_t.sf(t_abs, df=df_t)
    else:
        res["pvalue"] = np.nan
    return res[["source", "activity_score", "pvalue", "pvalue_dc"]]


def add_fdr(gsea_df, collection_map):
    g = gsea_df.copy()
    g["collection"] = g["source"].map(collection_map).fillna("NA")
    p = g["pvalue"].fillna(1.0).values
    g["FDR_pooled"] = false_discovery_control(p) if np.isfinite(p).any() else np.nan
    g["FDR"] = np.nan
    for coll, idx in g.groupby("collection").groups.items():
        sub = g.loc[idx, "pvalue"].fillna(1.0).values
        if len(sub):
            g.loc[idx, "FDR"] = false_discovery_control(sub)
    return g


def bh(p_arr):
    p = np.asarray(p_arr, dtype=float)
    ok = ~np.isnan(p)
    out = np.full(len(p), np.nan)
    if ok.sum() > 0:
        out[ok] = false_discovery_control(p[ok], method="bh")
    return out


# ============================================================================
def load_de():
    cmap = yaml.safe_load(MAP_YAML.read_text())["placenta_compartments"]
    rrho_comps = cmap["rrho_compartments"]
    m2h = dict(pd.read_csv(ORTHO, sep="\t")[["mouse_symbol", "human_symbol"]].values)
    print("[h09g] mouse E18.5 Late-vs-Relaxed DE")
    mo = sc.read_h5ad(MOUSE_H5AD)
    mo = mo[mo.obs["age"] == "E18.5"].copy()
    mo = mo[~unassigned_mask(mo.obs, ["celltype_majority"])].copy()
    mo.obs["compartment"] = mo.obs["celltype_majority"].astype(str).map(cmap["mouse"])
    mo = mo[mo.obs["compartment"].isin(rrho_comps)].copy()
    if "counts" not in mo.layers:
        mo.layers["counts"] = mo.X.copy()
    mouse_de = pseudobulk_de(mo, "compartment", "donor_id", "group",
                             "Relaxed", "Late_Stress", "sex", "mouse")
    del mo
    print("[h09g] human obese-vs-lean DE")
    hu = sc.read_h5ad(HUMAN_H5AD)
    hu = hu[hu.obs["compartment"].isin(rrho_comps)].copy()
    if "counts" not in hu.layers:
        hu.layers["counts"] = hu.X.copy()
    human_de = pseudobulk_de(hu, "compartment", "donor_id", "condition",
                             "lean", "obese", "side", "human")
    del hu
    return mouse_de, human_de, m2h, rrho_comps


def main():
    argparse.ArgumentParser().parse_args()
    mouse_de, human_de, m2h, rrho_comps = load_de()

    net = load_genesets_tsv(HUMAN_MSIGDB, COLLECTIONS, MIN_GENES)
    collection_map = dict(net[["source", "collection"]].drop_duplicates().values)
    print(f"[h09g] {net['source'].nunique()} human gene sets; loading CollecTRI human")
    collectri = load_collectri_human()

    pw_rows, tf_rows = [], []
    for comp in rrho_comps:
        if comp not in mouse_de or comp not in human_de:
            continue
        m_stat = mouse_de[comp].rename(index=m2h)
        m_stat = m_stat[~m_stat.index.duplicated()].dropna()
        h_stat = human_de[comp]
        print(f"  [{comp}] GSEA mouse+human")
        gm = add_fdr(run_gsea_on_ranks(m_stat, net, MIN_GENES), collection_map)
        gh = add_fdr(run_gsea_on_ranks(h_stat, net, MIN_GENES), collection_map)
        merged = gm.merge(gh, on="source", suffixes=("_mouse", "_human"))
        conc = merged[(merged["FDR_mouse"] < FDR_CUT) & (merged["FDR_human"] < FDR_CUT)
                      & (np.sign(merged["NES_mouse"]) == np.sign(merged["NES_human"]))].copy()
        conc["compartment"] = comp
        conc["direction"] = np.where(conc["NES_mouse"] > 0, "up_both", "down_both")
        pw_rows.append(conc)
        print(f"    {len(conc)} concordant pathways (FDR<{FDR_CUT} both, same sign)")

        print(f"  [{comp}] TF ULM mouse+human")
        tm = run_tf_ulm(m_stat, collectri); tm["FDR"] = bh(tm["pvalue"])
        th = run_tf_ulm(h_stat, collectri); th["FDR"] = bh(th["pvalue"])
        tmh = tm.merge(th, on="source", suffixes=("_mouse", "_human"))
        tconc = tmh[(tmh["FDR_mouse"] < FDR_CUT) & (tmh["FDR_human"] < FDR_CUT)
                    & (np.sign(tmh["activity_score_mouse"]) == np.sign(tmh["activity_score_human"]))].copy()
        tconc["compartment"] = comp
        tconc["direction"] = np.where(tconc["activity_score_mouse"] > 0, "up_both", "down_both")
        tf_rows.append(tconc)
        print(f"    {len(tconc)} concordant TFs")

    tab = OUT_DIR / "tables"; tab.mkdir(parents=True, exist_ok=True)
    pw = pd.concat(pw_rows, ignore_index=True) if pw_rows else pd.DataFrame()
    tf = pd.concat(tf_rows, ignore_index=True) if tf_rows else pd.DataFrame()
    pw_cols = ["compartment", "source", "collection_mouse", "direction",
               "NES_mouse", "FDR_mouse", "NES_human", "FDR_human"]
    pw = pw[[c for c in pw_cols if c in pw.columns]] if not pw.empty else pw
    pw.to_csv(tab / "h09g_concordant_pathways.csv", index=False)
    tf.to_csv(tab / "h09g_concordant_tfs.csv", index=False)

    print(f"\n[h09g] concordant pathways -> {tab/'h09g_concordant_pathways.csv'}")
    if not pw.empty:
        for comp in rrho_comps:
            sub = pw[pw["compartment"] == comp]
            if len(sub):
                print(f"\n  {comp} ({len(sub)} pathways), top by mouse FDR:")
                print(sub.nsmallest(8, "FDR_mouse")[["source", "direction", "NES_mouse", "NES_human"]]
                      .to_string(index=False))
    print(f"\n[h09g] concordant TFs -> {tab/'h09g_concordant_tfs.csv'}")
    if not tf.empty:
        print(tf.groupby("compartment").size().to_string())


if __name__ == "__main__":
    main()
