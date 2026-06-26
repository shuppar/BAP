#!/usr/bin/env python
"""h09h_leading_edge.py -- leading-edge genes for every concordant cross-species pathway.

For each pathway in h09g_concordant_pathways.csv, computes the leading-edge genes in
BOTH species (compute_leading_edge lifted VERBATIM from 08c) and reports each species'
leading edge separately + the overlap (the cross-species core driving that pathway).

Re-runs the same compartment DE (h09e pseudobulk_de); pulls pathway membership from the
human msigdb TSV. No fgsea re-run -- leading edge is deterministic from ranks + members.

Output:
  h09h_leading_edge.csv  -- one row per (compartment, pathway, gene) with columns:
    in_mouse_le, in_human_le, shared, mouse_stat, human_stat, mouse_rank_in_le, human_rank_in_le
  plus per-pathway overlap summary printed + h09h_le_overlap_summary.csv

Usage (from project root):
  uv run python scripts/h09h_leading_edge.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from _utils import unassigned_mask  # noqa: E402
from h09e_cross_species_rrho import (  # noqa: E402
    pseudobulk_de, MOUSE_H5AD, HUMAN_H5AD, ORTHO, MAP_YAML, OUT_DIR,
)

HUMAN_MSIGDB = Path("refs/msigdb_human.tsv")
TAB = OUT_DIR / "tables"


def compute_leading_edge(rank_series, members, nes):
    """VERBATIM from 08c_pathways.py. Members driving enrichment, sorted by |stat| desc."""
    present = [g for g in members if g in rank_series.index]
    if not present:
        return []
    sub = rank_series.loc[present]
    if nes is not None and not np.isnan(nes):
        sub = sub[sub > 0] if nes > 0 else sub[sub < 0]
    sub = sub.reindex(sub.abs().sort_values(ascending=False).index)
    return list(zip(sub.index.tolist(), sub.values.tolist()))


def load_de():
    cmap = yaml.safe_load(MAP_YAML.read_text())["placenta_compartments"]
    rrho_comps = cmap["rrho_compartments"]
    m2h = dict(pd.read_csv(ORTHO, sep="\t")[["mouse_symbol", "human_symbol"]].values)
    print("[h09h] mouse E18.5 Late-vs-Relaxed DE")
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
    print("[h09h] human obese-vs-lean DE")
    hu = sc.read_h5ad(HUMAN_H5AD)
    hu = hu[hu.obs["compartment"].isin(rrho_comps)].copy()
    if "counts" not in hu.layers:
        hu.layers["counts"] = hu.X.copy()
    human_de = pseudobulk_de(hu, "compartment", "donor_id", "condition",
                             "lean", "obese", "side", "human")
    del hu
    return mouse_de, human_de, m2h


def main():
    conc_path = TAB / "h09g_concordant_pathways.csv"
    if not conc_path.is_file():
        sys.exit(f"missing {conc_path} (run h09g first)")
    conc = pd.read_csv(conc_path)
    if conc.empty:
        sys.exit("no concordant pathways to process")

    # pathway -> member gene set (human symbols)
    gs = pd.read_csv(HUMAN_MSIGDB, sep="\t")
    members = gs.groupby("gs_name")["gene_symbol"].apply(set).to_dict()

    mouse_de, human_de, m2h = load_de()

    rows, summ = [], []
    for _, r in conc.iterrows():
        comp, pw = r["compartment"], r["source"]
        if comp not in mouse_de or comp not in human_de or pw not in members:
            continue
        m_stat = mouse_de[comp].rename(index=m2h)
        m_stat = m_stat[~m_stat.index.duplicated()].dropna()
        h_stat = human_de[comp]
        mem = members[pw]
        m_le = compute_leading_edge(m_stat, mem, r["NES_mouse"])
        h_le = compute_leading_edge(h_stat, mem, r["NES_human"])
        m_rank = {g: i + 1 for i, (g, _) in enumerate(m_le)}
        h_rank = {g: i + 1 for i, (g, _) in enumerate(h_le)}
        m_genes, h_genes = set(m_rank), set(h_rank)
        shared = m_genes & h_genes
        for g in sorted(m_genes | h_genes):
            rows.append({
                "compartment": comp, "pathway": pw, "direction": r["direction"], "gene": g,
                "in_mouse_le": g in m_genes, "in_human_le": g in h_genes,
                "shared": g in shared,
                "mouse_stat": round(float(m_stat.get(g, np.nan)), 3) if g in m_stat.index else np.nan,
                "human_stat": round(float(h_stat.get(g, np.nan)), 3) if g in h_stat.index else np.nan,
                "mouse_rank_in_le": m_rank.get(g, np.nan),
                "human_rank_in_le": h_rank.get(g, np.nan),
            })
        jacc = len(shared) / len(m_genes | h_genes) if (m_genes | h_genes) else 0.0
        summ.append({"compartment": comp, "pathway": pw, "direction": r["direction"],
                     "n_mouse_le": len(m_genes), "n_human_le": len(h_genes),
                     "n_shared": len(shared), "jaccard": round(jacc, 3),
                     "shared_genes": ",".join(sorted(shared))})

    le_df = pd.DataFrame(rows)
    sm_df = pd.DataFrame(summ).sort_values(["compartment", "n_shared"], ascending=[True, False])
    TAB.mkdir(parents=True, exist_ok=True)
    le_df.to_csv(TAB / "h09h_leading_edge.csv", index=False)
    sm_df.to_csv(TAB / "h09h_le_overlap_summary.csv", index=False)

    print(f"\n[h09h] leading-edge detail -> {TAB/'h09h_leading_edge.csv'} ({len(le_df)} rows)")
    print(f"[h09h] overlap summary -> {TAB/'h09h_le_overlap_summary.csv'}")
    for comp in sm_df["compartment"].unique():
        sub = sm_df[sm_df["compartment"] == comp]
        print(f"\n  {comp}: {len(sub)} pathways, mean Jaccard={sub['jaccard'].mean():.2f}")
        print(sub.head(6)[["pathway", "direction", "n_mouse_le", "n_human_le",
                           "n_shared", "jaccard"]].to_string(index=False))
    # spotlight hypoxia (the headline)
    hyp = sm_df[sm_df["pathway"] == "HALLMARK_HYPOXIA"]
    if not hyp.empty:
        print("\n[h09h] HALLMARK_HYPOXIA shared leading-edge genes:")
        for _, h in hyp.iterrows():
            print(f"  {h['compartment']}: {h['shared_genes']}")


if __name__ == "__main__":
    main()
