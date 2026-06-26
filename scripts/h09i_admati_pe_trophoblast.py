#!/usr/bin/env python
"""h09i_admati_pe_trophoblast.py -- second human placenta validation (PE, snRNA trophoblast).

Admati 2023 eoPE snRNA trophoblast (modality-matched to mouse + Gunter-Rahman). n=3 PE vs
2 control donors -> UNDERPOWERED for discovery, so the design is dual:
  (a) full pseudobulk DE eoPE-vs-control, flagged underpowered_exploratory (for the record),
  (b) HEADLINE = targeted confirmation: do the conserved trophoblast hypoxia leading-edge genes
      (mouse <-> Gunter-Rahman shared set, from h09h) move UP in eoPE trophoblast?
      -> sign test + targeted GSEA (HALLMARK_HYPOXIA) + scatter vs the conserved set.

Input file is TRANSPOSED: 22 metadata rows (cellID..donor_age) then gene rows; cells in columns.
Gene IDs are human symbols already (no Ensembl mapping).

Usage (from project root):
  uv run python scripts/h09i_admati_pe_trophoblast.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).parent))
from h09e_cross_species_rrho import OUT_DIR  # tables/plots root for the placenta arm
from h09g_pathways_tf import run_gsea_on_ranks, add_fdr, load_genesets_tsv, HUMAN_MSIGDB  # noqa

ADMATI = Path("data/human_validation/placenta/admati_2023_PE/"
              "sn_PE_TB_allcells_with_metadata_30-May-2023.txt")
META_ROWS = ["cellID", "celltype", "sample", "donorID", "total_molecules",
             "early_control", "early_PE", "female_fetus", "IUGR", "C-section_birth",
             "vaginal_birth", "induction", "non-induction", "magnesium",
             "spinal_anaesthesia", "epidural_anaesthesia", "general_anaesthesia",
             "delivery_week", "weight", "wieght_percentile-Dolberg", "donor_age"]
TAB = OUT_DIR / "tables"
PLOT = OUT_DIR / "plots" / "h09i_admati_pe"
MIN_CELLS = 10


def load_admati():
    """Parse the transposed file -> (counts cells x genes DataFrame, obs DataFrame)."""
    print(f"[h09i] reading {ADMATI}")
    n_meta = len(META_ROWS)
    # metadata rows: read first n_meta rows, transpose to per-cell obs
    meta = pd.read_csv(ADMATI, sep="\t", nrows=n_meta, header=None, index_col=0)
    meta.index = meta.index.str.strip()
    obs = meta.T
    obs.columns = [c.strip() for c in obs.columns]
    obs = obs.set_index("cellID")
    # gene rows: skip the metadata rows; first col = gene symbol, rest = per-cell counts
    expr = pd.read_csv(ADMATI, sep="\t", skiprows=range(1, n_meta + 1),
                       header=0, index_col=0, low_memory=False)
    expr.index = expr.index.str.strip()           # gene symbols
    # align cells
    expr = expr[obs.index]
    counts = sp.csr_matrix(expr.T.values.astype(np.float32))  # cells x genes
    print(f"  {counts.shape[0]} cells x {counts.shape[1]} genes; "
          f"donors: {obs['donorID'].nunique()}, celltypes: {obs['celltype'].nunique()}")
    return counts, expr.index.to_numpy(), obs


def pseudobulk_by_donor(counts, genes, obs):
    """Sum trophoblast counts per donor -> (count matrix donors x genes, donor meta)."""
    donors = obs["donorID"].astype(str).values
    rows, meta = [], []
    for d in sorted(set(donors)):
        m = donors == d
        if m.sum() < MIN_CELLS:
            continue
        rows.append(np.asarray(counts[m].sum(axis=0)).ravel())
        o = obs[obs["donorID"].astype(str) == d].iloc[0]
        # early_PE is a 0/1 indicator row
        cond = "PE" if str(o["early_PE"]).strip() in ("1", "1.0") else "control"
        meta.append({"donorID": d, "condition": cond,
                     "female_fetus": str(o.get("female_fetus", "")).strip()})
    cmat = pd.DataFrame(np.array(rows), index=[m["donorID"] for m in meta],
                        columns=genes).astype(int)
    md = pd.DataFrame(meta).set_index("donorID")
    return cmat, md


def run_de(cmat, md):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    md = md.copy()
    md["condition"] = pd.Categorical(md["condition"], categories=["control", "PE"])
    dds = DeseqDataSet(counts=cmat, metadata=md, design_factors=["condition"], quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "PE", "control"], quiet=True)
    st.summary()
    return st.results_df


def main():
    PLOT.mkdir(parents=True, exist_ok=True)
    counts, genes, obs = load_admati()

    # --- (a) full DE, flagged underpowered ---
    cmat, md = pseudobulk_by_donor(counts, genes, obs)
    print(f"[h09i] pseudobulk donors: {dict(md['condition'].value_counts())} "
          f"(UNDERPOWERED -- flagged)")
    res = run_de(cmat, md)
    res = res.dropna(subset=["stat"]).copy()
    res["flag"] = "underpowered_exploratory_3v2"
    res.sort_values("stat", ascending=False).to_csv(TAB / "h09i_admati_pe_de.csv")
    print(f"[h09i] DE -> {TAB/'h09i_admati_pe_de.csv'} ({len(res)} genes)")

    pe_stat = res["stat"]

    # --- (b) HEADLINE: targeted confirmation of the conserved hypoxia signature ---
    # conserved trophoblast hypoxia leading-edge genes from h09h (mouse <-> Gunter-Rahman shared)
    le = pd.read_csv(TAB / "h09h_leading_edge.csv")
    conserved = le[(le["pathway"] == "HALLMARK_HYPOXIA") & (le["shared"])
                   & (le["compartment"] == "trophoblast")]["gene"].tolist()
    present = [g for g in conserved if g in pe_stat.index]
    print(f"\n[h09i] conserved hypoxia genes: {len(conserved)} ({len(present)} testable in Admati)")

    sub = pe_stat.loc[present]
    n_up = int((sub > 0).sum()); n_tot = len(sub)
    bt = binomtest(n_up, n_tot, 0.5, alternative="greater")
    print(f"[h09i] sign test: {n_up}/{n_tot} conserved hypoxia genes UP in eoPE "
          f"(binom p={bt.pvalue:.2e})")

    # targeted GSEA: HALLMARK_HYPOXIA enrichment in the full eoPE ranking
    net = load_genesets_tsv(HUMAN_MSIGDB, ["H"], 15)
    cmap = dict(net[["source", "collection"]].drop_duplicates().values)
    g = add_fdr(run_gsea_on_ranks(pe_stat, net, 15), cmap)
    hyp = g[g["source"] == "HALLMARK_HYPOXIA"]
    if not hyp.empty:
        h = hyp.iloc[0]
        print(f"[h09i] HALLMARK_HYPOXIA in eoPE: NES={h['NES']:.2f}, FDR={h['FDR']:.2e}")

    # save targeted readout + per-gene table
    pd.DataFrame({"gene": present, "eoPE_stat": sub.values,
                  "up_in_eoPE": (sub > 0).values}).sort_values("eoPE_stat", ascending=False) \
        .to_csv(TAB / "h09i_conserved_hypoxia_in_eoPE.csv", index=False)

    # plot: conserved hypoxia genes' eoPE stat (bar), colored by direction
    s = sub.sort_values()
    fig, ax = plt.subplots(figsize=(6, max(3, 0.22 * len(s))))
    ax.barh(range(len(s)), s.values,
            color=np.where(s.values > 0, "#c0392b", "#2471a3"))
    ax.set_yticks(range(len(s))); ax.set_yticklabels(s.index, fontsize=6)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("eoPE-vs-control Wald stat (Admati trophoblast)")
    ax.set_title(f"Conserved hypoxia genes in eoPE: {n_up}/{n_tot} up "
                 f"(binom p={bt.pvalue:.1e})\n[underpowered 3v2 -- confirmatory]", fontsize=9)
    fig.tight_layout()
    fig.savefig(PLOT / "conserved_hypoxia_in_eoPE.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"\n[h09i] targeted readout -> {TAB/'h09i_conserved_hypoxia_in_eoPE.csv'}")
    print(f"[h09i] plot -> {PLOT/'conserved_hypoxia_in_eoPE.png'}")
    print("\nheadline: does a 3rd independent human stressor (eoPE) recover the "
          "conserved trophoblast hypoxia program? -> sign test + GSEA above")


if __name__ == "__main__":
    main()
