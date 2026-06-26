#!/usr/bin/env python
"""h09e_cross_species_rrho.py -- Phase 9 placenta cross-species validation.

Compartment-level RRHO: mouse E18.5 Late-vs-Relaxed  <->  human term obese-vs-lean
(Gunter-Rahman GSE271976). Both sides pseudobulked per compartment (the bridge level;
mouse/human trophoblast subtypes have no 1:1 homology). Genes bridged via the
mouse<->human 1:1 ortholog table. Ranking metric = signed DESeq2 Wald `stat`
(identical to 8f's RRHO, so this arm matches the brain arm method exactly).

Steps:
  1. mouse DE: subset E18.5 {Late_Stress, Relaxed}, celltype_majority->compartment,
     pseudobulk per (donor x compartment), PyDESeq2 ~ sex + group, Late vs Relaxed.
  2. human DE: h09c compartments, pseudobulk per (sample x compartment),
     PyDESeq2 ~ side + condition, obese vs lean.
  3. orthology bridge (refs/mouse_human_orthologs.tsv) -> shared human symbols.
  4. RRHO per matched compartment (rrho_matrix/classify lifted verbatim from 08f).

Usage (from project root):
  uv run python scripts/h09e_cross_species_rrho.py --n-jobs 8
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
from scipy.stats import hypergeom, spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from _utils import unassigned_mask  # noqa: E402

MOUSE_H5AD = Path("results/placenta/h5ad/08_annotated/all_samples.h5ad")
HUMAN_H5AD = Path("data/human_validation/placenta/gunter_rahman_2025_GSE271976/h5ad/h09c_integrated.h5ad")
ORTHO = Path("refs/mouse_human_orthologs.tsv")
MAP_YAML = Path("config/cross_species_celltype_map.yaml")
OUT_DIR = Path("data/human_validation/placenta/gunter_rahman_2025_GSE271976")
MIN_CELLS, MIN_DONORS = 10, 2   # pseudobulk inclusion (matches 8b spirit)


# ============================================================================
# RRHO -- lifted VERBATIM from 08f_cross_tissue.py (do not diverge; keeps the
# placenta and brain cross-species arms methodologically identical).
# ============================================================================
def rrho_matrix(stats_a, stats_b, step=100):
    common = stats_a.index.intersection(stats_b.index)
    if len(common) < 200:
        return None, None
    a = stats_a.loc[common]; b = stats_b.loc[common]
    rank_a = a.rank(ascending=False, method="first")
    rank_b = b.rank(ascending=False, method="first")
    n = len(common)
    cutoffs = np.arange(step, n, step)
    if len(cutoffs) < 3:
        return None, None
    if len(cutoffs) > 40:
        idx = np.linspace(0, len(cutoffs) - 1, 40).astype(int)
        cutoffs = cutoffs[idx]
    ra = rank_a.to_numpy(); rb = rank_b.to_numpy()
    cut = cutoffs[:, None]
    a_ind = (ra[None, :] <= cut).astype(np.float64)
    b_ind = (rb[None, :] <= cut).astype(np.float64)
    overlap = a_ind @ b_ind.T
    na = a_ind.sum(axis=1); nb = b_ind.sum(axis=1)
    na_grid, nb_grid = np.meshgrid(na, nb, indexing="ij")
    p = hypergeom.sf(overlap - 1, n, na_grid, nb_grid)
    mat = -np.log10(np.maximum(p, 1e-300))
    mat[overlap == 0] = 0.0
    return mat, cutoffs


def classify_rrho_concordance(mat, cutoffs):
    if mat is None:
        return "none", 0.0
    k = mat.shape[0]; h = k // 2
    quad_tl = mat[:h, :h].max(); quad_br = mat[h:, h:].max()
    quad_tr = mat[:h, h:].max(); quad_bl = mat[h:, :h].max()
    quads = {"concordant_up": quad_tl, "concordant_down": quad_br,
             "discordant": max(quad_tr, quad_bl)}
    best = max(quads, key=quads.get)
    if quads[best] < 2:
        return "none", float(quads[best])
    return best, float(quads[best])


# ============================================================================
# Pseudobulk + DE
# ============================================================================
def pseudobulk_de(adata, compartment_col, donor_col, group_col, ref_level, test_level,
                  covariate_col, label):
    """Per-compartment pseudobulk -> PyDESeq2 (~ covariate + group). Returns
    {compartment: Series(stat indexed by gene)}. Mirrors 8b: animal=unit, raw counts."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    adata = adata[adata.obs[group_col].isin([ref_level, test_level])].copy()
    out = {}
    for comp in sorted(adata.obs[compartment_col].dropna().unique()):
        sub = adata[adata.obs[compartment_col] == comp]
        donors = sub.obs[donor_col].unique()
        # build pseudobulk count matrix donors x genes
        rows, meta = [], []
        for d in donors:
            cells = sub[sub.obs[donor_col] == d]
            if cells.n_obs < MIN_CELLS:
                continue
            counts = np.asarray(cells.layers["counts"].sum(axis=0)).ravel() \
                if "counts" in cells.layers else np.asarray(cells.X.sum(axis=0)).ravel()
            rows.append(counts)
            o = cells.obs.iloc[0]
            meta.append({donor_col: d, group_col: o[group_col], covariate_col: o[covariate_col]})
        if len(rows) < 2 * MIN_DONORS:
            print(f"    [{label}/{comp}] only {len(rows)} donors -- skipped")
            continue
        md = pd.DataFrame(meta).set_index(donor_col)
        # need both groups with >=MIN_DONORS
        gc = md[group_col].value_counts()
        if gc.get(ref_level, 0) < MIN_DONORS or gc.get(test_level, 0) < MIN_DONORS:
            print(f"    [{label}/{comp}] groups too small {dict(gc)} -- skipped")
            continue
        cmat = pd.DataFrame(np.array(rows), index=md.index, columns=adata.var_names).astype(int)
        # drop covariate if constant or aliased (8b rule)
        design_factors = [group_col]
        if md[covariate_col].nunique() > 1:
            design_factors = [covariate_col, group_col]
        md[group_col] = pd.Categorical(md[group_col], categories=[ref_level, test_level])
        dds = DeseqDataSet(counts=cmat, metadata=md, design_factors=design_factors, quiet=True)
        dds.deseq2()
        st = DeseqStats(dds, contrast=[group_col, test_level, ref_level], quiet=True)
        st.summary()
        res = st.results_df.dropna(subset=["stat"])
        out[comp] = res["stat"]
        print(f"    [{label}/{comp}] {len(rows)} donors ({dict(gc)}), {len(res)} genes")
    return out


def _plot_rrho(mat, cutoffs, comp, klass, peak, rho, out_path):
    fig, ax = plt.subplots(figsize=(5, 4.3))
    im = ax.imshow(mat, origin="lower", aspect="auto", cmap="inferno")
    ax.set_xlabel("human obese-vs-lean rank cutoff")
    ax.set_ylabel("mouse Late-vs-Relaxed rank cutoff")
    ax.set_title(f"{comp}: {klass}\npeak -log10p={peak:.1f}  Spearman r={rho:.2f}", fontsize=9)
    fig.colorbar(im, ax=ax, label="-log10 hypergeom p")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-jobs", type=int, default=8)  # reserved; DESeq2 threads internally
    args = ap.parse_args()

    cmap = yaml.safe_load(MAP_YAML.read_text())["placenta_compartments"]
    mouse_map, human_map = cmap["mouse"], cmap["human"]
    rrho_comps = cmap["rrho_compartments"]
    ortho = pd.read_csv(ORTHO, sep="\t")
    m2h = dict(zip(ortho["mouse_symbol"], ortho["human_symbol"]))
    print(f"[h09e] {len(m2h)} mouse->human 1:1 orthologs")

    # --- mouse DE ---
    print("[h09e] mouse E18.5 Late-vs-Relaxed pseudobulk DE")
    mo = sc.read_h5ad(MOUSE_H5AD)
    mo = mo[mo.obs["age"] == "E18.5"].copy()
    mo = mo[~unassigned_mask(mo.obs, ["celltype_majority"])].copy()
    mo.obs["compartment"] = mo.obs["celltype_majority"].astype(str).map(mouse_map)
    mo = mo[mo.obs["compartment"].isin(rrho_comps)].copy()   # drop other_extraembryonic/NaN
    if "counts" not in mo.layers:
        mo.layers["counts"] = mo.X.copy()
    mouse_de = pseudobulk_de(mo, "compartment", "donor_id", "group",
                             ref_level="Relaxed", test_level="Late_Stress",
                             covariate_col="sex", label="mouse")
    del mo

    # --- human DE ---
    print("[h09e] human obese-vs-lean pseudobulk DE")
    hu = sc.read_h5ad(HUMAN_H5AD)
    hu = hu[hu.obs["compartment"].isin(rrho_comps)].copy()
    if "counts" not in hu.layers:
        hu.layers["counts"] = hu.X.copy()
    human_de = pseudobulk_de(hu, "compartment", "donor_id", "condition",
                             ref_level="lean", test_level="obese",
                             covariate_col="side", label="human")
    del hu

    # --- bridge mouse stats to human symbols, RRHO per compartment ---
    plot_dir = OUT_DIR / "plots" / "h09e"
    plot_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for comp in rrho_comps:
        if comp not in mouse_de or comp not in human_de:
            print(f"  [{comp}] missing in one species -- skipped")
            continue
        m_stat = mouse_de[comp].rename(index=m2h)        # mouse symbol -> human symbol
        m_stat = m_stat[~m_stat.index.duplicated()].dropna()
        h_stat = human_de[comp]
        mat, cutoffs = rrho_matrix(m_stat, h_stat)
        if mat is None:
            print(f"  [{comp}] <200 shared genes -- skipped")
            continue
        klass, peak = classify_rrho_concordance(mat, cutoffs)
        common = m_stat.index.intersection(h_stat.index)
        rho, rho_p = spearmanr(m_stat.loc[common], h_stat.loc[common])
        _plot_rrho(mat, cutoffs, comp, klass, peak, rho, plot_dir / f"rrho_{comp}.png")
        rows.append({"compartment": comp, "n_shared_genes": len(common),
                     "rrho_class": klass, "peak_neg_log10_p": round(peak, 2),
                     "spearman_r": round(rho, 3), "spearman_p": rho_p})
        print(f"  [{comp}] {klass}  peak={peak:.1f}  r={rho:.3f}  (n={len(common)})")

    summ = pd.DataFrame(rows)
    tab_dir = OUT_DIR / "tables"; tab_dir.mkdir(parents=True, exist_ok=True)
    summ.to_csv(tab_dir / "h09e_rrho_summary.csv", index=False)
    # also save the per-compartment ranked stats for overlap-gene follow-up
    pd.concat({c: mouse_de[c].rename(index=m2h) for c in mouse_de}, names=["compartment"]) \
        .rename("mouse_stat").to_csv(tab_dir / "h09e_mouse_stat_humansym.csv")
    pd.concat({c: human_de[c] for c in human_de}, names=["compartment"]) \
        .rename("human_stat").to_csv(tab_dir / "h09e_human_stat.csv")

    print(f"\n[h09e] RRHO summary -> {tab_dir/'h09e_rrho_summary.csv'}")
    print(summ.to_string(index=False))
    print(f"[h09e] plots -> {plot_dir}")
    print("\nheadline = trophoblast row (cleanest cross-species homology)")


if __name__ == "__main__":
    main()
