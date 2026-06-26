#!/usr/bin/env python
"""h09k_diagnostics.py -- interrogate the 2x2 trophoblast anomalies before interpreting.

Two anomalies in h09k:
  (1) E18.5 x loPE x trophoblast = DISCORDANT (peak 24.7) while E12.5 x loPE x trophoblast
      = concordant_up (peak 37) -- a sign flip on the MOUSE axis.
  (2) loPE concords far more strongly than eoPE everywhere -- possibly a power artifact
      (eoPE = 10 PE vs only 3 controls; loPE = 7 vs 6).

Diagnostics (trophoblast only):
  A. Spearman among the 4 rankings -- esp. mouse E12.5 vs E18.5 (is the flip mechanical?).
  B. stat distributions per contrast (is eoPE shrunken vs loPE?).
  C. CAUSAL: subsample loPE controls 6->3 (match eoPE), recompute, re-RRHO x5 seeds
     -- does loPE's peak collapse toward eoPE's?

Usage:  uv run python scripts/h09k_diagnostics.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from _utils import unassigned_mask  # noqa: E402
from h09e_cross_species_rrho import (  # noqa: E402
    pseudobulk_de, rrho_matrix, classify_rrho_concordance, MOUSE_H5AD, ORTHO, MAP_YAML, OUT_DIR,
)
from h09k_admati_2x2 import concordance_peak, PB_PATH, META_PATH, HUMAN_ARMS  # noqa: E402

COMP = "trophoblast"


def human_troph_de(pb, meta, ref, test, control_sub_n=None, seed=0):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    sel = meta[(meta["compartment"] == COMP) & (meta["condition"].isin([ref, test]))].copy()
    if control_sub_n is not None:
        ctrl = sel[sel["condition"] == ref]
        keep = ctrl.sample(n=control_sub_n, random_state=seed).index
        sel = sel[(sel["condition"] == test) | (sel.index.isin(keep))]
    cmat = pb.loc[sel.index].astype(int)
    cmat = cmat.loc[:, cmat.sum(axis=0) > 0]
    md = sel[["condition"]].copy()
    md["condition"] = pd.Categorical(md["condition"], categories=[ref, test])
    dds = DeseqDataSet(counts=cmat, metadata=md, design_factors=["condition"], quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", test, ref], quiet=True)
    st.summary()
    return st.results_df.dropna(subset=["stat"])


def main():
    m2h = dict(pd.read_csv(ORTHO, sep="\t")[["mouse_symbol", "human_symbol"]].values)
    cmap = yaml.safe_load(MAP_YAML.read_text())["placenta_compartments"]["mouse"]

    # --- mouse trophoblast rankings (E12.5 Early-vs-Relaxed, E18.5 Late-vs-Relaxed) ---
    mo = sc.read_h5ad(MOUSE_H5AD)
    mo = mo[~unassigned_mask(mo.obs, ["celltype_majority"])].copy()
    mo.obs["compartment"] = mo.obs["celltype_majority"].astype(str).map(cmap)
    if "counts" not in mo.layers:
        mo.layers["counts"] = mo.X.copy()
    def mouse_rank(age, test):
        sub = mo[(mo.obs["age"] == age) & (mo.obs["compartment"] == COMP)].copy()
        de = pseudobulk_de(sub, "compartment", "donor_id", "group", "Relaxed", test, "sex", f"m/{age}")
        s = de[COMP].rename(index=m2h)
        return s[~s.index.duplicated()].dropna()
    m_e125 = mouse_rank("E12.5", "Early_Stress")
    m_e185 = mouse_rank("E18.5", "Late_Stress")
    del mo

    # --- human trophoblast rankings ---
    pb = pd.read_parquet(PB_PATH); meta = pd.read_csv(META_PATH, index_col=0)
    eo = human_troph_de(pb, meta, **{k: HUMAN_ARMS["eoPE"][k] for k in ("ref", "test")})["stat"]
    lo = human_troph_de(pb, meta, **{k: HUMAN_ARMS["loPE"][k] for k in ("ref", "test")})["stat"]

    ranks = {"mouse_E12.5": m_e125, "mouse_E18.5": m_e185, "human_eoPE": eo, "human_loPE": lo}

    # === A. correlation matrix among the 4 rankings ===
    print("\n=== A. Spearman among trophoblast rankings (shared genes) ===")
    names = list(ranks)
    cor = pd.DataFrame(index=names, columns=names, dtype=float)
    for a in names:
        for b in names:
            common = ranks[a].index.intersection(ranks[b].index)
            cor.loc[a, b] = spearmanr(ranks[a].loc[common], ranks[b].loc[common])[0]
    print(cor.round(3).to_string())
    print(f"\n  >> mouse E12.5 vs E18.5 trophoblast: r = {cor.loc['mouse_E12.5','mouse_E18.5']:.3f}")
    print("     (negative => the discordant E18.5xloPE cell is MECHANICAL: the two mouse "
          "stress-timings disagree in sign, consistent with 8g age-specific effects)")

    # === B. stat distributions (power) ===
    print("\n=== B. stat distributions per contrast (is eoPE shrunken vs loPE?) ===")
    for nm, s in ranks.items():
        print(f"  {nm:14s} n={len(s):6d}  std={s.std():.2f}  max|stat|={s.abs().max():6.2f}  "
              f"q05={s.quantile(.05):6.2f}  q95={s.quantile(.95):6.2f}")
    # n significant if padj available — recompute quickly for human (approx via |stat|>3)
    for nm, s in [("human_eoPE", eo), ("human_loPE", lo)]:
        print(f"  {nm}: |stat|>3 genes = {(s.abs()>3).sum()}  (proxy for sig; eoPE has 3 controls)")

    # === C. CAUSAL subsample: loPE controls 6->3, recompute, re-RRHO vs both mouse ages ===
    print("\n=== C. loPE controls downsampled 6->3 (match eoPE), RRHO peak vs mouse ===")
    full_peaks = {}
    for age, mr in [("E12.5", m_e125), ("E18.5", m_e185)]:
        mat, cut = rrho_matrix(mr, lo)
        full_peaks[age] = concordance_peak(mat)
    print(f"  FULL loPE (6 ctrl): vs E12.5 peak={full_peaks['E12.5']:.1f}, "
          f"vs E18.5 peak={full_peaks['E18.5']:.1f}")
    print(f"  (reference: eoPE 3 ctrl gave E12.5 peak~4.8, E18.5 peak~7.2)")
    for age, mr in [("E12.5", m_e125), ("E18.5", m_e185)]:
        peaks = []
        for seed in range(5):
            lo_sub = human_troph_de(pb, meta, HUMAN_ARMS["loPE"]["ref"],
                                    HUMAN_ARMS["loPE"]["test"], control_sub_n=3, seed=seed)["stat"]
            mat, cut = rrho_matrix(mr, lo_sub)
            peaks.append(concordance_peak(mat))
        peaks = np.array(peaks)
        print(f"  loPE@3ctrl vs {age}: peak = {peaks.mean():.1f} ± {peaks.std():.1f} "
              f"(seeds: {np.round(peaks,1).tolist()})  [full was {full_peaks[age]:.1f}]")

    print("\n  INTERPRET: if loPE@3ctrl peak collapses toward eoPE's (~5-7) => the loPE>eoPE "
          "gap is a POWER artifact (control n). If it stays high => biology.")


if __name__ == "__main__":
    main()
