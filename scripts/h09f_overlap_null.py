#!/usr/bin/env python
"""h09f_overlap_null.py -- overlap gene lists + permutation null for the placenta RRHO.

Builds on h09e (re-runs the same compartment DE, no CSV round-trip). For each
compartment:
  1. extracts the concordant-UP leading-edge genes (top-up in BOTH species at the
     up-quadrant peak) -> the genes that name the cross-species signal for Fig 4.
  2. gene-label-shuffle null: permute human gene<->stat assignment N times, recompute
     the up-quadrant RRHO peak each time -> empirical p for the observed peak.

Reuses pseudobulk_de + rrho_matrix imported VERBATIM from h09e (same DE, same RRHO).

Usage (from project root):
  uv run python scripts/h09f_overlap_null.py --n-perm 10000 --n-jobs 24
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from _utils import unassigned_mask, parallel_map  # noqa: E402
from h09e_cross_species_rrho import (  # noqa: E402
    pseudobulk_de, rrho_matrix, MOUSE_H5AD, HUMAN_H5AD, ORTHO, MAP_YAML, OUT_DIR,
)


def up_peak(mat):
    h = mat.shape[0] // 2
    return float(mat[:h, :h].max())


def leading_edge_up(m_stat, h_stat):
    """Concordant-UP genes at the up-quadrant peak: top-up in both species."""
    mat, cutoffs = rrho_matrix(m_stat, h_stat)
    if mat is None:
        return None, pd.DataFrame()
    common = m_stat.index.intersection(h_stat.index)
    a = m_stat.loc[common]; b = h_stat.loc[common]
    ra = a.rank(ascending=False, method="first")
    rb = b.rank(ascending=False, method="first")
    h = mat.shape[0] // 2
    sub = mat[:h, :h]
    i, j = np.unravel_index(int(np.argmax(sub)), sub.shape)
    ci, cj = cutoffs[i], cutoffs[j]
    sel = common[(ra.values <= ci) & (rb.values <= cj)]
    df = pd.DataFrame({"gene": sel,
                       "mouse_stat": a.loc[sel].values,
                       "human_stat": b.loc[sel].values})
    df["mean_rank"] = (ra.loc[sel].values + rb.loc[sel].values) / 2
    return up_peak(mat), df.sort_values("mean_rank").reset_index(drop=True)


def _shuffle_chunk(args):
    """Top-level (picklable) worker: n_sh shuffles of human labels -> up-quadrant peaks."""
    m_vals, h_vals, n_sh, seed = args
    rng = np.random.default_rng(seed)
    n = len(m_vals)
    idx = pd.RangeIndex(n)
    m_ser = pd.Series(m_vals, index=idx)
    peaks = np.empty(n_sh)
    for s in range(n_sh):
        h_ser = pd.Series(rng.permutation(h_vals), index=idx)
        mat, _ = rrho_matrix(m_ser, h_ser)
        peaks[s] = up_peak(mat)
    return peaks


def permutation_null(m_stat, h_stat, obs_peak, n_perm, n_jobs):
    common = m_stat.index.intersection(h_stat.index)
    m_vals = m_stat.loc[common].to_numpy()
    h_vals = h_stat.loc[common].to_numpy()
    n_chunks = max(n_jobs, 1) * 2
    base = n_perm // n_chunks
    sizes = [base + (1 if k < n_perm - base * n_chunks else 0) for k in range(n_chunks)]
    jobs = [(m_vals, h_vals, sz, 1000 + k) for k, sz in enumerate(sizes) if sz > 0]
    peaks = []
    for _job, res, err in parallel_map(_shuffle_chunk, jobs, n_jobs=n_jobs,
                                       use_threads=False, desc="shuffle"):
        if err:
            print(f"    [warn] shuffle chunk failed: {err.splitlines()[-1]}")
            continue
        peaks.append(res)
    null = np.concatenate(peaks) if peaks else np.array([])
    n_ge = int((null >= obs_peak).sum())
    emp_p = (1 + n_ge) / (1 + len(null)) if len(null) else np.nan
    return emp_p, len(null), null


def load_de():
    """Re-run mouse + human compartment DE (mirrors h09e)."""
    cmap = yaml.safe_load(MAP_YAML.read_text())["placenta_compartments"]
    rrho_comps = cmap["rrho_compartments"]
    ortho = pd.read_csv(ORTHO, sep="\t")
    m2h = dict(zip(ortho["mouse_symbol"], ortho["human_symbol"]))

    print("[h09f] mouse E18.5 Late-vs-Relaxed DE")
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

    print("[h09f] human obese-vs-lean DE")
    hu = sc.read_h5ad(HUMAN_H5AD)
    hu = hu[hu.obs["compartment"].isin(rrho_comps)].copy()
    if "counts" not in hu.layers:
        hu.layers["counts"] = hu.X.copy()
    human_de = pseudobulk_de(hu, "compartment", "donor_id", "condition",
                             "lean", "obese", "side", "human")
    del hu
    return mouse_de, human_de, m2h, rrho_comps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--n-jobs", type=int, default=24)
    args = ap.parse_args()

    mouse_de, human_de, m2h, rrho_comps = load_de()
    tab_dir = OUT_DIR / "tables"; tab_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for comp in rrho_comps:
        if comp not in mouse_de or comp not in human_de:
            continue
        m_stat = mouse_de[comp].rename(index=m2h)
        m_stat = m_stat[~m_stat.index.duplicated()].dropna()
        h_stat = human_de[comp]
        obs_peak, le = leading_edge_up(m_stat, h_stat)
        if obs_peak is None:
            print(f"  [{comp}] <200 shared genes -- skipped")
            continue
        le.to_csv(tab_dir / f"h09f_overlap_{comp}.csv", index=False)
        print(f"  [{comp}] obs up-peak={obs_peak:.2f}, {len(le)} leading-edge genes; "
              f"running {args.n_perm} shuffles")
        emp_p, n_null, _ = permutation_null(m_stat, h_stat, obs_peak,
                                            args.n_perm, args.n_jobs)
        rows.append({"compartment": comp, "obs_up_peak_neg_log10_p": round(obs_peak, 2),
                     "n_leading_edge_genes": len(le), "empirical_p": emp_p,
                     "n_perm": n_null})
        print(f"  [{comp}] empirical p = {emp_p:.2e}  (top genes: "
              f"{', '.join(le['gene'].head(8))})")

    summ = pd.DataFrame(rows)
    summ.to_csv(tab_dir / "h09f_permutation_null.csv", index=False)
    print(f"\n[h09f] null summary -> {tab_dir/'h09f_permutation_null.csv'}")
    print(summ.to_string(index=False))
    print(f"[h09f] per-compartment overlap gene lists -> {tab_dir}/h09f_overlap_*.csv")


if __name__ == "__main__":
    main()
