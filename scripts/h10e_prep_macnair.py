#!/usr/bin/env python
"""h10e_prep_macnair.py -- pseudobulk the Macnair 2025 MS lesion snRNA-seq (discovery cohort).

Brain cross-species ARM-B (MS = STRESSED-GLIA REFERENCE, NOT etiology). Reported separately
from the psychiatric/neurodevelopmental datasets; the valid claim is "shared stressed-glia
program," never causation.

Input (data/human_validation/brain/macnair_2025_MS/):
  ms_lesions_snRNAseq_cleaned_counts_matrix_2023-09-12.mtx.gz -- 33939 genes x 632375 cells (gene-major)
  ms_lesions_snRNAseq_row_data_2023-09-12.txt.gz -- gene_id,ensembl,symbol (row order = mtx rows)
  ms_lesions_snRNAseq_col_data_2023-09-12.txt.gz -- per-cell metadata (row order = mtx cols);
       key cols: individual_id_anon (donor), sample_id_anon, type_broad, matter (GM/WM),
       lesion_type, diagnosis (CTR/SPMS/PPMS/RRMS), sex, exclude_pseudobulk.

Honors the authors' `exclude_pseudobulk` flag (mixed-cluster cells dropped, matching their
analysis). Diagnosis collapsed MS = {SPMS,PPMS,RRMS} vs CTR. Bridge = broad 7-class;
B/T cells dropped. Unit = sample_id_anon (individual x lesion), mirroring the mouse sample
unit; `matter` (GM/WM) recorded as the region axis. genes = `symbol` column.

Outputs (.../tables/):
  h10e_macnair_pseudobulk_primary.parquet -- (sample x broad) x gene counts (int)
  h10e_macnair_group_meta_primary.csv     -- group -> sample, donor, broad, diagnosis, sex,
                                             matter, lesion_type, n_cells

Usage (WS, from project root):
  uv run python scripts/h10e_prep_macnair.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

DATA = Path("data/human_validation/brain/macnair_2025_MS")
MTX = DATA / "ms_lesions_snRNAseq_cleaned_counts_matrix_2023-09-12.mtx.gz"
ROW = DATA / "ms_lesions_snRNAseq_row_data_2023-09-12.txt.gz"
COL = DATA / "ms_lesions_snRNAseq_col_data_2023-09-12.txt.gz"
TAB = DATA / "tables"
MIN_DONORS = 3

BROAD_MAP = {"Excitatory neurons": "ExN", "Inhibitory neurons": "InN",
             "Oligodendrocytes": "Oli", "Astrocytes": "Ast", "OPCs + COPs": "OPC",
             "Microglia": "Mic", "Endo + Peri": "Endo"}   # B cells / T cells -> dropped


def main():
    print(f"[h10e] reading col_data {COL}")
    obs = pd.read_csv(COL, low_memory=False)
    obs.columns = [c.strip() for c in obs.columns]
    n_cells = len(obs)
    print(f"  {n_cells} cells; cols: {list(obs.columns)[:8]}...")

    # exclude_pseudobulk is the authors' mixed-cluster flag -> honor it
    excl = obs["exclude_pseudobulk"].astype(str).str.upper().isin(["TRUE", "1"])
    print(f"  exclude_pseudobulk=TRUE: {excl.sum()} cells dropped (authors' QC)")

    obs["broad"] = obs["type_broad"].map(BROAD_MAP)
    diag_collapse = {"SPMS": "MS", "PPMS": "MS", "RRMS": "MS", "CTR": "Control"}
    obs["diagnosis"] = obs["diagnosis"].map(diag_collapse)   # unmapped -> NaN
    keep = (~excl) & obs["broad"].notna() & obs["diagnosis"].notna()
    print(f"  cells kept: {keep.sum()}/{n_cells}")
    print("  broad census (kept):\n", obs.loc[keep, "broad"].value_counts().to_string())
    print("  diagnosis census (kept cells):\n",
          obs.loc[keep, "diagnosis"].value_counts().to_string())
    print("  matter census (kept):\n", obs.loc[keep, "matter"].value_counts().to_string())

    # unit = donor (individual_id_anon). Mouse is whole-brain (no GM/WM analog), so GM/WM has
    # nothing to RRHO against -> donor is the clean unit; `matter`/`lesion_type` kept as
    # descriptive meta only (available for a within-Macnair secondary, not the bridge).
    grp = (obs["individual_id_anon"].astype(str) + "||" + obs["broad"].astype(str)) \
        .where(keep, other=np.nan)
    groups = sorted(grp.dropna().unique())
    gidx = {g: i for i, g in enumerate(groups)}
    rows_idx = [i for i, k in enumerate(grp) if isinstance(k, str)]
    cols_idx = [gidx[grp.iloc[i]] for i in rows_idx]
    G = sp.csr_matrix((np.ones(len(rows_idx)), (rows_idx, cols_idx)),
                      shape=(n_cells, len(groups)), dtype=np.float32)
    n_per = np.asarray(G.sum(axis=0)).ravel()
    print(f"  {len(groups)} (donor x broad) groups")

    print(f"[h10e] reading mtx {MTX} (genes x cells, gene-major) -- large, ~1-2 min")
    M = sio.mmread(MTX).tocsr()
    if M.shape[1] != n_cells:
        sys.exit(f"ERROR: mtx has {M.shape[1]} cells, expected {n_cells}")
    rd = pd.read_csv(ROW)
    rd.columns = [c.strip() for c in rd.columns]
    genes = rd["symbol"].astype(str).str.strip()
    if M.shape[0] != len(genes):
        sys.exit(f"ERROR: mtx has {M.shape[0]} genes, expected {len(genes)}")

    pb = (M @ G)
    pb = pd.DataFrame(np.asarray(pb.todense()), index=genes.values, columns=groups)
    pb = pb.groupby(level=0).sum().T
    pb = pb.round().astype(int)
    print(f"  pseudobulk: {pb.shape[0]} groups x {pb.shape[1]} genes")

    # donor -> diagnosis/sex luts (donor is one diagnosis, one sex). matter is per-donor
    # mixed (a donor has both GM and WM samples) -> record as comma-joined for reference.
    dlut = obs[keep].drop_duplicates("individual_id_anon").set_index("individual_id_anon")
    matter_by_donor = (obs[keep].groupby("individual_id_anon")["matter"]
                       .apply(lambda s: ",".join(sorted(s.unique()))))
    gmeta = pd.DataFrame({"group": groups})
    gmeta[["donor", "broad"]] = gmeta["group"].str.split(r"\|\|", expand=True)
    gmeta["diagnosis"] = gmeta["donor"].map(dlut["diagnosis"])
    gmeta["sex"] = gmeta["donor"].map(dlut["sex"])
    gmeta["matter"] = gmeta["donor"].map(matter_by_donor)
    gmeta["n_cells"] = n_per
    gmeta = gmeta.set_index("group")

    TAB.mkdir(parents=True, exist_ok=True)
    pb.to_parquet(TAB / "h10e_macnair_pseudobulk_primary.parquet")
    gmeta.to_csv(TAB / "h10e_macnair_group_meta_primary.csv")
    print(f"  -> {TAB / 'h10e_macnair_pseudobulk_primary.parquet'}")
    print(f"  -> {TAB / 'h10e_macnair_group_meta_primary.csv'}")

    print("\n[h10e] donors per broad x diagnosis (unit = donor):")
    tab = (gmeta.reset_index().groupby(["broad", "diagnosis"])["donor"].nunique()
           .unstack(fill_value=0))
    print(tab.to_string())
    for broad, row in tab.iterrows():
        thin = [d for d, n in row.items() if n < MIN_DONORS]
        if thin:
            print(f"  -- WARN {broad}: < {MIN_DONORS} donors in {thin} -> engine skips loudly")


if __name__ == "__main__":
    main()
