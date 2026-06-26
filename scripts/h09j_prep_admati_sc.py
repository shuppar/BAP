#!/usr/bin/env python
"""h09j_prep_admati_sc.py -- pseudobulk the Admati sc all-compartment PE matrix.

Input: sc_PE_allcells_with_metadata_29-May-2023.txt -- TRANSPOSED, 23 metadata rows then
gene rows; ~98k cells in columns (too wide to read whole -> STREAM gene rows in chunks and
accumulate pseudobulk via a sparse (cells x groups) indicator). Author-annotated, so no
clustering/SoupX -- we pseudobulk on their labels.

Compartment from celltype prefix (clean): TB_*->trophoblast, STROMAL_*->decidua_stromal,
VASCULAR_*->vascular, IMMUNE_*->immune. (No erythroid in this dataset -- already dropped
from the mouse RRHO.) Condition from the 4 indicator rows early/late_control/PE.

Output (reused by h09k):
  h09j_admati_pseudobulk.parquet   -- group (donor x compartment) x gene counts (int)
  h09j_admati_group_meta.csv       -- group -> donorID, compartment, condition, n_cells

Usage (from project root):
  uv run python scripts/h09j_prep_admati_sc.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).parent))
from h09e_cross_species_rrho import OUT_DIR  # noqa: E402

SC = Path("data/human_validation/placenta/admati_2023_PE/"
          "sc_PE_allcells_with_metadata_29-May-2023.txt")
N_META = 23   # cellID..donor_age (genes start row 24)
COND_ROWS = ["early_control", "late_control", "early_PE", "late_PE"]
PREFIX2COMP = {"TB": "trophoblast", "STROMAL": "decidua_stromal",
               "VASCULAR": "vascular", "IMMUNE": "immune"}
CHUNK = 2000
TAB = OUT_DIR / "tables"


def compartment_of(celltype):
    return PREFIX2COMP.get(str(celltype).split("_", 1)[0], None)


def main():
    print(f"[h09j] reading {N_META} metadata rows from {SC}")
    meta = pd.read_csv(SC, sep="\t", nrows=N_META, header=None, index_col=0, low_memory=False)
    meta.index = meta.index.astype(str).str.strip()
    obs = meta.T
    obs.columns = [c.strip() for c in obs.columns]
    n_cells = len(obs)
    print(f"  {n_cells} cells; metadata fields: {list(obs.columns)[:8]}...")

    # per-cell compartment + condition
    obs["compartment"] = obs["celltype"].map(compartment_of)
    def cond_of(row):
        for c in COND_ROWS:
            if str(row.get(c, "0")).strip() in ("1", "1.0"):
                return c
        return None
    obs["condition"] = obs.apply(cond_of, axis=1)
    obs["donorID"] = obs["donorID"].astype(str).str.strip()

    keep = obs["compartment"].notna() & obs["condition"].notna()
    print(f"  cells in kept compartments+conditions: {keep.sum()}/{n_cells}")
    print("  compartment census:\n",
          obs.loc[keep, "compartment"].value_counts().to_string())
    print("  condition census (cells):\n",
          obs.loc[keep, "condition"].value_counts().to_string())

    # group = (donor, compartment); build sparse indicator (cells x groups)
    grp_key = obs["donorID"] + "||" + obs["compartment"].astype(str)
    grp_key = grp_key.where(keep, other=np.nan)
    groups = sorted(grp_key.dropna().unique())
    gidx = {g: i for i, g in enumerate(groups)}
    rows_idx = [i for i, k in enumerate(grp_key) if isinstance(k, str)]
    cols_idx = [gidx[grp_key.iloc[i]] for i in rows_idx]
    G = sp.csr_matrix((np.ones(len(rows_idx)), (rows_idx, cols_idx)),
                      shape=(n_cells, len(groups)), dtype=np.float32)
    n_cells_per_group = np.asarray(G.sum(axis=0)).ravel()
    print(f"  {len(groups)} (donor x compartment) groups")

    # stream gene rows, accumulate genes x groups
    print(f"[h09j] streaming gene rows (chunksize={CHUNK})")
    parts, genes = [], []
    reader = pd.read_csv(SC, sep="\t", skiprows=N_META, header=None, index_col=0,
                         chunksize=CHUNK, low_memory=False)
    seen = 0
    for chunk in reader:
        vals = chunk.to_numpy(dtype=np.float32)
        if vals.shape[1] != n_cells:
            sys.exit(f"ERROR: gene chunk has {vals.shape[1]} cell cols, expected {n_cells}")
        np.nan_to_num(vals, copy=False)
        parts.append(vals @ G)                 # (chunk_genes x groups)
        genes.extend(chunk.index.astype(str).str.strip().tolist())
        seen += vals.shape[0]
        print(f"    {seen} genes processed", end="\r")
    print()
    pb = pd.DataFrame(np.vstack(parts), index=genes, columns=groups)
    pb = pb.groupby(level=0).sum()             # collapse duplicate gene symbols
    pb = pb.T                                   # groups x genes
    pb = pb.round().astype(int)
    print(f"  pseudobulk: {pb.shape[0]} groups x {pb.shape[1]} genes")

    # group metadata
    gmeta = pd.DataFrame({"group": groups})
    gmeta[["donorID", "compartment"]] = gmeta["group"].str.split(r"\|\|", expand=True)
    cond_map = obs.drop_duplicates("donorID").set_index("donorID")["condition"].to_dict()
    gmeta["condition"] = gmeta["donorID"].map(cond_map)
    gmeta["n_cells"] = n_cells_per_group
    gmeta = gmeta.set_index("group")

    TAB.mkdir(parents=True, exist_ok=True)
    pb.to_parquet(TAB / "h09j_admati_pseudobulk.parquet")
    gmeta.to_csv(TAB / "h09j_admati_group_meta.csv")
    print(f"\n[h09j] pseudobulk -> {TAB/'h09j_admati_pseudobulk.parquet'}")
    print(f"[h09j] group meta -> {TAB/'h09j_admati_group_meta.csv'}")

    # donor counts per condition x compartment (the powering check)
    print("\n[h09j] donors per condition x compartment:")
    tab = gmeta.reset_index().groupby(["compartment", "condition"])["donorID"].nunique().unstack(fill_value=0)
    print(tab.to_string())


if __name__ == "__main__":
    main()
