#!/usr/bin/env python
"""h10c_prep_maitra.py -- pseudobulk the Maitra 2023 FEMALE MDD dlPFC matrix.

Brain cross-species ARM, dataset 2 (MDD-female). The GSE213982 combined matrix holds both
sexes (F1-F38 female + M1-M34 male, the latter realigned from Nagy). We use the FEMALE
subset only; the male arm is Nagy standalone (h10d), which keeps condition inline and needs
no external crosswalk.

Input (data/human_validation/brain/maitra_2023_GSE213982/):
  GSE213982_combined_counts_matrix.mtx.gz            -- 36588 genes x 160711 cells (gene-major)
  GSE213982_combined_counts_matrix_genes_rows.csv.gz -- gene symbols (row order = mtx rows)
  GSE213982_combined_counts_matrix_cells_columns.csv.gz -- cell strings (col order = mtx cols),
       format: "{donor}.{barcode}.{broad}.{subtype}"  e.g. F1.AAACCCACACCTCTGT-1.Mic.Mic1
  maitra_donor_meta.csv  -- donor,gsm,sex,group,platform ; group in {Case, Control}

Unit = donor (all BA9 dlPFC, no region split -> donor is the clean unit, unlike Velmeshev).
Bridge = broad 7-class; Mix dropped.

Outputs (.../tables/):
  h10c_maitra_pseudobulk_primary.parquet   -- (donor x broad) x gene counts (int)
  h10c_maitra_group_meta_primary.csv       -- group -> donor, broad, diagnosis, sex, n_cells

Usage (WS, from project root):
  uv run python scripts/h10c_prep_maitra.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

DATA = Path("data/human_validation/brain/maitra_2023_GSE213982")
MTX = DATA / "GSE213982_combined_counts_matrix.mtx.gz"
GENES = DATA / "GSE213982_combined_counts_matrix_genes_rows.csv.gz"
CELLS = DATA / "GSE213982_combined_counts_matrix_cells_columns.csv.gz"
DONOR_META = DATA / "maitra_donor_meta.csv"
TAB = DATA / "tables"
MIN_DONORS = 3

BROAD_MAP = {"ExN": "ExN", "InN": "InN", "Oli": "Oli", "Ast": "Ast",
             "OPC": "OPC", "End": "Endo", "Mic": "Mic"}   # Mix -> dropped (absent here)
# condition label: 'Case' -> ASD-style 'test', 'Control' -> ref. Engine reads diagnosis col.
DIAG_MAP = {"Case": "MDD", "Control": "Control"}


def parse_cell(s):
    """'F1.AAACCCACACCTCTGT-1.Mic.Mic1' -> (donor='F1', broad='Mic')."""
    parts = s.strip().strip('"').split(".")
    return parts[0], parts[-2]   # donor, broad (3rd-from-? -> NF-1)


def main():
    print(f"[h10c] reading cell strings {CELLS}")
    cells = pd.read_csv(CELLS)["x"] if "x" in pd.read_csv(CELLS, nrows=1).columns \
        else pd.read_csv(CELLS, header=None).iloc[:, -1]
    cells = cells.astype(str).str.strip().str.strip('"')
    donors, broads = zip(*(parse_cell(s) for s in cells))
    obs = pd.DataFrame({"cell": cells.values, "donor": donors, "broad_raw": broads})
    n_cells = len(obs)
    print(f"  {n_cells} cells")

    # female donors only + condition join
    dm = pd.read_csv(DONOR_META)
    dm["donor"] = dm["donor"].astype(str).str.strip()
    fem = dm[dm["sex"].str.lower() == "female"].set_index("donor")
    obs["is_female"] = obs["donor"].isin(fem.index)
    obs["diagnosis"] = obs["donor"].map(fem["group"].map(DIAG_MAP))
    obs["sex"] = "F"
    obs["broad"] = obs["broad_raw"].map(BROAD_MAP)

    keep = obs["is_female"] & obs["broad"].notna() & obs["diagnosis"].notna()
    print(f"  female cells kept (mapped broad): {keep.sum()}/{n_cells}")
    print("  broad census (kept):\n", obs.loc[keep, "broad"].value_counts().to_string())
    print("  diagnosis census (kept cells):\n",
          obs.loc[keep, "diagnosis"].value_counts().to_string())

    # group = (donor, broad); sparse indicator (cells x groups)
    grp = (obs["donor"] + "||" + obs["broad"].astype(str)).where(keep, other=np.nan)
    groups = sorted(grp.dropna().unique())
    gidx = {g: i for i, g in enumerate(groups)}
    rows_idx = [i for i, k in enumerate(grp) if isinstance(k, str)]
    cols_idx = [gidx[grp.iloc[i]] for i in rows_idx]
    G = sp.csr_matrix((np.ones(len(rows_idx)), (rows_idx, cols_idx)),
                      shape=(n_cells, len(groups)), dtype=np.float32)
    n_per = np.asarray(G.sum(axis=0)).ravel()
    print(f"  {len(groups)} (donor x broad) groups")

    print(f"[h10c] reading mtx {MTX} (genes x cells, gene-major)")
    M = sio.mmread(MTX).tocsr()   # genes x cells
    if M.shape[1] != n_cells:
        sys.exit(f"ERROR: mtx has {M.shape[1]} cells, expected {n_cells}")
    genes = pd.read_csv(GENES)["x"] if "x" in pd.read_csv(GENES, nrows=1).columns \
        else pd.read_csv(GENES, header=None).iloc[:, -1]
    genes = genes.astype(str).str.strip().str.strip('"')
    if M.shape[0] != len(genes):
        sys.exit(f"ERROR: mtx has {M.shape[0]} genes, expected {len(genes)}")

    # pseudobulk = genes x groups  (M @ G), then transpose -> groups x genes
    pb = (M @ G)                                  # (genes x groups), sparse
    pb = pd.DataFrame(np.asarray(pb.todense()), index=genes.values, columns=groups)
    pb = pb.groupby(level=0).sum().T              # collapse dup symbols, groups x genes
    pb = pb.round().astype(int)
    print(f"  pseudobulk: {pb.shape[0]} groups x {pb.shape[1]} genes")

    gmeta = pd.DataFrame({"group": groups})
    gmeta[["donor", "broad"]] = gmeta["group"].str.split(r"\|\|", expand=True)
    gmeta["diagnosis"] = gmeta["donor"].map(fem["group"].map(DIAG_MAP))
    gmeta["sex"] = "F"
    gmeta["n_cells"] = n_per
    gmeta = gmeta.set_index("group")

    TAB.mkdir(parents=True, exist_ok=True)
    pb.to_parquet(TAB / "h10c_maitra_pseudobulk_primary.parquet")
    gmeta.to_csv(TAB / "h10c_maitra_group_meta_primary.csv")
    print(f"  -> {TAB / 'h10c_maitra_pseudobulk_primary.parquet'}")
    print(f"  -> {TAB / 'h10c_maitra_group_meta_primary.csv'}")

    print("\n[h10c] donors per broad x diagnosis:")
    tab = (gmeta.reset_index().groupby(["broad", "diagnosis"])["donor"].nunique()
           .unstack(fill_value=0))
    print(tab.to_string())
    for broad, row in tab.iterrows():
        thin = [d for d, n in row.items() if n < MIN_DONORS]
        if thin:
            print(f"  -- WARN {broad}: < {MIN_DONORS} donors in {thin} -> engine skips loudly")


if __name__ == "__main__":
    main()
