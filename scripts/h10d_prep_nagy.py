#!/usr/bin/env python
"""h10d_prep_nagy.py -- pseudobulk the Nagy 2020 MALE MDD dlPFC matrix.

Brain cross-species ARM, dataset 3 (MDD-male). Self-contained: condition + donor are encoded
in the barcode string, no external crosswalk. (This is the male half that GSE213982 folded in
as M1-M34; we use Nagy standalone to keep condition inline and avoid the unrecoverable
M#->Nagy-donor crosswalk.)

Input (data/human_validation/brain/nagy_2020_GSE144136/):
  GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz -- 30062 genes x 78886 cells (gene-major)
  GSE144136_GeneNames.csv.gz                    -- gene symbols (",x" header; row order = mtx)
  GSE144136_CellNames.csv.gz                    -- cell strings (",x" header; col order = mtx),
       format "{subtype}.{donor}_{condition}_{B#}_{barcode}"
       e.g. Ex_10_L2_4.3_Control_B3_AAACCTGAGGTAGCCA
       donor = the integer after the FIRST dot; condition in {Control, Suicide}; B# = batch.

Unit = donor (34 individuals, 17 Control / 17 Suicide; all BA9). Bridge = broad 7-class.
Broad = subtype prefix before first '_' or '.'  (Ex, Inhib, Astros, Oligos, OPCs,
Micro/Macro, Endo); Mix dropped.

Outputs (.../tables/):
  h10d_nagy_pseudobulk_primary.parquet   -- (donor x broad) x gene counts (int)
  h10d_nagy_group_meta_primary.csv       -- group -> donor, broad, diagnosis, sex, batch, n_cells

Usage (WS, from project root):
  uv run python scripts/h10d_prep_nagy.py
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

DATA = Path("data/human_validation/brain/nagy_2020_GSE144136")
MTX = DATA / "GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz"
GENES = DATA / "GSE144136_GeneNames.csv.gz"
CELLS = DATA / "GSE144136_CellNames.csv.gz"
TAB = DATA / "tables"
MIN_DONORS = 3

BROAD_MAP = {"Ex": "ExN", "Inhib": "InN", "Astros": "Ast", "Oligos": "Oli",
             "OPCs": "OPC", "Micro/Macro": "Mic", "Endo": "Endo"}   # Mix -> dropped
DIAG_MAP = {"Control": "Control", "Suicide": "MDD"}   # Suicide cohort = MDD cases

# cell string is "{prefix}[_subtype].{donor}_{cond}_{batch}_{bc}". The subtype segment is
# OPTIONAL: most types have it (Ex_10_L2_4.3_...), but Micro/Macro and Endo go straight to
# the donor-dot (Micro/Macro.3_Control_B3_..., Endo.3_...). Prefix = leading [A-Za-z/]+;
# donor = the integer immediately before _Control_/_Suicide_.
_CELL_RE = re.compile(r"^([A-Za-z/]+)(?:_\S*?)?\.(\d+)_(Control|Suicide)_(B\d+)_")


def parse_cell(s):
    """-> (broad_raw, donor, condition, batch) or (None,...) if unparseable."""
    m = _CELL_RE.match(s)
    if not m:
        return None, None, None, None
    return m.group(1), m.group(2), m.group(3), m.group(4)


def main():
    print(f"[h10d] reading cell strings {CELLS}")
    cdf = pd.read_csv(CELLS)
    cells = (cdf["x"] if "x" in cdf.columns else cdf.iloc[:, -1]).astype(str).str.strip()
    parsed = [parse_cell(s) for s in cells]
    obs = pd.DataFrame(parsed, columns=["broad_raw", "donor", "condition", "batch"])
    obs["cell"] = cells.values
    n_cells = len(obs)
    bad = obs["donor"].isna().sum()
    print(f"  {n_cells} cells; {bad} unparseable")
    if bad > 0.01 * n_cells:
        sys.exit(f"ERROR: {bad} cells failed the barcode regex (>1%) -- check format")

    obs["broad"] = obs["broad_raw"].map(BROAD_MAP)
    obs["diagnosis"] = obs["condition"].map(DIAG_MAP)
    obs["sex"] = "M"
    keep = obs["broad"].notna() & obs["diagnosis"].notna()
    print(f"  cells kept (mapped broad): {keep.sum()}/{n_cells}")
    print("  broad census (kept):\n", obs.loc[keep, "broad"].value_counts().to_string())
    print("  diagnosis census (kept cells):\n",
          obs.loc[keep, "diagnosis"].value_counts().to_string())

    grp = (obs["donor"] + "||" + obs["broad"].astype(str)).where(keep, other=np.nan)
    groups = sorted(grp.dropna().unique())
    gidx = {g: i for i, g in enumerate(groups)}
    rows_idx = [i for i, k in enumerate(grp) if isinstance(k, str)]
    cols_idx = [gidx[grp.iloc[i]] for i in rows_idx]
    G = sp.csr_matrix((np.ones(len(rows_idx)), (rows_idx, cols_idx)),
                      shape=(n_cells, len(groups)), dtype=np.float32)
    n_per = np.asarray(G.sum(axis=0)).ravel()
    print(f"  {len(groups)} (donor x broad) groups")

    print(f"[h10d] reading mtx {MTX} (genes x cells, gene-major)")
    M = sio.mmread(MTX).tocsr()
    if M.shape[1] != n_cells:
        sys.exit(f"ERROR: mtx has {M.shape[1]} cells, expected {n_cells}")
    gdf = pd.read_csv(GENES)
    genes = (gdf["x"] if "x" in gdf.columns else gdf.iloc[:, -1]).astype(str).str.strip()
    if M.shape[0] != len(genes):
        sys.exit(f"ERROR: mtx has {M.shape[0]} genes, expected {len(genes)}")

    pb = (M @ G)
    pb = pd.DataFrame(np.asarray(pb.todense()), index=genes.values, columns=groups)
    pb = pb.groupby(level=0).sum().T
    pb = pb.round().astype(int)
    print(f"  pseudobulk: {pb.shape[0]} groups x {pb.shape[1]} genes")

    # donor -> diagnosis/batch luts (a donor is one condition; batch is the modal batch)
    dlut = obs[keep].drop_duplicates("donor").set_index("donor")
    gmeta = pd.DataFrame({"group": groups})
    gmeta[["donor", "broad"]] = gmeta["group"].str.split(r"\|\|", expand=True)
    gmeta["diagnosis"] = gmeta["donor"].map(dlut["diagnosis"])
    gmeta["sex"] = "M"
    gmeta["batch"] = gmeta["donor"].map(dlut["batch"])
    gmeta["n_cells"] = n_per
    gmeta = gmeta.set_index("group")

    TAB.mkdir(parents=True, exist_ok=True)
    pb.to_parquet(TAB / "h10d_nagy_pseudobulk_primary.parquet")
    gmeta.to_csv(TAB / "h10d_nagy_group_meta_primary.csv")
    print(f"  -> {TAB / 'h10d_nagy_pseudobulk_primary.parquet'}")
    print(f"  -> {TAB / 'h10d_nagy_group_meta_primary.csv'}")

    print("\n[h10d] donors per broad x diagnosis:")
    tab = (gmeta.reset_index().groupby(["broad", "diagnosis"])["donor"].nunique()
           .unstack(fill_value=0))
    print(tab.to_string())
    for broad, row in tab.iterrows():
        thin = [d for d, n in row.items() if n < MIN_DONORS]
        if thin:
            print(f"  -- WARN {broad}: < {MIN_DONORS} donors in {thin} -> engine skips loudly")


if __name__ == "__main__":
    main()
