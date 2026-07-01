#!/usr/bin/env python
"""h10a_prep_velmeshev.py -- pseudobulk the Velmeshev 2019 ASD cortex matrix.

Brain cross-species ARM, dataset 1 (PATHFINDER). Builds per-(donor x broad-celltype)
pseudobulk from the UCSC Cell Browser dump, ready for the shared RRHO engine (h10*).

Input (data/human_validation/brain/velmeshev_2019_autism/):
  exprMatrix.tsv.gz  -- TRANSPOSED-wide: 1 header row of cell IDs, then GENE rows.
                        col 0 = gene_ID as "ENSG00000227232|WASH7P" (Ensembl|SYMBOL).
                        36502 genes x 104559 cells. Too wide to read whole -> STREAM
                        gene rows in chunks, accumulate via sparse (cells x groups) indicator
                        (same pattern as h09j).
  meta.tsv           -- per-cell metadata, ROW ORDER IDENTICAL to exprMatrix column order
                        (verified). Columns: cell, cluster, sample, individual, region, age,
                        sex, diagnosis, ...  -> align by POSITION, no per-cell join needed.

Bridge = broad 7-class {ExN, InN, Ast, Oli, OPC, Mic, Endo}. Velmeshev's 17 fine clusters
map via BROAD_PRIMARY. Ambiguous neurons (Neu-NRGN-I/-II, Neu-mat) are EXCLUDED from the
primary bridge and emitted as a separate SENSITIVITY variant (quarantined + README) per the
locked design decision.

Outputs (data/human_validation/brain/velmeshev_2019_autism/tables/):
  h10a_velmeshev_pseudobulk_primary.parquet     -- (donor x broad) x gene counts (int)
  h10a_velmeshev_pseudobulk_sensitivity.parquet -- same + Neu-NRGN/Neu-mat -> ExN
  h10a_velmeshev_group_meta_{primary,sensitivity}.csv
                                                -- group -> individual, diagnosis, sex, region, age, n_cells

Usage (from project root, WS):
  uv run python scripts/h10a_prep_velmeshev.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

DATA = Path("data/human_validation/brain/velmeshev_2019_autism")
EXPR = DATA / "exprMatrix.tsv.gz"
META = DATA / "meta.tsv"
TAB = DATA / "tables"
CHUNK = 2000

MIN_DONORS = 3  # per (broad, diagnosis); skip a broad type loudly if either arm is thinner

# --- cluster -> broad map (locked) -------------------------------------------
# Primary: drop the three ambiguous neuronal clusters (immature / NRGN).
BROAD_PRIMARY = {
    "L2/3": "ExN", "L4": "ExN", "L5/6": "ExN", "L5/6-CC": "ExN",
    "IN-PV": "InN", "IN-SST": "InN", "IN-VIP": "InN", "IN-SV2C": "InN",
    "AST-PP": "Ast", "AST-FB": "Ast",
    "Oligodendrocytes": "Oli",
    "OPC": "OPC",
    "Microglia": "Mic",
    "Endothelial": "Endo",
    # Neu-NRGN-I, Neu-NRGN-II, Neu-mat: intentionally absent (primary excludes)
}
# Sensitivity: same, plus the ambiguous neurons folded into ExN.
SENSITIVITY_EXTRA = {"Neu-NRGN-I": "ExN", "Neu-NRGN-II": "ExN", "Neu-mat": "ExN"}
BROAD_SENSITIVITY = {**BROAD_PRIMARY, **SENSITIVITY_EXTRA}

# Pseudobulk unit = `sample` (= individual x region, e.g. "1823_BA24"), mirroring the
# mouse organization (donor=sample, region recorded separately). `individual` (the person)
# is carried alongside so the person<->sample link is never lost (a donor contributing both
# PFC and ACC appears as two samples sharing one `individual`).
UNIT_COL = "sample"
# per-cell metadata columns carried to group meta (one value per UNIT)
DONOR_META_COLS = ["individual", "diagnosis", "sex", "region", "age"]


def symbol_of(gene_id: str) -> str:
    """'ENSG00000227232|WASH7P' -> 'WASH7P'. No pipe -> return as-is."""
    s = str(gene_id)
    return s.split("|", 1)[1] if "|" in s else s


def build_pseudobulk(obs: pd.DataFrame, cluster2broad: dict, n_cells: int):
    """Stream gene rows of EXPR, accumulate (donor x broad) x gene counts.

    Returns (pb_groups_x_genes_int, group_meta_df). Cells whose cluster isn't in
    cluster2broad are dropped (their indicator column is absent).
    """
    broad = obs["cluster"].map(cluster2broad)
    unit = obs[UNIT_COL].astype(str)
    grp_key = (unit + "||" + broad.astype(str))
    keep = broad.notna()
    grp_key = grp_key.where(keep, other=np.nan)

    groups = sorted(grp_key.dropna().unique())
    gidx = {g: i for i, g in enumerate(groups)}
    rows_idx = [i for i, k in enumerate(grp_key) if isinstance(k, str)]
    cols_idx = [gidx[grp_key.iloc[i]] for i in rows_idx]
    G = sp.csr_matrix((np.ones(len(rows_idx)), (rows_idx, cols_idx)),
                      shape=(n_cells, len(groups)), dtype=np.float32)
    n_cells_per_group = np.asarray(G.sum(axis=0)).ravel()
    print(f"    {keep.sum()}/{n_cells} cells kept; {len(groups)} (donor x broad) groups")

    # stream gene rows; first row is the cell-ID header (skiprows=1)
    parts, genes = [], []
    reader = pd.read_csv(EXPR, sep="\t", skiprows=1, header=None, index_col=0,
                         chunksize=CHUNK, low_memory=False)
    seen = 0
    for chunk in reader:
        vals = chunk.to_numpy(dtype=np.float32)
        if vals.shape[1] != n_cells:
            sys.exit(f"ERROR: gene chunk has {vals.shape[1]} cell cols, expected {n_cells} "
                     f"(ragged parse / header misalignment)")
        np.nan_to_num(vals, copy=False)
        parts.append(vals @ G)                       # (chunk_genes x groups)
        genes.extend(symbol_of(g) for g in chunk.index)
        seen += vals.shape[0]
        print(f"      {seen} genes processed", end="\r")
    print()

    pb = pd.DataFrame(np.vstack(parts), index=genes, columns=groups)
    pb = pb.groupby(level=0).sum()                   # collapse duplicate symbols
    pb = pb.T.round().astype(int)                    # groups x genes

    gmeta = pd.DataFrame({"group": groups})
    gmeta[[UNIT_COL, "broad"]] = gmeta["group"].str.split(r"\|\|", expand=True)
    donor_lut = obs.copy()
    donor_lut[UNIT_COL] = donor_lut[UNIT_COL].astype(str)
    donor_lut = donor_lut.drop_duplicates(UNIT_COL).set_index(UNIT_COL)
    for c in DONOR_META_COLS:
        gmeta[c] = gmeta[UNIT_COL].map(donor_lut[c])
    gmeta["n_cells"] = n_cells_per_group
    gmeta = gmeta.set_index("group")
    return pb, gmeta


def powering_report(gmeta: pd.DataFrame, tag: str):
    print(f"\n[h10a] ({tag}) {UNIT_COL}s per broad x diagnosis:")
    tab = (gmeta.reset_index()
           .groupby(["broad", "diagnosis"])[UNIT_COL].nunique()
           .unstack(fill_value=0))
    print(tab.to_string())
    # loud skip warning for any broad type too thin in either arm
    for broad, row in tab.iterrows():
        thin = [d for d, n in row.items() if n < MIN_DONORS]
        if thin:
            print(f"  -- WARN {broad}: < {MIN_DONORS} donors in {thin} "
                  f"-> RRHO engine will skip this compartment loudly")


def main():
    print(f"[h10a] reading metadata {META}")
    obs = pd.read_csv(META, sep="\t", low_memory=False)
    obs.columns = [c.strip() for c in obs.columns]
    n_cells = len(obs)
    print(f"  {n_cells} cells; columns: {list(obs.columns)[:8]}...")
    print("  cluster census:\n", obs["cluster"].value_counts().to_string())
    print("  diagnosis census:\n", obs["diagnosis"].value_counts().to_string())

    TAB.mkdir(parents=True, exist_ok=True)

    for tag, cmap in (("primary", BROAD_PRIMARY), ("sensitivity", BROAD_SENSITIVITY)):
        print(f"\n[h10a] === building {tag} pseudobulk ===")
        pb, gmeta = build_pseudobulk(obs, cmap, n_cells)
        print(f"    pseudobulk: {pb.shape[0]} groups x {pb.shape[1]} genes")
        pb.to_parquet(TAB / f"h10a_velmeshev_pseudobulk_{tag}.parquet")
        gmeta.to_csv(TAB / f"h10a_velmeshev_group_meta_{tag}.csv")
        print(f"    -> {TAB / f'h10a_velmeshev_pseudobulk_{tag}.parquet'}")
        print(f"    -> {TAB / f'h10a_velmeshev_group_meta_{tag}.csv'}")
        powering_report(gmeta, tag)


if __name__ == "__main__":
    main()
