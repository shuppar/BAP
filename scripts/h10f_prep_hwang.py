#!/usr/bin/env python
"""h10f_prep_hwang.py -- recluster + pseudobulk the Hwang/Girgenti 2025 PTSD/MDD dlPFC deposit.

Brain cross-species ARM 5 (PTSD + MDD, the closest adult TRAUMA analog). Unlike the other four
brain arms (h10a/c/d/e = pseudobulk straight off the authors' deposited celltype labels), the
Hwang per-cell obs (celltype `anno`, Condition, Channel) was NOT deposited -- it lives in the
undeposited RNA_FINAL.zarr. So this arm reclusters from the raw counts through the SAME pipeline
as the other arms + the mouse data (own scVI -> Leiden -> marker-annotate broad-7), then pseudobulks
per donor x broad. Diagnosis/Sex are recovered from Supplementary Table 1 by channel prefix.
The authors' Supp Tables are a cross-check layer only (h10f_validate_vs_tables.py), never input here.

Input (data/human_validation/brain/hwang_ptsd/):
  RNA_count_mat.npz   -- scipy sparse (COO), (935371 cells x 27982 genes), uint32 raw counts
  RNA_cellnames.txt   -- 935371 rows, row-aligned to matrix; "{MS####XX}-{16bp barcode}"; 105 channels
  RNA_genenames.txt   -- 9-col CSV, col1 = ENSG, 27982 rows = matrix columns
  supp/Supp_Data_Table_1_Sample_metadata.txt -- tab-sep; Sample "{MS####XX}/..." -> Condition, Sex
                                                 (111 rows: 39 CON / 36 MDD / 36 PTSD; 105 join to cells)

Bridge: ENSG -> human symbol via refs/mouse_human_orthologs.tsv (human_ensembl -> human_symbol).
  This is the SAME 1:1 ortholog table the RRHO uses, so restricting to ortholog-mapped genes loses
  nothing downstream (RRHO only compares ortholog-bridged genes anyway) and keeps var aligned
  between clustering and pseudobulk. Dup symbols collapsed at pseudobulk (groupby.sum), like h10c.

Unit = donor = channel (1 library per donor in this deposit). scVI batch_key = channel.
No SoupX (filtered counts only, no raw droplet matrix -- flagged, like Admati/Maitra/Nagy/Macnair).
No doublet removal (already the authors' post-QC discovery set; re-running would over-remove and
diverge from their Table-2 census, which h10f_validate_vs_tables.py checks).

Outputs (.../tables/):
  h10f_hwang_pseudobulk_primary.parquet   -- (donor x broad) x symbol counts (int)
  h10f_hwang_group_meta_primary.csv       -- group -> donor, broad, diagnosis, sex, n_cells
  (smoke mode appends _smoke to both, so it never clobbers the real run)

Usage (WS, from project root):
  uv run python scripts/h10f_prep_hwang.py --smoke                 # 2 channels, ~5-10 min
  uv run python scripts/h10f_prep_hwang.py                         # full, GPU, tmux (~3-5 hr)
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).parent))
from _utils import add_lognorm, select_accelerator  # noqa: E402

DATA = Path("data/human_validation/brain/hwang_ptsd")
NPZ = DATA / "RNA_count_mat.npz"
CELLNAMES = DATA / "RNA_cellnames.txt"
GENES = DATA / "RNA_genenames.txt"
T1 = DATA / "supp" / "Supp_Data_Table_1_Sample_metadata.txt"
ORTHO = Path("refs/mouse_human_orthologs.tsv")
TAB = DATA / "tables"

MIN_DONORS = 3          # human pseudobulk inclusion floor per broad x diagnosis (engine skips < this)
MIN_GENES_CELL = 200    # light floor only; this is already the authors' post-QC set
DIAG_LEVELS = ["CON", "MDD", "PTSD"]   # ref = CON (set in the engine's DATASETS entry)

# lake_genes (utils/lists.py) collapsed to a discriminative broad-7 panel.
# DUSP1 (stress IEG), STMN2/RBFOX3 (pan-neuronal) dropped -- not class-discriminative.
BROAD_MARKERS = {
    "Mic":  ["MRC1", "TMEM119", "CX3CR1", "APBB1IP"],
    "Endo": ["CLDN5", "FLT1", "COBLL1"],
    "Ast":  ["AQP4", "GJA1", "GFAP", "ALDH1L1", "SLC4A4", "NDRG2", "GLUL", "SOX9", "ALDH1A1", "VIM"],
    "OPC":  ["PDGFRA", "PCDH15", "OLIG1", "OLIG2"],
    "Oli":  ["PLP1", "MAG", "MOG", "MOBP", "MBP"],
    "ExN":  ["SATB2", "SLC17A7", "GRM4"],
    "InN":  ["GAD1", "GAD2", "SLC32A1", "SST", "PVALB"],
}
# HVG exclusion (symbol-based): mito, ribo, sex-linked -- kept OUT of scVI feature selection.
SEX_LINKED = {"XIST", "RPS4Y1", "DDX3Y", "UTY", "USP9Y", "ZFY", "KDM5D", "EIF1AY", "NLGN4Y", "TSIX"}


def load_table1():
    """channel -> (Condition, Sex) from Supp Table 1. Channel = first token of `Sample`."""
    t1 = pd.read_csv(T1, sep="\t")
    need = {"Sample", "Condition", "Sex"}
    miss = need - set(t1.columns)
    if miss:
        sys.exit(f"ERROR: Supp Table 1 missing columns {miss}; has {list(t1.columns)}")
    t1["channel"] = t1["Sample"].astype(str).str.split("/", n=1).str[0].str.strip()
    # no-silent-failure: a channel with two different Conditions must NOT be silently deduped
    conflict = t1.groupby("channel")["Condition"].nunique()
    conflict = conflict[conflict > 1]
    if len(conflict):
        sys.exit(f"ERROR: {len(conflict)} channels have conflicting Condition rows in Table 1 "
                 f"({list(conflict.index[:10])}) -- dedup would pick one arbitrarily")
    t1 = t1.drop_duplicates("channel").set_index("channel")
    print(f"[h10f] Table 1: {len(t1)} channels (full sample sheet; RNA is a subset); "
          f"Condition census:\n{t1['Condition'].value_counts().to_string()}")
    bad = set(t1["Condition"].unique()) - set(DIAG_LEVELS)
    if bad:
        sys.exit(f"ERROR: unexpected Condition values {bad}; expected {DIAG_LEVELS}")
    return t1["Condition"].to_dict(), t1["Sex"].astype(str).to_dict()


def build_adata(args):
    import anndata as ad

    print(f"[h10f] cellnames {CELLNAMES}")
    cellnames = pd.read_csv(CELLNAMES, header=None).iloc[:, 0].astype(str).str.strip()
    channel = cellnames.str.split("-", n=1).str[0]

    print(f"[h10f] genenames {GENES}")
    gtab = pd.read_csv(GENES, header=None, sep=None, engine="python")
    ensg = gtab.iloc[:, 0].astype(str).str.strip()

    cond_map, sex_map = load_table1()

    # --- channel resolution BEFORE loading the big matrix (fail fast, no silent drops) ---
    present = pd.Index(channel.unique())
    unresolved = [c for c in present if c not in cond_map]
    if unresolved:
        sys.exit(f"ERROR: {len(unresolved)} cell channels not in Table 1 (silent-drop guard): "
                 f"{unresolved[:10]}")
    print(f"[h10f] all {len(present)} cell channels resolve against Table 1")
    pres_cond = pd.Series({c: cond_map[c] for c in present})
    print(f"[h10f] diagnosis census of the {len(present)} RNA cell channels:\n"
          f"{pres_cond.value_counts().to_string()}")

    # --- smoke: pick 2 present channels spanning >=2 conditions ---
    if args.smoke:
        by_cond = {}
        for c in present:
            by_cond.setdefault(cond_map[c], []).append(c)
        pick = []
        for lvl in DIAG_LEVELS:                       # prefer CON + (MDD/PTSD)
            if by_cond.get(lvl):
                pick.append(by_cond[lvl][0])
            if len(pick) == 2:
                break
        if len(pick) < 2:                             # fallback: first two present
            pick = list(present[:2])
        print(f"[h10f] SMOKE channels: {pick} "
              f"({[cond_map[c] for c in pick]})")
        cell_keep = channel.isin(pick).values
    else:
        cell_keep = np.ones(len(channel), dtype=bool)

    print(f"[h10f] loading matrix {NPZ} (COO -> CSR); this is the big step")
    M = sp.load_npz(NPZ).tocsr()                       # (cells x genes)
    if M.shape[0] != len(cellnames):
        sys.exit(f"ERROR: matrix {M.shape[0]} cells != {len(cellnames)} cellnames")
    if M.shape[1] != len(ensg):
        sys.exit(f"ERROR: matrix {M.shape[1]} genes != {len(ensg)} genenames")
    print(f"  matrix {M.shape[0]:,} cells x {M.shape[1]:,} genes")

    if not cell_keep.all():
        M = M[cell_keep]
        cellnames = cellnames[cell_keep].reset_index(drop=True)
        channel = channel[cell_keep].reset_index(drop=True)
    print(f"  after cell subset: {M.shape[0]:,} cells")

    # --- restrict genes to ortholog-mapped, carry symbol ---
    o = pd.read_csv(ORTHO, sep="\t")
    e2s = dict(o[["human_ensembl", "human_symbol"]].dropna().drop_duplicates("human_ensembl").values)
    gene_sym = ensg.map(e2s)
    gene_keep = gene_sym.notna().values
    n_map = int(gene_keep.sum())
    print(f"  ENSG->symbol via ortholog table: {n_map}/{len(ensg)} genes mapped (kept)")
    if n_map < 5000:
        sys.exit(f"ERROR: only {n_map} genes mapped -- check ortholog human_ensembl join")
    M = M[:, gene_keep]
    ensg_keep = ensg[gene_keep].reset_index(drop=True)
    sym_keep = gene_sym[gene_keep].reset_index(drop=True)

    obs = pd.DataFrame({
        "channel": channel.values,
        "donor": channel.values,                       # 1 library per donor in this deposit
        "diagnosis": channel.map(cond_map).values,
        "sex": channel.map(sex_map).values,
    }, index=cellnames.values)
    var = pd.DataFrame({"symbol": sym_keep.values}, index=ensg_keep.values)   # var_names = ENSG (unique)
    var.index.name = "ensembl"

    A = ad.AnnData(X=M.astype(np.float32), obs=obs, var=var)
    print(f"[h10f] AnnData {A.n_obs:,} x {A.n_vars:,}")
    print("  diagnosis census (cells):\n" + A.obs["diagnosis"].value_counts().to_string())
    print("  sex census (cells):\n" + A.obs["sex"].value_counts().to_string())
    return A


def light_qc(A):
    import scanpy as sc
    A.var["mt"] = A.var["symbol"].str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    before = A.n_obs
    A = A[A.obs["n_genes_by_counts"] >= MIN_GENES_CELL].copy()
    print(f"[h10f] light QC (min {MIN_GENES_CELL} genes): {A.n_obs:,}/{before:,} cells kept "
          f"(post-QC deposit -> few expected)")
    keep_g = np.asarray((A.X > 0).sum(axis=0)).ravel() > 0
    A = A[:, keep_g].copy()
    print(f"  dropped all-zero genes -> {A.n_vars:,} genes")
    return A


def integrate(A, args):
    import scanpy as sc
    import scvi

    n_hvg = min(args.n_hvg, A.n_vars - 1)
    print(f"[h10f] HVG (seurat_v3, n_top={n_hvg}, batch=channel)")
    sc.pp.highly_variable_genes(A, flavor="seurat_v3", n_top_genes=n_hvg, batch_key="channel")
    excl = (A.var["mt"].values
            | A.var["symbol"].str.upper().str.startswith(("RPS", "RPL")).values
            | A.var["symbol"].str.upper().isin(SEX_LINKED).values)
    A.var["use_for_scvi"] = A.var["highly_variable"].values & ~excl
    print(f"  HVGs for scVI after mito/ribo/sex exclusion: {int(A.var['use_for_scvi'].sum())}")

    accel, prec = select_accelerator(force_cpu=args.cpu)
    max_epochs = args.max_epochs if args.max_epochs else (30 if args.smoke else 400)
    scvi.settings.seed = 42
    if accel == "gpu":
        import torch
        torch.set_float32_matmul_precision("high")   # Tensor-core hint (Lightning warned)
        # scVI is GPU-bound: a modest worker count keeps the loader ahead of the GPU without
        # eating the box. NOT Lightning's naive n_cores-1 (=55) -- that spawns overhead/RAM for
        # ~0 gain and blows the CPU margin. persistent_workers avoids per-epoch respawn cost.
        scvi.settings.dl_num_workers = args.n_workers
        scvi.settings.dl_persistent_workers = True
        print(f"[h10f] scVI dataloader: num_workers={args.n_workers}, persistent=True")
    Ah = A[:, A.var["use_for_scvi"]].copy()
    scvi.model.SCVI.setup_anndata(Ah, batch_key="channel",
                                  continuous_covariate_keys=["pct_counts_mt"])
    print(f"[h10f] scVI train (accel={accel}, prec={prec}, max_epochs={max_epochs}, "
          f"batch_size={args.batch_size})")
    model = scvi.model.SCVI(Ah, n_layers=2, n_latent=30)
    model.train(max_epochs=max_epochs, batch_size=args.batch_size,
                early_stopping=True, early_stopping_patience=30,
                accelerator=accel, devices=1, precision=prec)
    A.obsm["X_scVI"] = model.get_latent_representation()

    print(f"[h10f] neighbors + Leiden (igraph, res={args.resolution})")
    sc.pp.neighbors(A, use_rep="X_scVI", n_neighbors=30, random_state=42)
    sc.tl.leiden(A, flavor="igraph", n_iterations=2, resolution=args.resolution, random_state=42)
    print(f"  {A.obs['leiden'].nunique()} Leiden clusters")
    return A


def annotate(A):
    """Cluster-majority broad-7 via lake_genes marker scores (mouse-07 / h09d pattern)."""
    import scanpy as sc

    sym2ensg = {s: e for e, s in zip(A.var_names, A.var["symbol"])}
    marker_ensg = {}
    for cls, syms in BROAD_MARKERS.items():
        present = [sym2ensg[s] for s in syms if s in sym2ensg]
        marker_ensg[cls] = present
        print(f"  markers {cls}: {len(present)}/{len(syms)} present "
              f"({[s for s in syms if s not in sym2ensg]} missing)")
    ok = sum(1 for v in marker_ensg.values() if v)
    if ok < 5:
        sys.exit(f"ERROR: only {ok}/7 broad classes have any marker present -- panel/join broken")

    add_lognorm(A)
    Xraw = A.X
    A.X = A.layers["lognorm"]
    score_cols = []
    for cls, ens in marker_ensg.items():
        if not ens:
            continue
        sc.tl.score_genes(A, ens, score_name=f"score_{cls}", use_raw=False)
        score_cols.append(f"score_{cls}")
    A.X = Xraw
    del A.layers["lognorm"]

    cl = A.obs.groupby("leiden")[score_cols].mean()
    cl_broad = cl.idxmax(axis=1).str.replace("score_", "", regex=False)
    print("[h10f] cluster -> broad (mean marker score argmax):")
    print(pd.concat([cl.round(3), cl_broad.rename("broad")], axis=1).to_string())
    A.obs["broad"] = A.obs["leiden"].map(cl_broad).astype(str)
    print("[h10f] broad census (cells):\n" + A.obs["broad"].value_counts().to_string())
    return A


def pseudobulk(A, suffix):
    keep = A.obs["broad"].notna() & A.obs["diagnosis"].notna()
    grp = (A.obs["donor"].astype(str) + "||" + A.obs["broad"].astype(str)).where(keep)
    groups = sorted(grp.dropna().unique())
    gidx = {g: i for i, g in enumerate(groups)}
    rows = np.where(keep.values)[0]
    cols = grp.iloc[rows].map(gidx).values
    G = sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)),
                      shape=(A.n_obs, len(groups)))
    n_per = np.asarray(G.sum(axis=0)).ravel()
    print(f"[h10f] {len(groups)} (donor x broad) pseudobulk groups")

    pb = (G.T @ A.X)                                   # (groups x genes), sparse
    pb = pd.DataFrame(np.asarray(pb.todense()), index=groups, columns=A.var["symbol"].values)
    pb = pb.T.groupby(level=0).sum().T                 # collapse dup symbols
    pb = pb.round().astype(int)
    print(f"  pseudobulk {pb.shape[0]} groups x {pb.shape[1]} symbols")

    meta = pd.DataFrame({"group": groups})
    meta[["donor", "broad"]] = meta["group"].str.split(r"\|\|", expand=True)
    dmap = A.obs.drop_duplicates("donor").set_index("donor")
    meta["diagnosis"] = meta["donor"].map(dmap["diagnosis"])
    meta["sex"] = meta["donor"].map(dmap["sex"])
    meta["n_cells"] = n_per
    meta = meta.set_index("group")

    TAB.mkdir(parents=True, exist_ok=True)
    pb_path = TAB / f"h10f_hwang_pseudobulk_primary{suffix}.parquet"
    mt_path = TAB / f"h10f_hwang_group_meta_primary{suffix}.csv"
    pb.to_parquet(pb_path)
    meta.to_csv(mt_path)
    print(f"  -> {pb_path}")
    print(f"  -> {mt_path}")

    print("\n[h10f] donors per broad x diagnosis:")
    tab = (meta.reset_index().groupby(["broad", "diagnosis"])["donor"].nunique()
           .unstack(fill_value=0))
    print(tab.to_string())
    for broad, row in tab.iterrows():
        thin = [d for d, n in row.items() if 0 < n < MIN_DONORS]
        if thin:
            print(f"  -- WARN {broad}: < {MIN_DONORS} donors in {thin} -> engine skips loudly")
    if meta.empty:
        sys.exit("ERROR: empty pseudobulk -- no (donor x broad) groups produced")


def main():
    ap = argparse.ArgumentParser(description="Prep Hwang PTSD/MDD (recluster -> pseudobulk)")
    ap.add_argument("--smoke", action="store_true", help="2 channels spanning >=2 conditions")
    ap.add_argument("--cpu", action="store_true", help="force CPU (scVI)")
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--batch-size", type=int, default=2048,
                    help="scVI minibatch. VRAM-safe up to ~4096 on 2000 HVG / 24GB; 2048 ~halves "
                         "steps/epoch vs 1024. (scVI default 128 is ~16x slower at atlas scale.)")
    ap.add_argument("--n-workers", type=int, default=12,
                    help="scVI dataloader workers (GPU-bound: 8-12 saturates the feed; leaves "
                         ">40 cores free. NOT Lightning's naive 55.)")
    ap.add_argument("--resolution", type=float, default=1.0)
    ap.add_argument("--max-epochs", type=int, default=None, help="override (smoke default 30)")
    args = ap.parse_args()

    print(f"\n=== h10f prep Hwang PTSD/MDD (smoke={args.smoke}) ===")
    for p in (NPZ, CELLNAMES, GENES, T1, ORTHO):
        if not p.exists():
            sys.exit(f"ERROR: missing {p}")

    A = build_adata(args)
    A = light_qc(A)
    A = integrate(A, args)
    A = annotate(A)
    pseudobulk(A, suffix="_smoke" if args.smoke else "")
    print("\n[h10f] done." + ("  (SMOKE outputs suffixed _smoke)" if args.smoke else ""))


if __name__ == "__main__":
    main()
