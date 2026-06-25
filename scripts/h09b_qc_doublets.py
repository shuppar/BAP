#!/usr/bin/env python
"""h09b_qc_doublets.py -- concat 20 SoupX'd human placenta samples, per-sample QC, doublet removal.

Mirrors mouse Phase 2 (QC) + Phase 3 (doublets) but for the human cohort:
  concat -> per-sample MAD(n_mads=5) + hard floors/caps QC -> scDblFinder (one run,
  samples= preserves per-sample doublet simulation) -> drop doublets -> concatenated h5ad.

Human gene patterns: mito prefix 'MT-' (not mouse 'mt-'); hemoglobin = explicit symbol list.
Erythroid is a REAL placental compartment at term -> pct_hemo flags ambient, does NOT gate cells.

Usage (from project root):
  uv run python scripts/h09b_qc_doublets.py
  uv run python scripts/h09b_qc_doublets.py --sample-ids fs_lean_1 fs_mo_1   # smoke
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp

GSE_DIR = Path("data/human_validation/placenta/gunter_rahman_2025_GSE271976")
RWORKER = "scripts/run_scdblfinder.R"

# QC thresholds (from config/placenta.yaml qc:)
MIN_COUNTS, MIN_GENES = 500, 200
PCT_MT_MAX, PCT_HEMO_MAX, N_MADS = 1.0, 5.0, 5
HEMO_GENES = ["HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBM", "HBQ1", "HBZ"]
MAD_METRICS = ["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_in_top_20_genes"]


def mad_outlier(x, n_mads):
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    if mad == 0:
        return np.zeros(len(x), dtype=bool)
    return (x < med - n_mads * mad) | (x > med + n_mads * mad)


def run_scdblfinder(adata):
    """One scDblFinder call on the QC'd cohort; samples= keeps per-sample boundaries."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        sio.mmwrite(td / "m.mtx", sp.csr_matrix(adata.X).T)  # genes x cells for R
        (td / "bc.tsv").write_text("\n".join(adata.obs_names) + "\n")
        (td / "ft.tsv").write_text("\n".join(adata.var_names) + "\n")
        (td / "sm.tsv").write_text("\n".join(adata.obs["sample_id"].astype(str)) + "\n")
        cmd = ["Rscript", RWORKER, "--matrix", str(td / "m.mtx"),
               "--barcodes", str(td / "bc.tsv"), "--features", str(td / "ft.tsv"),
               "--samples", str(td / "sm.tsv"), "--output", str(td / "out.tsv")]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"scDblFinder failed:\n{r.stderr[-2000:]}")
        d = pd.read_csv(td / "out.tsv", sep="\t").set_index("barcode")
    return d.reindex(adata.obs_names)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-ids", nargs="*", help="subset (smoke test)")
    args = ap.parse_args()

    files = sorted((GSE_DIR / "h5ad").glob("*.h5ad"))
    files = [f for f in files if f.stem != "h09b_qc_doublets"]
    if args.sample_ids:
        files = [f for f in files if f.stem in set(args.sample_ids)]
    if not files:
        sys.exit("no per-sample h5ads found (run h09a first)")
    print(f"[h09b] concatenating {len(files)} samples")
    adata = sc.concat([sc.read_h5ad(f) for f in files], join="outer")
    adata.var_names_make_unique()
    n_in = adata.n_obs

    # --- QC metrics ---
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["hemo"] = adata.var_names.isin(HEMO_GENES)
    if adata.var["mt"].sum() == 0:
        sys.exit("no MT- genes matched -- check human gene symbols")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "hemo"], percent_top=(20,), inplace=True)

    # --- per-sample MAD + hard floors/caps ---
    mad = pd.Series(False, index=adata.obs_names)
    for sid, idx in adata.obs.groupby("sample_id", observed=True).groups.items():
        for m in MAD_METRICS:
            mad.loc[idx] |= mad_outlier(adata.obs.loc[idx, m].values, N_MADS)
    adata.obs["mad_outlier"] = mad.values
    floors = (adata.obs["total_counts"] >= MIN_COUNTS) & (adata.obs["n_genes_by_counts"] >= MIN_GENES)
    caps = (adata.obs["pct_counts_mt"] <= PCT_MT_MAX) & (adata.obs["pct_counts_hemo"] <= PCT_HEMO_MAX)
    keep = (~adata.obs["mad_outlier"].values) & floors.values & caps.values

    print(f"[h09b] QC: {n_in} -> {keep.sum()} "
          f"(MAD drop {mad.sum()}, floor fail {(~floors).sum()}, cap fail {(~caps).sum()})")
    adata = adata[keep].copy()

    # --- doublets ---
    print(f"[h09b] scDblFinder on {adata.n_obs} cells")
    d = run_scdblfinder(adata)
    adata.obs["doublet_score"] = d["doublet_score"].values
    adata.obs["doublet_class"] = d["doublet_class"].values
    n_dbl = (adata.obs["doublet_class"] == "doublet").sum()
    adata = adata[adata.obs["doublet_class"] == "singlet"].copy()
    print(f"[h09b] dropped {n_dbl} doublets -> {adata.n_obs} singlets")

    # --- save + manifest ---
    out = GSE_DIR / "h5ad" / "h09b_qc_doublets.h5ad"
    adata.write(out)
    man = (adata.obs.groupby(["condition", "side", "sample_id"], observed=True)
           .size().rename("n_final").reset_index())
    raw = pd.read_csv(GSE_DIR / "h09a_soupx_manifest.csv")[["sample_id", "n_cells"]]
    man = man.merge(raw.rename(columns={"n_cells": "n_soupx"}), on="sample_id", how="left")
    mpath = GSE_DIR / "h09b_qc_manifest.csv"
    man.to_csv(mpath, index=False)
    print(f"\n[h09b] wrote {out}\n[h09b] manifest -> {mpath}")
    print(man.to_string(index=False))
    print(f"\ntotal: {n_in} raw -> {adata.n_obs} final  (paper ~62,864)")


if __name__ == "__main__":
    main()
