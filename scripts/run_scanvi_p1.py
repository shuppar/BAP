#!/usr/bin/env python
"""
run_scanvi_p1.py — GPU subprocess worker: transfer Rosenberg labels onto P1 brain.

Called by 07_annotation.py (brain) as a subprocess, mirroring how 03_doublets.py
calls run_scdblfinder.R. Phase 7 stays CPU/in-process for 4W/3mo (ABC CellTypist)
and STAMP (placenta); only the P1 brain branch shells out to this GPU script.

Contract (the "API" between Phase 7 and this worker):
  INPUT  --query-h5ad : a TEMP h5ad containing ONLY P1 brain cells, raw counts
                        in .X, var_names = gene symbols, obs has 'pool'.
         --ref-h5ad   : refs/rosenberg_p2brain_reference.h5ad (raw counts,
                        obs['rosenberg_fine']).
  OUTPUT --out-tsv    : TSV [barcode, scanvi_fine, scanvi_conf], one row per
                        P1 query cell. Phase 7 reads this back and joins on
                        barcode, then derives class/region/broad from
                        scanvi_fine via the config CSVs.

Isolated as a subprocess because scANVI is a GPU torch job; keeping it out of
the main annotation process avoids loading torch/scvi into Phase 7 and lets a
crash here surface cleanly without killing the whole annotation run.

Usage (normally invoked by Phase 7, but runnable standalone for debugging):
  uv run python scripts/run_scanvi_p1.py \
      --query-h5ad /tmp/p1_query.h5ad \
      --ref-h5ad   refs/rosenberg_p2brain_reference.h5ad \
      --out-tsv    /tmp/p1_scanvi.tsv \
      --labels-key rosenberg_fine \
      --seed 42
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query-h5ad", required=True, type=Path)
    ap.add_argument("--ref-h5ad", required=True, type=Path)
    ap.add_argument("--out-tsv", required=True, type=Path)
    ap.add_argument("--labels-key", default="rosenberg_fine")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--max-epochs-scvi", type=int, default=200)
    ap.add_argument("--max-epochs-scanvi", type=int, default=100)
    ap.add_argument("--min-shared-genes", type=int, default=2000)
    ap.add_argument("--cpu", action="store_true", help="force CPU (debug only)")
    args = ap.parse_args()

    import torch
    import anndata as ad
    import scvi

    # perf: Tensor-Core matmul hint (RTX 4500 Ada). 'high' = TF32, safe for scVI.
    if not args.cpu:
        torch.set_float32_matmul_precision("high")

    for p in (args.query_h5ad, args.ref_h5ad):
        if not p.is_file():
            sys.exit(f"ERROR: input not found: {p}")

    UNLAB = "Unknown"
    REF_BATCH = "__reference__"

    print(f"[scanvi_p1] loading query {args.query_h5ad}")
    q = sc.read_h5ad(args.query_h5ad)
    print(f"  query: {q.n_obs:,} cells × {q.n_vars:,} genes")
    if q.n_obs == 0:
        sys.exit("ERROR: query has 0 cells")
    if "pool" not in q.obs.columns:
        sys.exit("ERROR: query missing 'pool' column (needed as query batch)")

    print(f"[scanvi_p1] loading reference {args.ref_h5ad}")
    ref = sc.read_h5ad(args.ref_h5ad)
    if args.labels_key not in ref.obs.columns:
        sys.exit(f"ERROR: labels_key '{args.labels_key}' not in reference obs. "
                 f"Have: {list(ref.obs.columns)}")
    print(f"  ref: {ref.n_obs:,} cells × {ref.n_vars:,} genes, "
          f"{ref.obs[args.labels_key].nunique()} labels")

    # --- gene intersection ---
    common = q.var_names.intersection(ref.var_names)
    print(f"[scanvi_p1] shared genes: {len(common):,}")
    if len(common) < args.min_shared_genes:
        sys.exit(f"ERROR: only {len(common)} shared genes "
                 f"(< {args.min_shared_genes}). Gene-naming mismatch?")
    q = q[:, common].copy()
    ref = ref[:, common].copy()

    # --- label + batch columns ---
    ref.obs["_lab"] = ref.obs[args.labels_key].astype(str).values
    q.obs["_lab"] = UNLAB
    ref.obs["_batch"] = REF_BATCH
    q.obs["_batch"] = q.obs["pool"].astype(str).values

    # --- concat (raw counts must be in .X) ---
    combined = ad.concat([ref, q], axis=0, label="_origin",
                         keys=["ref", "query"], index_unique="-")
    if combined.X is None:
        sys.exit("ERROR: combined .X is None — need raw counts in .X")
    print(f"[scanvi_p1] combined: {combined.n_obs:,} "
          f"({ref.n_obs:,} ref + {q.n_obs:,} query)")

    scvi.settings.seed = args.seed
    accel = "cpu" if args.cpu else "gpu"
    prec = "32" if args.cpu else "bf16-mixed"

    # --- scVI: integrate ref+query, correct platform/batch shift ---
    scvi.model.SCVI.setup_anndata(combined, batch_key="_batch")
    vae = scvi.model.SCVI(combined, n_layers=2, n_latent=30)
    print(f"[scanvi_p1] training scVI (max_epochs={args.max_epochs_scvi}) ...")
    vae.train(
        max_epochs=args.max_epochs_scvi,
        accelerator=accel, devices=1, precision=prec,
        early_stopping=True,
        datasplitter_kwargs={"num_workers": 4, "pin_memory": True},
    )

    # --- scANVI: learn labels in shared space, predict query ---
    print(f"[scanvi_p1] training scANVI (max_epochs={args.max_epochs_scanvi}) ...")
    lvae = scvi.model.SCANVI.from_scvi_model(
        vae, adata=combined, unlabeled_category=UNLAB, labels_key="_lab")
    lvae.train(
        max_epochs=args.max_epochs_scanvi, n_samples_per_label=100,
        accelerator=accel, devices=1, precision=prec,
        datasplitter_kwargs={"num_workers": 4, "pin_memory": True},
    )

    # --- predict query cells only ---
    pred = lvae.predict(combined)
    proba = lvae.predict(combined, soft=True)
    conf = proba.max(axis=1).values

    combined.obs["_pred"] = pred
    combined.obs["_conf"] = conf
    isq = (combined.obs["_origin"] == "query").values
    qpred = combined[isq].copy()
    # strip concat suffix to restore original query barcodes
    qpred.obs_names = [bc.rsplit("-query", 1)[0] for bc in qpred.obs_names]

    out = pd.DataFrame({
        "scanvi_fine": qpred.obs["_pred"].astype(str).values,
        "scanvi_conf": qpred.obs["_conf"].astype(float).values,
    }, index=qpred.obs_names)
    out.index.name = "barcode"

    # safety: one row per query cell, no NaN
    if out["scanvi_fine"].isna().any():
        sys.exit("ERROR: some query cells got no prediction")
    if len(out) != int(isq.sum()):
        sys.exit(f"ERROR: row count {len(out)} != query cells {int(isq.sum())}")

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t")
    print(f"[scanvi_p1] wrote {args.out_tsv}  ({len(out):,} cells)")
    print(f"[scanvi_p1] mean conf {out['scanvi_conf'].mean():.3f}, "
          f"median {out['scanvi_conf'].median():.3f}")


if __name__ == "__main__":
    main()
