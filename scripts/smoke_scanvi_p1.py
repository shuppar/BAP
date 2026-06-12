#!/usr/bin/env python
"""
smoke_scanvi_p1.py — standalone smoke test for P1 scANVI label transfer.

Validates the Rosenberg → P1 scANVI transfer in isolation BEFORE we wire it
into Phase 7. Checks the two things that decide viability:
  1. Gene overlap between Rosenberg (SPLiT-seq, ~26.9k genes) and our Flex
     query (~19k). scANVI needs enough shared informative genes.
  2. Whether transferred labels are biologically sane — specifically on the
     Leiden clusters Di Bella mislabeled as erythrocyte (2, 3, 51, 45).

Runs on a SUBSAMPLE of P1 (default 20k query cells) so it's fast (~5-10 min on
GPU). This is a throwaway diagnostic — it writes a CSV + prints a report, does
NOT modify the annotated object.

Usage (WS, GPU):
  uv run python scripts/smoke_scanvi_p1.py \
      --query results/brain/h5ad/08_annotated/all_samples.h5ad \
      --ref   refs/rosenberg_p2brain_reference.h5ad \
      --n-query 20000
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", required=True, type=Path)
    ap.add_argument("--ref", required=True, type=Path)
    ap.add_argument("--labels-key", default="rosenberg_fine")
    ap.add_argument("--n-query", type=int, default=20000,
                    help="subsample this many P1 query cells for the smoke test")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--cpu", action="store_true")
    ap.add_argument("--out", type=Path,
                    default=Path("results/brain/tables/smoke_scanvi_p1.csv"))
    args = ap.parse_args()

    import scvi
    import anndata as ad

    print("=== scANVI P1 smoke test ===")

    # --- load query, subset to P1 ---
    print(f"\n[1] loading query {args.query}")
    q = sc.read_h5ad(args.query)
    if "age" not in q.obs.columns:
        sys.exit("ERROR: query has no 'age' column")
    q = q[q.obs["age"] == "P1"].copy()
    print(f"  P1 query cells: {q.n_obs:,}")
    if q.n_obs == 0:
        sys.exit("ERROR: no P1 cells in query")

    # subsample for speed
    if q.n_obs > args.n_query:
        rng = np.random.default_rng(args.seed)
        idx = rng.choice(q.n_obs, args.n_query, replace=False)
        q = q[idx].copy()
        print(f"  subsampled to {q.n_obs:,}")

    # --- load reference ---
    print(f"\n[2] loading reference {args.ref}")
    ref = sc.read_h5ad(args.ref)
    print(f"  ref: {ref.n_obs:,} cells × {ref.n_vars:,} genes, "
          f"{ref.obs[args.labels_key].nunique()} labels")

    # --- GENE OVERLAP (the make-or-break check) ---
    common = q.var_names.intersection(ref.var_names)
    print(f"\n[3] GENE OVERLAP: {len(common):,} shared "
          f"(query {q.n_vars:,}, ref {ref.n_vars:,})")
    if len(common) < 200:
        sys.exit(f"FAIL: only {len(common)} shared genes — naming mismatch?")
    # check key markers survive the intersection
    key_markers = ["Snap25", "Rbfox3", "Aqp4", "Gja1", "Pdgfra", "Olig1",
                   "P2ry12", "Csf1r", "Cldn5", "Mbp", "Gad1", "Slc17a7",
                   "Eomes", "Sox2", "Dcx", "Meg3"]
    present = [g for g in key_markers if g in common]
    missing = [g for g in key_markers if g not in common]
    print(f"  key markers present: {len(present)}/{len(key_markers)}")
    if missing:
        print(f"  key markers MISSING from overlap: {missing}")

    q = q[:, common].copy()
    ref = ref[:, common].copy()

    # --- scANVI ---
    print(f"\n[4] scANVI transfer ...")
    UNLAB = "Unknown"
    ref.obs["_lab"] = ref.obs[args.labels_key].astype(str).values
    q.obs["_lab"] = UNLAB
    ref.obs["_batch"] = "__ref__"
    q.obs["_batch"] = q.obs["pool"].astype(str).values if "pool" in q.obs.columns else "query"

    comb = ad.concat([ref, q], axis=0, label="_origin", keys=["ref", "query"],
                     index_unique="-")
    if comb.X is None:
        sys.exit("ERROR: combined .X is None (need raw counts)")

    scvi.settings.seed = args.seed
    acc = "cpu" if args.cpu else "gpu"
    prec = "32" if args.cpu else "bf16-mixed"

    scvi.model.SCVI.setup_anndata(comb, batch_key="_batch")
    vae = scvi.model.SCVI(comb, n_layers=2, n_latent=30)
    print("  training scVI ...")
    vae.train(max_epochs=200, accelerator=acc, devices=1, precision=prec,
              early_stopping=True)

    print("  training scANVI ...")
    lvae = scvi.model.SCANVI.from_scvi_model(
        vae, adata=comb, unlabeled_category=UNLAB, labels_key="_lab")
    lvae.train(max_epochs=100, n_samples_per_label=100,
               accelerator=acc, devices=1, precision=prec)

    comb.obs["_pred"] = lvae.predict(comb)
    proba = lvae.predict(comb, soft=True)
    comb.obs["_conf"] = proba.max(axis=1).values

    isq = comb.obs["_origin"] == "query"
    qpred = comb[isq].copy()
    qpred.obs_names = [bc.rsplit("-query", 1)[0] for bc in qpred.obs_names]

    # --- sanity report ---
    print(f"\n[5] RESULTS")
    print(f"  mean confidence: {qpred.obs['_conf'].mean():.3f} "
          f"(median {qpred.obs['_conf'].median():.3f})")
    print(f"\n  transferred label distribution (top 25):")
    print(qpred.obs["_pred"].value_counts().head(25).to_string())

    # erythrocyte check: ZERO Rosenberg labels are erythroid, so any cell that
    # Di Bella called erythrocyte should now get a real neural/glial label
    print(f"\n  --- the Di-Bella-erythrocyte clusters (2,3,51,45) ---")
    if "leiden" in qpred.obs.columns:
        for cl in ["2", "3", "51", "45"]:
            sub = qpred.obs[qpred.obs["leiden"] == cl]
            if len(sub) == 0:
                continue
            top = sub["_pred"].value_counts().head(3)
            print(f"  cluster {cl} (n={len(sub)}):")
            for lab, n in top.items():
                print(f"      {n:5d}  {lab}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_cols = qpred.obs[["_pred", "_conf"]].copy()
    if "leiden" in qpred.obs.columns:
        out_cols["leiden"] = qpred.obs["leiden"].values
    out_cols.to_csv(args.out)
    print(f"\n  wrote {args.out}")
    print(f"\n=== smoke test done — inspect labels above ===")


if __name__ == "__main__":
    main()
