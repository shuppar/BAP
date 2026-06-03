#!/usr/bin/env python
"""
dev_split_h5.py — DEV-ONLY pre-processing. NOT part of the workstation pipeline.

Reads the dev subset's 10x Flex .h5 files and writes, for each one, N new
10x-format .h5 files whose cells are a random partition of the original. Then
emits a dev config (config/dev_split.yaml) listing the N×(n_samples) pseudo-
samples, each with its own id / donor_id (suffixed _ps1.._psN) but inheriting
the parent's group / age / sex / pool.

Why: dev has one real sample per group, so pseudobulk DE/composition (8a/8b/8c)
has no replicates. Splitting each sample into N pseudo-donors gives the
downstream stats n=N per group so those code paths run. The numbers are
MEANINGLESS (random cell partitions of one animal — no real between-animal
variance); this is a smoke test of the code paths only.

This runs ONCE, before the pipeline, and changes NOTHING in the phase scripts.
The pipeline then runs normally:  --config config/dev_split.yaml

Usage:
    uv run python dev_split_h5.py --config config/dev.yaml --n 3
    uv run python scripts/01_validate.py --config config/dev_split.yaml
    # ... then 02_qc.py etc. with config/dev_split.yaml

The written .h5 files are real 10x v3 feature-barcode HDF5 matrices, so
sc.read_10x_h5 (used by 02_qc) reads them with no code change. Each written
file is verified by reading it back and checking shape + total counts.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import scanpy as sc
import scipy.sparse as sp
import h5py
import yaml

sys.path.insert(0, "scripts")
from _utils import load_config


def write_10x_h5(adata, out_path: Path, genome: str = "mm10"):
    """Write an AnnData (cells × genes, integer counts in .X) as a 10x v3
    feature-barcode HDF5 matrix that sc.read_10x_h5 can read back.

    Schema (per 10x docs): /matrix/{barcodes,data,indices,indptr,shape} +
    /matrix/features/{_all_tag_keys,feature_type,genome,id,name}. Matrix stored
    CSC with shape [n_features, n_barcodes].
    """
    # CSC, genes × cells
    X = adata.X
    X = sp.csc_matrix(X.T) if not (sp.issparse(X) and X.format == "csc") else X.T.tocsc()
    X.eliminate_zeros()
    data = np.asarray(X.data, dtype=np.int32)
    indices = np.asarray(X.indices, dtype=np.int64)
    indptr = np.asarray(X.indptr, dtype=np.int64)
    n_genes, n_cells = X.shape

    # Feature ids/names: prefer Ensembl in var['gene_ids'] if present
    ids = (adata.var["gene_ids"].astype(str).values
           if "gene_ids" in adata.var.columns else adata.var_names.astype(str).values)
    names = adata.var_names.astype(str).values
    ftype = (adata.var["feature_types"].astype(str).values
             if "feature_types" in adata.var.columns
             else np.array(["Gene Expression"] * n_genes))
    barcodes = adata.obs_names.astype(str).values

    def _b(arr):  # encode to ascii bytes for h5py fixed-length strings
        return np.array([str(x).encode("ascii", "replace") for x in arr])

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out_path, "w") as f:
        m = f.create_group("matrix")
        m.create_dataset("barcodes", data=_b(barcodes))
        m.create_dataset("data", data=data)
        m.create_dataset("indices", data=indices)
        m.create_dataset("indptr", data=indptr)
        m.create_dataset("shape", data=np.array([n_genes, n_cells], dtype=np.int64))
        ft = m.create_group("features")
        ft.create_dataset("_all_tag_keys", data=_b(["genome"]))
        ft.create_dataset("feature_type", data=_b(ftype))
        ft.create_dataset("genome", data=_b([genome] * n_genes))
        ft.create_dataset("id", data=_b(ids))
        ft.create_dataset("name", data=_b(names))


def verify(out_path: Path, expect_cells: int, expect_total: int):
    """Read back with sc.read_10x_h5 and check shape + total counts match."""
    a = sc.read_10x_h5(out_path)
    got_cells = a.n_obs
    got_total = int(np.asarray(a.X.sum()))
    if got_cells != expect_cells:
        sys.exit(f"VERIFY FAIL {out_path.name}: {got_cells} cells, expected {expect_cells}")
    if got_total != expect_total:
        sys.exit(f"VERIFY FAIL {out_path.name}: total {got_total}, expected {expect_total}")
    print(f"    verified {out_path.name}: {got_cells} cells, total counts {got_total}")


def main():
    ap = argparse.ArgumentParser(description="DEV-ONLY: split dev h5 files into N pseudo-donors")
    ap.add_argument("--config", required=True, type=Path,
                    help="dev config to read the sample list from (e.g. config/dev.yaml)")
    ap.add_argument("--n", type=int, default=3, help="pseudo-donors per sample (default 3)")
    ap.add_argument("--outdir", type=Path, default=Path("data/dev_split"),
                    help="where to write split .h5 files (default data/dev_split)")
    ap.add_argument("--out-config", type=Path, default=Path("config/dev_split.yaml"),
                    help="dev config to write (default config/dev_split.yaml)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    if args.n < 2:
        sys.exit("--n must be >= 2")

    cfg = load_config(args.config)
    samples = cfg["samples"]
    print(f"=== DEV h5 split: {len(samples)} samples × {args.n} pseudo-donors ===")
    print(f"  (MEANINGLESS numbers — smoke test of stats code paths only)\n")

    rng = np.random.default_rng(args.seed)
    new_samples = []

    # Load the original raw config text to inherit top-level blocks verbatim.
    with open(args.config) as f:
        raw = yaml.safe_load(f)

    for s in samples:
        sid = s["id"]
        h5 = s["h5"]
        print(f"  [{sid}] reading {Path(h5).name}")
        adata = sc.read_10x_h5(h5)
        adata.var_names_make_unique()
        n = adata.n_obs
        parts = rng.integers(0, args.n, size=n)

        for k in range(args.n):
            sub = adata[parts == k].copy()
            new_id = f"{sid}_ps{k+1}"
            out_h5 = args.outdir / f"{new_id}.h5"
            write_10x_h5(sub, out_h5)
            verify(out_h5, sub.n_obs, int(np.asarray(sub.X.sum())))

            rec = dict(s)                      # inherit group/age/sex/pool/etc.
            rec["id"] = new_id
            rec["donor_id"] = f"{s.get('donor_id', sid)}_ps{k+1}"
            rec["h5"] = str(out_h5)
            rec.pop("raw_h5", None)            # no split raw matrix (CellBender skipped in dev)
            new_samples.append(rec)
        print(f"  [{sid}] -> {args.n} files, {n} cells split\n")

    # Build the split config: start from raw dev.yaml (so local edits win),
    # then ADD any blocks that load_config inherited from samples_from (brain.yaml)
    # but that aren't literally present in dev.yaml — contrasts, pathways,
    # composition, stress_focused_cell_types, reference. Without this the
    # dev_split config would lack those blocks and Phase 8 would fail with
    # "no 'contrasts:' block in config."
    out_cfg = {k: v for k, v in raw.items()
               if k not in ("samples_from", "subset", "samples")}
    INHERITED = ("contrasts", "stress_focused_cell_types", "composition",
                 "pathways", "reference", "annotation")
    for key in INHERITED:
        if key not in out_cfg and key in cfg:
            out_cfg[key] = cfg[key]
    out_cfg["samples"] = new_samples
    # Keep the dev cell cap if it was present, applied per pseudo-sample.
    cap = raw.get("subset", {}).get("max_cells_per_sample")
    if cap:
        out_cfg["subset"] = {"enabled": False, "max_cells_per_sample": cap}

    args.out_config.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_config, "w") as f:
        f.write("# GENERATED by dev_split_h5.py — DEV-ONLY, do not hand-edit.\n")
        f.write("# Pseudo-donor split for smoke-testing 8a/8b/8c. Numbers meaningless.\n\n")
        yaml.safe_dump(out_cfg, f, sort_keys=False)

    print(f"✓ Wrote {len(new_samples)} split samples to {args.outdir}/")
    print(f"✓ Wrote config: {args.out_config}")
    print(f"\nNext:\n  uv run python scripts/01_validate.py --config {args.out_config}")


if __name__ == "__main__":
    main()
