#!/usr/bin/env python
"""
02_soupx.py -- Phase 1: ambient RNA correction via SoupX (R subprocess).

Replaces the abandoned CellBender phase. For each sample in the config,
invokes scripts/run_soupx.R with cellranger filtered + raw counts, reads
the corrected MTX back, assembles an AnnData with sample metadata in obs,
writes to results/<tissue>/h5ad/02_soupx_corrected/<sample_id>.h5ad.

The raw counts referenced by raw_h5 may be:
  - a .h5 file (used directly)
  - a .tar.gz archive (extracted to a per-sample temp dir, MTX dir inside)

Per-sample is independent, so we parallelize with ProcessPoolExecutor.
Each R subprocess holds ~5-15 GB RAM (raw matrix). Default --n-jobs=4 is
safe on the 258 GB workstation; bump to 6-8 for production.

Outputs:
  results/<tissue>/h5ad/02_soupx_corrected/<sample_id>.h5ad
  results/<tissue>/tables/02_soupx/02_soupx_summary.csv
  results/<tissue>/logs/02_soupx/<sample_id>.log

Idempotent: re-running skips samples whose h5ad output already exists
(use --force to override).

Usage:
  # smoke test on one sample (sequential, fast log feedback)
  uv run python scripts/02_soupx.py --config config/brain.yaml \
      --sample-ids E1 --n-jobs 1

  # production
  uv run python scripts/02_soupx.py --config config/brain.yaml --n-jobs 6
  uv run python scripts/02_soupx.py --config config/placenta.yaml --n-jobs 6
"""

import argparse
import json
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io as scio
from scipy import sparse

from _utils import load_config


REPO_ROOT = Path(__file__).resolve().parent.parent
R_SCRIPT  = REPO_ROOT / "scripts" / "run_soupx.R"


def resolve_raw_path(raw_h5_path: str, tmp_dir: Path) -> str:
    """Return a path usable by R's read10xCounts.

    If raw_h5_path is .h5 -> return as-is.
    If raw_h5_path is .tar.gz -> extract into tmp_dir and return the
    directory containing matrix.mtx[.gz] / barcodes.tsv[.gz] / features.tsv[.gz].
    """
    raw_path = Path(raw_h5_path)
    if not raw_path.exists():
        raise FileNotFoundError(f"raw_h5 not found: {raw_path}")

    if raw_path.suffix == ".h5":
        return str(raw_path)

    if raw_path.name.endswith(".tar.gz") or raw_path.suffix in (".tar", ".gz", ".tgz"):
        extract_dir = tmp_dir / "raw_extracted"
        extract_dir.mkdir(exist_ok=True)
        with tarfile.open(raw_path, "r:*") as tar:
            # python 3.12+: filter='data' avoids tarfile security warnings
            try:
                tar.extractall(extract_dir, filter="data")
            except TypeError:
                tar.extractall(extract_dir)
        # Find a directory that looks like a 10X MTX dir. Some Cell Ranger
        # tar.gz archives have the three files at the archive root (no
        # enclosing directory), others put them in a subdir -- check both.
        candidates = [extract_dir] + [d for d in extract_dir.rglob("*") if d.is_dir()]
        for d in candidates:
            names = {p.name for p in d.iterdir() if p.is_file()}
            if {"matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"}.issubset(names):
                return str(d)
            if {"matrix.mtx", "barcodes.tsv", "features.tsv"}.issubset(names):
                return str(d)
        raise FileNotFoundError(
            f"After extracting {raw_path}, no 10X MTX directory found. "
            f"Looked in {extract_dir}"
        )

    raise ValueError(f"raw_h5 must be .h5 or .tar.gz, got: {raw_path}")


def run_soupx_one(sample: dict, out_h5ad: Path, log_path: Path,
                  rho: float | None = None, force: bool = False) -> dict:
    """Run SoupX on one sample, return a summary row dict.

    On success: row contains sample_id, rho_mean, pct_removed, n_cells, status.
    On skip:    row contains status='cached'.
    Raises on any failure (caught by main).
    """
    sample_id = sample["id"]
    if out_h5ad.exists() and not force:
        return {"sample_id": sample_id, "status": "cached"}

    tmp_dir = Path(tempfile.mkdtemp(prefix=f"soupx_{sample_id}_"))
    r_output_dir = tmp_dir / "r_output"
    r_output_dir.mkdir()

    try:
        # Resolve raw path (extract tar.gz if needed)
        raw_resolved = resolve_raw_path(sample["raw_h5"], tmp_dir)

        # Build R command
        cmd = [
            "Rscript", str(R_SCRIPT),
            "--filtered",   sample["h5"],
            "--raw",        raw_resolved,
            "--output_dir", str(r_output_dir),
            "--sample_id",  sample_id,
        ]
        if rho is not None:
            cmd.extend(["--rho", str(rho)])

        # Run, streaming stdout/stderr to per-sample log
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as logf:
            logf.write(f"CMD: {' '.join(cmd)}\n\n")
            logf.flush()
            result = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            raise RuntimeError(
                f"R subprocess for {sample_id} exited with code "
                f"{result.returncode}. See log: {log_path}"
            )

        # Read R outputs
        matrix_path   = r_output_dir / "matrix.mtx"
        barcodes_path = r_output_dir / "barcodes.tsv"
        features_path = r_output_dir / "features.tsv"
        summary_path  = r_output_dir / "soupx_summary.json"
        for p in (matrix_path, barcodes_path, features_path, summary_path):
            if not p.exists():
                raise FileNotFoundError(f"Expected R output missing: {p}")

        # Matrix Market is genes x cells; transpose to cells x genes for AnnData
        mtx = scio.mmread(str(matrix_path)).tocsr()
        mtx = mtx.T.tocsr()

        barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")[0].astype(str).tolist()
        features = pd.read_csv(features_path, sep="\t", dtype=str)

        if "ID" not in features.columns:
            raise ValueError(
                f"features.tsv from {sample_id} has no 'ID' column. "
                f"Columns found: {list(features.columns)}"
            )

        var = features.copy()
        var = var.set_index("ID")
        var.index.name = None
        # Standardise the symbol column name to 'symbol' (lowercase) for
        # consistency with the rest of the pipeline (project doc: var['symbol']).
        if "Symbol" in var.columns and "symbol" not in var.columns:
            var = var.rename(columns={"Symbol": "symbol"})

        obs = pd.DataFrame(index=pd.Index(barcodes, name=None))
        # Attach sample metadata to every cell (downstream phases use these)
        for key in ("id", "donor_id", "age", "group", "sex", "pool", "library"):
            if key in sample:
                # Rename 'id' to 'sample_id' in obs to match downstream convention
                obs_key = "sample_id" if key == "id" else key
                obs[obs_key] = sample[key]

        adata = ad.AnnData(X=mtx, obs=obs, var=var)

        # Pull the R summary into uns for traceability
        with open(summary_path) as f:
            soupx_summary = json.load(f)
        adata.uns["soupx"] = soupx_summary

        # Write
        out_h5ad.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(out_h5ad)

        return {
            "sample_id":   sample_id,
            "status":      "ok",
            "rho_mean":    soupx_summary["rho_mean"],
            "rho_min":     soupx_summary["rho_min"],
            "rho_max":     soupx_summary["rho_max"],
            "n_cells":     soupx_summary["n_cells"],
            "pct_removed": soupx_summary["pct_removed"],
            "n_clusters":  soupx_summary["n_clusters"],
            "elapsed_sec": soupx_summary["elapsed_sec"],
            "mode":        soupx_summary["mode"],
        }
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True,
                        help="Path to brain.yaml / placenta.yaml / dev.yaml")
    parser.add_argument("--n-jobs", type=int, default=4,
                        help="Parallel R subprocesses (default 4)")
    parser.add_argument("--rho", type=float, default=None,
                        help="Manual contamination fraction; overrides autoEst")
    parser.add_argument("--force", action="store_true",
                        help="Re-run even if output h5ad exists")
    parser.add_argument("--sample-ids", nargs="+", default=None,
                        help="Only process listed sample ids (smoke testing)")
    args = parser.parse_args()

    if not R_SCRIPT.exists():
        sys.exit(f"ERROR: R script not found at {R_SCRIPT}")

    cfg = load_config(args.config)
    tissue = cfg["tissue"]
    results_dir = Path(cfg.get("results_dir", f"results/{tissue}"))
    out_dir  = results_dir / "h5ad"  / "02_soupx_corrected"
    log_dir  = results_dir / "logs"  / "02_soupx"
    tbl_dir  = results_dir / "tables" / "02_soupx"
    for d in (out_dir, log_dir, tbl_dir):
        d.mkdir(parents=True, exist_ok=True)

    samples = cfg["samples"]
    if args.sample_ids:
        keep = set(args.sample_ids)
        samples = [s for s in samples if s["id"] in keep]
        missing = keep - {s["id"] for s in samples}
        if missing:
            sys.exit(f"ERROR: --sample-ids not found in config: {sorted(missing)}")

    print(f"[soupx] tissue={tissue}  samples={len(samples)}  n_jobs={args.n_jobs}")
    print(f"[soupx] R script:    {R_SCRIPT}")
    print(f"[soupx] output h5ad: {out_dir}")
    print(f"[soupx] logs:        {log_dir}")
    if args.rho is not None:
        print(f"[soupx] MANUAL rho = {args.rho}")
    else:
        print(f"[soupx] autoEst via scran::quickCluster + SoupX::autoEstCont")

    t_start = time.time()
    rows = []

    def submit_and_collect():
        with ProcessPoolExecutor(max_workers=args.n_jobs) as ex:
            futures = {
                ex.submit(
                    run_soupx_one,
                    s,
                    out_dir / f"{s['id']}.h5ad",
                    log_dir / f"{s['id']}.log",
                    args.rho,
                    args.force,
                ): s for s in samples
            }
            for fut in as_completed(futures):
                s = futures[fut]
                try:
                    row = fut.result()
                except Exception as e:
                    print(f"  [FAIL] {s['id']}: {e}")
                    rows.append({"sample_id": s["id"], "status": "failed",
                                 "error": str(e)})
                    continue
                if row["status"] == "cached":
                    print(f"  [skip] {row['sample_id']}: cached "
                          f"(use --force to re-run)")
                else:
                    print(f"  [done] {row['sample_id']}: "
                          f"rho={row['rho_mean']:.4f}  "
                          f"removed={row['pct_removed']:.2f}%  "
                          f"cells={row['n_cells']:,}  "
                          f"time={row['elapsed_sec']:.1f}s")
                rows.append(row)

    submit_and_collect()

    # Summary CSV
    df = pd.DataFrame(rows)
    summary_csv = tbl_dir / "02_soupx_summary.csv"
    df.to_csv(summary_csv, index=False)

    elapsed = (time.time() - t_start) / 60
    n_ok     = int((df["status"] == "ok").sum())     if len(df) else 0
    n_cache  = int((df["status"] == "cached").sum()) if len(df) else 0
    n_fail   = int((df["status"] == "failed").sum()) if len(df) else 0
    print(f"\n[soupx] Done. {n_ok} corrected, {n_cache} cached, "
          f"{n_fail} failed.  Wall {elapsed:.1f} min.")
    print(f"[soupx] Summary CSV: {summary_csv}")
    if n_fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
