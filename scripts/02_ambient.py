#!/usr/bin/env python
"""
02_ambient.py — Phase 1: Ambient RNA correction via CellBender.

Run from the CellBender sidecar venv, NOT the main uv venv:
    CUDA_VISIBLE_DEVICES=0 .venv-cellbender/bin/python scripts/02_ambient.py \\
        --config config/brain.yaml

What it does, per sample:
  1. Locates the raw h5 (from raw_h5 field in config — may be a .tar.gz, a
     directory, or a direct .h5 path). Extracts if needed.
  2. Runs CellBender remove-background (GPU, --epochs 150,
     --cells-posterior-reg 50, --cuda).
  3. Writes per-sample CellBender output to
       results/{tissue}/h5ad/01_cellbender/{sample_id}/
  4. Copies the _filtered.h5 output (CellBender's ambient-corrected cell
     matrix) alongside a barcode-rank plot and ambient-fraction summary.

Idempotent: samples whose _filtered.h5 output already exists are skipped.

GPU serialisation: CellBender uses the full GPU per sample. This script runs
samples sequentially (not parallel) by default. Pass --parallel N to run N
at once — only safe if VRAM per sample < 24 GB / N. With 34 brain samples at
typical snRNA-seq size, sequential is fine (~1-2h/sample → ~2 days total, but
see --parallel 2 option which is validated for the RTX 4500 Ada 24 GB).

Outputs:
  results/{tissue}/h5ad/01_cellbender/{sample_id}/
      {sample_id}_filtered.h5        — ambient-corrected cells (input to Phase 2)
      {sample_id}_cell_barcodes.csv  — barcodes called as cells by CellBender
      {sample_id}.log                — full CellBender stdout/stderr
  results/{tissue}/plots/01_ambient/
      {sample_id}_barcodes.png       — CellBender barcode rank curve
      summary_ambient_fraction.csv   — per-sample ambient fraction summary

Phase 2 (02_qc.py) reads from h5ad/01_cellbender/{sample_id}/{sample_id}_filtered.h5.
If you skip Phase 1 (not recommended), Phase 2 will fall back to the raw
filtered h5 from the config.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

# NOTE: this script runs under .venv-cellbender, which has cellbender + torch
# but NOT scanpy/anndata. Keep imports minimal.

# ============================================================================
# Helpers
# ============================================================================

def log(msg: str) -> None:
    print(f"[ambient] {msg}", flush=True)


def warn(msg: str) -> None:
    print(f"[ambient][WARN] {msg}", flush=True, file=sys.stderr)


def load_config_minimal(path: Path) -> dict:
    """Minimal YAML loader — avoids importing _utils (which needs scanpy)."""
    import yaml
    path = Path(path)
    with path.open() as f:
        cfg = yaml.safe_load(f)

    # Resolve relative paths to absolute (same logic as _utils.load_config)
    cwd = Path.cwd()
    for s in cfg.get("samples", []):
        for key in ("h5", "raw_h5"):
            if key in s and s[key]:
                p = Path(s[key])
                if not p.is_absolute():
                    s[key] = str((cwd / p).resolve())
    if "results_dir" in cfg:
        rp = Path(cfg["results_dir"])
        if not rp.is_absolute():
            cfg["results_dir"] = str((cwd / rp).resolve())
    return cfg


def resolve_raw_h5(raw_h5_field: str, sample_id: str, extract_dir: Path) -> Path:
    """Return path to a usable raw .h5 file.

    raw_h5 in the config can be:
      - a direct .h5 path   → return as-is
      - a .tar.gz           → extract to extract_dir/{sample_id}/, return the
                              sample_raw_feature_bc_matrix.h5 inside
      - a directory         → return directory/sample_raw_feature_bc_matrix.h5

    Raises FileNotFoundError if nothing usable is found.
    """
    p = Path(raw_h5_field)

    if p.suffix == ".h5" and p.exists():
        return p

    if p.name.endswith(".tar.gz") or p.name.endswith(".tgz"):
        dest = extract_dir / sample_id
        dest.mkdir(parents=True, exist_ok=True)
        # Check if already extracted
        candidate = dest / "sample_raw_feature_bc_matrix.h5"
        if candidate.exists():
            log(f"  [{sample_id}] raw h5 already extracted: {candidate}")
            return candidate
        log(f"  [{sample_id}] extracting {p.name} → {dest}/")
        if not p.exists():
            raise FileNotFoundError(f"raw_h5 tar.gz not found: {p}")
        with tarfile.open(p, "r:gz") as tf:
            tf.extractall(dest)
        # After extraction, find the .h5 file
        candidates = list(dest.rglob("sample_raw_feature_bc_matrix.h5"))
        if not candidates:
            candidates = list(dest.rglob("*.h5"))
        if not candidates:
            raise FileNotFoundError(
                f"No .h5 found after extracting {p} into {dest}"
            )
        return candidates[0]

    if p.is_dir():
        candidate = p / "sample_raw_feature_bc_matrix.h5"
        if candidate.exists():
            return candidate
        raise FileNotFoundError(f"raw_h5 dir exists but no h5 inside: {p}")

    raise FileNotFoundError(f"raw_h5 not found / unrecognised format: {p}")


def count_expected_cells(filtered_h5: str) -> int:
    """Estimate expected_cells from the Cell Ranger filtered matrix barcode count.

    CellBender's --expected-cells should be set to the Cell Ranger filtered
    cell count. We read the filtered h5 barcodes without scanpy (h5py only).
    Falls back to 5000 if anything goes wrong.
    """
    try:
        import h5py
        with h5py.File(filtered_h5, "r") as f:
            # Standard Cell Ranger v3 h5 layout
            for key in ("matrix/barcodes", "barcodes"):
                if key in f:
                    return len(f[key])
    except Exception as e:
        warn(f"Could not read expected_cells from {filtered_h5}: {e}")
    return 5000


def run_cellbender(
    sample_id: str,
    raw_h5: Path,
    expected_cells: int,
    out_dir: Path,
    epochs: int = 150,
    cuda: bool = True,
    fpr: float = 0.01,
) -> Path:
    """Run CellBender remove-background for one sample. Returns path to _filtered.h5."""
    out_dir.mkdir(parents=True, exist_ok=True)
    out_h5 = out_dir / f"{sample_id}.h5"
    filtered_h5 = out_dir / f"{sample_id}_filtered.h5"
    log_file = out_dir / f"{sample_id}.log"

    if filtered_h5.exists():
        log(f"  [{sample_id}] _filtered.h5 already exists — skipping CellBender.")
        return filtered_h5

    cmd = [
        str(Path(sys.executable).parent / "cellbender"), "remove-background",
        "--input", str(raw_h5),
        "--output", str(out_h5),
        "--expected-cells", str(expected_cells),
        "--epochs", str(epochs),
        "--posterior-regularization", "PRmu",
        "--fpr", str(fpr),
    ]
    if cuda:
        cmd.append("--cuda")

    log(f"  [{sample_id}] running CellBender (expected_cells={expected_cells}, "
        f"epochs={epochs}, cuda={cuda})")
    log(f"  [{sample_id}] command: {' '.join(cmd)}")
    log(f"  [{sample_id}] log: {log_file}")

    with open(log_file, "w") as lf:
        result = subprocess.run(
            cmd,
            stdout=lf,
            stderr=subprocess.STDOUT,
            cwd=str(out_dir),
        )

    if result.returncode != 0:
        # Print tail of log to help diagnose
        with open(log_file) as lf:
            lines = lf.readlines()
        tail = "".join(lines[-40:])
        warn(f"  [{sample_id}] CellBender FAILED (exit {result.returncode}). "
             f"Last 40 lines of log:\n{tail}")
        return None

    if not filtered_h5.exists():
        warn(f"  [{sample_id}] CellBender exited 0 but _filtered.h5 not found at {filtered_h5}")
        return None

    log(f"  [{sample_id}] done → {filtered_h5}")
    return filtered_h5


def collect_summary(samples_done: list[dict], out_csv: Path) -> None:
    """Write a per-sample summary CSV of ambient fraction estimates."""
    import csv
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "sample_id", "expected_cells", "status", "filtered_h5"
        ])
        writer.writeheader()
        writer.writerows(samples_done)
    log(f"Summary written → {out_csv}")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Phase 1: CellBender ambient RNA correction"
    )
    parser.add_argument("--config", required=True, type=Path,
                        help="Path to tissue YAML config (brain.yaml / placenta.yaml)")
    parser.add_argument("--epochs", type=int, default=150,
                        help="CellBender epochs (default 150)")
    parser.add_argument("--fpr", type=float, default=0.01,
                        help="CellBender false-positive rate (default 0.01)")
    parser.add_argument("--no-cuda", action="store_true",
                        help="Disable CUDA (CPU only, much slower)")
    parser.add_argument("--parallel", type=int, default=1,
                        help="Number of samples to run in parallel (default 1; "
                             "use 2 only if VRAM comfortably fits 2 samples)")
    parser.add_argument("--sample", type=str, default=None,
                        help="Run only this sample_id (for re-runs / debugging)")
    args = parser.parse_args()

    log(f"=== Phase 1: CellBender ambient RNA correction ===")
    log(f"Config: {args.config}")

    cfg = load_config_minimal(args.config)
    tissue = cfg.get("tissue", "unknown")
    samples = cfg["samples"]
    results_dir = Path(cfg["results_dir"])

    if args.sample:
        samples = [s for s in samples if s["id"] == args.sample]
        if not samples:
            sys.exit(f"ERROR: --sample {args.sample!r} not found in config")

    log(f"Tissue: {tissue} | Samples: {len(samples)}")

    # Output directories
    cb_dir = results_dir / "h5ad" / "01_cellbender"
    plots_dir = results_dir / "plots" / "01_ambient"
    extract_dir = results_dir / "h5ad" / "01_cellbender_raw_extracted"
    cb_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)
    extract_dir.mkdir(parents=True, exist_ok=True)

    cuda = not args.no_cuda

    # Check CUDA available
    if cuda:
        try:
            import torch
            if not torch.cuda.is_available():
                warn("CUDA not available — falling back to CPU. "
                     "Pass --no-cuda to suppress this warning.")
                cuda = False
            else:
                log(f"GPU: {torch.cuda.get_device_name(0)}")
        except ImportError:
            warn("torch not importable — running CPU mode.")
            cuda = False

    summary_rows = []
    failed = []

    if args.parallel > 1:
        # Parallel mode using multiprocessing
        import multiprocessing as mp
        import functools
        log(f"Parallel mode: {args.parallel} concurrent samples")

        _run_one = functools.partial(
            _process_sample,
            cb_dir=cb_dir, extract_dir=extract_dir,
            epochs=args.epochs, fpr=args.fpr, cuda=cuda
        )

        with mp.Pool(processes=args.parallel) as pool:
            results = pool.map(_run_one, samples)
        for row in results:
            summary_rows.append(row)
            if row["status"] != "done":
                failed.append(row["sample_id"])
    else:
        for s in samples:
            row = _process_sample(
                s, cb_dir, extract_dir, args.epochs,
                args.fpr, cuda
            )
            summary_rows.append(row)
            if row["status"] != "done":
                failed.append(s["id"])

    # Write summary
    collect_summary(summary_rows, plots_dir / "summary_ambient_fraction.csv")

    # Final status
    n_done = sum(1 for r in summary_rows if r["status"] == "done")
    n_skip = sum(1 for r in summary_rows if r["status"] == "skipped")
    n_fail = len(failed)

    log(f"\n=== Phase 1 complete ===")
    log(f"  Done:    {n_done}")
    log(f"  Skipped: {n_skip} (already existed)")
    log(f"  Failed:  {n_fail}")
    if failed:
        log(f"  Failed samples: {failed}")
        log(f"  Re-run with --sample <id> to retry individual samples.")
        sys.exit(1)

    log(f"\nOutputs in: {cb_dir}")
    log(f"Summary:    {plots_dir}/summary_ambient_fraction.csv")
    log(f"\nNext: run Phase 2 QC:")
    log(f"  uv run python scripts/02_qc.py --config {args.config}")


def _process_sample(
    s: dict,
    cb_dir: Path,
    extract_dir: Path,
    epochs: int,
    fpr: float,
    cuda: bool,
) -> dict:
    """Process one sample — resolve raw h5, run CellBender, return status row."""
    sid = s["id"]
    out_dir = cb_dir / sid

    # Idempotency check
    filtered_h5 = out_dir / f"{sid}_filtered.h5"
    if filtered_h5.exists():
        log(f"[{sid}] already done — skipping.")
        return {
            "sample_id": sid,
            "expected_cells": "?",
            "status": "skipped",
            "filtered_h5": str(filtered_h5),
        }

    # Resolve raw h5
    raw_h5_field = s.get("raw_h5")
    if not raw_h5_field:
        warn(f"[{sid}] no raw_h5 in config — skipping.")
        return {
            "sample_id": sid,
            "expected_cells": "?",
            "status": "no_raw_h5",
            "filtered_h5": "",
        }

    try:
        raw_h5 = resolve_raw_h5(raw_h5_field, sid, extract_dir)
    except FileNotFoundError as e:
        warn(f"[{sid}] {e}")
        return {
            "sample_id": sid,
            "expected_cells": "?",
            "status": "raw_h5_not_found",
            "filtered_h5": "",
        }

    # Expected cells from filtered h5
    expected_cells = count_expected_cells(s["h5"])

    # Run CellBender
    result = run_cellbender(
        sample_id=sid,
        raw_h5=raw_h5,
        expected_cells=expected_cells,
        out_dir=out_dir,
        epochs=epochs,
        cuda=cuda,
        fpr=fpr,
    )

    if result is None:
        return {
            "sample_id": sid,
            "expected_cells": expected_cells,
            "status": "failed",
            "filtered_h5": "",
        }

    return {
        "sample_id": sid,
        "expected_cells": expected_cells,
        "status": "done",
        "filtered_h5": str(result),
    }


if __name__ == "__main__":
    main()
