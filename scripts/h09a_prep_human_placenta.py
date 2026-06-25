#!/usr/bin/env python
"""h09a_prep_human_placenta.py -- SoupX prep for Gunter-Rahman GSE271976 (human term placenta).

GEO shipped raw_feature_bc_matrix only, so per sample we knee/inflection cell-call + SoupX
(via the h_run_soupx_from_raw.R subprocess), then assemble the corrected MTX trio into a
per-sample h5ad with parsed obs metadata. Mirrors mouse 02_soupx.py; parallel via parallel_map.

Filenames: GSM<digits>_<fs|ms>_<lean|mo>_<rep>_raw_feature_bc_matrix.h5
  fs/ms = fetal-/maternal-facing side; lean/mo = lean/maternal-obese; rep = replicate index.
Each GSM is treated as an independent sample (donor_id == sample_id); fs/ms pairing across
the same placenta is NOT assumed (mo has 6 fs vs 7 ms libraries -> not 1:1).

Usage (from project root):
  uv run python scripts/h09a_prep_human_placenta.py --n-jobs 12
  # smoke one sample:        --sample-ids fs_lean_1
  # re-assemble, skip SoupX: --skip-soupx
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from _utils import parallel_map  # noqa: E402

GSE_DIR = Path("data/human_validation/placenta/gunter_rahman_2025_GSE271976")
RWORKER = "scripts/h_run_soupx_from_raw.R"

FNAME_RE = re.compile(r"(GSM\d+)_(fs|ms)_(lean|mo)_(\d+)_raw_feature_bc_matrix\.h5$")
SIDE = {"fs": "fetal", "ms": "maternal"}
COND = {"lean": "lean", "mo": "obese"}


def parse_meta(h5: Path):
    m = FNAME_RE.search(h5.name)
    if not m:
        return None
    gsm, side, cond, rep = m.groups()
    return dict(
        sample_id=f"{side}_{cond}_{rep}",
        gsm=gsm,
        side=SIDE[side],
        condition=COND[cond],
        replicate=int(rep),
        h5=str(h5),
    )


def run_soupx(job):
    """Call the R worker for one sample. Raises on failure (captured by parallel_map)."""
    out_dir = GSE_DIR / "soupx" / job["sample_id"]
    cmd = [
        "Rscript", RWORKER,
        "--h5", job["h5"],
        "--out-dir", str(out_dir),
        "--sample-id", job["sample_id"],
        "--cutoff", job["cutoff"],
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{job['sample_id']} R worker failed:\n{r.stderr[-2000:]}")
    return r.stdout.strip()


def assemble(meta):
    """Read corrected MTX trio -> per-sample h5ad with obs metadata. Returns summary dict."""
    sid = meta["sample_id"]
    sample_dir = GSE_DIR / "soupx" / sid
    a = sc.read_10x_mtx(sample_dir, var_names="gene_symbols")
    a.var_names_make_unique()
    for k in ("sample_id", "gsm", "side", "condition", "replicate"):
        a.obs[k] = meta[k]
    a.obs["donor_id"] = sid
    a.obs_names = [f"{sid}_{bc}" for bc in a.obs_names]

    h5ad_dir = GSE_DIR / "h5ad"
    h5ad_dir.mkdir(parents=True, exist_ok=True)
    a.write(h5ad_dir / f"{sid}.h5ad")

    summ = json.loads((sample_dir / "summary.json").read_text())
    return dict(
        sample_id=sid, gsm=meta["gsm"], side=meta["side"], condition=meta["condition"],
        replicate=meta["replicate"], n_cells=a.n_obs,
        rho_mean=summ["rho_mean"], pct_removed=summ["pct_removed"],
        knee=summ["knee"], inflection=summ["inflection"], cutoff=summ["cutoff"],
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cutoff", default="inflection", choices=["knee", "inflection"])
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--sample-ids", nargs="*", help="subset of sample_ids (smoke test)")
    ap.add_argument("--skip-soupx", action="store_true", help="re-assemble h5ads without rerunning SoupX")
    args = ap.parse_args()

    metas = [m for h5 in sorted(GSE_DIR.glob("*_raw_feature_bc_matrix.h5"))
             if (m := parse_meta(h5)) is not None]
    if args.sample_ids:
        metas = [m for m in metas if m["sample_id"] in set(args.sample_ids)]
    if not metas:
        sys.exit("no matching h5 files found")
    for m in metas:
        m["cutoff"] = args.cutoff
    print(f"[h09a] {len(metas)} samples, cutoff={args.cutoff}, n_jobs={args.n_jobs}")

    if not args.skip_soupx:
        fails = []
        for job, out, err in parallel_map(run_soupx, metas, n_jobs=args.n_jobs, desc="SoupX"):
            if err:
                fails.append((job["sample_id"], err))
                print(f"  [FAIL] {job['sample_id']}: {err.splitlines()[-1]}")
            else:
                print(f"  {out}")
        if fails:
            sys.exit(f"{len(fails)} sample(s) failed SoupX; fix before assembling")

    rows = [assemble(m) for m in metas]  # mtx read is fast; serial is fine
    manifest = pd.DataFrame(rows).sort_values(["condition", "side", "replicate"])
    mpath = GSE_DIR / "h09a_soupx_manifest.csv"
    manifest.to_csv(mpath, index=False)
    print(f"\n[h09a] wrote {len(rows)} h5ads + manifest -> {mpath}")
    print(manifest.to_string(index=False))
    print(f"\ntotal cells: {manifest['n_cells'].sum()}  (paper reported ~62,864 across 20)")


if __name__ == "__main__":
    main()
