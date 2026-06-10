"""
05b_umap_sweep.py — Phase 5b: UMAP seed sweep on the integrated h5ad.

scVI training is locked. Only the UMAP projection is non-deterministic.
This script sweeps multiple seeds on the SAME scVI latent (no retraining,
no neighbor recomputation past the first seed), writes per-covariate
comparison plots, then optionally overwrites the chosen seed back to the
integrated h5ad.

Two modes:

  Sweep (default):
      uv run python scripts/05b_umap_sweep.py --config config/brain.yaml \\
          --seeds 42,0,7,123,2024
    - reads results/{tissue}/h5ad/06_integrated/all_samples.h5ad
    - computes neighbors ONCE on X_scVI
    - computes UMAP for each seed (~5-10 min per seed on 500K-cell scale)
    - saves N-panel comparison plots in results/{tissue}/plots/05b_umap_sweep/
    - DOES NOT modify the integrated h5ad

  Apply:
      uv run python scripts/05b_umap_sweep.py --config config/brain.yaml \\
          --apply 7
    - recomputes the chosen seed's UMAP
    - overwrites .obsm['X_umap'] in the integrated h5ad in place
    - records .uns['umap_seed_applied'] = 7
    - Phase 6 (clustering) inherits this UMAP automatically

UMAP hyperparameters locked at scanpy defaults: n_neighbors=15, min_dist=0.5,
spread=1.0, init_pos='spectral', metric='euclidean'.
"""

import argparse
import sys
import time
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import yaml


COVARIATES_CATEGORICAL = ["pool", "age", "group", "sex", "phase"]
COVARIATES_CONTINUOUS = ["pct_counts_mt"]


def log(msg: str) -> None:
    print(f"[umap_sweep] {msg}", flush=True)


def load_config(path: Path) -> dict:
    with path.open() as f:
        return yaml.safe_load(f)


def resolve_integrated_h5ad(cfg: dict) -> Path:
    results_dir = Path(cfg["results_dir"])
    return results_dir / "h5ad" / "06_integrated" / "all_samples.h5ad"


def compute_umaps(
    adata: ad.AnnData,
    seeds: list[int],
    n_neighbors: int = 15,
) -> None:
    """Compute UMAP for each seed; store each in obsm[f'X_umap_seed{seed}'].

    Mutates adata in place. Neighbors computed ONCE (k-NN is deterministic
    given the data + parameters; random_state only affects the connectivity
    smoothing, which we fix at the first seed for reproducibility).
    """
    log(f"Computing neighbors on X_scVI (n_neighbors={n_neighbors})...")
    t0 = time.time()
    sc.pp.neighbors(
        adata,
        use_rep="X_scVI",
        n_neighbors=n_neighbors,
        random_state=seeds[0],
    )
    log(f"  neighbors done in {time.time()-t0:.1f}s")

    for seed in seeds:
        log(f"Computing UMAP for seed={seed}...")
        t0 = time.time()
        sc.tl.umap(adata, random_state=seed)
        adata.obsm[f"X_umap_seed{seed}"] = adata.obsm["X_umap"].copy()
        log(f"  seed={seed} done in {time.time()-t0:.1f}s")


def plot_seed_comparison(
    adata: ad.AnnData,
    seeds: list[int],
    out_dir: Path,
) -> None:
    """One figure per covariate; N subplots side-by-side, one per seed."""
    out_dir.mkdir(parents=True, exist_ok=True)
    n_seeds = len(seeds)

    available_cat = [c for c in COVARIATES_CATEGORICAL if c in adata.obs.columns]
    available_cont = [c for c in COVARIATES_CONTINUOUS if c in adata.obs.columns]
    log(f"Plotting covariates: categorical={available_cat}, continuous={available_cont}")

    for covar in available_cat + available_cont:
        # Scale figure with number of seeds
        fig, axes = plt.subplots(1, n_seeds, figsize=(5 * n_seeds, 5),
                                  constrained_layout=True)
        if n_seeds == 1:
            axes = [axes]
        for ax, seed in zip(axes, seeds):
            # Point scanpy at the per-seed obsm key by swapping X_umap
            adata.obsm["X_umap"] = adata.obsm[f"X_umap_seed{seed}"]
            sc.pl.umap(
                adata,
                color=covar,
                ax=ax,
                show=False,
                frameon=False,
                title=f"seed={seed}",
                size=2,
                legend_loc="right margin" if covar in available_cat else None,
            )
        fig.suptitle(f"UMAP seed sweep — colored by {covar}", fontsize=14)
        out_path = out_dir / f"{covar}_seed_comparison.png"
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        log(f"  wrote {out_path}")


def apply_seed(
    adata: ad.AnnData,
    seed: int,
    h5ad_path: Path,
    n_neighbors: int = 15,
) -> None:
    """Recompute UMAP for the chosen seed and overwrite the integrated h5ad."""
    log(f"Apply mode — overwriting .obsm['X_umap'] with seed={seed}")
    log("Computing neighbors...")
    sc.pp.neighbors(
        adata,
        use_rep="X_scVI",
        n_neighbors=n_neighbors,
        random_state=seed,
    )
    log(f"Computing UMAP (random_state={seed})...")
    sc.tl.umap(adata, random_state=seed)
    adata.uns["umap_seed_applied"] = seed
    # Clean any per-seed obsm keys from prior sweep runs so we don't bloat
    for key in list(adata.obsm.keys()):
        if key.startswith("X_umap_seed"):
            del adata.obsm[key]
    log(f"Writing {h5ad_path}...")
    adata.write_h5ad(h5ad_path)
    log("Done. Phase 6 will inherit this UMAP.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Phase 5b: UMAP seed sweep")
    parser.add_argument("--config", required=True, type=Path,
                        help="Tissue YAML config (brain.yaml / placenta.yaml)")
    parser.add_argument("--seeds", type=str, default="42,0,7,123,2024",
                        help="Comma-separated seed list (sweep mode)")
    parser.add_argument("--apply", type=int, default=None,
                        help="Apply mode: recompute this single seed and "
                             "overwrite the integrated h5ad in place")
    parser.add_argument("--n-neighbors", type=int, default=15, dest="n_neighbors",
                        help="n_neighbors for sc.pp.neighbors (default 15)")
    args = parser.parse_args()

    cfg = load_config(args.config)
    tissue = cfg.get("tissue", "unknown")
    h5ad_path = resolve_integrated_h5ad(cfg)
    out_dir = Path(cfg["results_dir"]) / "plots" / "05b_umap_sweep"

    log(f"=== Phase 5b: UMAP seed sweep ({tissue}) ===")
    log(f"Integrated h5ad: {h5ad_path}")
    if not h5ad_path.exists():
        sys.exit(f"ERROR: integrated h5ad not found. Run Phase 5 first.")

    log(f"Loading h5ad...")
    t0 = time.time()
    adata = ad.read_h5ad(h5ad_path)
    log(f"  loaded {adata.n_obs:,} cells × {adata.n_vars:,} genes in {time.time()-t0:.1f}s")

    if "X_scVI" not in adata.obsm:
        sys.exit("ERROR: X_scVI not in .obsm. Phase 5 didn't write the latent.")

    if args.apply is not None:
        apply_seed(adata, args.apply, h5ad_path, n_neighbors=args.n_neighbors)
        return

    # Sweep mode
    seeds = [int(s.strip()) for s in args.seeds.split(",")]
    log(f"Sweep seeds: {seeds}")

    compute_umaps(adata, seeds, n_neighbors=args.n_neighbors)
    plot_seed_comparison(adata, seeds, out_dir)

    log("")
    log(f"✓ Sweep complete. Review plots in {out_dir}")
    log(f"  Pick the cleanest seed and re-run with --apply <seed>:")
    log(f"  uv run python scripts/05b_umap_sweep.py --config {args.config} --apply <seed>")


if __name__ == "__main__":
    main()
