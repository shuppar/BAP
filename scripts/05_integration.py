#!/usr/bin/env python
"""
05_integration.py — Phase 5: scVI integration.

Key changes from original:
  - Removed categorical_covariate_keys=[age, group, sex]. Those are nuisance
    removers in scVI — including them strips biological signal from the latent.
    batch_key=pool handles the technical correction; age/group/sex stay as
    real signal in the latent space.
  - Cell cycle: cc_difference is available in adata.obs from Phase 4. It is
    NOT conditioned on by default (proliferation is real biology at P1 and
    could be a stress effect). To enable: set scvi.condition_cell_cycle: true
    in the YAML — the script will add cc_difference to continuous_covariate_keys.
  - Post-integration UMAPs now include phase (cell cycle) as a coloring key,
    so you can immediately see whether cycle is driving any cluster structure
    before deciding whether to condition on it.

scVI setup:
  - batch_key = pool                      (technical batch — corrected)
  - continuous_covariate_keys = [pct_counts_mt]  (+ cc_difference if opted in)
  - No categorical covariates             (biology stays in latent)

Usage:
  uv run python scripts/05_integration.py --config config/dev.yaml
  uv run python scripts/05_integration.py --config config/brain.yaml

Inputs:
  {results_dir}/h5ad/05_integration_ready/all_samples.h5ad  (from Phase 4)

Outputs:
  {results_dir}/h5ad/06_integrated/all_samples.h5ad
  {results_dir}/h5ad/06_integrated/scvi_model/
  {results_dir}/plots/05_integration/
    - umap_pre_integration.png   : PCA-based UMAP before correction
    - umap_post_integration.png  : scVI latent UMAP after correction
    - umap_post_phase.png        : post-integration UMAP colored by cell cycle phase
    - scvi_loss_curve.png
  {results_dir}/tables/scvi_training_history.csv
"""

import argparse
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi

from _utils import load_config, phase_paths, select_accelerator


# ----------------------------------------------------------------------------
# UMAP helpers
# ----------------------------------------------------------------------------

def compute_pre_integration_umap(adata) -> None:
    """PCA on lognorm → neighbors → UMAP. Baseline before integration."""
    print("  Computing pre-integration UMAP (PCA on lognorm)...")
    hvg = adata.var["use_for_scvi"].values
    tmp = adata[:, hvg].copy()
    tmp.X = tmp.layers["lognorm"].copy()
    sc.pp.scale(tmp, max_value=10)
    sc.tl.pca(tmp, n_comps=min(50, tmp.n_vars - 1))
    sc.pp.neighbors(tmp, use_rep="X_pca")
    sc.tl.umap(tmp)
    adata.obsm["X_pca_pre"] = tmp.obsm["X_pca"]
    adata.obsm["X_umap_pre"] = tmp.obsm["X_umap"]


def compute_post_integration_umap(adata) -> None:
    """Neighbors → UMAP on scVI latent."""
    print("  Computing post-integration UMAP (on scVI latent)...")
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)


def plot_umap_panel(adata, basis: str, color_keys: list[str],
                    out: Path, title: str) -> None:
    """Multi-panel UMAP. basis is 'X_umap_pre' or 'X_umap'."""
    # Filter to keys that actually exist in obs
    color_keys = [k for k in color_keys if k in adata.obs.columns or k in adata.var_names]
    if not color_keys:
        return
    n = len(color_keys)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.5))
    if n == 1:
        axes = [axes]
    saved_umap = adata.obsm.get("X_umap")
    if basis != "X_umap":
        adata.obsm["X_umap"] = adata.obsm[basis].copy()
    try:
        for ax, key in zip(axes, color_keys):
            sc.pl.umap(adata, color=key, ax=ax, show=False, frameon=False,
                       legend_fontsize=7, size=8)
            ax.set_title(f"{title}: {key}")
    finally:
        if basis != "X_umap" and saved_umap is not None:
            adata.obsm["X_umap"] = saved_umap
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_loss_curve(history: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    candidates = [c for c in ("train_loss_epoch", "validation_loss",
                               "elbo_train", "elbo_validation",
                               "reconstruction_loss_train", "reconstruction_loss_validation")
                  if c in history.columns]
    for c in candidates:
        ax.plot(history.index, history[c], label=c, lw=1.2)
    ax.set_xlabel("epoch")
    ax.set_ylabel("loss")
    ax.set_title("scVI training history")
    if candidates:
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Phase 5: scVI integration")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--cpu", action="store_true", help="Force CPU")
    args = parser.parse_args()

    print(f"\n=== Phase 5: scVI integration ===")
    print(f"Config: {args.config}")

    cfg = load_config(args.config)
    prev_paths = phase_paths(cfg, "integration_prep")
    in_path = prev_paths["h5ad"] / "all_samples.h5ad"
    paths = phase_paths(cfg, "integration")
    out_dir = paths["h5ad"]
    plot_dir = paths["plots"]
    model_dir = out_dir / "scvi_model"

    if not in_path.is_file():
        sys.exit(f"ERROR: missing {in_path}. Run 04_integration_prep.py first.")

    print(f"\n[1/5] Loading {in_path}...")
    adata = sc.read_h5ad(in_path)
    print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    if "use_for_scvi" not in adata.var.columns:
        sys.exit("ERROR: adata.var['use_for_scvi'] not present. Re-run Phase 4.")
    n_hvg = int(adata.var["use_for_scvi"].sum())
    print(f"  HVGs for scVI: {n_hvg}")

    # Check Phase 4 cell cycle scores made it through
    has_cc = "S_score" in adata.obs.columns
    if not has_cc:
        print("  [warn] Cell cycle scores not found in adata.obs.")
        print("         Re-run Phase 4 to get S_score/G2M_score/phase/cc_difference.")

    # --- scVI config ---
    scvi_cfg = cfg.get("scvi", {})
    n_latent    = int(scvi_cfg.get("n_latent", 30))
    n_layers    = int(scvi_cfg.get("n_layers", 2))
    max_epochs  = int(scvi_cfg.get("max_epochs", 400))
    batch_size  = int(scvi_cfg.get("batch_size", 1024))
    es_patience = int(scvi_cfg.get("early_stopping_patience", 30))
    seed        = int(cfg.get("random_seed", 42))

    # Cell cycle conditioning: off by default (proliferation is real biology
    # at P1 and could itself be a stress effect). Enable in YAML if cycling
    # cells create unwanted cluster structure in the phase UMAP below.
    condition_cc = bool(scvi_cfg.get("condition_cell_cycle", False))
    continuous_covariates = ["pct_counts_mt"]
    if condition_cc and has_cc:
        continuous_covariates.append("cc_difference")
        print(f"  Cell cycle conditioning: ON (cc_difference added to continuous covariates)")
    else:
        print(f"  Cell cycle conditioning: OFF (default — check phase UMAP after training)")

    accelerator, precision = select_accelerator(force_cpu=args.cpu)
    if adata.n_obs < 5000:
        max_epochs = min(max_epochs, 50)
        print(f"  Small dataset (n={adata.n_obs}) — capping max_epochs at {max_epochs}")

    print(f"\n[2/5] scVI setup")
    print(f"  batch_key          = pool  (technical batch correction)")
    print(f"  categorical_covariates = none  (biology stays in latent)")
    print(f"  continuous_covariates  = {continuous_covariates}")
    print(f"  accelerator={accelerator}, precision={precision}, seed={seed}")

    scvi.settings.seed = seed
    adata_hvg = adata[:, adata.var["use_for_scvi"]].copy()

    # Check required obs columns exist
    missing_obs = [c for c in continuous_covariates if c not in adata_hvg.obs.columns]
    if missing_obs:
        print(f"  [warn] Missing obs columns for continuous covariates: {missing_obs}")
        continuous_covariates = [c for c in continuous_covariates
                                  if c in adata_hvg.obs.columns]

    scvi.model.SCVI.setup_anndata(
        adata_hvg,
        batch_key="pool",
        continuous_covariate_keys=continuous_covariates,
    )
    model = scvi.model.SCVI(adata_hvg, n_layers=n_layers, n_latent=n_latent)

    print(f"\n[3/5] Training scVI (max_epochs={max_epochs}, batch_size={batch_size}, "
          f"early_stopping_patience={es_patience})")
    model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        early_stopping=True,
        early_stopping_patience=es_patience,
        accelerator=accelerator,
        devices=1,
        precision=precision,
    )

    if model_dir.exists():
        shutil.rmtree(model_dir)
    model.save(str(model_dir), overwrite=True)
    history = pd.concat({k: pd.DataFrame(v) for k, v in model.history.items()}, axis=1)
    history.columns = history.columns.droplevel(1)
    history.to_csv(paths["tables"] / "scvi_training_history.csv", index=True)
    plot_loss_curve(history, plot_dir / "scvi_loss_curve.png")
    print(f"  Trained: actual epochs = {len(history)}")

    print(f"\n[4/5] Extracting latent + computing UMAPs...")
    adata.obsm["X_scVI"] = model.get_latent_representation()

    compute_pre_integration_umap(adata)
    compute_post_integration_umap(adata)

    # Pre-integration: pool, age, group, sex
    plot_umap_panel(adata, "X_umap_pre", ["pool", "age", "group", "sex"],
                    plot_dir / "umap_pre_integration.png", title="pre")

    # Post-integration: same biological keys
    plot_umap_panel(adata, "X_umap", ["pool", "age", "group", "sex"],
                    plot_dir / "umap_post_integration.png", title="post")

    # Cell cycle phase UMAP — key diagnostic: does phase drive clusters?
    # If yes and it looks like a technical artifact, enable condition_cell_cycle in YAML.
    if has_cc:
        plot_umap_panel(adata, "X_umap", ["phase", "S_score", "G2M_score"],
                        plot_dir / "umap_post_phase.png", title="cell cycle")

    print(f"\n[5/5] Writing integrated h5ad...")
    for layer in ("lognorm", "counts"):
        if layer in adata.layers:
            del adata.layers[layer]
    out_path = out_dir / "all_samples.h5ad"
    adata.write_h5ad(out_path)

    print(f"  Wrote {out_path}")
    print(f"  Model: {model_dir}")
    print(f"  Plots: {plot_dir}")
    print(f"\n✓ Phase 5 complete.")
    print(f"\nKey diagnostics:")
    print(f"  umap_pre_integration.png  → pool should cluster (batch effect visible)")
    print(f"  umap_post_integration.png → pool should mix; age/group/sex stay separated")
    print(f"  umap_post_phase.png       → does cell cycle phase drive any clusters?")
    print(f"    If yes: set scvi.condition_cell_cycle: true in YAML and re-run Phase 5.")
    print(f"\nNext step: Phase 6 clustering\n")


if __name__ == "__main__":
    main()
