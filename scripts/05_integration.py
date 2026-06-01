#!/usr/bin/env python
"""
05_integration.py — Phase 5: scVI integration.

Trains scVI on the concatenated post-doublet object from Phase 4, producing a
batch-corrected 30-dim latent embedding in adata.obsm["X_scVI"]. Computes neighbor
graphs + UMAPs both pre-integration (PCA on lognorm) and post-integration (on
the scVI latent) so you can see whether pool effects were corrected.

Key choices (per project doc §3, §5 Phase 5, §11):
  - batch_key=pool                 — the multiplexing/library structure
  - categorical_covariates=        — age, group, sex (preserved as biology)
  - continuous_covariates=         — pct_counts_mt
  - HVG subset for training        — uses adata.var["use_for_scvi"] from Phase 4
  - BF16 mixed precision           — Ada GPU; auto-falls-back to CPU
  - n_layers=2, n_latent=30, max_epochs=400, batch_size=1024, early stopping

Usage:
  uv run python scripts/05_integration.py --config config/dev.yaml
  uv run python scripts/05_integration.py --config config/brain.yaml

Inputs:
  Concatenated h5ad from Phase 4 at
    {results_dir}/h5ad/05_integration_ready/all_samples.h5ad

Outputs:
  {results_dir}/h5ad/06_integrated/all_samples.h5ad   — integrated AnnData
  {results_dir}/h5ad/06_integrated/scvi_model/        — trained scVI model dir
  {results_dir}/plots/05_integration/*.png            — pre/post UMAPs, loss
  {results_dir}/tables/scvi_training_history.csv      — per-epoch loss
"""

import argparse
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

def compute_pre_integration_umap(adata):
    """PCA on lognorm + neighbors + UMAP. Used as 'before integration' baseline.
    Writes to adata.obsm['X_pca_pre'] and adata.obsm['X_umap_pre']."""
    print("  Computing pre-integration UMAP (PCA on lognorm)...")
    # Subset to HVGs that survive exclusion (same gene set scVI will use), then
    # PCA on the log-normalized layer. Operate on a copy so we don't mutate
    # the input AnnData's main slots.
    hvg = adata.var["use_for_scvi"].values
    tmp = adata[:, hvg].copy()
    tmp.X = tmp.layers["lognorm"].copy()
    sc.pp.scale(tmp, max_value=10)
    sc.tl.pca(tmp, n_comps=min(50, tmp.n_vars - 1))
    sc.pp.neighbors(tmp, use_rep="X_pca")
    sc.tl.umap(tmp)
    adata.obsm["X_pca_pre"] = tmp.obsm["X_pca"]
    adata.obsm["X_umap_pre"] = tmp.obsm["X_umap"]


def compute_post_integration_umap(adata):
    """Neighbors + UMAP on scVI latent. Writes to adata.obsm['X_umap']."""
    print("  Computing post-integration UMAP (on scVI latent)...")
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)


def plot_umap_panel(adata, basis: str, color_keys: list[str], out: Path, title: str):
    """4-panel UMAP colored by pool, age, group, sex. basis is 'X_umap_pre' or 'X_umap'."""
    n = len(color_keys)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.5))
    if n == 1:
        axes = [axes]
    # scanpy's pl.umap reads obsm["X_umap"] by default — swap if needed
    if basis != "X_umap":
        adata.obsm["X_umap"], saved = adata.obsm[basis].copy(), adata.obsm.get("X_umap")
    try:
        for ax, key in zip(axes, color_keys):
            sc.pl.umap(adata, color=key, ax=ax, show=False, frameon=False,
                       legend_fontsize=7, size=8)
            ax.set_title(f"{title}: {key}")
    finally:
        if basis != "X_umap" and saved is not None:
            adata.obsm["X_umap"] = saved
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_loss_curve(history: pd.DataFrame, out: Path):
    """Loss curves from scVI's training history."""
    fig, ax = plt.subplots(figsize=(7, 4))
    # scVI history columns vary slightly by version; pick what's present
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
    parser.add_argument("--cpu", action="store_true", help="Force CPU (override GPU detection)")
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

    # scVI training config — pull from YAML if present, else use doc §6 defaults
    scvi_cfg = cfg.get("scvi", {})
    n_latent     = int(scvi_cfg.get("n_latent", 30))
    n_layers     = int(scvi_cfg.get("n_layers", 2))
    max_epochs   = int(scvi_cfg.get("max_epochs", 400))
    batch_size   = int(scvi_cfg.get("batch_size", 1024))
    es_patience  = int(scvi_cfg.get("early_stopping_patience", 30))
    seed         = int(cfg.get("random_seed", 42))

    accelerator, precision = select_accelerator(force_cpu=args.cpu)
    # For dev runs on tiny data, fewer epochs is enough and avoids overfitting
    if adata.n_obs < 5000:
        max_epochs = min(max_epochs, 50)
        print(f"  Small dataset (n={adata.n_obs}) — capping max_epochs at {max_epochs}")

    print(f"\n[2/5] scVI setup (accelerator={accelerator}, precision={precision}, seed={seed})")
    scvi.settings.seed = seed

    # Subset to the HVG set scVI should learn on. Raw counts are in .X already.
    adata_hvg = adata[:, adata.var["use_for_scvi"]].copy()

    scvi.model.SCVI.setup_anndata(
        adata_hvg,
        batch_key="pool",
        categorical_covariate_keys=["age", "group", "sex"],
        continuous_covariate_keys=["pct_counts_mt"],
    )
    model = scvi.model.SCVI(
        adata_hvg, n_layers=n_layers, n_latent=n_latent,
    )

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

    # Save trained model + training history
    if model_dir.exists():
        # SCVI.save with overwrite=True requires the model_dir to exist already
        # or it'll handle creation; we just need it to not error on stale dir
        import shutil
        shutil.rmtree(model_dir)
    model.save(str(model_dir), overwrite=True)
    history = pd.concat({k: pd.DataFrame(v) for k, v in model.history.items()}, axis=1)
    history.columns = history.columns.droplevel(1)
    history.to_csv(paths["tables"] / "scvi_training_history.csv", index=True)
    plot_loss_curve(history, plot_dir / "scvi_loss_curve.png")
    print(f"  Trained: actual epochs = {len(history)}")

    # Write latent back onto the FULL adata (not just HVG subset) so downstream
    # phases can use the full gene set for marker analysis etc.
    print(f"\n[4/5] Extracting latent + computing UMAPs...")
    adata.obsm["X_scVI"] = model.get_latent_representation()

    compute_pre_integration_umap(adata)
    compute_post_integration_umap(adata)

    color_keys = ["pool", "age", "group", "sex"]
    plot_umap_panel(adata, "X_umap_pre", color_keys,
                    plot_dir / "umap_pre_integration.png", title="pre")
    plot_umap_panel(adata, "X_umap", color_keys,
                    plot_dir / "umap_post_integration.png", title="post")

    print(f"\n[5/5] Writing integrated h5ad...")
    # Drop lognorm layer (cheap to recompute from raw counts when needed —
    # see project doc §3 on not carrying redundant layers at scale).
    # `counts` layer is also redundant with .X; drop too.
    for layer in ("lognorm", "counts"):
        if layer in adata.layers:
            del adata.layers[layer]
    out_path = out_dir / "all_samples.h5ad"
    adata.write_h5ad(out_path)
    print(f"  Wrote {out_path}")
    print(f"  Model: {model_dir}")
    print(f"  Plots: {plot_dir}")
    print(f"\n✓ Phase 5 complete.\n")
    print(f"Inspect umap_pre_integration.png vs umap_post_integration.png:")
    print(f"  - pool: should go from clustered (pre) to mixed (post)")
    print(f"  - age/group/sex: biology should stay clustered both pre and post")
    print(f"\nNote: lognorm layer dropped from saved h5ad. To re-add in a notebook:")
    print(f"  adata.layers['lognorm'] = adata.X.copy()")
    print(f"  sc.pp.normalize_total(adata, target_sum=1e4, layer='lognorm')")
    print(f"  sc.pp.log1p(adata, layer='lognorm')\n")


if __name__ == "__main__":
    main()
