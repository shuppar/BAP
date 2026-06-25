#!/usr/bin/env python
"""h09c_integrate.py -- scVI integration + Leiden clustering for the human placenta cohort.

Mirrors mouse Phase 5 (scVI) + Phase 6 (Leiden), compressed into one script since the
human path has no separate Phase 4/6. Adaptations:
  batch_key = sample_id (no pools); HVG computed in-script; placenta HVG exclusions
  translated to human symbols. continuous_covariate_keys=["pct_counts_mt"], no categorical
  covariates, n_latent=30/n_layers=2, BF16 on GPU via _utils.select_accelerator -- all as mouse.

HVG exclusions (integration only -- still scored at annotation): mito MT-, ribo RPS/RPL,
hemoglobin, sex-linked (XIST + Y genes), placenta hormone/secretory (PSG*/CGB*/CGA/CSH/GH2).
These mirror mouse mito/ribo/hemo/sex + placenta Prl*/Psg*/Cgb*/Cga.

Usage (from project root; pre-flight: nvidia-smi --query-gpu=memory.used):
  uv run python scripts/h09c_integrate.py
  uv run python scripts/h09c_integrate.py --cpu        # force CPU
  uv run python scripts/h09c_integrate.py --reuse-model
"""
import argparse
import re
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc
import scvi

sys.path.insert(0, str(Path(__file__).parent))
from _utils import select_accelerator  # noqa: E402

GSE_DIR = Path("data/human_validation/placenta/gunter_rahman_2025_GSE271976")

HEMO = {"HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBM", "HBQ1", "HBZ"}
Y_GENES = {"RPS4Y1", "RPS4Y2", "DDX3Y", "UTY", "EIF1AY", "KDM5D", "USP9Y",
           "NLGN4Y", "ZFY", "TXLNGY", "TMSB4Y", "PRKY"}
SEX = {"XIST"} | Y_GENES
# placenta hormone/secretory families (human analogs of mouse Prl*/Psg*/Cgb*/Cga)
HORMONE_EXPLICIT = {"CGA", "CSH1", "CSH2", "CSHL1", "GH2"}
HORMONE_PREFIX = ("PSG", "CGB", "PRL")
EXCLUDE_PREFIX = ("MT-", "RPS", "RPL") + HORMONE_PREFIX
SANITY_MARKERS = ["KRT7", "CGA", "HLA-G", "PAEP", "VIM", "PTPRC", "PECAM1", "HBB"]


def build_exclusion(var_names):
    excl = set(HEMO) | SEX | HORMONE_EXPLICIT
    pref = tuple(EXCLUDE_PREFIX)
    return {g for g in var_names if g in excl or g.startswith(pref)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpu", action="store_true")
    ap.add_argument("--reuse-model", action="store_true", dest="reuse_model")
    ap.add_argument("--n-top-hvg", type=int, default=2000)
    ap.add_argument("--n-latent", type=int, default=30)
    ap.add_argument("--resolution", type=float, default=1.0)
    ap.add_argument("--n-neighbors", type=int, default=30)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    in_path = GSE_DIR / "h5ad" / "h09b_qc_doublets.h5ad"
    if not in_path.is_file():
        sys.exit(f"missing {in_path} (run h09b first)")
    out_dir = GSE_DIR / "h5ad"
    plot_dir = GSE_DIR / "plots" / "h09c"
    plot_dir.mkdir(parents=True, exist_ok=True)
    model_dir = out_dir / "scvi_model"

    print(f"[h09c] loading {in_path}")
    adata = sc.read_h5ad(in_path)
    sc.pp.filter_genes(adata, min_cells=10)   # drop near-zero genes (seurat_v3 loess needs variance)
    adata.layers["counts"] = adata.X.copy()
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes (after min_cells=10 gene filter)")
    if "pct_counts_mt" not in adata.obs.columns:
        sys.exit("pct_counts_mt missing -- rerun h09b")

    # --- HVG: exclude hormone/mito/ribo/hemo/sex genes FIRST, then HVG on the rest ---
    # Placental hormone genes (CGA/CGB/CSH/PSG) have extreme mean -> singular seurat_v3 loess.
    # They're excluded from scVI anyway, so dropping them up front both fixes loess and is cleaner.
    excluded = build_exclusion(adata.var_names)
    adata.var["excluded_from_scvi"] = adata.var_names.isin(excluded)
    cand = adata[:, ~adata.var["excluded_from_scvi"]].copy()
    try:
        sc.pp.highly_variable_genes(cand, n_top_genes=args.n_top_hvg, flavor="seurat_v3",
                                    batch_key="sample_id", layer="counts")
        method = "seurat_v3 (batch-aware)"
    except Exception:
        try:
            sc.pp.highly_variable_genes(cand, n_top_genes=args.n_top_hvg, flavor="seurat_v3",
                                        layer="counts")
            method = "seurat_v3 global (batch-aware loess failed)"
        except Exception:
            tmp = cand.copy()                       # bin-based seurat flavor: no loess, always works
            sc.pp.normalize_total(tmp, target_sum=1e4); sc.pp.log1p(tmp)
            sc.pp.highly_variable_genes(tmp, n_top_genes=args.n_top_hvg, flavor="seurat")
            cand.var["highly_variable"] = tmp.var["highly_variable"].values
            method = "seurat-flavor on lognorm (seurat_v3 loess failed)"
    hv = set(cand.var_names[cand.var["highly_variable"]])
    adata.var["use_for_scvi"] = adata.var_names.isin(hv)
    n_hvg = int(adata.var["use_for_scvi"].sum())
    n_excl = int(adata.var["excluded_from_scvi"].sum())
    print(f"  HVG method: {method}")
    print(f"  HVGs for scVI: {n_hvg} ({n_excl} hormone/mito/ribo/hemo/sex genes excluded up front)")
    if n_hvg < 500:
        sys.exit(f"only {n_hvg} HVGs -- check exclusion list / gene symbols")

    accelerator, precision = select_accelerator(force_cpu=args.cpu)
    scvi.settings.seed = args.seed
    adata_hvg = adata[:, adata.var["use_for_scvi"]].copy()
    adata_hvg.X = adata_hvg.layers["counts"].copy()  # scVI models raw (SoupX-corrected) counts
    scvi.model.SCVI.setup_anndata(adata_hvg, batch_key="sample_id",
                                  continuous_covariate_keys=["pct_counts_mt"])

    if args.reuse_model:
        if not model_dir.exists():
            sys.exit(f"--reuse-model set but no model at {model_dir}")
        print(f"[h09c] loading trained model {model_dir} (skip training)")
        model = scvi.model.SCVI.load(str(model_dir), adata=adata_hvg)
    else:
        model = scvi.model.SCVI(adata_hvg, n_layers=2, n_latent=args.n_latent)
        print(f"[h09c] training scVI (accelerator={accelerator}, precision={precision})")
        model.train(max_epochs=400, batch_size=1024, early_stopping=True,
                    early_stopping_patience=30, accelerator=accelerator,
                    devices=1, precision=precision)
        if model_dir.exists():
            shutil.rmtree(model_dir)
        model.save(str(model_dir), overwrite=True)
        print(f"  trained epochs = {len(model.history['elbo_train'])}")

    # --- latent -> neighbors -> Leiden (igraph) -> UMAP ---
    adata.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=args.n_neighbors, random_state=args.seed)
    sc.tl.leiden(adata, resolution=args.resolution, flavor="igraph",
                 n_iterations=2, directed=False, random_state=args.seed)
    sc.tl.umap(adata, min_dist=0.3, spread=1.2, random_state=args.seed)
    print(f"[h09c] {adata.obs['leiden'].nunique()} Leiden clusters at res={args.resolution}")

    # --- diagnostic UMAPs ---
    adata.obs["leiden"] = adata.obs["leiden"].astype("category")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)  # lognorm in .X for the marker UMAP; counts preserved in layer
    for keys, fn, t in [
        (["sample_id", "condition", "side", "leiden"], "umap_batch.png", "integration"),
        ([m for m in SANITY_MARKERS if m in adata.var_names], "umap_markers.png", "markers"),
    ]:
        n = len(keys)
        fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 4))
        axes = [axes] if n == 1 else axes
        for ax, k in zip(axes, keys):
            sc.pl.umap(adata, color=k, ax=ax, show=False, frameon=False, size=5,
                       legend_fontsize=6, legend_loc="on data" if k == "leiden" else "right margin")
            for c in ax.collections:
                c.set_rasterized(True)
        fig.tight_layout()
        fig.savefig(plot_dir / fn, dpi=140, bbox_inches="tight")
        plt.close(fig)

    out_path = out_dir / "h09c_integrated.h5ad"
    adata.write_h5ad(out_path)
    print(f"[h09c] wrote {out_path}\n  plots -> {plot_dir}")
    print("  inspect umap_batch.png (samples should mix; condition/side stay real) "
          "+ umap_markers.png (compartments separate) before annotating in h09d")


if __name__ == "__main__":
    main()
