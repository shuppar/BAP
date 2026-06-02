#!/usr/bin/env python
"""
07b_subcluster.py — Phase 7b: subcluster a chosen cell type to resolve subtypes.

The joint scVI latent (Phase 5) is optimized to separate COARSE cell types, so
it lacks resolution for within-type structure (e.g. microglial states, oligo
lineage stages). This script subsets to one coarse cell type and RE-RUNS
HVG selection + scVI + Leiden on just those cells, so subtle subtypes that were
invisible in the joint embedding can separate.

Tissue-agnostic: the cell type to subcluster is a CLI arg; the label column and
batch key come from config / sensible defaults. Works for brain (microglia,
oligodendrocyte_lineage, ...) or placenta (trophoblast, ...) identically.

Which label column identifies the target cell type (priority):
  1. --label-key on the CLI (explicit override)
  2. manual_annotation   (if filled in and non-empty)
  3. scanvi_celltype     (from Phase 7c, if present)
  4. celltypist_majority (from Phase 7)
  5. provisional_celltype (marker-score fallback from Phase 7)
The chosen column and the exact --celltype value must match; the script lists
available labels and exits if the value isn't found (no silent empty subset).

Usage:
  # subcluster microglia from the dev brain object:
  uv run python scripts/07b_subcluster.py --config config/dev.yaml --celltype Microglia
  # explicit label column + custom resolution:
  uv run python scripts/07b_subcluster.py --config config/brain.yaml \\
      --celltype "Oligodendrocytes" --label-key manual_annotation --resolution 0.6

Inputs (first that exists, newest first):
  {results_dir}/h5ad/08b_label_transferred/all_samples.h5ad   (Phase 7c)
  {results_dir}/h5ad/08_annotated/all_samples.h5ad            (Phase 7)

Outputs (slug = sanitized celltype name):
  {results_dir}/h5ad/08c_subclustered/{slug}.h5ad
  {results_dir}/plots/07b_subcluster/{slug}/
    - umap_subclusters.png        : subclusters on the re-integrated subset
    - umap_by_group.png           : subset UMAP colored by stress group
    - umap_by_sample.png          : colored by sample (subcluster batch check)
    - subcluster_markers_dotplot.png
  {results_dir}/tables/
    - subcluster_{slug}_markers.csv
    - subcluster_{slug}_composition.csv   : subcluster × group/sample
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from _utils import load_config, add_lognorm, select_accelerator


LABEL_KEY_PRIORITY = [
    "manual_annotation", "scanvi_celltype", "celltypist_majority", "provisional_celltype",
]


def resolve_label_key(adata, explicit: str | None) -> str:
    """Pick the obs column that identifies cell types. Hard-fail if none usable."""
    if explicit:
        if explicit not in adata.obs.columns:
            sys.exit(f"ERROR: --label-key '{explicit}' not in adata.obs. "
                     f"Available: {list(adata.obs.columns)}")
        return explicit
    for key in LABEL_KEY_PRIORITY:
        if key in adata.obs.columns:
            # manual_annotation is only usable if actually filled in
            if key == "manual_annotation":
                vals = adata.obs[key].astype(str)
                if (vals.str.len() == 0).all() or vals.eq("").all():
                    continue
            print(f"  Using label column: '{key}'")
            return key
    sys.exit(
        "ERROR: no usable cell-type label column found.\n"
        f"  Looked for: {LABEL_KEY_PRIORITY}\n"
        "  Run Phase 7 (07_annotation.py) first, or pass --label-key explicitly."
    )


def slugify(name: str) -> str:
    return re.sub(r"[^0-9a-zA-Z]+", "_", name).strip("_").lower()


def reintegrate_subset(adata, seed, accelerator, precision, n_hvg, batch_key="pool"):
    """Re-run HVG + scVI on the subset to get a subtype-resolved latent.

    Mirrors Phase 4/5 logic at smaller scale: lognorm for HVG/markers, raw counts
    in .X for scVI, batch_key=pool, no categorical covariates.
    """
    import scvi

    # Need raw counts in .X. The annotated object stores raw counts in .X
    # (lognorm was dropped after Phase 5), so this holds — but verify loudly.
    Xmax = adata.X.max()
    if not np.isclose(Xmax, np.round(Xmax)) or Xmax < 0:
        sys.exit("ERROR: subset .X does not look like raw counts (scVI needs raw). "
                 f"max={Xmax}. Aborting rather than training on the wrong matrix.")

    # Only one batch level left after subsetting? scVI still runs, but warn.
    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        print(f"  [info] only {n_batches} '{batch_key}' level in this subset — "
              f"no batch correction needed within it.")

    # HVG on the subset (seurat_v3 on raw counts). Cap n_hvg at n_vars.
    n_hvg = min(n_hvg, adata.n_vars - 1)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat_v3",
                                batch_key=batch_key if n_batches >= 2 else None)
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    scvi.settings.seed = seed
    setup_kwargs = {} if n_batches < 2 else {"batch_key": batch_key}
    scvi.model.SCVI.setup_anndata(adata_hvg, **setup_kwargs)
    model = scvi.model.SCVI(adata_hvg, n_layers=2, n_latent=30)
    max_epochs = 50 if adata.n_obs < 5000 else 200
    print(f"  Training scVI on subset (n={adata.n_obs:,}, max_epochs={max_epochs})...")
    model.train(max_epochs=max_epochs, accelerator=accelerator, devices=1,
                precision=precision, early_stopping=True)

    adata.obsm["X_scVI_sub"] = model.get_latent_representation()
    return adata


def main():
    parser = argparse.ArgumentParser(description="Phase 7b: subcluster one cell type")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--celltype", required=True,
                        help="Cell type value to subcluster (must match the label column)")
    parser.add_argument("--label-key", default=None,
                        help="obs column holding cell type labels (default: auto-detect)")
    parser.add_argument("--resolution", type=float, default=0.6,
                        help="Leiden resolution for subclusters (default 0.6)")
    parser.add_argument("--cpu", action="store_true")
    args = parser.parse_args()

    print(f"\n=== Phase 7b: Subcluster '{args.celltype}' ===")
    cfg = load_config(args.config)
    seed = int(cfg.get("random_seed", 42))
    tissue = cfg["tissue"]
    n_hvg = int(cfg.get("integration", {}).get("n_hvg",
                2000 if tissue == "placenta" else 3000))

    # Input: prefer Phase 7c output, else Phase 7
    base = Path(cfg["results_dir"]) / "h5ad"
    candidates = [base / "08b_label_transferred" / "all_samples.h5ad",
                  base / "08_annotated" / "all_samples.h5ad"]
    in_path = next((p for p in candidates if p.is_file()), None)
    if in_path is None:
        sys.exit(f"ERROR: no annotated input found. Looked for:\n  " +
                 "\n  ".join(str(p) for p in candidates) +
                 "\n  Run 07_annotation.py (and optionally 07c) first.")
    print(f"  Input: {in_path}")

    slug = slugify(args.celltype)
    out_h5ad = base / "08c_subclustered"
    plot_dir = Path(cfg["results_dir"]) / "plots" / "07b_subcluster" / slug
    table_dir = Path(cfg["results_dir"]) / "tables"
    for d in (out_h5ad, plot_dir, table_dir):
        d.mkdir(parents=True, exist_ok=True)

    accelerator, precision = select_accelerator(force_cpu=args.cpu)

    print(f"\n[1/5] Loading + selecting cells...")
    adata = sc.read_h5ad(in_path)
    label_key = resolve_label_key(adata, args.label_key)

    available = adata.obs[label_key].astype(str).unique().tolist()
    if args.celltype not in available:
        sys.exit(
            f"ERROR: '{args.celltype}' not found in obs['{label_key}'].\n"
            f"  Available values: {sorted(available)}\n"
            f"  Pass an exact match (names are case-sensitive)."
        )

    mask = adata.obs[label_key].astype(str) == args.celltype
    sub = adata[mask].copy()
    print(f"  {sub.n_obs:,} cells labeled '{args.celltype}' (of {adata.n_obs:,})")
    MIN_CELLS = 50
    if sub.n_obs < MIN_CELLS:
        sys.exit(f"ERROR: only {sub.n_obs} cells — too few to subcluster reliably "
                 f"(need >= {MIN_CELLS}). Refusing to produce noise.")

    print(f"\n[2/5] Re-integrating subset (HVG + scVI)...")
    sub = reintegrate_subset(sub, seed, accelerator, precision, n_hvg)

    print(f"\n[3/5] Neighbors + Leiden (res={args.resolution}) + UMAP...")
    n_neighbors = min(15, max(5, sub.n_obs // 100))
    sc.pp.neighbors(sub, use_rep="X_scVI_sub", n_neighbors=n_neighbors, random_state=seed)
    sc.tl.leiden(sub, resolution=args.resolution, random_state=seed, key_added="subcluster")
    sc.tl.umap(sub, random_state=seed)
    n_sub = sub.obs["subcluster"].nunique()
    print(f"  {n_sub} subclusters")

    print(f"\n[4/5] Subcluster markers (on lognorm)...")
    add_lognorm(sub)
    sc.tl.rank_genes_groups(sub, groupby="subcluster", method="wilcoxon",
                            layer="lognorm", use_raw=False, key_added="sub_markers")
    markers = sc.get.rank_genes_groups_df(sub, group=None, key="sub_markers")
    top = (markers.sort_values("scores", ascending=False)
                  .groupby("group").head(20).reset_index(drop=True))
    top.to_csv(table_dir / f"subcluster_{slug}_markers.csv", index=False)

    print(f"\n[5/5] Composition + plots...")
    # Composition: subcluster × group and × sample (counts + fractions)
    comp = pd.crosstab(sub.obs["subcluster"], sub.obs["group"])
    comp_frac = comp.div(comp.sum(axis=1), axis=0)
    comp_out = comp.add_suffix("_n").join(comp_frac.add_suffix("_frac"))
    comp_out.to_csv(table_dir / f"subcluster_{slug}_composition.csv")

    def umap(color, title, fname, **kw):
        if color not in sub.obs.columns:
            return
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        sc.pl.umap(sub, color=color, ax=ax, show=False, frameon=False,
                   size=12, title=title, **kw)
        fig.tight_layout(); fig.savefig(plot_dir / fname, dpi=140, bbox_inches="tight")
        plt.close(fig)

    umap("subcluster", f"{args.celltype} subclusters (res={args.resolution})",
         "umap_subclusters.png", legend_loc="on data", legend_fontsize=7)
    umap("group", f"{args.celltype}: by stress group", "umap_by_group.png")
    umap("sample_id", f"{args.celltype}: by sample (batch check)", "umap_by_sample.png")

    top_genes = (top.groupby("group").head(5)["names"].unique().tolist())
    top_genes = [g for g in top_genes if g in sub.var_names]
    if top_genes:
        sub.obs["subcluster"] = sub.obs["subcluster"].astype("category")
        fig = sc.pl.dotplot(sub, var_names=top_genes, groupby="subcluster",
                            layer="lognorm", show=False, return_fig=True,
                            title=f"{args.celltype} subcluster markers")
        fig.savefig(plot_dir / "subcluster_markers_dotplot.png", dpi=130, bbox_inches="tight")
        plt.close()

    if "lognorm" in sub.layers:
        del sub.layers["lognorm"]
    sub.write_h5ad(out_h5ad / f"{slug}.h5ad")

    print(f"\n  Written: {out_h5ad / f'{slug}.h5ad'}")
    print(f"  Subclusters: {n_sub}  |  obs key: 'subcluster'  |  latent: 'X_scVI_sub'")
    print(f"  Plots: {plot_dir}")
    print(f"\n✓ Phase 7b complete for '{args.celltype}'.")
    print(f"  Re-run with a different --celltype for other lineages "
          f"(microglia, oligodendrocyte_lineage, ...).")
    print(f"  Check umap_by_sample.png — subclusters split by sample = batch, not biology.\n")


if __name__ == "__main__":
    main()
