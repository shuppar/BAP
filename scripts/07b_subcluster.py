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
  2. manual_annotation      (if filled in and non-empty)
  3. celltypist_broad       (brain Phase 7 output)
  4. celltype_majority      (placenta Phase 7 output)
  5. scanvi_celltype        (Phase 7c, if present)
  6. celltypist_majority    (legacy)
  7. provisional_celltype   (marker-score fallback)

The chosen column and the exact --celltype value must match; the script lists
available labels and exits if the value isn't found (no silent empty subset).

Usage:
  uv run python scripts/07b_subcluster.py --config config/brain.yaml \
      --celltype "Immune"
  uv run python scripts/07b_subcluster.py --config config/placenta.yaml \
      --celltype "Myeloid" --resolution 0.5
  # explicit label column override:
  uv run python scripts/07b_subcluster.py --config config/brain.yaml \
      --celltype "OPC/Oligodendrocytes" --label-key celltypist_broad --resolution 0.6

Inputs:
  {results_dir}/h5ad/08_annotated/all_samples.h5ad            (Phase 7)

Outputs (slug = sanitized celltype name):
  {results_dir}/h5ad/08c_subclustered/{slug}.h5ad
  {results_dir}/plots/07b_subcluster/{slug}/
    - umap_subclusters.png
    - umap_by_group.png
    - umap_by_age.png
    - umap_by_sample.png
    - subcluster_markers_dotplot.png
  {results_dir}/tables/07b_subcluster/
    - 07b_subcluster_{slug}_markers.csv
    - 07b_subcluster_{slug}_composition.csv
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

from _utils import load_config, add_lognorm, select_accelerator, phase_table_dir


# Priority order for auto-detecting the label column.
# celltypist_broad (brain) and celltype_majority (placenta) added first
# after manual_annotation so they are found before the legacy keys.
LABEL_KEY_PRIORITY = [
    "manual_annotation",
    "celltypist_broad",       # brain Phase 7
    "celltype_majority",      # placenta Phase 7
    "scanvi_celltype",        # Phase 7c if present
    "celltypist_majority",    # legacy
    "provisional_celltype",   # marker-score fallback
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

    # Verify .X looks like raw counts — scVI requires raw.
    Xmax = adata.X.max()
    if not np.isclose(Xmax, np.round(Xmax)) or Xmax < 0:
        sys.exit("ERROR: subset .X does not look like raw counts (scVI needs raw). "
                 f"max={Xmax}. Aborting rather than training on the wrong matrix.")

    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        print(f"  [info] only {n_batches} '{batch_key}' level in this subset — "
              f"no batch correction needed within it.")

    # HVG on the subset (seurat_v3 on raw counts). Cap at n_vars.
    n_hvg = min(n_hvg, adata.n_vars - 1)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_hvg, flavor="seurat_v3",
        batch_key=batch_key if n_batches >= 2 else None,
    )
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    scvi.settings.seed = seed
    setup_kwargs = {} if n_batches < 2 else {"batch_key": batch_key}
    scvi.model.SCVI.setup_anndata(adata_hvg, **setup_kwargs)
    model = scvi.model.SCVI(adata_hvg, n_layers=2, n_latent=30)
    # Small subsets converge fast; cap at 50 epochs below 5k cells.
    max_epochs = 50 if adata.n_obs < 5000 else 200
    print(f"  Training scVI on subset (n={adata.n_obs:,}, max_epochs={max_epochs})...")
    model.train(
        max_epochs=max_epochs, accelerator=accelerator, devices=1,
        precision=precision, early_stopping=True,
    )

    adata.obsm["X_scVI_sub"] = model.get_latent_representation()
    return adata


def main():
    parser = argparse.ArgumentParser(description="Phase 7b: subcluster one cell type")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--celltype", required=True,
                        help="Cell type value to subcluster (exact match to label column)")
    parser.add_argument("--label-key", default=None,
                        help="obs column holding cell type labels (default: auto-detect)")
    parser.add_argument("--resolution", type=float, default=0.6,
                        help="Leiden resolution for subclusters (default 0.6)")
    parser.add_argument("--cpu", action="store_true", help="Force CPU (no GPU)")
    args = parser.parse_args()

    print(f"\n=== Phase 7b: Subcluster '{args.celltype}' ===")
    cfg = load_config(args.config)
    seed = int(cfg.get("random_seed", 42))
    tissue = cfg["tissue"]
    n_hvg = int(cfg.get("integration", {}).get("n_hvg",
                2000 if tissue == "placenta" else 3000))

    # Input: prefer Phase 7c output, else Phase 7
    base = Path(cfg["results_dir"]) / "h5ad"
    candidates = [
        base / "08_annotated" / "all_samples.h5ad",
    ]
    in_path = next((p for p in candidates if p.is_file()), None)
    if in_path is None:
        sys.exit(
            "ERROR: no annotated input found. Looked for:\n  "
            + "\n  ".join(str(p) for p in candidates)
            + "\n  Run 07_annotation.py first."
        )
    print(f"  Input:   {in_path}")

    slug = slugify(args.celltype)
    out_h5ad = base / "08c_subclustered"
    plot_dir = Path(cfg["results_dir"]) / "plots" / "07b_subcluster" / slug
    table_dir = phase_table_dir(cfg, "07b_subcluster")
    for d in (out_h5ad, plot_dir, table_dir):
        d.mkdir(parents=True, exist_ok=True)

    accelerator, precision = select_accelerator(force_cpu=args.cpu)

    # ------------------------------------------------------------------ #
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
    MIN_CELLS = 2000  # agreed policy: don't subcluster below this — too few for stable scVI subtypes
    if sub.n_obs < MIN_CELLS:
        sys.exit(f"ERROR: only {sub.n_obs} cells — too few to subcluster reliably "
                 f"(need >= {MIN_CELLS}). Refusing to produce noise.")

    # ------------------------------------------------------------------ #
    print(f"\n[2/5] Re-integrating subset (HVG + scVI)...")
    sub = reintegrate_subset(sub, seed, accelerator, precision, n_hvg)

    # ------------------------------------------------------------------ #
    print(f"\n[3/5] Neighbors + Leiden (res={args.resolution}) + UMAP...")
    n_neighbors = min(15, max(5, sub.n_obs // 100))
    sc.pp.neighbors(sub, use_rep="X_scVI_sub", n_neighbors=n_neighbors,
                    random_state=seed)
    sc.tl.leiden(sub, resolution=args.resolution, random_state=seed,
                 key_added="subcluster", flavor="igraph",
                 n_iterations=2, directed=False)
    sc.tl.umap(sub, random_state=seed)
    n_sub = sub.obs["subcluster"].nunique()
    print(f"  {n_sub} subclusters found")

    # ------------------------------------------------------------------ #
    print(f"\n[4/5] Subcluster markers (on lognorm)...")
    add_lognorm(sub)
    sc.tl.rank_genes_groups(sub, groupby="subcluster", method="wilcoxon",
                            layer="lognorm", use_raw=False, key_added="sub_markers")
    markers = sc.get.rank_genes_groups_df(sub, group=None, key="sub_markers")
    top = (markers.sort_values("scores", ascending=False)
                  .groupby("group").head(20).reset_index(drop=True))
    # Map Ensembl var_names → symbols if needed
    if "symbol" in sub.var.columns:
        sym_map = sub.var["symbol"].to_dict()
        top["symbol"] = top["names"].map(sym_map).fillna(top["names"])
    top.to_csv(table_dir / f"07b_subcluster_{slug}_markers.csv", index=False)

    # ------------------------------------------------------------------ #
    print(f"\n[5/5] Composition + plots...")
    comp = pd.crosstab(sub.obs["subcluster"], sub.obs["group"])
    comp_frac = comp.div(comp.sum(axis=1), axis=0)
    comp_out = comp.add_suffix("_n").join(comp_frac.add_suffix("_frac"))
    comp_out.to_csv(table_dir / f"07b_subcluster_{slug}_composition.csv")

    # Ensure subcluster is categorical for on-data labels
    sub.obs["subcluster"] = sub.obs["subcluster"].astype("category")

    def save_umap(color, title, fname, **kw):
        if color not in sub.obs.columns:
            print(f"  [skip] '{color}' not in obs — skipping {fname}")
            return
        # Cast to category so on-data labels work
        if sub.obs[color].dtype.name != "category":
            sub.obs[color] = sub.obs[color].astype("category")
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        sc.pl.umap(sub, color=color, ax=ax, show=False, frameon=False,
                   size=12, title=title, **kw)
        fig.tight_layout()
        fig.savefig(plot_dir / fname, dpi=140, bbox_inches="tight")
        plt.close(fig)

    save_umap("subcluster",
              f"{args.celltype} subclusters (res={args.resolution})",
              "umap_subclusters.png",
              legend_loc="on data", legend_fontsize=7)
    save_umap("group",
              f"{args.celltype}: by stress group", "umap_by_group.png")
    save_umap("age",
              f"{args.celltype}: by age", "umap_by_age.png")
    save_umap("sample_id",
              f"{args.celltype}: by sample (batch check)", "umap_by_sample.png")

    # Dotplot: top 5 marker genes per subcluster.
    # Keep Ensembl IDs as var_names to avoid duplicate-symbol crash;
    # pass gene_symbols="symbol" so scanpy labels axes with readable names.
    top_genes = top.groupby("group", observed=True).head(5)["names"].unique().tolist()
    top_genes = [g for g in top_genes if g in sub.var_names]
    if top_genes:
        dotplot_kwargs = {"gene_symbols": "symbol"} if "symbol" in sub.var.columns else {}
        fig = sc.pl.dotplot(
            sub, var_names=top_genes, groupby="subcluster",
            layer="lognorm", use_raw=False, show=False, return_fig=True,
            title=f"{args.celltype} subcluster markers (top 5/cluster)",
            **dotplot_kwargs,
        )
        fig.savefig(plot_dir / "subcluster_markers_dotplot.png",
                    dpi=130, bbox_inches="tight")
        plt.close()

    # Drop lognorm layer before saving (raw counts only in .X)
    if "lognorm" in sub.layers:
        del sub.layers["lognorm"]
    out_path = out_h5ad / f"{slug}.h5ad"
    sub.write_h5ad(out_path)

    print(f"\n  Written:    {out_path}")
    print(f"  Subclusters: {n_sub}  |  obs key: 'subcluster'  |  latent: 'X_scVI_sub'")
    print(f"  Plots:      {plot_dir}")
    print(f"  Tables:     {table_dir}")
    print(f"\n✓ Phase 7b complete for '{args.celltype}'.")
    print(f"  Check umap_by_sample.png — subclusters split by sample = batch, not biology.\n")


if __name__ == "__main__":
    main()
