#!/usr/bin/env python
"""
replot_brain_annotation.py — regenerate the brain broad + per-age class UMAPs
from the (patched) annotated object, without re-running Phase 7.

Fixes two stale/ugly figures:
  - umap_celltypist_broad.png : was 29 region-tagged classes (pre-patch) +
    squished; now 9 region-free classes, widened so the legend doesn't compress
    the plot.
  - umap_celltypist_class_by_age.png : per-age panels had NO legend; now uses
    on-data labels so each panel is self-documenting.

Reads the saved UMAP coords already in the object (.obsm['X_umap']) — no
recompute. Fast (~1 min, mostly the h5ad read).

Usage:
  uv run python scripts/replot_brain_annotation.py \
      --h5ad results/brain/h5ad/08_annotated/all_samples.h5ad \
      --plot-dir results/brain/plots/07_annotation
"""

import argparse
from pathlib import Path

import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True, type=Path)
    ap.add_argument("--plot-dir", required=True, type=Path)
    args = ap.parse_args()

    a = sc.read_h5ad(args.h5ad)
    args.plot_dir.mkdir(parents=True, exist_ok=True)

    if "X_umap" not in a.obsm:
        raise SystemExit("ERROR: no X_umap in object — can't replot.")

    # --- broad UMAP (region-free, widened) ---
    if "celltypist_broad" in a.obs.columns:
        n = a.obs["celltypist_broad"].nunique()
        fig, ax = plt.subplots(figsize=(11, 6))
        sc.pl.umap(a, color="celltypist_broad", ax=ax, show=False, frameon=False,
                   size=6, legend_loc="right margin", legend_fontsize=8,
                   title=f"celltypist_broad ({n} classes; all ages aligned)")
        fig.tight_layout()
        out = args.plot_dir / "umap_celltypist_broad.png"
        fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
        print(f"  wrote {out}  ({n} region-free broad classes)")

    # --- per-age class UMAP (shared legend) ---
    if "celltypist_class" in a.obs.columns:
        a.obs["celltypist_class"] = a.obs["celltypist_class"].astype("category")
        cats = list(a.obs["celltypist_class"].cat.categories)
        # stable color per category (shared across panels). Use matplotlib
        # colormaps (guaranteed present) rather than scanpy palette attrs.
        import matplotlib as mpl
        import matplotlib.cm as cm
        base = cm.get_cmap("tab20", 20)
        if len(cats) <= 20:
            color_map = {c: base(i) for i, c in enumerate(cats)}
        else:
            big = cm.get_cmap("gist_ncar", len(cats))
            color_map = {c: big(i) for i, c in enumerate(cats)}

        age_col = a.obs["age"].astype(str)
        ages = sorted(age_col.unique())
        palette_list = [mpl.colors.to_hex(color_map[c]) for c in cats]
        ncols = len(ages)
        fig, axes = plt.subplots(1, ncols, figsize=(7 * ncols, 7))
        if ncols == 1:
            axes = [axes]
        for ax, age in zip(axes, ages):
            sub = a[age_col == age].copy()
            # keep the FULL category set so colors stay consistent per panel
            sub.obs["celltypist_class"] = sub.obs["celltypist_class"].cat.set_categories(cats)
            sc.pl.umap(sub, color="celltypist_class", ax=ax, show=False,
                       frameon=False, size=8, legend_loc=None,
                       palette=palette_list,
                       title=f"celltypist_class ({age})  n={sub.n_obs:,}")
        # one shared legend on the right
        handles = [mpl.lines.Line2D([0], [0], marker="o", linestyle="",
                                    markersize=6, markerfacecolor=color_map[c],
                                    markeredgewidth=0, label=c) for c in cats]
        fig.legend(handles=handles, loc="center left",
                   bbox_to_anchor=(1.0, 0.5), fontsize=7, frameon=False,
                   ncol=2 if len(cats) > 25 else 1, title="celltypist_class")
        fig.tight_layout()
        out = args.plot_dir / "umap_celltypist_class_by_age.png"
        fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
        print(f"  wrote {out}  (shared legend, {len(ages)} age panels, {len(cats)} classes)")


if __name__ == "__main__":
    main()
