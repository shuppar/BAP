#!/usr/bin/env python
"""
07c_label_transfer.py — Phase 7c: reference-based label transfer via scANVI.

Transfers cell type labels (and, where defensible, region labels) from a labeled
reference atlas onto the integrated query data, using scANVI.

Tissue-agnostic: all tissue-specific config comes from the YAML `reference:` block,
not hardcoded. Brain uses a spatially-resolved reference (e.g. ABC Atlas) and can
make regional claims; placenta has no spatial reference (region_key: null) and
gets cell-type labels only.

REGIONAL CLAIM POLICY (the key idea — see also project discussion):
  A region label is only transferred for cell types that are spatially RESTRICTED
  in the reference. For each reference cell type we compute the fraction of its
  cells in its single most common region; if that fraction >= the configured
  threshold (default 0.8), the type is "region-restricted" and its query cells
  get a region label. Otherwise the region label is SUPPRESSED (set to
  "region-ambiguous") — no fuzzy regional claims for broadly-distributed types
  like astrocytes, oligodendrocytes, microglia, vasculature.

  This split is DERIVED from the reference data, not hardcoded from memory.

scANVI workflow (verified against scvi-tools docs):
  1. Intersect genes between reference and query
  2. Concatenate; batch_key separates reference / query (and pools within query)
  3. Train scVI on the combined object
  4. SCANVI.from_scvi_model(..., unlabeled_category="Unknown", labels_key=...)
  5. scanvi.predict() → query cell type labels; .predict(soft=True) → confidence

Usage:
  uv run python scripts/07c_label_transfer.py --config config/brain.yaml
  uv run python scripts/07c_label_transfer.py --config config/dev.yaml
  uv run python scripts/07c_label_transfer.py --config config/placenta.yaml

Inputs:
  {results_dir}/h5ad/08_annotated/all_samples.h5ad   (from Phase 7)
  reference .h5ad path from the YAML `reference:` block

Outputs:
  {results_dir}/h5ad/08b_label_transferred/all_samples.h5ad
  {results_dir}/plots/07c_label_transfer/
    - umap_scanvi_celltype.png          : transferred cell type labels
    - umap_scanvi_confidence.png        : prediction confidence
    - umap_scanvi_region.png            : region labels (brain only)
  {results_dir}/tables/
    - scanvi_predictions.csv            : per-cell label + confidence (+ region)
    - region_restricted_celltypes.csv   : which types earned regional claims (brain)
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from _utils import load_config, add_lognorm


# ---------------------------------------------------------------------------
# Reference config validation
# ---------------------------------------------------------------------------

def read_reference_config(cfg: dict) -> dict:
    """Pull and validate the `reference:` block. Hard-fail on misconfiguration."""
    ref = cfg.get("reference")
    if ref is None:
        sys.exit(
            "ERROR: no 'reference:' block in config. Phase 7c needs one.\n"
            "  Re-run build_yaml.py (which now emits it) or add it by hand."
        )
    ref_h5ad = ref.get("ref_h5ad")
    if not ref_h5ad:
        sys.exit(
            "ERROR: reference.ref_h5ad is null. Set it to a labeled reference "
            ".h5ad path\n  (e.g. an ABC Atlas subset for brain) before running Phase 7c.\n"
            "  Until then, rely on the Phase 7 marker-based / CellTypist labels."
        )
    if not Path(ref_h5ad).is_file():
        sys.exit(f"ERROR: reference.ref_h5ad not found at: {ref_h5ad}")

    labels_key = ref.get("labels_key")
    if not labels_key:
        sys.exit("ERROR: reference.labels_key is required (cell type column in the reference).")

    return {
        "ref_h5ad": ref_h5ad,
        "labels_key": labels_key,
        "region_key": ref.get("region_key"),   # may be None → cell-type only
        "threshold": float(ref.get("region_concentration_threshold", 0.8)),
    }


# ---------------------------------------------------------------------------
# Region-restrictedness (derived from the reference, not hardcoded)
# ---------------------------------------------------------------------------

def compute_region_restricted_types(ref_adata, labels_key: str, region_key: str,
                                     threshold: float) -> pd.DataFrame:
    """For each reference cell type, find its dominant region and the fraction of
    its cells there. A type is 'region-restricted' if that fraction >= threshold.

    Returns a DataFrame: [cell_type, dominant_region, concentration, region_restricted].
    """
    if region_key not in ref_adata.obs.columns:
        sys.exit(
            f"ERROR: region_key '{region_key}' not in reference .obs. "
            f"Available: {list(ref_adata.obs.columns)}"
        )
    if labels_key not in ref_adata.obs.columns:
        sys.exit(
            f"ERROR: labels_key '{labels_key}' not in reference .obs. "
            f"Available: {list(ref_adata.obs.columns)}"
        )

    ct = pd.crosstab(ref_adata.obs[labels_key], ref_adata.obs[region_key])
    rows = []
    for cell_type in ct.index:
        counts = ct.loc[cell_type]
        total = counts.sum()
        if total == 0:
            continue
        dominant_region = counts.idxmax()
        concentration = float(counts.max() / total)
        rows.append({
            "cell_type": cell_type,
            "dominant_region": dominant_region,
            "concentration": round(concentration, 4),
            "n_ref_cells": int(total),
            "region_restricted": concentration >= threshold,
        })
    df = pd.DataFrame(rows).sort_values("concentration", ascending=False)
    return df


def assign_region_labels(query_obs, celltype_col: str,
                         restricted_df: pd.DataFrame) -> pd.Series:
    """Assign region labels to query cells based on their transferred cell type.

    A cell gets the dominant region of its type ONLY if that type is
    region-restricted; otherwise 'region-ambiguous'. Types not seen in the
    reference table get 'unknown'.
    """
    # Maps from cell type → region (only for restricted types)
    region_map = {
        r["cell_type"]: r["dominant_region"]
        for _, r in restricted_df.iterrows() if r["region_restricted"]
    }
    known_types = set(restricted_df["cell_type"])

    def lookup(ct):
        if ct in region_map:
            return region_map[ct]
        if ct in known_types:
            return "region-ambiguous"   # type exists in ref but is distributed
        return "unknown"                # type not in reference table
    return query_obs[celltype_col].map(lookup)


# ---------------------------------------------------------------------------
# scANVI label transfer
# ---------------------------------------------------------------------------

def run_scanvi_transfer(query, ref, labels_key: str, seed: int,
                        accelerator: str, precision: str):
    """Train scVI then scANVI on concatenated ref+query; predict query labels.

    Returns (query_with_predictions, predictions_df).
    API verified against scvi-tools docs (SCANVI.from_scvi_model workflow).
    """
    import anndata as ad
    import scvi

    UNLABELED = "Unknown"
    REF_BATCH = "__reference__"

    # --- Gene intersection (scANVI needs a shared feature space) ---
    common = query.var_names.intersection(ref.var_names)
    if len(common) < 200:
        sys.exit(
            f"ERROR: only {len(common)} genes shared between query and reference.\n"
            f"  Likely a gene-naming mismatch (symbols vs Ensembl IDs). scANVI needs\n"
            f"  a shared feature space; refusing to transfer on too few genes."
        )
    print(f"  Shared genes: {len(common):,}")
    query = query[:, common].copy()
    ref = ref[:, common].copy()

    # --- Label/batch columns ---
    # Reference carries true labels; query is all 'Unknown'.
    ref.obs["_scanvi_labels"] = ref.obs[labels_key].astype(str).values
    query.obs["_scanvi_labels"] = UNLABELED
    # batch_key separates reference from query so scANVI corrects the
    # technical/platform shift between them.
    ref.obs["_scanvi_batch"] = REF_BATCH
    query.obs["_scanvi_batch"] = query.obs["pool"].astype(str).values

    # --- Concatenate. Raw counts must be in .X for scVI. ---
    combined = ad.concat([ref, query], axis=0, join="inner",
                         label="_origin", keys=["ref", "query"],
                         index_unique="-")
    if combined.X is None:
        sys.exit("ERROR: combined .X is None — reference/query must have raw counts in .X.")

    print(f"  Combined: {combined.n_obs:,} cells "
          f"({ref.n_obs:,} ref + {query.n_obs:,} query)")

    scvi.settings.seed = seed

    # --- scVI ---
    scvi.model.SCVI.setup_anndata(combined, batch_key="_scanvi_batch")
    vae = scvi.model.SCVI(combined, n_layers=2, n_latent=30)
    max_epochs = 50 if combined.n_obs < 5000 else 200
    print(f"  Training scVI (max_epochs={max_epochs})...")
    vae.train(max_epochs=max_epochs, accelerator=accelerator,
              devices=1, precision=precision, early_stopping=True)

    # --- scANVI (verified API) ---
    print(f"  Training scANVI...")
    scanvi = scvi.model.SCANVI.from_scvi_model(
        vae, adata=combined,
        unlabeled_category="Unknown",
        labels_key="_scanvi_labels",
    )
    scanvi.train(max_epochs=min(max_epochs, 100), n_samples_per_label=100,
                 accelerator=accelerator, devices=1, precision=precision)

    # --- Predict on query cells only ---
    combined.obs["_pred"] = scanvi.predict(combined)
    combined.obsm["X_scANVI"] = scanvi.get_latent_representation(combined)
    # soft=True gives per-class probabilities; max = confidence
    proba = scanvi.predict(combined, soft=True)
    combined.obs["_pred_conf"] = proba.max(axis=1).values

    is_query = combined.obs["_origin"] == "query"
    q = combined[is_query].copy()

    # Strip the concat suffix to restore original query barcodes
    q.obs_names = [bc.rsplit("-query", 1)[0] for bc in q.obs_names]

    pred_df = pd.DataFrame({
        "scanvi_celltype": q.obs["_pred"].values,
        "scanvi_conf": q.obs["_pred_conf"].values,
    }, index=q.obs_names)

    return q, pred_df


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_umap(adata, color, title, out, cmap=None):
    if color not in adata.obs.columns:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    kw = {"color_map": cmap} if cmap else {}
    sc.pl.umap(adata, color=color, ax=ax, show=False, frameon=False,
               legend_fontsize=6, size=6, title=title, **kw)
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Phase 7c: scANVI label transfer")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--cpu", action="store_true", help="Force CPU")
    args = parser.parse_args()

    print(f"\n=== Phase 7c: Reference label transfer (scANVI) ===")
    print(f"Config: {args.config}")

    cfg = load_config(args.config)
    tissue = cfg["tissue"]
    seed = int(cfg.get("random_seed", 42))
    refcfg = read_reference_config(cfg)

    in_path = Path(cfg["results_dir"]) / "h5ad" / "08_annotated" / "all_samples.h5ad"
    if not in_path.is_file():
        sys.exit(f"ERROR: missing {in_path}. Run 07_annotation.py first.")

    out_dir = Path(cfg["results_dir"]) / "h5ad" / "08b_label_transferred"
    plot_dir = Path(cfg["results_dir"]) / "plots" / "07c_label_transfer"
    table_dir = Path(cfg["results_dir"]) / "tables"
    for d in (out_dir, plot_dir, table_dir):
        d.mkdir(parents=True, exist_ok=True)

    from _utils import select_accelerator
    accelerator, precision = select_accelerator(force_cpu=args.cpu)

    print(f"\n[1/5] Loading query: {in_path}")
    query = sc.read_h5ad(in_path)
    print(f"  Query: {query.n_obs:,} cells × {query.n_vars:,} genes")

    print(f"\n[2/5] Loading reference: {refcfg['ref_h5ad']}")
    ref = sc.read_h5ad(refcfg["ref_h5ad"])
    print(f"  Reference: {ref.n_obs:,} cells × {ref.n_vars:,} genes")
    if refcfg["labels_key"] not in ref.obs.columns:
        sys.exit(f"ERROR: labels_key '{refcfg['labels_key']}' not in reference .obs. "
                 f"Available: {list(ref.obs.columns)}")
    print(f"  Cell types in reference: {ref.obs[refcfg['labels_key']].nunique()}")

    # --- Region-restrictedness (only if region_key is set) ---
    do_region = refcfg["region_key"] is not None
    restricted_df = None
    if do_region:
        print(f"\n[3/5] Computing region-restricted cell types "
              f"(threshold={refcfg['threshold']})...")
        restricted_df = compute_region_restricted_types(
            ref, refcfg["labels_key"], refcfg["region_key"], refcfg["threshold"])
        restricted_df.to_csv(table_dir / "region_restricted_celltypes.csv", index=False)
        n_restricted = int(restricted_df["region_restricted"].sum())
        n_total = len(restricted_df)
        print(f"  {n_restricted}/{n_total} cell types are region-restricted "
              f"(>= {refcfg['threshold']:.0%} in one region):")
        for _, r in restricted_df[restricted_df["region_restricted"]].iterrows():
            print(f"    {r['cell_type']} → {r['dominant_region']} "
                  f"({r['concentration']:.0%})")
        print(f"  The other {n_total - n_restricted} types are distributed → "
              f"region label suppressed (no fuzzy claims).")
    else:
        print(f"\n[3/5] region_key is null ({tissue} has no spatial reference) — "
              f"cell-type transfer only, no regional claims.")

    # --- scANVI transfer ---
    print(f"\n[4/5] scANVI label transfer "
          f"(accelerator={accelerator}, precision={precision})...")
    q_pred, pred_df = run_scanvi_transfer(
        query, ref, refcfg["labels_key"], seed, accelerator, precision)

    # Align predictions back onto the query object (order-safe)
    pred_df = pred_df.reindex(query.obs_names)
    if pred_df["scanvi_celltype"].isna().any():
        n_missing = int(pred_df["scanvi_celltype"].isna().sum())
        sys.exit(f"ERROR: {n_missing} query cells got no prediction — barcode "
                 f"mismatch between transfer output and query. Aborting rather "
                 f"than writing partial labels.")
    query.obs["scanvi_celltype"] = pred_df["scanvi_celltype"].values
    query.obs["scanvi_conf"] = pred_df["scanvi_conf"].values

    # --- Region labels (only for restricted types) ---
    if do_region:
        query.obs["scanvi_region"] = assign_region_labels(
            query.obs, "scanvi_celltype", restricted_df).values
        n_regional = int((query.obs["scanvi_region"].isin(
            restricted_df.loc[restricted_df["region_restricted"], "dominant_region"]
        )).sum())
        print(f"  Regional labels assigned to {n_regional:,}/{query.n_obs:,} cells "
              f"(only region-restricted types).")

    # --- Save ---
    print(f"\n[5/5] Writing outputs + plots...")
    cols = ["scanvi_celltype", "scanvi_conf"] + (["scanvi_region"] if do_region else [])
    query.obs[cols].to_csv(table_dir / "scanvi_predictions.csv")

    plot_umap(query, "scanvi_celltype", "scANVI cell type (transferred)",
              plot_dir / "umap_scanvi_celltype.png")
    plot_umap(query, "scanvi_conf", "scANVI prediction confidence",
              plot_dir / "umap_scanvi_confidence.png", cmap="viridis")
    if do_region:
        plot_umap(query, "scanvi_region", "scANVI region (restricted types only)",
                  plot_dir / "umap_scanvi_region.png")

    query.write_h5ad(out_dir / "all_samples.h5ad")

    print(f"\n  Written: {out_dir / 'all_samples.h5ad'}")
    print(f"  obs columns added: {', '.join(cols)}")
    print(f"  Plots:  {plot_dir}")
    print(f"\n✓ Phase 7c complete.")
    if do_region:
        print(f"\n  Regional claims made ONLY for region-restricted types — see")
        print(f"  region_restricted_celltypes.csv for which types qualified.")
    else:
        print(f"\n  Cell-type labels only ({tissue} has no spatial reference).")
    print()


if __name__ == "__main__":
    main()
