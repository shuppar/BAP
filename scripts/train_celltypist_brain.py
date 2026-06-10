"""
train_celltypist_brain.py -- train CellTypist models from the ABC brain
reference, at TWO granularity tiers, and emit the mapping CSVs needed for
the derived columns.

Inputs (read from `reference:` block of config/brain.yaml):
  refs/abc_brain_ref.h5ad  -- labeled ABC reference (92K cells, prepared by
                              prepare_reference.py). Must have obs columns
                              `class`, `subclass`, `anatomical_division_label`.

Outputs:
  refs/celltypist_brain_adult_class.pkl     -- class-level model (34 categories)
  refs/celltypist_brain_adult_subclass.pkl  -- subclass-level model (~334 categories)
  refs/celltypist_brain_adult_region.pkl    -- region-level model on
                                              anatomical_division_label (~12 regions);
                                              gives per-cell region directly,
                                              avoiding the subclass->region collapse
                                              that buries cross-regional subclasses
                                              (e.g. cortical interneurons) in
                                              `multi_region`. Used at Phase 9.
  refs/abc_class_to_broad.csv               -- deterministic class -> broad map
  refs/abc_subclass_to_region.csv           -- subclass -> anatomical_division_label
                                              majority map (with fraction + n_cells);
                                              informational once the region model
                                              exists -- per-cell region comes from
                                              the model, not this CSV.

Design notes:
  - Trains on raw counts via CellTypist's standard pipeline
    (normalize -> log1p -> feature selection -> LogisticRegression).
  - GPU-backed via cuML (use_GPU=True). The first round of training (SGD-based
    feature selection) still runs on CPU; only the second-round LogReg fit is
    GPU-accelerated. Output classifiers are saved in sklearn-compatible format
    regardless of training backend.
  - n_jobs scales with available CPUs for the CPU prep + feature selection
    rounds. On workstation expect ~5-15 min for class, ~15-30 min for subclass
    (the CPU SGD round dominates; GPU LogReg fit is sub-minute), ~5-10 min for
    region. If skipping feature_selection becomes worthwhile, pass
    feature_selection=False inline (untested at scale).
  - Mapping CSVs are produced in the same script (single source of truth: the
    reference's own obs columns), so they can't drift from the trained models.

Idempotent: skips models/CSVs that already exist unless --force.

Usage:
  uv run python scripts/train_celltypist_brain.py --config config/brain.yaml
  # to rebuild from scratch:
  uv run python scripts/train_celltypist_brain.py --config config/brain.yaml --force
"""

import argparse
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml


# ---------------------------------------------------------------------------
# class -> broad collapse -- applied to ABC's 34 class labels.
# Each entry: ABC class name (exact string) -> broad biological class.
# Verified against Yao 2023 Nature Table 1 + the abc_atlas_access description
# (https://alleninstitute.github.io/abc_atlas_access/descriptions/WMB-taxonomy.html).
# Any ABC class not in this dict gets mapped to "Other" with a loud warning
# (catches reference-version drift).
# ---------------------------------------------------------------------------
CLASS_TO_BROAD = {
    # Cortical / pallial glutamatergic
    "01 IT-ET Glut":        "Excitatory neurons (cortical)",
    "02 NP-CT-L6b Glut":    "Excitatory neurons (cortical)",
    # Subpallial GABAergic
    "03 OB-CR Glut":        "Excitatory neurons (OB)",
    "04 DG-IMN Glut":       "Excitatory neurons (DG/IMN)",
    "05 OB-IMN GABA":       "Inhibitory neurons (OB)",
    "06 CTX-CGE GABA":      "Inhibitory neurons (cortical, CGE)",
    "07 CTX-MGE GABA":      "Inhibitory neurons (cortical, MGE)",
    "08 CNU-MGE GABA":      "Inhibitory neurons (CNU, MGE)",
    "09 CNU-LGE GABA":      "Inhibitory neurons (CNU, LGE)",
    "10 LSX GABA":          "Inhibitory neurons (LSX)",
    "11 CNU-HYa GABA":      "Inhibitory neurons (CNU/HYa)",
    "12 HY GABA":           "Inhibitory neurons (HY)",
    # Diencephalon / mesencephalon / hindbrain neurons
    "13 CNU-HYa Glut":      "Excitatory neurons (CNU/HYa)",
    "14 HY Glut":           "Excitatory neurons (HY)",
    "15 HY Gnrh1 Glut":     "Excitatory neurons (HY)",
    "16 HY MM Glut":        "Excitatory neurons (HY)",
    "17 MH-LH Glut":        "Excitatory neurons (MH/LH)",
    "18 TH Glut":           "Excitatory neurons (TH)",
    "19 MB Glut":           "Excitatory neurons (MB)",
    "20 MB GABA":           "Inhibitory neurons (MB)",
    "21 MB Dopa":           "Dopaminergic neurons",
    "22 MB-HB Sero":        "Serotonergic neurons",
    "23 P Glut":            "Excitatory neurons (P)",
    "24 MY Glut":           "Excitatory neurons (MY)",
    "25 Pineal Glut":       "Excitatory neurons (pineal)",
    "26 P GABA":            "Inhibitory neurons (P)",
    "27 MY GABA":           "Inhibitory neurons (MY)",
    "28 CB GABA":           "Inhibitory neurons (CB)",
    "29 CB Glut":           "Excitatory neurons (CB)",
    # Non-neuronal
    "30 Astro-Epen":        "Astrocytes/Ependymal",
    "31 OPC-Oligo":         "OPC/Oligodendrocytes",
    "32 OEC":               "Olfactory ensheathing cells",
    "33 Vascular":          "Vascular",
    "34 Immune":            "Immune",
}


def log(msg: str) -> None:
    print(f"[train] {msg}", flush=True)


def train_one_model(
    adata: ad.AnnData,
    labels_key: str,
    out_path: Path,
    n_jobs: int,
) -> None:
    """Train a CellTypist model at the given labels_key and save to out_path."""
    import celltypist  # heavy import -- defer until needed
    log(f"Training CellTypist at labels_key='{labels_key}' "
        f"({adata.obs[labels_key].nunique()} categories, "
        f"{adata.n_obs:,} cells)...")
    t0 = time.time()
    model = celltypist.train(
        adata,
        labels=labels_key,
        n_jobs=n_jobs,
        check_expression=False,   # reference already validated upstream
        # GPU path via cuML LogisticRegression. Confirmed working on RTX 4500 Ada
        # 2026-06-10. CPU LogReg was burning >18 hours on subclass tier (334
        # classes, lbfgs not parallelizable across classes); cuML does the same
        # fit in minutes. Output classifier is saved in sklearn-compatible
        # format (CellTypist's design) regardless of training backend, so the
        # produced .pkl files are interchangeable with CPU-trained pkls.
        use_GPU=True,
        feature_selection=True,
        top_genes=300,
    )
    elapsed = time.time() - t0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    model.write(str(out_path))
    log(f"  OK {labels_key} model written to {out_path} ({elapsed/60:.1f} min)")


def emit_class_to_broad_csv(adata: ad.AnnData, out_path: Path) -> None:
    """Emit the deterministic class->broad mapping with n_cells per class.

    Hard-fail if the reference has any class label we haven't categorized
    (catches drift between this script's hard-coded map and the ABC version
    actually loaded).
    """
    classes_in_ref = sorted(adata.obs["class"].astype(str).unique())
    unmapped = [c for c in classes_in_ref if c not in CLASS_TO_BROAD]
    if unmapped:
        sys.exit(
            f"ERROR: {len(unmapped)} class label(s) in the reference have no "
            f"entry in CLASS_TO_BROAD:\n  " + "\n  ".join(unmapped) +
            f"\nUpdate CLASS_TO_BROAD in train_celltypist_brain.py."
        )
    counts = adata.obs["class"].value_counts()
    rows = [
        {"class": c, "broad": CLASS_TO_BROAD[c], "n_cells_in_ref": int(counts.get(c, 0))}
        for c in classes_in_ref
    ]
    df = pd.DataFrame(rows).sort_values(["broad", "class"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    n_broad = df["broad"].nunique()
    log(f"  OK class->broad mapping: {len(df)} classes -> {n_broad} broad -> {out_path}")


def emit_subclass_to_region_csv(
    adata: ad.AnnData,
    out_path: Path,
    purity_threshold: float = 0.70,
) -> None:
    """Per-subclass majority map of anatomical_division_label.

    For each subclass, find the single most common anatomical_division. If that
    division contains >=purity_threshold of the subclass's cells, assign it.
    Otherwise label `multi_region` so the gap is explicit downstream.
    """
    region_col = "anatomical_division_label"
    if region_col not in adata.obs.columns:
        sys.exit(f"ERROR: '{region_col}' missing from reference .obs")

    sx_region = (
        adata.obs.groupby("subclass", observed=True)[region_col]
        .value_counts(normalize=False)
        .rename("n_cells")
        .reset_index()
    )
    rows = []
    for subclass, grp in sx_region.groupby("subclass", observed=True):
        grp_sorted = grp.sort_values("n_cells", ascending=False)
        total = grp_sorted["n_cells"].sum()
        top_div = grp_sorted.iloc[0][region_col]
        top_frac = grp_sorted.iloc[0]["n_cells"] / total
        if top_frac >= purity_threshold:
            region = top_div
        else:
            region = "multi_region"
        rows.append({
            "subclass": subclass,
            "region": region,
            "top_division": top_div,
            "top_division_fraction": round(top_frac, 3),
            "n_cells_in_ref": int(total),
        })
    df = pd.DataFrame(rows).sort_values("subclass")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    n_multi = (df["region"] == "multi_region").sum()
    log(f"  OK subclass->region mapping: {len(df)} subclasses "
        f"({n_multi} flagged multi_region at purity<{purity_threshold}) -> {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path,
                        help="config/brain.yaml")
    parser.add_argument("--force", action="store_true",
                        help="Retrain/rebuild even if outputs exist")
    parser.add_argument("--n-jobs", type=int, default=-1, dest="n_jobs",
                        help="Parallel jobs for CellTypist (default -1 = all cores)")
    args = parser.parse_args()

    cfg = yaml.safe_load(args.config.read_text())
    if cfg.get("tissue") != "brain":
        sys.exit(f"ERROR: this script is for brain only (tissue={cfg.get('tissue')})")

    ref_h5ad = Path(cfg["reference"]["ref_h5ad"])
    if not ref_h5ad.exists():
        sys.exit(f"ERROR: reference not found: {ref_h5ad}\n"
                 f"Run prepare_reference.py first.")

    out_dir = ref_h5ad.parent
    paths = {
        "class_pkl":          out_dir / "celltypist_brain_adult_class.pkl",
        "subclass_pkl":       out_dir / "celltypist_brain_adult_subclass.pkl",
        "region_pkl":         out_dir / "celltypist_brain_adult_region.pkl",
        "class_to_broad":     out_dir / "abc_class_to_broad.csv",
        "subclass_to_region": out_dir / "abc_subclass_to_region.csv",
    }

    log(f"=== Train CellTypist brain models ===")
    log(f"Reference: {ref_h5ad}")
    log(f"Outputs:")
    for k, p in paths.items():
        exists = "exists" if p.exists() else "to-build"
        log(f"  {k:20s} -> {p}  [{exists}]")

    # Load once, reused for all four outputs
    log(f"\nLoading reference...")
    t0 = time.time()
    adata = ad.read_h5ad(ref_h5ad)
    log(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes "
        f"({time.time()-t0:.1f}s)")

    for required in ("class", "subclass", "anatomical_division_label"):
        if required not in adata.obs.columns:
            sys.exit(f"ERROR: reference missing required obs column '{required}'")

    # ---- Mapping CSVs first (cheap, surface drift early) ----
    log(f"\n[1/5] Emitting class->broad mapping...")
    if paths["class_to_broad"].exists() and not args.force:
        log(f"  (exists, skipping -- pass --force to rebuild)")
    else:
        emit_class_to_broad_csv(adata, paths["class_to_broad"])

    log(f"\n[2/5] Emitting subclass->region mapping...")
    if paths["subclass_to_region"].exists() and not args.force:
        log(f"  (exists, skipping -- pass --force to rebuild)")
    else:
        emit_subclass_to_region_csv(adata, paths["subclass_to_region"])

    # ---- Class model ----
    log(f"\n[3/5] Training class model...")
    if paths["class_pkl"].exists() and not args.force:
        log(f"  (exists, skipping -- pass --force to rebuild)")
    else:
        train_one_model(adata, "class", paths["class_pkl"], n_jobs=args.n_jobs)

    # ---- Subclass model ----
    log(f"\n[4/5] Training subclass model...")
    if paths["subclass_pkl"].exists() and not args.force:
        log(f"  (exists, skipping -- pass --force to rebuild)")
    else:
        train_one_model(adata, "subclass", paths["subclass_pkl"], n_jobs=args.n_jobs)

    # ---- Region model ----
    # Trained on anatomical_division_label (~12 ABC regions: Isocortex, MB, MY,
    # HY, P, STR, HPF, TH, OLF, CB, CTXsp, PAL, ...). Per-cell region from this
    # model is the canonical Phase-9 region -- context-aware (e.g. a Pvalb-GABA
    # cell in cortical context predicts Isocortex), unlike the per-subclass
    # collapse which puts cross-regional subclasses in `multi_region`.
    log(f"\n[5/5] Training region model...")
    if paths["region_pkl"].exists() and not args.force:
        log(f"  (exists, skipping -- pass --force to rebuild)")
    else:
        train_one_model(adata, "anatomical_division_label",
                        paths["region_pkl"], n_jobs=args.n_jobs)

    log(f"\nOK Done. Update brain.yaml so Phase 7 picks up all three pkls:")
    log(f"  annotation:")
    log(f"    celltypist_models:")
    log(f"      P1:")
    log(f"        class:    Developing_Mouse_Brain.pkl   # broad-ish, derive class downstream")
    log(f"        subclass: null                          # not available for P1")
    log(f"        region:   null                          # not meaningful at P1 (regions not yet defined)")
    log(f"      4W:")
    log(f"        class:    refs/celltypist_brain_adult_class.pkl")
    log(f"        subclass: refs/celltypist_brain_adult_subclass.pkl")
    log(f"        region:   refs/celltypist_brain_adult_region.pkl")
    log(f"      3mo:")
    log(f"        class:    refs/celltypist_brain_adult_class.pkl")
    log(f"        subclass: refs/celltypist_brain_adult_subclass.pkl")
    log(f"        region:   refs/celltypist_brain_adult_region.pkl")
    log(f"    class_to_broad_csv:    refs/abc_class_to_broad.csv")
    log(f"    subclass_to_region_csv: refs/abc_subclass_to_region.csv  # informational; per-cell region comes from the region model")


if __name__ == "__main__":
    main()
