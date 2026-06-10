#!/usr/bin/env python
"""
07_annotation.py — Phase 7: cell type annotation (per-tissue primary track + marker cross-check).

BRAIN primary track: 3-tier CellTypist per age, locked by INSTRUCTIONS:

  class    -> canonical per-(Leiden cluster x age) label (celltypist_class).
              Assigned by MAJORITY VOTE of per-cell class predictions WITHIN
              each (Leiden cluster x age) group. Same value for all cells in
              the (cluster, age) group. Per-age native vocabulary:
                P1     -> Di Bella developing-brain labels
                4W/3mo -> ABC class labels (34 categories)
              8b/8c key off this column.

  subclass -> kept PER-CELL (celltypist_subclass), 4W/3mo only.
              P1 cells get the sentinel "no_subclass_model".
              Consumed at the subcluster level by 7b/7d.

  region   -> kept PER-CELL (celltypist_region), 4W/3mo only.
              P1 cells get "no_region_model".
              Consumed at Phase 9 for region matching against human atlases.

PLACENTA primary track: STAMP reference-based Spearman correlation
(Liu et al. 2024, eLife; 35 cell types covering E9.5-E18.5 incl. maternal
compartment). For each cell, log2FC-rank-correlate the per-cell expression
against each of 35 STAMP cell-type signatures; per-cluster majority -> label.
Then apply:
  - Tier 1+2 label-collapse map (35 -> ~21 labels) — same-lineage sub-states
    that can't be separated at main-pass Leiden resolution are merged.
  - Lineage-aware low-confidence flag: only fires when runner-up is a
    DIFFERENT lineage (not a sister state).
  - STRICT canonical-marker gates for noisy labels (Neutrophil / Lymphoid /
    Megakaryocyte):
        Neutrophil  requires Ly6g + ≥1 of (S100a8/S100a9/Csf3r); else demoted
                    to unassigned_immune.
        Lymphoid    requires ≥2 of (Ncr1/Klrb1c/Gzma/Gzmb/Cd3d/Cd3e/Cd19/Ms4a1);
                    else demoted to unassigned_immune.
        Megakaryocyte requires ≥2 of (Pf4/Itga2b/Gp1ba/Gp9/Ppbp); else
                    demoted to unassigned_blood.
  - Lymphoid sub-split after passing the gate: ≥2 markers per sub-class
    splits Lymphoid into NK / T cell / B cell. Multiple sub-classes →
    Lymphoid_mixed. Top gate passes but no sub-class clean → Lymphoid.
  - Xist + Y-gene compartment scoring. NOTE: Xist is absent from the 10x
    Flex probe panel — only fetal_male disambiguation (Y+ Xist-low) works
    in practice; maternal vs fetal-female cannot be resolved here.
  - EPC/TSC negative-control QC: Cdx2/Eomes/Elf5/Esrrb should NOT score
    positive at E12.5/E18.5; positive = contamination flag.
8b/8c key off celltype_majority (placenta) the same way they key off
celltypist_class (brain).

The marker track (rank_genes_groups + curated dotplot + score_genes UMAPs)
always runs for both tissues as a cross-check.

Outputs:
  {results_dir}/h5ad/08_annotated/all_samples.h5ad

  {results_dir}/plots/07_annotation/
    - umap_leiden_for_annotation.png                     (both tissues)
    - marker_dotplot.png                                  (both tissues)
    - marker_heatmap_top10.png                            (both tissues)
    - umap_marker_scores.png                              (both tissues)
    - cluster_composition_by_sample.png                   (both tissues)
    - celltype_composition_by_sample.png                  (both tissues)
    - celltype_composition_by_group.png                   (both tissues)
    BRAIN-only:
    - umap_celltypist_class.png                           (headline; per cluster x age majority)
    - umap_celltypist_class_confidence.png                (per-cell max-prob from class model)
    - umap_celltypist_subclass_by_age.png                 (faceted; 4W/3mo only)
    - umap_celltypist_region_by_age.png                   (faceted; 4W/3mo only)
    PLACENTA-only:
    - umap_celltype_majority.png                          (headline; STAMP-collapsed labels)
    - umap_celltype_gap.png                               (per-cell Spearman gap = confidence)
    - umap_celltype_compartment.png                       (Y+ / Xist+ / ambiguous)
    - cluster_fractions_heatmap.png                       (cluster x STAMP-collapsed type)
    - sanity_check_summary.png                            (per-cluster age + EPC QC)

  {results_dir}/tables/07_annotation/
    BRAIN:
    - 07_annotation_class_per_cluster_age.csv             *** the audit CSV ***
        rows = (leiden, age); cols = winner_class, purity, top1/2/3 label+frac,
        n_cells, low_purity flag (purity<0.6 OR runner-up>0.25),
        markers_checked, markers_present, gate_outcome
            (no_gate | passed | demoted), gate_label
    - 07_annotation_age_composition_sanity.csv            (developmentally
        implausible label-at-age combos -- e.g. erythrocyte at 4W/3mo, IPC
        outside P1. INFORMATIONAL only, does not modify labels.)
    - 07_annotation_celltypist_predictions.csv            (per-cell raw predictions)
    PLACENTA:
    - 07_annotation_cluster_summary.csv                   (per-cluster with flags)
    - 07_annotation_cluster_fractions.csv                 (top-5 STAMP-fine types per cluster)
    - 07_annotation_sanity_check.csv                      (age composition + EPC QC per cluster)
    BOTH:
    - 07_annotation_marker_genes_per_cluster.csv
    - 07_annotation_summary.csv                           (per-Leiden-cluster pooled summary)
    - 07_annotation_celltype_composition.csv              (per-sample × celltype fractions)
    - 07_annotation_celltype_by_age_wide.csv              (celltype × age mean/sem; one row per type)
    - 07_annotation_celltype_by_age_long.csv              (tidy: age × celltype with mean/sem/min/max/n_samples)

Obs columns written:
  BRAIN:
    celltypist_class_predicted, celltypist_class_conf     (per-cell raw)
    celltypist_class,           celltypist_class_purity   (per cluster x age majority — CANONICAL)
    celltypist_subclass,        celltypist_subclass_conf  (per-cell raw; "no_subclass_model" for P1)
    celltypist_region,          celltypist_region_conf    (per-cell raw; "no_region_model" for P1)
  PLACENTA:
    provisional_celltype_fine                             (per-cell STAMP fine label)
    provisional_celltype                                  (per-cell STAMP collapsed label)
    celltype_majority                                     (per-cluster collapsed label — CANONICAL)
    celltype_majority_flag                                (per-cluster flag string)
    celltype_gap                                          (per-cell Spearman gap)
    celltype_top_corr                                     (per-cell top Spearman correlation)
    compartment                                           (maternal / fetal_male / xist_positive / ambiguous)
  BOTH:
    provisional_celltype                                  (marker-based fallback; ALWAYS computed)
    manual_annotation                                     (empty string; for downstream editing)

Usage:
  uv run python scripts/07_annotation.py --config config/brain.yaml
  uv run python scripts/07_annotation.py --config config/placenta.yaml
  uv run python scripts/07_annotation.py --config config/dev.yaml
"""

import argparse
import sys
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import rankdata

from _utils import load_config, add_lognorm, phase_table_dir


# Sentinels — explicit strings rather than NaN so they never surface as "nan"
# categories in plots/tables and they're easy to filter on downstream.
NO_SUBCLASS = "no_subclass_model"
NO_REGION   = "no_region_model"

# Per-(cluster x age) audit-CSV low-purity criteria (BRAIN).
LOW_PURITY_THRESHOLD  = 0.60   # winner fraction below this -> flag
RUNNER_UP_THRESHOLD   = 0.25   # any runner-up above this   -> flag (cluster may need splitting)

# ---------------------------------------------------------------------------
# BRAIN -- STRICT canonical-marker gates for borderline CellTypist calls.
# Mirrors the placenta STAMP gates (Neutrophil / Lymphoid / Megakaryocyte).
# CellTypist's calibrated conf says "how sure is the LogReg among 34 trained
# classes" -- it cannot independently verify the cell expresses the biology
# the label requires. A cluster can score high CellTypist conf for a glial
# label while expressing only the generic gene-weight pattern, not the
# canonical markers. This gate catches that.
#
# Gating logic: if the winner_class label matches any of `keywords` (case-
# insensitive substring), check the cluster's marker expression. A marker is
# "present" if a fraction >= MARKER_PRESENCE_THRESHOLD of cells in that
# (cluster, age) group have lognorm > 0 for that gene. If fewer than
# min_present markers pass, the label is DEMOTED to `demoted_to` (and the
# audit CSV records the outcome). All cells in that (cluster, age) group
# share the demoted label.
# ---------------------------------------------------------------------------

MARKER_PRESENCE_THRESHOLD = 0.20   # >=20% of (cluster x age) cells must express (lognorm > 0)

BRAIN_GATE_CONFIG = {
    "microglia": {
        "keywords":   ["microglia", "microglial"],
        "markers":    ["Cx3cr1", "P2ry12", "Tmem119", "Csf1r", "Aif1"],
        "min_present": 2,
        "demoted_to": "unassigned_immune",
    },
    "astrocyte": {
        "keywords":   ["astrocyte", "astro"],
        "markers":    ["Aqp4", "Gja1", "Slc1a3", "Aldh1l1"],
        "min_present": 2,
        "demoted_to": "unassigned_glia",
    },
    "ol_lineage": {
        # Matches "Oligodendrocyte", "OPC", "Mature OL", "31 OPC-Oligo" etc.
        # Word boundary handled by checking " ol " or label endings; substrings
        # of "oligodendrocyte" / "opc" / " oligo" suffice for ABC + Di Bella vocabs.
        "keywords":   ["oligodendrocyte", "opc", "oligo"],
        "markers":    ["Mbp", "Mog", "Plp1", "Mag"],
        "min_present": 1,
        "demoted_to": "unassigned_glia",
    },
    "endothelial": {
        "keywords":   ["endothelial", "endothelium"],
        "markers":    ["Cldn5", "Pecam1", "Cdh5"],
        "min_present": 2,
        "demoted_to": "unassigned_vascular",
    },
}

# Age-composition expectations. Pure DIAGNOSTIC -- writes a separate CSV of
# (cluster x age) rows whose CellTypist winner is developmentally implausible.
# Does NOT enforce / demote. The point is to give reviewers a short
# "investigate me" list rather than re-reading the full 143-row audit.
BRAIN_AGE_EXPECTATIONS = [
    # (substring_in_winner_lower, expected_ages_set, flag_name)
    ("radial glia",          {"P1"},        "radial_glia_outside_P1"),
    ("intermediate progenitor", {"P1"},     "ipc_outside_P1"),
    ("ipc",                  {"P1"},        "ipc_outside_P1"),
    ("neuroblast",           {"P1"},        "neuroblast_outside_P1"),
    ("glioblast",            {"P1"},        "glioblast_outside_P1"),
    ("erythrocyte",          {"P1"},        "erythrocyte_outside_P1"),
    ("erythroid progenitor", {"P1"},        "erythroid_progenitor_outside_P1"),
]


# ---------------------------------------------------------------------------
# PLACENTA constants — STAMP correlation labelling
# ---------------------------------------------------------------------------

# Tier 1 + 2 collapse map (locked 2026-06-10).
# Key = collapsed label, value = list of STAMP fine labels to map to it.
# Same-lineage sub-states that can't be separated at main-pass Leiden
# resolution are merged. Phase 7b subclustering can split them on a focal
# subset if needed.
STAMP_LABEL_COLLAPSE = {
    "DSC":         ["Angiogenic DSC", "Nourishing DSC", "DSC precursor", "Endometrium stroma"],
    "Endothelium": ["Angiogenic EC", "Venous EC", "Proliferating EC"],
    "SynTI":       ["SynTI", "SynTI precursor"],
    "SynTII":      ["SynTII", "SynTII precursor"],
    "S-TGC":       ["S-TGC", "S-TGC precursor"],
    "SpT":         ["SpT", "SpT precursor"],
    "GC":          ["GC", "GC-1", "GC-2", "GC precursor"],
    "LaTP":        ["LaTP", "LaTP2"],
    "Myeloid":     ["Macrophage", "Monocyte"],
    "Lymphoid":    ["NK", "T cell", "B cell"],
    # The rest stay as their own labels:
    # P-TGC (literature has clean Prl3d1 marker; keep separate)
    # JZP, Lymphatic EC, Mesenchyme, DC, Erythrocyte, Megakaryocyte,
    # Neutrophil, Yolk sac epithelial, Parietal endoderm, Epithelial
}

# Lineage groups for the lineage-aware low-confidence flag.
# If a cell's top and runner-up are in the same lineage group, the low-gap
# flag is suppressed (they're sister states, not different cell types).
STAMP_LINEAGE_GROUPS = {
    "dsc":      ["Angiogenic DSC", "Nourishing DSC", "DSC precursor", "Endometrium stroma"],
    "ec":       ["Angiogenic EC", "Venous EC", "Proliferating EC"],
    "tgc_S":    ["S-TGC", "S-TGC precursor"],
    "synti":    ["SynTI", "SynTI precursor"],
    "syntii":   ["SynTII", "SynTII precursor"],
    "spt":      ["SpT", "SpT precursor"],
    "gc":       ["GC", "GC-1", "GC-2", "GC precursor", "JZP"],
    "latp":     ["LaTP", "LaTP2"],
    "myeloid":  ["Macrophage", "Monocyte", "DC"],
    "lymphoid": ["NK", "T cell", "B cell"],
}

# Strict canonical-marker gates for noisy labels.
# Two-stage check per label:
#   required: ALL of these markers must be present (expr_threshold met).
#             Empty list = no required marker.
#   confirming: at least N of these must also be present (in addition to required).
# Threshold lognorm mean > 0.1 in the cluster counts as 'present' (see _check_markers).
# A cluster failing the gate is demoted to the `demoted_to` label.
STAMP_GATE_CONFIG = {
    "Neutrophil": {
        "required":   ["Ly6g"],                          # neutrophil-specific
        "confirming": ["S100a8", "S100a9", "Csf3r"],     # supporting
        "n_confirming": 1,                               # need at least 1 of confirming
        "demoted_to": "unassigned_immune",
    },
    "Lymphoid": {
        # Top-level Lymphoid gate: must look immune-lymphoid at all.
        "required":   [],
        "confirming": ["Ncr1", "Klrb1c", "Gzma", "Gzmb", "Cd3d", "Cd3e", "Cd19", "Ms4a1"],
        "n_confirming": 2,
        "demoted_to": "unassigned_immune",
    },
    "Megakaryocyte": {
        "required":   [],
        "confirming": ["Pf4", "Itga2b", "Gp1ba", "Gp9", "Ppbp"],
        "n_confirming": 2,
        "demoted_to": "unassigned_blood",
    },
}

# Lymphoid sub-class canonicals — applied AFTER a cluster passes the top-level
# Lymphoid gate. Each sub-class requires ≥2 of its markers. If multiple
# sub-classes pass, the cluster becomes "Lymphoid_mixed".
STAMP_LYMPHOID_SUBSPLIT = {
    "NK":     ["Ncr1", "Klrb1c", "Gzma", "Gzmb", "Eomes"],
    "T cell": ["Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a", "Cd8b1"],
    "B cell": ["Cd19", "Ms4a1", "Cd79a", "Cd79b"],
}
STAMP_LYMPHOID_SUBSPLIT_MIN = 2  # require ≥2 markers per sub-class to assign

# EPC/TSC negative-control markers (should NOT score positive at E12.5/E18.5).
STAMP_EPC_QC_MARKERS = ["Cdx2", "Eomes", "Elf5", "Esrrb"]

# Y-chromosome markers for compartment scoring (fetal-male disambiguation).
STAMP_Y_GENES = ["Ddx3y", "Eif2s3y", "Uty", "Kdm5d"]

# Contamination labels (kept in main analysis but flagged for Phase 8 filtering).
STAMP_CONTAMINATION = {"Yolk sac epithelial", "Parietal endoderm"}


# ---------------------------------------------------------------------------
# Built-in fallback markers (placenta, or brain when no CellTypist).
# These mostly serve the cross-check / provisional-label path.
# ---------------------------------------------------------------------------

BRAIN_MARKERS = {
    "Excitatory neurons":  ["Slc17a7", "Slc17a6", "Neurod2", "Neurod6", "Tbr1"],
    "Inhibitory neurons":  ["Gad1", "Gad2", "Slc32a1", "Lhx6", "Sst", "Pvalb", "Vip"],
    "Astrocytes":          ["Aqp4", "Gfap", "Aldh1l1", "Sox9", "Slc1a3"],
    "Oligodendrocytes":    ["Mbp", "Plp1", "Mog", "Cnp", "Mal"],
    "OPC":                 ["Pdgfra", "Cspg4", "Olig2", "Sox10"],
    "Microglia":           ["Cx3cr1", "P2ry12", "Tmem119", "Hexb"],
    "Endothelial":         ["Cldn5", "Pecam1", "Ly6c1"],
    "Pericytes":           ["Pdgfrb", "Rgs5", "Vtn"],
    "Radial glia / NPCs":  ["Nes", "Sox2", "Pax6", "Vim", "Fabp7"],
    "Choroid plexus":      ["Ttr", "Folr1", "Clic6"],
}

PLACENTA_MARKERS = {
    "Trophoblast (SpT)":     ["Prl3d1", "Prl8a8", "Tpbpa"],
    "Trophoblast (LaT)":     ["Gcm1", "Syna", "Synb"],
    "Trophoblast (TGC)":     ["Prl3b1", "Hand1"],
    "Decidual stromal":      ["Igfbp1", "Foxo1"],
    "Endothelial":           ["Cldn5", "Pecam1"],
    "Hofbauer / Macrophage": ["Cd68", "Adgre1", "Cx3cr1"],
    "NK cells":              ["Ncr1", "Klrb1c", "Gzma"],
    "Erythroblasts":         ["Hbb-bt", "Hba-a1"],
}


def get_markers(cfg: dict) -> dict:
    custom = cfg.get("annotation", {}).get("markers")
    if custom:
        return custom
    return BRAIN_MARKERS if cfg["tissue"] == "brain" else PLACENTA_MARKERS


# ===========================================================================
# CellTypist tiers
# ===========================================================================

def run_celltypist_tier(adata_age, model_path, tier_label: str):
    """Run CellTypist on one (age) subset for one tier (class/subclass/region).

    Returns a DataFrame indexed by adata_age.obs_names with columns
      [predicted, conf]
    or None if model_path is None or the model can't be loaded.

    Uses majority_voting=False — we do our OWN per-(cluster x age) majority on
    the class tier downstream; subclass/region stay raw per-cell.
    """
    if model_path is None:
        return None
    try:
        import celltypist
        from celltypist import models
    except ImportError:
        print(f"    [skip {tier_label}] celltypist not installed.")
        return None

    # If model_path isn't a local file, treat as built-in name and verify.
    is_local = Path(model_path).is_file()
    if not is_local:
        try:
            desc = models.models_description()
            available = set(desc["model"].astype(str))
            cand = model_path if model_path.endswith(".pkl") else model_path + ".pkl"
            if cand not in available:
                print(f"    [skip {tier_label}] '{model_path}' not local and not a known built-in.")
                return None
            model_path = cand
        except Exception as e:
            print(f"    [warn {tier_label}] registry check failed: {e}; trying download anyway.")

    try:
        if Path(model_path).is_file():
            model = models.Model.load(model_path)
            print(f"    [{tier_label}] loaded local: {model_path}")
        else:
            print(f"    [{tier_label}] downloading built-in: {model_path}")
            models.download_models(model=model_path, force_update=False)
            model = models.Model.load(model_path)
    except Exception as e:
        print(f"    [skip {tier_label}] failed to load '{model_path}': {e}")
        return None

    # CellTypist expects log1p-normalized counts in .X. Our adata holds raw in .X
    # and lognorm in layers["lognorm"] (computed at the top of main()).
    tmp = adata_age.copy()
    tmp.X = tmp.layers["lognorm"].copy()

    predictions = celltypist.annotate(tmp, model=model, majority_voting=False)
    pred = predictions.predicted_labels
    # With majority_voting=False, predicted_labels is a Series of strings
    # (one per cell). Be defensive across CellTypist versions:
    if isinstance(pred, pd.DataFrame):
        col = "predicted_labels" if "predicted_labels" in pred.columns else pred.columns[0]
        labels = pred[col].astype(str).values
    else:
        labels = np.asarray(pred).astype(str)
    conf = np.asarray(predictions.probability_matrix.max(axis=1)).astype(float)

    return pd.DataFrame({"predicted": labels, "conf": conf}, index=adata_age.obs_names)


def run_all_celltypist_tiers(adata, per_age_models: dict, table_dir: Path) -> None:
    """Run class/subclass/region tiers per age; mutate adata.obs IN PLACE.

    Hard-fail if a class model is missing for any age in the data — class is
    the canonical label, no silent gaps allowed.

    Writes obs columns:
      celltypist_class_predicted, celltypist_class_conf
      celltypist_subclass,        celltypist_subclass_conf
      celltypist_region,          celltypist_region_conf
    P1 (no subclass/region model) -> sentinel strings + NaN conf.
    """
    ages_in_data = sorted(adata.obs["age"].astype(str).unique())
    print(f"  Ages in data:     {ages_in_data}")
    print(f"  Ages configured:  {sorted(per_age_models)}")

    missing_class = [a for a in ages_in_data
                     if not per_age_models.get(a, {}).get("class")]
    if missing_class:
        sys.exit(
            f"ERROR: no class model configured for age(s) {missing_class}.\n"
            f"  Add annotation.celltypist_models.<age>.class to the YAML.\n"
            f"  Refusing to leave cells without a class label."
        )

    # Init columns up front (sentinels + NaN) so cells in a tier that's skipped
    # for their age have a visible value, not NaN that pandas renders as 'nan'.
    adata.obs["celltypist_class_predicted"] = "unset"   # must be overwritten by every age
    adata.obs["celltypist_class_conf"]      = np.nan
    adata.obs["celltypist_subclass"]        = NO_SUBCLASS
    adata.obs["celltypist_subclass_conf"]   = np.nan
    adata.obs["celltypist_region"]          = NO_REGION
    adata.obs["celltypist_region_conf"]     = np.nan

    for age in ages_in_data:
        print(f"\n  --- age={age} ---")
        models_for_age = per_age_models.get(age, {})
        age_mask = adata.obs["age"].astype(str) == age
        adata_age = adata[age_mask].copy()
        print(f"    {adata_age.n_obs:,} cells")

        for tier, (label_col, conf_col) in (
            ("class",    ("celltypist_class_predicted", "celltypist_class_conf")),
            ("subclass", ("celltypist_subclass",        "celltypist_subclass_conf")),
            ("region",   ("celltypist_region",          "celltypist_region_conf")),
        ):
            mp = models_for_age.get(tier)
            res = run_celltypist_tier(adata_age, mp, tier_label=f"{age}/{tier}")
            if res is None:
                continue
            # pandas.loc with the result's index handles alignment by barcode.
            adata.obs.loc[res.index, label_col] = res["predicted"].values
            adata.obs.loc[res.index, conf_col]  = res["conf"].values
            n_types = res["predicted"].nunique()
            print(f"      -> {len(res):,} cells labeled across {n_types} {tier} categories")

    # Sanity gate — every cell must have a class label.
    n_unset = int((adata.obs["celltypist_class_predicted"] == "unset").sum())
    if n_unset:
        sys.exit(f"ERROR: {n_unset} cells have class_predicted='unset' "
                 f"— a tier silently failed. Re-check logs above.")

    # Per-cell raw predictions CSV.
    cols = ["celltypist_class_predicted", "celltypist_class_conf",
            "celltypist_subclass", "celltypist_subclass_conf",
            "celltypist_region", "celltypist_region_conf"]
    adata.obs[cols].to_csv(table_dir / "07_annotation_celltypist_predictions.csv")
    print(f"\n  Per-cell raw predictions: {table_dir / '07_annotation_celltypist_predictions.csv'}")


# ===========================================================================
# Per-(cluster x age) majority -> canonical celltypist_class
# ===========================================================================

def assign_class_per_cluster_age(adata, audit_csv: Path) -> None:
    """Vote per (Leiden cluster x age), assign winner as celltypist_class to
    every cell in that group, write the full audit CSV with the top-3 breakdown
    and a low_purity flag.
    """
    for required in ("celltypist_class_predicted", "leiden", "age"):
        if required not in adata.obs.columns:
            sys.exit(f"ERROR: '{required}' missing from obs — cannot vote.")

    obs = adata.obs[["leiden", "age", "celltypist_class_predicted"]].copy()
    obs["leiden"] = obs["leiden"].astype(str)
    obs["age"]    = obs["age"].astype(str)

    rows = []
    for (leiden, age), grp in obs.groupby(["leiden", "age"], observed=True):
        counts = grp["celltypist_class_predicted"].astype(str).value_counts()
        n = int(counts.sum())
        top = counts.head(3)
        winner = str(top.index[0])
        win_frac = float(top.iloc[0] / n)
        rup_frac = float(top.iloc[1] / n) if len(top) > 1 else 0.0
        low_purity = (win_frac < LOW_PURITY_THRESHOLD) or (rup_frac > RUNNER_UP_THRESHOLD)
        rows.append({
            "leiden":       leiden,
            "age":          age,
            "n_cells":      n,
            "winner_class": winner,
            "purity":       round(win_frac, 4),
            "top1_label":   winner,
            "top1_frac":    round(win_frac, 4),
            "top2_label":   str(top.index[1]) if len(top) > 1 else "",
            "top2_frac":    round(rup_frac, 4),
            "top3_label":   str(top.index[2]) if len(top) > 2 else "",
            "top3_frac":    round(float(top.iloc[2] / n), 4) if len(top) > 2 else 0.0,
            "low_purity":   bool(low_purity),
        })

    audit = pd.DataFrame(rows).sort_values(["leiden", "age"]).reset_index(drop=True)
    audit.to_csv(audit_csv, index=False)

    # Map winner/purity back per cell via merge on (leiden, age).
    winner_map = audit.set_index(["leiden", "age"])[["winner_class", "purity"]]
    keyed = obs[["leiden", "age"]].merge(
        winner_map, left_on=["leiden", "age"], right_index=True, how="left")
    adata.obs["celltypist_class"]        = keyed["winner_class"].values
    adata.obs["celltypist_class_purity"] = keyed["purity"].values

    n_groups = len(audit)
    n_low = int(audit["low_purity"].sum())
    print(f"\n  Class assignment: {n_groups} (cluster x age) groups; "
          f"{n_low} flagged low_purity")
    if n_low:
        print(f"  Examples to review (lowest purity first):")
        for _, r in audit[audit["low_purity"]].sort_values("purity").head(5).iterrows():
            print(f"    leiden={r['leiden']:>3s} age={r['age']:<4s} n={r['n_cells']:>6d} "
                  f"winner={r['winner_class']!r} purity={r['purity']:.0%} "
                  f"runner-up={r['top2_label']!r} ({r['top2_frac']:.0%})")
    print(f"  Audit CSV: {audit_csv}")


def apply_brain_marker_gate(adata, audit_csv: Path) -> None:
    """STRICT canonical-marker gate for borderline brain CellTypist calls.

    For each (cluster x age) row whose winner_class matches a key in
    BRAIN_GATE_CONFIG, check the cluster's expression of canonical markers.
    A marker is "present" if >=MARKER_PRESENCE_THRESHOLD fraction of the
    cells in that (cluster, age) group have lognorm > 0 for that gene.
    If fewer than min_present markers pass, DEMOTE the label to
    `demoted_to` for all cells in that group.

    Augments the audit CSV (in place) with columns:
        markers_checked  -- comma-separated list of canonical markers checked
        markers_present  -- comma-separated subset that passed the threshold
        gate_outcome     -- one of: no_gate, passed, demoted
        gate_label       -- which gate fired (microglia / astrocyte / ...) or ""
    Re-writes audit_csv at the end.
    Requires lognorm layer present on adata (Phase 5 policy: rebuilt at top
    of main()).
    """
    if "lognorm" not in adata.layers:
        sys.exit("ERROR: lognorm layer missing. apply_brain_marker_gate "
                 "requires lognorm; this should have been rebuilt at the "
                 "top of main(). Aborting rather than checking against raw counts.")

    audit = pd.read_csv(audit_csv)
    # NEW columns initialised so every row has a visible value (no NaN-strings)
    audit["markers_checked"] = ""
    audit["markers_present"] = ""
    audit["gate_outcome"]    = "no_gate"
    audit["gate_label"]      = ""

    # Build a symbol -> var_name lookup. var_names may be Ensembl IDs with
    # symbols in var['symbol']; fall back to identity if no symbol col.
    if "symbol" in adata.var.columns:
        # symbol column might have duplicates; first wins (matches scanpy default).
        symbol_to_var = (
            adata.var.reset_index()
            .drop_duplicates(subset="symbol", keep="first")
            .set_index("symbol")[adata.var.index.name or "index"]
            .to_dict()
        )
    else:
        symbol_to_var = {g: g for g in adata.var_names}

    X_log = adata.layers["lognorm"]
    var_name_to_idx = {v: i for i, v in enumerate(adata.var_names)}

    obs_leiden = adata.obs["leiden"].astype(str).values
    obs_age    = adata.obs["age"].astype(str).values

    demote_count = 0
    pass_count = 0
    skip_genes_total = []

    for i, row in audit.iterrows():
        winner_lc = str(row["winner_class"]).lower()
        gate_label, gate_cfg = None, None
        for name, cfg in BRAIN_GATE_CONFIG.items():
            if any(kw in winner_lc for kw in cfg["keywords"]):
                gate_label, gate_cfg = name, cfg
                break
        if gate_label is None:
            continue

        mask = (obs_leiden == str(row["leiden"])) & (obs_age == str(row["age"]))
        if mask.sum() == 0:
            continue
        cell_idx = np.where(mask)[0]

        markers_checked = gate_cfg["markers"]
        markers_present = []
        for sym in markers_checked:
            var_name = symbol_to_var.get(sym)
            if var_name is None or var_name not in var_name_to_idx:
                skip_genes_total.append(sym)
                continue
            gi = var_name_to_idx[var_name]
            col = X_log[cell_idx, gi]
            if sparse.issparse(col):
                col = col.toarray().ravel()
            else:
                col = np.asarray(col).ravel()
            frac = float((col > 0).mean())
            if frac >= MARKER_PRESENCE_THRESHOLD:
                markers_present.append(sym)

        n_present = len(markers_present)
        passed = n_present >= gate_cfg["min_present"]

        audit.at[i, "markers_checked"] = ",".join(markers_checked)
        audit.at[i, "markers_present"] = ",".join(markers_present)
        audit.at[i, "gate_outcome"]    = "passed" if passed else "demoted"
        audit.at[i, "gate_label"]      = gate_label

        if passed:
            pass_count += 1
        else:
            demote_count += 1
            # Demote on obs. Cast to object to avoid Categorical assignment errors.
            if pd.api.types.is_categorical_dtype(adata.obs["celltypist_class"]):
                adata.obs["celltypist_class"] = adata.obs["celltypist_class"].astype("object")
            adata.obs.loc[mask, "celltypist_class"] = gate_cfg["demoted_to"]

    # Re-write audit
    audit.to_csv(audit_csv, index=False)

    n_gated = pass_count + demote_count
    print(f"\n  Marker gate: {n_gated} (cluster x age) rows matched a gated label "
          f"({pass_count} passed, {demote_count} demoted)")
    if skip_genes_total:
        uniq = sorted(set(skip_genes_total))
        print(f"  WARN: {len(uniq)} canonical marker(s) not found in var_names: {uniq}")
    if demote_count:
        demoted = audit[audit["gate_outcome"] == "demoted"][
            ["leiden", "age", "n_cells", "winner_class", "gate_label",
             "markers_present", "markers_checked"]
        ]
        print(f"  Demoted rows:")
        for _, r in demoted.iterrows():
            present = r["markers_present"] if r["markers_present"] else "(none)"
            print(f"    leiden={str(r['leiden']):>3s} age={r['age']:<4s} "
                  f"n={int(r['n_cells']):>6d} winner={r['winner_class']!r} "
                  f"-> demoted ({r['gate_label']}); markers present: {present}")


def check_brain_age_composition(audit_csv: Path, sanity_csv: Path) -> None:
    """Diagnostic age-composition sanity check (informational only).

    Reads the audit CSV produced by assign_class_per_cluster_age, scans each
    row's winner_class for developmentally-implausible label-at-age combos
    defined in BRAIN_AGE_EXPECTATIONS (e.g. mature erythrocyte at 4W/3mo),
    writes a short CSV. Does NOT modify obs or the audit. Empty file written
    when no flags fire so reviewers can see the check ran.
    """
    audit = pd.read_csv(audit_csv)
    out_rows = []
    for _, r in audit.iterrows():
        winner_lc = str(r["winner_class"]).lower()
        for substr, expected_ages, flag_name in BRAIN_AGE_EXPECTATIONS:
            if substr in winner_lc and str(r["age"]) not in expected_ages:
                out_rows.append({
                    "leiden":         r["leiden"],
                    "age":            r["age"],
                    "n_cells":        int(r["n_cells"]),
                    "winner_class":   r["winner_class"],
                    "flag":           flag_name,
                    "expected_ages":  "/".join(sorted(expected_ages)),
                })
                break  # one flag per row is enough

    cols = ["leiden", "age", "n_cells", "winner_class", "flag", "expected_ages"]
    out = pd.DataFrame(out_rows, columns=cols)
    out.to_csv(sanity_csv, index=False)
    if len(out):
        print(f"\n  Age-composition sanity: {len(out)} developmentally-implausible "
              f"(cluster x age) row(s) -> {sanity_csv}")
        for _, r in out.iterrows():
            print(f"    leiden={str(r['leiden']):>3s} age={r['age']:<4s} "
                  f"n={int(r['n_cells']):>6d} winner={r['winner_class']!r} "
                  f"flag={r['flag']} (expected: {r['expected_ages']})")
    else:
        print(f"\n  Age-composition sanity: no developmentally-implausible "
              f"labels (empty CSV at {sanity_csv})")


# ===========================================================================
# PLACENTA primary track — STAMP reference Spearman correlation
# ===========================================================================

def _load_stamp_reference(path: Path) -> pd.DataFrame:
    """Load STAMP reference matrix (ref_types x genes of avg_log2FC)."""
    with h5py.File(path, "r") as f:
        mat = f["log2fc"][:]
        types = [s.decode() if isinstance(s, bytes) else s for s in f["cell_types"][:]]
        genes = [s.decode() if isinstance(s, bytes) else s for s in f["genes"][:]]
    return pd.DataFrame(mat, index=pd.Index(types, name="cell_type"), columns=genes)


def _rank_along_axis(X: np.ndarray, axis: int = 1) -> np.ndarray:
    return rankdata(X, method="average", axis=axis).astype(np.float32)


def _spearman_matrix(query_ranked: np.ndarray, ref_ranked: np.ndarray) -> np.ndarray:
    """Spearman = Pearson on ranks. Centered + normalized matmul."""
    def standardize(R):
        R = R - R.mean(axis=1, keepdims=True)
        denom = np.linalg.norm(R, axis=1, keepdims=True)
        denom[denom == 0] = 1.0
        return (R / denom).astype(np.float32)
    return standardize(query_ranked) @ standardize(ref_ranked).T


def _stamp_apply_collapse(fine_label: str) -> str:
    for parent, children in STAMP_LABEL_COLLAPSE.items():
        if fine_label in children:
            return parent
    return fine_label


def _stamp_label_to_lineage(label: str):
    for lineage, members in STAMP_LINEAGE_GROUPS.items():
        if label in members:
            return lineage
    return None


def _stamp_canonical_check(adata, cluster_id: str, leiden_key: str,
                           markers: list, sym_to_idx: pd.Series,
                           min_present: int = 1, expr_threshold: float = 0.1):
    """Check if a cluster expresses any of the canonical markers.

    Returns (passed, present_markers). Passed = at least min_present markers
    have mean lognorm > expr_threshold in this cluster.
    """
    mask = adata.obs[leiden_key].astype(str).values == str(cluster_id)
    if mask.sum() == 0:
        return False, []
    X = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
    present = []
    for g in markers:
        if g not in sym_to_idx.index:
            continue
        idx = int(sym_to_idx[g])
        col = X[mask, idx]
        if sparse.issparse(col):
            mean_expr = float(col.mean())
        else:
            mean_expr = float(np.asarray(col).mean())
        if mean_expr > expr_threshold:
            present.append(g)
    return len(present) >= min_present, present


def _stamp_resolve_noisy_label(adata, cluster_id: str, leiden_key: str,
                               majority_collapsed: str, sym_to_idx: pd.Series,
                               expr_threshold: float = 0.1):
    """Strict canonical-marker gate + Lymphoid split.

    Returns (final_label, flags, marker_audit).
      - For labels not in STAMP_GATE_CONFIG: returns (majority_collapsed, [], "")
      - For gated labels that PASS:
          * Lymphoid additionally runs the sub-split: returns NK / T cell /
            B cell / Lymphoid_mixed depending on which sub-canonicals score.
          * Other passing labels return as-is.
      - For gated labels that FAIL: returns the configured demoted_to label
        (e.g. unassigned_immune) and a marker_fail flag.

    The marker_audit string lists which canonicals were present, for the
    sanity-check CSV.
    """
    if majority_collapsed not in STAMP_GATE_CONFIG:
        return majority_collapsed, [], ""

    cfg = STAMP_GATE_CONFIG[majority_collapsed]

    # Required-marker check (e.g., Ly6g for Neutrophil)
    required = cfg.get("required", [])
    required_pass = True
    required_present = []
    if required:
        passed, present = _stamp_canonical_check(
            adata, cluster_id, leiden_key, required, sym_to_idx,
            min_present=len(required), expr_threshold=expr_threshold)
        required_pass = passed
        required_present = present

    # Confirming-marker check
    confirming = cfg.get("confirming", [])
    n_confirming_req = cfg.get("n_confirming", 1)
    confirming_pass, confirming_present = _stamp_canonical_check(
        adata, cluster_id, leiden_key, confirming, sym_to_idx,
        min_present=n_confirming_req, expr_threshold=expr_threshold)

    audit_parts = []
    if required_present:
        audit_parts.append("req=" + ",".join(required_present))
    if confirming_present:
        audit_parts.append("conf=" + ",".join(confirming_present))
    marker_audit = ";".join(audit_parts)

    if not (required_pass and confirming_pass):
        # Gate fail → demote
        return cfg["demoted_to"], [f"marker_gate_fail:{majority_collapsed}"], marker_audit

    # Gate passed. For Lymphoid, run sub-split.
    if majority_collapsed == "Lymphoid":
        sub_hits = {}
        for sub, sub_markers in STAMP_LYMPHOID_SUBSPLIT.items():
            passed, present = _stamp_canonical_check(
                adata, cluster_id, leiden_key, sub_markers, sym_to_idx,
                min_present=STAMP_LYMPHOID_SUBSPLIT_MIN,
                expr_threshold=expr_threshold)
            if passed:
                sub_hits[sub] = present
        if len(sub_hits) == 1:
            sub_label = next(iter(sub_hits))
            audit_parts.append(f"sub={sub_label}({','.join(sub_hits[sub_label])})")
            return sub_label, [], ";".join(audit_parts)
        elif len(sub_hits) > 1:
            audit_parts.append("sub=mixed:" + "+".join(sub_hits.keys()))
            return "Lymphoid_mixed", ["lymphoid_subsplit_mixed"], ";".join(audit_parts)
        else:
            # Top gate passed but no sub-class is clean — keep as Lymphoid
            audit_parts.append("sub=none_clean")
            return "Lymphoid", ["lymphoid_subsplit_unresolved"], ";".join(audit_parts)

    return majority_collapsed, [], marker_audit


def run_stamp_correlation(adata, reference_path: Path, table_dir: Path,
                          leiden_key: str = "leiden",
                          age_key: str = "age",
                          gap_threshold: float = 0.05,
                          purity_threshold: float = 0.5,
                          min_cluster_size: int = 50,
                          chunk: int = 5000) -> None:
    """Placenta primary track: STAMP Spearman correlation + label collapse +
    sanity checks. Writes obs columns and CSV tables. Does NOT plot.

    Obs columns written:
      provisional_celltype_fine   (per-cell STAMP fine label, 35-type)
      provisional_celltype        (per-cell STAMP collapsed label, ~21-type)
      celltype_majority           (per-cluster majority of collapsed; CANONICAL)
      celltype_majority_flag      (per-cluster flag string)
      celltype_gap                (per-cell Spearman gap)
      celltype_top_corr           (per-cell top Spearman correlation)
      compartment                 (fetal_male / xist_positive / ambiguous)

    Tables written (in table_dir):
      07_annotation_cluster_summary.csv
      07_annotation_cluster_fractions.csv
      07_annotation_sanity_check.csv
    """
    print(f"  Loading STAMP reference: {reference_path}")
    if not reference_path.is_file():
        sys.exit(f"ERROR: STAMP reference not found at {reference_path}.\n"
                 f"  Build it with scripts/build_placenta_reference.py first.")
    ref = _load_stamp_reference(reference_path)
    print(f"  Reference: {ref.shape[0]} cell types x {ref.shape[1]:,} genes")

    # Symbol lookup (placenta config has var_names = symbols; no Ensembl mapping)
    if adata.var_names[0].startswith("ENSMUS"):
        if "symbol" not in adata.var.columns:
            sys.exit("ERROR: adata uses Ensembl IDs but no 'symbol' column in .var")
        sym = adata.var["symbol"].astype(str)
        sym_to_idx = pd.Series(np.arange(adata.n_vars), index=sym.values)
    else:
        sym_to_idx = pd.Series(np.arange(adata.n_vars), index=adata.var_names)
    sym_to_idx = sym_to_idx[~sym_to_idx.index.duplicated(keep="first")]

    # Intersect reference genes with query
    common = ref.columns.intersection(sym_to_idx.index)
    print(f"  Common genes (ref ∩ query): {len(common):,}/{ref.shape[1]:,}")
    if len(common) < 500:
        print("  WARN: <500 common genes — annotation will be unreliable",
              file=sys.stderr)

    ref_common = ref[common].values.astype(np.float32)
    query_idx = sym_to_idx.loc[common].values

    # Query expression (lognorm) on common genes — densify (needed for ranking)
    X_lognorm = adata.layers["lognorm"]
    if sparse.issparse(X_lognorm):
        Xs = X_lognorm[:, query_idx].toarray().astype(np.float32)
    else:
        Xs = X_lognorm[:, query_idx].astype(np.float32)
    print(f"  Query matrix on common genes: {Xs.shape}")

    # Chunked Spearman correlation
    print("  Computing Spearman correlations (chunked)...")
    ref_ranked = _rank_along_axis(ref_common, axis=1)
    n_cells = Xs.shape[0]
    corrs = np.empty((n_cells, ref.shape[0]), dtype=np.float32)
    for start in range(0, n_cells, chunk):
        end = min(start + chunk, n_cells)
        q_ranked = _rank_along_axis(Xs[start:end], axis=1)
        corrs[start:end] = _spearman_matrix(q_ranked, ref_ranked)
    del Xs

    # Per-cell top match + gap
    type_names = np.array(list(ref.index), dtype=object)
    top_idx = corrs.argmax(axis=1)
    top_corr = corrs.max(axis=1)
    corrs_gap = corrs.copy()
    corrs_gap[np.arange(n_cells), top_idx] = -np.inf
    runner_corr = corrs_gap.max(axis=1)
    gap = (top_corr - runner_corr).astype(np.float32)
    del corrs, corrs_gap
    top_label_fine = type_names[top_idx]
    top_label_collapsed = np.array(
        [_stamp_apply_collapse(lab) for lab in top_label_fine], dtype=object)

    # Compartment scoring (Xist + Y genes)
    print("  Computing compartment (Xist + Y genes)...")
    compartment = np.full(n_cells, "ambiguous", dtype=object)
    xist_present = "Xist" in sym_to_idx.index
    y_present = [g for g in STAMP_Y_GENES if g in sym_to_idx.index]
    if xist_present and y_present:
        xist_idx = int(sym_to_idx["Xist"])
        x_col = X_lognorm[:, xist_idx]
        xist_score = (np.asarray(x_col.todense()).ravel() if sparse.issparse(x_col)
                      else np.asarray(x_col).ravel())
        y_idx = [int(sym_to_idx[g]) for g in y_present]
        y_cols = X_lognorm[:, y_idx]
        y_score = (np.asarray(y_cols.todense()).mean(axis=1) if sparse.issparse(y_cols)
                   else np.asarray(y_cols).mean(axis=1))
        y_score = np.asarray(y_score).ravel()
        # fetal-male: Y+ Xist-low.  xist_positive: Xist+ Y-low.  rest: ambiguous.
        compartment[(y_score > 0.05) & (xist_score < 0.5)] = "fetal_male"
        compartment[(y_score < 0.01) & (xist_score > 0.5)] = "xist_positive"
        n_male = int((compartment == "fetal_male").sum())
        n_xist = int((compartment == "xist_positive").sum())
        n_amb  = int((compartment == "ambiguous").sum())
        print(f"    fetal_male: {n_male:,}   xist_positive: {n_xist:,}   ambiguous: {n_amb:,}")
    else:
        print(f"    WARN: Xist={xist_present}, Y genes={y_present} — compartment defaulted to ambiguous")

    # Per-cluster aggregation
    print("  Aggregating per cluster...")
    clusters_arr = adata.obs[leiden_key].astype(str).values
    df_cell = pd.DataFrame({
        "cluster": clusters_arr,
        "top_fine": top_label_fine,
        "top_collapsed": top_label_collapsed,
        "top_corr": top_corr,
        "gap": gap,
    })

    summary_rows = []
    fractions_rows = []
    cluster_to_label: dict = {}
    cluster_to_flag: dict = {}

    for cl, sub in df_cell.groupby("cluster"):
        n = len(sub)
        counts_fine = sub["top_fine"].value_counts()
        counts_coll = sub["top_collapsed"].value_counts()
        majority_coll = counts_coll.index[0]
        purity_coll = float(counts_coll.iloc[0] / n)
        majority_fine = counts_fine.index[0]
        purity_fine = float(counts_fine.iloc[0] / n)

        runner_fine = counts_fine.index[1] if len(counts_fine) > 1 else ""
        runner_fine_frac = float(counts_fine.iloc[1] / n) if len(counts_fine) > 1 else 0.0

        mean_gap = float(sub["gap"].mean())
        median_top_corr = float(sub["top_corr"].median())

        # Lineage-aware low-confidence flag
        lin_top = _stamp_label_to_lineage(majority_fine)
        lin_runner = _stamp_label_to_lineage(runner_fine) if runner_fine else None
        same_lineage = lin_top is not None and lin_top == lin_runner

        flags = []
        if n < min_cluster_size:
            flags.append("under_populated")
        elif purity_coll < purity_threshold:
            flags.append("low_purity_mixed")
        elif mean_gap < gap_threshold and not same_lineage:
            flags.append("low_confidence_different_lineage")

        if majority_coll in STAMP_CONTAMINATION:
            flags.append("contamination_putative")

        # Noisy-label resolver: strict gate + Lymphoid split. Applies only to
        # clusters above min_cluster_size — under-populated clusters keep their
        # majority label so we don't accidentally demote tiny but real types.
        marker_audit = ""
        if majority_coll in STAMP_GATE_CONFIG and n >= min_cluster_size:
            new_label, gate_flags, marker_audit = _stamp_resolve_noisy_label(
                adata, cl, leiden_key, majority_coll, sym_to_idx,
                expr_threshold=0.1)
            if new_label != majority_coll:
                majority_coll = new_label
            flags.extend(gate_flags)

        cluster_to_label[cl] = majority_coll
        cluster_to_flag[cl] = ";".join(flags)

        summary_rows.append({
            "cluster": cl, "n_cells": n,
            "majority_label": majority_coll,
            "majority_fine": majority_fine,
            "purity_collapsed": round(purity_coll, 3),
            "purity_fine": round(purity_fine, 3),
            "runner_up_fine": runner_fine,
            "runner_up_fine_frac": round(runner_fine_frac, 3),
            "same_lineage_runner": same_lineage,
            "mean_gap": round(mean_gap, 4),
            "median_top_corr": round(median_top_corr, 3),
            "flags": ";".join(flags),
            "marker_audit": marker_audit,
        })
        for ct, c in counts_fine.head(5).items():
            fractions_rows.append({
                "cluster": cl, "cell_type_fine": ct,
                "fraction": round(float(c) / n, 4),
                "n_cells": int(c),
            })

    summary = (pd.DataFrame(summary_rows)
                 .sort_values("n_cells", ascending=False)
                 .reset_index(drop=True))
    fractions = pd.DataFrame(fractions_rows)
    summary.to_csv(table_dir / "07_annotation_cluster_summary.csv", index=False)
    fractions.to_csv(table_dir / "07_annotation_cluster_fractions.csv", index=False)

    # Sanity check (per-cluster age + EPC)
    print("  Sanity checks (age composition, EPC QC)...")
    sanity_rows = []
    age_col = adata.obs[age_key].astype(str).values if age_key in adata.obs else None
    for cl in summary["cluster"]:
        lab = cluster_to_label[cl]
        mask = clusters_arr == cl
        n = int(mask.sum())
        if age_col is not None:
            age_counts = pd.Series(age_col[mask]).value_counts(normalize=True).to_dict()
            pct_e125 = round(age_counts.get("E12.5", 0.0) * 100, 1)
            pct_e185 = round(age_counts.get("E18.5", 0.0) * 100, 1)
        else:
            pct_e125 = pct_e185 = np.nan
        # EPC QC: should NOT be positive
        epc_present_pass, epc_present = _stamp_canonical_check(
            adata, cl, leiden_key, STAMP_EPC_QC_MARKERS, sym_to_idx,
            min_present=2, expr_threshold=0.1)
        epc_qc_pass = not epc_present_pass  # QC passes if EPC markers NOT expressed
        sanity_rows.append({
            "cluster": cl, "n_cells": n, "assigned_label": lab,
            "pct_E12.5": pct_e125, "pct_E18.5": pct_e185,
            "epc_qc_pass": epc_qc_pass,
            "epc_markers_present": ";".join(epc_present),
        })
    sanity = pd.DataFrame(sanity_rows)
    sanity.to_csv(table_dir / "07_annotation_sanity_check.csv", index=False)

    # Write obs columns
    adata.obs["provisional_celltype_fine"] = pd.Categorical(top_label_fine)
    # NOTE: provisional_celltype is also written by the marker track later;
    # the STAMP collapsed label is more informative for placenta so we use it.
    # The marker-based provisional gets overwritten — kept as a cross-check
    # via the marker-score UMAPs and rank_genes_groups CSV.
    adata.obs["provisional_celltype_stamp_collapsed"] = pd.Categorical(top_label_collapsed)
    adata.obs["celltype_majority"] = pd.Categorical(
        pd.Series(clusters_arr).map(cluster_to_label).values)
    adata.obs["celltype_majority_flag"] = pd.Categorical(
        pd.Series(clusters_arr).map(cluster_to_flag).fillna("").values)
    adata.obs["celltype_gap"] = gap.astype(np.float32)
    adata.obs["celltype_top_corr"] = top_corr.astype(np.float32)
    adata.obs["compartment"] = pd.Categorical(compartment)

    # Stash sanity table on adata.uns for the plotting step
    adata.uns["_stamp_sanity"] = sanity
    adata.uns["_stamp_summary"] = summary

    # Print top-line summary
    n_unique = summary["majority_label"].nunique()
    n_flagged = (summary["flags"] != "").sum()
    print(f"\n  STAMP correlation labels:")
    print(f"    unique majority labels: {n_unique}")
    print(f"    clusters flagged:       {n_flagged}/{len(summary)}")
    n_epc_fail = int((~sanity["epc_qc_pass"]).sum())
    if n_epc_fail:
        print(f"    EPC QC fails:           {n_epc_fail} clusters "
              f"(Cdx2/Eomes/Elf5/Esrrb positive — review)")


# ===========================================================================
# Marker track
# ===========================================================================

def run_marker_genes(adata, obs_key="leiden", n_jobs: int = -1) -> pd.DataFrame:
    """rank_genes_groups (Wilcoxon).

    scanpy >=1.10 supports n_jobs for Wilcoxon and parallelizes per-cluster.
    Older scanpy silently ignores n_jobs (no error). Either way the call works.
    """
    try:
        sc.tl.rank_genes_groups(adata, groupby=obs_key, method="wilcoxon",
                                layer="lognorm", use_raw=False,
                                key_added="rank_genes_groups",
                                n_jobs=n_jobs)
    except TypeError:
        # scanpy < 1.10: no n_jobs kwarg
        sc.tl.rank_genes_groups(adata, groupby=obs_key, method="wilcoxon",
                                layer="lognorm", use_raw=False,
                                key_added="rank_genes_groups")
    result = sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_groups")
    return (result.sort_values("scores", ascending=False)
                  .groupby("group").head(20).reset_index(drop=True))


def score_marker_sets(adata, markers: dict, n_jobs: int = -1) -> None:
    """Score each marker set in parallel via threading.

    score_genes releases the GIL during numpy/scipy ops, so threading
    (rather than multiprocessing) is correct — adata is shared, no copy needed.
    n_jobs=-1 uses all cores; cap explicitly if workstation is shared.
    """
    from joblib import Parallel, delayed

    def _score_one(ct, genes):
        present = [g for g in genes if g in adata.var_names]
        if not present:
            return None
        key = "score_" + ct.replace(" ", "_").replace("/", "_") \
                          .replace("(", "").replace(")", "")
        # copy=False writes in-place to adata.obs[key]
        sc.tl.score_genes(adata, present, score_name=key,
                          layer="lognorm", copy=False)
        return key

    Parallel(n_jobs=n_jobs, backend="threading")(
        delayed(_score_one)(ct, genes) for ct, genes in markers.items()
    )


def assign_provisional_celltype(adata, markers, cluster_key="leiden") -> str:
    """Marker-based per-cluster fallback (used when no celltypist_class)."""
    score_cols, ct_names = [], []
    for ct in markers:
        key = "score_" + ct.replace(" ", "_").replace("/", "_") \
                          .replace("(", "").replace(")", "")
        if key in adata.obs.columns:
            score_cols.append(key); ct_names.append(ct)
    if not score_cols:
        raise ValueError(
            "No marker scores found — cannot assign provisional cell types.\n"
            f"  var_names sample: {list(adata.var_names[:5])}\n"
            "  Marker lists use mouse gene symbols. If var_names are Ensembl IDs,\n"
            "  pass symbol-based marker lists via YAML annotation.markers."
        )
    scores = adata.obs[score_cols].values
    per_cell = np.array([ct_names[i] for i in scores.argmax(axis=1)])
    cluster = adata.obs[cluster_key].astype(str).values
    cluster_label = {}
    for c in np.unique(cluster):
        m = cluster == c
        vals, counts = np.unique(per_cell[m], return_counts=True)
        cluster_label[c] = (vals[counts.argmax()], counts.max() / counts.sum())
    adata.obs["provisional_celltype"] = pd.Categorical(
        [cluster_label[c][0] for c in cluster])
    print(f"  Provisional cell types: {adata.obs['provisional_celltype'].nunique()} types")
    return "provisional_celltype"


# ===========================================================================
# Plots
# ===========================================================================

def _umap_save(adata, color, title, out, figsize=(8, 6), size=6, **kw):
    fig, ax = plt.subplots(figsize=figsize)
    sc.pl.umap(adata, color=color, ax=ax, show=False, frameon=False,
               size=size, title=title, **kw)
    fig.tight_layout(); fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)


def plot_leiden_umap(adata, out):
    n = adata.obs["leiden"].nunique()
    _umap_save(adata, "leiden", f"Leiden clusters ({n})", out,
               figsize=(7, 6), legend_loc="on data", legend_fontsize=7)


def plot_class_umap(adata, out_label, out_conf):
    if "celltypist_class" not in adata.obs.columns:
        return
    n = adata.obs["celltypist_class"].nunique()
    _umap_save(adata, "celltypist_class",
               f"celltypist_class (per cluster x age majority; {n} unique labels)",
               out_label, legend_loc="right margin", legend_fontsize=6)
    if "celltypist_class_conf" in adata.obs.columns:
        _umap_save(adata, "celltypist_class_conf",
                   "celltypist class confidence (per-cell max prob)",
                   out_conf, figsize=(7, 5), color_map="viridis")


def plot_per_age_umap(adata, col, out, title_prefix, sentinel):
    """Faceted UMAP, one panel per age, showing only ages where the model
    actually ran (P1 cells = sentinel -> skip that age). For many-category
    columns (subclass) the legend is suppressed.
    """
    if col not in adata.obs.columns:
        return
    age_col = adata.obs["age"].astype(str)
    ages = []
    for a in sorted(age_col.unique()):
        vals = adata.obs.loc[age_col == a, col].astype(str)
        if not (vals == sentinel).all():
            ages.append(a)
    if not ages:
        return
    n_categories = adata.obs[col].astype(str).nunique()
    legend = None if n_categories > 30 else "right margin"
    ncols = len(ages)
    fig, axes = plt.subplots(1, ncols, figsize=(8 * ncols, 6))
    if ncols == 1:
        axes = [axes]
    for ax, age in zip(axes, ages):
        sub = adata[age_col == age].copy()
        sc.pl.umap(sub, color=col, ax=ax, show=False, frameon=False,
                   size=6, legend_loc=legend, legend_fontsize=5,
                   title=f"{title_prefix} ({age})  n={sub.n_obs:,}")
    fig.tight_layout(); fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)


def plot_marker_dotplot(adata, markers, obs_key, out):
    present = {ct: [g for g in genes if g in adata.var_names]
               for ct, genes in markers.items()}
    present = {ct: g for ct, g in present.items() if g}
    if not present:
        print("  [skip] marker_dotplot: no genes match var_names")
        return
    seen, gene_list = set(), []
    for gs in present.values():
        for g in gs:
            if g not in seen:
                gene_list.append(g); seen.add(g)
    adata.obs[obs_key] = adata.obs[obs_key].astype("category")
    fig = sc.pl.dotplot(adata, var_names=gene_list, groupby=obs_key,
                        layer="lognorm", show=False, return_fig=True,
                        title=f"Curated markers x {obs_key}")
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close()


def plot_top_marker_heatmap(adata, top_markers, obs_key, out, n=10):
    top_genes = (top_markers.groupby("group")
                 .apply(lambda x: x.nlargest(n, "scores"))
                 .reset_index(drop=True)["names"].unique().tolist())
    top_genes = [g for g in top_genes if g in adata.var_names]
    if not top_genes:
        return
    adata.obs[obs_key] = adata.obs[obs_key].astype("category")
    fig = sc.pl.matrixplot(adata, var_names=top_genes, groupby=obs_key,
                           layer="lognorm", show=False, return_fig=True,
                           title=f"Top {n} markers per cluster")
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close()


def plot_marker_score_umaps(adata, markers, out):
    score_keys = []
    for ct in markers:
        key = "score_" + ct.replace(" ", "_").replace("/", "_") \
                          .replace("(", "").replace(")", "")
        if key in adata.obs.columns:
            score_keys.append((ct, key))
    if not score_keys:
        return
    n = len(score_keys); ncols = min(4, n); nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows))
    axes = np.array(axes).flatten()
    for ax, (ct, key) in zip(axes, score_keys):
        sc.pl.umap(adata, color=key, ax=ax, show=False, frameon=False,
                   color_map="Reds", size=6, title=ct)
    for ax in axes[len(score_keys):]:
        ax.set_visible(False)
    fig.tight_layout(); fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)


# ----- placenta-specific plots -----

def plot_celltype_majority_umap(adata, out):
    """Headline UMAP coloured by STAMP-collapsed celltype_majority."""
    if "celltype_majority" not in adata.obs.columns:
        return
    n = adata.obs["celltype_majority"].nunique()
    _umap_save(adata, "celltype_majority",
               f"celltype_majority (STAMP collapsed; {n} unique labels)",
               out, legend_loc="right margin", legend_fontsize=6)


def plot_celltype_gap_umap(adata, out):
    """Per-cell Spearman gap (top - runner-up) = per-cell confidence."""
    if "celltype_gap" not in adata.obs.columns:
        return
    _umap_save(adata, "celltype_gap",
               "Per-cell Spearman gap (top - runner-up) — higher = more confident",
               out, figsize=(7, 5), color_map="viridis")


def plot_compartment_umap(adata, out):
    """Maternal (xist_positive) vs fetal_male vs ambiguous."""
    if "compartment" not in adata.obs.columns:
        return
    _umap_save(adata, "compartment",
               "Compartment (Y+ = fetal_male; Xist+ Y- = xist_positive)",
               out, figsize=(7, 5))


def plot_cluster_fractions_heatmap(adata, out, n_top_types: int = 20):
    """Cluster x cell-type-fraction heatmap (STAMP-collapsed labels)."""
    if "provisional_celltype_stamp_collapsed" not in adata.obs.columns:
        return
    df = pd.DataFrame({
        "cluster": adata.obs["leiden"].astype(str).values,
        "collapsed": adata.obs["provisional_celltype_stamp_collapsed"].astype(str).values,
    })
    top_types = df["collapsed"].value_counts().head(n_top_types).index.tolist()
    mat_rows = []
    for cl, sub in df.groupby("cluster"):
        n = len(sub)
        row = {t: float((sub["collapsed"] == t).sum()) / n for t in top_types}
        row["cluster"] = cl
        row["n_cells"] = n
        mat_rows.append(row)
    mat = pd.DataFrame(mat_rows).set_index("cluster")
    mat = mat.sort_values("n_cells", ascending=False)
    n_cells_col = mat.pop("n_cells")
    fig, ax = plt.subplots(figsize=(max(8, len(top_types) * 0.4),
                                    max(6, len(mat) * 0.2)))
    im = ax.imshow(mat.values, aspect="auto", cmap="magma", vmin=0, vmax=1)
    ax.set_xticks(range(len(top_types)))
    ax.set_xticklabels(top_types, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels([f"{cl} (n={int(n_cells_col[cl])})" for cl in mat.index],
                       fontsize=6)
    ax.set_xlabel("STAMP cell type (collapsed)")
    ax.set_ylabel("Leiden cluster")
    ax.set_title("Cluster x cell-type composition (fraction per cluster)")
    fig.colorbar(im, ax=ax, label="fraction of cluster")
    fig.tight_layout(); fig.savefig(out, dpi=300); plt.close(fig)


def plot_sanity_check_summary(adata, out):
    """Per-cluster age composition + EPC QC red/green marker."""
    sanity = adata.uns.get("_stamp_sanity")
    if sanity is None or not len(sanity):
        return
    sanity_sorted = sanity.sort_values("n_cells", ascending=True)
    fig, ax = plt.subplots(figsize=(10, max(4, len(sanity_sorted) * 0.18)))
    y = list(range(len(sanity_sorted)))
    pass_color = ["#3aa856" if v else "#d33" for v in sanity_sorted["epc_qc_pass"]]
    ax.barh(y, sanity_sorted["pct_E12.5"], color="#4c72b0",
            label="% E12.5", height=0.4)
    ax.barh(y, -sanity_sorted["pct_E18.5"], color="#dd8452",
            label="% E18.5", height=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{cl}: {lab}" for cl, lab in zip(
            sanity_sorted["cluster"], sanity_sorted["assigned_label"])],
        fontsize=6)
    for i, c in enumerate(pass_color):
        ax.scatter([-105], [i], c=c, s=20, marker="s")
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlim(-110, 110)
    ax.set_xlabel("E12.5 ← % cells → E18.5    (red square = EPC marker QC fail)")
    ax.set_title("Per-cluster age composition + EPC QC")
    ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout(); fig.savefig(out, dpi=300); plt.close(fig)


def _stacked_bar(ct_table, title, xlabel, out):
    n = len(ct_table)
    width = max(7.0, 0.4 * n)
    fontsize = 9 if n <= 15 else (8 if n <= 30 else 6)
    cmap = plt.get_cmap("tab20")
    colors = [cmap(i % 20) for i in range(len(ct_table.columns))]
    fig, ax = plt.subplots(figsize=(width, 5))
    ct_table.plot(kind="bar", stacked=True, ax=ax, width=0.8,
                  color=colors, edgecolor="none", legend=True)
    ax.set_ylabel("fraction of cells"); ax.set_xlabel(xlabel); ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=6, ncol=2)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=fontsize)
    fig.tight_layout(); fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)


def plot_cluster_composition_by_sample(adata, out):
    ct = pd.crosstab(adata.obs["leiden"], adata.obs["sample_id"], normalize="index")
    _stacked_bar(ct,
        "Sample composition per Leiden cluster — single-sample bars = potential batch artifact",
        "leiden cluster", out)


def plot_celltype_composition_by_sample(adata, key, out):
    ct = pd.crosstab(adata.obs[key], adata.obs["sample_id"], normalize="index")
    _stacked_bar(ct,
        f"Sample composition per cell type ({key})\n"
        f"Tells you if any cell type is dominated by one sample.",
        "cell type", out)


def plot_celltype_composition_by_group(adata, key, out):
    if "group" not in adata.obs.columns:
        return
    ct = pd.crosstab(adata.obs["sample_id"], adata.obs[key], normalize="index")
    sample_meta = (adata.obs[["sample_id", "group", "age"]]
                   .drop_duplicates().set_index("sample_id").reindex(ct.index))
    ct = ct.loc[sample_meta.sort_values(["group", "age"]).index]
    n = len(ct); width = max(8, 0.4 * n)
    fontsize = 9 if n <= 15 else (8 if n <= 30 else 6)
    cmap = plt.get_cmap("tab20")
    colors = [cmap(i % 20) for i in range(len(ct.columns))]
    fig, ax = plt.subplots(figsize=(width, 5))
    ct.plot(kind="bar", stacked=True, ax=ax, width=0.8, color=colors,
            edgecolor="none", legend=True)
    ax.set_ylabel("fraction of cells")
    ax.set_xlabel("sample (sorted by group x age)")
    ax.set_title(f"Cell type composition per sample [{key}] — quantitative test in 8a")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=6, ncol=2)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=fontsize)
    fig.tight_layout(); fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)


# ===========================================================================
# Per-Leiden-cluster summary CSV (pooled across ages — quick-look)
# ===========================================================================

def build_annotation_summary(adata, top_markers, obs_key="leiden") -> pd.DataFrame:
    clusters = sorted(adata.obs[obs_key].unique(), key=lambda x: int(x))
    has_class = "celltypist_class" in adata.obs.columns
    rows = []
    for c in clusters:
        mask = adata.obs[obs_key] == c
        row = {"cluster": c, "n_cells": int(mask.sum())}
        if has_class:
            row["celltypist_class_majority_pooled"] = \
                adata.obs.loc[mask, "celltypist_class"].astype(str).mode().iloc[0]
            row["celltypist_class_conf_median"] = round(
                float(adata.obs.loc[mask, "celltypist_class_conf"].median()), 3)
        top3 = top_markers[top_markers["group"] == c].head(3)["names"].tolist()
        row["top_markers"] = ", ".join(top3)
        row["manual_annotation"] = ""
        rows.append(row)
    return pd.DataFrame(rows)


def save_composition_table(adata, key, out):
    counts = pd.crosstab(adata.obs["sample_id"], adata.obs[key])
    counts.columns = counts.columns.astype(str)
    fracs = counts.div(counts.sum(axis=1), axis=0)
    meta = (adata.obs[["sample_id", "group", "age", "sex", "pool"]]
            .drop_duplicates().set_index("sample_id"))
    meta.join(fracs).to_csv(out)


def save_composition_by_age(adata, key, out_wide, out_long):
    """Per-age cell-type composition summary.

    Aggregates the per-sample fractions saved by `save_composition_table` up
    to (age × celltype). Writes two views:

      out_wide: rows = celltype, columns = mean/sem per age. Ordered by
                E12.5-or-P1 mean descending (so the most abundant types are
                at the top). One file you can scan visually.

      out_long: rows = (age, celltype). Columns = mean, sem, n_samples,
                min_fraction, max_fraction. Tidy/long form for plotting.

    Both tables compute fractions per sample first, then average ACROSS
    samples within an age — so each pup contributes equally regardless of
    how many cells it has. This is the right summary for biology, not for
    cell-count totals.
    """
    counts = pd.crosstab(adata.obs["sample_id"], adata.obs[key])
    counts.columns = counts.columns.astype(str)
    fracs = counts.div(counts.sum(axis=1), axis=0)

    # Map sample_id -> age (one-to-one)
    sample_age = (adata.obs[["sample_id", "age"]]
                  .drop_duplicates()
                  .set_index("sample_id")["age"]
                  .astype(str))

    df = fracs.join(sample_age).reset_index()
    long_rows = []
    for age, sub in df.groupby("age"):
        n = len(sub)
        for ct in fracs.columns:
            v = sub[ct].astype(float)
            long_rows.append({
                "age": age,
                "celltype": ct,
                "n_samples": int(n),
                "mean_fraction": float(v.mean()),
                "sem_fraction": float(v.sem()) if n > 1 else 0.0,
                "min_fraction": float(v.min()),
                "max_fraction": float(v.max()),
            })
    long = pd.DataFrame(long_rows).sort_values(["age", "mean_fraction"],
                                                ascending=[True, False])
    long.to_csv(out_long, index=False)

    # Wide pivot: rows = celltype, columns = (age, stat)
    wide = long.pivot(index="celltype", columns="age",
                      values=["mean_fraction", "sem_fraction", "n_samples"])
    # Flatten MultiIndex columns so the CSV reads cleanly
    wide.columns = [f"{stat}_{age}" for stat, age in wide.columns]
    # Sort by first age's mean (heuristic: youngest age has the most types active)
    first_age = sorted(set(long["age"]))[0]
    wide = wide.sort_values(f"mean_fraction_{first_age}", ascending=False)
    wide.to_csv(out_wide)


# ===========================================================================
# Main
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(description="Phase 7: annotation (tissue-aware primary track + marker cross-check)")
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()

    print(f"\n=== Phase 7: Annotation ===")
    print(f"Config: {args.config}")
    cfg = load_config(args.config)
    tissue = cfg["tissue"]

    in_path = Path(cfg["results_dir"]) / "h5ad" / "07_clustered" / "all_samples.h5ad"
    if not in_path.is_file():
        sys.exit(f"ERROR: missing {in_path}. Run 06_clustering.py first.")

    out_dir   = Path(cfg["results_dir"]) / "h5ad" / "08_annotated"
    plot_dir  = Path(cfg["results_dir"]) / "plots" / "07_annotation"
    table_dir = phase_table_dir(cfg, "07_annotation")
    for d in (out_dir, plot_dir, table_dir):
        d.mkdir(parents=True, exist_ok=True)

    print(f"\n[1/4] Loading {in_path}...")
    adata = sc.read_h5ad(in_path)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"  Tissue: {tissue}")
    if "leiden" not in adata.obs.columns:
        sys.exit("ERROR: 'leiden' not in obs. Run 06_clustering.py first.")
    print(f"  Leiden clusters: {adata.obs['leiden'].nunique()}")

    # Drop any stale columns from an older 07 run (don't shadow the new schema).
    for stale in ("celltypist_majority", "celltypist_predicted", "celltypist_conf_score"):
        if stale in adata.obs.columns:
            print(f"  Removing stale obs column: {stale}")
            del adata.obs[stale]

    print(f"  Recomputing lognorm layer (dropped after Phase 5)...")
    add_lognorm(adata)

    # --- Tissue-specific primary track ---
    annot_cfg = cfg.get("annotation", {})

    if tissue == "brain":
        print(f"\n[2/4] BRAIN primary track: CellTypist 3-tier annotation...")
        per_age_models = annot_cfg.get("celltypist_models", {}) or {}
        if not per_age_models:
            print(f"  No CellTypist models configured — skipping CellTypist track.")
            print(f"  Provisional marker-based labels will be assigned in step [3/4].")
        else:
            run_all_celltypist_tiers(adata, per_age_models, table_dir)
            audit_path = table_dir / "07_annotation_class_per_cluster_age.csv"
            assign_class_per_cluster_age(adata, audit_csv=audit_path)
            apply_brain_marker_gate(adata, audit_csv=audit_path)
            check_brain_age_composition(
                audit_csv=audit_path,
                sanity_csv=table_dir / "07_annotation_age_composition_sanity.csv",
            )

    elif tissue == "placenta":
        print(f"\n[2/4] PLACENTA primary track: STAMP Spearman correlation...")
        ref_path_str = annot_cfg.get("stamp_reference",
                                     "refs/stamp/stamp_ref_allcells.h5")
        ref_path = Path(ref_path_str)
        if not ref_path.is_absolute():
            # Resolve relative to repo root (script lives in scripts/, repo is parent)
            ref_path = (Path(__file__).parent.parent / ref_path).resolve()
        run_stamp_correlation(
            adata, reference_path=ref_path, table_dir=table_dir,
            leiden_key="leiden", age_key="age",
            gap_threshold=annot_cfg.get("stamp_gap_threshold", 0.05),
            purity_threshold=annot_cfg.get("stamp_purity_threshold", 0.5),
            min_cluster_size=annot_cfg.get("stamp_min_cluster_size", 50),
        )

    else:
        print(f"  WARN: unknown tissue '{tissue}' — skipping primary track.")

    # --- Marker track (cross-check; always runs) ---
    print(f"\n[3/4] Marker track (cross-check)...")
    markers = get_markers(cfg)
    print(f"  Marker sets: {list(markers.keys())}")
    score_marker_sets(adata, markers)
    top_markers = run_marker_genes(adata, "leiden")
    top_markers.to_csv(table_dir / "07_annotation_marker_genes_per_cluster.csv", index=False)
    print(f"  Top markers: {len(top_markers)} rows -> "
          f"{table_dir / '07_annotation_marker_genes_per_cluster.csv'}")

    # Marker-based provisional labels (independent cross-check).
    # For placenta, STAMP run already wrote provisional_celltype_stamp_collapsed
    # — assign_provisional_celltype here writes provisional_celltype using
    # the curated marker dict, which serves as a separate cross-check.
    assign_provisional_celltype(adata, markers)

    summary = build_annotation_summary(adata, top_markers)
    summary.to_csv(table_dir / "07_annotation_summary.csv", index=False)
    print(f"\n  Per-Leiden-cluster summary:")
    print(summary.to_string(index=False))

    # Pick the canonical cell-type key for composition plots / downstream phases.
    #   BRAIN:    celltypist_class (if CellTypist ran) else provisional_celltype
    #   PLACENTA: celltype_majority (STAMP collapsed) else provisional_celltype
    if "celltypist_class" in adata.obs.columns:
        celltype_key = "celltypist_class"
    elif "celltype_majority" in adata.obs.columns:
        celltype_key = "celltype_majority"
    else:
        celltype_key = "provisional_celltype"
    print(f"\n  Canonical celltype key: {celltype_key}")

    # --- Plots ---
    print(f"\n[4/4] Plots...")
    plot_leiden_umap(adata, plot_dir / "umap_leiden_for_annotation.png")

    if tissue == "brain":
        plot_class_umap(adata,
                        plot_dir / "umap_celltypist_class.png",
                        plot_dir / "umap_celltypist_class_confidence.png")
        plot_per_age_umap(adata, "celltypist_subclass",
                          plot_dir / "umap_celltypist_subclass_by_age.png",
                          "celltypist_subclass", NO_SUBCLASS)
        plot_per_age_umap(adata, "celltypist_region",
                          plot_dir / "umap_celltypist_region_by_age.png",
                          "celltypist_region", NO_REGION)
    elif tissue == "placenta":
        plot_celltype_majority_umap(adata,
            plot_dir / "umap_celltype_majority.png")
        plot_celltype_gap_umap(adata,
            plot_dir / "umap_celltype_gap.png")
        plot_compartment_umap(adata,
            plot_dir / "umap_celltype_compartment.png")
        plot_cluster_fractions_heatmap(adata,
            plot_dir / "cluster_fractions_heatmap.png")
        plot_sanity_check_summary(adata,
            plot_dir / "sanity_check_summary.png")

    # Shared marker + composition plots (both tissues)
    plot_marker_dotplot(adata, markers, "leiden", plot_dir / "marker_dotplot.png")
    plot_top_marker_heatmap(adata, top_markers, "leiden",
                             plot_dir / "marker_heatmap_top10.png")
    plot_marker_score_umaps(adata, markers, plot_dir / "umap_marker_scores.png")
    plot_cluster_composition_by_sample(adata,
        plot_dir / "cluster_composition_by_sample.png")
    plot_celltype_composition_by_sample(adata, celltype_key,
        plot_dir / "celltype_composition_by_sample.png")
    plot_celltype_composition_by_group(adata, celltype_key,
        plot_dir / "celltype_composition_by_group.png")
    save_composition_table(adata, celltype_key,
        table_dir / "07_annotation_celltype_composition.csv")
    save_composition_by_age(adata, celltype_key,
        out_wide=table_dir / "07_annotation_celltype_by_age_wide.csv",
        out_long=table_dir / "07_annotation_celltype_by_age_long.csv")

    # Initialize manual_annotation (preserve across reruns).
    if "manual_annotation" not in adata.obs.columns:
        adata.obs["manual_annotation"] = ""

    # Strip uns scratch entries before saving
    for k in ("_stamp_sanity", "_stamp_summary"):
        if k in adata.uns:
            del adata.uns[k]

    # Drop lognorm before saving (Phase 5 policy).
    if "lognorm" in adata.layers:
        del adata.layers["lognorm"]

    adata.write_h5ad(out_dir / "all_samples.h5ad")
    print(f"\n  Written: {out_dir / 'all_samples.h5ad'}")
    print(f"  Plots:   {plot_dir}")
    print(f"  Tables:  {table_dir}")
    print(f"\n✓ Phase 7 complete (canonical celltype key: {celltype_key}).")

    print(f"\nReview before 7b:")
    if tissue == "brain":
        print(f"  - tables/07_annotation/07_annotation_class_per_cluster_age.csv "
              f"(low_purity flags + gate_outcome + markers_present)")
        print(f"  - tables/07_annotation/07_annotation_age_composition_sanity.csv "
              f"(developmentally-implausible label-at-age combos)")
        print(f"  - plots/07_annotation/umap_celltypist_class.png")
        print(f"  - plots/07_annotation/umap_celltypist_subclass_by_age.png "
              f"(4W/3mo facets only)")
    elif tissue == "placenta":
        print(f"  - tables/07_annotation/07_annotation_cluster_summary.csv "
              f"(flags column)")
        print(f"  - tables/07_annotation/07_annotation_sanity_check.csv "
              f"(EPC QC + age composition)")
        print(f"  - plots/07_annotation/umap_celltype_majority.png "
              f"(STAMP-collapsed labels)")
        print(f"  - plots/07_annotation/cluster_fractions_heatmap.png "
              f"(fine-type composition per cluster)")
    print(f"  - plots/07_annotation/cluster_composition_by_sample.png "
          f"(any single-sample clusters?)\n")


if __name__ == "__main__":
    main()
