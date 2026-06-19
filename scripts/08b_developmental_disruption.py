#!/usr/bin/env python
"""
08b_developmental_disruption.py — post-hoc on within_group_across_age.

Tests whether prenatal stress attenuates normal age-dependent gene regulation
("developmental disruption") and/or induces new age-dependent programs
("stress-induced developmental remodeling"). Operates entirely on the Phase 8b
master DE CSV — no re-fitting needed.

For each (tissue, sex stratum, level, region, cell type), partitions genes
into FIVE direction-classes by significance pattern across the three group
columns (Relaxed / Early_Stress / Late_Stress) of the `within_group_across_age`
contrast:

    universal       sig age-DE in ALL three groups (developmental, unaffected by stress)
    relaxed_only    sig in Relaxed, NOT in Early or Late (trajectory LOST under stress)
    stress_shared   sig in BOTH Early AND Late, NOT in Relaxed (trajectory GAINED under stress)
    early_only      sig in Early only
    late_only       sig in Late only

For each class within a slice, computes n_genes, mean and median |log2FC| in
each group (so the user can verify effect-size collapse vs. just-below-threshold
significance shifts).

Writes:
    tables/08b_de/08b_developmental_disruption_summary.csv
    tables/08b_de/08b_developmental_disruption_genes.csv

Sig thresholds (LOCKED, match 08b_de.py):
    padj < 0.05  AND  |log2FC| > 1.0

Usage:
    uv run python scripts/08b_developmental_disruption.py --config config/brain.yaml
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from _utils import load_config, phase_table_dir


PADJ_THR = 0.05
LFC_THR = 1.0

GROUPS = ["Relaxed", "Early_Stress", "Late_Stress"]
NEEDED_COLS = ["contrast", "test_method", "sex", "group_level", "level",
               "celltype", "gene", "log2FC", "padj"]


def _is_sig(sub):
    """Returns boolean Series: padj<0.05 & |log2FC|>1, NaN-safe."""
    return (sub["padj"].notna() & (sub["padj"] < PADJ_THR)
            & sub["log2FC"].notna() & (sub["log2FC"].abs() > LFC_THR))


def classify_slice(rel, early, late):
    """Given three gene-indexed frames (each with log2FC + padj columns) for
    one (sex × level × celltype × region) slice, return dict mapping
    direction-class -> set(gene). Uses index union with NaN sig => False."""
    def sig_set(df):
        s = _is_sig(df)
        return set(df.index[s])

    R, E, L = sig_set(rel), sig_set(early), sig_set(late)
    return {
        "universal":     R & E & L,
        "relaxed_only":  R - E - L,
        "stress_shared": (E & L) - R,
        "early_only":    E - R - L,
        "late_only":     L - R - E,
    }


def per_class_stats(direction, genes, rel, early, late):
    """For a direction-class gene set, compute mean/median |log2FC| in each
    group (using the union of all 3 group frames, NaN-safe)."""
    if not genes:
        return {
            "direction": direction, "n_genes": 0,
            "mean_abs_lfc_Relaxed": np.nan,
            "mean_abs_lfc_Early_Stress": np.nan,
            "mean_abs_lfc_Late_Stress": np.nan,
            "median_abs_lfc_Relaxed": np.nan,
            "median_abs_lfc_Early_Stress": np.nan,
            "median_abs_lfc_Late_Stress": np.nan,
        }
    genes = list(genes)
    out = {"direction": direction, "n_genes": len(genes)}
    for grp_lbl, frame in [("Relaxed", rel), ("Early_Stress", early),
                            ("Late_Stress", late)]:
        present = frame.index.intersection(genes)
        if len(present) == 0:
            out[f"mean_abs_lfc_{grp_lbl}"] = np.nan
            out[f"median_abs_lfc_{grp_lbl}"] = np.nan
            continue
        v = frame.loc[present, "log2FC"].abs()
        out[f"mean_abs_lfc_{grp_lbl}"] = float(v.mean())
        out[f"median_abs_lfc_{grp_lbl}"] = float(v.median())
    return out


def per_gene_long(direction, genes, rel, early, late):
    """Long-form rows: one per (gene × direction) carrying log2FC + padj in
    all three groups for downstream filtering / GO / pathway use."""
    if not genes:
        return []
    rows = []
    for g in genes:
        rec = {"direction": direction, "gene": g}
        for grp_lbl, frame in [("Relaxed", rel), ("Early_Stress", early),
                                ("Late_Stress", late)]:
            if g in frame.index:
                rec[f"log2FC_{grp_lbl}"] = float(frame.loc[g, "log2FC"])
                rec[f"padj_{grp_lbl}"] = float(frame.loc[g, "padj"])
            else:
                rec[f"log2FC_{grp_lbl}"] = np.nan
                rec[f"padj_{grp_lbl}"] = np.nan
        rows.append(rec)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--subcluster", default=None,
                    help="Read the subcluster master CSV instead of the main one.")
    args = ap.parse_args()

    print("\n=== 08b developmental disruption ===")
    cfg = load_config(args.config)
    tissue = cfg.get("tissue")

    table_dir = phase_table_dir(cfg, "08b_de")
    suffix = f"_subcluster_{args.subcluster}" if args.subcluster else ""
    csv_path = table_dir / f"08b_de_results{suffix}.csv"
    if not csv_path.is_file():
        sys.exit(f"ERROR: master CSV not found: {csv_path}")
    print(f"Reading {csv_path} ({csv_path.stat().st_size / 1e6:.1f} MB)...")

    df = pd.read_csv(csv_path, usecols=lambda c: c in NEEDED_COLS,
                     low_memory=False)
    print(f"  {len(df):,} rows loaded.")

    w = df[(df["test_method"] == "Wald")
           & (df["contrast"] == "within_group_across_age")]
    if w.empty:
        sys.exit("No within_group_across_age Wald rows. Nothing to do.")
    print(f"  within_group_across_age Wald rows: {len(w):,}")
    # Each gene appears multiple times per (slice × group_level) due to the
    # pairwise age comparisons (P1/4W, 4W/3mo, P1/3mo). Collapse by keeping
    # the most-significant age-pair as the gene's representative — otherwise
    # frame.loc[g] downstream returns a Series, not a scalar.
    w = (w.sort_values("padj")
          .drop_duplicates(
              ["sex", "level", "celltype", "group_level", "gene"],
              keep="first"))
    print(f"  After collapsing age-pairs (keep smallest padj): {len(w):,} rows")

    summary_rows = []
    gene_rows = []

    # Iterate sex strata × level (whole + regions for brain; whole for placenta)
    # × celltype. Within each (sex, level, celltype), partition by group_level.
    keys = ["sex", "level", "celltype"]
    grouped = w.groupby(keys, observed=True)
    n_slices = 0
    n_slices_skipped = 0
    for (sex_label, level, ct), g in grouped:
        per_grp = {grp: gg.set_index("gene")[["log2FC", "padj"]]
                   for grp, gg in g.groupby("group_level", observed=True)}
        if not all(k in per_grp for k in GROUPS):
            n_slices_skipped += 1
            continue
        rel, early, late = per_grp["Relaxed"], per_grp["Early_Stress"], per_grp["Late_Stress"]
        classes = classify_slice(rel, early, late)

        for direction, genes in classes.items():
            stats = per_class_stats(direction, genes, rel, early, late)
            summary_rows.append({
                "tissue": tissue, "sex": sex_label, "level": level,
                "celltype": ct, **stats,
            })
            for r in per_gene_long(direction, genes, rel, early, late):
                gene_rows.append({
                    "tissue": tissue, "sex": sex_label, "level": level,
                    "celltype": ct, **r,
                })
        n_slices += 1

    print(f"  Slices analyzed: {n_slices:,} "
          f"(skipped {n_slices_skipped} slices missing one of the 3 groups)")
    if n_slices == 0:
        sys.exit("No slice had all three groups present — nothing to write.")

    summary = pd.DataFrame(summary_rows)
    # Order direction column for readability
    dir_order = ["universal", "relaxed_only", "stress_shared",
                 "early_only", "late_only"]
    summary["direction"] = pd.Categorical(summary["direction"],
                                           categories=dir_order, ordered=True)
    summary = summary.sort_values(["sex", "level", "celltype", "direction"]).reset_index(drop=True)

    genes_df = pd.DataFrame(gene_rows)
    genes_df["direction"] = pd.Categorical(genes_df["direction"],
                                            categories=dir_order, ordered=True)
    genes_df = genes_df.sort_values(["sex", "level", "celltype", "direction",
                                      "gene"]).reset_index(drop=True)

    out_sum = table_dir / f"08b_developmental_disruption_summary{suffix}.csv"
    out_genes = table_dir / f"08b_developmental_disruption_genes{suffix}.csv"
    summary.to_csv(out_sum, index=False)
    genes_df.to_csv(out_genes, index=False)
    print(f"\n  Wrote {len(summary):,} summary rows -> {out_sum}")
    print(f"  Wrote {len(genes_df):,} gene rows    -> {out_genes}")

    # Quick view: combined sex × whole level, the headline slice
    print("\nHeadline (sex=combined, level=whole):")
    headline = summary[(summary["sex"] == "combined") & (summary["level"] == "whole")]
    pivot = (headline.set_index(["celltype", "direction"])
                     [["n_genes", "mean_abs_lfc_Relaxed",
                        "mean_abs_lfc_Early_Stress", "mean_abs_lfc_Late_Stress"]]
                     .round(2))
    print(pivot.to_string())

    # Asymmetry ratio (relaxed_only / stress_shared) — the headline number
    print("\nAsymmetry (relaxed_only / stress_shared, sex=combined, level=whole):")
    pivot2 = (headline.pivot_table(index="celltype", columns="direction",
                                    values="n_genes", observed=True, aggfunc="first")
                       .fillna(0))
    if "relaxed_only" in pivot2.columns and "stress_shared" in pivot2.columns:
        pivot2["asym_ratio"] = (pivot2["relaxed_only"]
                                / pivot2["stress_shared"].replace(0, np.nan))
        print(pivot2[["relaxed_only", "stress_shared", "asym_ratio"]]
              .sort_values("asym_ratio", ascending=False).round(2).to_string())


if __name__ == "__main__":
    main()
