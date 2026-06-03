#!/usr/bin/env python
"""
08a_composition.py — Phase 8a: cell type composition analysis (propeller / speckle).

Tests whether stress changes cell type PROPORTIONS, iterating over the
declarative contrasts spec (load_contrasts kind="de").

Statistical unit is the ANIMAL (donor_id): composition is per-donor cell-type
counts (one row per pup), the compositional analog of pseudobulk. Treating cells
as replicates would be wrong (project doc §2). No dam ID -> each pup independent
-> anti-conservative; caveat carried in output via `flag`.

METHOD: propeller (speckle + limma), run via R subprocess (scripts/run_propeller.R),
same clean pattern as scDblFinder. Replaces scCODA, whose TF/arviz/numpy/setuptools
stack proved unworkable. propeller transforms proportions (logit) then runs limma
empirical-Bayes moderated tests — the variance moderation borrows information
across cell types, which suits the small n here. 2-group contrasts use a t-test,
3-group omnibus uses ANOVA F-test; confounders (sex, pool) enter as extra design
columns.

Usage:
  uv run python scripts/08a_composition.py --config config/dev.yaml
  uv run python scripts/08a_composition.py --config config/brain.yaml --min-donors 2
  uv run python scripts/08a_composition.py --config config/brain.yaml --rscript /usr/bin/Rscript

Inputs (first that exists):
  {results_dir}/h5ad/08b_label_transferred/all_samples.h5ad   (Phase 7c)
  {results_dir}/h5ad/08_annotated/all_samples.h5ad            (Phase 7)

Outputs:
  {results_dir}/plots/08a_composition/{contrast}/{level}/
    - stacked_bar.png          : per-donor composition, ordered by group
    - propeller_effects.png    : -log10 FDR per cell type (significant highlighted)
  {results_dir}/tables/composition_results.csv
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from _utils import load_config, load_contrasts, phase_table_dir


LABEL_KEY_PRIORITY = [
    "manual_annotation", "scanvi_celltype", "celltypist_majority", "provisional_celltype",
]


def resolve_celltype_key(adata, explicit):
    if explicit:
        if explicit not in adata.obs.columns:
            sys.exit(f"ERROR: --celltype-key '{explicit}' not in adata.obs. "
                     f"Available: {list(adata.obs.columns)}")
        return explicit
    for key in LABEL_KEY_PRIORITY:
        if key in adata.obs.columns:
            if key == "manual_annotation" and adata.obs[key].astype(str).eq("").all():
                continue
            return key
    sys.exit("ERROR: no usable cell-type label column. Run Phase 7 first or pass "
             "--celltype-key.")


def build_count_matrix(adata, celltype_key, covariates):
    """Per-donor cell-type count matrix + sample-level covariates.

    Returns (df, celltype_cols). df indexed by donor_id, with `covariates`
    columns then one integer-count column per cell type.
    """
    obs = adata.obs
    counts = pd.crosstab(obs["donor_id"], obs[celltype_key])
    counts.columns = [str(c) for c in counts.columns]
    cov = obs[["donor_id"] + covariates].drop_duplicates("donor_id").set_index("donor_id")
    dup = cov.index[cov.index.duplicated()]
    if len(dup):
        sys.exit(f"ERROR: covariate(s) {covariates} vary within donor(s) {list(dup)}.")
    df = cov.join(counts)
    df[counts.columns] = df[counts.columns].fillna(0).astype(int)
    return df, list(counts.columns)


def run_propeller(df, celltype_cols, covariates, test_factor, levels, rscript, transform="logit"):
    """Write composition CSV, call run_propeller.R, read results back.

    levels: [test, ref] for pairwise; None for omnibus ANOVA.
    Returns a results DataFrame. Raises on subprocess failure.
    """
    with tempfile.TemporaryDirectory(prefix="propeller_") as td:
        td = Path(td)
        in_csv, out_csv = td / "counts.csv", td / "res.csv"
        df.to_csv(in_csv)  # donor_id is the index
        cmd = [
            rscript, "scripts/run_propeller.R",
            "--counts", str(in_csv),
            "--celltypes", ",".join(celltype_cols),
            "--covariates", ",".join(covariates),
            "--test", test_factor,
            "--levels", (",".join(levels) if levels else ""),
            "--transform", transform,
            "--out", str(out_csv),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError("propeller R subprocess failed:\n"
                               f"  stdout: {proc.stdout.strip()}\n"
                               f"  stderr: {proc.stderr.strip()[-900:]}")
        return pd.read_csv(out_csv)


def plot_stacked_bar(df, celltype_cols, group_col, out):
    frac = df[celltype_cols].copy()
    frac = frac.div(frac.sum(axis=1), axis=0)
    order = df.sort_values(group_col).index
    frac = frac.loc[order]
    fig, ax = plt.subplots(figsize=(max(7, 0.5 * len(frac)), 5))
    frac.plot(kind="bar", stacked=True, ax=ax, width=0.85, colormap="tab20", edgecolor="none")
    ax.set_ylabel("fraction of cells")
    ax.set_title(f"Per-donor composition (ordered by {group_col})")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=6, ncol=2)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)
    fig.tight_layout(); fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)


def plot_effects(res, out, fdr_thr=0.05):
    """Horizontal bar of -log10(FDR) per cell type; significant highlighted."""
    fdr_col = next((c for c in res.columns if c.upper() == "FDR"), None)
    if fdr_col is None or res.empty:
        return
    r = res.copy().sort_values(fdr_col)
    sig = r[fdr_col] < fdr_thr
    colors = ["salmon" if s else "lightgray" for s in sig]
    fig, ax = plt.subplots(figsize=(7, max(3, 0.4 * len(r))))
    ax.barh(r["celltype"].astype(str), -np.log10(r[fdr_col].clip(lower=1e-300)),
            color=colors, edgecolor="k")
    ax.axvline(-np.log10(fdr_thr), color="k", lw=0.8, ls="--")
    ax.set_xlabel("-log10(FDR)")
    ax.set_title("propeller — differential composition (salmon = FDR<0.05)")
    fig.tight_layout(); fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Phase 8a: composition (propeller)")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--celltype-key", default=None)
    parser.add_argument("--min-donors", type=int, default=None,
                        help="Min donors/group. CLI > YAML > 3. Set 2 to attempt n=2 "
                             "(flagged unreliable).")
    parser.add_argument("--rscript", default=None,
                        help="Path to Rscript (default: found on PATH)")
    args = parser.parse_args()

    print(f"\n=== Phase 8a: Composition analysis (propeller / speckle) ===")
    cfg = load_config(args.config)
    contrasts = load_contrasts(cfg, kind="de")
    min_donors = (args.min_donors if args.min_donors is not None
                  else int(cfg.get("composition", {}).get("min_donors", 3)))
    print(f"  min_donors per group: {min_donors}"
          + ("  (n=2 attempts run FLAGGED unreliable)" if min_donors <= 2 else ""))

    rscript = args.rscript or shutil.which("Rscript")
    if not rscript:
        sys.exit("ERROR: Rscript not found on PATH. Install R + speckle, or pass --rscript.\n"
                 "  R deps:  BiocManager::install(c('speckle','limma'))  + optparse")
    print(f"  Rscript: {rscript}")

    base = Path(cfg["results_dir"]) / "h5ad"
    candidates = [base / "08b_label_transferred" / "all_samples.h5ad",
                  base / "08_annotated" / "all_samples.h5ad"]
    in_path = next((p for p in candidates if p.is_file()), None)
    if in_path is None:
        sys.exit("ERROR: no annotated input found. Checked:\n  " +
                 "\n  ".join(str(p) for p in candidates))
    print(f"  Input: {in_path}")

    plot_root = Path(cfg["results_dir"]) / "plots" / "08a_composition"
    table_dir = phase_table_dir(cfg, "08a_composition")
    plot_root.mkdir(parents=True, exist_ok=True); table_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(in_path)
    celltype_key = resolve_celltype_key(adata, args.celltype_key)
    print(f"  Cell type column: '{celltype_key}' ({adata.obs[celltype_key].nunique()} types)")

    rows = []
    for cname, spec in contrasts.items():
        test = spec.get("test"); flag = spec.get("flag")
        levels = spec.get("levels"); group_by = spec.get("group_by")
        confound = spec.get("confound_warnings", {})

        # propeller handles pairwise (t-test) and 3-group omnibus (ANOVA).
        # Interaction terms aren't a propeller test -> skip, announced.
        if test and ":" in test:
            print(f"\n  [skip] {cname}: interaction not a propeller test.")
            continue
        is_omnibus = (test == "group_omnibus")
        print(f"\n  Contrast: {cname} (test={test}, flag={flag})")

        if isinstance(group_by, list):
            combos = adata.obs[group_by].drop_duplicates().itertuples(index=False)
            group_iter = [dict(zip(group_by, c)) for c in combos]
        elif group_by:
            group_iter = [{group_by: lvl} for lvl in adata.obs[group_by].unique()]
        else:
            group_iter = [{}]

        test_factor = "group" if is_omnibus else test
        design = spec.get("design", f"~ {test_factor}")
        cov_terms = [t.strip() for t in design.replace("~", "").split("+")]
        cov_terms = [t for t in cov_terms if t and t != test_factor and "*" not in t]

        for sel in group_iter:
            sub = adata
            label = "_".join(f"{k}-{v}" for k, v in sel.items()) or "all"
            for k, v in sel.items():
                sub = sub[sub.obs[k] == v]
            if not is_omnibus and levels:
                sub = sub[sub.obs[test].isin(levels)]
            if sub.n_obs == 0:
                continue

            n_donors = sub.obs["donor_id"].nunique()
            if n_donors < max(min_donors, 2) * (1 if is_omnibus else 1):
                # need at least min_donors per group; quick check on total then per-group below
                pass
            # per-group donor counts
            gd = sub.obs.groupby(test_factor, observed=True)["donor_id"].nunique()
            if (gd < min_donors).any() or len(gd) < (3 if is_omnibus else 2):
                print(f"    [skip] {label}: donors/group {dict(gd)} < min_donors={min_donors}.")
                continue
            reliability = "unreliable_n<3" if gd.min() < 3 else "ok"

            covs = [c for c in cov_terms if c in sub.obs.columns and sub.obs[c].nunique() > 1]
            covariates = [test_factor] + covs
            df, ct_cols = build_count_matrix(sub, celltype_key, covariates)
            note = confound.get(sel.get(group_by) if not isinstance(group_by, list) else None, "")

            print(f"    {label}: donors/group {dict(gd)}, {len(ct_cols)} cell types")
            try:
                res = run_propeller(df, ct_cols, covariates, test_factor,
                                    None if is_omnibus else levels, rscript)
            except Exception as e:
                print(f"      [warn] {e}")
                rows.append({"contrast": cname, "flag": flag, "group_level": label,
                             "celltype": None, "statistic": None, "pvalue": None,
                             "fdr": None, "reliability": reliability,
                             "note": f"propeller failed; {note}"})
                continue

            pdir = plot_root / cname / label
            pdir.mkdir(parents=True, exist_ok=True)
            plot_stacked_bar(df, ct_cols, test_factor, pdir / "stacked_bar.png")
            plot_effects(res, pdir / "propeller_effects.png")

            # normalize column names from propeller output
            fdr_col = next((c for c in res.columns if c.upper() == "FDR"), None)
            p_col = next((c for c in res.columns if c.lower() in ("p.value", "pvalue", "p_value")), None)
            stat_col = next((c for c in res.columns if c.lower() in
                             ("tstatistic", "fstatistic", "statistic", "t", "f")), None)
            for _, r in res.iterrows():
                rows.append({
                    "contrast": cname, "flag": flag, "group_level": label,
                    "celltype": r.get("celltype"),
                    "statistic": r.get(stat_col) if stat_col else None,
                    "pvalue": r.get(p_col) if p_col else None,
                    "fdr": r.get(fdr_col) if fdr_col else None,
                    "test_type": r.get("test_type"),
                    "reliability": reliability, "note": note,
                })

    out_csv = table_dir / "08a_composition_results.csv"
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    n_sig = 0
    if rows:
        dfh = pd.DataFrame(rows)
        n_sig = int((dfh["fdr"] < 0.05).sum()) if "fdr" in dfh else 0
    print(f"\n  Master table: {out_csv}  ({len(rows)} rows, {n_sig} at FDR<0.05)")
    print(f"  Plots: {plot_root}")
    print(f"\n✓ Phase 8a complete.")
    print(f"  propeller uses limma moderation (good for small n). Read 'fdr' with")
    print(f"  'reliability' + 'flag' + 'note'. n small — trust only clear shifts.\n")


if __name__ == "__main__":
    main()
