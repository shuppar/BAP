#!/usr/bin/env python
"""
08b_de.py — Phase 8b: pseudobulk differential expression.

THE statistical core. Pseudobulk only — never single-cell-level (cells are not
independent replicates; treating them as such inflates significance and reviewers
will flag it — project doc §2). Statistical unit is the ANIMAL (donor_id).

Plan (locked 2026-06-16, supersedes earlier celltypist_class scheme):

  L1 (sex strata):       combined / M / F  (iter_strata, applied to every contrast)
  L2 (broad cell types): celltypist_broad (brain, ~9 types) / celltype_majority (placenta)
                         × {whole + each celltypist_region}  brain
                         × {whole only}                       placenta
  L3 (subclusters):      subcluster_name (focal coarse types only; loaded via
                         --subcluster <slug>)
                         × {whole + each celltypist_region}  brain
                         × {whole only}                       placenta

`celltypist_class` is NO LONGER used in 8b: it has different vocabularies at P1
(Rosenberg) vs 4W/3mo (ABC), so the same column does not map cleanly across ages.
celltypist_broad is harmonized across ages; region is added as its own axis so
the (cell type × region) grid recovers the regional resolution that class
embedded.

Contrasts handled:
  - per-age Wald            : early_vs_relaxed_per_age, late_vs_relaxed_per_age,
                              early_vs_late_per_age (brain only)
                              -> PyDESeq2, `~ sex + pool + group`
  - cross-age Wald          : within_group_across_age (pairwise age pairs within
                              each group)
                              -> PyDESeq2, `~ sex + age` (pool aliased with age
                              -> dropped + flagged confounded_with_pool)
  - per-age 3-group omnibus : omnibus_3group_per_age   (LRT, df=2)
                              -> R DESeq2 subprocess, full=~sex+pool+group
                              vs reduced=~sex+pool
  - cross-age interaction   : group_x_age_interaction  (LRT, df=4)
                              -> R DESeq2 subprocess, full=~sex+pool+group*age
                              vs reduced=~sex+pool+group+age

Why R for LRT: PyDESeq2 exposes only Wald in its public API; doing a JOINT LRT
(testing several coefficients at once via a chi-squared) requires reaching into
internal arrays whose layout has shifted between 0.4 and 0.5. R DESeq2 has done
this since 2014 in one line: `DESeq(dds, test="LRT", reduced=~...)`. Matches our
SoupX / propeller / scDblFinder R-subprocess pattern.

Phase-8 standard (shared across 8a-8g via _utils, see INSTRUCTIONS):
  - SEX STRATA (iter_strata): combined / M / F on every contrast; every output
    row carries a `sex` column.
  - DROP non-cell-types up front: contaminants (Contamination_* / unresolved in
    subcluster_name) + gate `unassigned*` labels (_utils.unassigned_mask).
    Removed from BOTH numerator and the .X matrix; per-donor dropped counts
    written to a diagnostic CSV.
  - POOL ALIASING: any covariate perfectly aliased with the test factor in a
    slice is dropped from the design and the row flagged confounded_with_pool.
    Canonical case: P1 Late_Stress = Pool3 only -> pool == group.
  - min_donors FLOOR = 2: a group needs >=2 donors to RUN; <reliable_donors (3)
    in any group -> reliability=low_n. Gene-detection floor also = 2 (a gene
    must show >=10 counts in >=2 donors in the slice).
  - PARALLELISM: Wald jobs run in a ProcessPoolExecutor (CPU-bound DESeq2 in
    pure Python); LRT jobs run in a ThreadPoolExecutor (R subprocesses release
    the GIL while running). Both via _utils.parallel_map with --n-jobs.

Significance thresholds (LOCKED, project-wide in 8b):
  padj < 0.05  AND  |log2FC| > 1  (= 2x fold change)
Matches Maitra 2023 / Hwang/Girgenti 2025 / typical psych-genomics convention,
so Phase 9 RRHO and Phase 8f cross-tissue comparisons line up directly.

Caveats carried into output (never silently dropped):
  - No dam ID -> each pup independent. Anti-conservative for litter-aggregated
    traits. Every row carries the contrast `flag`.
  - n=2-4 per group -> only large effects (|log2FC|>1) are trustworthy.
  - Sex strata (M/F) halve n -> low_n. Sex × group INTERACTION not tested at
    this n; per-sex strata describe per-sex effects.
  - LRT uses R DESeq2's dispersion estimator, Wald uses PyDESeq2's. They are
    NOT directly comparable in effect size — never compare a Wald log2FC to an
    LRT p-value on the same gene. Downstream phases (8f / 8g / 9) use Wald
    rows; the LRT rows are tagged `test_method=LRT` and `log2FC=NaN`.

Usage:
  uv run python scripts/08b_de.py --config config/brain.yaml --n-jobs 16
  uv run python scripts/08b_de.py --config config/placenta.yaml --n-jobs 16
  uv run python scripts/08b_de.py --config config/brain.yaml --subcluster immune
  uv run python scripts/08b_de.py --config config/dev_split.yaml --n-jobs 4

Inputs:
  {results_dir}/h5ad/08_annotated/all_samples.h5ad      (main mode)
  {results_dir}/h5ad/08c_subclustered/{slug}.h5ad       (--subcluster mode; also
                                                          consulted in main mode
                                                          for the contam join)

Outputs:
  {results_dir}/tables/08b_de/
    08b_de_results.csv                       [master; columns below]
    08b_dropped_cells_per_donor.csv          [per-donor contam + unassigned mass]
    08b_de_gene_expression_per_sample.csv    [offline audit; sig DEGs only]
    08b_sample_metadata.csv                  [companion to above]
  {results_dir}/plots/08b_de/{sex}/{contrast}/{group_level}/{level}/{ct}/volcano.png

Master CSV columns:
  contrast, flag, test_method (Wald|LRT), sex (combined|M|F),
  group_level (age or group), pair (e.g. [Early_Stress, Relaxed]; NaN for LRT),
  level (whole|<region_name>), celltype, gene,
  log2FC, lfcSE, stat, pvalue, padj, direction (up|down|NaN for LRT),
  n_donors_test, n_donors_ref, n_donors_total,
  reliability (ok|low_n), note
"""

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp

from _utils import (load_config, load_contrasts, phase_table_dir, iter_strata,
                    parallel_map, unassigned_mask)


# ---------------------------------------------------------------------------
# Constants — single source of truth
# ---------------------------------------------------------------------------

CONTAM_PREFIX = "Contamination"

# Significance thresholds, project-wide in 8b. Volcanoes, n_DEGs counts, the
# summary plots in 08b_de_summary.py, and the Venn diagrams all read from here.
PADJ_THR = 0.05
LFC_THR = 1.0   # log2 scale; |log2FC|>1 == >=2x fold change

# Per-tissue tier (L2 main mode) + level grid + focal types for subcluster mode.
# subcluster mode uses subcluster_name as the cell-type column.
TISSUE_CFG = {
    "brain": {
        "tier_l2": "celltypist_broad",         # ~9 broad types, region-free
        "region_key": "celltypist_region",     # whole + each region
        # focal coarse types subclustered in 7b; their value strings must match
        # what's in celltypist_broad. used for the contam join (main mode).
        "focal": ["Immune", "OPC/Oligodendrocytes", "Astrocytes/Ependymal"],
        # columns checked for `unassigned*` gate labels (any match -> drop the cell)
        "ua_keys": ["celltypist_broad", "celltypist_class"],
    },
    "placenta": {
        "tier_l2": "celltype_majority",
        "region_key": None,                    # whole only
        "focal": ["DSC", "Endothelium", "Myeloid", "NK"],
        "ua_keys": ["celltype_majority"],
    },
}


# ---------------------------------------------------------------------------
# Small shared helpers
# ---------------------------------------------------------------------------

def slugify(name: str) -> str:
    return re.sub(r"[^0-9a-zA-Z]+", "_", str(name)).strip("_").lower()


def is_contam(name) -> bool:
    s = str(name)
    return s.startswith(CONTAM_PREFIX) or s == "unresolved"


def read_obs(path, cols=None):
    """Read obs only (backed); restrict to `cols` if given. Used for the contam
    join from 08c subcluster objects without loading their full .X."""
    a = ad.read_h5ad(path, backed="r")
    obs = a.obs if cols is None else a.obs[[c for c in cols if c in a.obs.columns]]
    df = obs.copy()
    a.file.close()
    return df


def aliased_with(df, factor, cov):
    """True if `cov` is perfectly collinear with `factor` in this slice (either
    direction), which makes the design rank-deficient. Classic case: at P1,
    Late_Stress is Pool3-only -> `pool` functionally determined by `group`.
    We drop such covariates and flag the contrast confounded_with_pool, rather
    than crash. Evaluated per donor (covariates are donor-level constants)."""
    dd = df[["donor_id", factor, cov]].drop_duplicates()
    if dd[factor].nunique() < 2 or dd[cov].nunique() < 2:
        return False
    f_determines_c = dd.groupby(factor, observed=True)[cov].nunique().max() <= 1
    c_determines_f = dd.groupby(cov, observed=True)[factor].nunique().max() <= 1
    return bool(f_determines_c or c_determines_f)


def design_terms_from_formula(formula: str, drop_test: str | None = None) -> list:
    """Parse '~ a + b + c' -> ['a','b','c']. Drops `*`-bearing terms (those are
    handled separately by the LRT path) and an optional explicit factor name."""
    terms = [t.strip() for t in formula.replace("~", "").split("+")]
    terms = [t for t in terms if t and "*" not in t]
    if drop_test:
        terms = [t for t in terms if t != drop_test]
    return terms


def informative_error(err_str: str, fallback_chars: int = 160) -> str:
    """Pull the diagnostic line from a multi-line error message.

    R errors print 'ERROR: ...' from our message() call, then quit(status=1)
    appends 'Execution halted' as the final stderr line — which is what the
    naive `splitlines()[-1]` returns and which is useless. Walk back from the
    end looking for ERROR:/Error in/Error: lines first; fall back to the last
    non-empty line if none matches.
    """
    if not err_str or not err_str.strip():
        return err_str
    lines = [ln.strip() for ln in err_str.strip().splitlines() if ln.strip()]
    for ln in reversed(lines):
        low = ln.lower()
        if "error:" in low or low.startswith("error in") or "error in " in low:
            return ln[:fallback_chars]
    return lines[-1][:fallback_chars]


# ---------------------------------------------------------------------------
# Pseudobulk (sum raw counts per donor for one slice)
# ---------------------------------------------------------------------------

def make_pseudobulk(adata_slice, covariates):
    """Sum raw counts per donor across the cells in adata_slice. Returns
    (counts_df [donor x gene, int], metadata_df [donor x covariates],
     n_cells_per_donor Series). Raw counts must be in .X."""
    donors = adata_slice.obs["donor_id"].astype(str)
    uniq = sorted(donors.unique())
    X = adata_slice.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)

    rows = []
    n_cells = {}
    for d in uniq:
        m = (donors == d).values
        n_cells[d] = int(m.sum())
        rows.append(np.asarray(X[m].sum(axis=0)).ravel())
    counts = pd.DataFrame(np.vstack(rows), index=uniq, columns=adata_slice.var_names)
    counts = counts.round().astype(int)   # DESeq2 needs integer counts

    meta = (adata_slice.obs[["donor_id"] + covariates]
            .astype({"donor_id": str})
            .drop_duplicates("donor_id").set_index("donor_id").loc[uniq])
    return counts, meta, pd.Series(n_cells)


def build_expression_matrix(adata, ct_key, genes, samples_meta_cols):
    """Per-sample mean lognorm of `genes` per cell type. Long format. Memory-
    safe: log-normalize only the gene subset using the FULL per-cell library
    size so we never copy the 30K-wide matrix. Mirrors scanpy's
    normalize_total(1e4)+log1p applied to the subset."""
    genes = [g for g in genes if g in adata.var_names]
    if not genes:
        return (pd.DataFrame(),
                adata.obs[["sample_id"] + samples_meta_cols]
                .astype({"sample_id": str}).drop_duplicates("sample_id")
                .set_index("sample_id"))
    gi = [adata.var_names.get_loc(g) for g in genes]

    Xraw = adata.X
    if not sp.issparse(Xraw):
        Xraw = sp.csr_matrix(Xraw)
    tot = np.asarray(Xraw.sum(axis=1)).ravel().astype(float)
    tot[tot == 0] = 1.0
    sub = Xraw[:, gi].astype(float)
    L = sub.multiply((1e4 / tot)[:, None]).tocsr()
    L.data = np.log1p(L.data)

    samples = sorted(adata.obs["sample_id"].astype(str).unique())
    ct_labels = adata.obs[ct_key].astype(str)
    sid = adata.obs["sample_id"].astype(str)

    rows = []
    for ct in sorted(ct_labels.unique()):
        for s in samples:
            m = ((ct_labels == ct) & (sid == s)).values
            n = int(m.sum())
            if n == 0:
                continue
            subL = L[m]
            means = (np.asarray(subL.mean(axis=0)).ravel() if sp.issparse(subL)
                     else np.asarray(subL).mean(axis=0))
            for gene, val in zip(genes, means):
                rows.append({"celltype": ct, "gene": gene, "sample_id": s,
                             "mean_lognorm": round(float(val), 5), "n_cells": n})

    long_df = pd.DataFrame(rows)
    meta = (adata.obs[["sample_id"] + samples_meta_cols]
            .astype({"sample_id": str}).drop_duplicates("sample_id")
            .set_index("sample_id"))
    return long_df, meta


# ---------------------------------------------------------------------------
# PyDESeq2 Wald — fits the full model and returns the Wald results
# ---------------------------------------------------------------------------

def run_pydeseq2_wald(counts, meta, design_terms, contrast_levels, test_factor,
                     n_cpus=1):
    """Fit DESeq2 + Wald test for one pairwise contrast.

    design_terms : list of design factors (intercept implied).
    contrast_levels : [test_level, ref_level] for test_factor.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # Keep only design terms that actually vary (DESeq2 errors on constants).
    design_terms = [t for t in design_terms if meta[t].nunique() > 1]
    if test_factor not in design_terms:
        design_terms.append(test_factor)

    counts = counts.loc[:, counts.sum(axis=0) > 0]
    formula = "~ " + " + ".join(design_terms)

    try:
        from pydeseq2.default_inference import DefaultInference
        inference = DefaultInference(n_cpus=n_cpus)
        dds = DeseqDataSet(counts=counts, metadata=meta, design=formula,
                           refit_cooks=True, inference=inference, quiet=True)
    except (TypeError, ImportError):
        dds = DeseqDataSet(counts=counts, metadata=meta,
                           design_factors=design_terms, refit_cooks=True,
                           n_cpus=n_cpus)
    dds.deseq2()

    contrast = [test_factor, contrast_levels[0], contrast_levels[1]]
    try:
        stat = DeseqStats(dds, contrast=contrast,
                          inference=DefaultInference(n_cpus=n_cpus), quiet=True)
    except (TypeError, NameError):
        stat = DeseqStats(dds, contrast=contrast, n_cpus=n_cpus)
    stat.summary()
    return stat.results_df.copy()


# ---------------------------------------------------------------------------
# R DESeq2 LRT — calls scripts/run_deseq2_lrt.R via subprocess
# ---------------------------------------------------------------------------

def run_r_deseq2_lrt(counts, meta, full_terms, reduced_terms, rscript,
                    test_label="LRT"):
    """Call R DESeq2 LRT for one slice. Writes counts.csv + meta.csv to a temp
    dir, invokes the R worker, reads results.csv back. Returns a results_df
    indexed by gene with columns: baseMean, stat (LRT chi-squared), pvalue, padj.
    log2FC is intentionally NOT returned (LRT tests multiple coefficients
    jointly; any single LFC is ambiguous - downstream callers set log2FC=NaN).
    """
    counts = counts.loc[:, counts.sum(axis=0) > 0]
    full_terms = [t for t in full_terms if meta[t].nunique() > 1]
    reduced_terms = [t for t in reduced_terms if meta[t].nunique() > 1]

    # Sanity: full must add a term beyond reduced; otherwise LRT is undefined.
    extra = set(full_terms) - set(reduced_terms)
    interaction_in_full = any("*" in t or ":" in t for t in full_terms)
    if not extra and not interaction_in_full:
        raise ValueError("LRT requires the full design to add a term beyond the "
                         f"reduced one; full={full_terms} reduced={reduced_terms}")

    with tempfile.TemporaryDirectory(prefix="deseq2_lrt_") as td:
        td = Path(td)
        counts_csv = td / "counts.csv"
        meta_csv = td / "meta.csv"
        out_csv = td / "results.csv"
        # write counts (genes as rows, donors as columns - R convention) and meta
        counts.T.to_csv(counts_csv)         # rows=gene, cols=donor
        meta.to_csv(meta_csv)               # rows=donor, cols=covariate
        # Empty terms list -> '~1' (intercept-only). Triggers when every
        # nuisance covariate becomes constant within the slice (e.g. at 3mo all
        # samples are Pool1; subset to M/F also makes sex constant; reduced
        # becomes ['] -> '~' which R can't parse).
        full_str    = "~" + "+".join(full_terms)    if full_terms    else "~1"
        reduced_str = "~" + "+".join(reduced_terms) if reduced_terms else "~1"
        cmd = [
            rscript, "scripts/run_deseq2_lrt.R",
            "--counts", str(counts_csv),
            "--meta", str(meta_csv),
            "--full", full_str,
            "--reduced", reduced_str,
            "--out", str(out_csv),
            "--label", test_label,
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError("R DESeq2 LRT subprocess failed:\n"
                               f"  stdout: {proc.stdout.strip()[-600:]}\n"
                               f"  stderr: {proc.stderr.strip()[-900:]}")
        res = pd.read_csv(out_csv, index_col=0)
        return res


# ---------------------------------------------------------------------------
# Volcano (Wald only — LRT has no per-gene log2FC)
# ---------------------------------------------------------------------------

def plot_volcano(res, title, out, padj_thr=PADJ_THR, lfc_thr=LFC_THR,
                max_labels=25, symbol_map=None):
    """Volcano with top significant genes labeled (symbols, not Ensembl IDs)."""
    df = res.dropna(subset=["padj", "log2FoldChange"]).copy()
    if df.empty:
        return
    df["gene"] = df.index.astype(str)
    if symbol_map:
        df["gene"] = df["gene"].map(lambda g: symbol_map.get(g, g))
    df["neglog10padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    sig = (df["padj"] < padj_thr) & (df["log2FoldChange"].abs() > lfc_thr)

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.scatter(df.loc[~sig, "log2FoldChange"], df.loc[~sig, "neglog10padj"],
               s=6, color="lightgray", rasterized=True)
    ax.scatter(df.loc[sig, "log2FoldChange"], df.loc[sig, "neglog10padj"],
               s=10, color="salmon", edgecolor="none", rasterized=True)
    ax.axhline(-np.log10(padj_thr), color="k", lw=0.6, ls="--")
    ax.axvline(lfc_thr, color="k", lw=0.6, ls="--")
    ax.axvline(-lfc_thr, color="k", lw=0.6, ls="--")
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 padj")
    ax.set_title(title, fontsize=9)

    to_label = (df[sig]
                .sort_values(["padj", "log2FoldChange"],
                             key=lambda s: s.abs() if s.name == "log2FoldChange" else s,
                             ascending=[True, False])
                .head(max_labels))
    texts = []
    for _, r in to_label.iterrows():
        texts.append(ax.text(r["log2FoldChange"], r["neglog10padj"], r["gene"],
                             fontsize=6, ha="left", va="bottom"))
    if texts:
        try:
            from adjustText import adjust_text
            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle="-", color="gray", lw=0.4))
        except ImportError:
            for t in texts:
                x, y = t.get_position()
                t.set_position((x + 0.05, y + 0.05))

    n_sig = int(sig.sum())
    ax.text(0.02, 0.98, f"{n_sig} sig (padj<{padj_thr}, |LFC|>{lfc_thr})",
            transform=ax.transAxes, fontsize=6, va="top", color="gray")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Parallel workers — top-level for picklability (ProcessPoolExecutor / threads).
# ---------------------------------------------------------------------------

def _wald_worker(job):
    """Run one Wald test (PyDESeq2) + draw its volcano. Top-level so
    ProcessPoolExecutor can pickle it. DESeq2 internal threads pinned to
    job['deseq_cpus'] so N processes don't oversubscribe the cores."""
    res = run_pydeseq2_wald(job["counts"], job["meta"], job["covariates"],
                            job["levels"], job["test_factor"],
                            n_cpus=job["deseq_cpus"])
    plot_volcano(res, job["title"], job["plot_out"], symbol_map=job["symbols"])
    return res


def _lrt_worker(job):
    """Run one LRT (R DESeq2 subprocess). Top-level for ThreadPoolExecutor;
    threads are fine because R is the CPU-bound side and runs in a separate
    process (releases the Python GIL). No volcano — LRT has no per-gene log2FC."""
    res = run_r_deseq2_lrt(job["counts"], job["meta"],
                           job["full_terms"], job["reduced_terms"],
                           job["rscript"], test_label=job["test_label"])
    return res


# ---------------------------------------------------------------------------
# Contaminant join (main mode) — pull subcluster_name from each focal type's
# 08c object so the contaminant mask matches 08a exactly.
# ---------------------------------------------------------------------------

def join_subcluster_name(adata, h5_base, tissue):
    """Return a Series (aligned to adata.obs.index) of subcluster_name pulled
    from the per-focal-type 08c objects, plus the list of focal types found.
    Mirrors 08a so 8a/8b drop the same contaminants."""
    info = TISSUE_CFG[tissue]
    base_tier = info["tier_l2"]
    sub_base = h5_base / "08c_subclustered"
    subname = pd.Series(index=adata.obs.index, dtype="object")
    focal_present = []
    if base_tier not in adata.obs.columns:
        print(f"  [note] base tier '{base_tier}' absent -> no contaminant join.")
        return subname, focal_present
    base_vals = set(adata.obs[base_tier].astype(str))
    for fl in info["focal"]:
        if fl not in base_vals:
            continue
        p = sub_base / f"{slugify(fl)}.h5ad"
        if not p.is_file():
            print(f"  [note] no 08c object for focal '{fl}' ({p.name}); "
                  f"contaminants can't be flagged.")
            continue
        sobs = read_obs(p, cols=["subcluster_name"])
        if "subcluster_name" not in sobs.columns:
            print(f"  [note] {p.name} lacks subcluster_name (run 07d); "
                  f"skipping '{fl}'.")
            continue
        s = sobs["subcluster_name"].reindex(adata.obs.index)
        subname = subname.where(s.isna(), s)
        focal_present.append(fl)
    return subname, focal_present


# ---------------------------------------------------------------------------
# Contrast classification
# ---------------------------------------------------------------------------

def classify_contrast(spec):
    """Map a contrast spec to one of: 'pairwise_age_wald' (per-age pairwise),
    'within_group_wald' (cross-age pairwise within group), 'omnibus_lrt' (per-
    age 3-group F-test as LRT), 'interaction_lrt' (cross-age interaction LRT),
    or 'skip' (anything else — derived contrasts handled by 8g)."""
    test = spec.get("test", "")
    group_by = spec.get("group_by")
    if test == "group_omnibus":
        return "omnibus_lrt"
    if ":" in test:
        return "interaction_lrt"
    if "pairwise" in spec:
        # within_group_across_age: pairwise list, group_by typically 'group'
        return "within_group_wald"
    if spec.get("levels"):
        return "pairwise_age_wald"
    return "skip"


# ---------------------------------------------------------------------------
# Job-building helpers — one per contrast kind. Each appends to (jobs_wald,
# jobs_lrt) and prints a one-line per-slice summary. They share many of the
# same gates so a small per-slice context object is passed in.
# ---------------------------------------------------------------------------

def _slice_donors_ok(m_ct, test, levels, min_donors):
    """Return (n_per_level dict, ok_bool). Each `level` value (e.g. age,
    group) must have >=min_donors donors for the slice to run."""
    g = m_ct.groupby(test, observed=True)["donor_id"].nunique()
    ok = all(int(g.get(lv, 0)) >= min_donors for lv in levels)
    return {lv: int(g.get(lv, 0)) for lv in levels}, ok


def _resolve_design(m_ct, design_terms, test_factor):
    """Build the design term list for one slice, dropping anything that would
    make the model matrix rank-deficient.

    Three layers of cleanup, in order:
      1. Drop covariates that are CONSTANT in this slice (DESeq2 errors on
         single-level factors).
      2. Drop covariates that are perfectly aliased with the TEST FACTOR
         (e.g. P1: Late_Stress is Pool3-only -> pool == group). These get
         flagged confounded_with_pool.
      3. Drop covariates that are perfectly aliased with EACH OTHER (e.g. at
         4W in the combined sex stratum: every M is Pool1, every F is Pool2,
         so sex == pool -> design has two identical columns -> Singular
         matrix). We keep the FIRST one in declared order; the second is
         dropped and flagged.

    Returns (covariates_in_order, dropped_with_test, dropped_inter_aliased)
    where the two dropped lists are reported in the slice `note`.
    """
    # Layer 1: drop constants
    covs = [t for t in design_terms
            if t != test_factor and t in m_ct.columns and m_ct[t].nunique() > 1]
    # Layer 2: drop covariates aliased with the test factor
    aliased_test = [c for c in covs if aliased_with(m_ct, test_factor, c)]
    covs = [c for c in covs if c not in aliased_test]
    # Layer 3: drop covariates aliased with each other (keep first in order)
    keep = []
    dropped_inter = []
    for c in covs:
        if any(aliased_with(m_ct, kept, c) for kept in keep):
            dropped_inter.append(c)
        else:
            keep.append(c)
    return [test_factor] + keep, aliased_test, dropped_inter


def _pseudobulk_or_none(adata, m_ct, covariates, min_samples_expr=2):
    """Build pseudobulk and apply the gene-detection floor. Returns
    (counts, meta) or (None, None) if too few genes survive."""
    ct_ad = adata[adata.obs.index.isin(m_ct.index)]
    counts, meta, _ = make_pseudobulk(ct_ad, covariates)
    keep = (counts >= 10).sum(axis=0) >= min_samples_expr
    counts = counts.loc[:, keep]
    if counts.shape[1] < 10:
        return None, None
    return counts, meta


def queue_pairwise_age_wald(ctx, adata, m_obs, cname, spec, jobs_wald):
    """early_vs_relaxed_per_age / late_vs_relaxed_per_age / early_vs_late_per_age.

    For each (sex × level × age × cell type), test 2 group levels via Wald."""
    test = spec.get("test")                       # 'group'
    levels = list(map(str, spec["levels"]))       # [test_level, ref_level]
    group_by = spec.get("group_by", "age")
    flag = spec.get("flag")
    design_terms = design_terms_from_formula(spec.get("design", f"~ {test}"))
    confound = spec.get("confound_warnings", {}) or {}

    # Optional contrast-level subset (placenta age-specific contrasts).
    c_subset = spec.get("subset", {}) or {}
    for k, v in c_subset.items():
        if k in m_obs.columns:
            m_obs = m_obs[m_obs[k].astype(str) == str(v)]

    m_obs = m_obs[m_obs[test].astype(str).isin(levels)]
    if m_obs.empty:
        return

    for age in sorted(m_obs[group_by].astype(str).unique()):
        m_age = m_obs[m_obs[group_by].astype(str) == age]
        if m_age.empty:
            continue
        slice_note_base = confound.get(age, "")

        for level_label, level_filter in ctx["level_grid"]:
            m_lvl = (m_age if level_filter is None
                     else m_age[m_age[ctx["region_key"]].astype(str) == level_filter])
            if m_lvl.empty:
                continue

            for ct in sorted(m_lvl[ctx["ct_key"]].astype(str).unique()):
                m_ct = m_lvl[m_lvl[ctx["ct_key"]].astype(str) == ct]
                # donors with enough cells of this type
                per_donor = m_ct["donor_id"].value_counts()
                good = per_donor[per_donor >= ctx["min_cells"]].index
                m_ct = m_ct[m_ct["donor_id"].isin(good)]
                if m_ct.empty:
                    continue
                n_per, ok = _slice_donors_ok(m_ct, test, levels, ctx["min_donors"])
                if not ok:
                    continue
                reliability = ("ok" if min(n_per.values()) >= ctx["reliable_donors"]
                               else "low_n")

                covariates, aliased_test, aliased_inter = _resolve_design(
                    m_ct, design_terms, test)
                note = slice_note_base
                if ctx["sex_val"] is not None:
                    note = (note + "; " if note else "") + "sex-specific (low power)"
                if aliased_test:
                    amsg = (f"{'/'.join(aliased_test)} aliased with {test} -> dropped "
                            f"from design (confounded_with_pool)")
                    note = (note + "; " if note else "") + amsg
                if aliased_inter:
                    # e.g. 4W combined: every M->Pool1, every F->Pool2 so sex==pool.
                    # We kept the first in declared order; the second is dropped.
                    imsg = (f"{'/'.join(aliased_inter)} aliased with kept covariate(s)"
                            f" -> dropped from design")
                    note = (note + "; " if note else "") + imsg

                counts, meta = _pseudobulk_or_none(adata, m_ct, covariates,
                                                  min_samples_expr=ctx["min_samples_expr"])
                if counts is None:
                    continue

                ct_slug = slugify(ct)
                plot_out = (ctx["plot_root"] / ctx["sex_label"] / cname / age
                            / level_label / ct_slug / "volcano.png")
                job_syms = ({g: ctx["symbol_map"][g] for g in counts.columns
                             if g in ctx["symbol_map"]} if ctx["symbol_map"] else None)

                jobs_wald.append(dict(
                    contrast=cname, flag=flag, sex=ctx["sex_label"],
                    group_level=age, pair=levels, level=level_label,
                    celltype=ct, reliability=reliability, note=note,
                    n_test=n_per[levels[0]], n_ref=n_per[levels[1]],
                    n_total=int(sum(n_per.values())),
                    counts=counts, meta=meta, covariates=covariates,
                    levels=levels, test_factor=test,
                    symbols=job_syms, deseq_cpus=ctx["deseq_cpus"],
                    title=f"{cname} | sex={ctx['sex_label']}\n"
                          f"{age} | {level_label} | {ct}",
                    plot_out=plot_out))
                print(f"    {ctx['sex_label']}|{age}|{level_label}|{ct}: "
                      f"{n_per[levels[0]]} vs {n_per[levels[1]]} donors, "
                      f"{counts.shape[1]} genes [{reliability}]")


def queue_within_group_wald(ctx, adata, m_obs, cname, spec, jobs_wald):
    """within_group_across_age: for each group (group_by='group'), test pairwise
    age pairs. Design = '~ sex + age' (pool fully aliased with age -> dropped).
    Tier is celltypist_broad (harmonized across ages) per the L2 plan."""
    flag = spec.get("flag")
    pairwise = spec.get("pairwise", [])
    test = spec.get("test")               # 'age'
    group_by = spec.get("group_by", "group")
    design_terms = design_terms_from_formula(spec.get("design", f"~ {test}"))

    for grp in sorted(m_obs[group_by].astype(str).unique()):
        m_grp = m_obs[m_obs[group_by].astype(str) == grp]
        if m_grp.empty:
            continue
        for pair in pairwise:
            pair = list(map(str, pair))
            m_pair = m_grp[m_grp[test].astype(str).isin(pair)]
            if m_pair.empty:
                continue

            for level_label, level_filter in ctx["level_grid"]:
                m_lvl = (m_pair if level_filter is None
                         else m_pair[m_pair[ctx["region_key"]].astype(str) == level_filter])
                if m_lvl.empty:
                    continue

                for ct in sorted(m_lvl[ctx["ct_key"]].astype(str).unique()):
                    m_ct = m_lvl[m_lvl[ctx["ct_key"]].astype(str) == ct]
                    per_donor = m_ct["donor_id"].value_counts()
                    good = per_donor[per_donor >= ctx["min_cells"]].index
                    m_ct = m_ct[m_ct["donor_id"].isin(good)]
                    if m_ct.empty:
                        continue
                    n_per, ok = _slice_donors_ok(m_ct, test, pair, ctx["min_donors"])
                    if not ok:
                        continue
                    reliability = ("ok" if min(n_per.values()) >= ctx["reliable_donors"]
                                   else "low_n")

                    covariates, aliased_test, aliased_inter = _resolve_design(
                        m_ct, design_terms, test)
                    # within_group_across_age is flagged confounded_with_pool by
                    # design (pool determined by age) — record explicitly.
                    note = "pool-age confounded (cross-age within group)"
                    if ctx["sex_val"] is not None:
                        note += "; sex-specific (low power)"
                    if aliased_test:
                        amsg = (f"{'/'.join(aliased_test)} aliased with {test} -> dropped")
                        note += "; " + amsg
                    if aliased_inter:
                        imsg = (f"{'/'.join(aliased_inter)} aliased with kept covariate(s)"
                                f" -> dropped from design")
                        note += "; " + imsg

                    counts, meta = _pseudobulk_or_none(adata, m_ct, covariates,
                                                      min_samples_expr=ctx["min_samples_expr"])
                    if counts is None:
                        continue

                    ct_slug = slugify(ct)
                    pair_slug = f"{pair[0]}_vs_{pair[1]}"
                    grp_slug = slugify(grp)
                    plot_out = (ctx["plot_root"] / ctx["sex_label"] / cname
                                / grp_slug / pair_slug / level_label / ct_slug
                                / "volcano.png")
                    job_syms = ({g: ctx["symbol_map"][g] for g in counts.columns
                                 if g in ctx["symbol_map"]} if ctx["symbol_map"] else None)

                    jobs_wald.append(dict(
                        contrast=cname, flag=flag, sex=ctx["sex_label"],
                        group_level=grp, pair=pair, level=level_label,
                        celltype=ct, reliability=reliability, note=note,
                        n_test=n_per[pair[0]], n_ref=n_per[pair[1]],
                        n_total=int(sum(n_per.values())),
                        counts=counts, meta=meta, covariates=covariates,
                        levels=pair, test_factor=test,
                        symbols=job_syms, deseq_cpus=ctx["deseq_cpus"],
                        title=f"{cname} | sex={ctx['sex_label']} | {grp}\n"
                              f"{pair[0]} vs {pair[1]} | {level_label} | {ct}",
                        plot_out=plot_out))
                    print(f"    {ctx['sex_label']}|{grp}|{pair[0]}v{pair[1]}|"
                          f"{level_label}|{ct}: {n_per[pair[0]]} vs "
                          f"{n_per[pair[1]]} donors, {counts.shape[1]} genes "
                          f"[{reliability}]")


def queue_omnibus_lrt(ctx, adata, m_obs, cname, spec, jobs_lrt):
    """omnibus_3group_per_age: per age, LRT of full=~sex+pool+group vs
    reduced=~sex+pool. df=2 (Early, Late vs Relaxed). One LRT per gene."""
    flag = spec.get("flag")
    group_by = spec.get("group_by", "age")
    full_terms = design_terms_from_formula(spec.get("design", "~ sex + pool + group"))
    reduced_terms = [t for t in full_terms if t != "group"]
    confound = spec.get("confound_warnings", {}) or {}

    for age in sorted(m_obs[group_by].astype(str).unique()):
        m_age = m_obs[m_obs[group_by].astype(str) == age]
        if m_age["group"].nunique() < 2:
            continue
        slice_note_base = confound.get(age, "")

        for level_label, level_filter in ctx["level_grid"]:
            m_lvl = (m_age if level_filter is None
                     else m_age[m_age[ctx["region_key"]].astype(str) == level_filter])
            if m_lvl.empty:
                continue

            for ct in sorted(m_lvl[ctx["ct_key"]].astype(str).unique()):
                m_ct = m_lvl[m_lvl[ctx["ct_key"]].astype(str) == ct]
                per_donor = m_ct["donor_id"].value_counts()
                good = per_donor[per_donor >= ctx["min_cells"]].index
                m_ct = m_ct[m_ct["donor_id"].isin(good)]
                if m_ct.empty or m_ct["group"].nunique() < 2:
                    continue
                # require each present group to have >=min_donors donors
                gd = m_ct.groupby("group", observed=True)["donor_id"].nunique()
                if (gd < ctx["min_donors"]).any() or len(gd) < 2:
                    continue
                reliability = ("ok" if gd.min() >= ctx["reliable_donors"]
                               else "low_n")

                # drop covariates aliased with group; keep full-vs-reduced diff
                full_use = [t for t in full_terms
                            if t in m_ct.columns and m_ct[t].nunique() > 1]
                aliased_test = [c for c in full_use
                                if c != "group" and aliased_with(m_ct, "group", c)]
                full_use = [t for t in full_use if t not in aliased_test]
                # NEW: also drop covariates aliased with each other (e.g. 4W
                # combined: sex==pool -> rank-deficient design). Keep first in
                # declared order.
                kept_in_full = []
                aliased_inter = []
                for t in full_use:
                    if t == "group":
                        kept_in_full.append(t)
                        continue
                    if any(aliased_with(m_ct, kept, t) for kept in kept_in_full
                           if kept != "group"):
                        aliased_inter.append(t)
                    else:
                        kept_in_full.append(t)
                full_use = kept_in_full
                reduced_use = [t for t in full_use if t != "group"]
                if "group" not in full_use:
                    continue          # group itself became constant
                # all covariates kept by m_ct should be in pseudobulk meta too
                covariates = full_use[:]   # for the pseudobulk step

                note = slice_note_base
                if ctx["sex_val"] is not None:
                    note = (note + "; " if note else "") + "sex-specific (low power)"
                if aliased_test:
                    amsg = (f"{'/'.join(aliased_test)} aliased with group -> dropped "
                            f"(confounded_with_pool)")
                    note = (note + "; " if note else "") + amsg
                if aliased_inter:
                    imsg = (f"{'/'.join(aliased_inter)} aliased with kept covariate(s)"
                            f" -> dropped from design")
                    note = (note + "; " if note else "") + imsg

                counts, meta = _pseudobulk_or_none(adata, m_ct, covariates,
                                                  min_samples_expr=ctx["min_samples_expr"])
                if counts is None:
                    continue

                jobs_lrt.append(dict(
                    contrast=cname, flag=flag, sex=ctx["sex_label"],
                    group_level=age, pair=None, level=level_label,
                    celltype=ct, reliability=reliability, note=note,
                    n_total=int(m_ct["donor_id"].nunique()),
                    n_per_group=dict(gd),
                    counts=counts, meta=meta,
                    full_terms=full_use, reduced_terms=reduced_use,
                    test_label="omnibus_3group",
                    rscript=ctx["rscript"]))
                print(f"    {ctx['sex_label']}|{age}|{level_label}|{ct} [LRT/omni]: "
                      f"groups {dict(gd)}, {counts.shape[1]} genes [{reliability}]")


def queue_interaction_lrt(ctx, adata, m_obs, cname, spec, jobs_lrt):
    """group_x_age_interaction: LRT of full=~sex+pool+group*age vs
    reduced=~sex+pool+group+age. df=4 (2 non-ref groups × 2 non-ref ages).
    One LRT per gene per (sex × level × cell type)."""
    flag = spec.get("flag")
    # build full / reduced term lists from formula manually since 'group*age'
    # produces non-atomic terms
    formula = spec.get("design", "~ sex + pool + group * age")
    base_terms = [t.strip() for t in formula.replace("~", "").split("+")]
    base_terms = [t for t in base_terms if t]
    full_terms = base_terms[:]                              # includes 'group*age'
    reduced_terms = [t for t in base_terms if "*" not in t]
    reduced_terms.extend(["group", "age"])                  # main effects only
    reduced_terms = list(dict.fromkeys(reduced_terms))      # dedup

    # interaction is meaningful only with >=2 groups AND >=2 ages
    if m_obs["group"].nunique() < 2 or m_obs["age"].nunique() < 2:
        return

    for level_label, level_filter in ctx["level_grid"]:
        m_lvl = (m_obs if level_filter is None
                 else m_obs[m_obs[ctx["region_key"]].astype(str) == level_filter])
        if m_lvl.empty:
            continue

        for ct in sorted(m_lvl[ctx["ct_key"]].astype(str).unique()):
            m_ct = m_lvl[m_lvl[ctx["ct_key"]].astype(str) == ct]
            per_donor = m_ct["donor_id"].value_counts()
            good = per_donor[per_donor >= ctx["min_cells"]].index
            m_ct = m_ct[m_ct["donor_id"].isin(good)]
            if m_ct.empty:
                continue
            if m_ct["group"].nunique() < 2 or m_ct["age"].nunique() < 2:
                continue
            # each (group, age) cell should have >=1 donor (relax min_donors here
            # since interaction tests across cells of the table)
            cell_table = m_ct.groupby(["group", "age"], observed=True)["donor_id"].nunique()
            if cell_table.empty or cell_table.min() < 1:
                continue
            n_total = int(m_ct["donor_id"].nunique())
            # tag low_n if any cell has <reliable_donors
            reliability = ("ok" if cell_table.min() >= ctx["reliable_donors"]
                           else "low_n")

            # filter constants & resolve aliasing
            atom_full = [t for t in base_terms if "*" not in t]
            atom_full_use = [t for t in atom_full
                             if t in m_ct.columns and m_ct[t].nunique() > 1]
            # Inter-covariate aliasing (rare here since spanning multiple
            # ages usually breaks the M-Pool1 / F-Pool2 alias, but cheap to
            # check). Keep first in declared order; protect group and age.
            kept = []
            inter_aliased = []
            for t in atom_full_use:
                if t in ("group", "age"):
                    kept.append(t)
                    continue
                if any(aliased_with(m_ct, k, t) for k in kept
                       if k not in ("group", "age")):
                    inter_aliased.append(t)
                else:
                    kept.append(t)
            atom_full_use = kept
            covariates = atom_full_use[:]
            # construct full / reduced from what survived
            inter_term = "group*age"
            full_use = atom_full_use + [inter_term]
            reduced_use = atom_full_use[:]
            if "group" not in atom_full_use or "age" not in atom_full_use:
                continue   # need both to test interaction

            note = "interaction (underpowered_exploratory)"
            if ctx["sex_val"] is not None:
                note += "; sex-specific (low power)"
            if inter_aliased:
                note += (f"; {'/'.join(inter_aliased)} aliased with kept "
                         f"covariate(s) -> dropped from design")

            counts, meta = _pseudobulk_or_none(adata, m_ct, covariates,
                                              min_samples_expr=ctx["min_samples_expr"])
            if counts is None:
                continue

            jobs_lrt.append(dict(
                contrast=cname, flag=flag, sex=ctx["sex_label"],
                group_level="cross_age", pair=None, level=level_label,
                celltype=ct, reliability=reliability, note=note,
                n_total=n_total, n_per_group=dict(cell_table),
                counts=counts, meta=meta,
                full_terms=full_use, reduced_terms=reduced_use,
                test_label="interaction",
                rscript=ctx["rscript"]))
            print(f"    {ctx['sex_label']}|cross_age|{level_label}|{ct} "
                  f"[LRT/interaction]: cells {dict(cell_table)}, "
                  f"{counts.shape[1]} genes [{reliability}]")


# ---------------------------------------------------------------------------
# Result collection (Wald + LRT -> unified row schema)
# ---------------------------------------------------------------------------

MASTER_COLS = [
    "contrast", "flag", "test_method", "sex", "group_level", "pair", "level",
    "celltype", "gene",
    "log2FC", "lfcSE", "stat", "pvalue", "padj", "direction",
    "n_donors_test", "n_donors_ref", "n_donors_total",
    "reliability", "note",
]


def _wald_rows(job, res, err):
    """Turn a Wald (job, results_df, err) into MASTER_COLS rows."""
    base = dict(contrast=job["contrast"], flag=job["flag"], test_method="Wald",
                sex=job["sex"], group_level=job["group_level"],
                pair=str(job["pair"]), level=job["level"],
                celltype=job["celltype"],
                n_donors_test=job["n_test"], n_donors_ref=job["n_ref"],
                n_donors_total=job["n_total"],
                reliability=job["reliability"], note=job["note"])
    if err:
        last = informative_error(err)
        r = dict(base, gene=None, log2FC=None, lfcSE=None, stat=None,
                 pvalue=None, padj=None, direction=None)
        r["note"] = (r["note"] + "; " if r["note"] else "") + f"DESeq2 failed: {last[:120]}"
        return [r]
    out = []
    res = res.reset_index().rename(columns={"index": "gene"})
    for _, rr in res.iterrows():
        lfc = rr.get("log2FoldChange")
        out.append(dict(
            base, gene=rr["gene"], log2FC=lfc, lfcSE=rr.get("lfcSE"),
            stat=rr.get("stat"), pvalue=rr.get("pvalue"), padj=rr.get("padj"),
            direction=(None if pd.isna(lfc) else ("up" if lfc > 0 else "down"))))
    return out


def _lrt_rows(job, res, err):
    """Turn an LRT (job, results_df, err) into MASTER_COLS rows. log2FC + lfcSE
    + direction are intentionally NaN — LRT tests multiple coefficients jointly,
    a single LFC would be misleading."""
    base = dict(contrast=job["contrast"], flag=job["flag"], test_method="LRT",
                sex=job["sex"], group_level=job["group_level"],
                pair=None, level=job["level"], celltype=job["celltype"],
                n_donors_test=None, n_donors_ref=None,
                n_donors_total=job["n_total"],
                reliability=job["reliability"], note=job["note"])
    if err:
        last = informative_error(err)
        r = dict(base, gene=None, log2FC=None, lfcSE=None, stat=None,
                 pvalue=None, padj=None, direction=None)
        r["note"] = (r["note"] + "; " if r["note"] else "") + f"LRT failed: {last[:120]}"
        return [r]
    out = []
    res = res.reset_index().rename(columns={res.index.name or "index": "gene"})
    for _, rr in res.iterrows():
        out.append(dict(
            base, gene=rr["gene"], log2FC=np.nan, lfcSE=np.nan,
            stat=rr.get("stat"), pvalue=rr.get("pvalue"), padj=rr.get("padj"),
            direction=None))
    return out


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Phase 8b: pseudobulk DE (Wald + LRT)")
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--celltype-key", default=None,
                    help="Override the per-tissue L2 column (default: "
                         "celltypist_broad / celltype_majority). Use this only "
                         "to force a different label column on the master h5ad.")
    ap.add_argument("--min-cells", type=int, default=10,
                    help="Min cells of a type per donor for that donor to count (default 10).")
    ap.add_argument("--min-donors", type=int, default=None,
                    help="Min donors/group to RUN (floored at 2; CLI > "
                         "composition.min_donors > 2).")
    ap.add_argument("--reliable-donors", type=int, default=None,
                    help="Donors/group at/above which a slice is 'ok' else 'low_n' "
                         "(CLI > composition.reliable_donors > 3).")
    ap.add_argument("--n-jobs", type=int, default=8,
                    help="Concurrent workers (default 8; ~16-24 on the WS). "
                         "Wald uses processes, LRT uses threads (R subprocess "
                         "releases the GIL).")
    ap.add_argument("--deseq-cpus", type=int, default=1,
                    help="DESeq2 internal n_cpus PER worker (default 1; keep low "
                         "so n_jobs processes don't oversubscribe cores).")
    ap.add_argument("--expr-matrix", action="store_true", default=True,
                    help="Write per-sample expression matrix of DE genes (default on; "
                         "significant genes only unless --expr-all-genes).")
    ap.add_argument("--no-expr-matrix", dest="expr_matrix", action="store_false",
                    help="Skip the per-sample expression matrix.")
    ap.add_argument("--expr-all-genes", action="store_true",
                    help="Expression matrix over ALL tested genes (~25K) instead "
                         "of just significant ones. Memory-heavier, >1 GB CSV.")
    ap.add_argument("--subcluster", default=None,
                    help="Run on a 7b subcluster object (slug, e.g. 'immune'): "
                         "reads h5ad/08c_subclustered/{slug}.h5ad and uses "
                         "subcluster_name as the label. Writes *_subcluster_{slug} "
                         "outputs.")
    ap.add_argument("--rscript", default=None,
                    help="Path to Rscript (default: shutil.which('Rscript')). "
                         "Required for LRT contrasts.")
    ap.add_argument("--skip-wald", action="store_true",
                    help="Skip Wald (PyDESeq2) contrasts; only re-run LRT. If "
                         "the output CSV already exists, Wald rows are preserved "
                         "from it (single-command iteration on LRT bugs without "
                         "re-fitting ~hours of Wald).")
    ap.add_argument("--skip-lrt", action="store_true",
                    help="Skip LRT (R DESeq2) contrasts; only re-run Wald. If "
                         "the output CSV already exists, LRT rows are preserved "
                         "from it.")
    args = ap.parse_args()
    if args.skip_wald and args.skip_lrt:
        sys.exit("ERROR: --skip-wald and --skip-lrt set together leaves nothing to run.")

    print("\n=== Phase 8b: pseudobulk DE (PyDESeq2 Wald + R DESeq2 LRT) ===")
    cfg = load_config(args.config)
    tissue = cfg.get("tissue")
    if tissue not in TISSUE_CFG:
        sys.exit(f"ERROR: unknown tissue '{tissue}' (expected {list(TISSUE_CFG)}).")
    info = TISSUE_CFG[tissue]
    contrasts = load_contrasts(cfg, kind="de")

    # tuning knobs (CLI > config > defaults)
    comp_cfg = cfg.get("composition", {})
    min_donors = (args.min_donors if args.min_donors is not None
                  else int(comp_cfg.get("min_donors", 2)))
    min_donors = max(2, min_donors)
    reliable = (args.reliable_donors if args.reliable_donors is not None
                else int(comp_cfg.get("reliable_donors", 3)))
    min_samples_expr = 2
    n_jobs = max(1, args.n_jobs)
    print(f"  tissue={tissue}  min_cells/donor={args.min_cells}  "
          f"min_donors/group={min_donors}  reliable>={reliable}  workers={n_jobs}")
    print(f"  thresholds: padj<{PADJ_THR}, |log2FC|>{LFC_THR}")

    rscript = args.rscript or shutil.which("Rscript")
    if not rscript:
        print("  [warn] Rscript not found on PATH; LRT contrasts (omnibus, "
              "interaction) will be skipped. Install R + DESeq2 or pass --rscript.")
    else:
        print(f"  Rscript: {rscript}")

    # ---- locate input ----
    base = Path(cfg["results_dir"]) / "h5ad"
    if args.subcluster:
        in_path = base / "08c_subclustered" / f"{args.subcluster}.h5ad"
        if not in_path.is_file():
            sys.exit(f"ERROR: subcluster object not found: {in_path}\n"
                     f"  Run 07b_subcluster.py --celltype ... first.")
        out_suffix = f"_subcluster_{args.subcluster}"
        print(f"  SUBCLUSTER mode: {args.subcluster}")
    else:
        in_path = base / "08_annotated" / "all_samples.h5ad"
        if not in_path.is_file():
            sys.exit(f"ERROR: annotated input not found: {in_path}")
        out_suffix = ""
    print(f"  Input: {in_path}")

    plot_root = Path(cfg["results_dir"]) / "plots" / ("08b_de" + out_suffix)
    table_dir = phase_table_dir(cfg, "08b_de")
    plot_root.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(in_path)

    # ---- resolve cell-type key ----
    if args.subcluster:
        ct_key = ("subcluster_name" if "subcluster_name" in adata.obs.columns
                  else "subcluster")
    elif args.celltype_key:
        ct_key = args.celltype_key
    else:
        ct_key = info["tier_l2"]
    if ct_key not in adata.obs.columns:
        sys.exit(f"ERROR: cell-type column '{ct_key}' not in adata.obs. "
                 f"Available: {sorted(adata.obs.columns)}")
    print(f"  Cell-type column: '{ct_key}' ({adata.obs[ct_key].nunique()} types)")

    # gene symbols for readable volcano labels
    symbol_map = None
    for sym_col in ("symbol", "gene_symbol", "gene_symbols", "Symbol"):
        if sym_col in adata.var.columns:
            symbol_map = dict(zip(adata.var_names.astype(str),
                                  adata.var[sym_col].astype(str)))
            print(f"  Gene symbols for labels: var['{sym_col}']")
            break
    symbol_map = symbol_map or {}

    # .X must be raw counts
    xmax = adata.X.max()
    if not np.isclose(xmax, round(float(xmax))) or xmax < 0:
        sys.exit(f"ERROR: .X doesn't look like raw counts (max={xmax}). "
                 f"Pseudobulk needs raw counts.")

    # ---- ensure required obs columns are strings ----
    if "donor_id" not in adata.obs.columns:
        sys.exit("ERROR: obs has no 'donor_id' (the statistical unit).")
    for c in ("donor_id", "group", "age", "sex", "pool"):
        if c in adata.obs.columns:
            adata.obs[c] = adata.obs[c].astype(str)

    # ---- DROP non-cell-types: contaminants + unassigned (record per-donor) ----
    if args.subcluster:
        # subcluster_name lives directly on this object
        contam_mask = adata.obs[ct_key].astype(str).map(is_contam)
        ua_keys = [ct_key]
    else:
        subname, focal_present = join_subcluster_name(adata, base, tissue)
        adata.obs["subcluster_name"] = subname.values
        print(f"  focal types with subclusters: {focal_present or '(none)'} "
              f"({int(subname.notna().sum()):,} cells carry subcluster_name)")
        contam_mask = (adata.obs["subcluster_name"].notna()
                       & adata.obs["subcluster_name"].map(is_contam))
        ua_keys = info["ua_keys"]

    ua = unassigned_mask(adata.obs, ua_keys)
    grp_keys = ["donor_id", "age", "group"] + (
        ["sex"] if "sex" in adata.obs.columns else [])
    diag = (adata.obs.assign(_contam=contam_mask.values, _unassigned=ua.values)
            .groupby(grp_keys, observed=True)
            .agg(n_total=("_contam", "size"),
                 n_contaminant=("_contam", "sum"),
                 n_unassigned=("_unassigned", "sum")).reset_index())
    diag["frac_contaminant"] = diag.n_contaminant / diag.n_total
    diag["frac_unassigned"] = diag.n_unassigned / diag.n_total
    diag.to_csv(table_dir / f"08b_dropped_cells_per_donor{out_suffix}.csv", index=False)

    drop = (contam_mask | ua).values
    adata = adata[~drop].copy()
    print(f"  dropped {int(contam_mask.sum()):,} contaminant + "
          f"{int(ua.sum()):,} unassigned cells -> {adata.n_obs:,} cells remain")
    print(f"  dropped-cell diagnostic -> {table_dir / ('08b_dropped_cells_per_donor' + out_suffix + '.csv')}")

    # ---- sex strata ----
    sex_strata = iter_strata(cfg, axis="sex")
    if "sex" not in adata.obs.columns:
        sex_strata = [("combined", None)]
        print("  [note] no 'sex' column -> only the combined stratum.")
    print(f"  sex strata: {[s for s, _ in sex_strata]}")

    # ---- level grid (whole + each region) ----
    region_key = info["region_key"]
    has_region = bool(region_key) and region_key in adata.obs.columns \
        and adata.obs[region_key].notna().any()
    if has_region:
        region_vals = sorted(adata.obs[region_key].dropna().astype(str).unique())
        level_grid = [("whole", None)] + [(r, r) for r in region_vals]
        print(f"  level grid: whole + {len(region_vals)} regions "
              f"({', '.join(region_vals[:6])}{'...' if len(region_vals)>6 else ''})")
    else:
        level_grid = [("whole", None)]
        print(f"  level grid: whole only (no region key for {tissue})")

    # ---- BUILD JOB LIST ----
    jobs_wald, jobs_lrt = [], []
    obs = adata.obs

    for cname, spec in contrasts.items():
        kind = classify_contrast(spec)
        if kind == "skip":
            print(f"\n  [skip] {cname}: derived/non-DE contrast (handled by 8g).")
            continue
        if kind in ("omnibus_lrt", "interaction_lrt") and not rscript:
            print(f"\n  [skip] {cname}: needs Rscript for LRT.")
            continue
        if kind in ("pairwise_age_wald", "within_group_wald") and args.skip_wald:
            print(f"\n  [skip] {cname}: --skip-wald set (kind={kind}).")
            continue
        if kind in ("omnibus_lrt", "interaction_lrt") and args.skip_lrt:
            print(f"\n  [skip] {cname}: --skip-lrt set (kind={kind}).")
            continue

        flag = spec.get("flag")
        levels_disp = spec.get("levels") or spec.get("pairwise") or ("omnibus"
                       if kind == "omnibus_lrt" else "interaction")
        print(f"\n  Contrast {cname}: kind={kind} flag={flag} levels={levels_disp}")

        # Early-vs-Late only makes sense for brain (placenta arms don't share an age)
        if kind == "pairwise_age_wald" and tissue != "brain":
            lv = set(map(str, spec.get("levels", [])))
            if lv == {"Early_Stress", "Late_Stress"}:
                print(f"    [skip] Early-vs-Late impossible for {tissue}.")
                continue

        for sex_label, sex_val in sex_strata:
            m_sex = obs if sex_val is None else obs[obs["sex"] == sex_val]
            if m_sex.empty:
                continue
            ctx = dict(
                sex_label=sex_label, sex_val=sex_val,
                ct_key=ct_key, region_key=region_key, level_grid=level_grid,
                plot_root=plot_root, symbol_map=symbol_map,
                min_cells=args.min_cells, min_donors=min_donors,
                reliable_donors=reliable, min_samples_expr=min_samples_expr,
                deseq_cpus=args.deseq_cpus, rscript=rscript,
            )
            if kind == "pairwise_age_wald":
                queue_pairwise_age_wald(ctx, adata, m_sex, cname, spec, jobs_wald)
            elif kind == "within_group_wald":
                queue_within_group_wald(ctx, adata, m_sex, cname, spec, jobs_wald)
            elif kind == "omnibus_lrt":
                queue_omnibus_lrt(ctx, adata, m_sex, cname, spec, jobs_lrt)
            elif kind == "interaction_lrt":
                queue_interaction_lrt(ctx, adata, m_sex, cname, spec, jobs_lrt)

    print(f"\n  Queued: {len(jobs_wald)} Wald + {len(jobs_lrt)} LRT jobs.")

    # ---- RUN Wald (processes) ----
    rows = []
    if jobs_wald:
        print(f"\n  Running {len(jobs_wald)} Wald jobs (ProcessPool, {n_jobs} workers)...")
        for job, res, err in parallel_map(_wald_worker, jobs_wald, n_jobs=n_jobs,
                                          use_threads=False, desc="wald"):
            if err:
                last = informative_error(err)
                print(f"    [warn] {job['contrast']}|{job['sex']}|{job['group_level']}|"
                      f"{job['level']}|{job['celltype']}: DESeq2 failed: {last[:160]}")
            rows.extend(_wald_rows(job, res, err))

    # ---- RUN LRT (threads; R subprocesses are CPU-bound but release GIL) ----
    if jobs_lrt:
        print(f"\n  Running {len(jobs_lrt)} LRT jobs (ThreadPool, {n_jobs} workers)...")
        for job, res, err in parallel_map(_lrt_worker, jobs_lrt, n_jobs=n_jobs,
                                          use_threads=True, desc="lrt"):
            if err:
                last = informative_error(err)
                print(f"    [warn] {job['contrast']}|{job['sex']}|{job['group_level']}|"
                      f"{job['level']}|{job['celltype']}: LRT failed: {last[:160]}")
            rows.extend(_lrt_rows(job, res, err))

    # ---- MASTER CSV ----
    df_out = pd.DataFrame(rows)
    if len(df_out):
        df_out = df_out[[c for c in MASTER_COLS if c in df_out.columns]]
    out_csv = table_dir / f"08b_de_results{out_suffix}.csv"

    # --skip-wald / --skip-lrt: AUTO-SPLICE with the existing output CSV.
    # If the output already exists from a previous full run, preserve the rows
    # of the test method we just skipped (they're already correct) and replace
    # only the rows of the method we just re-ran. Single command, no separate
    # merge utility needed.
    if (args.skip_wald or args.skip_lrt) and out_csv.is_file():
        keep_method = "Wald" if args.skip_wald else "LRT"
        try:
            prev = pd.read_csv(out_csv)
            preserved = prev[prev["test_method"] == keep_method]
            print(f"  [--skip-{'wald' if args.skip_wald else 'lrt'}] "
                  f"preserving {len(preserved):,} {keep_method} rows from existing "
                  f"{out_csv.name}.")
            # align columns (use master order; tolerate older CSVs missing cols)
            all_cols = [c for c in MASTER_COLS if c in df_out.columns
                        or c in preserved.columns]
            df_out = pd.concat([preserved.reindex(columns=all_cols),
                                df_out.reindex(columns=all_cols)],
                               ignore_index=True)
        except Exception as e:
            print(f"  [warn] could not splice with existing {out_csv.name} "
                  f"({type(e).__name__}: {e}); writing only this run's rows.")
    df_out.to_csv(out_csv, index=False)
    n_sig = 0
    if len(df_out):
        n_sig = int(((df_out["padj"] < PADJ_THR) & (df_out["log2FC"].abs() > LFC_THR)).sum())
        n_lrt_sig = int(((df_out["test_method"] == "LRT") & (df_out["padj"] < PADJ_THR)).sum())
        print(f"\n  Master table: {out_csv}")
        print(f"    {len(df_out)} rows, {n_sig} Wald hits at padj<{PADJ_THR} & "
              f"|log2FC|>{LFC_THR}, {n_lrt_sig} LRT hits at padj<{PADJ_THR}")
    print(f"  Volcanoes: {plot_root}")

    # ---- PER-SAMPLE EXPRESSION MATRIX (offline audit; same as before) ----
    # Skip when --skip-wald is set: the gene universe would be picked from LRT
    # rows only, whose log2FC is NaN by design (LRT tests multiple coefs
    # jointly), so the |log2FC|>1 filter rejects everything. Re-use the
    # expression matrix from the previous full run instead.
    if args.expr_matrix and len(df_out) and not args.skip_wald:
        dfh = df_out.dropna(subset=["gene"])
        if args.expr_all_genes:
            gene_universe = sorted(dfh["gene"].astype(str).unique())
            gtag = "all tested genes"
        else:
            sel = dfh[(dfh["padj"] < PADJ_THR) & (dfh["log2FC"].abs() > LFC_THR)]
            gene_universe = sorted(sel["gene"].astype(str).unique())
            gtag = f"significant Wald DEGs (padj<{PADJ_THR}, |LFC|>{LFC_THR})"
        meta_cols = [c for c in ("group", "age", "sex", "pool")
                     if c in adata.obs.columns]
        if gene_universe:
            long_df, sample_meta = build_expression_matrix(
                adata, ct_key, gene_universe, meta_cols)
            if symbol_map:
                long_df["gene_symbol"] = long_df["gene"].map(
                    lambda g: symbol_map.get(g, g))
            expr_csv = table_dir / f"08b_de_gene_expression_per_sample{out_suffix}.csv"
            meta_csv = table_dir / f"08b_sample_metadata{out_suffix}.csv"
            long_df.to_csv(expr_csv, index=False)
            sample_meta.to_csv(meta_csv)
            print(f"  Expression matrix: {expr_csv}")
            print(f"    {len(gene_universe)} {gtag} × {long_df['celltype'].nunique()} "
                  f"cell types × {sample_meta.shape[0]} samples ({len(long_df)} rows)")
            print(f"  Sample metadata:   {meta_csv}")

    print("\n✓ Phase 8b complete.")
    print("  Pup = statistical unit; no dam ID (anti-conservative); n small. "
          f"Trust only large effects (|log2FC|>{LFC_THR}, padj<{PADJ_THR}).")
    print("  Read 'test_method' (Wald|LRT), 'sex', 'reliability', 'flag', 'note' "
          "on every row. LRT rows have NaN log2FC by design.")
    print(f"  Cross-slice summary plots are produced by 08b_de_summary.py "
          f"(reads {out_csv.name}).\n")


if __name__ == "__main__":
    main()
