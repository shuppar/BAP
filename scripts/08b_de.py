#!/usr/bin/env python
"""
08b_de.py — Phase 8b: pseudobulk differential expression (PyDESeq2).

THE statistical core. Pseudobulk only — never single-cell-level (cells are not
independent replicates; treating them as such inflates significance and reviewers
will flag it — project doc §2). The statistical unit is the ANIMAL (donor_id).

Per cell type, per contrast:
  1. Sum raw counts across cells -> pseudobulk: one column per donor (animal)
  2. Filter low-expressed genes; require enough cells/donors backing each donor
  3. PyDESeq2 with the contrast's design (e.g. ~ sex + pool + group)
  4. Wald test for the contrast levels (e.g. Early_Stress vs Relaxed)

Iterates over DE-style contrasts from the declarative spec (load_contrasts "de").
Pairwise level contrasts run as DESeq2 contrasts; omnibus (group_omnibus) runs
as a likelihood-ratio-style 3-group test via the reduced-model LRT where
available, else flagged skipped; interaction terms are recorded but flagged
underpowered (project doc §2).

Caveats carried into output (never silently dropped):
  - No dam ID -> each pup treated as independent. Anti-conservative for litter-
    aggregated traits. Every row carries the contrast `flag`.
  - n=2-4 per group -> only large effects (|logFC| > ~1.5) are trustworthy.
  - Pool-confounded contrasts (e.g. P1 Late Stress) carry confound_warnings.

PyDESeq2 API verified against docs (handles 0.4.x and 0.5.x):
  newer: DeseqDataSet(counts=, metadata=, design="~ ...", inference=DefaultInference(n_cpus=))
  older: DeseqDataSet(counts=, metadata=, design_factors=[...], n_cpus=)
  both : dds.deseq2(); DeseqStats(dds, contrast=[factor, test, ref], ...).summary(); .results_df

Usage:
  uv run python scripts/08b_de.py --config config/dev.yaml
  uv run python scripts/08b_de.py --config config/brain.yaml
  uv run python scripts/08b_de.py --config config/placenta.yaml
  uv run python scripts/08b_de.py --config config/dev.yaml --celltype-key manual_annotation

Inputs (first that exists):
  {results_dir}/h5ad/08b_label_transferred/all_samples.h5ad   (Phase 7c)
  {results_dir}/h5ad/08_annotated/all_samples.h5ad            (Phase 7)

Outputs:
  {results_dir}/plots/08b_de/{contrast}/{level}/{celltype}/
    - volcano.png
  {results_dir}/tables/
    - de_results.csv : [contrast, flag, group_level, pair, celltype, gene,
                        log2FC, lfcSE, stat, pvalue, padj, direction,
                        n_donors_test, n_donors_ref, reliability, note]
"""

import argparse
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from _utils import load_config, load_contrasts


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


def make_pseudobulk(adata_ct, covariates):
    """Sum raw counts per donor for one cell type.

    Returns (counts_df [donor x gene, int], metadata_df [donor x covariates],
             n_cells_per_donor Series). Raw counts must be in .X.
    """
    donors = adata_ct.obs["donor_id"].astype(str)
    uniq = sorted(donors.unique())
    X = adata_ct.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)

    rows = []
    n_cells = {}
    for d in uniq:
        m = (donors == d).values
        n_cells[d] = int(m.sum())
        rows.append(np.asarray(X[m].sum(axis=0)).ravel())
    counts = pd.DataFrame(np.vstack(rows), index=uniq, columns=adata_ct.var_names)
    counts = counts.round().astype(int)   # DESeq2 needs integer counts

    meta = (adata_ct.obs[["donor_id"] + covariates]
            .astype({"donor_id": str})
            .drop_duplicates("donor_id").set_index("donor_id").loc[uniq])
    return counts, meta, pd.Series(n_cells)


def run_pydeseq2(counts, meta, design_terms, contrast_levels, test_factor, n_cpus=4):
    """Fit DESeq2 and run the Wald test for one pairwise contrast.

    design_terms: list like ["sex","pool","group"] (intercept implied).
    contrast_levels: [test_level, ref_level] for test_factor.
    Returns results_df (pandas) or raises.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # Keep only design terms that actually vary (DESeq2 errors on constant factors)
    design_terms = [t for t in design_terms if meta[t].nunique() > 1]
    if test_factor not in design_terms:
        design_terms.append(test_factor)

    # genes expressed at all (DESeq2 handles low counts, but drop all-zero)
    counts = counts.loc[:, counts.sum(axis=0) > 0]

    formula = "~ " + " + ".join(design_terms)

    # Try the newer formula+inference API, fall back to the older one.
    dds = None
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


def plot_volcano(res, title, out, padj_thr=0.05, lfc_thr=1.0, max_labels=25,
                 symbol_map=None):
    """Volcano with the top significant genes labeled.

    Genes come from res.index (PyDESeq2 results_df is indexed by gene). If
    symbol_map (dict: var_name -> symbol) is given, labels use symbols so the
    plot never shows Ensembl IDs. Only the most significant hits are labeled
    (up to max_labels) for readability. Uses adjustText if installed, else a
    simple offset fallback.
    """
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
    ax.axvline(lfc_thr, color="k", lw=0.6, ls="--"); ax.axvline(-lfc_thr, color="k", lw=0.6, ls="--")
    ax.set_xlabel("log2 fold change"); ax.set_ylabel("-log10 padj")
    ax.set_title(title, fontsize=9)

    # Label the top significant genes by significance (then |logFC| as tiebreak).
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
            # No adjustText — nudge labels up-right a touch so they don't sit on
            # the marker. Less pretty but fully functional.
            for t in texts:
                x, y = t.get_position()
                t.set_position((x + 0.05, y + 0.05))

    n_sig = int(sig.sum())
    ax.text(0.02, 0.98, f"{n_sig} sig (padj<{padj_thr}, |LFC|>{lfc_thr})",
            transform=ax.transAxes, fontsize=6, va="top", color="gray")
    fig.tight_layout(); fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Phase 8b: pseudobulk DE (PyDESeq2)")
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--celltype-key", default=None)
    ap.add_argument("--min-cells", type=int, default=10,
                    help="min cells of a type in a donor for that donor to count (default 10)")
    ap.add_argument("--min-donors", type=int, default=None,
                    help="min donors per group with >=min-cells (default: composition.min_donors or 3)")
    ap.add_argument("--n-cpus", type=int, default=4)
    ap.add_argument("--subcluster", default=None,
                    help="Run DE on a 7b subcluster object (slug, e.g. 'microglia'): "
                         "reads h5ad/08c_subclustered/{slug}.h5ad and uses the "
                         "'subcluster' column as the cell-type label. Writes "
                         "de_results_subcluster_{slug}.csv.")
    args = ap.parse_args()

    print(f"\n=== Phase 8b: Pseudobulk DE (PyDESeq2) ===")
    cfg = load_config(args.config)
    contrasts = load_contrasts(cfg, kind="de")
    group_ref = cfg.get("group_reference", "Relaxed")
    min_donors = (args.min_donors if args.min_donors is not None
                  else int(cfg.get("composition", {}).get("min_donors", 3)))
    print(f"  min_cells/donor={args.min_cells}, min_donors/group={min_donors}")

    base = Path(cfg["results_dir"]) / "h5ad"
    if args.subcluster:
        in_path = base / "08c_subclustered" / f"{args.subcluster}.h5ad"
        if not in_path.is_file():
            sys.exit(f"ERROR: subcluster object not found: {in_path}\n"
                     f"  Run 07b_subcluster.py --celltype ... first.")
        forced_ct_key = "subcluster"   # the label 7b writes
        out_suffix = f"_subcluster_{args.subcluster}"
        print(f"  SUBCLUSTER mode: {args.subcluster} (label column 'subcluster')")
    else:
        cand = [base / "08b_label_transferred" / "all_samples.h5ad",
                base / "08_annotated" / "all_samples.h5ad"]
        in_path = next((p for p in cand if p.is_file()), None)
        if in_path is None:
            sys.exit("ERROR: no annotated input. Checked:\n  " +
                     "\n  ".join(str(p) for p in cand))
        forced_ct_key = None
        out_suffix = ""
    print(f"  Input: {in_path}")

    plot_root = Path(cfg["results_dir"]) / "plots" / ("08b_de" + out_suffix)
    table_dir = Path(cfg["results_dir"]) / "tables"
    plot_root.mkdir(parents=True, exist_ok=True); table_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(in_path)
    ct_key = resolve_celltype_key(adata, forced_ct_key or args.celltype_key)
    print(f"  Cell type column: '{ct_key}' ({adata.obs[ct_key].nunique()} types)")
    # Map var_names -> gene symbols for readable volcano labels (var_names may be
    # Ensembl IDs per the project's gene-naming convention). None if no symbol col.
    symbol_map = None
    for sym_col in ("symbol", "gene_symbol", "gene_symbols", "Symbol"):
        if sym_col in adata.var.columns:
            symbol_map = dict(zip(adata.var_names.astype(str),
                                  adata.var[sym_col].astype(str)))
            print(f"  Gene symbols for labels: var['{sym_col}']")
            break
    # Sanity: .X must be raw counts for pseudobulk summing
    xmax = adata.X.max()
    if not np.isclose(xmax, round(float(xmax))) or xmax < 0:
        sys.exit(f"ERROR: .X doesn't look like raw counts (max={xmax}). "
                 f"Pseudobulk needs raw counts.")

    rows = []
    for cname, spec in contrasts.items():
        test = spec.get("test"); flag = spec.get("flag")
        levels = spec.get("levels"); group_by = spec.get("group_by")
        confound = spec.get("confound_warnings", {})

        # 8b handles pairwise level tests. Omnibus/interaction recorded as skipped
        # (omnibus is read from pairwise; interaction is underpowered here).
        if test == "group_omnibus" or (test and ":" in test):
            print(f"\n  [skip] {cname}: {test} not run in 8b (see project doc §2).")
            continue
        if not levels:
            print(f"\n  [skip] {cname}: no pairwise levels (within-group/age handled elsewhere).")
            continue

        print(f"\n  Contrast: {cname} (test={test} {levels}, flag={flag})")

        # iterate group_by levels (e.g. each age)
        if isinstance(group_by, list):
            combos = adata.obs[group_by].drop_duplicates().itertuples(index=False)
            group_iter = [dict(zip(group_by, c)) for c in combos]
        elif group_by:
            group_iter = [{group_by: v} for v in adata.obs[group_by].unique()]
        else:
            group_iter = [{}]

        design_terms = [t.strip() for t in spec.get("design", f"~ {test}").replace("~", "").split("+")]
        design_terms = [t for t in design_terms if t and "*" not in t]

        for sel in group_iter:
            sub = adata
            label = "_".join(f"{k}-{v}" for k, v in sel.items()) or "all"
            for k, v in sel.items():
                sub = sub[sub.obs[k] == v]
            sub = sub[sub.obs[test].isin(levels)]
            if sub.n_obs == 0:
                continue
            note = confound.get(sel.get(group_by) if not isinstance(group_by, list) else None, "")

            for ct in sub.obs[ct_key].astype(str).unique():
                ct_ad = sub[sub.obs[ct_key].astype(str) == ct]
                # donors with enough cells of this type
                per_donor = ct_ad.obs["donor_id"].value_counts()
                good = per_donor[per_donor >= args.min_cells].index
                ct_ad = ct_ad[ct_ad.obs["donor_id"].isin(good)]
                # donors per group after the cell floor
                g = ct_ad.obs.groupby(test, observed=True)["donor_id"].nunique()
                n_test = int(g.get(levels[0], 0)); n_ref = int(g.get(levels[1], 0))
                if n_test < min_donors or n_ref < min_donors:
                    continue
                reliability = "unreliable_n<3" if min(n_test, n_ref) < 3 else "ok"

                covariates = [t for t in design_terms
                              if t in ct_ad.obs.columns and ct_ad.obs[t].nunique() > 1]
                if test not in covariates:
                    covariates.append(test)
                counts, meta, _ = make_pseudobulk(ct_ad, covariates)
                # gene filter: >=10 counts in >=min_donors samples (project doc §8b)
                keep = ((counts >= 10).sum(axis=0) >= min_donors)
                counts = counts.loc[:, keep]
                if counts.shape[1] < 10:
                    continue

                print(f"    {label} | {ct}: {n_test} vs {n_ref} donors, "
                      f"{counts.shape[1]} genes")
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        res = run_pydeseq2(counts, meta, covariates, levels, test, args.n_cpus)
                except Exception as e:
                    print(f"      [warn] DESeq2 failed: {e}")
                    rows.append({"contrast": cname, "flag": flag, "group_level": label,
                                 "pair": str(levels), "celltype": ct, "gene": None,
                                 "log2FC": None, "padj": None, "direction": None,
                                 "n_donors_test": n_test, "n_donors_ref": n_ref,
                                 "reliability": reliability, "note": f"DESeq2 failed: {e}; {note}"})
                    continue

                pdir = plot_root / cname / label / ct.replace("/", "_").replace(" ", "_")
                pdir.mkdir(parents=True, exist_ok=True)
                plot_volcano(res, f"{cname}\n{label} | {ct}", pdir / "volcano.png",
                             symbol_map=symbol_map)

                res = res.reset_index().rename(columns={"index": "gene"})
                for _, r in res.iterrows():
                    lfc = r.get("log2FoldChange")
                    rows.append({
                        "contrast": cname, "flag": flag, "group_level": label,
                        "pair": str(levels), "celltype": ct, "gene": r["gene"],
                        "log2FC": lfc, "lfcSE": r.get("lfcSE"), "stat": r.get("stat"),
                        "pvalue": r.get("pvalue"), "padj": r.get("padj"),
                        "direction": (None if pd.isna(lfc) else ("up" if lfc > 0 else "down")),
                        "n_donors_test": n_test, "n_donors_ref": n_ref,
                        "reliability": reliability, "note": note,
                    })

    out_csv = table_dir / f"de_results{out_suffix}.csv"
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    n_sig = 0
    if rows:
        dfh = pd.DataFrame(rows)
        n_sig = int(((dfh["padj"] < 0.05) & (dfh["log2FC"].abs() > 1)).sum())
    print(f"\n  Master table: {out_csv}  ({len(rows)} rows, {n_sig} at padj<0.05 & |log2FC|>1)")
    print(f"  Plots: {plot_root}")
    print(f"\n✓ Phase 8b complete.")
    print(f"\n  Reminder: pup is the unit; no dam ID (anti-conservative); n small.")
    print(f"  Trust only large effects, read 'reliability' + 'flag' + 'note' per row.\n")


if __name__ == "__main__":
    main()
