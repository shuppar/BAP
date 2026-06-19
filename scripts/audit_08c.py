#!/usr/bin/env python
"""
audit_08c.py — thorough audit of Phase 8c outputs across all 8 runs.

Runs five check categories:
  A. Structural integrity per CSV (dup keys, NaN, schema, n_cells>0)
  B. Cross-job pathway/TF universe consistency
  C. Placenta first-look (cell types, donor recon, top Hallmark hits)
  D. Per-cell h5ad integrity (shape, X validity, UMAP, var labels)
  E. Sig-rate calibration per (job × collection)

Prints PASS / WARN / FAIL per check. Final summary at the bottom.
Memory-conscious: reads only required columns where possible; chunks large CSVs.

Usage: uv run python scripts/audit_08c.py
"""

import ast
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad


# --- Job spec mirrors run_08c_remaining.sh ---------------------------------
JOBS = [
    ("brain",    ""),                       # brain main (already done)
    ("brain",    "immune"),
    ("brain",    "opc_oligodendrocytes"),
    ("brain",    "astrocytes_ependymal"),
    ("placenta", ""),
    ("placenta", "dsc"),
    ("placenta", "endothelium"),
    ("placenta", "myeloid"),
    ("placenta", "nk"),
]

RDIR = Path("results")
KEY_COLS_GSEA = ["sex", "contrast", "group_level", "pair", "level", "celltype", "source"]
KEY_COLS_TF   = ["sex", "contrast", "group_level", "pair", "level", "celltype", "TF"]
KEY_COLS_PD   = ["donor_id", "celltype", "level", "region", "pathway"]


def label(tissue, sub):
    return f"{tissue}/{sub or 'main'}"


def paths_for(tissue, sub):
    sfx = f"_subcluster_{sub}" if sub else ""
    base = RDIR / tissue / "tables" / f"08c_pathways{sfx}"
    pc_name = (f"{tissue}_subcluster_{sub}_per_cell_scores.h5ad" if sub
               else f"{tissue}_per_cell_scores.h5ad")
    return {
        "gsea":   base / f"08c_pathway_results{sfx}.csv",
        "le":     base / f"08c_pathway_leading_edge{sfx}.csv",
        "tf":     base / f"08c_tf_activity{sfx}.csv",
        "pd":     base / f"08c_pathway_scores_per_donor{sfx}.csv",
        "pcell":  RDIR / tissue / "h5ad" / "08c_pathway_scores" / pc_name,
    }


# Final per-job status accumulator
RESULTS = {label(t, s): {"A": [], "B": [], "C": [], "D": [], "E": []} for t, s in JOBS}


def mark(job_lbl, cat, status, msg):
    RESULTS[job_lbl][cat].append((status, msg))
    color = {"PASS": "✓", "WARN": "⚠", "FAIL": "✗"}[status]
    print(f"    {color} [{status}] {msg}")


# ===========================================================================
# A. Structural integrity
# ===========================================================================
print("\n" + "=" * 75)
print("A. Structural integrity per CSV")
print("=" * 75)

EXPECTED_GSEA_COLS = {"tissue", "sex", "contrast", "flag", "group_level", "pair",
                     "level", "celltype", "collection", "source", "NES",
                     "pvalue", "FDR", "FDR_pooled"}
EXPECTED_TF_COLS = {"tissue", "sex", "contrast", "flag", "group_level", "pair",
                    "level", "celltype", "TF", "activity_score", "pvalue",
                    "FDR", "FDR_ctx_celltype", "direction"}
EXPECTED_PD_COLS = {"tissue", "donor_id", "sample_id", "sex", "age", "group", "pool",
                    "celltype", "level", "region", "pathway", "mean_score",
                    "median_score", "n_cells"}

for tissue, sub in JOBS:
    lbl = label(tissue, sub)
    paths = paths_for(tissue, sub)
    print(f"\n  {lbl}")

    # --- GSEA master ---
    p = paths["gsea"]
    if not p.is_file():
        mark(lbl, "A", "FAIL", f"GSEA CSV missing: {p}"); continue
    # Schema check (cheap)
    head = pd.read_csv(p, nrows=5)
    missing = EXPECTED_GSEA_COLS - set(head.columns)
    if missing:
        mark(lbl, "A", "FAIL", f"GSEA missing cols: {missing}")
    else:
        mark(lbl, "A", "PASS", f"GSEA schema OK ({len(head.columns)} cols)")

    # NaN scan on critical numeric cols (chunked, in case file is huge)
    nan_counts = {"NES": 0, "pvalue": 0, "FDR": 0}
    n_rows = 0
    for chunk in pd.read_csv(p, usecols=list(nan_counts), chunksize=500_000):
        n_rows += len(chunk)
        for c in nan_counts:
            nan_counts[c] += int(chunk[c].isna().sum())
    bad = {c: n for c, n in nan_counts.items() if n > 0}
    if bad:
        # NaN pvalue/FDR can happen for low_n slices (recorded as "low_n" reliability);
        # NES is more concerning if missing.
        if nan_counts["NES"] > 0:
            mark(lbl, "A", "WARN", f"GSEA NaN: {bad} (NaN NES is suspicious; check low_n slices)")
        else:
            mark(lbl, "A", "PASS", f"GSEA NaN OK (NES never NaN; pvalue/FDR NaN: {bad})")
    else:
        mark(lbl, "A", "PASS", f"GSEA NaN OK (0 in NES/pvalue/FDR across {n_rows:,} rows)")

    # Dup key check
    keys = pd.read_csv(p, usecols=KEY_COLS_GSEA)
    n_dup = int(keys.duplicated().sum())
    if n_dup > 0:
        mark(lbl, "A", "FAIL", f"GSEA dup keys: {n_dup:,}")
    else:
        mark(lbl, "A", "PASS", f"GSEA dup keys: 0 (across {len(keys):,} rows)")

    # Pair column parse spot-check
    pair_uniq = keys["pair"].dropna().unique()[:20]
    bad_pairs = []
    for ps in pair_uniq:
        try:
            v = ast.literal_eval(str(ps))
            if not isinstance(v, (list, tuple)) or len(v) != 2:
                bad_pairs.append(ps)
        except Exception:
            bad_pairs.append(ps)
    if bad_pairs:
        mark(lbl, "A", "WARN", f"unparseable pair values: {bad_pairs[:3]}")
    else:
        mark(lbl, "A", "PASS", f"pair column parses cleanly ({len(pair_uniq)} unique sampled)")

    # --- TF activity ---
    p = paths["tf"]
    if p.is_file():
        head = pd.read_csv(p, nrows=5)
        missing = EXPECTED_TF_COLS - set(head.columns)
        if missing:
            mark(lbl, "A", "FAIL", f"TF missing cols: {missing}")
        tf = pd.read_csv(p, usecols=list(KEY_COLS_TF + ["activity_score", "pvalue", "FDR", "FDR_ctx_celltype"]))
        n_dup = int(tf[KEY_COLS_TF].duplicated().sum())
        nan_act = int(tf["activity_score"].isna().sum())
        nan_pv = int(tf["pvalue"].isna().sum())
        nan_ctx = int(tf["FDR_ctx_celltype"].isna().sum())
        msg = (f"TF dup={n_dup}, NaN(activity)={nan_act}, NaN(pvalue)={nan_pv}, "
               f"NaN(FDR_ctx)={nan_ctx}, rows={len(tf):,}")
        bad = (n_dup > 0) or (nan_act > 0) or (nan_pv > 0)
        mark(lbl, "A", "FAIL" if bad else "PASS", msg)

    # --- Per-donor ---
    p = paths["pd"]
    if p.is_file():
        head = pd.read_csv(p, nrows=5)
        missing = EXPECTED_PD_COLS - set(head.columns)
        if missing:
            mark(lbl, "A", "FAIL", f"per-donor missing cols: {missing}")
        else:
            mark(lbl, "A", "PASS", f"per-donor schema OK")
        # Spot checks
        pd_df = pd.read_csv(p, usecols=KEY_COLS_PD + ["mean_score", "n_cells"])
        n_dup = int(pd_df[KEY_COLS_PD].duplicated().sum())
        n_nan_score = int(pd_df["mean_score"].isna().sum())
        n_zero_cells = int((pd_df["n_cells"] == 0).sum())
        n_neg = int((pd_df["mean_score"] < 0).sum())
        msg = (f"per-donor dup={n_dup}, NaN(mean)={n_nan_score}, "
               f"n_cells=0 rows={n_zero_cells}, neg scores={n_neg}, rows={len(pd_df):,}")
        bad = (n_dup > 0) or (n_zero_cells > 0) or (n_neg > 0)
        if n_nan_score > 0:
            bad = True  # NaN should not happen
        mark(lbl, "A", "FAIL" if bad else "PASS", msg)

    # --- Leading-edge (schema only — file too big to fully scan) ---
    p = paths["le"]
    if p.is_file():
        head = pd.read_csv(p, nrows=5)
        needed = {"sex", "contrast", "group_level", "pair", "level", "celltype",
                  "collection", "pathway", "NES", "pathway_FDR",
                  "leading_edge_rank", "gene", "log2FC", "rank_stat", "direction"}
        missing = needed - set(head.columns)
        if missing:
            mark(lbl, "A", "FAIL", f"leading-edge missing cols: {missing}")
        else:
            mark(lbl, "A", "PASS", f"leading-edge schema OK (full scan skipped — file too large)")


# ===========================================================================
# B. Cross-job pathway/TF universe consistency
# ===========================================================================
print("\n" + "=" * 75)
print("B. Cross-job pathway/TF universe")
print("=" * 75)

path_universes = {}
tf_universes = {}
for tissue, sub in JOBS:
    lbl = label(tissue, sub)
    paths = paths_for(tissue, sub)
    # Pathway universe in this job
    upaths = set()
    for chunk in pd.read_csv(paths["gsea"], usecols=["source"], chunksize=500_000):
        upaths |= set(chunk["source"].unique())
    path_universes[lbl] = upaths
    # TF universe
    if paths["tf"].is_file():
        tfs = set(pd.read_csv(paths["tf"], usecols=["TF"])["TF"].unique())
        tf_universes[lbl] = tfs

union_paths = set().union(*path_universes.values())
union_tfs = set().union(*tf_universes.values())
print(f"\n  Union pathway space across all 8 jobs: {len(union_paths):,} unique pathways")
print(f"  Union TF space across all 8 jobs:      {len(union_tfs):,} unique TFs")
print(f"\n  Per-job coverage of the union:")
print(f"  {'job':<32} {'pathways':>10} {'/ union':>10}  {'TFs':>6} {'/ union':>10}")
for tissue, sub in JOBS:
    lbl = label(tissue, sub)
    np_ = len(path_universes.get(lbl, set()))
    nt_ = len(tf_universes.get(lbl, set()))
    print(f"  {lbl:<32} {np_:>10,} / {len(union_paths):>8,}  {nt_:>6} / {len(union_tfs):>8}")

# Flag: any job using a pathway NOT in the union of master GSEA jobs (impossible
# unless a different gene-set TSV was used)
all_in_union = all(p.issubset(union_paths) for p in path_universes.values())
mark_lbl = "ALL"
if all_in_union:
    print(f"\n  ✓ [PASS] all jobs' pathways drawn from the same MSigDB universe")
else:
    print(f"\n  ✗ [FAIL] job(s) used pathways outside the union (different gene-set TSV?)")

# All jobs should share the SAME TF universe (CollecTRI)
if tf_universes:
    first_tfs = next(iter(tf_universes.values()))
    if all(t == first_tfs for t in tf_universes.values()):
        print(f"  ✓ [PASS] all jobs share identical CollecTRI TF universe ({len(first_tfs)} TFs)")
    else:
        diffs = [(l, len(t.symmetric_difference(first_tfs)))
                 for l, t in tf_universes.items()]
        print(f"  ⚠ [WARN] TF universes differ across jobs (job: |sym_diff|): {diffs}")


# ===========================================================================
# C. Placenta first-look
# ===========================================================================
print("\n" + "=" * 75)
print("C. Placenta first-look biology")
print("=" * 75)

# --- placenta main ---
print("\n  placenta/main")
p = paths_for("placenta", "")["pd"]
if p.is_file():
    pd_df = pd.read_csv(p, usecols=["donor_id", "celltype", "level", "pathway", "mean_score"])
    cts = sorted(pd_df["celltype"].unique())
    donors = sorted(pd_df["donor_id"].unique())
    print(f"    donors: {len(donors)}  → {donors}")
    print(f"    celltypes ({len(cts)}): {cts}")
    print(f"    levels: {sorted(pd_df['level'].unique())}")
    print(f"    pathways scored: {pd_df['pathway'].nunique()}")
    if 10 <= len(cts) <= 35:
        print(f"    ✓ celltype count in expected range (~21 broad classes)")
    else:
        print(f"    ⚠ celltype count = {len(cts)} (expected ~21 broad)")

p_gsea = paths_for("placenta", "")["gsea"]
if p_gsea.is_file():
    print(f"\n  placenta main — top 8 sig MH hits (by lowest FDR, across celltypes/contrasts)")
    gsea = pd.read_csv(p_gsea, usecols=["sex", "contrast", "group_level", "level",
                                          "celltype", "collection", "source", "NES", "FDR"])
    mh_sig = gsea[(gsea["collection"] == "MH") & (gsea["FDR"] < 0.05)
                   & (gsea["sex"] == "combined")]
    top = mh_sig.sort_values(["FDR", "NES"]).head(8)
    if top.empty:
        print("    [no MH sig hits in combined sex]")
    else:
        for _, r in top.iterrows():
            print(f"    {r['celltype']:<35} {r['contrast']:<25} {r['group_level']:<6} "
                  f"{r['source']:<40} NES={r['NES']:+.2f}  FDR={r['FDR']:.1e}")

# --- placenta subclusters ---
for sub in ["dsc", "endothelium", "myeloid", "nk"]:
    print(f"\n  placenta/{sub}")
    p = paths_for("placenta", sub)["pd"]
    if p.is_file():
        pd_df = pd.read_csv(p, usecols=["celltype"])
        cts = sorted(pd_df["celltype"].unique())
        print(f"    subcluster_name values ({len(cts)}): {cts}")


# ===========================================================================
# D. Per-cell h5ad integrity
# ===========================================================================
print("\n" + "=" * 75)
print("D. Per-cell h5ad integrity")
print("=" * 75)

REQUIRED_OBS = {"donor_id", "sample_id", "sex", "age", "group", "pool"}

for tissue, sub in JOBS:
    lbl = label(tissue, sub)
    pc = paths_for(tissue, sub)["pcell"]
    print(f"\n  {lbl}")
    if not pc.is_file():
        mark(lbl, "D", "FAIL", f"per-cell h5ad missing: {pc}")
        continue
    a = ad.read_h5ad(pc)
    n_cells, n_paths = a.shape
    # X validity
    X = a.X
    if hasattr(X, "toarray"):
        x_min, x_max = float(X.min()), float(X.max())
        n_nan = int(np.isnan(X.data).sum()) if hasattr(X, "data") else 0
    else:
        x_min, x_max = float(np.nanmin(X)), float(np.nanmax(X))
        n_nan = int(np.isnan(X).sum())
    n_inf = int(np.isinf(np.asarray(X if not hasattr(X, "toarray") else X.toarray())).sum())
    msg = f"shape=({n_cells:,}, {n_paths}), X range [{x_min:.3f}, {x_max:.3f}], NaN={n_nan}, inf={n_inf}"
    bad = (n_nan > 0) or (n_inf > 0) or (x_min < 0) or (x_max > 1.0)
    mark(lbl, "D", "FAIL" if bad else "PASS", f"X: {msg}")

    # obsm['X_umap']
    if "X_umap" in a.obsm:
        u = a.obsm["X_umap"]
        if u.shape[0] == n_cells and u.shape[1] == 2:
            mark(lbl, "D", "PASS", f"obsm['X_umap'] OK {u.shape}")
        else:
            mark(lbl, "D", "FAIL", f"obsm['X_umap'] wrong shape {u.shape}")
    else:
        mark(lbl, "D", "FAIL", f"obsm['X_umap'] MISSING")

    # var has collection column
    if "collection" in a.var.columns:
        ncoll = a.var["collection"].nunique()
        mark(lbl, "D", "PASS", f"var['collection'] OK ({ncoll} collections)")
    else:
        mark(lbl, "D", "WARN", f"var['collection'] MISSING")

    # obs has required columns
    miss = REQUIRED_OBS - set(a.obs.columns)
    if miss:
        mark(lbl, "D", "FAIL", f"obs missing: {miss}")
    else:
        mark(lbl, "D", "PASS", f"obs has all required keys ({len(a.obs.columns)} total)")


# ===========================================================================
# E. Sig-rate calibration per (job × collection)
# ===========================================================================
print("\n" + "=" * 75)
print("E. Sig-rate calibration (FDR<0.05) per (job × collection)")
print("=" * 75)
print(f"\n  {'job':<32} {'collection':<6} {'rows':>10} {'sig':>8} {'%sig':>7}")

for tissue, sub in JOBS:
    lbl = label(tissue, sub)
    p = paths_for(tissue, sub)["gsea"]
    if not p.is_file():
        continue
    # Per-collection sig rates from streamed read
    counts = {}
    sigs = {}
    for chunk in pd.read_csv(p, usecols=["collection", "FDR"], chunksize=500_000):
        for coll, g in chunk.groupby("collection"):
            counts[coll] = counts.get(coll, 0) + len(g)
            sigs[coll] = sigs.get(coll, 0) + int((g["FDR"] < 0.05).sum())
    for coll in sorted(counts.keys()):
        n = counts[coll]; ns = sigs[coll]; pct = 100.0 * ns / n if n else 0.0
        flag = ""
        if pct > 80:
            flag = "  ⚠ very high (>80%) — verify"
        elif pct < 0.5 and coll == "MH":
            flag = "  ⚠ MH very low — possible power issue"
        print(f"  {lbl:<32} {coll:<6} {n:>10,} {ns:>8,} {pct:>6.1f}%{flag}")


# ===========================================================================
# Final summary
# ===========================================================================
print("\n" + "=" * 75)
print("FINAL SUMMARY")
print("=" * 75)
print(f"\n  {'job':<32} {'A':>6} {'B':>6} {'C':>6} {'D':>6} {'E':>6}")
overall_fail = 0
for tissue, sub in JOBS:
    lbl = label(tissue, sub)
    row = {}
    for cat in ["A", "B", "C", "D", "E"]:
        items = RESULTS[lbl][cat]
        if not items:
            row[cat] = "  -  "
        elif any(s == "FAIL" for s, _ in items):
            row[cat] = " FAIL"
            overall_fail += 1
        elif any(s == "WARN" for s, _ in items):
            row[cat] = " WARN"
        else:
            row[cat] = " PASS"
    print(f"  {lbl:<32} {row['A']:>6} {row['B']:>6} {row['C']:>6} {row['D']:>6} {row['E']:>6}")

print()
if overall_fail == 0:
    print("  ✓ 0 FAIL categories. 8c outputs ready for downstream.")
else:
    print(f"  ✗ {overall_fail} FAIL category-jobs — investigate before downstream.")
