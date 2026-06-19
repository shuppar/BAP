#!/usr/bin/env python
"""
08c_pathways.py — Phase 8c: pathway / TF activity / per-cell scoring on 8b DE.

Operates entirely on the 8b master DE CSV (and the subcluster variant under
--subcluster). Produces FOUR CSVs + ONE per-cell h5ad — NO plots. Plots live
in scripts/08c_pathways_summary.py (mirrors the 8b / 8b_summary split).

Iteration matches 8b's iteration. Work unit (one GSEA + one TF ULM call):
    (sex, contrast, group_level, pair, level, celltype)   on Wald rows
LRT rows (test_method=="LRT", log2FC=NaN) are skipped — no Wald stat to rank.

Parallelism is mandatory (project standard). Threads via _utils.parallel_map —
GSEA goes through an R subprocess (fgsea-multilevel; see run_gsea_on_ranks);
ULM regression is numpy/numba (GIL released).

Outputs (under {results_dir}/tables/08c_pathways{suffix}/):
  08c_pathway_results.csv            (GSEA master)
  08c_pathway_results_<COLLECTION>.csv     (paper-panel splits: MH/M2/M5/M8)
  08c_pathway_leading_edge.csv       (genes driving sig pathways)
  08c_tf_activity.csv                (CollecTRI ULM, BH within slice)
  08c_pathway_scores_per_donor.csv   (per-donor pathway means; what 8f/8g eat)

And one h5ad ({results_dir}/h5ad/08c_pathway_scores/):
  {tissue}{suffix}_per_cell_scores.h5ad  (cells × pathways AUCell, + UMAP)

Master CSV row schema (mirrors 8b keys so 8f/8g join cleanly):
  tissue, sex, contrast, flag, group_level, pair, level, celltype,
  collection, source, NES, pvalue, FDR, FDR_pooled, n_donors_total,
  reliability, note

Per-collection BH FDR is the headline FDR (keeps the 50 Hallmark sets from
being buried under thousands of GO:BP sets). FDR_pooled (BH across all sets
in the slice) is kept as a reference column only.

TF activity (CollecTRI mouse) is REQUIRED — it gates 8f view 5 and 8g view 3.
Default ON; --no-tf only for QA reruns.

Per-cell pathway scoring uses decoupler AUCell (rank-based; insensitive to
normalization; consistent with the rank-based GSEA above). Scope: HALLMARK
(always) + union of pathways sig at FDR < le_fdr in the pseudobulk GSEA.
Bounded by what 8b/8c actually flagged — not the whole MSigDB universe.

Smoke-test mode (--smoke-test) filters the work list to ONE slice. Hard-fails
if the requested slice doesn't exist (no silent auto-pick).

Usage:
  uv run python scripts/08c_pathways.py --config config/brain.yaml --n-jobs 16
  uv run python scripts/08c_pathways.py --config config/placenta.yaml --n-jobs 16
  uv run python scripts/08c_pathways.py --config config/brain.yaml \
      --subcluster immune --n-jobs 16
  # smoke test:
  uv run python scripts/08c_pathways.py --config config/brain.yaml \
      --smoke-test --smoke-celltype Immune \
      --smoke-contrast early_vs_relaxed_per_age \
      --smoke-sex combined --smoke-group-level 4W --smoke-level whole \
      --no-tf --no-per-cell --n-jobs 1
"""

import argparse
import ast
import re
import sys
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from _utils import (load_config, phase_table_dir, parallel_map, unassigned_mask)


# ---------------------------------------------------------------------------
# Constants — single source of truth
# ---------------------------------------------------------------------------

CONTAM_PREFIX = "Contamination"

# Per-tissue obs columns used to identify unassigned cells (matches 8b).
TISSUE_UA_KEYS = {
    "brain": ["celltypist_broad", "celltypist_class"],
    "placenta": ["celltype_majority"],
}
# Focal subclusters per tissue (for the contam join in main mode).
TISSUE_FOCAL = {
    "brain":    ["Immune", "OPC/Oligodendrocytes", "Astrocytes/Ependymal"],
    "placenta": ["DSC", "Endothelium", "Myeloid", "NK"],
}
# Region key per tissue (None == whole only).
TISSUE_REGION_KEY = {"brain": "celltypist_region", "placenta": None}
# Per-tissue cell-type label column (for per-cell aggregation).
TISSUE_CELLTYPE_COL = {"brain": "celltypist_broad", "placenta": "celltype_majority"}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _safe(name) -> str:
    """Filesystem-safe slug."""
    return re.sub(r"[^0-9A-Za-z._-]+", "_", str(name)).strip("_")


def slugify(name) -> str:
    """Match 8b's slugify (lowercase, _-separated)."""
    return re.sub(r"[^0-9a-zA-Z]+", "_", str(name)).strip("_").lower()


def parse_pair(pair_repr) -> list:
    """8b stores `pair` as str(list), e.g. "['Early_Stress', 'Relaxed']".
    Parse safely back to a list of strings. Returns [] on failure / NaN."""
    if pair_repr is None or (isinstance(pair_repr, float) and np.isnan(pair_repr)):
        return []
    s = str(pair_repr).strip()
    if not s or s.lower() == "nan":
        return []
    try:
        v = ast.literal_eval(s)
        return [str(x) for x in v] if isinstance(v, (list, tuple)) else []
    except (ValueError, SyntaxError):
        return []


def pair_slug(pair_repr) -> str:
    """Turn pair into a folder/file-safe slug. '['Early_Stress','Relaxed']'
    -> 'Early_Stress_vs_Relaxed'. Empty pair -> 'none'."""
    p = parse_pair(pair_repr)
    if not p:
        return "none"
    return _safe("_vs_".join(p))


def is_contam(name) -> bool:
    s = str(name)
    return s.startswith(CONTAM_PREFIX) or s == "unresolved"


def get_symbol_map(adata):
    """Find a symbol column in adata.var and return {var_name: symbol} + col."""
    for col in ("symbol", "gene_symbol", "gene_symbols", "Symbol"):
        if col in adata.var.columns:
            return dict(zip(adata.var_names.astype(str), adata.var[col].astype(str))), col
    return None, None


def load_genesets_tsv(path: Path, collections: list, min_genes: int) -> pd.DataFrame:
    """Load the msigdbr TSV (collection / subcollection / gs_name / gene_symbol)
    into a decoupler net DataFrame [source, target, collection]. Filters to the
    requested collections and drops tiny sets."""
    df = pd.read_csv(path, sep="\t")
    needed = {"collection", "gs_name", "gene_symbol"}
    if not needed.issubset(df.columns):
        sys.exit(f"ERROR: {path} missing columns {needed - set(df.columns)}. "
                 f"Regenerate with fetch_genesets.R.")
    if collections:
        df = df[df["collection"].isin(collections)]
    net = (df.rename(columns={"gs_name": "source", "gene_symbol": "target"})
             [["source", "target", "collection"]].drop_duplicates())
    sizes = net.groupby("source").size()
    net = net[net["source"].isin(sizes[sizes >= min_genes].index)]
    return net.reset_index(drop=True)


# ---------------------------------------------------------------------------
# GSEA + TF (version-robust across decoupler 2.0 and 1.9)
# ---------------------------------------------------------------------------

def load_collectri():
    """Fetch CollecTRI mouse TF–target network."""
    import decoupler as dc
    if hasattr(dc, "op") and hasattr(dc.op, "collectri"):
        return dc.op.collectri(organism="mouse")
    return dc.get_collectri(organism="mouse", split_complexes=False)


def run_gsea_on_ranks(rank_series, net, min_genes, times=None, seed=42,
                       max_genes=500):
    """One GSEA call on a single ranking vector. Returns [source, NES, pvalue].
    FDR added separately by add_fdr() (per-collection BH + pooled).

    Implementation: fgsea-multilevel via R subprocess (Korotkevich et al. 2021,
    bioRxiv 060012). Replaces the prior decoupler.mt.gsea call — decoupler uses
    basic permutation with a 1/nperm p-value floor that produced bimodal FDRs
    (everything either FDR=0 or FDR>>0.05). fgsea-multilevel estimates
    continuous p-values down to ~1e-50 via adaptive multilevel splitting.

    Args:
      rank_series : pd.Series indexed by gene symbol; values = continuous rank
                    statistic (Wald stat from 8b).
      net         : long-form pathway DataFrame with columns `source` (pathway
                    name) and `target` (gene symbol). Same `net` previously
                    passed to decoupler. May contain multiple collections —
                    fgsea is called once on the union; per-collection FDR is
                    applied downstream by add_fdr() using collection_map.
      min_genes   : pathway size lower bound (passed to fgsea minSize).
      times       : ignored. Accepted for back-compat with the previous
                    decoupler-based signature; fgsea-multilevel adapts
                    internally and does not use a fixed permutation count.
      seed        : RNG seed for fgsea's internal multilevel splitting.
      max_genes   : pathway size upper bound (default 500; fgsea convention).

    Returns:
      pd.DataFrame with columns [source, NES, pvalue]. Same schema as before.
    """
    import subprocess
    import tempfile

    if rank_series is None or rank_series.empty:
        return pd.DataFrame(columns=["source", "NES", "pvalue"])

    worker = Path(__file__).resolve().parent / "run_fgsea.R"
    if not worker.is_file():
        raise FileNotFoundError(
            f"fgsea worker not found at {worker}. "
            f"Expected sibling of 08c_pathways.py.")

    with tempfile.TemporaryDirectory(prefix="fgsea_") as td:
        td = Path(td)
        ranks_path    = td / "ranks.tsv"
        pathways_path = td / "pathways.tsv"
        out_path      = td / "out.tsv"

        # ranks: gene\tstat
        pd.DataFrame({
            "gene": rank_series.index.astype(str),
            "stat": pd.to_numeric(rank_series.values, errors="coerce"),
        }).to_csv(ranks_path, sep="\t", index=False)

        # pathways: pathway_name\tgene  (from net.source/net.target)
        if not {"source", "target"}.issubset(net.columns):
            raise ValueError("net must have columns 'source' and 'target'")
        (net[["source", "target"]]
            .rename(columns={"source": "pathway_name", "target": "gene"})
            .to_csv(pathways_path, sep="\t", index=False))

        cmd = ["Rscript", str(worker),
               str(ranks_path), str(pathways_path), str(out_path),
               str(int(min_genes)), str(int(max_genes)), str(int(seed))]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
        if proc.returncode != 0:
            raise RuntimeError(
                f"fgsea worker failed (returncode={proc.returncode})\n"
                f"  cmd: {' '.join(cmd)}\n"
                f"  stderr (tail): {proc.stderr[-1500:]}\n"
                f"  stdout (tail): {proc.stdout[-500:]}")
        if (not out_path.exists()) or out_path.stat().st_size == 0:
            return pd.DataFrame(columns=["source", "NES", "pvalue"])

        out = pd.read_csv(out_path, sep="\t")

    if out.empty:
        return pd.DataFrame(columns=["source", "NES", "pvalue"])

    # fgsea schema -> our schema. add_fdr() recomputes FDR per collection.
    # fgsea's extra columns (ES, padj, size, leading_edge, log2err) are
    # intentionally dropped here — leading edges are computed downstream in
    # _gsea_worker by compute_leading_edge() to stay consistent with the
    # existing leading-edge CSV layout.
    return pd.DataFrame({
        "source": out["pathway"].astype(str).values,
        "NES":    pd.to_numeric(out["NES"], errors="coerce").values,
        "pvalue": pd.to_numeric(out["pval"], errors="coerce").values,
    })


def run_tf_ulm(rank_series, collectri, min_targets=5):
    """TF activity via ULM on the ranked DE-stat vector.
    Returns [source, activity_score, pvalue]. Version-robust.

    NOTE on p-values: decoupler 2.0's mt.ulm returns a permutation-based p-value
    that defaults to very few permutations -> badly underpowered (we observed
    activity_score |t|≈3.7 yielding p≈0.06 and ties across distinct TFs). The
    activity_score IS a regression t-statistic, so we OVERRIDE decoupler's
    p-value with the analytical t-distribution p-value (df = n_genes - 2).
    This matches what every TF-activity paper using ULM does (Holland 2020,
    Müller-Dott 2023, decoupler paper). Decoupler's raw permutation p is kept
    as `pvalue_dc` for reference / sensitivity QA.
    """
    import decoupler as dc
    from scipy.stats import t as student_t
    mat = rank_series.to_frame().T
    mat.index = ["contrast"]
    n_genes = len(rank_series)
    df_t = max(n_genes - 2, 1)

    if hasattr(dc, "mt") and hasattr(dc.mt, "ulm"):
        out = dc.mt.ulm(data=mat, net=collectri, tmin=min_targets)
        if isinstance(out, tuple):
            est, pval = out[0], out[1]
            res = pd.DataFrame({"source": est.columns,
                                "activity_score": est.iloc[0].values,
                                "pvalue_dc": pval.iloc[0].values})
        else:
            # AnnData-style fallback
            try:
                est = dc.pp.get_obsm(out, key="score_ulm")
                pv = dc.pp.get_obsm(out, key="padj_ulm")
                res = pd.DataFrame({"source": est.var_names,
                                    "activity_score": est.X[0],
                                    "pvalue_dc": pv.X[0] if pv is not None else np.nan})
            except Exception:
                res = pd.DataFrame(columns=["source", "activity_score", "pvalue_dc"])
    else:
        # decoupler 1.9 fallback
        acts, pvals = dc.run_ulm(mat=mat, net=collectri, min_n=min_targets)
        res = pd.DataFrame({"source": acts.columns,
                            "activity_score": acts.iloc[0].values,
                            "pvalue_dc": pvals.iloc[0].values})

    # Analytical 2-sided p-value from the t-statistic — overrides decoupler's
    if not res.empty:
        t_abs = np.abs(res["activity_score"].astype(float).values)
        res["pvalue"] = 2.0 * student_t.sf(t_abs, df=df_t)
    else:
        res["pvalue"] = np.nan
    return res[["source", "activity_score", "pvalue", "pvalue_dc"]]


def add_fdr(gsea_df, collection_map):
    """Attach collection + BH-correct two ways: within each collection (FDR)
    and pooled across all sets (FDR_pooled)."""
    from scipy.stats import false_discovery_control
    g = gsea_df.copy()
    g["collection"] = g["source"].map(collection_map).fillna("NA")
    p = g["pvalue"].fillna(1.0).values
    if np.isfinite(p).any():
        g["FDR_pooled"] = false_discovery_control(p)
    else:
        g["FDR_pooled"] = np.nan
    g["FDR"] = np.nan
    for coll, idx in g.groupby("collection").groups.items():
        sub = g.loc[idx, "pvalue"].fillna(1.0).values
        if len(sub):
            g.loc[idx, "FDR"] = false_discovery_control(sub)
    return g


def bh(p_arr):
    """BH FDR with NaN-safety. Returns array same length as p_arr."""
    from scipy.stats import false_discovery_control
    p = np.asarray(p_arr, dtype=float)
    ok = ~np.isnan(p)
    out = np.full(len(p), np.nan)
    if ok.sum() > 0:
        out[ok] = false_discovery_control(p[ok], method="bh")
    return out


def compute_leading_edge(rank_series, members, nes):
    """Members driving the enrichment, sorted by |stat| descending. Sign of
    NES determines whether we pick stat>0 or stat<0 members."""
    present = [g for g in members if g in rank_series.index]
    if not present:
        return []
    sub = rank_series.loc[present]
    if nes is not None and not np.isnan(nes):
        sub = sub[sub > 0] if nes > 0 else sub[sub < 0]
    sub = sub.reindex(sub.abs().sort_values(ascending=False).index)
    return list(zip(sub.index.tolist(), sub.values.tolist()))


# ---------------------------------------------------------------------------
# Worker — runs in parallel via parallel_map
# ---------------------------------------------------------------------------

def _gsea_worker(job, net=None, collectri=None, min_genes=5, times=1000,
                 le_fdr=0.05, stat_col="stat", run_tf=True):
    """Process one slice. `job` is a dict with: sex, contrast, flag,
    group_level, pair, level, celltype, n_donors_total, reliability, note,
    and `de_df` (the gene-level rows pre-ranked into a Series under 'rank').

    Returns dict with three lists: gsea_rows, le_rows, tf_rows.
    """
    rank_series = job["rank_series"]
    lfc_map = job["lfc_map"]
    collection_map = job["collection_map"]
    set_members = job["set_members"]   # dict: pathway -> set(target)

    base = {
        "tissue": job["tissue"],
        "sex": job["sex"], "contrast": job["contrast"], "flag": job["flag"],
        "group_level": job["group_level"], "pair": job["pair"],
        "level": job["level"], "celltype": job["celltype"],
        "n_donors_total": job["n_donors_total"],
        "reliability": job["reliability"], "note": job["note"],
    }

    gsea_rows, le_rows, tf_rows = [], [], []

    # GSEA
    try:
        gsea = run_gsea_on_ranks(rank_series, net, min_genes, times)
    except Exception as e:
        return {"gsea_rows": [], "le_rows": [], "tf_rows": [],
                "error": f"GSEA failed: {e}"}
    if gsea is None or gsea.empty:
        return {"gsea_rows": [], "le_rows": [], "tf_rows": [], "error": None}

    gsea = add_fdr(gsea, collection_map)

    for _, g in gsea.iterrows():
        gsea_rows.append({
            **base, "collection": g["collection"], "source": g["source"],
            "NES": g["NES"], "pvalue": g["pvalue"],
            "FDR": g["FDR"], "FDR_pooled": g["FDR_pooled"],
        })

    # Leading-edge for pathways at FDR<le_fdr
    sig = gsea.dropna(subset=["FDR"])
    sig = sig[sig["FDR"] < le_fdr]
    for _, g in sig.iterrows():
        members = set_members.get(g["source"], set())
        le = compute_leading_edge(rank_series, members, g["NES"])
        for rank_i, (gene, stat_val) in enumerate(le, start=1):
            le_rows.append({
                **base, "collection": g["collection"], "pathway": g["source"],
                "NES": round(float(g["NES"]), 4),
                "pathway_FDR": round(float(g["FDR"]), 6),
                "leading_edge_rank": rank_i, "gene": gene,
                "log2FC": (round(float(lfc_map[gene]), 4)
                           if gene in lfc_map and not np.isnan(lfc_map[gene])
                           else np.nan),
                "rank_stat": round(float(stat_val), 4),
                "direction": "up" if stat_val > 0 else "down",
            })

    # TF activity (CollecTRI ULM); BH within this slice
    if run_tf and collectri is not None:
        try:
            tf = run_tf_ulm(rank_series, collectri)
        except Exception as e:
            return {"gsea_rows": gsea_rows, "le_rows": le_rows, "tf_rows": [],
                    "error": f"TF ULM failed: {e}"}
        if tf is not None and not tf.empty:
            tf["FDR"] = bh(tf["pvalue"].values)
            for _, t in tf.iterrows():
                tf_rows.append({
                    **base, "TF": t["source"],
                    "activity_score": round(float(t["activity_score"]), 4),
                    "pvalue": float(t["pvalue"]),
                    "pvalue_dc": (float(t["pvalue_dc"])
                                  if "pvalue_dc" in t and pd.notna(t["pvalue_dc"])
                                  else np.nan),
                    "FDR": (round(float(t["FDR"]), 6)
                            if not np.isnan(t["FDR"]) else np.nan),
                    "direction": ("activated" if t["activity_score"] > 0
                                  else "repressed"),
                })

    return {"gsea_rows": gsea_rows, "le_rows": le_rows, "tf_rows": tf_rows,
            "error": None}


# ---------------------------------------------------------------------------
# Build work list
# ---------------------------------------------------------------------------

GROUP_KEYS = ["sex", "contrast", "flag", "group_level", "pair", "level", "celltype"]


def build_work_jobs(de_df, net, collection_map, symbol_map, stat_col,
                    tissue, min_genes_per_slice=10):
    """Return a list of work dicts. One per Wald slice with >=min_genes rows.

    Wald-only (LRT has no Wald stat). Each job carries its own ranking Series
    and lfc map keyed by gene symbol (already resolved from Ensembl if needed).
    `set_members` is the shared per-pathway membership dict — passed by reference
    so we don't rebuild it for every job.
    """
    wald = de_df[de_df["test_method"] == "Wald"].copy()
    if wald.empty:
        return []
    # Drop rows without a gene or stat
    wald = wald.dropna(subset=["gene", stat_col])
    if wald.empty:
        return []

    # gene_symbol resolution (Ensembl -> symbol if a map is present)
    if symbol_map:
        wald["gene_sym"] = wald["gene"].astype(str).map(lambda g: symbol_map.get(g, g))
    else:
        wald["gene_sym"] = wald["gene"].astype(str)

    set_members = (net.groupby("source")["target"]
                      .apply(lambda s: set(s.values)).to_dict())

    jobs = []
    for keys, sub in wald.groupby(GROUP_KEYS, observed=True, dropna=False):
        if sub.shape[0] < min_genes_per_slice:
            continue
        sex, contrast, flag, group_level, pair, level, celltype = keys
        # Drop gene-symbol duplicates by keeping the larger |stat| (e.g. when
        # multiple Ensembl IDs map to the same symbol — rare but happens).
        s = (sub.reindex(sub[stat_col].abs().sort_values(ascending=False).index)
                .drop_duplicates("gene_sym"))
        rank_series = s.set_index("gene_sym")[stat_col].astype(float)
        lfc_map = (s.set_index("gene_sym")["log2FC"].astype(float).to_dict()
                   if "log2FC" in s.columns else {})
        n_total = (int(s["n_donors_total"].iloc[0])
                   if "n_donors_total" in s.columns and pd.notna(s["n_donors_total"].iloc[0])
                   else -1)
        reliability = (str(s["reliability"].iloc[0])
                       if "reliability" in s.columns else "")
        note = str(s["note"].iloc[0]) if "note" in s.columns else ""
        jobs.append({
            "tissue": tissue,
            "sex": sex, "contrast": contrast, "flag": flag,
            "group_level": group_level, "pair": pair, "level": level,
            "celltype": celltype,
            "n_donors_total": n_total, "reliability": reliability, "note": note,
            "rank_series": rank_series, "lfc_map": lfc_map,
            "collection_map": collection_map, "set_members": set_members,
        })
    return jobs


def apply_smoke_filter(jobs, args):
    """Filter the job list down to the matching slice for --smoke-test. Hard-
    fails if no match (no silent auto-pick).

    --smoke-pair accepts EITHER the repr ("['Early_Stress', 'Relaxed']") OR
    the slug ("Early_Stress_vs_Relaxed") — the latter is easier on a CLI.
    """
    direct = {
        "celltype": args.smoke_celltype, "contrast": args.smoke_contrast,
        "sex": args.smoke_sex, "group_level": args.smoke_group_level,
        "level": args.smoke_level,
    }
    pair_filter = args.smoke_pair

    def _match(j):
        for k, v in direct.items():
            if v is not None and str(j.get(k)) != str(v):
                return False
        if pair_filter is not None:
            slug = pair_slug(j.get("pair"))
            if str(pair_filter) != slug and str(pair_filter) != str(j.get("pair")):
                return False
        return True

    filtered = [j for j in jobs if _match(j)]
    if not filtered:
        spec = ", ".join(f"{k}={v}" for k, v in direct.items() if v is not None)
        if pair_filter is not None:
            spec += f", pair={pair_filter}"
        sys.exit(f"ERROR: --smoke-test matched 0 slices ({spec}). "
                 f"Check that 8b emitted that exact slice.")
    return filtered


# ---------------------------------------------------------------------------
# Per-cell AUCell + per-donor aggregation
# ---------------------------------------------------------------------------

def _collect_contam_barcodes(rdir: Path, focal_list: list) -> set:
    """In main mode the integrated h5ad has no subcluster_name. Walk the focal
    subcluster h5ads under h5ad/08c_subclustered/ to pull contam barcodes."""
    subc_dir = rdir / "h5ad" / "08c_subclustered"
    bc = set()
    for focal in focal_list:
        slug = slugify(focal)
        path = subc_dir / f"{slug}.h5ad"
        if not path.is_file():
            print(f"  [warn] subcluster h5ad not found: {path.name} — skipping contam join")
            continue
        a = ad.read_h5ad(path, backed="r")
        if "subcluster_name" in a.obs.columns:
            mask = a.obs["subcluster_name"].apply(is_contam)
            bc |= set(a.obs_names[mask].astype(str))
        a.file.close()
    return bc


def run_per_cell_aucell(adata, pathway_scope: list, net: pd.DataFrame,
                        min_genes_per_set=5):
    """Run AUCell. Modifies adata in place (adds obsm key). Returns the
    cells × pathways score DataFrame [index=cell, columns=pathway]."""
    import decoupler as dc
    net_sub = net[net["source"].isin(pathway_scope)].copy()
    sizes = net_sub.groupby("source").size()
    net_sub = net_sub[net_sub["source"].isin(sizes[sizes >= min_genes_per_set].index)]
    if net_sub.empty:
        raise RuntimeError("AUCell: no pathways in scope have enough genes.")

    # decoupler 2.0 (mt.aucell) — modifies AnnData in place
    if hasattr(dc, "mt") and hasattr(dc.mt, "aucell"):
        dc.mt.aucell(data=adata, net=net_sub)
        # decoupler stores scores in obsm under a versioned key — find it
        key = next((k for k in adata.obsm if "aucell" in str(k).lower()), None)
        if key is None:
            raise RuntimeError("AUCell ran but no obsm key with 'aucell' found.")
        scores = adata.obsm[key]
        if isinstance(scores, pd.DataFrame):
            return scores
        # If it's an AnnData (newer decoupler), extract X + var_names
        if hasattr(scores, "X") and hasattr(scores, "var_names"):
            arr = scores.X.toarray() if sp.issparse(scores.X) else np.asarray(scores.X)
            return pd.DataFrame(arr, index=adata.obs_names, columns=scores.var_names)
        # Otherwise assume ndarray
        return pd.DataFrame(np.asarray(scores), index=adata.obs_names,
                            columns=sorted(net_sub["source"].unique()))

    # decoupler 1.9
    dc.run_aucell(mat=adata, net=net_sub, source="source", target="target",
                  use_raw=False)
    key = next((k for k in adata.obsm if "aucell" in str(k).lower()), None)
    if key is None:
        raise RuntimeError("AUCell (1.9) ran but no obsm key with 'aucell' found.")
    return adata.obsm[key]


def aggregate_per_donor(scores_df: pd.DataFrame, obs: pd.DataFrame,
                        tissue: str, celltype_col: str,
                        region_key: str | None) -> pd.DataFrame:
    """Long-form per-donor pathway means. For each (donor_id, celltype,
    level, pathway): mean_score, median_score, n_cells.

    level = 'whole' aggregates all the donor's cells of that celltype;
    level = <region_name> aggregates the donor's cells in that region.
    """
    meta_cols = ["donor_id", "sample_id", "sex", "age", "group", "pool"]
    meta_cols = [c for c in meta_cols if c in obs.columns]
    if not meta_cols or celltype_col not in obs.columns:
        return pd.DataFrame()

    # Align scores and obs to common barcodes
    common = obs.index.intersection(scores_df.index)
    obs = obs.loc[common]
    scores = scores_df.loc[common].copy()
    scores.columns.name = "pathway"       # named level so stack() produces it cleanly

    # Build the grouping frame: meta_cols + celltype, indexed by cell barcode
    gframe = obs[meta_cols + [celltype_col]].copy()
    gframe = gframe.rename(columns={celltype_col: "celltype"})
    gb_cols = meta_cols + ["celltype"]

    # Region lookup (None for placenta)
    if region_key and region_key in obs.columns:
        region_series = obs[region_key].astype(str)
    else:
        region_series = None

    def _agg_one_level(level_label: str, mask) -> pd.DataFrame:
        if not mask.any():
            return None
        g_sub = gframe.loc[mask]
        s_sub = scores.loc[mask]
        # Group on real columns of a combined frame for clarity (no list-of-Series trick)
        combined = pd.concat([g_sub.reset_index(drop=True),
                              s_sub.reset_index(drop=True)], axis=1)
        grp = combined.groupby(gb_cols, observed=True, sort=False)
        means = grp[s_sub.columns].mean()
        medians = grp[s_sub.columns].median()
        sizes = grp.size().rename("n_cells")
        # The pathway-level name is lost across pd.concat -> grp -> mean — re-set
        # explicitly so .stack() yields a named "pathway" level (otherwise we'd
        # get a default "level_N" and the downstream merge fails).
        means.columns.name = "pathway"
        medians.columns.name = "pathway"
        mean_long = means.stack().rename("mean_score").reset_index()
        med_long = medians.stack().rename("median_score").reset_index()
        merged = mean_long.merge(med_long, on=gb_cols + ["pathway"])
        merged = merged.merge(sizes.reset_index(), on=gb_cols)
        merged["tissue"] = tissue
        merged["level"] = level_label
        merged["region"] = level_label if level_label != "whole" else "all"
        return merged

    pieces = []
    whole = _agg_one_level("whole", np.ones(len(obs), dtype=bool))
    if whole is not None:
        pieces.append(whole)
    if region_series is not None:
        for region in sorted(region_series.dropna().unique()):
            m = (region_series == region).values
            piece = _agg_one_level(region, m)
            if piece is not None:
                pieces.append(piece)

    if not pieces:
        return pd.DataFrame()
    out = pd.concat(pieces, ignore_index=True)
    head = ["tissue", "donor_id", "sample_id", "sex", "age", "group", "pool",
            "celltype", "level", "region", "pathway",
            "mean_score", "median_score", "n_cells"]
    head = [c for c in head if c in out.columns]
    out = out[head + [c for c in out.columns if c not in head]]
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Phase 8c: pathway/TF/per-cell on 8b DE")
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--n-jobs", type=int, default=16,
                    help="Parallel workers (threads). 16 is fine on the WS.")
    ap.add_argument("--stat-col", default="stat",
                    help="Column to rank genes by (default 'stat' = DESeq2 Wald).")
    ap.add_argument("--times", type=int, default=1000,
                    help="DEPRECATED: GSEA permutations per slice. "
                         "Ignored since switch to fgsea-multilevel which "
                         "adapts internally. Kept for back-compat with "
                         "existing runner scripts.")
    ap.add_argument("--le-fdr", type=float, default=0.05,
                    help="FDR cutoff for leading-edge rows AND per-cell pathway scope.")
    ap.add_argument("--tf", dest="tf", action="store_true", default=True,
                    help="TF activity ON (default).")
    ap.add_argument("--no-tf", dest="tf", action="store_false",
                    help="TF activity OFF (QA reruns only).")
    ap.add_argument("--per-cell", dest="per_cell", action="store_true", default=True,
                    help="Per-cell AUCell ON (default).")
    ap.add_argument("--no-per-cell", dest="per_cell", action="store_false",
                    help="Skip per-cell AUCell (saves ~10-30 min for QA reruns).")
    ap.add_argument("--per-cell-only", action="store_true",
                    help="Skip GSEA+TF; load existing master CSV from disk and "
                         "run per-cell AUCell only. Use to resume after a crash "
                         "in the AUCell step.")
    ap.add_argument("--per-cell-cap", type=int, default=1000,
                    help="Max pathways to score per cell (Hallmark + top-N "
                         "non-MH by sig-frequency across slices). 1000 keeps "
                         "the per-cell matrix at ~2.7 GB for 670k cells.")
    ap.add_argument("--subcluster", default=None,
                    help="Run on the 7b subcluster DE table for {slug} "
                         "(reads 08b_de_results_subcluster_{slug}.csv + "
                         "08c_subclustered/{slug}.h5ad).")
    # Smoke-test args (filter to one slice; no auto-select)
    ap.add_argument("--smoke-test", action="store_true",
                    help="Filter work list to ONE slice via --smoke-* args.")
    ap.add_argument("--smoke-celltype", default=None)
    ap.add_argument("--smoke-contrast", default=None)
    ap.add_argument("--smoke-sex", default=None)
    ap.add_argument("--smoke-group-level", default=None)
    ap.add_argument("--smoke-level", default=None)
    ap.add_argument("--smoke-pair", default=None,
                    help="Pair string repr, e.g. \"['Early_Stress', 'Relaxed']\".")
    args = ap.parse_args()

    print(f"\n=== Phase 8c: pathway / TF / per-cell ===")
    cfg = load_config(args.config)
    tissue = cfg.get("tissue")
    if tissue not in TISSUE_UA_KEYS:
        sys.exit(f"ERROR: unknown tissue '{tissue}' (expected {list(TISSUE_UA_KEYS)}).")

    rdir = Path(cfg["results_dir"])
    pcfg = cfg.get("pathways", {})
    geneset_tsv = Path(pcfg.get("geneset_tsv", "refs/msigdb_mouse.tsv"))
    collections = pcfg.get("collections", ["MH", "M2", "M5", "M8"])
    min_genes = int(pcfg.get("min_genes_per_set", 5))

    # ---- Inputs -------------------------------------------------------------
    suffix = f"_subcluster_{args.subcluster}" if args.subcluster else ""
    de_dir = rdir / "tables" / "08b_de"
    de_path = de_dir / f"08b_de_results{suffix}.csv"
    if not de_path.is_file():
        sys.exit(f"ERROR: {de_path} not found. Run 08b_de.py"
                 + (f" --subcluster {args.subcluster}" if args.subcluster else ""))

    # Output table dir is needed for both branches
    out_table_dir = phase_table_dir(cfg, f"08c_pathways{suffix}")
    out_table_dir.mkdir(parents=True, exist_ok=True)

    # Gene-set network is needed for both branches (per-cell uses it for AUCell)
    geneset_tsv = Path(pcfg.get("geneset_tsv", "refs/msigdb_mouse.tsv"))
    collections = pcfg.get("collections", ["MH", "M2", "M5", "M8"])
    min_genes = int(pcfg.get("min_genes_per_set", 5))
    if not geneset_tsv.is_file():
        sys.exit(f"ERROR: gene-set TSV not found: {geneset_tsv}\n"
                 f"  Generate once: Rscript scripts/fetch_genesets.R --out {geneset_tsv}")
    net = load_genesets_tsv(geneset_tsv, collections, min_genes)
    print(f"  Gene sets: {net['source'].nunique()} from {geneset_tsv} "
          f"(collections: {', '.join(collections)})")
    collection_map = dict(net.drop_duplicates("source")
                             .set_index("source")["collection"])

    # Integrated h5ad for symbol_map + per-cell
    if args.subcluster:
        h5_path = rdir / "h5ad" / "08c_subclustered" / f"{args.subcluster}.h5ad"
    else:
        h5_path = rdir / "h5ad" / "08_annotated" / "all_samples.h5ad"

    # ---- --per-cell-only RESUME BRANCH -------------------------------------
    if args.per_cell_only:
        gsea_csv = out_table_dir / f"08c_pathway_results{suffix}.csv"
        if not gsea_csv.is_file():
            sys.exit(f"ERROR: --per-cell-only requires {gsea_csv}; not found.\n"
                     f"  Run without --per-cell-only first to generate it.")
        print(f"\n  RESUME (--per-cell-only): loading {gsea_csv}")
        gsea_df = pd.read_csv(gsea_csv, low_memory=False)
        print(f"  Master GSEA: {len(gsea_df):,} rows loaded.")
        run_per_cell(rdir, cfg, args, tissue, gsea_df, net,
                     out_table_dir, suffix, h5_path)
        print(f"\n✓ Phase 8c (per-cell-only) complete.")
        return

    print(f"  DE table: {de_path}  ({de_path.stat().st_size/1e6:.1f} MB)")
    de = pd.read_csv(de_path, low_memory=False)
    if de.empty:
        sys.exit(f"ERROR: {de_path.name} has no rows.")

    symbol_map = None
    sym_col = None
    if h5_path.is_file():
        ann = ad.read_h5ad(h5_path, backed="r")
        symbol_map, sym_col = get_symbol_map(ann)
        ann.file.close()
        if symbol_map:
            print(f"  Symbol map from var['{sym_col}']: {len(symbol_map)} entries")
    else:
        print(f"  [warn] h5ad not found: {h5_path}  (per-cell AUCell will be skipped)")

    # Overlap sanity (fail loudly on Ensembl/symbol mismatch)
    de_genes = de["gene"].dropna().astype(str).unique()
    mapped = ([symbol_map.get(g, g) for g in de_genes] if symbol_map
              else list(de_genes))
    overlap = len(set(mapped) & set(net["target"]))
    print(f"  Gene overlap (DE symbols ∩ gene sets): {overlap}")
    if overlap < 5:
        sys.exit(f"ERROR: only {overlap} DE genes overlap gene sets.\n"
                 f"  Likely Ensembl-vs-symbol mismatch or wrong organism GMT.\n"
                 f"  DE examples: {list(de_genes[:5])}\n"
                 f"  Set targets: {list(net['target'].unique()[:5])}")

    # CollecTRI (TF activity)
    collectri = None
    if args.tf:
        try:
            collectri = load_collectri()
            n_tf = (collectri["source"].nunique()
                    if "source" in collectri.columns
                    else collectri.iloc[:, 0].nunique())
            print(f"  CollecTRI mouse: {n_tf} TFs (ULM per slice)")
        except Exception as e:
            print(f"  [warn] CollecTRI fetch failed ({e}); TF activity SKIPPED.")
            collectri = None

    # ---- Work list ---------------------------------------------------------
    n_total_rows = len(de)
    n_lrt = int((de["test_method"] == "LRT").sum()) if "test_method" in de.columns else 0
    print(f"  DE rows: {n_total_rows:,}  (LRT rows skipped: {n_lrt:,})")
    jobs = build_work_jobs(de, net, collection_map, symbol_map,
                           args.stat_col, tissue)
    print(f"  Work jobs (Wald slices, n_genes>=10): {len(jobs):,}")
    if args.smoke_test:
        jobs = apply_smoke_filter(jobs, args)
        print(f"  SMOKE-TEST mode: filtered to {len(jobs)} job(s)")
    if not jobs:
        sys.exit("ERROR: no jobs to run.")

    # ---- Parallel run ------------------------------------------------------
    def _fn(job):
        return _gsea_worker(
            job, net=net, collectri=collectri, min_genes=min_genes,
            times=args.times, le_fdr=args.le_fdr,
            stat_col=args.stat_col, run_tf=(args.tf and collectri is not None))

    gsea_rows, le_rows, tf_rows = [], [], []
    n_err = 0
    print(f"\n  Running parallel_map ({args.n_jobs} threads)...")
    for job, res, err in parallel_map(_fn, jobs, n_jobs=args.n_jobs,
                                       use_threads=True,
                                       desc="08c GSEA+TF"):
        if err:
            n_err += 1
            print(f"  [err] {job.get('contrast')}|{job.get('sex')}|"
                  f"{job.get('group_level')}|{job.get('level')}|"
                  f"{job.get('celltype')}: {str(err)[:160]}")
            continue
        if res is None:
            continue
        if res.get("error"):
            n_err += 1
            print(f"  [err] {job.get('celltype')}: {res['error'][:160]}")
        gsea_rows.extend(res["gsea_rows"])
        le_rows.extend(res["le_rows"])
        tf_rows.extend(res["tf_rows"])
    print(f"  Completed: {len(jobs)-n_err}/{len(jobs)} slices "
          f"({n_err} errors).")

    # ---- Write CSVs --------------------------------------------------------
    gsea_df = pd.DataFrame(gsea_rows)
    gsea_csv = out_table_dir / f"08c_pathway_results{suffix}.csv"
    gsea_df.to_csv(gsea_csv, index=False)
    print(f"\n  Master GSEA: {gsea_csv}  ({len(gsea_df):,} rows)")
    if not gsea_df.empty:
        for coll, sub in gsea_df.groupby("collection"):
            sub_path = out_table_dir / f"08c_pathway_results{suffix}_{_safe(str(coll))}.csv"
            sub.to_csv(sub_path, index=False)
        print(f"  Per-collection splits: {gsea_df['collection'].nunique()} files")

    le_df = pd.DataFrame(le_rows)
    le_csv = out_table_dir / f"08c_pathway_leading_edge{suffix}.csv"
    le_df.to_csv(le_csv, index=False)
    if le_rows:
        n_paths = le_df.groupby(["contrast", "group_level", "level",
                                  "celltype", "pathway"]).ngroups
        print(f"  Leading-edge: {le_csv}  ({len(le_df):,} gene rows, "
              f"{n_paths} sig pathways at FDR<{args.le_fdr})")
    else:
        print(f"  Leading-edge: {le_csv}  (empty — no pathways below FDR<{args.le_fdr})")

    if args.tf:
        tf_df = pd.DataFrame(tf_rows)
        # Per project doc: BH within (contrast × celltype) pools across
        # {sex, group_level, pair, level} for the same (contrast, celltype).
        # Stronger than within-slice when a TF is active in several
        # co-occurring slices. Within-slice FDR (`FDR`) is kept as the strict
        # column; this is the preferred sig call for downstream (8f / 8g).
        if not tf_df.empty:
            tf_df["FDR_ctx_celltype"] = np.nan
            for keys, idx in tf_df.groupby(["contrast", "celltype"],
                                            observed=True).groups.items():
                p = tf_df.loc[idx, "pvalue"].values
                tf_df.loc[idx, "FDR_ctx_celltype"] = bh(p)
        tf_csv = out_table_dir / f"08c_tf_activity{suffix}.csv"
        tf_df.to_csv(tf_csv, index=False)
        n_sig_tf = int((tf_df["FDR"] < 0.05).sum()) if not tf_df.empty else 0
        n_sig_tf_pooled = (int((tf_df["FDR_ctx_celltype"] < 0.05).sum())
                           if not tf_df.empty else 0)
        print(f"  TF activity: {tf_csv}  ({len(tf_df):,} rows, "
              f"{n_sig_tf} at FDR<0.05 within-slice, "
              f"{n_sig_tf_pooled} at FDR_ctx_celltype<0.05)")

    # ---- Per-cell AUCell ----------------------------------------------------
    if args.per_cell and not args.smoke_test:
        run_per_cell(rdir, cfg, args, tissue, gsea_df, net,
                     out_table_dir, suffix, h5_path)
    elif args.smoke_test:
        print(f"\n  Per-cell AUCell SKIPPED (--smoke-test mode).")
    elif not args.per_cell:
        print(f"\n  Per-cell AUCell SKIPPED (--no-per-cell).")

    print(f"\n✓ Phase 8c complete.")


def run_per_cell(rdir, cfg, args, tissue, gsea_df, net, out_table_dir,
                 suffix, h5_path):
    """Run AUCell on all cells, write per-cell h5ad + per-donor CSV.

    RESUME: if the per-cell h5ad already exists, load it and skip AUCell —
    only the per-donor aggregation re-runs. Delete the h5ad to force a
    fresh AUCell run.
    """
    print(f"\n=== Per-cell AUCell ===")

    # Resume path — if per-cell h5ad already exists, skip AUCell entirely.
    pc_dir = rdir / "h5ad" / "08c_pathway_scores"
    pc_dir.mkdir(parents=True, exist_ok=True)
    pc_path = pc_dir / f"{tissue}{suffix}_per_cell_scores.h5ad"
    if pc_path.is_file():
        print(f"  Per-cell h5ad exists -> SKIPPING AUCell.")
        print(f"    {pc_path}  ({pc_path.stat().st_size/1e6:.0f} MB)")
        print(f"    (delete the file to force a fresh AUCell run)")
        pc_adata = ad.read_h5ad(pc_path)
        X = pc_adata.X
        if sp.issparse(X):
            X = X.toarray()
        scores_df = pd.DataFrame(np.asarray(X), index=pc_adata.obs_names,
                                  columns=pc_adata.var_names)
        celltype_col = TISSUE_CELLTYPE_COL[tissue]
        region_key = TISSUE_REGION_KEY[tissue]
        if args.subcluster and "subcluster_name" in pc_adata.obs.columns:
            celltype_col = "subcluster_name"
        per_donor = aggregate_per_donor(scores_df, pc_adata.obs, tissue,
                                         celltype_col, region_key)
        pd_csv = out_table_dir / f"08c_pathway_scores_per_donor{suffix}.csv"
        per_donor.to_csv(pd_csv, index=False)
        print(f"  Per-donor pathway means: {pd_csv}  ({len(per_donor):,} rows)")
        return

    if not h5_path.is_file():
        print(f"  [skip] h5ad not found: {h5_path}")
        return

    # Pathway scope: Hallmark (always) + top-N non-Hallmark by sig-frequency
    # across slices (the most consistently-flagged pathways are the most
    # robust biology to score at single-cell level). Capping keeps the
    # cells × pathways matrix tractable: at 670k cells, 1000 pathways =
    # ~2.7 GB float32 — fine. 8000+ pathways = 21+ GB — unusable.
    hallmark_paths = sorted(set(net.loc[net["collection"] == "MH", "source"]))
    if not gsea_df.empty:
        sig_mask = gsea_df["FDR"].notna() & (gsea_df["FDR"] < args.le_fdr)
        sig_df = gsea_df.loc[sig_mask]
        # Rank non-Hallmark sig pathways by frequency of significance
        sig_counts = sig_df["source"].value_counts()
        non_mh_counts = sig_counts[~sig_counts.index.isin(hallmark_paths)]
        cap_non_mh = max(args.per_cell_cap - len(hallmark_paths), 0)
        top_non_mh = non_mh_counts.head(cap_non_mh).index.tolist()
        pathway_scope = sorted(set(hallmark_paths) | set(top_non_mh))
        print(f"  Pathway scope: {len(pathway_scope)} sets  "
              f"(MH: {len(hallmark_paths)} anchor + top {len(top_non_mh)} "
              f"non-MH by sig-frequency; cap={args.per_cell_cap}, "
              f"{len(non_mh_counts)} non-MH sig before cap)")
    else:
        pathway_scope = sorted(set(hallmark_paths))
        print(f"  Pathway scope: {len(pathway_scope)} sets (Hallmark only — no sig)")
    if not pathway_scope:
        print(f"  [skip] no pathways in scope.")
        return

    # Load h5ad (full, not backed — AUCell needs to score)
    print(f"  Loading {h5_path}...")
    adata = sc.read_h5ad(h5_path)
    n_start = adata.n_obs
    print(f"  Cells loaded: {n_start:,}")

    # Drop unassigned + contaminants (mirror 8b's drop logic)
    ua_keys = [k for k in TISSUE_UA_KEYS[tissue] if k in adata.obs.columns]
    if ua_keys:
        ua = unassigned_mask(adata.obs, ua_keys)
        adata = adata[~ua].copy()
        print(f"  Dropped unassigned: {int(ua.sum()):,}  -> {adata.n_obs:,} cells")

    # Contam: subcluster mode = directly from subcluster_name; main mode = join
    if args.subcluster:
        if "subcluster_name" in adata.obs.columns:
            contam = adata.obs["subcluster_name"].apply(is_contam)
            adata = adata[~contam].copy()
            print(f"  Dropped contam (subcluster_name): {int(contam.sum()):,}  "
                  f"-> {adata.n_obs:,} cells")
    else:
        contam_bc = _collect_contam_barcodes(rdir, TISSUE_FOCAL[tissue])
        if contam_bc:
            keep = ~adata.obs_names.astype(str).isin(contam_bc)
            n_dropped = int((~keep).sum())
            adata = adata[keep].copy()
            print(f"  Dropped contam (focal subcluster join): {n_dropped:,}  "
                  f"-> {adata.n_obs:,} cells")

    if adata.n_obs == 0:
        print(f"  [skip] all cells dropped.")
        return

    # Symbol resolution on var (AUCell needs to match net's symbol targets)
    if "symbol" in adata.var.columns:
        adata.var_names = adata.var["symbol"].astype(str).values
        # Drop duplicate symbols (keep first)
        keep_mask = ~adata.var_names.duplicated()
        if (~keep_mask).any():
            print(f"  Dropping {(~keep_mask).sum()} duplicate symbol var rows")
            adata = adata[:, keep_mask].copy()

    # Run AUCell
    print(f"  Running AUCell on {adata.n_obs:,} cells × {len(pathway_scope)} pathways...")
    try:
        scores_df = run_per_cell_aucell(adata, pathway_scope, net,
                                         min_genes_per_set=5)
    except Exception as e:
        print(f"  [err] AUCell failed: {e}")
        return
    print(f"  AUCell done. Scores shape: {scores_df.shape}")

    # Build the per-cell AnnData (cells × pathways)
    keep_obs = [c for c in ["donor_id", "sample_id", "sex", "age", "group",
                             "pool", "celltypist_broad", "celltype_majority",
                             "celltypist_class", "celltypist_subclass",
                             "celltypist_region", "subcluster_name",
                             "leiden", "n_genes", "n_counts"]
                if c in adata.obs.columns]
    obs_small = adata.obs[keep_obs].copy()
    var_small = pd.DataFrame({
        "collection": [net.loc[net["source"] == p, "collection"].iloc[0]
                       if (net["source"] == p).any() else "NA"
                       for p in scores_df.columns],
    }, index=scores_df.columns)
    pc_adata = ad.AnnData(
        X=scores_df.astype(np.float32).values,
        obs=obs_small.loc[scores_df.index],
        var=var_small,
    )
    if "X_umap" in adata.obsm:
        pc_adata.obsm["X_umap"] = adata.obsm["X_umap"][
            adata.obs_names.get_indexer(scores_df.index)]
    pc_adata.write_h5ad(pc_path)
    print(f"  Per-cell scores: {pc_path}  "
          f"({pc_adata.shape[0]:,} cells × {pc_adata.shape[1]} pathways, "
          f"{pc_path.stat().st_size/1e6:.0f} MB)")

    # Per-donor aggregation
    celltype_col = TISSUE_CELLTYPE_COL[tissue]
    region_key = TISSUE_REGION_KEY[tissue]
    # For subcluster mode, prefer subcluster_name as the celltype column
    if args.subcluster and "subcluster_name" in adata.obs.columns:
        celltype_col = "subcluster_name"
    per_donor = aggregate_per_donor(scores_df, adata.obs, tissue,
                                     celltype_col, region_key)
    pd_csv = out_table_dir / f"08c_pathway_scores_per_donor{suffix}.csv"
    per_donor.to_csv(pd_csv, index=False)
    print(f"  Per-donor pathway means: {pd_csv}  ({len(per_donor):,} rows)")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    main()
