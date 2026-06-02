#!/usr/bin/env python
"""
08c_pathways.py — Phase 8c: pathway / gene-set enrichment on the 8b DE results.

Runs GSEA on the ranked DE statistics (decoupler get_gsea_df), per contrast x
cell type, against gene sets. Optionally infers TF activity (CollecTRI).
Reads de_results.csv from 8b — does NOT re-run DE.

Ranking metric: the DESeq2 Wald 'stat' column. decoupler docs: contrast-level
stats (Wald / logFC) need no transformation before GSEA — ideal as the ranking.

GENE SETS — read from LOCAL .gmt files, not decoupler's get_resource('MSigDB',
organism='mouse'), which is broken (open issues: pypath ortholog-translation /
decompression errors). Download mouse-native GMTs once (MSigDB provides mouse
symbol GMTs directly) and point the YAML at them. Plus a built-in dict of the
stress-relevant gene sets the project specifies (GR targets, HPA, neuroinflam,
synaptic, mito, ER-stress, OXPHOS) — small, curated, no download.

Gene IDs: gene sets use SYMBOLS. If the DE genes are Ensembl IDs, they are
mapped to symbols via var['symbol'] from the annotated h5ad BEFORE enrichment.
If the symbol overlap with the gene sets is near-zero, the script HARD-FAILS
(an Ensembl-vs-symbol mismatch would otherwise silently return no pathways).

Config (YAML `pathways:` block, optional):
  pathways:
    gmt_files:                 # local GMT paths (mouse symbols)
      - /refs/mh.all.v2024.1.Mm.symbols.gmt    # MSigDB mouse hallmark
      # - /refs/m2.cp.reactome.v2024.1.Mm.symbols.gmt
    run_tf_activity: false     # CollecTRI TF activity (needs network); off by default
    min_genes_per_set: 5
  # If gmt_files is empty/absent, only the built-in stress gene sets are used.

Usage:
  uv run python scripts/08c_pathways.py --config config/dev.yaml
  uv run python scripts/08c_pathways.py --config config/brain.yaml

Inputs:
  {results_dir}/tables/de_results.csv           (from 8b)
  {results_dir}/h5ad/08_annotated/all_samples.h5ad   (for var['symbol'] map only)
  GMT files from the YAML pathways.gmt_files (optional)

Outputs:
  {results_dir}/plots/08c_pathways/{contrast}/{level}/{celltype}/
    - gsea_top_pathways.png     : top enriched pathways NAMED, colored by NES
  {results_dir}/tables/pathway_results.csv
    [contrast, flag, group_level, celltype, source(geneset), NES, pvalue, FDR, note]
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

from _utils import load_config


# Optional small supplement of niche stress sets not well captured by MSigDB
# collections. OFF by default (use_builtin_stress_sets=false). The PRIMARY source
# is the MSigDB TSV from fetch_genesets.R (MH + M2/Reactome + M5/GO:BP + M8).
# These are representative members only — flagged UNVERIFIED; enable only if you
# have a vetted list and know why MSigDB doesn't already cover it.
SUPPLEMENT_GENE_SETS = {  # UNVERIFIED — optional, off by default
    "GR_target_genes_custom": ["Fkbp5", "Tsc22d3", "Sgk1", "Zbtb16", "Ddit4", "Per1"],
}


def load_genesets_tsv(path, collections, min_genes):
    """Load the msigdbr export TSV (cols: collection, subcollection, gs_name,
    gene_symbol) into a decoupler net DataFrame [source=gs_name, target=gene,
    collection]. Filters to the requested collections and drops tiny sets."""
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
    return net


def get_symbol_map(adata):
    for col in ("symbol", "gene_symbol", "gene_symbols", "Symbol"):
        if col in adata.var.columns:
            return dict(zip(adata.var_names.astype(str), adata.var[col].astype(str))), col
    return None, None


def run_gsea_on_ranks(rank_series, net, min_genes, times, seed=42):
    """Run GSEA on a single ranking vector (index=gene symbol, values=stat).

    Version-robust across decoupler 2.0 (dc.mt.gsea, matrix/df input) and
    1.9 (dc.get_gsea_df, long-format DE table). Returns a tidy DataFrame with
    columns [source, NES, pvalue, FDR] (names normalized by the caller).
    """
    import decoupler as dc

    # decoupler 2.0: dc.mt.gsea on a 1-row DataFrame (1 "sample" x genes).
    if hasattr(dc, "mt") and hasattr(dc.mt, "gsea"):
        mat = rank_series.to_frame().T          # 1 x n_genes
        mat.index = ["contrast"]
        # 2.0 returns results into the object / as a tuple depending on input type.
        # For a DataFrame, dc.mt.gsea returns (estimate_df, pvalue_df).
        out = dc.mt.gsea(data=mat, net=net, tmin=min_genes, times=times, seed=seed)
        if isinstance(out, tuple):
            est, pval = out[0], out[1]
            res = pd.DataFrame({
                "source": est.columns,
                "NES": est.iloc[0].values,
                "pvalue": pval.iloc[0].values,
            })
        else:  # AnnData-like fallback
            est = out
            res = pd.DataFrame({"source": est.columns, "NES": est.iloc[0].values})
        # BH-adjust
        if "pvalue" in res.columns:
            from scipy.stats import false_discovery_control
            res["FDR"] = false_discovery_control(res["pvalue"].fillna(1.0))
        return res

    # decoupler 1.9 fallback: long-format wrapper.
    df = rank_series.to_frame("stat")
    gsea = dc.get_gsea_df(df, stat="stat", net=net,
                          source="source", target="target",
                          times=times, min_n=min_genes)
    gsea = gsea.reset_index().rename(columns={"index": "source"})
    return gsea


def plot_top_pathways(res, title, out, n=20):
    """Horizontal bar of top pathways by |NES|, NAMED, colored by sign,
    significant (FDR<0.05) outlined."""
    if res.empty:
        return
    nes_col = "NES" if "NES" in res.columns else next((c for c in res.columns if "nes" in c.lower()), None)
    fdr_col = next((c for c in res.columns if c.lower() in ("fdr", "padj", "adj_pvalue")), None)
    if nes_col is None:
        return
    r = res.reindex(res[nes_col].abs().sort_values(ascending=False).index).head(n)
    colors = ["salmon" if v > 0 else "steelblue" for v in r[nes_col]]
    edges = ["black" if (fdr_col and r.loc[i, fdr_col] < 0.05) else "none" for i in r.index]
    fig, ax = plt.subplots(figsize=(7.5, max(3, 0.4 * len(r))))
    ax.barh(r["source"].astype(str), r[nes_col], color=colors, edgecolor=edges, linewidth=1.2)
    ax.axvline(0, color="k", lw=0.8)
    ax.invert_yaxis()
    ax.set_xlabel("NES (red=up in stress, blue=down; black outline=FDR<0.05)")
    ax.set_title(title, fontsize=9)
    fig.tight_layout(); fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Phase 8c: pathway/GSEA on DE results")
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--stat-col", default="stat",
                    help="DE column to rank genes by (default: Wald 'stat')")
    ap.add_argument("--times", type=int, default=1000, help="GSEA permutations")
    ap.add_argument("--subcluster", default=None,
                    help="Run on a 7b subcluster DE table instead of the main one. "
                         "Pass the cell-type slug (e.g. 'microglia') — reads "
                         "de_results_subcluster_{slug}.csv. See 08b --subcluster.")
    args = ap.parse_args()

    print(f"\n=== Phase 8c: Pathway / GSEA (decoupler) ===")
    cfg = load_config(args.config)
    pcfg = cfg.get("pathways", {})
    geneset_tsv = pcfg.get("geneset_tsv", "refs/msigdb_mouse.tsv")
    collections = pcfg.get("collections", ["MH", "M2", "M5", "M8"])
    min_genes = int(pcfg.get("min_genes_per_set", 5))
    run_tf = bool(pcfg.get("run_tf_activity", False))
    use_supplement = bool(pcfg.get("use_builtin_stress_sets", False))

    rdir = Path(cfg["results_dir"])
    # Main DE table, or a 7b-subcluster DE table when --subcluster is given.
    if args.subcluster:
        de_path = rdir / "tables" / f"de_results_subcluster_{args.subcluster}.csv"
        out_suffix = f"_subcluster_{args.subcluster}"
        print(f"  SUBCLUSTER mode: {args.subcluster}")
    else:
        de_path = rdir / "tables" / "de_results.csv"
        out_suffix = ""
    if not de_path.is_file():
        sys.exit(f"ERROR: {de_path} not found. Run 08b_de.py"
                 + (f" --subcluster {args.subcluster}" if args.subcluster else "") + " first.")
    de = pd.read_csv(de_path)
    if de.empty or de["gene"].isna().all():
        sys.exit(f"ERROR: {de_path.name} has no gene-level rows to rank.")

    # Symbol map (gene sets are symbols; DE genes may be Ensembl)
    ann = rdir / "h5ad" / "08_annotated" / "all_samples.h5ad"
    symbol_map = None
    if ann.is_file():
        adata = sc.read_h5ad(ann, backed="r")
        symbol_map, sym_col = get_symbol_map(adata)
        if symbol_map:
            print(f"  Mapping genes Ensembl->symbol via var['{sym_col}']")

    import decoupler as dc
    # Primary gene sets: MSigDB export from fetch_genesets.R.
    if not Path(geneset_tsv).is_file():
        sys.exit(
            f"ERROR: gene-set file not found: {geneset_tsv}\n"
            f"  Generate it once with:  Rscript scripts/fetch_genesets.R --out {geneset_tsv}\n"
            f"  (exports mouse MSigDB {collections} via msigdbr)."
        )
    net = load_genesets_tsv(geneset_tsv, collections, min_genes)
    # Optional niche supplement (off by default)
    if use_supplement:
        sup = pd.DataFrame(
            [(n, g, "SUPPLEMENT") for n, gs in SUPPLEMENT_GENE_SETS.items() for g in gs],
            columns=["source", "target", "collection"])
        net = pd.concat([net, sup], ignore_index=True).drop_duplicates()
    print(f"  Gene sets: {net['source'].nunique()} from {geneset_tsv} "
          f"(collections: {', '.join(collections)})"
          + ("  + builtin supplement" if use_supplement else ""))

    plot_root = rdir / "plots" / ("08c_pathways" + out_suffix)
    plot_root.mkdir(parents=True, exist_ok=True)

    # Overlap sanity: map a sample of DE genes to symbols, check intersection
    # with the gene sets. Near-zero overlap => ID mismatch => hard fail.
    de_genes = de["gene"].dropna().astype(str).unique()
    mapped = ([symbol_map.get(g, g) for g in de_genes] if symbol_map else list(de_genes))
    overlap = len(set(mapped) & set(net["target"]))
    print(f"  Gene overlap (DE symbols ∩ gene sets): {overlap}")
    if overlap < 5:
        sys.exit(
            f"ERROR: only {overlap} DE genes overlap the gene sets.\n"
            f"  Likely a gene-ID mismatch (Ensembl vs symbol) or wrong organism GMT.\n"
            f"  DE gene examples: {list(de_genes[:5])}\n"
            f"  Gene-set target examples: {list(net['target'].unique()[:5])}\n"
            f"  Refusing to run GSEA that would silently return no pathways."
        )

    rows = []
    # iterate per contrast x group_level x celltype
    keys = ["contrast", "group_level", "celltype"]
    for (contrast, level, ct), sub in de.groupby(keys, observed=True):
        sub = sub.dropna(subset=[args.stat_col, "gene"]).copy()
        if sub.shape[0] < 10:
            continue
        sub["gene_sym"] = (sub["gene"].astype(str).map(lambda g: symbol_map.get(g, g))
                           if symbol_map else sub["gene"].astype(str))
        # collapse duplicate symbols (keep max |stat|)
        sub = (sub.reindex(sub[args.stat_col].abs().sort_values(ascending=False).index)
                  .drop_duplicates("gene_sym"))
        rank_series = sub.set_index("gene_sym")[args.stat_col]
        flag = sub["flag"].iloc[0] if "flag" in sub.columns else None
        note = sub["note"].iloc[0] if "note" in sub.columns else ""

        try:
            gsea = run_gsea_on_ranks(rank_series, net, min_genes, args.times)
        except Exception as e:
            print(f"  [warn] GSEA failed for {contrast}|{level}|{ct}: {e}")
            continue
        if gsea is None or gsea.empty:
            continue

        pdir = plot_root / str(contrast) / str(level) / str(ct).replace("/", "_").replace(" ", "_")
        pdir.mkdir(parents=True, exist_ok=True)
        plot_top_pathways(gsea, f"{contrast}\n{level} | {ct}", pdir / "gsea_top_pathways.png")

        nes_col = "NES" if "NES" in gsea.columns else next(
            (c for c in gsea.columns if "nes" in c.lower()), None)
        fdr_col = next((c for c in gsea.columns if c.lower() in ("fdr", "padj")), None)
        p_col = next((c for c in gsea.columns if c.lower() in ("pval", "pvalue", "p_value")), None)
        for _, g in gsea.iterrows():
            rows.append({
                "contrast": contrast, "flag": flag, "group_level": level, "celltype": ct,
                "source": g.get("source"),
                "NES": g.get(nes_col) if nes_col else None,
                "pvalue": g.get(p_col) if p_col else None,
                "FDR": g.get(fdr_col) if fdr_col else None,
                "note": note,
            })

    out_csv = rdir / "tables" / f"pathway_results{out_suffix}.csv"
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"\n  Master table: {out_csv}  ({len(rows)} rows)")

    # Optional TF activity (CollecTRI) — needs network; guarded.
    if run_tf:
        print(f"\n  TF activity (CollecTRI)...")
        try:
            # decoupler 2.0: dc.op.collectri; 1.9: dc.get_collectri
            if hasattr(dc, "op") and hasattr(dc.op, "collectri"):
                collectri = dc.op.collectri(organism="mouse")
            else:
                collectri = dc.get_collectri(organism="mouse", split_complexes=False)
            src_col = "source" if "source" in collectri.columns else collectri.columns[0]
            print(f"    CollecTRI: {collectri[src_col].nunique()} TFs. "
                  f"(per-contrast ULM TF activity wired here in a follow-up.)")
            # NOTE: per-contrast TF activity via dc.run_ulm on the stat vector can
            # be added once gene sets are validated; left as a deliberate stub so
            # this run stays focused on GSEA. Not a silent skip — announced.
        except Exception as e:
            print(f"    [skip] CollecTRI fetch failed (network?): {e}")

    print(f"  Plots: {plot_root}")
    print(f"\n✓ Phase 8c complete.")
    print(f"\n  NOTE: built-in stress gene sets are a SCAFFOLD (representative members,")
    print(f"  flagged UNVERIFIED). Add mouse MSigDB GMTs via pathways.gmt_files and")
    print(f"  refine the stress sets with literature lists before the real run.\n")


if __name__ == "__main__":
    main()
