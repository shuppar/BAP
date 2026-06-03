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
    - gsea_dotplot_panels.png   : dot plots, ONE PANEL PER COLLECTION side by side
    - gsea_volcano_panels.png   : pathway volcanoes, one panel per collection
    - running_<coll>_<pathway>.png : GSEA running-enrichment for top hits per collection
    {contrast}/{level}/celltype_pathway_heatmap_panels.png :
        cell-type x pathway NES heatmap, one panel per collection (shared cell-type axis)
  All panel figures scale width with the number of collections so panels stay
  full-size (no clipping); within-collection FDR throughout.
  {results_dir}/tables/pathway_results.csv
    [contrast, flag, group_level, celltype, source, collection, NES, pvalue,
     FDR(per-collection), FDR_pooled, note]
  FDR is BH-corrected WITHIN each collection (MH/M2/M5/M8) by default — keeps the
  50 Hallmark sets from being buried under thousands of GO:BP sets. FDR_pooled
  (BH across all sets) is kept as a reference column.
"""

import argparse
import re
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


def _safe(name):
    """Filesystem-safe token for folder/file names."""
    return re.sub(r"[^0-9A-Za-z._-]+", "_", str(name)).strip("_")


def get_symbol_map(adata):
    for col in ("symbol", "gene_symbol", "gene_symbols", "Symbol"):
        if col in adata.var.columns:
            return dict(zip(adata.var_names.astype(str), adata.var[col].astype(str))), col
    return None, None


def run_gsea_on_ranks(rank_series, net, min_genes, times, seed=42):
    """Run GSEA on a single ranking vector (index=gene symbol, values=stat).

    Returns a tidy DataFrame [source, NES, pvalue]. FDR is intentionally NOT
    computed here — the caller applies BH correction (both per-collection and
    pooled) once the collection of each set is known, so multiple-testing is
    handled consistently regardless of decoupler version.
    Version-robust across decoupler 2.0 (dc.mt.gsea) and 1.9 (get_gsea_df).
    """
    import decoupler as dc

    if hasattr(dc, "mt") and hasattr(dc.mt, "gsea"):
        mat = rank_series.to_frame().T
        mat.index = ["contrast"]
        out = dc.mt.gsea(data=mat, net=net, tmin=min_genes, times=times, seed=seed)
        if isinstance(out, tuple):
            est, pval = out[0], out[1]
            return pd.DataFrame({"source": est.columns, "NES": est.iloc[0].values,
                                 "pvalue": pval.iloc[0].values})
        est = out
        return pd.DataFrame({"source": est.columns, "NES": est.iloc[0].values,
                             "pvalue": np.nan})

    # decoupler 1.9 fallback
    df = rank_series.to_frame("stat")
    gsea = dc.get_gsea_df(df, stat="stat", net=net, source="source", target="target",
                          times=times, min_n=min_genes).reset_index()
    nes = next((c for c in gsea.columns if "nes" in c.lower()), "NES")
    pv = next((c for c in gsea.columns if c.lower() in ("pval", "pvalue", "p_value")), None)
    src = "source" if "source" in gsea.columns else gsea.columns[0]
    return pd.DataFrame({"source": gsea[src], "NES": gsea[nes],
                         "pvalue": gsea[pv] if pv else np.nan})


def compute_leading_edge(rank_series, members, nes):
    """Leading-edge genes for one pathway: the member genes that drive the
    enrichment, with their ranking stat (DE Wald stat = magnitude + direction).

    GSEA leading edge = members up to the running-sum peak. We approximate it
    directionally from the NES sign (decoupler doesn't return the ES position):
      NES > 0 (enriched at top, upregulated)  -> members with stat > 0
      NES < 0 (enriched at bottom, down)       -> members with stat < 0
    Returns a list of (gene, stat) tuples sorted by |stat| desc.

    The 'stat' here is whatever 8b ranked on (default DESeq2 Wald 'stat'), so
    sign = direction of regulation in the contrast, magnitude = strength.
    """
    present = [g for g in members if g in rank_series.index]
    if not present:
        return []
    sub = rank_series.loc[present]
    if nes is not None and not np.isnan(nes):
        sub = sub[sub > 0] if nes > 0 else sub[sub < 0]
    sub = sub.reindex(sub.abs().sort_values(ascending=False).index)
    return list(zip(sub.index.tolist(), sub.values.tolist()))


def add_fdr(gsea, collection_map):
    """Attach collection, then BH-correct two ways: within each collection
    (the kosher default — keeps small high-quality collections from being
    buried by large redundant ones) and pooled across all sets (reference)."""
    from scipy.stats import false_discovery_control
    g = gsea.copy()
    g["collection"] = g["source"].map(collection_map).fillna("NA")
    p = g["pvalue"].fillna(1.0).values
    # pooled
    g["FDR_pooled"] = false_discovery_control(p) if np.isfinite(p).any() else np.nan
    # per-collection
    g["FDR"] = np.nan
    for coll, idx in g.groupby("collection").groups.items():
        sub = g.loc[idx, "pvalue"].fillna(1.0).values
        g.loc[idx, "FDR"] = false_discovery_control(sub) if len(sub) else np.nan
    return g


def _draw_dotplot(ax, gsea, set_sizes, title, n=20, fdr_thr=0.05):
    """Draw a dot plot onto a given axes. Returns the scatter handle (for colorbar)."""
    g = gsea.dropna(subset=["NES"]).copy()
    if g.empty:
        ax.text(0.5, 0.5, "no enriched sets", ha="center", va="center",
                transform=ax.transAxes, fontsize=8, color="gray")
        ax.set_title(title, fontsize=9); ax.axis("off")
        return None
    g["neglog10fdr"] = -np.log10(g["FDR"].clip(lower=1e-300))
    g = g.reindex(g["FDR"].fillna(1).sort_values().index).head(n)
    sizes = np.array([set_sizes.get(s, 10) for s in g["source"]], float)
    sizes = 30 + 220 * (sizes - sizes.min()) / (np.ptp(sizes) + 1e-9)
    sc_ = ax.scatter(g["NES"], range(len(g)), s=sizes, c=g["neglog10fdr"],
                     cmap="viridis", edgecolor="k", linewidth=0.4, zorder=3)
    ax.set_yticks(range(len(g))); ax.set_yticklabels(g["source"].astype(str), fontsize=7)
    ax.invert_yaxis(); ax.axvline(0, color="k", lw=0.8, zorder=1)
    ax.set_xlabel("NES (>0 up in stress)")
    ax.set_title(title, fontsize=9)
    return sc_


def _draw_volcano(ax, gsea, title, fdr_thr=0.05, max_labels=20):
    """Draw a pathway volcano onto a given axes."""
    g = gsea.dropna(subset=["NES", "FDR"]).copy()
    if g.empty:
        ax.text(0.5, 0.5, "no enriched sets", ha="center", va="center",
                transform=ax.transAxes, fontsize=8, color="gray")
        ax.set_title(title, fontsize=9); ax.axis("off")
        return
    g["neglog10fdr"] = -np.log10(g["FDR"].clip(lower=1e-300))
    sig = g["FDR"] < fdr_thr
    ax.scatter(g.loc[~sig, "NES"], g.loc[~sig, "neglog10fdr"], s=10,
               color="lightgray", rasterized=True)
    ax.scatter(g.loc[sig, "NES"], g.loc[sig, "neglog10fdr"], s=16,
               c=["salmon" if v > 0 else "steelblue" for v in g.loc[sig, "NES"]])
    ax.axhline(-np.log10(fdr_thr), color="k", lw=0.6, ls="--")
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("NES (red=up, blue=down)"); ax.set_ylabel("-log10 FDR")
    ax.set_title(title, fontsize=9)
    lab = g[sig].reindex(g[sig]["FDR"].sort_values().index).head(max_labels)
    texts = [ax.text(r["NES"], r["neglog10fdr"], str(r["source"])[:45], fontsize=6)
             for _, r in lab.iterrows()]
    if texts:
        try:
            from adjustText import adjust_text
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.3))
        except ImportError:
            pass
    ax.text(0.02, 0.98, f"{int(sig.sum())} sig", transform=ax.transAxes,
            fontsize=6, va="top", color="gray")


def _draw_heatmap(ax, per_ct_coll, title, n_paths=25, fdr_thr=0.05):
    """Draw a cell-type x pathway NES heatmap onto a given axes. Returns image handle."""
    frames = []
    for ct, g in per_ct_coll.items():
        gg = g[["source", "NES", "FDR"]].copy(); gg["celltype"] = ct
        frames.append(gg)
    if not frames:
        ax.text(0.5, 0.5, "no enriched sets", ha="center", va="center",
                transform=ax.transAxes, fontsize=8, color="gray")
        ax.set_title(title, fontsize=9); ax.axis("off")
        return None
    long = pd.concat(frames, ignore_index=True)
    top_paths = (long.sort_values("FDR").drop_duplicates("source").head(n_paths)["source"].tolist())
    long = long[long["source"].isin(top_paths)]
    nes_mat = long.pivot_table(index="source", columns="celltype", values="NES").reindex(top_paths)
    fdr_mat = long.pivot_table(index="source", columns="celltype", values="FDR").reindex(top_paths)
    if nes_mat.empty:
        ax.text(0.5, 0.5, "no enriched sets", ha="center", va="center",
                transform=ax.transAxes, fontsize=8, color="gray")
        ax.set_title(title, fontsize=9); ax.axis("off")
        return None
    vmax = np.nanmax(np.abs(nes_mat.values)) or 1
    im = ax.imshow(nes_mat.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(nes_mat.shape[1]))
    ax.set_xticklabels(nes_mat.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(nes_mat.shape[0])); ax.set_yticklabels(nes_mat.index, fontsize=6)
    for i, path in enumerate(nes_mat.index):
        for j, ct in enumerate(nes_mat.columns):
            f = (fdr_mat.loc[path, ct] if (path in fdr_mat.index and ct in fdr_mat.columns)
                 else np.nan)
            if pd.notna(f) and f < fdr_thr:
                ax.text(j, i, "•", ha="center", va="center", fontsize=9, color="black")
    ax.set_title(title, fontsize=9)
    return im


def panel_by_collection(gsea_full, kind, title, out, **kw):
    """One wide figure with a subplot per collection (paper-panel style).

    kind: 'dotplot' | 'volcano'. gsea_full has a 'collection' column. The figure
    width scales with the number of collections present, so each panel stays
    full-size (no shrinking, no clipping). FDR is already per-collection.
    """
    colls = [c for c in ["MH", "M2", "M5", "M8"] if c in set(gsea_full["collection"])]
    colls += [c for c in sorted(set(gsea_full["collection"])) if c not in colls]
    if not colls:
        return
    per_w = 6.5 if kind == "volcano" else 7.0
    fig, axes = plt.subplots(1, len(colls), figsize=(per_w * len(colls), 6.0),
                             squeeze=False, constrained_layout=True)
    axes = axes[0]
    last_handle = None
    for ax, coll in zip(axes, colls):
        gc = gsea_full[gsea_full["collection"] == coll]
        if kind == "dotplot":
            last_handle = _draw_dotplot(ax, gc, kw["set_sizes"], coll) or last_handle
        else:
            _draw_volcano(ax, gc, coll)
    if kind == "dotplot" and last_handle is not None:
        cb = fig.colorbar(last_handle, ax=axes, pad=0.01, fraction=0.025)
        cb.set_label("-log10 FDR", fontsize=8)
    fig.suptitle(f"{title}  —  per-collection panels (within-collection FDR)", fontsize=10)
    # constrained_layout handles suptitle/colorbar spacing without clipping;
    # do NOT also pass bbox_inches='tight' (the two fight each other).
    fig.savefig(out, dpi=140); plt.close(fig)


def panel_heatmap_by_collection(per_ct, title, out, n_paths=25, fdr_thr=0.05):
    """One wide figure: a cell-type x pathway heatmap subplot per collection."""
    colls_all = sorted({c for g in per_ct.values() for c in g["collection"].unique()})
    colls = [c for c in ["MH", "M2", "M5", "M8"] if c in colls_all]
    colls += [c for c in colls_all if c not in colls]
    if not colls:
        return
    fig, axes = plt.subplots(1, len(colls), figsize=(7.0 * len(colls), 8.0),
                             squeeze=False, constrained_layout=True)
    axes = axes[0]
    last_im = None
    for ax, coll in zip(axes, colls):
        per_ct_coll = {ct: g[g["collection"] == coll] for ct, g in per_ct.items()}
        per_ct_coll = {ct: g for ct, g in per_ct_coll.items() if not g.empty}
        im = _draw_heatmap(ax, per_ct_coll, coll, n_paths=n_paths, fdr_thr=fdr_thr)
        last_im = im or last_im
    if last_im is not None:
        cb = fig.colorbar(last_im, ax=axes, pad=0.01, fraction=0.02)
        cb.set_label("NES", fontsize=8)
    fig.suptitle(f"{title}  —  NES by cell type x pathway (• = FDR<{fdr_thr}, per-collection)",
                 fontsize=10)
    fig.savefig(out, dpi=140); plt.close(fig)


def plot_running_enrichment(rank_series, members, title, out):
    """Classic GSEA running-enrichment 'mountain' plot for ONE pathway.
    (Single-pathway by nature — not panelled.)"""
    r = rank_series.sort_values(ascending=False)
    genes = r.index.to_numpy()
    in_set = np.isin(genes, list(members))
    if in_set.sum() < 2:
        return
    scores = np.abs(r.values) ** 1.0
    hit_norm = scores[in_set].sum()
    miss_norm = (~in_set).sum()
    inc = np.where(in_set, scores / (hit_norm + 1e-12), -1.0 / (miss_norm + 1e-12))
    running = np.cumsum(inc)
    es_idx = np.argmax(np.abs(running))
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7, 5), height_ratios=[3, 1], sharex=True)
    a1.plot(running, color="green", lw=1.3)
    a1.axhline(0, color="k", lw=0.6)
    a1.scatter([es_idx], [running[es_idx]], color="red", zorder=5,
               label=f"ES={running[es_idx]:.2f}")
    a1.set_ylabel("running enrichment"); a1.legend(fontsize=7, loc="best")
    a1.set_title(title, fontsize=9)
    a2.vlines(np.where(in_set)[0], 0, 1, color="black", lw=0.5)
    a2.set_yticks([]); a2.set_xlabel("gene rank (high stat -> low)")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Phase 8c: pathway/GSEA on DE results")
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--stat-col", default="stat",
                    help="DE column to rank genes by (default: Wald 'stat')")
    ap.add_argument("--times", type=int, default=1000, help="GSEA permutations")
    ap.add_argument("--le-fdr", type=float, default=0.05,
                    help="FDR cutoff for which pathways get leading-edge gene rows "
                         "(default 0.05). Leading-edge written to pathway_leading_edge.csv.")
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

    set_sizes = net.groupby("source").size().to_dict()
    collection_map = dict(net.drop_duplicates("source").set_index("source")["collection"])
    n_top_running = int(pcfg.get("n_running_enrichment", 3))

    rows = []
    le_rows = []   # leading-edge: one row per contrast×level×celltype×pathway×gene
    # Group by contrast x group_level so we can build a cross-cell-type heatmap
    # per contrast, while still emitting per-cell-type dot/volcano/mountain plots.
    for (contrast, level), block in de.groupby(["contrast", "group_level"], observed=True):
        per_ct = {}   # celltype -> gsea df (with FDR) for the heatmap
        for ct, sub in block.groupby("celltype", observed=True):
            sub = sub.dropna(subset=[args.stat_col, "gene"]).copy()
            if sub.shape[0] < 10:
                continue
            sub["gene_sym"] = (sub["gene"].astype(str).map(lambda g: symbol_map.get(g, g))
                               if symbol_map else sub["gene"].astype(str))
            sub = (sub.reindex(sub[args.stat_col].abs().sort_values(ascending=False).index)
                      .drop_duplicates("gene_sym"))
            rank_series = sub.set_index("gene_sym")[args.stat_col].astype(float)
            # gene_sym -> log2FC map (for leading-edge magnitude/direction)
            lfc_map = (sub.set_index("gene_sym")["log2FC"].astype(float).to_dict()
                       if "log2FC" in sub.columns else {})
            flag = sub["flag"].iloc[0] if "flag" in sub.columns else None
            note = sub["note"].iloc[0] if "note" in sub.columns else ""

            try:
                gsea = run_gsea_on_ranks(rank_series, net, min_genes, args.times)
            except Exception as e:
                print(f"  [warn] GSEA failed for {contrast}|{level}|{ct}: {e}")
                continue
            if gsea is None or gsea.empty:
                continue
            gsea = add_fdr(gsea, collection_map)   # per-collection + pooled FDR
            per_ct[ct] = gsea

            cts = str(ct).replace("/", "_").replace(" ", "_")
            pdir = plot_root / str(contrast) / str(level) / cts
            pdir.mkdir(parents=True, exist_ok=True)
            ttl = f"{contrast} | {level} | {ct}"
            # Multi-panel figures: one subplot per collection, side by side.
            # Figure width scales with #collections so panels stay full-size.
            panel_by_collection(gsea, "dotplot", ttl, pdir / "gsea_dotplot_panels.png",
                                set_sizes=set_sizes)
            panel_by_collection(gsea, "volcano", ttl, pdir / "gsea_volcano_panels.png")
            # Running-enrichment: top significant per collection (single-pathway plots).
            for coll, gcoll in gsea.groupby("collection", observed=True):
                top_sig = gcoll.dropna(subset=["FDR"]).sort_values("FDR").head(n_top_running)
                for _, tp in top_sig.iterrows():
                    members = set(net.loc[net["source"] == tp["source"], "target"])
                    plot_running_enrichment(
                        rank_series, members,
                        f"{ct} | {coll} | {str(tp['source'])[:45]} (FDR={tp['FDR']:.2g})",
                        pdir / f"running_{_safe(coll)}_{_safe(str(tp['source'])[:35])}.png")

            for _, g in gsea.iterrows():
                rows.append({
                    "contrast": contrast, "flag": flag, "group_level": level, "celltype": ct,
                    "source": g["source"], "collection": g["collection"],
                    "NES": g["NES"], "pvalue": g["pvalue"],
                    "FDR": g["FDR"], "FDR_pooled": g["FDR_pooled"], "note": note,
                })

            # Leading-edge genes for significant pathways (FDR < le_fdr).
            # One row per gene driving each enriched pathway, with log2FC + stat.
            sig_gsea = gsea.dropna(subset=["FDR"])
            sig_gsea = sig_gsea[sig_gsea["FDR"] < args.le_fdr]
            for _, g in sig_gsea.iterrows():
                members = set(net.loc[net["source"] == g["source"], "target"])
                le = compute_leading_edge(rank_series, members, g["NES"])
                for rank_i, (gene, stat_val) in enumerate(le, start=1):
                    le_rows.append({
                        "contrast": contrast, "group_level": level, "celltype": ct,
                        "collection": g["collection"], "pathway": g["source"],
                        "NES": round(float(g["NES"]), 4), "pathway_FDR": g["FDR"],
                        "leading_edge_rank": rank_i,
                        "gene": gene,
                        "log2FC": round(lfc_map.get(gene, np.nan), 4),
                        "rank_stat": round(float(stat_val), 4),
                        "direction": "up" if stat_val > 0 else "down",
                        "flag": flag, "note": note,
                    })

        # Cross-cell-type heatmap per contrast x level — one panel per collection.
        if per_ct:
            ldir = plot_root / str(contrast) / str(level)
            ldir.mkdir(parents=True, exist_ok=True)
            panel_heatmap_by_collection(per_ct, f"{contrast} | {level}",
                                        ldir / "celltype_pathway_heatmap_panels.png")

    df_rows = pd.DataFrame(rows)
    out_csv = rdir / "tables" / f"pathway_results{out_suffix}.csv"
    df_rows.to_csv(out_csv, index=False)
    print(f"\n  Master table: {out_csv}  ({len(rows)} rows)")
    # Per-collection sub-tables (paper-panel style)
    if not df_rows.empty:
        for coll, sub in df_rows.groupby("collection"):
            sub_path = rdir / "tables" / f"pathway_results{out_suffix}_{_safe(str(coll))}.csv"
            sub.to_csv(sub_path, index=False)
        print(f"  Per-collection tables: pathway_results{out_suffix}_<collection>.csv "
              f"({df_rows['collection'].nunique()} collections)")

    # Leading-edge genes per significant pathway (the genes driving each result,
    # with log2FC magnitude + direction). One row per gene; join to the
    # 08b expression matrix on (celltype, gene) for per-sample levels.
    le_df = pd.DataFrame(le_rows)
    le_csv = rdir / "tables" / f"pathway_leading_edge{out_suffix}.csv"
    le_df.to_csv(le_csv, index=False)
    if le_rows:
        n_paths = le_df.groupby(["contrast", "group_level", "celltype", "pathway"]).ngroups
        print(f"  Leading-edge table: {le_csv}  ({len(le_rows)} gene rows, "
              f"{n_paths} significant pathways at FDR<{args.le_fdr})")
    else:
        print(f"  Leading-edge table: {le_csv}  (empty — no pathways below FDR<{args.le_fdr})")

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
