#!/usr/bin/env python
"""
09_cross_species_validation.py — Phase 9: cross-species validation.

Pseudobulk-DE-vs-DE RRHO2 comparison between your mouse pseudobulk per-cell-type
DE table (from 08b_de_results.csv) and human pseudobulk per-cell-type DE tables
derived from publicly available human snRNA-seq cohorts.

Pipeline:
  1. Load mouse DE table (08b output; one or more contrasts per cell type).
  2. For each human dataset, load its raw counts + cell metadata, pseudobulk
     by (donor × cell_type), run pyDESeq2 between cases and controls,
     produce a per-cell-type human DE table with the SAME schema as 08b.
  3. Map mouse cell types → human cell types via config/cross_species_celltype_map.yaml.
  4. For each matching cell-type pair, run RRHO2 (rank-rank hypergeometric
     overlap) on signed Wald stats. Save heatmap PNGs + summary CSV.

Output (under results/<tissue>/tables/09_cross_species/ and plots/09_cross_species/):
  09_human_de_<dataset>.csv       — DE table per human dataset (pseudobulk)
  09_rrho_summary.csv             — one row per (mouse_contrast, mouse_ct,
                                    human_dataset, human_ct) with concordance
                                    class + peak −log10(p)
  rrho_<...>.png                  — heatmap per significant pair

Usage:
  uv run python scripts/09_cross_species_validation.py \\
      --config config/brain.yaml \\
      --celltype-map config/cross_species_celltype_map.yaml

This script is a SCAFFOLD. The per-human-dataset loaders are stubs — each
needs to be filled in once you've seen the actual downloaded file structure.
Each loader function says exactly what it needs to return; the orchestration
around them is complete.

Honest dependencies on upstream:
  - Mouse DE results from 08b_de_results.csv MUST exist (Phase 8b done).
  - Human datasets downloaded via scripts/download_human_validation.sh.
  - Cell-type map YAML refined to actually match downloaded labels
    (run scripts/preview_human_dataset.py first when available — TBW).
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import scanpy as sc
import anndata as ad
from scipy.stats import hypergeom

import matplotlib.pyplot as plt

from _utils import load_config, phase_table_dir


# ============================================================================
# Mouse DE loader — uses 08b output schema (fixed)
# ============================================================================
def load_mouse_de(cfg: dict) -> pd.DataFrame:
    """Load 08b_de_results.csv produced by 08b_de.py.

    Expected columns (per 08b_de.py):
      contrast_name, level, group_by_level, celltype, gene, log2FoldChange,
      stat (Wald statistic, signed), pvalue, padj, flag, confound_note

    Returns the full table; callers filter by (contrast, celltype) as needed.
    """
    de_path = phase_table_dir(cfg, "08b_de") / "08b_de_results.csv"
    if not de_path.is_file():
        sys.exit(
            f"ERROR: mouse DE table not found at {de_path}.\n"
            f"  Run Phase 8b first (`uv run python scripts/08b_de.py --config {cfg.get('_config_path','...')}`)."
        )
    df = pd.read_csv(de_path)
    required = ["contrast_name", "celltype", "gene", "stat", "padj"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: 08b DE table missing required columns: {missing}")
    return df


# ============================================================================
# Generic pseudobulk DE for a human snRNA-seq dataset
# ============================================================================
def pseudobulk_de_human(
    adata: ad.AnnData,
    donor_col: str,
    celltype_col: str,
    condition_col: str,
    case_label: str,
    control_label: str,
    covariates: list[str] = None,
    min_cells_per_donor_ct: int = 10,
    min_donors_per_group: int = 3,
) -> pd.DataFrame:
    """Pseudobulk by donor × celltype, then run pyDESeq2 per celltype.

    Returns long-form DE table with same schema as mouse 08b:
      contrast_name, level, celltype, gene, log2FoldChange, stat, pvalue, padj, flag

    `flag` is set to 'human_validation_<celltype>' so downstream filtering by
    flag works the same way as mouse.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        sys.exit("ERROR: pydeseq2 not installed. `uv add pydeseq2`")

    if covariates is None:
        covariates = []

    rows = []
    for ct in adata.obs[celltype_col].astype(str).unique():
        sub = adata[adata.obs[celltype_col].astype(str) == ct]
        if sub.n_obs == 0:
            continue
        # Donors with enough cells in this celltype
        per_donor = sub.obs[donor_col].value_counts()
        good_donors = per_donor[per_donor >= min_cells_per_donor_ct].index
        sub = sub[sub.obs[donor_col].isin(good_donors)]
        if sub.n_obs == 0:
            continue
        # Restrict to two conditions
        sub = sub[sub.obs[condition_col].isin([case_label, control_label])]
        # Group sizes
        g = sub.obs.groupby(condition_col, observed=True)[donor_col].nunique()
        if g.get(case_label, 0) < min_donors_per_group or \
           g.get(control_label, 0) < min_donors_per_group:
            print(f"  [skip] {ct}: donor counts {dict(g)} < {min_donors_per_group}")
            continue

        # Build pseudobulk matrix: rows = donors, cols = genes
        donor_ids = sub.obs[donor_col].astype(str).unique()
        pb = np.zeros((len(donor_ids), sub.n_vars), dtype=np.float32)
        for i, d in enumerate(donor_ids):
            cells = sub[sub.obs[donor_col].astype(str) == d]
            X = cells.X
            pb[i] = np.asarray(X.sum(axis=0)).ravel() if hasattr(X, "toarray") \
                else X.sum(axis=0)
        pb_df = pd.DataFrame(pb.astype(int), index=donor_ids, columns=sub.var_names)

        # Donor-level metadata
        meta = (sub.obs[[donor_col, condition_col] + covariates]
                .drop_duplicates(subset=[donor_col])
                .set_index(donor_col)
                .loc[donor_ids])

        # Design: ~ covariate1 + covariate2 + condition
        design_factors = covariates + [condition_col]
        try:
            dds = DeseqDataSet(counts=pb_df, metadata=meta, design_factors=design_factors)
            dds.deseq2()
            stats = DeseqStats(dds, contrast=[condition_col, case_label, control_label])
            stats.summary()
            res = stats.results_df.reset_index().rename(columns={"index": "gene"})
            res["contrast_name"] = f"{case_label}_vs_{control_label}"
            res["level"] = case_label
            res["celltype"] = ct
            res["flag"] = "human_validation"
            rows.append(res[[
                "contrast_name", "level", "celltype", "gene",
                "log2FoldChange", "stat", "pvalue", "padj", "flag"
            ]])
        except Exception as e:
            print(f"  [error] {ct}: {type(e).__name__}: {e}")

    if not rows:
        return pd.DataFrame(columns=[
            "contrast_name", "level", "celltype", "gene",
            "log2FoldChange", "stat", "pvalue", "padj", "flag"
        ])
    return pd.concat(rows, ignore_index=True)


# ============================================================================
# Per-dataset loaders — STUBS, fill in after data lands
# ============================================================================
# Each loader returns:
#   adata (raw counts in .X), donor_col, celltype_col, condition_col,
#   case_label, control_label, covariates
# Then pseudobulk_de_human() handles the rest.
# ----------------------------------------------------------------------------

def load_nagy_2020(data_dir: Path):
    """Nagy 2020 — MDD dlPFC male, GSE144136.

    TODO after download: inspect data_dir for exact files. Expected layout
    (from GEO supplementary): counts matrix as 10X MTX or .csv.gz, cell
    metadata with columns like Sample, Cluster, Condition.

    Returns (adata, donor_col, celltype_col, condition_col, case, control, covariates)
    """
    raise NotImplementedError(
        f"Nagy 2020 loader is a stub. Inspect {data_dir} and fill in:\n"
        f"  - how to read counts (sc.read_mtx / sc.read_csv / sc.read_h5ad)\n"
        f"  - donor_col (likely 'Sample' or 'subject_id')\n"
        f"  - celltype_col (cluster label column)\n"
        f"  - condition_col + case/control labels ('Suicide' vs 'Control'?)\n"
        f"  - sex/age covariates if present"
    )


def load_maitra_2023(data_dir: Path):
    """Maitra 2023 — MDD dlPFC F+M re-analysed with Nagy 2020, GSE213982."""
    raise NotImplementedError(
        f"Maitra 2023 loader is a stub. Same pattern as Nagy 2020 loader.\n"
        f"  Key difference: Maitra 2023 stratifies by sex — covariates should\n"
        f"  include 'sex', and you may want to run sex-stratified DE if your\n"
        f"  mouse contrast is also sex-stratified."
    )


def load_velmeshev_2019(data_dir: Path):
    """Velmeshev 2019 — ASD PFC+ACC, ages 4-22, UCSC autism dataset.

    Downloaded files (per cells.ucsc.edu/autism/downloads.html):
      rawMatrix.zip → matrix.mtx + barcodes.tsv + genes.tsv
      meta.tsv      → cell metadata with columns:
                      cluster, sample, individual, region, age, sex,
                      diagnosis (ASD vs Control), Capbatch, Seqbatch, ...

    Returns (adata, "individual", "cluster", "diagnosis", "ASD", "Control",
             ["sex", "Seqbatch"])
    """
    raise NotImplementedError(
        f"Velmeshev 2019 loader is a stub. Most informed of the four because\n"
        f"  file layout is documented:\n"
        f"    1. Unzip {data_dir}/rawMatrix.zip\n"
        f"    2. sc.read_mtx + attach barcodes/genes\n"
        f"    3. Read meta.tsv and join onto adata.obs\n"
        f"    4. Filter to region == 'PFC' OR keep both PFC+ACC\n"
        f"  Then return:\n"
        f"    adata, 'individual', 'cluster', 'diagnosis', 'ASD', 'Control',\n"
        f"    covariates=['sex', 'Seqbatch']"
    )


def load_herring_2022(data_dir: Path):
    """Herring 2022 — developmental PFC, GSE168408.

    NOTE: Herring 2022 is the dataset to use for cross-AGE comparison of your
    mouse P1/4W/3mo timepoints. The validation isn't case-vs-control but
    "does mouse 4W look like human child/adolescent? does 3mo look like adult?"
    Cross-age validation = different analytic approach (correlate per-age
    mouse signatures vs per-age human signatures across the lifespan).
    Defer this dataset's loader until brain Phase 8b is done and you can
    decide the comparison framing.
    """
    raise NotImplementedError(
        "Herring 2022 loader is a stub. Use only after Phase 8b is done and\n"
        "you've decided the cross-age comparison framing. Defer."
    )


def load_marsh_2022(data_dir: Path):
    """Marsh 2022 — mid-gestation placenta, GSE198373.

    Mid-gestation = closest match for mouse E12. The dataset is normal
    (smooth vs villous chorion); there is NO case/control contrast in this
    dataset, so it's a CELL-TYPE REFERENCE rather than a DE-vs-DE comparator.

    For RRHO comparison against your mouse E12 trophoblast DE, you'd use this
    only as a cell-type-label reference. The real DE comparator for stress
    on placenta is ECHO-PATHWAYS (Cao-Lei 2024, dbGaP controlled).
    """
    raise NotImplementedError(
        "Marsh 2022 loader is a stub. The dataset is normal-only — no\n"
        "case/control contrast. Use for cell-type reference; DE-vs-DE\n"
        "stress comparison must wait for ECHO-PATHWAYS dbGaP approval."
    )


LOADERS = {
    "nagy_2020_GSE144136":     load_nagy_2020,
    "maitra_2023_GSE213982":   load_maitra_2023,
    "velmeshev_2019_autism":   load_velmeshev_2019,
    "herring_2022_GSE168408":  load_herring_2022,
    "marsh_2022_GSE198373":    load_marsh_2022,
}


# ============================================================================
# RRHO2 — adapted from 08f_cross_tissue.py
# ============================================================================
def rrho_matrix(stats_a: pd.Series, stats_b: pd.Series, step: int = 100):
    """Rank-rank hypergeometric overlap on signed Wald stats.
    Returns (matrix of -log10(p) shape [k, k], cutoffs)."""
    common = stats_a.index.intersection(stats_b.index)
    if len(common) < 200:
        return None, None
    a = stats_a.loc[common]
    b = stats_b.loc[common]
    rank_a = a.rank(ascending=False, method="first")
    rank_b = b.rank(ascending=False, method="first")
    n = len(common)

    cutoffs = np.arange(step, n, step)
    if len(cutoffs) < 3:
        return None, None
    if len(cutoffs) > 40:
        idx = np.linspace(0, len(cutoffs) - 1, 40).astype(int)
        cutoffs = cutoffs[idx]

    mat = np.zeros((len(cutoffs), len(cutoffs)))
    for i, ci in enumerate(cutoffs):
        top_a = set(common[rank_a <= ci])
        for j, cj in enumerate(cutoffs):
            top_b = set(common[rank_b <= cj])
            overlap = len(top_a & top_b)
            if overlap == 0:
                continue
            p = hypergeom.sf(overlap - 1, n, len(top_a), len(top_b))
            mat[i, j] = -np.log10(max(p, 1e-300))
    return mat, cutoffs


def classify_rrho(mat: np.ndarray) -> tuple[str, float]:
    """Classify RRHO heatmap as concordant_up, concordant_down, discordant, or none."""
    if mat is None:
        return "none", 0.0
    k = mat.shape[0]
    h = k // 2
    quad_tl = mat[:h, :h].max()   # both top: concordant up
    quad_br = mat[h:, h:].max()   # both bottom: concordant down
    quad_tr = mat[:h, h:].max()   # discordant
    quad_bl = mat[h:, :h].max()
    quads = {"concordant_up": quad_tl, "concordant_down": quad_br,
             "discordant": max(quad_tr, quad_bl)}
    best = max(quads, key=quads.get)
    if quads[best] < 2:           # less than ~p=0.01
        return "none", float(quads[best])
    return best, float(quads[best])


def plot_rrho(mat: np.ndarray, cutoffs: np.ndarray, title: str, out_path: Path):
    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(mat, origin="lower", cmap="viridis", aspect="auto",
                   extent=[cutoffs[0], cutoffs[-1], cutoffs[0], cutoffs[-1]])
    ax.set_xlabel("Rank in human DE (top → up in case)")
    ax.set_ylabel("Rank in mouse DE (top → up in stress)")
    ax.set_title(title, fontsize=9)
    plt.colorbar(im, ax=ax, label="−log10(p)")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ============================================================================
# Cell-type label matching
# ============================================================================
def match_human_labels(human_labels: pd.Series, mapping_entry: dict) -> pd.Series:
    """Given a Series of human cell-type labels and a mapping entry like
    {labels: [...]} or {labels_pattern: '^Ex'}, return boolean mask."""
    if mapping_entry is None:
        return pd.Series(False, index=human_labels.index)
    if "labels" in mapping_entry and mapping_entry["labels"]:
        return human_labels.astype(str).isin(mapping_entry["labels"])
    if "labels_pattern" in mapping_entry:
        return human_labels.astype(str).str.match(mapping_entry["labels_pattern"])
    return pd.Series(False, index=human_labels.index)


# ============================================================================
# Mouse → human gene symbol mapping
# ============================================================================
def mouse_to_human_symbols(mouse_genes: pd.Index, ortholog_table: Path = None) -> pd.Series:
    """Map mouse gene symbols to human orthologs.

    Naive default: uppercase the mouse symbol. ~85% of mouse-human orthologs
    differ only by case (Mecp2 ↔ MECP2, Gfap ↔ GFAP), so this catches the
    majority. For unambiguous accuracy use HGNC's homologene/HCOP table —
    pass `ortholog_table` pointing at a TSV with columns [mouse, human].

    Returns a Series indexed by mouse_genes with human symbol values
    (or NaN if not mappable).
    """
    if ortholog_table is not None and ortholog_table.is_file():
        tbl = pd.read_csv(ortholog_table, sep="\t")
        m2h = dict(zip(tbl["mouse"], tbl["human"]))
        return pd.Series([m2h.get(g, np.nan) for g in mouse_genes], index=mouse_genes)
    # Naive uppercase fallback
    print("  [warn] Using naive uppercase mouse→human mapping. For accuracy,")
    print("         provide a curated ortholog table (HGNC HCOP TSV).")
    return pd.Series(mouse_genes.astype(str).str.upper().values, index=mouse_genes)


# ============================================================================
# Main
# ============================================================================
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--celltype-map", type=Path,
                   default=Path("config/cross_species_celltype_map.yaml"))
    p.add_argument("--human-data-root", type=Path,
                   default=Path("data/human_validation"))
    p.add_argument("--datasets", nargs="+", default=None,
                   help="Subset of human datasets to run (default: all in map).")
    p.add_argument("--ortholog-table", type=Path, default=None,
                   help="TSV with [mouse, human] columns. Falls back to uppercase.")
    p.add_argument("--mouse-contrast", default=None,
                   help="Mouse contrast_name to validate (default: all primary).")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)
    cfg["_config_path"] = str(args.config)
    tissue = cfg["tissue"]

    print(f"\n{'='*60}")
    print(f"Phase 9: cross-species validation — tissue: {tissue}")
    print(f"{'='*60}")

    # -- Mouse DE --
    print(f"\n[1/4] Loading mouse DE table from 08b...")
    mouse_de = load_mouse_de(cfg)
    print(f"  {len(mouse_de):,} mouse DE rows, "
          f"{mouse_de['celltype'].nunique()} cell types, "
          f"{mouse_de['contrast_name'].nunique()} contrasts.")
    if args.mouse_contrast:
        mouse_de = mouse_de[mouse_de["contrast_name"] == args.mouse_contrast]

    # -- Cell-type map --
    with open(args.celltype_map) as f:
        ct_map_all = yaml.safe_load(f)
    ct_map = ct_map_all.get(tissue, {})
    if not ct_map:
        sys.exit(f"ERROR: no celltype map entries for tissue '{tissue}' in {args.celltype_map}")

    # -- Output dirs --
    out_tab = phase_table_dir(cfg, "09_cross_species")
    out_plot = Path(cfg["results_dir"]) / "plots" / "09_cross_species"
    out_plot.mkdir(parents=True, exist_ok=True)

    # -- Datasets to run --
    datasets = args.datasets or sorted(LOADERS.keys())

    summary_rows = []
    for dataset in datasets:
        print(f"\n[2/4] Processing human dataset: {dataset}")
        loader = LOADERS.get(dataset)
        if loader is None:
            print(f"  [skip] no loader registered for {dataset}")
            continue
        data_dir = args.human_data_root / ("brain" if tissue == "brain" else "placenta") / dataset
        if not data_dir.is_dir():
            print(f"  [skip] data not downloaded yet: {data_dir}")
            continue

        try:
            adata, donor_col, ct_col, cond_col, case, ctrl, covs = loader(data_dir)
        except NotImplementedError as e:
            print(f"  [skip] loader is a stub: {e}")
            continue

        # Pseudobulk DE
        print(f"  Running pseudobulk DE: {case} vs {ctrl} ...")
        human_de = pseudobulk_de_human(
            adata, donor_col=donor_col, celltype_col=ct_col,
            condition_col=cond_col, case_label=case, control_label=ctrl,
            covariates=covs,
        )
        if human_de.empty:
            print(f"  [skip] no cell types passed donor-count filters.")
            continue
        # Persist per-dataset DE table
        out_de = out_tab / f"09_human_de_{dataset}.csv"
        human_de.to_csv(out_de, index=False)
        print(f"  Wrote {out_de}")

        # -- RRHO per mouse-celltype × human-celltype --
        print(f"\n[3/4] RRHO2: mouse cell types vs {dataset} ...")
        for mouse_ct, datasets_map in ct_map.items():
            entry = datasets_map.get(dataset)
            if not entry:
                continue
            # Subset mouse DE to this celltype
            m_sub = mouse_de[mouse_de["celltype"] == mouse_ct]
            if m_sub.empty:
                continue

            # Map mouse genes → human symbols once per pair
            m_sub = m_sub.copy()
            m_sub["human_gene"] = mouse_to_human_symbols(
                pd.Index(m_sub["gene"]), args.ortholog_table).values
            m_sub = m_sub.dropna(subset=["human_gene"])

            # Subset human DE to the labels matching this mouse_ct
            h_labels = human_de["celltype"].astype(str)
            mask = match_human_labels(h_labels, entry)
            h_sub = human_de[mask]
            if h_sub.empty:
                continue

            for m_contrast in m_sub["contrast_name"].unique():
                m_stats = (m_sub[m_sub["contrast_name"] == m_contrast]
                           .set_index("human_gene")["stat"])
                for h_ct in h_sub["celltype"].unique():
                    h_stats = h_sub[h_sub["celltype"] == h_ct].set_index("gene")["stat"]
                    mat, cutoffs = rrho_matrix(m_stats, h_stats)
                    cls, peak = classify_rrho(mat)
                    summary_rows.append({
                        "dataset": dataset,
                        "mouse_contrast": m_contrast,
                        "mouse_celltype": mouse_ct,
                        "human_celltype": h_ct,
                        "n_common_genes": int(len(m_stats.index.intersection(h_stats.index))),
                        "concordance": cls,
                        "peak_neg_log10_p": peak,
                    })
                    if cls != "none" and mat is not None:
                        fname = f"rrho_{dataset}_{mouse_ct}_vs_{h_ct}_{m_contrast}.png"
                        fname = re.sub(r"[^A-Za-z0-9._-]", "_", fname)
                        plot_rrho(
                            mat, cutoffs,
                            title=f"{mouse_ct} ({m_contrast}) ↔ {h_ct}\n"
                                  f"{dataset} | {cls} | peak −log10(p)={peak:.2f}",
                            out_path=out_plot / fname,
                        )

    # -- Summary --
    print(f"\n[4/4] Writing summary...")
    if summary_rows:
        summary = pd.DataFrame(summary_rows).sort_values(
            ["mouse_celltype", "dataset", "peak_neg_log10_p"], ascending=[True, True, False]
        )
        summary_path = out_tab / "09_rrho_summary.csv"
        summary.to_csv(summary_path, index=False)
        print(f"  Wrote {summary_path}")
        # Headline
        sig = summary[summary["concordance"] != "none"]
        print(f"\n  {len(sig)}/{len(summary)} cell-type-pair comparisons with concordance ≠ none.")
        if not sig.empty:
            print(f"\n  Top 10 by peak −log10(p):")
            print(sig.head(10)[["dataset", "mouse_celltype", "human_celltype",
                                "concordance", "peak_neg_log10_p"]].to_string(index=False))
    else:
        print("  No RRHO comparisons produced (likely all loaders are still stubs).")

    print(f"\n✓ Phase 9 cross-species validation complete.\n")


if __name__ == "__main__":
    main()
