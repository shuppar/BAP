#!/usr/bin/env python
"""
08f_cross_tissue.py — Phase 8f: Placenta → Brain transcriptional cascade.

The headline cross-tissue analysis. Two biologically aligned arms (project doc §8f):
  - Early arm: E12.5 placenta (Early-vs-Relaxed) → P1/4W/3mo brain (Early-vs-Relaxed)
  - Late  arm: E18.5 placenta (Late-vs-Relaxed)  → P1/4W/3mo brain (Late-vs-Relaxed)
    (Brain P1 Late carries 'confounded_with_pool' flag — propagated to output)

Six analytical views (all operate on completed 8b/8c tables — no re-running DE).
Join is always placenta-WHOLE × brain-{whole + 13 regions} (placenta is whole-only),
so every output row/path carries a `brain_level`:
  1. DEG overlap         — Hypergeometric test of overlapping DEGs per ct_pair
  2. RRHO                — Rank-rank hypergeometric overlap (custom NumPy impl)
  3. Pathway concordance — Same-direction pathway hits per ct_pair (from 8c)
  4. LR cross-tissue     — Placental ligand × brain receptor mechanistic
                           hypotheses (liana mouseconsensus) — KEY publication table
  5. TF concordance      — Same-direction TF activity per ct_pair (from 8c --tf)
  6. Overlap-gene ORA    — MSigDB enrichment of the concordant overlap gene set

Filters (all views): sex stratum = --sex (default combined); placenta level = whole;
brain level iterates whole + regions. DEG cutoff |log2FC| > --logfc-cutoff (default 1.0,
matches 8b) & padj < --padj-cutoff.

A 7th view (bulk-deconvolved sample-level concordance) is deferred to a later
session per user decision; placeholder marked with `# TODO 7th view`.

Inputs:
  results/brain/tables/08b_de/08b_de_results.csv
  results/brain/tables/08c_pathways/08c_pathway_results.csv (+ 08c_tf_activity.csv)
  results/placenta/tables/08b_de/08b_de_results.csv
  results/placenta/tables/08c_pathways/08c_pathway_results.csv (+ 08c_tf_activity.csv)

Outputs:
  plots/08f_cross_tissue/
    01_overview/                                  — overlap heatmaps per arm×age×level
    02_deg_overlap/{arm}_{brain_age}/{brain_level}/
    03_rrho/{arm}_{brain_age}/{brain_level}/
    04_pathway_concordance/{arm}_{brain_age}/{brain_level}/
    05_lr_cross_tissue/{arm}_{brain_age}/{brain_level}/
    06_tf_concordance/{arm}_{brain_age}/{brain_level}/
    07_overlap_enrichment/{arm}_{brain_age}/{brain_level}/

  tables/08f_cross_tissue/
    08f_deg_overlap.csv / 08f_rrho_summary.csv / 08f_pathway_concordance.csv
    08f_lr_cross_tissue.csv (KEY) / 08f_tf_concordance.csv / 08f_overlap_enrichment.csv

Usage:
  uv run python scripts/08f_cross_tissue.py \\
      --brain-config config/brain.yaml \\
      --placenta-config config/placenta.yaml

  # Dev smoke test (brain-only — placenta arms cleanly skip with warning):
  uv run python scripts/08f_cross_tissue.py \\
      --brain-config config/dev_split.yaml \\
      --placenta-config config/dev_split.yaml
"""

import argparse
import sys
import warnings
from itertools import product
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from _utils import load_config, phase_table_dir

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

PHASE = "08f_cross_tissue"


# ============================================================================
# Cross-tissue arm definitions (project doc §8f)
# ============================================================================
# Each arm pairs ONE placenta contrast×age with ONE brain contrast.
# Brain contrast yields multiple ages; we iterate over them.
# Confound flags propagated to output rows.

ARMS = [
    {
        "arm": "Early",
        "placenta_contrast": "early_vs_relaxed_E12.5",
        "placenta_age": "E12.5",
        "brain_contrast": "early_vs_relaxed_per_age",
        "brain_ages": ["P1", "4W", "3mo"],
        "confound_flags": {},  # no per-age confounds for Early arm
    },
    {
        "arm": "Late",
        "placenta_contrast": "late_vs_relaxed_E18.5",
        "placenta_age": "E18.5",
        "brain_contrast": "late_vs_relaxed_per_age",
        "brain_ages": ["P1", "4W", "3mo"],
        "confound_flags": {
            "P1": "Late stress at P1 is Pool3 only; group fully confounded with pool.",
        },
    },
]


# ============================================================================
# CLI
# ============================================================================

def parse_args():
    p = argparse.ArgumentParser(description="Phase 8f: cross-tissue placenta→brain")
    p.add_argument("--brain-config", required=True, type=Path)
    p.add_argument("--placenta-config", required=True, type=Path)
    p.add_argument("--out-results-dir", type=Path, default=None,
                   help="Output dir for tables/plots (default: brain config's results_dir)")
    p.add_argument("--padj-cutoff", type=float, default=0.05,
                   help="FDR cutoff for calling a gene DEG (default 0.05)")
    p.add_argument("--logfc-cutoff", type=float, default=1.0,
                   help="|log2FC| cutoff for DEG calling (default 1.0, matches 8b)")
    p.add_argument("--pathway-fdr-cutoff", type=float, default=0.1,
                   help="FDR cutoff for calling a pathway hit (default 0.1)")
    p.add_argument("--sex", type=str, default="combined",
                   choices=["combined", "M", "F"],
                   help="Sex stratum to use on both tissues (default combined; "
                        "M/F are low_n exploratory per Phase 8 conventions)")
    p.add_argument("--top-n-label", type=int, default=20,
                   help="Hard backstop: max items to plot after the quantile floor (default 20)")
    p.add_argument("--plot-quantile", type=float, default=0.75,
                   help="Plot-only density floor: keep items whose effect statistic is "
                        ">= this quantile WITHIN each plot's own slice (default 0.75 = top "
                        "quartile). Tables are unaffected. Floors the LR bar + pathway/TF "
                        "dotplots only; heatmaps and distribution scatters are never floored.")
    p.add_argument("--dev-test", action="store_true",
                   help="Smoke test mode: treat the 'placenta' input as duplicated "
                        "brain data, remap arms to brain contrasts so joins exercise. "
                        "DO NOT USE on real data.")
    return p.parse_args()


# ============================================================================
# Loading helpers
# ============================================================================

def _slug(s: str) -> str:
    return str(s).replace(" ", "_").replace("/", "-").replace(".", "")


def load_tissue_tables(cfg, tissue_label):
    """Load 08b DE + 08c pathway + 08c TF tables for a tissue.
    Returns (de_df, pw_df, tf_df). Any can be empty if file is missing."""
    base = Path(cfg["results_dir"]) / "tables"
    de_path = base / "08b_de" / "08b_de_results.csv"
    pw_path = base / "08c_pathways" / "08c_pathway_results.csv"
    tf_path = base / "08c_pathways" / "08c_tf_activity.csv"

    de_df = pd.DataFrame()
    pw_df = pd.DataFrame()
    tf_df = pd.DataFrame()

    if de_path.is_file():
        de_df = pd.read_csv(de_path, low_memory=False)
        print(f"  [{tissue_label}] 08b DE: {len(de_df):,} rows from {de_path.name}")
    else:
        print(f"  [{tissue_label}] [warn] 08b DE not found at {de_path}")

    if pw_path.is_file():
        pw_df = pd.read_csv(pw_path, low_memory=False)
        print(f"  [{tissue_label}] 08c pathways: {len(pw_df):,} rows from {pw_path.name}")
    else:
        print(f"  [{tissue_label}] [warn] 08c pathways not found at {pw_path}")

    if tf_path.is_file():
        tf_df = pd.read_csv(tf_path, low_memory=False)
        print(f"  [{tissue_label}] 08c TF activity: {len(tf_df):,} rows from {tf_path.name}")
    else:
        print(f"  [{tissue_label}] [info] 08c TF activity not found "
              f"(run 8c with --tf to enable TF concordance view)")

    return de_df, pw_df, tf_df


def _extract_age(group_level):
    """For the per-age contrasts 8f uses, group_level holds the age DIRECTLY
    (P1 / 4W / 3mo / E12.5 / E18.5) — not an 'age-4W' encoding. Return as-is."""
    return group_level if isinstance(group_level, str) else None


def _apply_strata(sub, sex="combined", level="whole"):
    """Filter a DE/pathway/TF frame to one sex stratum and one level.
    sex/level == None means 'do not filter on that axis' (used for the
    pre-check that asks whether ANY rows exist for a contrast×age)."""
    if sex is not None and "sex" in sub.columns:
        sub = sub[sub["sex"] == sex]
    if level is not None and "level" in sub.columns:
        sub = sub[sub["level"] == level]
    return sub


def prep_de_df(de_df, contrast_name, age_filter=None, sex="combined", level="whole"):
    """Filter to one contrast, sex stratum, level, and optionally one age. Returns
    gene-level rows only with normalised columns: celltype, gene, stat, padj,
    log2FC, age. sex/level=None skips that filter."""
    if de_df.empty:
        return pd.DataFrame()
    if "contrast" not in de_df.columns:
        return pd.DataFrame()

    sub = de_df[de_df["contrast"] == contrast_name].copy()
    if sub.empty:
        return pd.DataFrame()
    sub = _apply_strata(sub, sex=sex, level=level)
    if sub.empty:
        return pd.DataFrame()
    sub = sub.dropna(subset=["gene", "stat"])
    if sub.empty:
        return pd.DataFrame()
    sub["age"] = sub["group_level"].map(_extract_age)
    if age_filter:
        sub = sub[sub["age"] == age_filter]
    return sub[["celltype", "gene", "stat", "padj", "log2FC", "age"]].copy()


def prep_pw_df(pw_df, contrast_name, age_filter=None, sex="combined", level="whole"):
    """Same idea for 08c pathway table.
    08c actual schema: tissue, sex, contrast, flag, group_level, pair, level,
    celltype, collection, source (= pathway NAME), NES, pvalue, FDR, FDR_pooled.
    We rename source→pathway and FDR→padj for consistency in 8f."""
    if pw_df.empty or "contrast" not in pw_df.columns:
        return pd.DataFrame()
    sub = pw_df[pw_df["contrast"] == contrast_name].copy()
    if sub.empty:
        return pd.DataFrame()
    sub = _apply_strata(sub, sex=sex, level=level)
    if sub.empty:
        return pd.DataFrame()
    sub["age"] = sub["group_level"].map(_extract_age)
    if age_filter:
        sub = sub[sub["age"] == age_filter]
    # Normalise column names: 08c uses 'source' for pathway name, 'FDR' for padj
    rename = {}
    if "source" in sub.columns and "pathway" not in sub.columns:
        rename["source"] = "pathway"
    if "FDR" in sub.columns and "padj" not in sub.columns:
        rename["FDR"] = "padj"
    if rename:
        sub = sub.rename(columns=rename)
    return sub


def prep_tf_df(tf_df, contrast_name, age_filter=None, sex="combined", level="whole"):
    """Filter 08c TF activity table.
    08c schema: tissue, sex, contrast, flag, group_level, pair, level, celltype,
    TF, activity_score, pvalue, FDR, direction, ... Rename FDR→padj."""
    if tf_df.empty or "contrast" not in tf_df.columns:
        return pd.DataFrame()
    sub = tf_df[tf_df["contrast"] == contrast_name].copy()
    if sub.empty:
        return pd.DataFrame()
    sub = _apply_strata(sub, sex=sex, level=level)
    if sub.empty:
        return pd.DataFrame()
    sub["age"] = sub["group_level"].map(_extract_age)
    if age_filter:
        sub = sub[sub["age"] == age_filter]
    if "FDR" in sub.columns and "padj" not in sub.columns:
        sub = sub.rename(columns={"FDR": "padj"})
    return sub


# ============================================================================
# Stress-axis gene curation — HPA, glucocorticoid, neuro-inflammation
# Used to flag rows in the LR cross-tissue table that touch the canonical
# prenatal-stress neuroendocrine pathways. Reviewers always ask about these.
# ============================================================================
STRESS_AXIS_GENES = {
    # Glucocorticoid / mineralocorticoid receptor signalling
    "Nr3c1", "Nr3c2", "Fkbp5", "Fkbp4", "Hsp90aa1", "Hsp90ab1",
    # Corticotropin axis (CRH/CRHR)
    "Crh", "Crhr1", "Crhr2", "Crhbp", "Ucn", "Ucn2", "Ucn3",
    # HPA terminal hormones
    "Pomc", "Mc2r", "Mc3r", "Mc4r", "Avp", "Avpr1a", "Avpr1b",
    # Inflammatory cytokines linked to maternal/fetal stress
    "Il6", "Il6r", "Il1b", "Il1r1", "Tnf", "Tnfrsf1a", "Tnfrsf1b",
    "Il17a", "Il17ra", "Ifng", "Ifngr1", "Ifngr2",
    # Steroidogenic / neurosteroid
    "Cyp11a1", "Cyp17a1", "Hsd11b1", "Hsd11b2", "Star",
    # Allopregnanolone / GABA-A modulation (Vacher 2021 axis)
    "Gabra1", "Gabra2", "Gabra3", "Gabra4", "Gabra5", "Gabrb1", "Gabrb2", "Gabrb3",
    # Serotonin (Bonnin 2011, Goeden 2016 axis)
    "Tph1", "Tph2", "Htr1a", "Htr2a", "Htr2c", "Slc6a4",
    # BDNF / neurotrophin signalling
    "Bdnf", "Ntrk2", "Ngf", "Ntrk1", "Ntf3", "Ntrk3",
}


# ============================================================================
# View 1: DEG overlap (hypergeometric per cell-type pair)
# ============================================================================

def _brain_levels(brain_de):
    """All brain levels present (whole + regions), whole first. Placenta is
    whole-only, so every cross-tissue join is placenta-whole × brain-<level>."""
    if brain_de.empty or "level" not in brain_de.columns:
        return ["whole"]
    levs = sorted(brain_de["level"].dropna().astype(str).unique())
    return (["whole"] if "whole" in levs else []) + [l for l in levs if l != "whole"]


def deg_overlap_test(genes_a, genes_b, gene_universe):
    """Hypergeometric test for overlap of two gene sets in a common universe.
    Returns (n_overlap, n_a, n_b, n_universe, p_value)."""
    A = set(genes_a) & set(gene_universe)
    B = set(genes_b) & set(gene_universe)
    overlap = A & B
    n_universe = len(gene_universe)
    n_a, n_b, n_overlap = len(A), len(B), len(overlap)
    if n_a == 0 or n_b == 0 or n_overlap == 0:
        return n_overlap, n_a, n_b, n_universe, 1.0
    # P(>= n_overlap successes) = 1 - CDF(n_overlap - 1)
    p = float(hypergeom.sf(n_overlap - 1, n_universe, n_a, n_b))
    return n_overlap, n_a, n_b, n_universe, p


def run_view1_deg_overlap(brain_de, placenta_de, args):
    """For each arm × brain_age × placenta_celltype × brain_celltype: overlap test
    of DEGs (up and down separately, and pooled). BH-FDR within arm × brain_age.
    Returns master DataFrame.
    """
    from statsmodels.stats.multitest import multipletests
    print("\n[View 1] DEG overlap (hypergeometric)")
    rows = []
    brain_levels = _brain_levels(brain_de)

    for arm in ARMS:
        pl_de = prep_de_df(placenta_de, arm["placenta_contrast"], arm["placenta_age"],
                           sex=args.sex, level="whole")
        if pl_de.empty:
            print(f"  [warn] Arm {arm['arm']}: no placenta DE for "
                  f"{arm['placenta_contrast']} / {arm['placenta_age']} "
                  f"(sex={args.sex}, level=whole) — arm skipped")
            continue
        for br_age in arm["brain_ages"]:
            # Pre-check: does this brain contrast×age exist at all (any level)?
            br_any = prep_de_df(brain_de, arm["brain_contrast"], br_age,
                                sex=args.sex, level=None)
            if br_any.empty:
                print(f"  [warn] Arm {arm['arm']} brain {br_age}: NO DE rows for "
                      f"contrast={arm['brain_contrast']} sex={args.sex} "
                      f"(check contrast name / age encoding)")
                continue

            for br_level in brain_levels:
                br_de = prep_de_df(brain_de, arm["brain_contrast"], br_age,
                                   sex=args.sex, level=br_level)
                if br_de.empty:
                    continue  # region legitimately absent for this slice

                # Universe = genes tested in either tissue (intersection)
                gene_universe = set(pl_de["gene"]) & set(br_de["gene"])
                if len(gene_universe) < 100:
                    continue

                confound = arm["confound_flags"].get(br_age, "")

                for pl_ct, br_ct in product(pl_de["celltype"].unique(),
                                            br_de["celltype"].unique()):
                    pl_sub = pl_de[pl_de["celltype"] == pl_ct]
                    br_sub = br_de[br_de["celltype"] == br_ct]

                    pl_sig = pl_sub[(pl_sub["padj"] < args.padj_cutoff) &
                                    (pl_sub["log2FC"].abs() > args.logfc_cutoff)]
                    br_sig = br_sub[(br_sub["padj"] < args.padj_cutoff) &
                                    (br_sub["log2FC"].abs() > args.logfc_cutoff)]

                    pl_up = set(pl_sig.loc[pl_sig["log2FC"] > 0, "gene"])
                    pl_dn = set(pl_sig.loc[pl_sig["log2FC"] < 0, "gene"])
                    br_up = set(br_sig.loc[br_sig["log2FC"] > 0, "gene"])
                    br_dn = set(br_sig.loc[br_sig["log2FC"] < 0, "gene"])
                    pl_all = pl_up | pl_dn
                    br_all = br_up | br_dn

                    # Three directional tests + one pooled
                    tests = {
                        "concordant_up": (pl_up, br_up),
                        "concordant_down": (pl_dn, br_dn),
                        "discordant_pl_up_br_dn": (pl_up, br_dn),
                        "any_overlap": (pl_all, br_all),
                    }
                    for direction, (a_set, b_set) in tests.items():
                        n_o, n_a, n_b, n_u, p = deg_overlap_test(a_set, b_set,
                                                                  gene_universe)
                        rows.append({
                            "arm": arm["arm"],
                            "brain_age": br_age,
                            "brain_level": br_level,
                            "placenta_celltype": pl_ct,
                            "brain_celltype": br_ct,
                            "direction": direction,
                            "n_overlap": n_o,
                            "n_placenta": n_a,
                            "n_brain": n_b,
                            "n_universe": n_u,
                            "pvalue": p,
                            "overlap_genes": ";".join(sorted(set(a_set) & set(b_set))[:50]),
                            "confound_note": confound,
                        })

    if not rows:
        print("  No DEG overlap tests produced.")
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    # BH-FDR within (arm, brain_age, brain_level, direction)
    fdr_list = []
    for _, g in df.groupby(["arm", "brain_age", "brain_level", "direction"]):
        _, fdr, _, _ = multipletests(g["pvalue"].fillna(1.0), method="fdr_bh")
        fdr_list.append(pd.Series(fdr, index=g.index))
    df["fdr"] = pd.concat(fdr_list)
    df["neg_log10_p"] = -np.log10(df["pvalue"].clip(lower=1e-300))
    print(f"  Total tests: {len(df):,}; FDR<0.05: {(df['fdr'] < 0.05).sum():,}")
    return df


# ============================================================================
# View 2: RRHO (rank-rank hypergeometric overlap)
# ============================================================================

def rrho_matrix(stats_a, stats_b, step=100):
    """Custom NumPy RRHO. Sort both ranked lists by signed Wald stat (descending)
    and scan rank cutoffs in steps; at each (i, j) compute hypergeometric overlap
    of top-i from list A and top-j from list B.

    stats_a, stats_b: pd.Series indexed by gene with signed Wald stats.
    Returns (matrix of -log10(p) shape [k, k], cutoffs).
    Aligns on common gene set first; assumes genes already filtered to common universe.
    """
    common = stats_a.index.intersection(stats_b.index)
    if len(common) < 200:
        return None, None
    a = stats_a.loc[common]
    b = stats_b.loc[common]

    # Rank descending (largest positive = rank 1)
    rank_a = a.rank(ascending=False, method="first")
    rank_b = b.rank(ascending=False, method="first")
    n = len(common)

    cutoffs = np.arange(step, n, step)
    if len(cutoffs) < 3:
        return None, None
    if len(cutoffs) > 40:
        # Cap matrix at 40×40 for plotting performance
        idx = np.linspace(0, len(cutoffs) - 1, 40).astype(int)
        cutoffs = cutoffs[idx]

    # Vectorized RRHO: overlap[i,j] = #genes with rank_a<=cutoffs[i] AND rank_b<=cutoffs[j].
    # Membership matrices (k×n), then a single matmul for all overlaps, then ONE
    # array-valued hypergeom.sf over the whole k×k grid. Identical output to the
    # naive double loop; ~50-100× faster (no per-cell set rebuild, no per-cell scipy call).
    ra = rank_a.to_numpy()
    rb = rank_b.to_numpy()
    cut = cutoffs[:, None]                              # (k, 1)
    a_ind = (ra[None, :] <= cut).astype(np.float64)     # (k, n): top-ci membership of A
    b_ind = (rb[None, :] <= cut).astype(np.float64)     # (k, n): top-cj membership of B
    overlap = a_ind @ b_ind.T                           # (k, k): |top_a_i ∩ top_b_j|
    na = a_ind.sum(axis=1)                              # (k,) = len(top_a) per cutoff
    nb = b_ind.sum(axis=1)                              # (k,) = len(top_b) per cutoff
    na_grid, nb_grid = np.meshgrid(na, nb, indexing="ij")
    p = hypergeom.sf(overlap - 1, n, na_grid, nb_grid)  # array-valued; broadcasts
    mat = -np.log10(np.maximum(p, 1e-300))
    mat[overlap == 0] = 0.0                             # sf(-1)=1 already → 0, set explicit
    return mat, cutoffs


def classify_rrho_concordance(mat, cutoffs):
    """Classify the RRHO result as 'concordant_up', 'concordant_down',
    'discordant', or 'none' based on which quadrant has the strongest signal.
    Concordant up = signal in top-left (top genes both lists).
    Concordant down = signal in bottom-right (bottom genes both lists).
    Discordant = signal in off-diagonal quadrants.
    """
    if mat is None:
        return "none", 0.0
    k = mat.shape[0]
    h = k // 2
    quad_tl = mat[:h, :h].max()         # both top: concordant up
    quad_br = mat[h:, h:].max()         # both bottom: concordant down
    quad_tr = mat[:h, h:].max()         # off: discordant
    quad_bl = mat[h:, :h].max()         # off: discordant
    quads = {"concordant_up": quad_tl, "concordant_down": quad_br,
             "discordant": max(quad_tr, quad_bl)}
    best = max(quads, key=quads.get)
    if quads[best] < 2:  # less than ~p=0.01
        return "none", float(quads[best])
    return best, float(quads[best])


def run_view2_rrho(brain_de, placenta_de, args, plots_root):
    """For each arm × brain_age × placenta_ct × brain_ct: compute RRHO matrix,
    classify concordance, save plot. Returns summary table."""
    print("\n[View 2] RRHO (rank-rank hypergeometric overlap)")
    rows = []
    brain_levels = _brain_levels(brain_de)

    for arm in ARMS:
        pl_de = prep_de_df(placenta_de, arm["placenta_contrast"], arm["placenta_age"],
                           sex=args.sex, level="whole")
        if pl_de.empty:
            print(f"  [warn] Arm {arm['arm']}: no placenta DE "
                  f"(sex={args.sex}, level=whole) — arm skipped")
            continue
        for br_age in arm["brain_ages"]:
            br_any = prep_de_df(brain_de, arm["brain_contrast"], br_age,
                                sex=args.sex, level=None)
            if br_any.empty:
                print(f"  [warn] Arm {arm['arm']} brain {br_age}: NO DE rows "
                      f"(contrast={arm['brain_contrast']}, sex={args.sex})")
                continue
            confound = arm["confound_flags"].get(br_age, "")

            for br_level in brain_levels:
                br_de = prep_de_df(brain_de, arm["brain_contrast"], br_age,
                                   sex=args.sex, level=br_level)
                if br_de.empty:
                    continue
                pdir = plots_root / f"03_rrho/{arm['arm']}_{_slug(br_age)}/{_slug(br_level)}"
                pdir.mkdir(parents=True, exist_ok=True)

                for pl_ct in pl_de["celltype"].unique():
                    pl_stats = (pl_de[pl_de["celltype"] == pl_ct]
                                .set_index("gene")["stat"])
                    if len(pl_stats) < 200:
                        continue
                    for br_ct in br_de["celltype"].unique():
                        br_stats = (br_de[br_de["celltype"] == br_ct]
                                    .set_index("gene")["stat"])
                        if len(br_stats) < 200:
                            continue

                        mat, cutoffs = rrho_matrix(pl_stats, br_stats)
                        if mat is None:
                            continue
                        klass, peak = classify_rrho_concordance(mat, cutoffs)
                        common = pl_stats.index.intersection(br_stats.index)
                        rho, rho_p = spearmanr(pl_stats.loc[common], br_stats.loc[common])

                        rows.append({
                            "arm": arm["arm"], "brain_age": br_age,
                            "brain_level": br_level,
                            "placenta_celltype": pl_ct, "brain_celltype": br_ct,
                            "concordance_class": klass,
                            "peak_neg_log10_p": peak,
                            "spearman_rho": float(rho),
                            "spearman_p": float(rho_p),
                            "n_common_genes": int(len(common)),
                            "confound_note": confound,
                        })

                        if klass != "none":
                            _plot_rrho(mat, cutoffs, pl_ct, br_ct, arm["arm"],
                                       br_age, klass, peak, rho, pdir)

    if not rows:
        print("  No RRHO pairs produced.")
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    print(f"  Total ct_pair RRHOs: {len(df):,}; concordant (any direction): "
          f"{(df['concordance_class'].isin(['concordant_up', 'concordant_down'])).sum():,}")
    return df


def _plot_rrho(mat, cutoffs, pl_ct, br_ct, arm, br_age, klass, peak, rho, pdir):
    fig, ax = plt.subplots(figsize=(5.5, 5))
    im = ax.imshow(mat, origin="upper", cmap="Reds", aspect="auto")
    plt.colorbar(im, ax=ax, label="−log10(p)")
    ax.set_xlabel(f"Brain rank cutoff ({br_ct})")
    ax.set_ylabel(f"Placenta rank cutoff ({pl_ct})")
    ax.set_title(f"RRHO — {arm} arm | brain {br_age}\n"
                 f"{pl_ct} ↔ {br_ct}\n"
                 f"class: {klass}  peak={peak:.1f}  ρ={rho:.2f}", fontsize=9)
    fig.tight_layout()
    fname = f"rrho_{_slug(pl_ct)}_vs_{_slug(br_ct)}.png"
    fig.savefig(pdir / fname, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ============================================================================
# View 3: Pathway concordance
# ============================================================================

def run_view3_pathway_concordance(brain_pw, placenta_pw, args):
    """For each pathway hit (FDR<cutoff) in both tissues: classify directional
    concordance based on NES sign. Returns master table."""
    print("\n[View 3] Pathway concordance")
    rows = []
    brain_levels = _brain_levels(brain_pw)

    for arm in ARMS:
        pl_pw = prep_pw_df(placenta_pw, arm["placenta_contrast"], arm["placenta_age"],
                           sex=args.sex, level="whole")
        if pl_pw.empty:
            continue
        for br_age in arm["brain_ages"]:
            for br_level in brain_levels:
                br_pw = prep_pw_df(brain_pw, arm["brain_contrast"], br_age,
                                   sex=args.sex, level=br_level)
                if br_pw.empty:
                    continue

                confound = arm["confound_flags"].get(br_age, "")

                # Significant pathways per side
                pl_sig = pl_pw[pl_pw["padj"] < args.pathway_fdr_cutoff].copy()
                br_sig = br_pw[br_pw["padj"] < args.pathway_fdr_cutoff].copy()
                if pl_sig.empty or br_sig.empty:
                    continue

                pl_sig["direction"] = np.sign(pl_sig.get("NES", pd.Series(dtype=float))).fillna(0)
                br_sig["direction"] = np.sign(br_sig.get("NES", pd.Series(dtype=float))).fillna(0)

                # Join on pathway × celltype pair
                for pl_ct, br_ct in product(pl_sig["celltype"].unique(),
                                            br_sig["celltype"].unique()):
                    pl_ct_sig = pl_sig[pl_sig["celltype"] == pl_ct]
                    br_ct_sig = br_sig[br_sig["celltype"] == br_ct]
                    common_pw = set(pl_ct_sig["pathway"]) & set(br_ct_sig["pathway"])
                    for pw in common_pw:
                        pl_row = pl_ct_sig[pl_ct_sig["pathway"] == pw].iloc[0]
                        br_row = br_ct_sig[br_ct_sig["pathway"] == pw].iloc[0]
                        direction_a = pl_row["direction"]
                        direction_b = br_row["direction"]
                        if direction_a == 0 or direction_b == 0:
                            klass = "unknown"
                        elif direction_a > 0 and direction_b > 0:
                            klass = "concordant_up"
                        elif direction_a < 0 and direction_b < 0:
                            klass = "concordant_down"
                        else:
                            klass = "discordant"
                        rows.append({
                            "arm": arm["arm"], "brain_age": br_age,
                            "brain_level": br_level,
                            "placenta_celltype": pl_ct, "brain_celltype": br_ct,
                            "pathway": pw,
                            "placenta_NES": pl_row.get("NES"),
                            "brain_NES": br_row.get("NES"),
                            "placenta_padj": pl_row["padj"],
                            "brain_padj": br_row["padj"],
                            "concordance_class": klass,
                            "collection": pl_row.get("collection"),
                            "confound_note": confound,
                        })

    if not rows:
        print("  No pathway concordances produced.")
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    print(f"  Total pathway hits in both tissues: {len(df):,}")
    print(f"  Concordance breakdown:")
    print(df["concordance_class"].value_counts().to_string().replace("\n", "\n    "))
    return df


# ============================================================================
# View 4: LR cross-tissue (mechanistic hypotheses) — the KEY view
# ============================================================================

def run_view4_lr_cross_tissue(brain_de, placenta_de, args):
    """For each LR pair (L, R) in liana mouseconsensus:
      - is L deregulated in placenta (any celltype)?
      - is R deregulated in brain (any celltype)?
    If both pass FDR cutoff with same-direction effect (up/down in BOTH), output
    a mechanistic hypothesis row. This is the publication table.
    """
    print("\n[View 4] LR cross-tissue mechanism hypotheses")
    try:
        import liana as li
        resource = li.rs.select_resource("mouseconsensus")
    except Exception as e:
        print(f"  [skip] could not load liana mouseconsensus: {e}")
        return pd.DataFrame()

    # Resource columns
    l_col = next((c for c in resource.columns if "ligand" in c.lower()), None)
    r_col = next((c for c in resource.columns if "receptor" in c.lower()), None)
    if not l_col or not r_col:
        print(f"  [skip] no ligand/receptor cols in resource. Has: {list(resource.columns)}")
        return pd.DataFrame()
    lr_pairs = resource[[l_col, r_col]].drop_duplicates()
    lr_pairs.columns = ["ligand", "receptor"]
    print(f"  Loaded {len(lr_pairs):,} LR pairs from mouseconsensus")

    rows = []
    brain_levels = _brain_levels(brain_de)

    for arm in ARMS:
        pl_de = prep_de_df(placenta_de, arm["placenta_contrast"], arm["placenta_age"],
                           sex=args.sex, level="whole")
        if pl_de.empty:
            continue
        # Ligand candidates: genes DE in placenta
        pl_sig = pl_de[(pl_de["padj"] < args.padj_cutoff) &
                       (pl_de["log2FC"].abs() > args.logfc_cutoff)]
        if pl_sig.empty:
            continue

        for br_age in arm["brain_ages"]:
            for br_level in brain_levels:
                br_de = prep_de_df(brain_de, arm["brain_contrast"], br_age,
                                   sex=args.sex, level=br_level)
                if br_de.empty:
                    continue
                br_sig = br_de[(br_de["padj"] < args.padj_cutoff) &
                               (br_de["log2FC"].abs() > args.logfc_cutoff)]
                if br_sig.empty:
                    continue

                confound = arm["confound_flags"].get(br_age, "")

                # Match LR pairs where ligand is in pl_sig AND receptor in br_sig
                pl_genes = set(pl_sig["gene"])
                br_genes = set(br_sig["gene"])
                cand = lr_pairs[lr_pairs["ligand"].isin(pl_genes) &
                                lr_pairs["receptor"].isin(br_genes)]
                if cand.empty:
                    continue

                for _, lr in cand.iterrows():
                    lg, rg = lr["ligand"], lr["receptor"]
                    # All (pl_celltype, br_celltype) combinations where both pass
                    lg_rows = pl_sig[pl_sig["gene"] == lg]
                    rg_rows = br_sig[br_sig["gene"] == rg]
                    for _, lr_row in lg_rows.iterrows():
                        for _, rr_row in rg_rows.iterrows():
                            lfc_l = lr_row["log2FC"]
                            lfc_r = rr_row["log2FC"]
                            direction = (
                                "concordant_up" if (lfc_l > 0 and lfc_r > 0) else
                                "concordant_down" if (lfc_l < 0 and lfc_r < 0) else
                                "discordant"
                            )
                            rows.append({
                                "arm": arm["arm"], "brain_age": br_age,
                                "brain_level": br_level,
                                "ligand": lg, "receptor": rg,
                                "placenta_celltype": lr_row["celltype"],
                                "brain_celltype": rr_row["celltype"],
                                "placenta_log2FC": lfc_l,
                                "brain_log2FC": lfc_r,
                                "placenta_padj": lr_row["padj"],
                                "brain_padj": rr_row["padj"],
                                "placenta_stat": lr_row["stat"],
                                "brain_stat": rr_row["stat"],
                                "direction": direction,
                                "stress_axis": (
                                    "ligand+receptor"
                                    if lg in STRESS_AXIS_GENES and rg in STRESS_AXIS_GENES
                                    else "ligand"
                                    if lg in STRESS_AXIS_GENES
                                    else "receptor"
                                    if rg in STRESS_AXIS_GENES
                                    else ""
                                ),
                                "confound_note": confound,
                            })

    if not rows:
        print("  No LR cross-tissue hypotheses found.")
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    n_concord = (df["direction"].isin(["concordant_up", "concordant_down"])).sum()
    print(f"  Total LR cross-tissue rows: {len(df):,}  ({n_concord:,} concordant)")
    return df


# ============================================================================
# View 5: TF concordance (mirror of view 3, on 08c_tf_activity.csv)
# ============================================================================

def run_view5_tf_concordance(brain_tf, placenta_tf, args):
    """For each {arm × brain_age × placenta_ct × brain_ct × TF}: classify
    same-direction TF activation/repression. The intracellular analog of
    pathway concordance — robust at small n; doesn't require ligand transport
    assumptions (TFs are within-cell, so cross-tissue concordance reflects
    shared upstream signalling rather than direct communication)."""
    print("\n[View 5] TF concordance")
    if brain_tf.empty or placenta_tf.empty:
        print("  [skip] TF activity table missing for at least one tissue "
              "(run 8c with --tf to enable).")
        return pd.DataFrame()

    rows = []
    brain_levels = _brain_levels(brain_tf)
    for arm in ARMS:
        pl_tf = prep_tf_df(placenta_tf, arm["placenta_contrast"], arm["placenta_age"],
                           sex=args.sex, level="whole")
        if pl_tf.empty:
            continue
        for br_age in arm["brain_ages"]:
            for br_level in brain_levels:
                br_tf = prep_tf_df(brain_tf, arm["brain_contrast"], br_age,
                                   sex=args.sex, level=br_level)
                if br_tf.empty:
                    continue
                confound = arm["confound_flags"].get(br_age, "")

                pl_sig = pl_tf[pl_tf["padj"] < args.pathway_fdr_cutoff].copy()
                br_sig = br_tf[br_tf["padj"] < args.pathway_fdr_cutoff].copy()
                if pl_sig.empty or br_sig.empty:
                    continue

                for pl_ct, br_ct in product(pl_sig["celltype"].unique(),
                                            br_sig["celltype"].unique()):
                    pl_ct_sig = pl_sig[pl_sig["celltype"] == pl_ct]
                    br_ct_sig = br_sig[br_sig["celltype"] == br_ct]
                    common_tfs = set(pl_ct_sig["TF"]) & set(br_ct_sig["TF"])
                    for tf_name in common_tfs:
                        pl_row = pl_ct_sig[pl_ct_sig["TF"] == tf_name].iloc[0]
                        br_row = br_ct_sig[br_ct_sig["TF"] == tf_name].iloc[0]
                        sa = float(pl_row["activity_score"])
                        sb = float(br_row["activity_score"])
                        if sa > 0 and sb > 0:
                            klass = "concordant_activated"
                        elif sa < 0 and sb < 0:
                            klass = "concordant_repressed"
                        else:
                            klass = "discordant"
                        rows.append({
                            "arm": arm["arm"], "brain_age": br_age,
                            "brain_level": br_level,
                            "placenta_celltype": pl_ct, "brain_celltype": br_ct,
                            "TF": tf_name,
                            "placenta_activity_score": sa,
                            "brain_activity_score": sb,
                            "placenta_padj": float(pl_row["padj"]),
                            "brain_padj": float(br_row["padj"]),
                            "concordance_class": klass,
                            "stress_axis": "TF" if tf_name in STRESS_AXIS_GENES else "",
                            "confound_note": confound,
                        })

    if not rows:
        print("  No TF concordances found.")
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    print(f"  Total TF concordances: {len(df):,}")
    print(f"  Class breakdown:")
    print(df["concordance_class"].value_counts().to_string().replace("\n", "\n    "))
    return df


# ============================================================================
# View 6: Pathway over-representation enrichment of the OVERLAP gene set
# Takes the concordant-up and concordant-down DEGs (genes deregulated in BOTH
# placenta and brain in the same direction) and runs ORA against MSigDB.
# Answers the discussion-figure question: "what biology is shared?"
# ============================================================================

def run_view6_overlap_enrichment(brain_de, placenta_de, args):
    """Pathway ORA on cross-tissue overlap gene sets per arm × brain_age × direction.

    Uses MSigDB mouse via msigdbr if available; else falls back to liana's
    resource pathway annotations as a small reference set."""
    print("\n[View 6] Pathway enrichment of cross-tissue overlap genes")

    # Load gene sets from refs/msigdb_mouse.tsv if present (built by 8c).
    # Schema: gs_name (pathway), gene_symbol, collection
    msigdb_path = Path("refs/msigdb_mouse.tsv")
    if not msigdb_path.is_file():
        print(f"  [skip] {msigdb_path} not found (built by 8c's fetch_genesets.R).")
        print(f"         Run 8c first, or symlink the file from a brain results dir.")
        return pd.DataFrame()

    print(f"  Loading gene sets from {msigdb_path}...")
    msig = pd.read_csv(msigdb_path, sep="\t")
    # Tolerate column name variants
    cols = {c.lower(): c for c in msig.columns}
    gs_col = cols.get("gs_name") or cols.get("pathway") or cols.get("source")
    gene_col = cols.get("gene_symbol") or cols.get("target") or cols.get("gene")
    coll_col = cols.get("collection") or cols.get("gs_cat")
    if not gs_col or not gene_col:
        print(f"  [skip] could not find gs_name/gene_symbol columns. Has: {list(msig.columns)}")
        return pd.DataFrame()
    msig = msig.rename(columns={gs_col: "pathway", gene_col: "gene"})
    if coll_col and coll_col != "collection":
        msig = msig.rename(columns={coll_col: "collection"})
    if "collection" not in msig.columns:
        msig["collection"] = "all"

    # Filter to reasonably-sized gene sets (5–500 genes)
    gs_sizes = msig.groupby("pathway").size()
    keep_gs = gs_sizes[(gs_sizes >= 5) & (gs_sizes <= 500)].index
    msig = msig[msig["pathway"].isin(keep_gs)]
    print(f"  {msig['pathway'].nunique():,} pathways after size filtering")

    from statsmodels.stats.multitest import multipletests
    rows = []
    brain_levels = _brain_levels(brain_de)

    for arm in ARMS:
        pl_de = prep_de_df(placenta_de, arm["placenta_contrast"], arm["placenta_age"],
                           sex=args.sex, level="whole")
        if pl_de.empty:
            continue
        for br_age in arm["brain_ages"]:
            for br_level in brain_levels:
                br_de = prep_de_df(brain_de, arm["brain_contrast"], br_age,
                                   sex=args.sex, level=br_level)
                if br_de.empty:
                    continue
                confound = arm["confound_flags"].get(br_age, "")

                # Universe = genes tested in both tissues (any celltype)
                universe = set(pl_de["gene"]) & set(br_de["gene"])
                if len(universe) < 100:
                    continue
                n_universe = len(universe)

                # Build cross-tissue overlap gene sets, pooled across all celltype pairs
                pl_sig = pl_de[(pl_de["padj"] < args.padj_cutoff) &
                               (pl_de["log2FC"].abs() > args.logfc_cutoff)]
                br_sig = br_de[(br_de["padj"] < args.padj_cutoff) &
                               (br_de["log2FC"].abs() > args.logfc_cutoff)]

                pl_up = set(pl_sig.loc[pl_sig["log2FC"] > 0, "gene"])
                pl_dn = set(pl_sig.loc[pl_sig["log2FC"] < 0, "gene"])
                br_up = set(br_sig.loc[br_sig["log2FC"] > 0, "gene"])
                br_dn = set(br_sig.loc[br_sig["log2FC"] < 0, "gene"])

                overlap_sets = {
                    "concordant_up": (pl_up & br_up) & universe,
                    "concordant_down": (pl_dn & br_dn) & universe,
                    "any_overlap": ((pl_up | pl_dn) & (br_up | br_dn)) & universe,
                }

                for direction, ovset in overlap_sets.items():
                    if len(ovset) < 5:
                        continue
                    # ORA against each pathway: hypergeometric
                    slice_rows = []
                    for pw, pw_block in msig.groupby(["pathway", "collection"]):
                        pathway_name, collection = pw
                        pw_genes = set(pw_block["gene"]) & universe
                        if len(pw_genes) < 5:
                            continue
                        overlap = ovset & pw_genes
                        if len(overlap) == 0:
                            continue
                        p = float(hypergeom.sf(len(overlap) - 1, n_universe,
                                                len(pw_genes), len(ovset)))
                        slice_rows.append({
                            "arm": arm["arm"], "brain_age": br_age,
                            "brain_level": br_level,
                            "direction": direction,
                            "pathway": pathway_name, "collection": collection,
                            "n_overlap": len(overlap),
                            "n_pathway": len(pw_genes),
                            "n_overlap_set": len(ovset),
                            "n_universe": n_universe,
                            "pvalue": p,
                            "overlap_genes": ";".join(sorted(overlap)[:50]),
                            "confound_note": confound,
                        })
                    if not slice_rows:
                        continue
                    slice_df = pd.DataFrame(slice_rows)
                    # BH-FDR within (arm, brain_age, brain_level, direction, collection)
                    fdr_parts = []
                    for _, gg in slice_df.groupby("collection"):
                        _, fdr, _, _ = multipletests(gg["pvalue"].fillna(1.0), method="fdr_bh")
                        fdr_parts.append(pd.Series(fdr, index=gg.index))
                    slice_df["fdr"] = pd.concat(fdr_parts)
                    rows.extend(slice_df.to_dict("records"))

    if not rows:
        print("  No overlap enrichments produced.")
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    n_sig = (df["fdr"] < 0.05).sum()
    print(f"  Total enrichment tests: {len(df):,};  FDR<0.05: {n_sig:,}")
    return df




def _save_fig(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {path.name}")


def _quantile_floor(values, q):
    """Slice-specific adaptive plot floor (mirrors 8e). Given a 1-D array/Series
    of effect magnitudes for ONE plot's slice, return the threshold = q-quantile.
    Items >= threshold are kept. Returns -inf (keep all) when too few to rank."""
    v = pd.Series(values).dropna()
    if len(v) < 4:
        return float("-inf")
    return float(v.quantile(q))


def plot_overview_overlap(deg_df, plots_root):
    """One heatmap per arm: brain_age × concordant_overlap_score.
    Score = max -log10(p) across all ct pairs for direction='any_overlap'."""
    if deg_df.empty:
        return
    sub = deg_df[deg_df["direction"] == "any_overlap"]
    if sub.empty:
        return
    pdir = plots_root / "01_overview"
    pdir.mkdir(parents=True, exist_ok=True)

    for arm in sub["arm"].unique():
        arm_df = sub[sub["arm"] == arm]
        pivot = (arm_df.groupby(["placenta_celltype", "brain_celltype",
                                 "brain_age", "brain_level"])
                 ["neg_log10_p"].max().reset_index())
        for br_age, br_level in pivot[["brain_age", "brain_level"]].drop_duplicates().itertuples(index=False):
            slice_df = pivot[(pivot["brain_age"] == br_age) &
                             (pivot["brain_level"] == br_level)]
            if slice_df.empty:
                continue
            mat = (slice_df.pivot(index="placenta_celltype",
                                   columns="brain_celltype",
                                   values="neg_log10_p").fillna(0))
            if mat.empty:
                continue
            fig, ax = plt.subplots(
                figsize=(max(5, mat.shape[1] * 0.6 + 2),
                         max(4, mat.shape[0] * 0.5 + 2)))
            im = ax.imshow(mat.values, cmap="Reds", aspect="auto")
            plt.colorbar(im, ax=ax, label="max −log10(p), any-direction overlap")
            ax.set_xticks(range(mat.shape[1]))
            ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=8)
            ax.set_yticks(range(mat.shape[0]))
            ax.set_yticklabels(mat.index, fontsize=8)
            ax.set_xlabel("Brain celltype")
            ax.set_ylabel("Placenta celltype")
            ax.set_title(f"DEG overlap: {arm} arm → brain {br_age} [{br_level}]", fontsize=10)
            _save_fig(fig, pdir / f"overlap_overview_{arm}_{_slug(br_age)}_{_slug(br_level)}.png")


def plot_deg_overlap_per_pair(deg_df, plots_root, top_n_label):
    """Per arm×brain_age, a 4-panel heatmap (one per direction)."""
    if deg_df.empty:
        return
    for (arm, br_age, br_level), grp in deg_df.groupby(["arm", "brain_age", "brain_level"]):
        pdir = plots_root / f"02_deg_overlap/{arm}_{_slug(br_age)}/{_slug(br_level)}"
        pdir.mkdir(parents=True, exist_ok=True)
        directions = ["concordant_up", "concordant_down",
                      "discordant_pl_up_br_dn", "any_overlap"]
        for direction in directions:
            sub = grp[grp["direction"] == direction]
            if sub.empty:
                continue
            mat = (sub.pivot(index="placenta_celltype",
                              columns="brain_celltype",
                              values="neg_log10_p").fillna(0))
            fdr_mat = (sub.pivot(index="placenta_celltype",
                                  columns="brain_celltype",
                                  values="fdr").fillna(1))
            fig, ax = plt.subplots(
                figsize=(max(5, mat.shape[1] * 0.6 + 2),
                         max(4, mat.shape[0] * 0.5 + 2)))
            im = ax.imshow(mat.values, cmap="Reds", aspect="auto")
            plt.colorbar(im, ax=ax, label="−log10(p)")
            # Star FDR<0.05 cells
            for i in range(mat.shape[0]):
                for j in range(mat.shape[1]):
                    if fdr_mat.values[i, j] < 0.05:
                        ax.text(j, i, "*", ha="center", va="center",
                                color="white", fontsize=10, fontweight="bold")
            ax.set_xticks(range(mat.shape[1]))
            ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=8)
            ax.set_yticks(range(mat.shape[0]))
            ax.set_yticklabels(mat.index, fontsize=8)
            ax.set_xlabel("Brain celltype")
            ax.set_ylabel("Placenta celltype")
            ax.set_title(f"{direction}  |  {arm} arm → brain {br_age} [{br_level}]\n"
                         f"(* = FDR<0.05)", fontsize=9)
            _save_fig(fig, pdir / f"deg_overlap_{direction}.png")


def plot_pathway_concordance(pw_df, plots_root, top_n_label, plot_q=0.75):
    """Dotplot: pathway × ct_pair, color = direction, size = combined significance.
    ct_pair columns quantile-floored within slice (by max sig_score) + top_n backstop."""
    if pw_df.empty:
        return
    for (arm, br_age, br_level), grp in pw_df.groupby(["arm", "brain_age", "brain_level"]):
        pdir = plots_root / f"04_pathway_concordance/{arm}_{_slug(br_age)}/{_slug(br_level)}"
        pdir.mkdir(parents=True, exist_ok=True)

        # Top pathways by concordance abundance
        concord = grp[grp["concordance_class"].isin(["concordant_up", "concordant_down"])]
        if concord.empty:
            continue
        top_pw = (concord["pathway"].value_counts()
                  .head(top_n_label).index.tolist())
        sub = grp[grp["pathway"].isin(top_pw)].copy()
        sub["ct_pair"] = sub["placenta_celltype"] + " → " + sub["brain_celltype"]
        sub["sig_score"] = -np.log10(sub["placenta_padj"].clip(lower=1e-300)) \
                           - np.log10(sub["brain_padj"].clip(lower=1e-300))
        color_map = {"concordant_up": "#d73027", "concordant_down": "#4575b4",
                     "discordant": "#fdae61", "unknown": "lightgray"}
        sub["color"] = sub["concordance_class"].map(color_map)

        if sub.empty:
            continue
        # Plot-only cap: keep ct_pair columns whose max sig_score >= q-quantile within
        # this slice, then top_n_label backstop. Full set stays in the CSV.
        ct_rank = (sub.groupby("ct_pair")["sig_score"].max()
                   .sort_values(ascending=False))
        n_ct_total = len(ct_rank)
        thr = _quantile_floor(ct_rank, plot_q)
        keep_ct = ct_rank[ct_rank >= thr].head(top_n_label).index.tolist()
        sub = sub[sub["ct_pair"].isin(keep_ct)]
        n_ct_hidden = n_ct_total - len(keep_ct)
        ct_pairs = sorted(sub["ct_pair"].unique())
        pw_list = top_pw
        ct_idx = {c: i for i, c in enumerate(ct_pairs)}
        pw_idx = {p: i for i, p in enumerate(pw_list)}

        fig, ax = plt.subplots(
            figsize=(max(7, len(ct_pairs) * 0.55 + 2),
                     max(6, len(pw_list) * 0.4 + 2)))
        ax.scatter(sub["ct_pair"].map(ct_idx), sub["pathway"].map(pw_idx),
                   s=np.clip(sub["sig_score"] * 10, 30, 400),
                   c=sub["color"], alpha=0.85, edgecolors="k", linewidths=0.3)
        ax.set_xticks(range(len(ct_pairs)))
        ax.set_xticklabels(ct_pairs, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(pw_list)))
        ax.set_yticklabels(pw_list, fontsize=7)
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color=v, label=k) for k, v in color_map.items()],
                  fontsize=7, loc="best")
        ct_note = f"; {n_ct_hidden} more ct_pairs in CSV" if n_ct_hidden else ""
        ax.set_title(f"Pathway concordance — {arm} arm → brain {br_age} [{br_level}]\n"
                     f"top {len(pw_list)} pathways × top {len(ct_pairs)} ct_pairs{ct_note}",
                     fontsize=9)
        _save_fig(fig, pdir / "pathway_concordance.png")


def plot_lr_cross_tissue(lr_df, plots_root, top_n_label, plot_q=0.75):
    """Per arm × brain_age × brain_level:
       1. Scatter (ALL matched pairs — distribution plot, never floored; labels capped)
       2. Top concordant LR bar — quantile-floored within slice + top_n_label backstop
    """
    if lr_df.empty:
        return
    print("\n[Plots] LR cross-tissue mechanism scatters")
    color_map = {"concordant_up": "#d73027", "concordant_down": "#4575b4",
                 "discordant": "#fdae61"}

    for (arm, br_age, br_level), grp in lr_df.groupby(["arm", "brain_age", "brain_level"]):
        pdir = plots_root / f"05_lr_cross_tissue/{arm}_{_slug(br_age)}/{_slug(br_level)}"
        pdir.mkdir(parents=True, exist_ok=True)
        grp = grp.copy()
        grp["color"] = grp["direction"].map(color_map).fillna("lightgray")

        # ---- 1. Scatter (distribution plot: ALL matched pairs, per 8e principle) ----
        grp["prod"] = grp["placenta_log2FC"].abs() * grp["brain_log2FC"].abs()
        draw = grp

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.axhline(0, color="k", lw=0.5, alpha=0.4)
        ax.axvline(0, color="k", lw=0.5, alpha=0.4)
        lim = float(max(draw[["placenta_log2FC", "brain_log2FC"]].abs().max().max(), 1))
        ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.5, alpha=0.3)
        ax.scatter(draw["placenta_log2FC"], draw["brain_log2FC"],
                   c=draw["color"], s=40, alpha=0.7,
                   edgecolors="k", linewidths=0.3)
        # Labels (not points) are capped — privilege pairs deregulated strongly in both
        for _, r in draw.nlargest(top_n_label, "prod").iterrows():
            stress_marker = " ★" if r.get("stress_axis", "") else ""
            ax.annotate(
                f"{r['ligand']}→{r['receptor']}{stress_marker}\n"
                f"[{r['placenta_celltype']} → {r['brain_celltype']}]",
                (r["placenta_log2FC"], r["brain_log2FC"]),
                fontsize=5.5, ha="center", va="bottom",
                xytext=(0, 3), textcoords="offset points",
                arrowprops=dict(arrowstyle="-", color="gray", lw=0.4))
        ax.set_xlabel("Placenta ligand log2FC")
        ax.set_ylabel("Brain receptor log2FC")
        from matplotlib.patches import Patch
        handles = [Patch(color=v, label=k) for k, v in color_map.items()]
        handles.append(Patch(color="white", label="★ touches stress axis",
                              ec="black"))
        ax.legend(handles=handles, fontsize=7)
        ax.set_title(f"LR cross-tissue mechanism — {arm} arm → brain {br_age} [{br_level}]\n"
                     f"(all {len(draw):,} matched pairs; top {top_n_label} labelled "
                     f"by |L log2FC × R log2FC|)", fontsize=9)
        _save_fig(fig, pdir / "lr_cross_tissue_scatter.png")

        # ---- 2. Top concordant bar chart (quantile-floored within slice) ----
        concord = grp[grp["direction"].isin(["concordant_up", "concordant_down"])]
        if concord.empty:
            continue
        thr = _quantile_floor(concord["prod"], plot_q)
        kept = concord[concord["prod"] >= thr]
        top = kept.nlargest(top_n_label, "prod").copy()
        n_more = len(concord) - len(top)
        top["label"] = (
            top["ligand"] + "→" + top["receptor"]
            + top.get("stress_axis", pd.Series("", index=top.index))
              .map(lambda s: " ★" if s else "")
            + "  [" + top["placenta_celltype"] + "→" + top["brain_celltype"] + "]"
        )
        fig, ax = plt.subplots(figsize=(8, max(4, len(top) * 0.35 + 1)))
        ax.barh(range(len(top)), top["prod"], color=top["color"], alpha=0.85)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top["label"].values, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("|Placenta log2FC × Brain log2FC|")
        ax.set_title(f"Top concordant LR mechanisms — {arm} arm → brain {br_age} [{br_level}]\n"
                     f"(q≥{plot_q} floor; {n_more} more concordant pairs in CSV; "
                     f"★ = ligand/receptor in stress axis)", fontsize=9)
        _save_fig(fig, pdir / "lr_cross_tissue_top_concordant.png")


def plot_overlap_effect_size_scatter(brain_de, placenta_de, args, plots_root):
    """Per arm × brain_age: scatter of placenta log2FC vs brain log2FC for genes
    DE in either tissue. RRHO confirms overlap; this shows the *quantitative*
    relationship. Coloured by directional concordance, top genes labelled.
    One plot per arm × brain_age (pooled across celltype pairs)."""
    print("\n[Plots] Effect-size scatters (DEG overlap)")
    brain_levels = _brain_levels(brain_de)
    for arm in ARMS:
        pl_de = prep_de_df(placenta_de, arm["placenta_contrast"], arm["placenta_age"],
                           sex=args.sex, level="whole")
        if pl_de.empty:
            continue
        for br_age in arm["brain_ages"]:
            for br_level in brain_levels:
                br_de = prep_de_df(brain_de, arm["brain_contrast"], br_age,
                                   sex=args.sex, level=br_level)
                if br_de.empty:
                    continue

                # Aggregate to gene-level: take the strongest |log2FC| per gene
                # across celltypes (a gene's "tissue-level" effect).
                def _gene_max(df):
                    return (df.assign(abs_lfc=df["log2FC"].abs())
                            .sort_values("abs_lfc", ascending=False)
                            .drop_duplicates("gene")[["gene", "log2FC", "padj"]])
                pl_g = _gene_max(pl_de).set_index("gene")
                br_g = _gene_max(br_de).set_index("gene")
                common = pl_g.index.intersection(br_g.index)
                if len(common) < 20:
                    continue
                pl_g = pl_g.loc[common].rename(columns=lambda c: "placenta_" + c)
                br_g = br_g.loc[common].rename(columns=lambda c: "brain_" + c)
                merged = pl_g.join(br_g)
                # Only label genes DE in at least one tissue
                merged["sig_either"] = (
                    (merged["placenta_padj"] < args.padj_cutoff) |
                    (merged["brain_padj"] < args.padj_cutoff)
                )
                sig = merged[merged["sig_either"]]
                if len(sig) < 5:
                    continue

                rho, rho_p = spearmanr(sig["placenta_log2FC"], sig["brain_log2FC"])

                def _col(row):
                    pl, br = row["placenta_log2FC"], row["brain_log2FC"]
                    if pl > 0 and br > 0: return "#d73027"
                    if pl < 0 and br < 0: return "#4575b4"
                    return "#fdae61"
                sig = sig.copy()
                sig["color"] = sig.apply(_col, axis=1)

                fig, ax = plt.subplots(figsize=(7, 6))
                ax.axhline(0, color="k", lw=0.5, alpha=0.4)
                ax.axvline(0, color="k", lw=0.5, alpha=0.4)
                lim = float(max(sig[["placenta_log2FC", "brain_log2FC"]].abs().max().max(), 1))
                ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.5, alpha=0.3)
                ax.scatter(sig["placenta_log2FC"], sig["brain_log2FC"],
                           c=sig["color"], s=22, alpha=0.7, edgecolors="k", linewidths=0.2)
                # Label top by product of |logFC|
                sig["prod"] = sig["placenta_log2FC"].abs() * sig["brain_log2FC"].abs()
                top = sig.nlargest(args.top_n_label, "prod")
                for gene, r in top.iterrows():
                    ax.annotate(gene, (r["placenta_log2FC"], r["brain_log2FC"]),
                                fontsize=6, ha="center", va="bottom",
                                xytext=(0, 3), textcoords="offset points",
                                arrowprops=dict(arrowstyle="-", color="gray", lw=0.3))
                ax.set_xlabel("Placenta log2FC (per-gene max |log2FC| across celltypes)")
                ax.set_ylabel("Brain log2FC (per-gene max |log2FC| across celltypes)")
                ax.set_title(
                    f"DEG effect-size concordance — {arm['arm']} arm → brain {br_age} [{br_level}]\n"
                    f"Spearman ρ = {rho:.2f} (p={rho_p:.1e}; n={len(sig):,} genes DE in either)",
                    fontsize=9)
                pdir = plots_root / f"02_deg_overlap/{arm['arm']}_{_slug(br_age)}/{_slug(br_level)}"
                _save_fig(fig, pdir / "effect_size_scatter.png")


def plot_pathway_concordance_scatter(pw_df, plots_root, top_n_label):
    """Per arm × brain_age: scatter of placenta NES vs brain NES for pathways
    significant in both tissues (any cell-type pair). Spearman ρ on plot."""
    if pw_df.empty:
        return
    print("\n[Plots] Effect-size scatters (pathway concordance)")
    for (arm, br_age, br_level), grp in pw_df.groupby(["arm", "brain_age", "brain_level"]):
        # Aggregate: pathway-level mean NES across celltype pairs in this slice
        agg = (grp.groupby("pathway")
               .agg(placenta_NES=("placenta_NES", "mean"),
                    brain_NES=("brain_NES", "mean"),
                    n_ct_pairs=("placenta_NES", "size"))
               .dropna())
        if len(agg) < 5:
            continue
        rho, rho_p = spearmanr(agg["placenta_NES"], agg["brain_NES"])

        def _col(row):
            pl, br = row["placenta_NES"], row["brain_NES"]
            if pl > 0 and br > 0: return "#d73027"
            if pl < 0 and br < 0: return "#4575b4"
            return "#fdae61"
        agg["color"] = agg.apply(_col, axis=1)

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.axhline(0, color="k", lw=0.5, alpha=0.4)
        ax.axvline(0, color="k", lw=0.5, alpha=0.4)
        lim = float(max(agg[["placenta_NES", "brain_NES"]].abs().max().max(), 1))
        ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.5, alpha=0.3)
        ax.scatter(agg["placenta_NES"], agg["brain_NES"],
                   s=20 + agg["n_ct_pairs"] * 4,
                   c=agg["color"], alpha=0.75, edgecolors="k", linewidths=0.2)
        agg["prod"] = agg["placenta_NES"].abs() * agg["brain_NES"].abs()
        for pw, r in agg.nlargest(top_n_label, "prod").iterrows():
            ax.annotate(str(pw)[:60],
                        (r["placenta_NES"], r["brain_NES"]),
                        fontsize=5.5, ha="center", va="bottom",
                        xytext=(0, 3), textcoords="offset points",
                        arrowprops=dict(arrowstyle="-", color="gray", lw=0.3))
        ax.set_xlabel("Placenta NES (mean across ct_pairs)")
        ax.set_ylabel("Brain NES (mean across ct_pairs)")
        ax.set_title(
            f"Pathway NES concordance — {arm} arm → brain {br_age} [{br_level}]\n"
            f"Spearman ρ = {rho:.2f} (p={rho_p:.1e}; n={len(agg):,} pathways)",
            fontsize=9)
        pdir = plots_root / f"04_pathway_concordance/{arm}_{_slug(br_age)}/{_slug(br_level)}"
        _save_fig(fig, pdir / "pathway_nes_scatter.png")


def plot_tf_concordance(tf_df, plots_root, top_n_label, plot_q=0.75):
    """Dotplot of concordant TFs per arm × brain_age × brain_level, plus per-pair
    activity scatter. ct_pair columns quantile-floored within slice + top_n backstop."""
    if tf_df.empty:
        return
    print("\n[Plots] TF concordance")
    color_map = {"concordant_activated": "#d73027",
                 "concordant_repressed": "#4575b4",
                 "discordant": "#fdae61"}
    for (arm, br_age, br_level), grp in tf_df.groupby(["arm", "brain_age", "brain_level"]):
        pdir = plots_root / f"06_tf_concordance/{arm}_{_slug(br_age)}/{_slug(br_level)}"

        # 1. Dotplot: TF × ct_pair, colour = class
        sub = grp.copy()
        concord = sub[sub["concordance_class"].isin(
            ["concordant_activated", "concordant_repressed"])]
        if concord.empty:
            continue
        top_tfs = (concord["TF"].value_counts().head(top_n_label).index.tolist())
        sub = sub[sub["TF"].isin(top_tfs)].copy()
        sub["ct_pair"] = sub["placenta_celltype"] + " → " + sub["brain_celltype"]
        sub["sig_score"] = -np.log10(sub["placenta_padj"].clip(lower=1e-300)) \
                           - np.log10(sub["brain_padj"].clip(lower=1e-300))
        sub["color"] = sub["concordance_class"].map(color_map).fillna("lightgray")

        # Plot-only cap: keep ct_pair columns whose max sig_score >= q-quantile within
        # this slice, then top_n_label backstop. Full set stays in the CSV.
        ct_rank = (sub.groupby("ct_pair")["sig_score"].max()
                   .sort_values(ascending=False))
        n_ct_total = len(ct_rank)
        thr = _quantile_floor(ct_rank, plot_q)
        keep_ct = ct_rank[ct_rank >= thr].head(top_n_label).index.tolist()
        sub = sub[sub["ct_pair"].isin(keep_ct)]
        n_ct_hidden = n_ct_total - len(keep_ct)

        ct_pairs = sorted(sub["ct_pair"].unique())
        tf_list = top_tfs
        ct_idx = {c: i for i, c in enumerate(ct_pairs)}
        tf_idx = {t: i for i, t in enumerate(tf_list)}
        fig, ax = plt.subplots(
            figsize=(max(7, len(ct_pairs) * 0.55 + 2),
                     max(6, len(tf_list) * 0.35 + 2)))
        ax.scatter(sub["ct_pair"].map(ct_idx), sub["TF"].map(tf_idx),
                   s=np.clip(sub["sig_score"] * 10, 30, 400),
                   c=sub["color"], alpha=0.85, edgecolors="k", linewidths=0.3)
        ax.set_xticks(range(len(ct_pairs)))
        ax.set_xticklabels(ct_pairs, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(tf_list)))
        ax.set_yticklabels(tf_list, fontsize=8)
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color=v, label=k) for k, v in color_map.items()],
                  fontsize=7, loc="best")
        ct_note = f"; {n_ct_hidden} more ct_pairs in CSV" if n_ct_hidden else ""
        ax.set_title(f"TF concordance — {arm} arm → brain {br_age} [{br_level}]\n"
                     f"top {len(tf_list)} TFs × top {len(ct_pairs)} ct_pairs{ct_note}",
                     fontsize=9)
        _save_fig(fig, pdir / "tf_concordance_dotplot.png")

        # 2. Activity-score scatter (placenta vs brain), TF-level mean
        agg = (grp.groupby("TF")
               .agg(placenta_act=("placenta_activity_score", "mean"),
                    brain_act=("brain_activity_score", "mean"),
                    n_pairs=("placenta_activity_score", "size"))
               .dropna())
        if len(agg) < 5:
            continue
        rho, rho_p = spearmanr(agg["placenta_act"], agg["brain_act"])
        def _col(row):
            pa, ba = row["placenta_act"], row["brain_act"]
            if pa > 0 and ba > 0: return "#d73027"
            if pa < 0 and ba < 0: return "#4575b4"
            return "#fdae61"
        agg["color"] = agg.apply(_col, axis=1)

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.axhline(0, color="k", lw=0.5, alpha=0.4)
        ax.axvline(0, color="k", lw=0.5, alpha=0.4)
        lim = float(max(agg[["placenta_act", "brain_act"]].abs().max().max(), 1))
        ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.5, alpha=0.3)
        ax.scatter(agg["placenta_act"], agg["brain_act"],
                   s=20 + agg["n_pairs"] * 4,
                   c=agg["color"], alpha=0.75, edgecolors="k", linewidths=0.2)
        agg["prod"] = agg["placenta_act"].abs() * agg["brain_act"].abs()
        for tf_name, r in agg.nlargest(top_n_label, "prod").iterrows():
            ax.annotate(tf_name, (r["placenta_act"], r["brain_act"]),
                        fontsize=6, ha="center", va="bottom",
                        xytext=(0, 3), textcoords="offset points",
                        arrowprops=dict(arrowstyle="-", color="gray", lw=0.3))
        ax.set_xlabel("Placenta TF activity score (mean across ct_pairs)")
        ax.set_ylabel("Brain TF activity score (mean across ct_pairs)")
        ax.set_title(
            f"TF activity concordance — {arm} arm → brain {br_age}\n"
            f"Spearman ρ = {rho:.2f} (p={rho_p:.1e}; n={len(agg):,} TFs)", fontsize=9)
        _save_fig(fig, pdir / "tf_activity_scatter.png")


def plot_overlap_enrichment(enr_df, plots_root, top_n_label):
    """Bar chart of top enriched pathways from the cross-tissue overlap gene
    set, per arm × brain_age × direction × collection. The 'discussion-figure'
    answer to 'what biology is shared'."""
    if enr_df.empty:
        return
    print("\n[Plots] Overlap-gene-set enrichment")
    for (arm, br_age, br_level, direction), grp in enr_df.groupby(
            ["arm", "brain_age", "brain_level", "direction"]):
        sig = grp[grp["fdr"] < 0.05]
        # If nothing FDR-significant in this slice, still draw top-by-p so smoke-test
        # has output to inspect (annotate that nothing was significant).
        if sig.empty:
            top = grp.nsmallest(min(top_n_label, len(grp)), "pvalue").copy()
            sig_marker = " (none FDR<0.05; showing top by p)"
        else:
            top = sig.nsmallest(top_n_label, "fdr").copy()
            sig_marker = ""
        if top.empty:
            continue
        top["neg_log10_fdr"] = -np.log10(top["fdr"].clip(lower=1e-300))
        top = top.sort_values("neg_log10_fdr")
        # Group bars by collection (colour)
        colls = sorted(top["collection"].unique())
        coll_colors = {c: plt.cm.Set2(i / max(len(colls) - 1, 1))
                       for i, c in enumerate(colls)}
        colors = [coll_colors[c] for c in top["collection"]]

        fig, ax = plt.subplots(figsize=(9, max(4, len(top) * 0.35 + 1)))
        ax.barh(range(len(top)), top["neg_log10_fdr"], color=colors, alpha=0.85)
        ax.set_yticks(range(len(top)))
        labels = [f"{r['pathway'][:60]}  ({r['n_overlap']}/{r['n_pathway']})"
                  for _, r in top.iterrows()]
        ax.set_yticklabels(labels, fontsize=7)
        ax.axvline(-np.log10(0.05), color="k", lw=0.6, ls="--")
        ax.set_xlabel("−log10(FDR)")
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color=coll_colors[c], label=c) for c in colls],
                  fontsize=7, loc="best")
        ax.set_title(
            f"Pathways enriched in {direction} cross-tissue DEGs{sig_marker}\n"
            f"{arm} arm → brain {br_age} [{br_level}]  "
            f"(annot = n_overlap / n_pathway_genes)", fontsize=9)
        pdir = plots_root / f"07_overlap_enrichment/{arm}_{_slug(br_age)}/{_slug(br_level)}"
        _save_fig(fig, pdir / f"enrichment_{direction}.png")


# ============================================================================
# Main
# ============================================================================

def main():
    args = parse_args()
    print(f"\n{'='*60}")
    print("Phase 8f: Cross-tissue placenta → brain")
    print(f"{'='*60}")

    global ARMS
    if args.dev_test:
        print("  ⚠ DEV-TEST MODE: arms remapped to brain contrasts; treating "
              "'placenta' input as duplicated brain data. Smoke test only.")
        ARMS = [
            {
                "arm": "Early_devtest",
                "placenta_contrast": "early_vs_relaxed_per_age",
                "placenta_age": "4W",
                "brain_contrast": "early_vs_relaxed_per_age",
                "brain_ages": ["4W"],
                "confound_flags": {},
            },
            {
                "arm": "Late_devtest",
                "placenta_contrast": "late_vs_relaxed_per_age",
                "placenta_age": "4W",
                "brain_contrast": "late_vs_relaxed_per_age",
                "brain_ages": ["4W"],
                "confound_flags": {},
            },
        ]

    print(f"  Brain config:    {args.brain_config}")
    print(f"  Placenta config: {args.placenta_config}")
    brain_cfg = load_config(args.brain_config)
    placenta_cfg = load_config(args.placenta_config)

    # Output location: brain config's results_dir by default
    out_cfg = brain_cfg if args.out_results_dir is None else {
        **brain_cfg, "results_dir": str(args.out_results_dir),
    }
    tdir = phase_table_dir(out_cfg, PHASE)
    pdir_root = Path(out_cfg["results_dir"]) / "plots" / PHASE
    pdir_root.mkdir(parents=True, exist_ok=True)

    print("\n[Loading] Brain + placenta 8b/8c tables")
    brain_de, brain_pw, brain_tf = load_tissue_tables(brain_cfg, "brain")
    placenta_de, placenta_pw, placenta_tf = load_tissue_tables(placenta_cfg, "placenta")

    if brain_de.empty and placenta_de.empty:
        sys.exit("ERROR: neither tissue has 08b DE results. Nothing to do.")
    if brain_de.empty or placenta_de.empty:
        print(f"\n  [info] Only one tissue has DE results — output will be partial.")
        print(f"         (This is expected during dev / before placenta has run.)")

    # ---- Run all six views ----
    deg_df = run_view1_deg_overlap(brain_de, placenta_de, args)
    rrho_df = run_view2_rrho(brain_de, placenta_de, args, pdir_root)
    pw_df = run_view3_pathway_concordance(brain_pw, placenta_pw, args)
    lr_df = run_view4_lr_cross_tissue(brain_de, placenta_de, args)
    tf_df = run_view5_tf_concordance(brain_tf, placenta_tf, args)
    enr_df = run_view6_overlap_enrichment(brain_de, placenta_de, args)

    # TODO 7th view: bulk-deconvolved sample-level concordance.
    # Deferred per user; needs sample-level deconvolution scaffolding.

    # ---- Persist tables ----
    print("\n[Tables]")
    if not deg_df.empty:
        deg_df.to_csv(tdir / "08f_deg_overlap.csv", index=False)
        print(f"  Saved: 08f_deg_overlap.csv  ({len(deg_df):,} rows)")
    if not rrho_df.empty:
        rrho_df.to_csv(tdir / "08f_rrho_summary.csv", index=False)
        print(f"  Saved: 08f_rrho_summary.csv  ({len(rrho_df):,} rows)")
    if not pw_df.empty:
        pw_df.to_csv(tdir / "08f_pathway_concordance.csv", index=False)
        print(f"  Saved: 08f_pathway_concordance.csv  ({len(pw_df):,} rows)")
    if not lr_df.empty:
        lr_df.to_csv(tdir / "08f_lr_cross_tissue.csv", index=False)
        n_axis = (lr_df["stress_axis"] != "").sum() if "stress_axis" in lr_df.columns else 0
        print(f"  Saved: 08f_lr_cross_tissue.csv  ({len(lr_df):,} rows; "
              f"{n_axis:,} touch the stress axis) ← KEY FILE")
    if not tf_df.empty:
        tf_df.to_csv(tdir / "08f_tf_concordance.csv", index=False)
        print(f"  Saved: 08f_tf_concordance.csv  ({len(tf_df):,} rows)")
    if not enr_df.empty:
        enr_df.to_csv(tdir / "08f_overlap_enrichment.csv", index=False)
        n_sig = (enr_df["fdr"] < 0.05).sum()
        print(f"  Saved: 08f_overlap_enrichment.csv  ({len(enr_df):,} rows; "
              f"{n_sig:,} FDR<0.05)")

    # ---- Plots ----
    print("\n[Plots]")
    plot_overview_overlap(deg_df, pdir_root)
    plot_deg_overlap_per_pair(deg_df, pdir_root, args.top_n_label)
    plot_overlap_effect_size_scatter(brain_de, placenta_de, args, pdir_root)
    plot_pathway_concordance(pw_df, pdir_root, args.top_n_label, plot_q=args.plot_quantile)
    plot_pathway_concordance_scatter(pw_df, pdir_root, args.top_n_label)
    plot_lr_cross_tissue(lr_df, pdir_root, args.top_n_label, plot_q=args.plot_quantile)
    plot_tf_concordance(tf_df, pdir_root, args.top_n_label, plot_q=args.plot_quantile)
    plot_overlap_enrichment(enr_df, pdir_root, args.top_n_label)

    # ---- Summary ----
    print(f"\n{'='*60}")
    print("Phase 8f complete.")
    print(f"  Tables: {tdir}")
    print(f"  Plots:  {pdir_root}")
    if not lr_df.empty:
        n_concord = (lr_df["direction"].isin(["concordant_up", "concordant_down"])).sum()
        n_axis = (lr_df["stress_axis"] != "").sum() if "stress_axis" in lr_df.columns else 0
        print(f"\n  KEY RESULT: {n_concord:,} concordant LR cross-tissue mechanisms")
        print(f"             {n_axis:,} touch the canonical stress axis (GR/CRH/cytokines/...)")
        print(f"  See: tables/08f_cross_tissue/08f_lr_cross_tissue.csv")
        for arm in lr_df["arm"].unique():
            for br_age in lr_df[lr_df["arm"] == arm]["brain_age"].unique():
                n = ((lr_df["arm"] == arm) & (lr_df["brain_age"] == br_age) &
                     (lr_df["direction"].isin(["concordant_up", "concordant_down"]))).sum()
                print(f"    {arm} arm → brain {br_age}: {n:,} concordant LR hypotheses")
    print()


if __name__ == "__main__":
    main()
