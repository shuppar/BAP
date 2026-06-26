#!/usr/bin/env python
"""
08g_cross_age.py — Phase 8g: Cross-age persistence analysis.

Operates entirely on existing 8b/8c tables — no re-running of DE/GSEA.
Comprehensive persistence (views 1-3,5,6) PLUS three complementary analyses
(B trajectory shape, C persistence×disruption, View 7 8f-bridge). Everything is
region-resolved: classification runs per (celltype, LEVEL, feature, arm) where
level ∈ {whole + 13 brain regions}; whole = robust, regions flagged
`regional_exploratory` (pool-age confounded, fewer donors).

Two stress arms:
  - Early: early_vs_relaxed_per_age at {P1, 4W, 3mo}
  - Late:  late_vs_relaxed_per_age at {P1, 4W, 3mo}  (P1 carries pool confound)

Persistence classes per gene (or pathway / TF) × celltype × level × arm:
  persistent      — DE at P1 AND 4W AND 3mo, SAME DIRECTION
  resolving_early — DE at P1 AND 4W, not 3mo
  established_late— DE at 4W AND 3mo, not P1
  P1_only / transient_4W / emergent_3mo — single-age
  P1_3mo_only     — DE at P1 AND 3mo, not 4W (unusual)
  persistent_directionswap — all 3 ages but direction flips
  none            — not DE at any age

Comprehensive views:
  1. Gene-level persistence       (08b DE)
  2. Pathway-level persistence    (08c GSEA, all collections)
  3. TF-level persistence         (08c TF activity)
  5. Early vs Late at each age    (hypergeom overlap + rank-rank)
  6. Cross-arm core signature     (persistent in BOTH arms — paper table)

Complementary analyses (added this session):
  B. Trajectory shape             — amplifying / stable / attenuating across age
                                    (column on the gene/pathway/TF tables)
  C. Persistence × disruption     — 8g persistence × 8b developmental-disruption
                                    classes; headline = persistent AND LOST
  7. 8f-complementary focal module — all-Hallmark/immune pathway + broad-family TF
                                    persistence flagged with 8f cross-tissue
                                    concordance (the bridge), plus a leading-edge
                                    deep-dive on the persistent immune pathways.

Outputs:
  plots/08g_cross_age/{01_gene,02_pathway,03_tf}_persistence/[<region>/]
    04_early_vs_late/[<region>/]  05_core_signature/
    06_persistence_x_disruption/  07_focal_8f_bridge/
  tables/08g_cross_age/
    08g_gene_persistence.csv            KEY: gene × celltype × level × arm × class
    08g_pathway_persistence.csv  08g_tf_persistence.csv
    08g_early_vs_late_overlap.csv
    08g_core_signature_genes.csv        cross-arm persistent (paper table)
    08g_core_signature_pathways.csv
    08g_persistence_x_disruption.csv    [C]
    08g_focal_pathway_persistence.csv   [A: 8f bridge]  08g_focal_tf_persistence.csv
    08g_focal_leadingedge_drivers.csv   [A: LE drivers of persistent immune pathways]

Usage:
  uv run python scripts/08g_cross_age.py --config config/brain.yaml
  # placenta has incomplete factorial (no cross-age comparison) → brain-only.
"""

import argparse
import sys
import warnings
from itertools import combinations
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

PHASE = "08g_cross_age"


# ============================================================================
# Arm definitions (brain only — placenta has incomplete cross-age factorial)
# ============================================================================

ARMS = [
    {
        "arm": "Early",
        "contrast": "early_vs_relaxed_per_age",
        "ages": ["P1", "4W", "3mo"],
        "confound_flags": {},
    },
    {
        "arm": "Late",
        "contrast": "late_vs_relaxed_per_age",
        "ages": ["P1", "4W", "3mo"],
        "confound_flags": {
            "P1": "Late stress at P1 is Pool3 only; group fully confounded with pool.",
        },
    },
]

# Persistence-class color map (used across all plots)
CLASS_COLORS = {
    "persistent":         "#a50026",  # darkest red — paper-quality signal
    "resolving_early":    "#f46d43",
    "established_late":   "#fdae61",
    "P1_only":            "#fee090",
    "transient_4W":       "#abd9e9",
    "emergent_3mo":       "#4575b4",
    "P1_3mo_only":        "#74add1",
    "none":               "#cccccc",
}

PERSISTENCE_RANK = {
    "persistent": 0, "resolving_early": 1, "established_late": 2,
    "P1_only": 3, "transient_4W": 4, "emergent_3mo": 5,
    "P1_3mo_only": 6, "none": 7,
}


# ============================================================================
# CLI
# ============================================================================

def parse_args():
    p = argparse.ArgumentParser(description="Phase 8g: cross-age persistence")
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--padj-cutoff", type=float, default=0.05,
                   help="FDR cutoff for calling a gene DEG (default 0.05)")
    p.add_argument("--logfc-cutoff", type=float, default=0.5,
                   help="|log2FC| cutoff for DEG (default 0.5)")
    p.add_argument("--pathway-fdr-cutoff", type=float, default=0.1,
                   help="FDR cutoff for pathway/TF hits (default 0.1)")
    p.add_argument("--top-n-label", type=int, default=20,
                   help="Top N to label in trajectory plots (default 20)")
    p.add_argument("--top-n-plot", type=int, default=10,
                   help="Top N celltypes/categories per faceted plot (default 10)")
    p.add_argument("--sex", type=str, default="combined",
                   choices=["combined", "M", "F"],
                   help="Sex stratum (default combined; M/F persistence is "
                        "meaningless at n~2 and is not recommended)")
    p.add_argument("--amplify-ratio", type=float, default=1.2,
                   help="Trajectory shape: |effect_last|/|effect_first| >= this "
                        "=> amplifying (default 1.2)")
    p.add_argument("--attenuate-ratio", type=float, default=0.8,
                   help="Trajectory shape: ratio <= this => attenuating (default 0.8)")
    return p.parse_args()


# ============================================================================
# Loaders
# ============================================================================

def _slug(s: str) -> str:
    return str(s).replace(" ", "_").replace("/", "-").replace(".", "")


def _extract_age(group_level):
    """For the per-age stress contrasts 8g uses (early/late_vs_relaxed_per_age),
    group_level holds the age DIRECTLY (P1 / 4W / 3mo) — not an 'age-4W' encoding.
    Return as-is."""
    return group_level if isinstance(group_level, str) else None


def _filter_sex(sub, sex="combined"):
    """Pin to one sex stratum. M/F persistence on n~2/group is meaningless, so
    8g uses combined only by default. Mixing strata corrupts the per-age direction
    map in classify_dataframe (last row wins), so this filter is mandatory."""
    if sex is not None and "sex" in sub.columns:
        sub = sub[sub["sex"] == sex]
    return sub


def load_tables(cfg):
    base = Path(cfg["results_dir"]) / "tables"
    paths = {
        "de": base / "08b_de" / "08b_de_results.csv",
        "pw": base / "08c_pathways" / "08c_pathway_results.csv",
        "tf": base / "08c_pathways" / "08c_tf_activity.csv",
        "disruption": base / "08b_de" / "08b_developmental_disruption_genes.csv",
    }
    out = {}
    for kind, p in paths.items():
        if p.is_file():
            df = pd.read_csv(p, low_memory=False)
            print(f"  {kind}: {len(df):,} rows from {p.name}")
            out[kind] = df
        else:
            print(f"  [info] {kind} table missing: {p.name}")
            out[kind] = pd.DataFrame()
    return out


def prep_de(df, contrast, sex="combined"):
    """Filter 08b DE to one contrast + sex stratum, add 'age', keep 'level',
    drop missing genes. Level is RETAINED (not filtered) so persistence can be
    classified per (celltype, level, gene) — whole + every brain region."""
    if df.empty or "contrast" not in df.columns:
        return pd.DataFrame()
    sub = df[df["contrast"] == contrast].copy()
    if sub.empty:
        return pd.DataFrame()
    sub = _filter_sex(sub, sex)
    if sub.empty:
        return pd.DataFrame()
    sub = sub.dropna(subset=["gene", "stat"])
    sub["age"] = sub["group_level"].map(_extract_age)
    if "level" not in sub.columns:
        sub["level"] = "whole"
    return sub.dropna(subset=["age"])[
        ["celltype", "level", "gene", "stat", "padj", "log2FC", "age"]
    ]


def prep_pw(df, contrast, sex="combined"):
    """Filter 08c pathway. 08c writes pathway name as 'source', FDR as 'FDR'.
    Sex-filtered; 'level' retained for per-region persistence."""
    if df.empty or "contrast" not in df.columns:
        return pd.DataFrame()
    sub = df[df["contrast"] == contrast].copy()
    if sub.empty:
        return pd.DataFrame()
    sub = _filter_sex(sub, sex)
    if sub.empty:
        return pd.DataFrame()
    sub["age"] = sub["group_level"].map(_extract_age)
    if "level" not in sub.columns:
        sub["level"] = "whole"
    rename = {}
    if "source" in sub.columns and "pathway" not in sub.columns:
        rename["source"] = "pathway"
    if "FDR" in sub.columns and "padj" not in sub.columns:
        rename["FDR"] = "padj"
    if rename:
        sub = sub.rename(columns=rename)
    return sub.dropna(subset=["age"])


def prep_tf(df, contrast, sex="combined"):
    """Filter 08c TF activity. Columns: contrast, group_level, celltype, TF,
    activity_score, pvalue, FDR, direction. Sex-filtered; 'level' retained."""
    if df.empty or "contrast" not in df.columns:
        return pd.DataFrame()
    sub = df[df["contrast"] == contrast].copy()
    if sub.empty:
        return pd.DataFrame()
    sub = _filter_sex(sub, sex)
    if sub.empty:
        return pd.DataFrame()
    sub["age"] = sub["group_level"].map(_extract_age)
    if "level" not in sub.columns:
        sub["level"] = "whole"
    if "FDR" in sub.columns and "padj" not in sub.columns:
        sub = sub.rename(columns={"FDR": "padj"})
    return sub.dropna(subset=["age"])


# ============================================================================
# Persistence classification (shared logic for gene, pathway, TF)
# ============================================================================

def classify_persistence(per_age_status):
    """Given dict {age: 'up' | 'down' | 'none'}, return persistence class.

    'persistent' requires same-direction DE at all three ages.
    Direction switches → classify by the multi-age presence pattern without
    direction sensitivity, since direction-flipping isn't 'persistent' by any
    reasonable definition.
    """
    p1, w4, m3 = per_age_status["P1"], per_age_status["4W"], per_age_status["3mo"]
    sig_p1 = p1 != "none"
    sig_4w = w4 != "none"
    sig_3m = m3 != "none"
    n_sig = sum([sig_p1, sig_4w, sig_3m])

    if n_sig == 0:
        return "none"
    if n_sig == 3:
        # All three ages — require same direction for "persistent"
        if p1 == w4 == m3:
            return "persistent"
        # Direction inconsistent — fall through to "P1_3mo_only" style
        # (rare; still informative). Classify by majority age subset structure.
        return "persistent"  # still all 3 ages; flag direction in a separate col
    if n_sig == 2:
        if sig_p1 and sig_4w:
            return "resolving_early"
        if sig_4w and sig_3m:
            return "established_late"
        if sig_p1 and sig_3m:
            return "P1_3mo_only"
    # n_sig == 1
    if sig_p1: return "P1_only"
    if sig_4w: return "transient_4W"
    if sig_3m: return "emergent_3mo"
    return "none"


def classify_dataframe(per_age_df, age_col="age", direction_col="direction",
                       feature_cols=("celltype", "feature")):
    """Take long-format per-age status and return classified rows.

    per_age_df schema: feature_cols + [age, direction] where direction is
    'up' / 'down' / 'none'.
    Output: feature_cols + [persistence_class, direction_consistent,
        P1_dir, 4W_dir, 3mo_dir, n_sig_ages].
    """
    rows = []
    for keys, grp in per_age_df.groupby(list(feature_cols), observed=True):
        per_age = {a: "none" for a in ("P1", "4W", "3mo")}
        for _, r in grp.iterrows():
            per_age[r[age_col]] = r[direction_col]
        klass = classify_persistence(per_age)
        n_sig = sum(1 for v in per_age.values() if v != "none")
        sig_dirs = [v for v in per_age.values() if v != "none"]
        direction_consistent = len(set(sig_dirs)) <= 1
        # If three sig ages but direction inconsistent, downgrade label
        if klass == "persistent" and not direction_consistent:
            klass = "persistent_directionswap"
        row = dict(zip(feature_cols, keys if isinstance(keys, tuple) else (keys,)))
        row.update({
            "persistence_class": klass,
            "direction_consistent": direction_consistent,
            "P1_dir": per_age["P1"],
            "4W_dir": per_age["4W"],
            "3mo_dir": per_age["3mo"],
            "n_sig_ages": n_sig,
        })
        rows.append(row)
    return pd.DataFrame(rows)


# ============================================================================
# View 1: Gene-level persistence
# ============================================================================

def run_view1_gene_persistence(de_df, args):
    """Classify each (celltype, gene, arm) into persistence class."""
    print("\n[View 1] Gene-level persistence")
    if de_df.empty:
        print("  [skip] no DE table")
        return pd.DataFrame()
    all_rows = []
    for arm in ARMS:
        sub = prep_de(de_df, arm["contrast"], sex=args.sex)
        if sub.empty:
            print(f"  [skip] {arm['arm']} arm: no DE rows")
            continue
        sub = sub[sub["age"].isin(arm["ages"])].copy()
        # Per (celltype, level, gene, age) — direction column
        sub["direction"] = "none"
        sig = (sub["padj"] < args.padj_cutoff) & \
              (sub["log2FC"].abs() > args.logfc_cutoff)
        sub.loc[sig & (sub["log2FC"] > 0), "direction"] = "up"
        sub.loc[sig & (sub["log2FC"] < 0), "direction"] = "down"

        # Keep only (celltype, level, gene) that are DE at >=1 age.
        any_sig = sub.groupby(["celltype", "level", "gene"])["direction"].apply(
            lambda x: (x != "none").any())
        keep_keys = any_sig[any_sig].index
        idx = sub.set_index(["celltype", "level", "gene"]).index
        sub = sub.set_index(["celltype", "level", "gene"]).loc[
            idx.isin(keep_keys)].reset_index()

        if sub.empty:
            continue

        classified = classify_dataframe(
            sub, age_col="age", direction_col="direction",
            feature_cols=("celltype", "level", "gene"))
        classified["arm"] = arm["arm"]
        classified["confound_note"] = arm["confound_flags"].get("P1", "") \
            if arm["arm"] == "Late" else ""
        # Attach the actual log2FC per age (useful in output)
        wide = (sub.pivot_table(index=["celltype", "level", "gene"],
                                 columns="age", values="log2FC",
                                 aggfunc="first")
                 .reset_index())
        wide.columns = ["celltype", "level", "gene"] + [f"{c}_log2FC" for c in wide.columns[3:]]
        classified = classified.merge(wide, on=["celltype", "level", "gene"], how="left")
        all_rows.append(classified)

    if not all_rows:
        return pd.DataFrame()
    df = pd.concat(all_rows, ignore_index=True)
    print(f"  Total (gene, celltype, arm) rows: {len(df):,}")
    cls_counts = df.groupby(["arm", "persistence_class"]).size().unstack(fill_value=0)
    print(f"  Class counts:\n{cls_counts.to_string()}")
    return df


# ============================================================================
# View 2: Pathway-level persistence
# ============================================================================

def run_view2_pathway_persistence(pw_df, args):
    print("\n[View 2] Pathway-level persistence")
    if pw_df.empty:
        print("  [skip] no pathway table")
        return pd.DataFrame()
    all_rows = []
    for arm in ARMS:
        sub = prep_pw(pw_df, arm["contrast"], sex=args.sex)
        if sub.empty:
            continue
        sub = sub[sub["age"].isin(arm["ages"])].copy()
        sub["direction"] = "none"
        sig = sub["padj"] < args.pathway_fdr_cutoff
        sub.loc[sig & (sub["NES"] > 0), "direction"] = "up"
        sub.loc[sig & (sub["NES"] < 0), "direction"] = "down"

        any_sig = sub.groupby(["celltype", "level", "pathway"])["direction"].apply(
            lambda x: (x != "none").any())
        keep = any_sig[any_sig].index
        idx = sub.set_index(["celltype", "level", "pathway"]).index
        sub = sub.set_index(["celltype", "level", "pathway"]).loc[
            idx.isin(keep)].reset_index()
        if sub.empty:
            continue

        classified = classify_dataframe(
            sub, age_col="age", direction_col="direction",
            feature_cols=("celltype", "level", "pathway"))
        classified["arm"] = arm["arm"]
        classified["confound_note"] = arm["confound_flags"].get("P1", "") \
            if arm["arm"] == "Late" else ""
        wide = (sub.pivot_table(index=["celltype", "level", "pathway"],
                                 columns="age", values="NES", aggfunc="first")
                 .reset_index())
        wide.columns = ["celltype", "level", "pathway"] + [f"{c}_NES" for c in wide.columns[3:]]
        # Also carry collection column if present
        if "collection" in sub.columns:
            coll = sub.groupby(["celltype", "level", "pathway"])["collection"].first().reset_index()
            classified = classified.merge(coll, on=["celltype", "level", "pathway"], how="left")
        classified = classified.merge(wide, on=["celltype", "level", "pathway"], how="left")
        all_rows.append(classified)

    if not all_rows:
        return pd.DataFrame()
    df = pd.concat(all_rows, ignore_index=True)
    print(f"  Total (pathway, celltype, arm) rows: {len(df):,}")
    cls_counts = df.groupby(["arm", "persistence_class"]).size().unstack(fill_value=0)
    print(f"  Class counts:\n{cls_counts.to_string()}")
    return df


# ============================================================================
# View 3: TF-level persistence
# ============================================================================

def run_view3_tf_persistence(tf_df, args):
    print("\n[View 3] TF-level persistence")
    if tf_df.empty:
        print("  [skip] no TF activity table "
              "(run 8c with --tf to enable)")
        return pd.DataFrame()
    all_rows = []
    for arm in ARMS:
        sub = prep_tf(tf_df, arm["contrast"], sex=args.sex)
        if sub.empty:
            continue
        sub = sub[sub["age"].isin(arm["ages"])].copy()
        sub["direction"] = "none"
        sig = sub["padj"] < args.pathway_fdr_cutoff
        sub.loc[sig & (sub["activity_score"] > 0), "direction"] = "up"
        sub.loc[sig & (sub["activity_score"] < 0), "direction"] = "down"

        any_sig = sub.groupby(["celltype", "level", "TF"])["direction"].apply(
            lambda x: (x != "none").any())
        keep = any_sig[any_sig].index
        idx = sub.set_index(["celltype", "level", "TF"]).index
        sub = sub.set_index(["celltype", "level", "TF"]).loc[
            idx.isin(keep)].reset_index()
        if sub.empty:
            continue

        classified = classify_dataframe(
            sub, age_col="age", direction_col="direction",
            feature_cols=("celltype", "level", "TF"))
        classified["arm"] = arm["arm"]
        classified["confound_note"] = arm["confound_flags"].get("P1", "") \
            if arm["arm"] == "Late" else ""
        wide = (sub.pivot_table(index=["celltype", "level", "TF"],
                                 columns="age", values="activity_score",
                                 aggfunc="first").reset_index())
        wide.columns = ["celltype", "level", "TF"] + [f"{c}_activity" for c in wide.columns[3:]]
        classified = classified.merge(wide, on=["celltype", "level", "TF"], how="left")
        all_rows.append(classified)

    if not all_rows:
        return pd.DataFrame()
    df = pd.concat(all_rows, ignore_index=True)
    print(f"  Total (TF, celltype, arm) rows: {len(df):,}")
    cls_counts = df.groupby(["arm", "persistence_class"]).size().unstack(fill_value=0)
    print(f"  Class counts:\n{cls_counts.to_string()}")
    return df


# ============================================================================
# View 5: Early vs Late at each age (overlap + rank-rank)
# ============================================================================

def run_view5_early_vs_late(de_df, args):
    """Per age × celltype: hypergeometric overlap of ES-vs-Rel DEGs and
    LS-vs-Rel DEGs, plus full-list Spearman ρ for signature concordance."""
    from statsmodels.stats.multitest import multipletests
    print("\n[View 5] Early vs Late at each age")
    if de_df.empty:
        return pd.DataFrame()

    es = prep_de(de_df, "early_vs_relaxed_per_age", sex=args.sex)
    ls = prep_de(de_df, "late_vs_relaxed_per_age", sex=args.sex)
    if es.empty or ls.empty:
        print("  [skip] missing ES or LS DE rows")
        return pd.DataFrame()

    levels = sorted(set(es["level"]).union(ls["level"]))
    rows = []
    for lvl in levels:
        es_l = es[es["level"] == lvl]
        ls_l = ls[ls["level"] == lvl]
        for age in ("P1", "4W", "3mo"):
            es_age = es_l[es_l["age"] == age]
            ls_age = ls_l[ls_l["age"] == age]
            if es_age.empty or ls_age.empty:
                continue

            common_cts = set(es_age["celltype"]) & set(ls_age["celltype"])
            for ct in sorted(common_cts):
                es_ct = es_age[es_age["celltype"] == ct]
                ls_ct = ls_age[ls_age["celltype"] == ct]
                universe = set(es_ct["gene"]) & set(ls_ct["gene"])
                if len(universe) < 100:
                    continue
                n = len(universe)

                es_sig = es_ct[(es_ct["padj"] < args.padj_cutoff) &
                               (es_ct["log2FC"].abs() > args.logfc_cutoff)]
                ls_sig = ls_ct[(ls_ct["padj"] < args.padj_cutoff) &
                               (ls_ct["log2FC"].abs() > args.logfc_cutoff)]

                es_up = set(es_sig.loc[es_sig["log2FC"] > 0, "gene"]) & universe
                es_dn = set(es_sig.loc[es_sig["log2FC"] < 0, "gene"]) & universe
                ls_up = set(ls_sig.loc[ls_sig["log2FC"] > 0, "gene"]) & universe
                ls_dn = set(ls_sig.loc[ls_sig["log2FC"] < 0, "gene"]) & universe

                for direction, (a_set, b_set) in [
                    ("concordant_up", (es_up, ls_up)),
                    ("concordant_down", (es_dn, ls_dn)),
                    ("discordant_es_up_ls_dn", (es_up, ls_dn)),
                    ("any_overlap", (es_up | es_dn, ls_up | ls_dn)),
                ]:
                    overlap = a_set & b_set
                    if len(a_set) == 0 or len(b_set) == 0 or len(overlap) == 0:
                        p = 1.0
                    else:
                        p = float(hypergeom.sf(len(overlap) - 1, n,
                                                len(a_set), len(b_set)))

                    # Full-list spearman ρ (signed Wald stats)
                    es_stats = es_ct.set_index("gene")["stat"]
                    ls_stats = ls_ct.set_index("gene")["stat"]
                    comm = es_stats.index.intersection(ls_stats.index)
                    rho, rho_p = (float("nan"), float("nan"))
                    if len(comm) >= 50:
                        rho_val, rho_pv = spearmanr(es_stats.loc[comm], ls_stats.loc[comm])
                        rho, rho_p = float(rho_val), float(rho_pv)

                    rows.append({
                        "level": lvl, "age": age, "celltype": ct,
                        "direction": direction,
                        "n_overlap": len(overlap),
                        "n_ES": len(a_set), "n_LS": len(b_set),
                        "n_universe": n,
                        "pvalue": p,
                        "spearman_rho_full_stats": rho,
                        "spearman_p_full_stats": rho_p,
                        "overlap_genes": ";".join(sorted(overlap)[:50]),
                    })

    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    # BH-FDR within (level, age, direction)
    fdr_parts = []
    for _, g in df.groupby(["level", "age", "direction"]):
        _, fdr, _, _ = multipletests(g["pvalue"].fillna(1.0), method="fdr_bh")
        fdr_parts.append(pd.Series(fdr, index=g.index))
    df["fdr"] = pd.concat(fdr_parts)
    df["neg_log10_p"] = -np.log10(df["pvalue"].clip(lower=1e-300))
    print(f"  Total ES vs LS overlap tests: {len(df):,}; "
          f"FDR<0.05: {(df['fdr'] < 0.05).sum():,}")
    return df


# ============================================================================
# View 6: Cross-arm core signature (intersection of persistent in both arms)
# ============================================================================

def run_view6_core_signature(gene_persist, pw_persist):
    """Cross-arm intersection of 'persistent' classification.

    A gene/pathway in this table is one that's deregulated at all three ages
    in BOTH the Early arm AND the Late arm in the same direction. This is the
    most robust 'core stress signature' — independent of stress timing,
    persistent across development.

    The publication-worthy distillation.
    """
    print("\n[View 6] Cross-arm core signature (intersection of persistent calls)")
    results = {}

    for kind, df, feature_col in (
        ("genes", gene_persist, "gene"),
        ("pathways", pw_persist, "pathway"),
    ):
        if df.empty:
            results[kind] = pd.DataFrame()
            continue
        early = df[(df["arm"] == "Early") &
                   (df["persistence_class"] == "persistent")]
        late = df[(df["arm"] == "Late") &
                  (df["persistence_class"] == "persistent")]
        if early.empty or late.empty:
            print(f"  {kind}: no persistent calls in one or both arms (empty)")
            results[kind] = pd.DataFrame()
            continue

        # Join on (celltype, level, feature) and require same direction in both arms
        key_cols = ["celltype", "level", feature_col]
        early_keyed = early[key_cols + ["P1_dir", "4W_dir", "3mo_dir"]].rename(
            columns={c: f"Early_{c}" for c in ("P1_dir", "4W_dir", "3mo_dir")})
        late_keyed = late[key_cols + ["P1_dir", "4W_dir", "3mo_dir"]].rename(
            columns={c: f"Late_{c}" for c in ("P1_dir", "4W_dir", "3mo_dir")})
        merged = early_keyed.merge(late_keyed, on=key_cols, how="inner")
        if merged.empty:
            print(f"  {kind}: 0 shared persistent features")
            results[kind] = pd.DataFrame()
            continue
        # Same direction in both arms (P1 dir as canonical, or all three)
        merged["both_arms_same_direction"] = (
            (merged["Early_P1_dir"] == merged["Late_P1_dir"]) &
            (merged["Early_4W_dir"] == merged["Late_4W_dir"]) &
            (merged["Early_3mo_dir"] == merged["Late_3mo_dir"])
        )
        merged["core_signature"] = merged["both_arms_same_direction"]

        # Attach effect-size columns from the source tables
        es_col = "log2FC" if kind == "genes" else "NES"
        for arm_name, arm_df in (("Early", early), ("Late", late)):
            cols = key_cols + [f"P1_{es_col}", f"4W_{es_col}", f"3mo_{es_col}"]
            cols = [c for c in cols if c in arm_df.columns]
            arm_es = arm_df[cols].rename(columns={
                c: f"{arm_name}_{c}" for c in cols if c not in key_cols
            })
            merged = merged.merge(arm_es, on=key_cols, how="left")

        merged["P1_confound_note"] = (
            "Late P1 pool-confounded; core signature robust if P1 dir matches "
            "downstream ages."
        )
        n_core = merged["core_signature"].sum()
        print(f"  {kind}: {len(merged):,} shared-persistent features; "
              f"{n_core:,} are full-direction-concordant core signature")
        results[kind] = merged

    return results.get("genes", pd.DataFrame()), results.get("pathways", pd.DataFrame())


# ============================================================================
# Plots
# ============================================================================

def _save_fig(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot: {path.name}")


def plot_persistence_class_barplot(persist_df, kind_label, plots_root):
    """Stacked barchart of persistence classes per celltype, faceted by arm."""
    if persist_df.empty:
        return
    arms = sorted(persist_df["arm"].unique())
    fig, axes = plt.subplots(1, len(arms), figsize=(8 * len(arms), 6),
                             sharey=True, constrained_layout=True)
    if len(arms) == 1:
        axes = [axes]
    classes_present = [c for c in PERSISTENCE_RANK.keys()
                       if c in persist_df["persistence_class"].unique()]
    for ax, arm in zip(axes, arms):
        sub = persist_df[persist_df["arm"] == arm]
        counts = (sub.groupby(["celltype", "persistence_class"])
                  .size().unstack(fill_value=0))
        counts = counts.reindex(columns=classes_present, fill_value=0)
        # Sort celltypes by total non-'none' count (most-affected first)
        order = (counts.drop(columns=["none"], errors="ignore")
                 .sum(axis=1).sort_values(ascending=False).index)
        counts = counts.loc[order]
        if counts.empty:
            continue
        bottom = np.zeros(len(counts))
        for klass in classes_present:
            vals = counts[klass].values
            ax.barh(range(len(counts)), vals, left=bottom,
                    color=CLASS_COLORS.get(klass, "#999999"),
                    label=klass)
            bottom = bottom + vals
        ax.set_yticks(range(len(counts)))
        ax.set_yticklabels(counts.index, fontsize=8)
        ax.set_xlabel(f"# {kind_label}")
        ax.set_title(f"{arm} arm", fontsize=10)
        ax.invert_yaxis()
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower right", fontsize=7, ncol=1)
    fig.suptitle(f"{kind_label} persistence classes per cell type", fontsize=12)
    _save_fig(fig, plots_root / f"{kind_label.replace(' ', '_').lower()}_persistence_class_barplot.png")


def plot_persistence_heatmap(persist_df, feature_col, kind_label, plots_root,
                              top_n_label=50):
    """Heatmap: top features × celltypes, colored by persistence class.
    One per arm."""
    if persist_df.empty:
        return
    for arm in sorted(persist_df["arm"].unique()):
        sub = persist_df[persist_df["arm"] == arm].copy()
        # Prioritize persistent > resolving > established > others
        sub["sort_rank"] = sub["persistence_class"].map(PERSISTENCE_RANK).fillna(99)
        # Pick features by sort_rank (persistent first), then by celltype coverage
        feat_score = (sub.groupby(feature_col)
                      .agg(min_rank=("sort_rank", "min"),
                           n_cts=("celltype", "nunique")))
        feat_score = feat_score.sort_values(["min_rank", "n_cts"],
                                             ascending=[True, False])
        top_feats = feat_score.head(top_n_label).index
        if len(top_feats) == 0:
            continue
        plot_df = sub[sub[feature_col].isin(top_feats)]
        pivot = (plot_df.pivot_table(index=feature_col, columns="celltype",
                                      values="sort_rank", aggfunc="min")
                 .reindex(top_feats))
        if pivot.empty:
            continue
        # Build a discrete colormap by mapping rank → color
        from matplotlib.colors import ListedColormap, BoundaryNorm
        classes_ordered = list(PERSISTENCE_RANK.keys())
        colors = [CLASS_COLORS[c] for c in classes_ordered]
        cmap = ListedColormap(colors)
        bounds = list(range(len(classes_ordered) + 1))
        norm = BoundaryNorm(bounds, cmap.N)

        fig, ax = plt.subplots(
            figsize=(max(6, pivot.shape[1] * 0.55 + 2),
                     max(8, pivot.shape[0] * 0.28 + 2)))
        plot_mat = pivot.fillna(PERSISTENCE_RANK["none"]).values
        im = ax.imshow(plot_mat, aspect="auto", cmap=cmap, norm=norm,
                       interpolation="none")
        ax.set_xticks(range(pivot.shape[1]))
        ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(pivot.shape[0]))
        ax.set_yticklabels(pivot.index, fontsize=7)
        ax.set_title(
            f"{kind_label} persistence classes: top {len(top_feats)} features\n"
            f"{arm} arm  (rows sorted by best class then celltype coverage)",
            fontsize=9)
        from matplotlib.patches import Patch
        legend_handles = [Patch(color=CLASS_COLORS[c], label=c)
                          for c in classes_ordered]
        ax.legend(handles=legend_handles, fontsize=6, loc="center left",
                  bbox_to_anchor=(1.02, 0.5))
        _save_fig(fig, plots_root / f"{kind_label.replace(' ', '_').lower()}_persistence_heatmap_{arm}.png")


def plot_effect_size_trajectory(persist_df, feature_col, value_prefix,
                                kind_label, plots_root, top_n_label=20):
    """Trajectory: P1 → 4W → 3mo effect size for top persistent features.
    Faceted by celltype × arm."""
    if persist_df.empty:
        return
    # Pick only persistent + resolving + established classes (the multi-age ones)
    multi_classes = {"persistent", "resolving_early", "established_late",
                     "P1_3mo_only"}
    multi = persist_df[persist_df["persistence_class"].isin(multi_classes)].copy()
    if multi.empty:
        return

    age_cols = [f"P1_{value_prefix}", f"4W_{value_prefix}", f"3mo_{value_prefix}"]
    if not all(c in multi.columns for c in age_cols):
        return

    for arm in sorted(multi["arm"].unique()):
        arm_df = multi[multi["arm"] == arm]
        celltypes = arm_df["celltype"].value_counts().head(8).index
        if len(celltypes) == 0:
            continue
        ncols = min(3, len(celltypes))
        nrows = int(np.ceil(len(celltypes) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows),
                                 constrained_layout=True)
        axes_flat = np.array(axes).flatten() if len(celltypes) > 1 else [axes]

        for ax, ct in zip(axes_flat, celltypes):
            ct_df = arm_df[arm_df["celltype"] == ct].copy()
            ct_df["max_abs"] = ct_df[age_cols].abs().max(axis=1)
            top = ct_df.nlargest(top_n_label, "max_abs")
            if top.empty:
                ax.set_visible(False)
                continue
            ages = ["P1", "4W", "3mo"]
            for _, r in top.iterrows():
                vals = [r[c] if not pd.isna(r[c]) else 0 for c in age_cols]
                color = CLASS_COLORS.get(r["persistence_class"], "gray")
                ax.plot(ages, vals, color=color, alpha=0.55, lw=1.2, marker="o",
                        markersize=3)
                # Label endpoint
                end_x = ages[-1] if abs(vals[-1]) >= abs(vals[0]) else ages[0]
                end_y = vals[-1] if abs(vals[-1]) >= abs(vals[0]) else vals[0]
                ax.annotate(str(r[feature_col])[:30], (end_x, end_y),
                            fontsize=5.5, alpha=0.7,
                            xytext=(3, 0), textcoords="offset points")
            ax.axhline(0, color="k", lw=0.5, alpha=0.4)
            ax.set_title(ct, fontsize=8)
            ax.set_ylabel(f"{value_prefix}", fontsize=7)

        for ax in axes_flat[len(celltypes):]:
            ax.set_visible(False)
        fig.suptitle(
            f"{kind_label} effect-size trajectory across ages — {arm} arm\n"
            f"(top {top_n_label} per cell type by max |{value_prefix}|; "
            f"lines coloured by persistence class)", fontsize=10)
        _save_fig(fig,
                  plots_root / f"{kind_label.replace(' ', '_').lower()}_trajectory_{arm}.png")


def plot_early_vs_late_overlap(ovr_df, plots_root, top_n_label=20):
    """Venn-ish bar plot of overlap per age × celltype, plus rank-rank scatter
    summary table per age."""
    if ovr_df.empty:
        return
    # 1. Bar chart: -log10(p) for any_overlap per age × celltype
    for age in sorted(ovr_df["age"].unique()):
        sub = ovr_df[(ovr_df["age"] == age) &
                     (ovr_df["direction"] == "any_overlap")]
        if sub.empty:
            continue
        sub = sub.sort_values("neg_log10_p", ascending=False).head(top_n_label)
        colors = ["#a50026" if f < 0.05 else "#fdae61"
                  for f in sub["fdr"]]
        fig, ax = plt.subplots(figsize=(8, max(4, len(sub) * 0.35 + 1)))
        ax.barh(range(len(sub)), sub["neg_log10_p"], color=colors, alpha=0.85)
        ax.set_yticks(range(len(sub)))
        labels = [f"{r['celltype']}  "
                  f"(n_overlap={r['n_overlap']}, ES={r['n_ES']}, LS={r['n_LS']})"
                  for _, r in sub.iterrows()]
        ax.set_yticklabels(labels, fontsize=7)
        ax.axvline(-np.log10(0.05), color="k", lw=0.6, ls="--")
        ax.set_xlabel("−log10(p) of ES ∩ LS overlap")
        ax.set_title(
            f"ES-vs-Rel ∩ LS-vs-Rel DEG overlap at {age}\n"
            f"(red = FDR<0.05; any-direction overlap)", fontsize=9)
        _save_fig(fig, plots_root / f"overlap_es_vs_ls_{age}.png")

    # 2. Spearman ρ summary heatmap: age × celltype
    rho_df = ovr_df.drop_duplicates(subset=["age", "celltype"])
    if not rho_df.empty and "spearman_rho_full_stats" in rho_df.columns:
        pivot = rho_df.pivot(index="celltype", columns="age",
                             values="spearman_rho_full_stats")
        if not pivot.empty:
            max_abs = max(np.abs(pivot.fillna(0).values).max(), 0.3)
            fig, ax = plt.subplots(figsize=(max(5, pivot.shape[1] + 2),
                                            max(4, pivot.shape[0] * 0.5 + 2)))
            im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r",
                           vmin=-max_abs, vmax=max_abs)
            plt.colorbar(im, ax=ax,
                         label="Spearman ρ of ES vs LS Wald stats")
            ax.set_xticks(range(pivot.shape[1]))
            ax.set_xticklabels(pivot.columns, fontsize=8)
            ax.set_yticks(range(pivot.shape[0]))
            ax.set_yticklabels(pivot.index, fontsize=8)
            ax.set_title("ES vs LS signature concordance across ages\n"
                         "(positive ρ = same direction; full-list Spearman)",
                         fontsize=9)
            for i in range(pivot.shape[0]):
                for j in range(pivot.shape[1]):
                    v = pivot.values[i, j]
                    if not np.isnan(v):
                        ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                                fontsize=7,
                                color="white" if abs(v) > max_abs * 0.6 else "black")
            _save_fig(fig, plots_root / "es_vs_ls_signature_concordance_heatmap.png")


def plot_core_signature(core_genes, core_pw, plots_root, top_n_label=30):
    """Tables-as-plots for the publication-quality cross-arm core signature."""
    for label, df, value_prefix in (
        ("genes", core_genes, "log2FC"),
        ("pathways", core_pw, "NES"),
    ):
        if df.empty:
            continue
        core = df[df.get("core_signature", False)] if "core_signature" in df.columns else df
        if core.empty:
            continue
        # Effect-size product summed across the 6 (arm × age) cells
        es_cols = [c for c in core.columns if c.endswith(f"_{value_prefix}")]
        if not es_cols:
            continue
        core = core.copy()
        core["mean_abs_effect"] = core[es_cols].abs().mean(axis=1)
        top = core.nlargest(top_n_label, "mean_abs_effect").copy()
        feature_col = "gene" if label == "genes" else "pathway"
        top["row_label"] = top["celltype"] + " | " + top[feature_col].astype(str)

        # Heatmap of effect sizes
        plot_mat = top[es_cols]
        max_abs = max(float(np.abs(plot_mat.fillna(0).values).max()), 1e-6)
        fig, ax = plt.subplots(figsize=(max(7, len(es_cols) * 1.2 + 2),
                                        max(6, len(top) * 0.28 + 2)))
        im = ax.imshow(plot_mat.values, aspect="auto", cmap="RdBu_r",
                       vmin=-max_abs, vmax=max_abs)
        plt.colorbar(im, ax=ax, label=f"{value_prefix}")
        ax.set_xticks(range(len(es_cols)))
        ax.set_xticklabels([c.replace(f"_{value_prefix}", "") for c in es_cols],
                           rotation=30, ha="right", fontsize=7)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top["row_label"].values, fontsize=7)
        ax.set_title(
            f"Cross-arm core signature: top {len(top)} {label}\n"
            f"(persistent in BOTH Early and Late arms, same direction throughout)",
            fontsize=9)
        _save_fig(fig, plots_root / f"core_signature_{label}_heatmap.png")


# ============================================================================
# B: Trajectory shape (amplifying / attenuating / stable)
# ============================================================================

def add_trajectory_shape(df, effect_prefix, args):
    """Add a 'trajectory_shape' column to a persistence table.

    For features significant at >=2 ages with CONSISTENT direction, compare the
    |effect| at the first vs last significant age (age order P1 < 4W < 3mo):
        amplifying  ratio >= amplify_ratio (effect grows across development)
        attenuating ratio <= attenuate_ratio (effect shrinks / resolving)
        stable      otherwise
    Others: 'single_age' (<2 sig ages), 'directionswap', 'undetermined'.
    """
    if df.empty:
        return df
    ages = ["P1", "4W", "3mo"]
    eff = {a: f"{a}_{effect_prefix}" for a in ages}
    for a in ages:
        if eff[a] not in df.columns:
            df[eff[a]] = np.nan
    shapes = []
    for _, r in df.iterrows():
        sig_ages = [a for a in ages if str(r.get(f"{a}_dir", "none")) != "none"]
        if len(sig_ages) < 2:
            shapes.append("single_age"); continue
        if len({r[f"{a}_dir"] for a in sig_ages}) > 1:
            shapes.append("directionswap"); continue
        e0 = abs(float(r.get(eff[sig_ages[0]], np.nan)))
        e1 = abs(float(r.get(eff[sig_ages[-1]], np.nan)))
        if not (np.isfinite(e0) and np.isfinite(e1)) or e0 == 0:
            shapes.append("undetermined"); continue
        ratio = e1 / e0
        if ratio >= args.amplify_ratio:
            shapes.append("amplifying")
        elif ratio <= args.attenuate_ratio:
            shapes.append("attenuating")
        else:
            shapes.append("stable")
    df = df.copy()
    df["trajectory_shape"] = shapes
    return df


# ============================================================================
# C: Persistence x developmental-disruption cross-reference
# ============================================================================

def run_persistence_x_disruption(gene_df, disruption_df, args):
    """Cross-reference 8g stress-persistence classes against the 8b developmental
    disruption classes (LOST / GAINED / universal / ...).

    NOTE: links two DIFFERENT contrasts — 8g uses the per-age stress contrasts
    (early/late_vs_relaxed_per_age); 8b disruption uses within_group_across_age
    (developmental trajectory). The interesting intersection is genes that are
    BOTH persistently stress-DE (8g 'persistent') AND developmentally LOST under
    stress (8b 'relaxed_only'). Matched on (celltype, level, gene), sex=combined.
    """
    print("\n[View C] Persistence x developmental disruption")
    if gene_df.empty or disruption_df.empty:
        print("  [skip] need both gene persistence and 08b disruption table")
        return pd.DataFrame()

    d = disruption_df.copy()
    if "sex" in d.columns:
        d = d[d["sex"] == "combined"]
    # 08b disruption stores the LOST/GAINED class in the 'direction' column
    # (values: universal / relaxed_only=LOST / stress_shared=GAINED / early_only /
    # late_only). Note this 'direction' is the disruption CLASS, not up/down.
    class_col = next((c for c in ["direction", "direction_class", "disruption_class",
                                   "class", "category"] if c in d.columns), None)
    if class_col is None:
        raise ValueError(
            "08b disruption table has no recognizable class column "
            f"(looked for direction/direction_class/...); got {list(d.columns)}. "
            "Cannot run persistence x disruption — fix the column mapping rather "
            "than skipping silently.")
    expected = {"universal", "relaxed_only", "stress_shared", "early_only", "late_only"}
    if not (set(d[class_col].dropna().unique()) & expected):
        raise ValueError(
            f"08b disruption '{class_col}' values {sorted(set(d[class_col].dropna().unique()))[:6]} "
            f"don't match expected disruption classes {sorted(expected)}. Refusing to "
            "proceed on a mismatched schema.")
    if "level" not in d.columns:
        d["level"] = "whole"
    keep = [c for c in ["celltype", "level", "gene", class_col] if c in d.columns]
    d = d[keep].rename(columns={class_col: "disruption_class"})

    on = [c for c in ["celltype", "level", "gene"] if c in d.columns]
    g = gene_df[["arm", "celltype", "level", "gene",
                 "persistence_class", "direction_consistent"]]
    merged = g.merge(d, on=on, how="inner")
    if merged.empty:
        print("  [warn] no (celltype, level, gene) overlap between 8g and 8b disruption")
        return pd.DataFrame()

    merged["level_flag"] = np.where(merged["level"] == "whole",
                                    "robust", "regional_exploratory")
    # Headline subset: persistent stress effect AND developmentally LOST
    merged["persistent_and_LOST"] = (
        (merged["persistence_class"] == "persistent") &
        (merged["disruption_class"] == "relaxed_only")
    )
    n_pl = int(merged["persistent_and_LOST"].sum())
    print(f"  {len(merged):,} matched (gene,celltype,level,arm) rows; "
          f"{n_pl:,} persistent-AND-developmentally-LOST")
    ct = (merged.groupby(["arm", "persistence_class", "disruption_class"])
          .size().reset_index(name="n"))
    print("  Top persistence x disruption cells:")
    print(ct.sort_values("n", ascending=False).head(8).to_string(index=False))
    return merged


# ============================================================================
# A: View 7 — 8f-complementary focal module
# ============================================================================

# Broad regulator families (regex on TF symbol) — immune / IEG-AP1 / stress-GR
FOCAL_TF_FAMILIES = ["IRF", "STAT", "NFKB", "REL", "JUN", "FOS", "FOSL",
                     "EGR", "ATF", "CEBP", "SPI", "RUNX", "NFIL3", "BCL3", "NR3C"]
# Immune/inflammatory pathway keywords for the leading-edge deep-dive
IMMUNE_PATHWAY_KEYWORDS = ["INTERFERON", "IFN", "COMPLEMENT", "INFLAMMAT",
                           "CYTOKINE", "TNF", "NFKB", "IL2_STAT", "IL6", "JAK_STAT",
                           "INNATE", "ALLOGRAFT", "CHEMOKINE", "ANTIGEN", "ISG"]


def _load_8f_concordant(cfg, kind):
    """Return the 8f concordance table (pathway or TF) if it exists, else empty.
    kind in {'pathway','tf'}."""
    fname = ("08f_pathway_concordance.csv" if kind == "pathway"
             else "08f_tf_concordance.csv")
    p = Path(cfg["results_dir"]) / "tables" / "08f_cross_tissue" / fname
    if not p.is_file():
        print(f"  [info] 8f table not found ({fname}); cross-tissue flag skipped")
        return pd.DataFrame()
    df = pd.read_csv(p, low_memory=False)
    return df[df["concordance_class"].astype(str).str.startswith("concordant")]


def run_view7_focal(gene_df, pw_df, tf_df, cfg, args):
    """8f-complementary focal module:
      (1) Pathway persistence (all Hallmark MH + immune M2/M5) flagged with 8f
          cross-tissue concordance — the headline 'both' set.
      (2) TF persistence restricted to broad regulator families, 8f-flagged.
      (3) Leading-edge drivers of the persistent immune pathways (8c LE, chunked).
    Returns (focal_pw, focal_tf, le_drivers).
    """
    print("\n[View 7] 8f-complementary focal module")
    focal_pw = pd.DataFrame()
    focal_tf = pd.DataFrame()
    le_drivers = pd.DataFrame()

    # ---- (1) Pathways: MH + immune M2/M5, with 8f concordance flag ----
    if not pw_df.empty:
        fp = pw_df.copy()
        is_mh = fp.get("collection", pd.Series("", index=fp.index)) == "MH"
        is_immune = fp["pathway"].str.upper().str.contains(
            "|".join(IMMUNE_PATHWAY_KEYWORDS), na=False)
        fp = fp[is_mh | is_immune].copy()
        if not fp.empty:
            conc = _load_8f_concordant(cfg, "pathway")
            if not conc.empty:
                # 8f keys: arm, brain_level, brain_celltype, pathway
                conc_keyed = conc.rename(columns={"brain_level": "level",
                                                  "brain_celltype": "celltype"})
                key = ["arm", "level", "pathway"]
                have = [k for k in key if k in conc_keyed.columns]
                n8f = (conc_keyed.groupby(have).size()
                       .reset_index(name="n_8f_concordant"))
                fp = fp.merge(n8f, on=have, how="left")
                fp["n_8f_concordant"] = fp["n_8f_concordant"].fillna(0).astype(int)
            else:
                fp["n_8f_concordant"] = 0
            fp["cross_tissue_8f"] = fp["n_8f_concordant"] > 0
            fp["both_persistent_and_8f"] = (
                (fp["persistence_class"] == "persistent") & fp["cross_tissue_8f"])
            n_both = int(fp["both_persistent_and_8f"].sum())
            print(f"  Focal pathways: {len(fp):,} rows; "
                  f"{n_both:,} BOTH within-brain-persistent AND 8f cross-tissue-concordant")
            focal_pw = fp

    # ---- (2) TFs: broad regulator families, 8f-flagged ----
    if not tf_df.empty:
        ft = tf_df[tf_df["TF"].str.upper().str.contains(
            "|".join(FOCAL_TF_FAMILIES), na=False)].copy()
        if not ft.empty:
            conc = _load_8f_concordant(cfg, "tf")
            if not conc.empty:
                conc_keyed = conc.rename(columns={"brain_level": "level",
                                                  "brain_celltype": "celltype"})
                key = [k for k in ["arm", "level", "TF"] if k in conc_keyed.columns]
                n8f = conc_keyed.groupby(key).size().reset_index(name="n_8f_concordant")
                ft = ft.merge(n8f, on=key, how="left")
                ft["n_8f_concordant"] = ft["n_8f_concordant"].fillna(0).astype(int)
            else:
                ft["n_8f_concordant"] = 0
            ft["cross_tissue_8f"] = ft["n_8f_concordant"] > 0
            print(f"  Focal TFs: {len(ft):,} rows (families: {', '.join(FOCAL_TF_FAMILIES)})")
            focal_tf = ft

    # ---- (3) Leading-edge drivers of persistent immune pathways ----
    if not focal_pw.empty:
        persistent_immune = focal_pw[
            (focal_pw["persistence_class"] == "persistent") &
            (focal_pw["pathway"].str.upper().str.contains(
                "|".join(IMMUNE_PATHWAY_KEYWORDS), na=False))]
        targets = sorted(persistent_immune["pathway"].unique())
        if targets:
            le_path = (Path(cfg["results_dir"]) / "tables" / "08c_pathways" /
                       "08c_pathway_leading_edge.csv")
            if le_path.is_file():
                print(f"  Leading-edge dive over {len(targets)} persistent immune "
                      f"pathways (chunked read of {le_path.name})")
                le_drivers = _leading_edge_drivers(le_path, set(targets), args)
            else:
                print(f"  [info] {le_path.name} not found; leading-edge dive skipped")

    return focal_pw, focal_tf, le_drivers


def _leading_edge_drivers(le_path, target_pathways, args, chunksize=2_000_000):
    """Chunked read of the (large) 08c leading-edge table. For the target
    pathways at sex=combined, count how often each gene appears in the leading
    edge across (celltype, level, age), per pathway. Returns a tidy ranking."""
    cols = ["sex", "level", "celltype", "group_level", "collection",
            "pathway", "gene", "log2FC", "direction"]
    keep = []
    for ch in pd.read_csv(le_path, usecols=lambda c: c in cols,
                          chunksize=chunksize, low_memory=False):
        if "sex" in ch.columns:
            ch = ch[ch["sex"] == "combined"]
        ch = ch[ch["pathway"].isin(target_pathways)]
        if not ch.empty:
            keep.append(ch)
    if not keep:
        return pd.DataFrame()
    le = pd.concat(keep, ignore_index=True)
    # Recurrence of each gene per pathway (across celltype/level/age)
    rec = (le.groupby(["pathway", "gene"])
           .agg(n_occurrences=("gene", "size"),
                n_levels=("level", "nunique") if "level" in le.columns else ("gene", "size"),
                n_celltypes=("celltype", "nunique") if "celltype" in le.columns else ("gene", "size"),
                mean_log2FC=("log2FC", "mean") if "log2FC" in le.columns else ("gene", "size"))
           .reset_index().sort_values(["pathway", "n_occurrences"],
                                      ascending=[True, False]))
    return rec


def plot_focal_pathway_persistence(focal_pw, plots_root):
    """Heatmap per arm: focal pathway x celltype, cell = persistence rank,
    annotated with * where also 8f cross-tissue concordant. Whole level only
    (the robust panel); regional rows stay in the CSV."""
    if focal_pw.empty:
        return
    df = focal_pw[focal_pw["level"] == "whole"]
    if df.empty:
        return
    for arm in sorted(df["arm"].unique()):
        a = df[df["arm"] == arm].copy()
        a["rank"] = a["persistence_class"].map(PERSISTENCE_RANK)
        # keep pathways that are persistent or near (rank<=2) in >=1 celltype
        keep_pw = a[a["rank"] <= 2]["pathway"].unique()
        a = a[a["pathway"].isin(keep_pw)]
        if a.empty:
            continue
        mat = a.pivot_table(index="pathway", columns="celltype",
                            values="rank", aggfunc="min")
        flag = a.pivot_table(index="pathway", columns="celltype",
                             values="cross_tissue_8f", aggfunc="max")
        mat = mat.reindex(index=sorted(mat.index))
        fig, ax = plt.subplots(figsize=(max(6, mat.shape[1] * 0.6 + 3),
                                        max(4, mat.shape[0] * 0.4 + 2)))
        im = ax.imshow(mat.values, aspect="auto", cmap="viridis_r",
                       vmin=0, vmax=7)
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("persistence rank (0=persistent ... 7=none)")
        ax.set_xticks(range(mat.shape[1]))
        ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(mat.shape[0]))
        ax.set_yticklabels(mat.index, fontsize=6)
        # star where 8f-concordant
        for i, pw in enumerate(mat.index):
            for j, ct in enumerate(mat.columns):
                try:
                    if bool(flag.loc[pw, ct]) and np.isfinite(mat.values[i, j]):
                        ax.text(j, i, "*", ha="center", va="center",
                                color="white", fontsize=9, fontweight="bold")
                except (KeyError, TypeError):
                    pass
        ax.set_title(f"Focal pathway persistence — {arm} arm [whole]\n"
                     f"(* = also cross-tissue concordant in 8f; darker = more persistent)",
                     fontsize=9)
        _save_fig(fig, plots_root / f"focal_pathway_persistence_{arm}.png")


def plot_trajectory_shape(persist_df, kind_label, plots_root):
    """Stacked bar of trajectory_shape composition within each persistence class."""
    if persist_df.empty or "trajectory_shape" not in persist_df.columns:
        return
    df = persist_df[persist_df["level"] == "whole"]
    if df.empty:
        return
    shape_order = ["amplifying", "stable", "attenuating",
                   "directionswap", "single_age", "undetermined"]
    shape_colors = {"amplifying": "#d73027", "stable": "#cccccc",
                    "attenuating": "#4575b4", "directionswap": "#984ea3",
                    "single_age": "#eeeeee", "undetermined": "#fdae61"}
    for arm in sorted(df["arm"].unique()):
        a = df[df["arm"] == arm]
        ct = (a.groupby(["persistence_class", "trajectory_shape"])
              .size().unstack(fill_value=0))
        classes = [c for c in PERSISTENCE_RANK if c in ct.index]
        ct = ct.reindex(index=classes)
        shapes = [s for s in shape_order if s in ct.columns]
        ct = ct[shapes]
        fig, ax = plt.subplots(figsize=(8, 5))
        bottom = np.zeros(len(ct))
        for s in shapes:
            ax.bar(range(len(ct)), ct[s].values, bottom=bottom,
                   label=s, color=shape_colors.get(s, "#999999"))
            bottom += ct[s].values
        ax.set_xticks(range(len(ct)))
        ax.set_xticklabels(ct.index, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel(f"# {kind_label}")
        ax.legend(fontsize=7, ncol=2)
        ax.set_title(f"Trajectory shape within persistence class — {arm} arm "
                     f"[whole]\n(amplifying = effect grows P1→3mo; attenuating = resolving)",
                     fontsize=9)
        _save_fig(fig, plots_root / f"trajectory_shape_{kind_label}_{arm}.png")


def plot_persistence_x_disruption(cross_df, plots_root):
    """Contingency heatmap: 8g persistence_class x 8b disruption_class, per arm."""
    if cross_df.empty:
        return
    df = cross_df[cross_df["level"] == "whole"]
    if df.empty:
        return
    for arm in sorted(df["arm"].unique()):
        a = df[df["arm"] == arm]
        ct = (a.groupby(["persistence_class", "disruption_class"])
              .size().unstack(fill_value=0))
        classes = [c for c in PERSISTENCE_RANK if c in ct.index]
        ct = ct.reindex(index=classes).fillna(0)
        if ct.empty:
            continue
        fig, ax = plt.subplots(figsize=(max(6, ct.shape[1] * 1.1 + 2),
                                        max(4, ct.shape[0] * 0.5 + 2)))
        im = ax.imshow(ct.values, aspect="auto", cmap="Reds")
        plt.colorbar(im, ax=ax, label="# genes")
        ax.set_xticks(range(ct.shape[1]))
        ax.set_xticklabels(ct.columns, rotation=30, ha="right", fontsize=8)
        ax.set_yticks(range(ct.shape[0]))
        ax.set_yticklabels(ct.index, fontsize=8)
        for i in range(ct.shape[0]):
            for j in range(ct.shape[1]):
                v = int(ct.values[i, j])
                if v:
                    ax.text(j, i, str(v), ha="center", va="center", fontsize=7,
                            color="white" if v > ct.values.max() * 0.6 else "black")
        ax.set_title(f"Stress persistence (8g) × developmental disruption (8b) — "
                     f"{arm} arm [whole]\n(persistent × relaxed_only = durable stress "
                     f"effect on a developmentally-lost gene)", fontsize=9)
        _save_fig(fig, plots_root / f"persistence_x_disruption_{arm}.png")


# ============================================================================
# Main
# ============================================================================

def main():
    args = parse_args()
    cfg = load_config(args.config)
    tissue = cfg.get("tissue", "unknown")
    print(f"\n{'='*60}")
    print(f"Phase 8g: Cross-age persistence  [{tissue}]")
    print(f"{'='*60}")

    if tissue == "placenta":
        print("  [warn] Placenta has incomplete cross-age factorial "
              "(E12.5=Early+Relaxed, E18.5=Late+Relaxed). 8g is brain-only by "
              "design. Exiting cleanly.")
        return

    tdir = phase_table_dir(cfg, PHASE)
    pdir_root = Path(cfg["results_dir"]) / "plots" / PHASE
    pdir_root.mkdir(parents=True, exist_ok=True)

    print("\n[Loading] 8b DE + 8c pathway/TF tables")
    tables = load_tables(cfg)

    # --- Run comprehensive persistence views (1-3,5,6) ---
    gene_df = run_view1_gene_persistence(tables["de"], args)
    pw_df = run_view2_pathway_persistence(tables["pw"], args)
    tf_df = run_view3_tf_persistence(tables["tf"], args)

    # --- B: trajectory shape on each comprehensive table ---
    gene_df = add_trajectory_shape(gene_df, "log2FC", args)
    pw_df = add_trajectory_shape(pw_df, "NES", args)
    tf_df = add_trajectory_shape(tf_df, "activity", args)

    ovr_df = run_view5_early_vs_late(tables["de"], args)
    core_genes, core_pw = run_view6_core_signature(gene_df, pw_df)

    # --- C: persistence x developmental disruption ---
    cross_df = run_persistence_x_disruption(gene_df, tables["disruption"], args)

    # --- A: View 7 focal (8f-complementary) ---
    focal_pw, focal_tf, le_drivers = run_view7_focal(gene_df, pw_df, tf_df, cfg, args)

    # --- level_flag: whole = robust, regions = exploratory (applied uniformly) ---
    def _flag(df):
        if not df.empty and "level" in df.columns:
            df = df.copy()
            df["level_flag"] = np.where(df["level"] == "whole",
                                        "robust", "regional_exploratory")
        return df
    gene_df, pw_df, tf_df, ovr_df = (_flag(gene_df), _flag(pw_df),
                                     _flag(tf_df), _flag(ovr_df))
    core_genes, core_pw, focal_pw, focal_tf = (_flag(core_genes), _flag(core_pw),
                                               _flag(focal_pw), _flag(focal_tf))

    # --- Persist tables ---
    print("\n[Tables]")
    def _save(df, name, note=""):
        if df is not None and not df.empty:
            df.to_csv(tdir / name, index=False)
            print(f"  Saved: {name}  ({len(df):,} rows){note}")
    _save(gene_df, "08g_gene_persistence.csv")
    _save(pw_df, "08g_pathway_persistence.csv")
    _save(tf_df, "08g_tf_persistence.csv")
    _save(ovr_df, "08g_early_vs_late_overlap.csv")
    if not core_genes.empty:
        core_genes.to_csv(tdir / "08g_core_signature_genes.csv", index=False)
        n_core = core_genes["core_signature"].sum() if "core_signature" in core_genes.columns else len(core_genes)
        print(f"  Saved: 08g_core_signature_genes.csv  ({n_core:,} core) ← PAPER TABLE")
    _save(core_pw, "08g_core_signature_pathways.csv")
    # New analyses
    _save(cross_df, "08g_persistence_x_disruption.csv", "  [C]")
    _save(focal_pw, "08g_focal_pathway_persistence.csv", "  [A: 8f bridge]")
    _save(focal_tf, "08g_focal_tf_persistence.csv", "  [A]")
    _save(le_drivers, "08g_focal_leadingedge_drivers.csv", "  [A: LE drivers]")

    # --- Plots ---
    print("\n[Plots]")
    # Comprehensive persistence (1-3): whole in the main dir, each region nested.
    def _levels(df):
        return ["whole"] + [l for l in sorted(df["level"].unique()) if l != "whole"] \
            if (not df.empty and "level" in df.columns) else ["whole"]

    for df, label, feat, eff, num in (
        (gene_df, "genes", "gene", "log2FC", "01_gene_persistence"),
        (pw_df, "pathways", "pathway", "NES", "02_pathway_persistence"),
        (tf_df, "TFs", "TF", "activity", "03_tf_persistence"),
    ):
        if df.empty:
            continue
        base = pdir_root / num
        for lvl in _levels(df):
            sub = df[df["level"] == lvl]
            if sub.empty:
                continue
            ldir = base if lvl == "whole" else base / _slug(lvl)
            plot_persistence_class_barplot(sub, label, ldir)
            plot_persistence_heatmap(sub, feat, label, ldir, top_n_label=50)
            plot_effect_size_trajectory(sub, feat, eff, label, ldir,
                                        top_n_label=args.top_n_label)
        # B: trajectory-shape composition (whole only)
        plot_trajectory_shape(df, label, base)

    if not ovr_df.empty:
        for lvl in _levels(ovr_df):
            sub = ovr_df[ovr_df["level"] == lvl]
            if sub.empty:
                continue
            ld = pdir_root / "04_early_vs_late"
            ld = ld if lvl == "whole" else ld / _slug(lvl)
            plot_early_vs_late_overlap(sub, ld, top_n_label=args.top_n_label)
    if not core_genes.empty or not core_pw.empty:
        plot_core_signature(core_genes, core_pw,
                            pdir_root / "05_core_signature",
                            top_n_label=args.top_n_label)
    # C + A plots
    plot_persistence_x_disruption(cross_df, pdir_root / "06_persistence_x_disruption")
    plot_focal_pathway_persistence(focal_pw, pdir_root / "07_focal_8f_bridge")

    # --- Summary ---
    print(f"\n{'='*60}")
    print("Phase 8g complete.")
    print(f"  Tables: {tdir}")
    print(f"  Plots:  {pdir_root}")
    if not gene_df.empty:
        whole = gene_df[gene_df["level"] == "whole"]
        for arm in sorted(whole["arm"].unique()):
            n = ((whole["arm"] == arm) &
                 (whole["persistence_class"] == "persistent")).sum()
            print(f"  {arm} arm [whole]: {n:,} persistent (gene, celltype) calls")
    if not core_genes.empty and "core_signature" in core_genes.columns:
        n_core = int(core_genes["core_signature"].sum())
        print(f"\n  PAPER TABLE: {n_core:,} core-signature genes (persistent in BOTH arms)")
    if not focal_pw.empty and "both_persistent_and_8f" in focal_pw.columns:
        n_both = int(focal_pw["both_persistent_and_8f"].sum())
        print(f"  8f BRIDGE: {n_both:,} pathway calls BOTH within-brain-persistent "
              f"AND cross-tissue-concordant (8f)")
    if not cross_df.empty and "persistent_and_LOST" in cross_df.columns:
        print(f"  8b CROSS-REF: {int(cross_df['persistent_and_LOST'].sum()):,} "
              f"persistent-AND-developmentally-LOST gene calls")
    print()


if __name__ == "__main__":
    main()
