#!/usr/bin/env python
"""
08g_cross_age.py — Phase 8g: Cross-age persistence analysis.

Operates entirely on existing 8b/8c tables — no re-running of DE/GSEA.
Six analytical views answering: "of the signals deregulated by prenatal stress,
which persist across development, which resolve, which emerge later?"

Two stress arms:
  - Early: early_vs_relaxed_per_age at {P1, 4W, 3mo}
  - Late:  late_vs_relaxed_per_age at {P1, 4W, 3mo}  (P1 carries pool confound)

Persistence classes per gene (or pathway / TF) × celltype × arm:
  persistent      — DE at P1 AND 4W AND 3mo, SAME DIRECTION
  resolving_early — DE at P1 AND 4W, not 3mo
  established_late— DE at 4W AND 3mo, not P1
  P1_only         — DE at P1 only
  transient_4W    — DE at 4W only
  emergent_3mo    — DE at 3mo only
  P1_3mo_only     — DE at P1 AND 3mo, not 4W (unusual)
  none            — not DE at any age

Six views:
  1. Gene-level persistence       — 08b DE → classification + trajectory plots
  2. Pathway-level persistence    — 08c GSEA → same classification
  3. TF-level persistence         — 08c TF activity → same
  4. Effect-size trajectories     — top persistent features, log2FC/NES vs age
  5. Early vs Late at each age    — hypergeometric overlap + rank-rank scatter
  6. Cross-arm core signature     — features persistent in BOTH arms (paper table)

Outputs:
  plots/08g_cross_age/{01_gene,02_pathway,03_tf,04_early_vs_late,05_core_signature}/
  tables/08g_cross_age/
    08g_gene_persistence.csv            KEY: gene × celltype × arm × class
    08g_pathway_persistence.csv
    08g_tf_persistence.csv
    08g_early_vs_late_overlap.csv
    08g_core_signature_genes.csv        cross-arm persistent (paper table)
    08g_core_signature_pathways.csv

Usage:
  uv run python scripts/08g_cross_age.py --config config/dev_split.yaml
  uv run python scripts/08g_cross_age.py --config config/brain.yaml
  # placenta has incomplete factorial (no cross-age comparison possible) so
  # this script is brain-only by design.
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
    return p.parse_args()


# ============================================================================
# Loaders
# ============================================================================

def _slug(s: str) -> str:
    return str(s).replace(" ", "_").replace("/", "-").replace(".", "")


def _extract_age(group_level):
    if not isinstance(group_level, str):
        return None
    for part in group_level.split("_"):
        if part.startswith("age-"):
            return part[4:]
    return None


def load_tables(cfg):
    base = Path(cfg["results_dir"]) / "tables"
    paths = {
        "de": base / "08b_de" / "08b_de_results.csv",
        "pw": base / "08c_pathways" / "08c_pathway_results.csv",
        "tf": base / "08c_pathways" / "08c_tf_activity.csv",
    }
    out = {}
    for kind, p in paths.items():
        if p.is_file():
            df = pd.read_csv(p)
            print(f"  {kind}: {len(df):,} rows from {p.name}")
            out[kind] = df
        else:
            print(f"  [info] {kind} table missing: {p.name}")
            out[kind] = pd.DataFrame()
    return out


def prep_de(df, contrast):
    """Filter 08b DE to one contrast, add 'age', drop missing genes."""
    if df.empty or "contrast" not in df.columns:
        return pd.DataFrame()
    sub = df[df["contrast"] == contrast].copy()
    if sub.empty:
        return pd.DataFrame()
    sub = sub.dropna(subset=["gene", "stat"])
    sub["age"] = sub["group_level"].map(_extract_age)
    return sub.dropna(subset=["age"])[
        ["celltype", "gene", "stat", "padj", "log2FC", "age"]
    ]


def prep_pw(df, contrast):
    """Filter 08c pathway. 08c writes pathway name as 'source', FDR as 'FDR'."""
    if df.empty or "contrast" not in df.columns:
        return pd.DataFrame()
    sub = df[df["contrast"] == contrast].copy()
    if sub.empty:
        return pd.DataFrame()
    sub["age"] = sub["group_level"].map(_extract_age)
    rename = {}
    if "source" in sub.columns and "pathway" not in sub.columns:
        rename["source"] = "pathway"
    if "FDR" in sub.columns and "padj" not in sub.columns:
        rename["FDR"] = "padj"
    if rename:
        sub = sub.rename(columns=rename)
    return sub.dropna(subset=["age"])


def prep_tf(df, contrast):
    """Filter 08c TF activity. Columns: contrast, group_level, celltype, TF,
    activity_score, pvalue, FDR, direction."""
    if df.empty or "contrast" not in df.columns:
        return pd.DataFrame()
    sub = df[df["contrast"] == contrast].copy()
    if sub.empty:
        return pd.DataFrame()
    sub["age"] = sub["group_level"].map(_extract_age)
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
        sub = prep_de(de_df, arm["contrast"])
        if sub.empty:
            print(f"  [skip] {arm['arm']} arm: no DE rows")
            continue
        sub = sub[sub["age"].isin(arm["ages"])].copy()
        # Per (celltype, gene, age) — direction column
        sub["direction"] = "none"
        sig = (sub["padj"] < args.padj_cutoff) & \
              (sub["log2FC"].abs() > args.logfc_cutoff)
        sub.loc[sig & (sub["log2FC"] > 0), "direction"] = "up"
        sub.loc[sig & (sub["log2FC"] < 0), "direction"] = "down"

        # Keep only genes that are DE at ≥1 age (otherwise the classification
        # table is just 'none' for every untested gene — useless and huge).
        any_sig = sub.groupby(["celltype", "gene"])["direction"].apply(
            lambda x: (x != "none").any())
        keep_keys = any_sig[any_sig].index
        sub = sub.set_index(["celltype", "gene"]).loc[
            sub.set_index(["celltype", "gene"]).index.isin(keep_keys)].reset_index()

        if sub.empty:
            continue

        classified = classify_dataframe(
            sub, age_col="age", direction_col="direction",
            feature_cols=("celltype", "gene"))
        classified["arm"] = arm["arm"]
        # Carry P1 confound flag for the Late arm (since P1 row entered)
        classified["confound_note"] = arm["confound_flags"].get("P1", "") \
            if arm["arm"] == "Late" else ""
        # Attach the actual log2FC per age (useful in output)
        wide = (sub.pivot_table(index=["celltype", "gene"],
                                 columns="age", values="log2FC",
                                 aggfunc="first")
                 .reset_index())
        wide.columns = ["celltype", "gene"] + [f"{c}_log2FC" for c in wide.columns[2:]]
        classified = classified.merge(wide, on=["celltype", "gene"], how="left")
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
        sub = prep_pw(pw_df, arm["contrast"])
        if sub.empty:
            continue
        sub = sub[sub["age"].isin(arm["ages"])].copy()
        sub["direction"] = "none"
        sig = sub["padj"] < args.pathway_fdr_cutoff
        sub.loc[sig & (sub["NES"] > 0), "direction"] = "up"
        sub.loc[sig & (sub["NES"] < 0), "direction"] = "down"

        any_sig = sub.groupby(["celltype", "pathway"])["direction"].apply(
            lambda x: (x != "none").any())
        keep = any_sig[any_sig].index
        sub = sub.set_index(["celltype", "pathway"]).loc[
            sub.set_index(["celltype", "pathway"]).index.isin(keep)].reset_index()
        if sub.empty:
            continue

        classified = classify_dataframe(
            sub, age_col="age", direction_col="direction",
            feature_cols=("celltype", "pathway"))
        classified["arm"] = arm["arm"]
        classified["confound_note"] = arm["confound_flags"].get("P1", "") \
            if arm["arm"] == "Late" else ""
        wide = (sub.pivot_table(index=["celltype", "pathway"],
                                 columns="age", values="NES", aggfunc="first")
                 .reset_index())
        wide.columns = ["celltype", "pathway"] + [f"{c}_NES" for c in wide.columns[2:]]
        # Also carry collection column if present
        if "collection" in sub.columns:
            coll = sub.groupby(["celltype", "pathway"])["collection"].first().reset_index()
            classified = classified.merge(coll, on=["celltype", "pathway"], how="left")
        classified = classified.merge(wide, on=["celltype", "pathway"], how="left")
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
        sub = prep_tf(tf_df, arm["contrast"])
        if sub.empty:
            continue
        sub = sub[sub["age"].isin(arm["ages"])].copy()
        sub["direction"] = "none"
        sig = sub["padj"] < args.pathway_fdr_cutoff
        sub.loc[sig & (sub["activity_score"] > 0), "direction"] = "up"
        sub.loc[sig & (sub["activity_score"] < 0), "direction"] = "down"

        any_sig = sub.groupby(["celltype", "TF"])["direction"].apply(
            lambda x: (x != "none").any())
        keep = any_sig[any_sig].index
        sub = sub.set_index(["celltype", "TF"]).loc[
            sub.set_index(["celltype", "TF"]).index.isin(keep)].reset_index()
        if sub.empty:
            continue

        classified = classify_dataframe(
            sub, age_col="age", direction_col="direction",
            feature_cols=("celltype", "TF"))
        classified["arm"] = arm["arm"]
        classified["confound_note"] = arm["confound_flags"].get("P1", "") \
            if arm["arm"] == "Late" else ""
        wide = (sub.pivot_table(index=["celltype", "TF"],
                                 columns="age", values="activity_score",
                                 aggfunc="first").reset_index())
        wide.columns = ["celltype", "TF"] + [f"{c}_activity" for c in wide.columns[2:]]
        classified = classified.merge(wide, on=["celltype", "TF"], how="left")
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

    es = prep_de(de_df, "early_vs_relaxed_per_age")
    ls = prep_de(de_df, "late_vs_relaxed_per_age")
    if es.empty or ls.empty:
        print("  [skip] missing ES or LS DE rows")
        return pd.DataFrame()

    rows = []
    for age in ("P1", "4W", "3mo"):
        es_age = es[es["age"] == age]
        ls_age = ls[ls["age"] == age]
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
                    "age": age, "celltype": ct, "direction": direction,
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
    # BH-FDR within (age, direction)
    fdr_parts = []
    for _, g in df.groupby(["age", "direction"]):
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

        # Join on (celltype, feature) and require same direction in both arms
        key_cols = ["celltype", feature_col]
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

    # --- Run six views ---
    gene_df = run_view1_gene_persistence(tables["de"], args)
    pw_df = run_view2_pathway_persistence(tables["pw"], args)
    tf_df = run_view3_tf_persistence(tables["tf"], args)
    ovr_df = run_view5_early_vs_late(tables["de"], args)
    core_genes, core_pw = run_view6_core_signature(gene_df, pw_df)

    # --- Persist tables ---
    print("\n[Tables]")
    if not gene_df.empty:
        gene_df.to_csv(tdir / "08g_gene_persistence.csv", index=False)
        print(f"  Saved: 08g_gene_persistence.csv  ({len(gene_df):,} rows)")
    if not pw_df.empty:
        pw_df.to_csv(tdir / "08g_pathway_persistence.csv", index=False)
        print(f"  Saved: 08g_pathway_persistence.csv  ({len(pw_df):,} rows)")
    if not tf_df.empty:
        tf_df.to_csv(tdir / "08g_tf_persistence.csv", index=False)
        print(f"  Saved: 08g_tf_persistence.csv  ({len(tf_df):,} rows)")
    if not ovr_df.empty:
        ovr_df.to_csv(tdir / "08g_early_vs_late_overlap.csv", index=False)
        print(f"  Saved: 08g_early_vs_late_overlap.csv  ({len(ovr_df):,} rows)")
    if not core_genes.empty:
        core_genes.to_csv(tdir / "08g_core_signature_genes.csv", index=False)
        n_core = core_genes["core_signature"].sum() if "core_signature" in core_genes.columns else len(core_genes)
        print(f"  Saved: 08g_core_signature_genes.csv  ({n_core:,} core) ← PAPER TABLE")
    if not core_pw.empty:
        core_pw.to_csv(tdir / "08g_core_signature_pathways.csv", index=False)
        n_core = core_pw["core_signature"].sum() if "core_signature" in core_pw.columns else len(core_pw)
        print(f"  Saved: 08g_core_signature_pathways.csv  ({n_core:,} core)")

    # --- Plots ---
    print("\n[Plots]")
    if not gene_df.empty:
        gene_dir = pdir_root / "01_gene_persistence"
        plot_persistence_class_barplot(gene_df, "genes", gene_dir)
        plot_persistence_heatmap(gene_df, "gene", "genes", gene_dir,
                                  top_n_label=50)
        plot_effect_size_trajectory(gene_df, "gene", "log2FC", "genes",
                                     gene_dir, top_n_label=args.top_n_label)
    if not pw_df.empty:
        pw_dir = pdir_root / "02_pathway_persistence"
        plot_persistence_class_barplot(pw_df, "pathways", pw_dir)
        plot_persistence_heatmap(pw_df, "pathway", "pathways", pw_dir,
                                  top_n_label=40)
        plot_effect_size_trajectory(pw_df, "pathway", "NES", "pathways",
                                     pw_dir, top_n_label=args.top_n_label)
    if not tf_df.empty:
        tf_dir = pdir_root / "03_tf_persistence"
        plot_persistence_class_barplot(tf_df, "TFs", tf_dir)
        plot_persistence_heatmap(tf_df, "TF", "TFs", tf_dir, top_n_label=40)
        plot_effect_size_trajectory(tf_df, "TF", "activity", "TFs",
                                     tf_dir, top_n_label=args.top_n_label)
    if not ovr_df.empty:
        plot_early_vs_late_overlap(ovr_df, pdir_root / "04_early_vs_late",
                                    top_n_label=args.top_n_label)
    if not core_genes.empty or not core_pw.empty:
        plot_core_signature(core_genes, core_pw,
                             pdir_root / "05_core_signature",
                             top_n_label=args.top_n_label)

    # --- Summary ---
    print(f"\n{'='*60}")
    print("Phase 8g complete.")
    print(f"  Tables: {tdir}")
    print(f"  Plots:  {pdir_root}")
    if not gene_df.empty:
        for arm in sorted(gene_df["arm"].unique()):
            persist = (gene_df["arm"] == arm) & \
                      (gene_df["persistence_class"] == "persistent")
            print(f"  {arm} arm: {persist.sum():,} persistent (gene, celltype) calls")
    if not core_genes.empty and "core_signature" in core_genes.columns:
        n_core = int(core_genes["core_signature"].sum())
        print(f"\n  PAPER TABLE: {n_core:,} core-signature genes (persistent in BOTH arms)")
        print(f"  See: tables/08g_cross_age/08g_core_signature_genes.csv")
    print()


if __name__ == "__main__":
    main()
