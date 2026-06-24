#!/usr/bin/env python
"""
08e_comms_summary.py — Phase 8e PLOTTING. Reads ONLY the 08e CSVs (mirrors the
8b/8c compute/summary split). No LIANA compute, no adata access (except the
static LR→pathway resource via liana.rs, fetched once for heatmap annotations).

Reads (from tables/08e_communication{,_subcluster_<slug>}/):
  08e_lr_baseline.csv       — arm 1, DESCRIPTIVE landscape (pooled cells)
  08e_lr_differential.csv   — arm 2, PRIMARY inferential (df_to_lr on 8b stats)
  08e_lr_per_donor.csv      — arm 3 raw (per-donor activity)
  08e_lr_quantified.csv     — arm 3 stats (MW-U across donors, FDR)
  08e_sender_receiver.csv   — derived sender/receiver summary

ARM-HIERARCHY FRAMING (locked):
  - baseline (02_baseline/, 04_sender_receiver/ landscape figs): DESCRIPTIVE only,
    NOT a stress test (pooled cells). Tagged in titles/folder README.
  - differential (03_differential/): PRIMARY inferential arm. Stress claims here.
    Shares signal with 8b by construction (transcriptional proxy for signalling).
  - per-donor (05_per_donor/): independent-design corroboration, expected null at
    n≈4. Read effect sizes, not p-values.

Output: plots/08e_communication{,_subcluster_<slug>}/{01_overview,02_baseline,
        03_differential,04_sender_receiver,05_per_donor}/

Usage:
  uv run python scripts/08e_comms_summary.py --config config/placenta.yaml
  uv run python scripts/08e_comms_summary.py --config config/brain.yaml \\
      --top-lr 25 --magnitude-cutoff 0.05
"""

import argparse
import sys
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _utils import load_config

import _08e_plots_baseline as pb
import _08e_plots_differential as pd_diff
import _08e_plots_perdonor as pp
import _08e_plots_stats as ps          # statistically-grounded volcano/dotplot/network/SR

PHASE = "08e_communication"
GROUP_ORDER = ["Relaxed", "Early_Stress", "Late_Stress"]

# Folders that hold DESCRIPTIVE (not-a-stress-test) figures — a README is dropped
# into each so the distinction survives outside this script.
DESCRIPTIVE_NOTE = (
    "Figures here are a DESCRIPTIVE communication landscape from the baseline arm "
    "(rank_aggregate on POOLED cells per group×age). They are NOT a statistical "
    "stress test — pooling cells discards donor structure. Stress claims come from "
    "03_differential/ (df_to_lr on 8b pseudobulk Wald stats) and 05_per_donor/ "
    "(Mann-Whitney U across donors).\n")


def parse_args():
    p = argparse.ArgumentParser(description="Phase 8e plotting (reads CSVs only)")
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--subcluster", default=None,
                   help="Plot a subcluster run's CSVs (slug; e.g. opc_oligodendrocytes).")
    p.add_argument("--top-lr", type=int, default=25,
                   help="Top-N LR pairs in dotplots/volcanos (default 25).")
    p.add_argument("--top-lr-large", type=int, default=50,
                   help="Top-N for the large LR dotplot (default 50).")
    p.add_argument("--magnitude-cutoff", type=float, default=0.05,
                   help="magnitude_rank cutoff defining an 'active' interaction (default 0.05).")
    p.add_argument("--focus-celltypes", default=None,
                   help="Comma-separated cell types to focus Δ heatmaps on "
                        "(default: stress_focused_cell_types from config, else all).")
    p.add_argument("--zscore-rows", action="store_true",
                   help="Z-score rows of the Δ LR heatmap.")
    return p.parse_args()


# ---------------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------------

def read_csv(tdir, name):
    p = tdir / name
    if not p.is_file():
        print(f"  [info] {name} not present — related plots skipped.")
        return pd.DataFrame()
    df = pd.read_csv(p, low_memory=False)
    print(f"  Loaded {name}: {len(df):,} rows")
    return df


def resolve_focus_celltypes(baseline_df, cfg, explicit):
    """Cell-type universe comes from the baseline source/target union (NO adata).
    If --focus-celltypes given, use that; else fall back to the config stress list
    intersected with observed types; else None (no focus filtering)."""
    if baseline_df.empty:
        observed = set()
    else:
        observed = set(baseline_df["source"].astype(str)) | set(baseline_df["target"].astype(str))
    if explicit:
        req = [c.strip() for c in explicit.split(",") if c.strip()]
        miss = [c for c in req if c not in observed]
        if miss:
            print(f"  [warn] --focus-celltypes not observed: {miss}")
        keep = [c for c in req if c in observed]
        return keep or None
    stress = cfg.get("stress_focused_cell_types", [])
    keep = [c for c in stress if c in observed]
    if keep:
        print(f"  Focus cell types (from config ∩ observed): {keep}")
        return keep
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    cfg = load_config(args.config)
    tissue = cfg.get("tissue", "unknown")
    ref_group = cfg.get("group_reference", "Relaxed")
    label = PHASE if not args.subcluster else f"{PHASE}_subcluster_{args.subcluster}"

    tdir = Path(cfg["results_dir"]) / "tables" / label
    if not tdir.is_dir():
        sys.exit(f"ERROR: {tdir} not found — run 08e_comms.py first.")
    pdir_root = Path(cfg["results_dir"]) / "plots" / label

    print(f"\n{'='*60}\nPhase 8e summary plots  [{tissue}]  ({label})\n{'='*60}")

    baseline = read_csv(tdir, "08e_lr_baseline.csv")
    differential = read_csv(tdir, "08e_lr_differential.csv")
    per_donor = read_csv(tdir, "08e_lr_per_donor.csv")
    quant = read_csv(tdir, "08e_lr_quantified.csv")
    sr = read_csv(tdir, "08e_sender_receiver.csv")

    # static LR→pathway map for heatmap annotations (liana resource, not adata)
    try:
        pathway_map = pd_diff.get_pathway_map()
    except Exception as e:
        print(f"  [warn] pathway map unavailable ({e}) — heatmap row colors skipped.")
        pathway_map = None

    focus = resolve_focus_celltypes(baseline, cfg, args.focus_celltypes)

    d_overview = pdir_root / "01_overview"
    d_base = pdir_root / "02_baseline"
    d_diff = pdir_root / "03_differential"
    d_sr = pdir_root / "04_sender_receiver"
    d_pd = pdir_root / "05_per_donor"
    for d in (d_overview, d_base, d_diff, d_sr, d_pd):
        d.mkdir(parents=True, exist_ok=True)
    (d_base / "README.txt").write_text(DESCRIPTIVE_NOTE)
    (d_sr / "README.txt").write_text(DESCRIPTIVE_NOTE)

    cutoff = args.magnitude_cutoff

    # ---- group×age and age iterators from what's actually present -----------
    if not baseline.empty:
        ages = sorted(baseline["age"].astype(str).unique())
        ga_present = (baseline[["group", "age"]].drop_duplicates()
                      .itertuples(index=False))
        ga_present = [(g, a) for g, a in ga_present]
    else:
        ages, ga_present = [], []

    def groups_at(age):
        gs = baseline.loc[baseline["age"].astype(str) == age, "group"].astype(str).unique()
        return [g for g in GROUP_ORDER if g in gs]

    def group_pairs_at(age):
        """Stress-vs-ref pairs present at an age (ref on the right)."""
        gs = groups_at(age)
        pairs = []
        for a, b in combinations(gs, 2):
            if b == ref_group:
                pairs.append((a, b))
            elif a == ref_group:
                pairs.append((b, a))
            else:
                pairs.append((a, b))
        return pairs

    # =======================================================================
    # 02 BASELINE — descriptive landscape (per group×age, per age, global)
    # =======================================================================
    if not baseline.empty:
        print("\n[02_baseline] descriptive landscape")
        for group, age in ga_present:
            try:
                pb.plot_chord_diagram(baseline, group, age, cutoff, d_base)
                pb.plot_network_graph(baseline, group, age, cutoff, d_base)
                pb.plot_interaction_count_heatmap(baseline, group, age, cutoff, d_base)
                pb.plot_baseline_dotplot(baseline, group, age, args.top_lr, d_base)
                pb.plot_large_lr_dotplot(baseline, group, age, args.top_lr_large, d_base)
                pb.plot_top_lr_per_celltype_pair(baseline, group, age, args.top_lr, d_base)
            except Exception as e:
                print(f"  [warn] baseline {group}/{age}: {e}")
        for age in ages:
            try:
                pb.plot_delta_chord_diagram(baseline, age, cutoff, d_base)
            except Exception as e:
                print(f"  [warn] delta_chord {age}: {e}")
        try:
            pb.plot_interaction_counts_barplot(baseline, ref_group, d_overview)
            pb.plot_pathway_activity_heatmap(baseline, cutoff, d_overview)
            pb.plot_lr_persistence_across_ages(baseline, cutoff, args.top_lr, d_overview)
        except Exception as e:
            print(f"  [warn] overview baseline plots: {e}")

    # =======================================================================
    # 03 DIFFERENTIAL — PRIMARY inferential arm
    # =======================================================================
    if not differential.empty:
        print("\n[03_differential] PRIMARY inferential arm (FDR from interaction_padj)")
        for cname in sorted(differential["contrast_name"].astype(str).unique()):
            for age in sorted(differential.loc[differential["contrast_name"] == cname,
                                               "age"].astype(str).unique()):
                try:
                    ps.plot_differential_volcano(differential, cname, age, args.top_lr, d_diff)
                    ps.plot_differential_dotplot(differential, cname, age, args.top_lr, d_diff)
                    ps.plot_delta_network(differential, cname, age, d_diff, focus=focus)
                    ps.plot_delta_chord(differential, cname, age, d_diff, focus=focus)
                    ps.plot_delta_celltype_heatmap(differential, cname, age, d_diff)
                except Exception as e:
                    print(f"  [warn] differential {cname}/{age}: {e}")
        try:
            # writes 08e_lr_stress_signature_pivot.csv into tdir + the heatmap
            pd_diff.plot_stress_signature_heatmap(differential, tdir, d_diff, top_n=100)
        except Exception as e:
            print(f"  [warn] stress_signature_heatmap: {e}")

        # Δ LR heatmaps + rank-rank scatters are built from the BASELINE landscape
        # (group-pair differences), kept under 03 since they motivate stress contrasts.
        if not baseline.empty:
            for age in ages:
                for grp_a, grp_b in group_pairs_at(age):
                    try:
                        pd_diff.plot_delta_lr_heatmap(
                            baseline, grp_a, grp_b, age, cutoff, d_diff,
                            focus_celltypes=focus, zscore_rows=args.zscore_rows,
                            pathway_map=pathway_map)
                    except Exception as e:
                        print(f"  [warn] delta_lr_heatmap {grp_a}-{grp_b}/{age}: {e}")
                pairs = group_pairs_at(age)
                for ca, cb in combinations(pairs, 2):
                    try:
                        pd_diff.plot_rank_rank_scatter(baseline, ca, cb, age, cutoff,
                                                       d_diff, focus_celltypes=focus)
                    except Exception as e:
                        print(f"  [warn] rank_rank {ca}|{cb}/{age}: {e}")

    # =======================================================================
    # 04 SENDER / RECEIVER — descriptive
    # =======================================================================
    if not baseline.empty:
        print("\n[04_sender_receiver] descriptive")
        for age in ages:
            try:
                ps.plot_sender_receiver_bubble(baseline, age, cutoff, d_sr)
            except Exception as e:
                print(f"  [warn] sr_bubble {age}: {e}")
            if not sr.empty:
                try:
                    pd_diff.plot_sender_receiver_heatmap(sr, age, ref_group, d_sr)
                except Exception as e:
                    print(f"  [warn] sr_heatmap {age}: {e}")
                for grp_a, grp_b in group_pairs_at(age):
                    try:
                        pd_diff.plot_delta_sender_receiver_heatmap(sr, grp_a, grp_b, age, d_sr)
                    except Exception as e:
                        print(f"  [warn] delta_sr {grp_a}-{grp_b}/{age}: {e}")

    # =======================================================================
    # 05 PER-DONOR — corroboration (per age; covers all contrasts internally)
    # =======================================================================
    if not per_donor.empty and not quant.empty:
        print("\n[05_per_donor] corroboration (animal as unit)")
        for age in sorted(quant["age"].astype(str).unique()):
            pdf_a = per_donor[per_donor["age"].astype(str) == age]
            q_a = quant[quant["age"].astype(str) == age]
            if pdf_a.empty or q_a.empty:
                continue
            try:
                pp.plot_per_donor_quantification(pdf_a, q_a, ref_group, args.top_lr, d_pd)
            except Exception as e:
                print(f"  [warn] per_donor {age}: {e}")

    print(f"\n{'='*60}\n✓ Phase 8e summary done. Plots in {pdir_root}/")
    for d in (d_overview, d_base, d_diff, d_sr, d_pd):
        n = len(list(d.glob("*.png")))
        print(f"    {d.name}/: {n} figures")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
