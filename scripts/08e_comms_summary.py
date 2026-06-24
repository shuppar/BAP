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
import _08e_plots_pathway as pw_plot   # per-pathway CCC graphs

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
    p.add_argument("--node-scheme", default="broad", choices=["broad", "subtype"],
                   help="Which compute output to plot: 'broad' → 08e_communication/, "
                        "'subtype' → 08e_communication_subtype/.")
    p.add_argument("--by-pathway", action="store_true",
                   help="Draw per-pathway CCC graphs (one chord+network per stress "
                        "pathway) into 06_by_pathway/. Uses config/stress_pathways_8e.yaml.")
    p.add_argument("--stress-spec", default="config/stress_pathways_8e.yaml", type=Path,
                   help="Pathway whitelist + focal map for --by-pathway.")
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
    if args.node_scheme == "subtype" and not args.subcluster:
        label = f"{PHASE}_subtype"

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

    def groups_at_in(df_slice, age):
        """Groups present at an age within an arbitrary baseline slice (ordered)."""
        gs = df_slice.loc[df_slice["age"].astype(str) == str(age), "group"].astype(str).unique()
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
                    ps.plot_sender_receiver_updown_bars(differential, cname, age, d_diff)
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

    # =======================================================================
    # 06 BY-PATHWAY — one chord+network per stress pathway (noise-reduced)
    # =======================================================================
    if args.by_pathway:
        print("\n[06_by_pathway] per-pathway CCC graphs")
        import yaml
        if not args.stress_spec.is_file():
            print(f"  [warn] {args.stress_spec} not found — skipping by-pathway.")
        else:
            spec = yaml.safe_load(args.stress_spec.read_text())
            le_path = tdir.parent.parent / "tables" / "08c_pathways" / "08c_pathway_leading_edge.csv"
            # le lives under the MAIN 8c dir regardless of node scheme/subcluster
            le_path = Path(cfg["results_dir"]) / "tables" / "08c_pathways" / "08c_pathway_leading_edge.csv"
            if not le_path.is_file():
                print(f"  [warn] 8c leading-edge CSV not found: {le_path} — skipping.")
            else:
                genesets = pw_plot.load_pathway_genesets(cfg, tissue, spec, le_path)
                print(f"  loaded {len(genesets)} pathway gene sets: {sorted(_slugmap(genesets))}")
                # Build contrast specs from what's present in baseline (group/age) +
                # the canonical primary contrasts. test vs ref(=Relaxed).
                contrasts = _build_pathway_contrasts(baseline, differential, ref_group)
                d_pw = pdir_root / "06_by_pathway"
                d_pw.mkdir(parents=True, exist_ok=True)
                pw_plot.plot_by_pathway(
                    genesets, baseline,
                    differential if not differential.empty else None,
                    contrasts, d_pw)

                # ---- cross-scheme companion (broad vs subtype, whole level) ----
                # Load the OTHER scheme's baseline CSV. This run's `label` is either
                # 08e_communication (broad) or 08e_communication_subtype.
                base_tables = Path(cfg["results_dir"]) / "tables"
                broad_csv = base_tables / "08e_communication" / "08e_lr_baseline.csv"
                sub_csv = base_tables / "08e_communication_subtype" / "08e_lr_baseline.csv"
                if broad_csv.is_file() and sub_csv.is_file():
                    print("  cross-scheme companion (broad vs subtype)...")
                    bb = pd.read_csv(broad_csv, low_memory=False)
                    bs = pd.read_csv(sub_csv, low_memory=False)
                    pw_plot.plot_cross_scheme_companion(
                        genesets, bb, bs, contrasts, d_pw)
                else:
                    print("  [info] cross-scheme companion needs BOTH broad + subtype "
                          "baseline CSVs — run both schemes first; skipping.")

    # =======================================================================
    # 07 FOCAL-FAN GRIDS — readable per-cell-type fan layout (per-group + Δ)
    #     across whole + per-pathway + regional. Subtype nodes handled
    #     automatically when run with --node-scheme subtype.
    # =======================================================================
    if not baseline.empty:
        print("\n[07_focal_grids] focal-fan network grids (per-group + Δ)")
        levels_present = (sorted(baseline["level"].astype(str).unique())
                          if "level" in baseline.columns else ["whole"])

        def _grids_for_slice(df_slice, age, out_dir, tag, level=None,
                             diff_slice=None):
            out_dir.mkdir(parents=True, exist_ok=True)
            # per-group descriptive grids (reuse proven plot_network_graph)
            for g in groups_at_in(df_slice, age):
                try:
                    pb.plot_network_graph(df_slice, g, age, cutoff, out_dir, level=level)
                except Exception as e:
                    print(f"  [warn] grid {tag} {g}/{age}: {e}")
            # Δ grids — baseline arm (both tissues): count + magnitude
            for grp_a, grp_b in group_pairs_at(age):
                for metric in ("count", "magnitude"):
                    try:
                        pb.plot_delta_network_grid(
                            df_slice, grp_a, grp_b, age, cutoff, out_dir,
                            level=level, metric=metric, arm="baseline")
                    except Exception as e:
                        print(f"  [warn] delta_grid baseline/{metric} {tag} "
                              f"{grp_a}-{grp_b}/{age}: {e}")
            # Δ grids — differential arm (placenta): count + magnitude
            if diff_slice is not None and not diff_slice.empty:
                for cname in sorted(diff_slice.loc[diff_slice["age"].astype(str) == str(age),
                                                   "contrast_name"].astype(str).unique()):
                    # test/ref are encoded in the contrast; pass through for title only
                    tg, rg = _contrast_groups(cname, ref_group)
                    for metric in ("count", "magnitude"):
                        try:
                            pb.plot_delta_network_grid(
                                df_slice, tg, rg, age, cutoff, out_dir,
                                level=level, metric=metric, arm="differential",
                                differential_df=diff_slice, contrast_name=cname)
                        except Exception as e:
                            print(f"  [warn] delta_grid diff/{metric} {tag} "
                                  f"{cname}/{age}: {e}")

        # whole + regional levels
        for lvl in levels_present:
            df_lvl = baseline[baseline.get("level", "whole").astype(str) == lvl] \
                if "level" in baseline.columns else baseline
            # differential is whole-only (no level column) → only feed it at whole
            diff_for_lvl = differential if (lvl == "whole" and not differential.empty) else None
            for age in sorted(df_lvl["age"].astype(str).unique()):
                _grids_for_slice(df_lvl, age, pdir_root / "07_focal_grids" / lvl, lvl,
                                 level=lvl, diff_slice=diff_for_lvl)

        # per-pathway (whole level), if the spec + 8c leading edge are available
        if args.by_pathway:
            import yaml as _yaml
            if args.stress_spec.is_file():
                _spec = _yaml.safe_load(args.stress_spec.read_text())
                _le = Path(cfg["results_dir"]) / "tables" / "08c_pathways" / "08c_pathway_leading_edge.csv"
                if _le.is_file():
                    _gs = pw_plot.load_pathway_genesets(cfg, tissue, _spec, _le)
                    bw = baseline[baseline.get("level", "whole").astype(str) == "whole"] \
                        if "level" in baseline.columns else baseline
                    for pw, genes in _gs.items():
                        mask = pw_plot._in_pathway_mask(bw, genes)
                        df_pw = bw[mask]
                        if df_pw.empty:
                            continue
                        # differential slice restricted to this pathway's LR pairs
                        diff_pw = None
                        if not differential.empty:
                            dmask = pw_plot._in_pathway_mask(differential, genes)
                            diff_pw = differential[dmask]
                            if diff_pw.empty:
                                diff_pw = None
                        for age in sorted(df_pw["age"].astype(str).unique()):
                            _grids_for_slice(df_pw, age,
                                             pdir_root / "07_focal_grids" / "by_pathway" / pw_plot._slug(pw),
                                             f"pw:{pw}", diff_slice=diff_pw)

    print(f"\n{'='*60}\n✓ Phase 8e summary done. Plots in {pdir_root}/")
    for d in (d_overview, d_base, d_diff, d_sr, d_pd):
        n = len(list(d.glob("*.png")))
        print(f"    {d.name}/: {n} figures")
    if args.by_pathway:
        npw = len(list((pdir_root / "06_by_pathway").rglob("*.png")))
        print(f"    06_by_pathway/: {npw} figures")
    nfg = len(list((pdir_root / "07_focal_grids").rglob("*.png")))
    print(f"    07_focal_grids/: {nfg} figures")
    print("=" * 60 + "\n")


def _slugmap(genesets):
    return [k.replace("HALLMARK_", "") for k in genesets]


def _contrast_groups(contrast_name, ref_group):
    """Map a differential contrast name → (test_group, ref_group) for titles.
    Brain/placenta both encode early/late vs relaxed; fall back to (name, ref)."""
    cl = contrast_name.lower()
    if "early_vs_late" in cl or "early_vs_late" in cl:
        return ("Early_Stress", "Late_Stress")
    if "early" in cl:
        return ("Early_Stress", ref_group)
    if "late" in cl:
        return ("Late_Stress", ref_group)
    return (contrast_name, ref_group)


def _build_pathway_contrasts(baseline, differential, ref_group):
    """Primary stress contrasts present in the data: (Early|Late)_Stress vs Relaxed,
    per age where both groups are present. label = EVR/LVR for filenames; name =
    the differential contrast_name (for the differential arm join)."""
    out = []
    if baseline.empty:
        return out
    ga = baseline[["group", "age"]].drop_duplicates()
    by_age = ga.groupby("age")["group"].apply(lambda s: set(s.astype(str)))
    # map test group -> (label, differential contrast prefix)
    specs = [("Early_Stress", "EVR", "early_vs_relaxed"),
             ("Late_Stress", "LVR", "late_vs_relaxed")]
    diff_names = set(differential["contrast_name"].astype(str).unique()) \
        if (differential is not None and not differential.empty) else set()
    for age, groups in by_age.items():
        if ref_group not in groups:
            continue
        for test_g, label, prefix in specs:
            if test_g not in groups:
                continue
            # find the matching differential contrast_name for this age (brain:
            # *_per_age; placenta: *_E12.5/*_E18.5) by prefix
            cname = next((n for n in diff_names
                          if n.startswith(prefix) and (str(age) in n or n.endswith("per_age"))),
                         f"{prefix}_per_age")
            out.append(dict(name=cname, test_group=test_g, ref_group=ref_group,
                            age=str(age), label=label))
    return out


if __name__ == "__main__":
    main()
