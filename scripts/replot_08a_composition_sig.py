#!/usr/bin/env python
"""
replot_08a_composition_sig.py
─────────────────────────────
Re-render 8a makeup stacked bars WITH per-cell-type propeller FDR markers,
without re-running propeller. Reads only existing 8a outputs.

What it does:
  * Loads obs from the annotated h5ad and joins subcluster_name from 08c
    subcluster objects (same prep as 08a_composition.py).
  * Drops contaminants + unassigned cells (same prep as 08a).
  * Reads 08a_composition_results.csv to look up FDR per cell type per age
    per (sex×level×granularity).
  * Renders the same stacked-bar makeup figures, with markers in the stress
    bars only (Relaxed is the reference):
      * inline `*` / `**` / `***` at slab y-center if slab height ≥ INLINE_THR
      * leader-line callout to the right of the bar if slab is thinner
      * collision avoidance: sweep callouts top-to-bottom, push by MIN_SEP if
        they would overlap
  * Writes new PNGs alongside the originals with a configurable suffix
    (default `_sig`) so the originals stay intact.

Marker → contrast mapping:
  Early_Stress bar  ← contrast name containing "early_vs_relaxed"
  Late_Stress bar   ← contrast name containing "late_vs_relaxed"
  (omnibus + early_vs_late are not used for inline markers — single-stress-vs-
  Relaxed is the cleanest story for a descriptive makeup bar.)

Usage (WS, from /home/poller/BAP-BrainPlacenta/):
  uv run python scripts/replot_08a_composition_sig.py --config config/placenta.yaml
  uv run python scripts/replot_08a_composition_sig.py --config config/brain.yaml

The script is read-only with respect to propeller results — it only writes
new PNGs.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import numpy as np
import pandas as pd
import anndata as ad

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _utils import load_config, iter_strata, unassigned_mask  # noqa: E402

# ───────────────────────────────────────────────────────────────────────────
# Constants mirrored from 08a_composition.py (kept in sync manually; if 08a
# changes these, update here too)
# ───────────────────────────────────────────────────────────────────────────

GROUP_ORDER = ["Relaxed", "Early_Stress", "Late_Stress"]
AGE_ORDER = ["P1", "4W", "3mo", "E12.5", "E18.5"]
CONTAM_PREFIX = "Contamination"

TISSUE_TIERS = {
    "brain": {
        "granularities": {"broad": "celltypist_broad", "class": "celltypist_class"},
        "subtype_base": "celltypist_broad",
        "region_key": "celltypist_region",
        "focal": ["Immune", "OPC/Oligodendrocytes", "Astrocytes/Ependymal"],
    },
    "placenta": {
        "granularities": {"broad": "celltype_majority"},
        "subtype_base": "celltype_majority",
        "region_key": None,
        "focal": ["DSC", "Endothelium", "Myeloid", "NK"],
    },
}

# Marker layout tunables
INLINE_THR = 0.025      # slab height below which we use a leader line
MIN_SEP = 0.030         # min y-spacing between adjacent callouts
LEADER_OFFSET = 0.22    # how far past the bar right edge the callout sits
BAR_WIDTH = 0.7


# ───────────────────────────────────────────────────────────────────────────
# helpers
# ───────────────────────────────────────────────────────────────────────────

def slugify(name: str) -> str:
    return re.sub(r"[^0-9a-zA-Z]+", "_", str(name)).strip("_").lower()


def ordered(values, order):
    vals = list(dict.fromkeys(values))
    return [v for v in order if v in vals] + sorted(v for v in vals if v not in order)


def is_contam(name) -> bool:
    s = str(name)
    return s.startswith(CONTAM_PREFIX) or s == "unresolved"


def read_obs(path, cols=None):
    a = ad.read_h5ad(path, backed="r")
    obs = a.obs if cols is None else a.obs[[c for c in cols if c in a.obs.columns]]
    df = obs.copy()
    a.file.close()
    return df


def stars(fdr):
    """FDR -> significance stars."""
    if pd.isna(fdr):
        return ""
    if fdr < 0.001:
        return "***"
    if fdr < 0.01:
        return "**"
    if fdr < 0.05:
        return "*"
    return ""


def stress_group_for_contrast(cname: str):
    """Map a propeller contrast name to the stress group its bar should mark.

    Only `early_vs_relaxed*` (-> Early_Stress) and `late_vs_relaxed*` (->
    Late_Stress) yield markers. Omnibus 3-group and Early-vs-Late tests are
    skipped — their stars don't have a single bar to land on.
    """
    lc = cname.lower()
    if "early_vs_relaxed" in lc:
        return "Early_Stress"
    if "late_vs_relaxed" in lc:
        return "Late_Stress"
    return None


def build_sig_lookup(res_df: pd.DataFrame, tissue: str, sex_label: str,
                    granularity: str) -> dict:
    """Return dict[(age, stress_group, level, category)] -> fdr."""
    sub = res_df[
        (res_df["tissue"] == tissue)
        & (res_df["sex"].astype(str) == sex_label)
        & (res_df["granularity"].astype(str) == granularity)
    ]
    out = {}
    for _, r in sub.iterrows():
        sg = stress_group_for_contrast(str(r.get("contrast", "")))
        if sg is None or pd.isna(r.get("fdr")):
            continue
        key = (str(r["age"]), sg, str(r["level"]), str(r["category"]))
        # Keep smallest FDR if multiple rows for the same key (shouldn't happen
        # but defensive — propeller writes one row per category per slice)
        cur = out.get(key)
        new = float(r["fdr"])
        if cur is None or new < cur:
            out[key] = new
    return out


# ───────────────────────────────────────────────────────────────────────────
# makeup with significance overlay
# ───────────────────────────────────────────────────────────────────────────

def plot_makeup_with_sig(meta_slice, label_col, title, footnote, out,
                         sig_lookup, level_name="whole"):
    """Render the stacked bar makeup and overlay propeller FDR markers."""
    d = meta_slice
    if d.empty or d[label_col].nunique() == 0:
        return

    ages = ordered(d["age"].astype(str).unique(), AGE_ORDER)
    cats = sorted(d[label_col].astype(str).unique())
    cmap = plt.get_cmap("tab20")
    colors = {c: cmap(i % 20) for i, c in enumerate(cats)}

    fig, axes = plt.subplots(
        1, len(ages),
        figsize=(max(4.0, 3.2 * len(ages)), 5),
        squeeze=False,
    )

    for ax, age in zip(axes[0], ages):
        sub = d[d["age"].astype(str) == age]
        groups = ordered(sub["group"].astype(str).unique(), GROUP_ORDER)
        ct = (pd.crosstab(sub["group"].astype(str), sub[label_col].astype(str))
                .reindex(index=groups, columns=cats, fill_value=0))
        frac = ct.div(ct.sum(axis=1), axis=0).fillna(0)
        bottom = np.zeros(len(groups))
        x = np.arange(len(groups))

        # Track each slab's y-center + height as we draw, keyed (group_idx, cat).
        slab_info = {}
        for c in cats:
            heights = frac[c].values
            for gi in range(len(groups)):
                slab_info[(gi, c)] = (bottom[gi] + heights[gi] / 2.0, heights[gi])
            ax.bar(x, heights, bottom=bottom, color=colors[c], label=c,
                   width=BAR_WIDTH, edgecolor="white", linewidth=0.3)
            bottom += heights

        for i, g in enumerate(groups):
            ax.text(i, 1.01, f"n={int(ct.loc[g].sum()):,}", ha="center",
                    va="bottom", fontsize=7, color="0.3")
        ax.set_xticks(x)
        ax.set_xticklabels([g.replace("_Stress", "") for g in groups])
        ax.set_title(age, fontsize=11, pad=18)  # pad clears the n=NNN labels
        ax.set_ylim(0, 1)
        ax.set_ylabel("fraction of cells" if age == ages[0] else "")
        ax.spines[["top", "right"]].set_visible(False)

        # ────────── significance overlay ──────────
        # Only stress group bars get markers — Relaxed is the reference.
        for gi, g in enumerate(groups):
            if g == "Relaxed":
                continue

            # Gather sig hits in this bar
            hits = []
            for c in cats:
                key = (age, g, level_name, c)
                fdr = sig_lookup.get(key, np.nan)
                star = stars(fdr)
                if not star:
                    continue
                y_center, h = slab_info[(gi, c)]
                if h <= 0:  # category absent in this group
                    continue
                hits.append({"cat": c, "y": y_center, "h": h,
                             "star": star, "color": colors[c]})

            if not hits:
                continue

            # Sort by y for callout collision sweep
            hits.sort(key=lambda r: r["y"])

            x_bar_right = gi + BAR_WIDTH / 2.0
            x_leader_end = gi + BAR_WIDTH / 2.0 + LEADER_OFFSET

            placed_callout_y = []
            for hit in hits:
                if hit["h"] >= INLINE_THR:
                    # Inline marker in the slab — black with white halo
                    txt = ax.text(gi, hit["y"], hit["star"],
                                  ha="center", va="center",
                                  fontsize=9, fontweight="bold", color="black",
                                  zorder=6)
                    txt.set_path_effects(
                        [PathEffects.withStroke(linewidth=2.2, foreground="white")]
                    )
                else:
                    # Leader-line callout, color-matched to slab
                    y_pos = hit["y"]
                    for prev_y in reversed(placed_callout_y):
                        if abs(y_pos - prev_y) < MIN_SEP:
                            y_pos = prev_y + MIN_SEP
                            break
                    placed_callout_y.append(y_pos)
                    ax.plot([x_bar_right, x_leader_end],
                            [hit["y"], y_pos],
                            color=hit["color"], lw=0.6, alpha=0.9,
                            clip_on=False, zorder=4)
                    ax.text(x_leader_end + 0.02, y_pos, hit["star"],
                            ha="left", va="center",
                            fontsize=8.5, fontweight="bold", color="black",
                            clip_on=False, zorder=5)

    # Legend (right of figure, like the original)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[c]) for c in cats]
    fig.legend(handles, cats, loc="center left", bbox_to_anchor=(1.0, 0.5),
               fontsize=7, frameon=False, ncol=1 if len(cats) <= 16 else 2)
    fig.suptitle(title, y=1.03, fontsize=12)

    sig_legend = "(* FDR<0.05  ** <0.01  *** <0.001 — propeller per-donor; markers in stress bars only)"
    full_foot = f"{footnote}\n{sig_legend}"
    fig.text(0.5, -0.08, full_foot, ha="center", fontsize=7, style="italic")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ───────────────────────────────────────────────────────────────────────────
# main
# ───────────────────────────────────────────────────────────────────────────

def main():
    global INLINE_THR
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--results-csv", type=Path, default=None,
                    help="Override 08a_composition_results.csv (default: from results_dir)")
    ap.add_argument("--suffix", default="_sig",
                    help="Suffix appended to output filenames (default '_sig'). "
                         "Pass '' to overwrite originals.")
    ap.add_argument("--inline-threshold", type=float, default=INLINE_THR,
                    help=f"Slab-height cutoff for inline vs leader-line marker "
                         f"(default {INLINE_THR}).")
    args = ap.parse_args()
    INLINE_THR = float(args.inline_threshold)

    cfg = load_config(args.config)
    tissue = cfg["tissue"]
    if tissue not in TISSUE_TIERS:
        sys.exit(f"ERROR: unknown tissue '{tissue}'.")
    tiers = TISSUE_TIERS[tissue]

    results_dir = Path(cfg["results_dir"])
    h5 = results_dir / "h5ad"
    annotated = h5 / "08_annotated" / "all_samples.h5ad"
    if not annotated.is_file():
        sys.exit(f"ERROR: annotated input not found: {annotated}")
    res_csv = (args.results_csv
               or results_dir / "tables" / "08a_composition" / "08a_composition_results.csv")
    if not res_csv.is_file():
        sys.exit(f"ERROR: 8a results CSV not found: {res_csv}")
    print(f"  annotated h5ad : {annotated}")
    print(f"  8a results CSV : {res_csv}")
    print(f"  output suffix  : '{args.suffix}'  (originals preserved)")
    print(f"  inline thr     : {INLINE_THR}")

    # ── Mirror 8a's prep: read obs, join subcluster_name, drop contaminants/unassigned ──
    base_keys = list(dict.fromkeys(list(tiers["granularities"].values())
                                   + [tiers["subtype_base"]]))
    region_key = tiers["region_key"]
    want = ["donor_id", "group", "age", "sex", "pool"] + base_keys + (
        [region_key] if region_key else [])
    meta = read_obs(annotated, cols=want)
    for c in ("donor_id", "group", "age", "sex", "pool"):
        if c in meta.columns:
            meta[c] = meta[c].astype(str)

    has_region = bool(region_key) and region_key in meta.columns \
        and meta[region_key].notna().any()

    sub_base = h5 / "08c_subclustered"
    subname = pd.Series(index=meta.index, dtype="object")
    base_for_focal = tiers["subtype_base"]
    for fl in tiers["focal"]:
        if fl not in set(meta[base_for_focal].astype(str)):
            continue
        p = sub_base / f"{slugify(fl)}.h5ad"
        if not p.is_file():
            continue
        sobs = read_obs(p, cols=["subcluster_name"])
        if "subcluster_name" not in sobs.columns:
            continue
        s = sobs["subcluster_name"].reindex(meta.index)
        subname = subname.where(s.isna(), s)
    meta["subcluster_name"] = subname

    ua_keys = base_keys
    contam_mask = meta["subcluster_name"].notna() & meta["subcluster_name"].map(is_contam)
    ua_mask = unassigned_mask(meta, ua_keys)
    meta = meta.loc[~(contam_mask | ua_mask)].copy()
    print(f"  cleaned obs: {len(meta):,} cells "
          f"(dropped {int(contam_mask.sum()):,} contaminant + {int(ua_mask.sum()):,} unassigned)")

    # ── 8a results ──
    res_df = pd.read_csv(res_csv)
    print(f"  loaded {len(res_df):,} propeller rows  "
          f"(unique contrasts: {sorted(res_df['contrast'].dropna().unique())})")

    sex_strata = iter_strata(cfg, axis="sex")
    if "sex" not in meta.columns:
        sex_strata = [("combined", None)]

    plot_root = results_dir / "plots" / "08a_composition"
    foot = ("Pooled cells, descriptive — propeller (per-donor) does the test. "
            "Contaminants dropped.")

    n_plots = 0
    for sex_label, sex_val in sex_strata:
        m_sex = meta if sex_val is None else meta[meta["sex"] == sex_val]
        if m_sex.empty:
            continue

        sig_broad = build_sig_lookup(res_df, tissue, sex_label, "broad")
        sig_subtype = build_sig_lookup(res_df, tissue, sex_label, "subtype")

        # 1. whole / all_cells (broad)
        out = plot_root / sex_label / "whole" / "all_cells" / f"makeup{args.suffix}.png"
        plot_makeup_with_sig(
            m_sex.assign(_l=m_sex[base_for_focal].astype(str)), "_l",
            f"{tissue} whole — all cells (broad), sex={sex_label}",
            foot, out, sig_broad, level_name="whole")
        n_plots += 1

        # 2. whole / <focal>/<subtype>
        for fl in tiers["focal"]:
            sub = m_sex[(m_sex[base_for_focal].astype(str) == fl)
                        & m_sex["subcluster_name"].notna()]
            if sub.empty:
                continue
            out = plot_root / sex_label / "whole" / slugify(fl) / f"makeup{args.suffix}.png"
            plot_makeup_with_sig(
                sub.assign(_l=sub["subcluster_name"].astype(str)), "_l",
                f"{tissue} whole — {fl} subtypes, sex={sex_label}",
                foot, out, sig_subtype, level_name="whole")
            n_plots += 1

        # 3. region / all_cells (brain only)
        if has_region:
            for r in ordered(m_sex[region_key].dropna().astype(str).unique(), []):
                mr = m_sex[m_sex[region_key].astype(str) == r]
                out = (plot_root / sex_label / "region" / slugify(r)
                       / "all_cells" / f"makeup{args.suffix}.png")
                plot_makeup_with_sig(
                    mr.assign(_l=mr[base_for_focal].astype(str)), "_l",
                    f"{tissue} {r} — all cells (broad), sex={sex_label}",
                    foot, out, sig_broad, level_name=str(r))
                n_plots += 1

    print(f"\n✓ Re-rendered {n_plots} makeup figures with FDR markers.")
    print(f"  Output suffix: '{args.suffix}'  (originals preserved alongside)")
    print(f"  Plot root: {plot_root}")


if __name__ == "__main__":
    main()
