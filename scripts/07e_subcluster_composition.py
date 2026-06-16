#!/usr/bin/env python
"""
07e_subcluster_composition.py — subcluster (07d-named) composition across stress
groups and ages, two complementary views, plus the tidy source CSV.

Two fractions, two denominators:
  A) frac_within_type  — subtype / all (real) cells of THIS coarse type, per
     (group, age). Stacked bars, sum to 1. "Of the microglia, what's PAM vs homeo?"
  B) frac_of_tissue    — subtype / ALL tissue cells (from 08_annotated), per
     (group, age). Grouped bars. "What % of the whole tissue is this subtype?"

POOLED CELLS, DESCRIPTIVE ONLY. The animal (donor) is the statistical unit; 8a
(propeller, per-donor) does the actual test. These bars are for eyeballing.

Both figures are fully reproducible from the saved CSV (UMAPs / marker dotplots
are NOT — those need the per-cell h5ad).

Run on the workstation (reads 08c objects + the full 08_annotated for totals),
then rsync the PNGs + CSV to Mac.

Usage:
  uv run python scripts/07e_subcluster_composition.py --config config/brain.yaml
  uv run python scripts/07e_subcluster_composition.py --config config/placenta.yaml

Outputs:
  {results_dir}/plots/07b_subcluster/{slug}/subcluster_composition_by_group.png
  {results_dir}/plots/07b_subcluster/{slug}/subcluster_fraction_of_tissue.png
  {results_dir}/tables/07e_subcluster_composition/07e_composition_{tissue}.csv
"""

import argparse
import glob
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata as ad

from _utils import load_config, phase_table_dir

AGE_ORDER = ["P1", "4W", "3mo", "E12.5", "E18.5"]
GROUP_ORDER = ["Relaxed", "Early_Stress", "Late_Stress"]


def slugify(name: str) -> str:
    return re.sub(r"[^0-9a-zA-Z]+", "_", name).strip("_").lower()


def ordered(values, order):
    vals = set(values)
    return [v for v in order if v in vals] + sorted(vals - set(order))


def resolve_cols(obs_cols):
    grp = "group" if "group" in obs_cols else ("condition" if "condition" in obs_cols else None)
    age = "age" if "age" in obs_cols else ("stage" if "stage" in obs_cols else None)
    return grp, age


def load_obs(path, cols):
    a = ad.read_h5ad(path, backed="r")
    grp, age = resolve_cols(a.obs.columns)
    if grp is None or age is None:
        a.file.close()
        return None, False
    extra = [c for c in cols if c in a.obs.columns]
    df = a.obs[extra + [grp, age]].copy()
    a.file.close()
    df = df.rename(columns={grp: "group", age: "age"})
    return df, True


def tissue_totals(annotated_path):
    df, ok = load_obs(annotated_path, cols=[])
    if not ok:
        sys.exit(f"ERROR: {annotated_path} missing group/age columns for tissue totals.")
    return df.groupby(["group", "age"], observed=True).size().rename("n_tissue_total")


def is_contam(name: str) -> bool:
    return name.startswith("Contamination") or name == "unresolved"


def plot_within_type(df, slug, plot_dir):
    d = df[~df["subcluster_name"].map(is_contam)]
    if d.empty:
        print(f"  [skip A] {slug}: only contamination/unresolved"); return
    ages = ordered(d["age"].unique(), AGE_ORDER)
    subtypes = sorted(d["subcluster_name"].unique())
    cmap = plt.get_cmap("tab20"); colors = {s: cmap(i % 20) for i, s in enumerate(subtypes)}
    fig, axes = plt.subplots(1, len(ages), figsize=(max(4.0, 3.2 * len(ages)), 5), squeeze=False)
    for ax, age in zip(axes[0], ages):
        sub = d[d["age"] == age]
        groups = ordered(sub["group"].unique(), GROUP_ORDER)
        ct = (sub.groupby(["group", "subcluster_name"], observed=True)["n_cells"].sum()
                 .unstack(fill_value=0).reindex(index=groups, columns=subtypes, fill_value=0))
        frac = ct.div(ct.sum(axis=1), axis=0).fillna(0)
        bottom = np.zeros(len(groups)); x = np.arange(len(groups))
        for s in subtypes:
            ax.bar(x, frac[s].values, bottom=bottom, color=colors[s], label=s,
                   width=0.7, edgecolor="white", linewidth=0.3)
            bottom += frac[s].values
        for i, g in enumerate(groups):
            ax.text(i, 1.01, f"n={int(ct.loc[g].sum()):,}", ha="center", va="bottom",
                    fontsize=7, color="0.3")
        ax.set_xticks(x); ax.set_xticklabels([g.replace("_Stress", "") for g in groups])
        ax.set_title(age, fontsize=11); ax.set_ylim(0, 1)
        ax.set_ylabel("fraction within type" if age == ages[0] else "")
        ax.spines[["top", "right"]].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[s]) for s in subtypes]
    fig.legend(handles, subtypes, loc="center left", bbox_to_anchor=(1.0, 0.5),
               fontsize=8, frameon=False, title="subtype")
    fig.suptitle(f"{slug}: subtype composition WITHIN type, by group", y=1.03, fontsize=12)
    fig.text(0.5, -0.07, "Pooled cells, descriptive — 8a (per-donor propeller) does the test.",
             ha="center", fontsize=7, style="italic")
    fig.tight_layout()
    fig.savefig(plot_dir / "subcluster_composition_by_group.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {plot_dir/'subcluster_composition_by_group.png'}")


def plot_fraction_of_tissue(df, slug, plot_dir):
    d = df[~df["subcluster_name"].map(is_contam)]
    if d.empty:
        print(f"  [skip B] {slug}: only contamination/unresolved"); return
    ages = ordered(d["age"].unique(), AGE_ORDER)
    subtypes = sorted(d["subcluster_name"].unique())
    cmap = plt.get_cmap("tab20"); colors = {s: cmap(i % 20) for i, s in enumerate(subtypes)}
    fig, axes = plt.subplots(1, len(ages), figsize=(max(4.5, 3.6 * len(ages)), 5),
                             squeeze=False, sharey=True)
    for ax, age in zip(axes[0], ages):
        sub = d[d["age"] == age]
        groups = ordered(sub["group"].unique(), GROUP_ORDER)
        pivot = (sub.groupby(["group", "subcluster_name"], observed=True)["frac_of_tissue"].sum()
                    .unstack(fill_value=0).reindex(index=groups, columns=subtypes, fill_value=0))
        x = np.arange(len(groups)); w = 0.8 / max(1, len(subtypes))
        for j, s in enumerate(subtypes):
            ax.bar(x + (j - (len(subtypes) - 1) / 2) * w, pivot[s].values, width=w,
                   color=colors[s], label=s, edgecolor="white", linewidth=0.2)
        ax.set_xticks(x); ax.set_xticklabels([g.replace("_Stress", "") for g in groups])
        ax.set_title(age, fontsize=11)
        ax.set_ylabel("fraction of all tissue cells" if age == ages[0] else "")
        ax.spines[["top", "right"]].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[s]) for s in subtypes]
    fig.legend(handles, subtypes, loc="center left", bbox_to_anchor=(1.0, 0.5),
               fontsize=8, frameon=False, title="subtype")
    fig.suptitle(f"{slug}: subtype fraction of WHOLE tissue, by group", y=1.03, fontsize=12)
    fig.text(0.5, -0.07, "Pooled cells, descriptive — 8a (per-donor propeller) does the test.",
             ha="center", fontsize=7, style="italic")
    fig.tight_layout()
    fig.savefig(plot_dir / "subcluster_fraction_of_tissue.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {plot_dir/'subcluster_fraction_of_tissue.png'}")


def main():
    ap = argparse.ArgumentParser(description="Phase 7e: subcluster composition + tissue fraction")
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--celltype", default=None, help="One cell type (default: all 08c objects)")
    args = ap.parse_args()

    cfg = load_config(args.config)
    tissue = cfg["tissue"]
    base = Path(cfg["results_dir"]) / "h5ad" / "08c_subclustered"

    annotated = Path(cfg["results_dir"]) / "h5ad" / "08_annotated" / "all_samples.h5ad"
    if not annotated.is_file():
        sys.exit(f"ERROR: full annotated object not found for tissue totals: {annotated}")
    print(f"Tissue totals from {annotated}")
    totals = tissue_totals(annotated)

    if args.celltype:
        paths = [base / f"{slugify(args.celltype)}.h5ad"]
    else:
        paths = sorted(Path(p) for p in glob.glob(str(base / "*.h5ad")))
    if not paths:
        sys.exit(f"No 08c subcluster h5ads found in {base}")

    all_rows = []
    for p in paths:
        if not p.is_file():
            print(f"  [skip] missing {p}"); continue
        slug = p.stem
        obs, ok = load_obs(p, cols=["subcluster_name"])
        if obs is None or "subcluster_name" not in obs.columns:
            print(f"  [skip] {slug}: no subcluster_name / group / age (run 07d first)"); continue

        cnt = (obs.groupby(["subcluster_name", "group", "age"], observed=True)
                  .size().rename("n_cells").reset_index())
        real = cnt[~cnt["subcluster_name"].map(is_contam)]
        type_tot = (real.groupby(["group", "age"], observed=True)["n_cells"].sum()
                        .rename("n_type_total"))
        cnt = cnt.merge(type_tot, on=["group", "age"], how="left")
        cnt = cnt.merge(totals, on=["group", "age"], how="left")
        cnt["frac_within_type"] = cnt["n_cells"] / cnt["n_type_total"]
        cnt["frac_of_tissue"] = cnt["n_cells"] / cnt["n_tissue_total"]
        cnt.insert(0, "coarse_celltype", slug)
        cnt.insert(0, "tissue", tissue)
        all_rows.append(cnt)

        plot_dir = Path(cfg["results_dir"]) / "plots" / "07b_subcluster" / slug
        plot_dir.mkdir(parents=True, exist_ok=True)
        print(f"{slug}:")
        plot_within_type(cnt, slug, plot_dir)
        plot_fraction_of_tissue(cnt, slug, plot_dir)

    if not all_rows:
        sys.exit("No cell types produced rows — check 07d ran.")

    out = pd.concat(all_rows, ignore_index=True)
    col_order = ["tissue", "coarse_celltype", "subcluster_name", "age", "group",
                 "n_cells", "n_type_total", "n_tissue_total",
                 "frac_within_type", "frac_of_tissue"]
    out = out[col_order].sort_values(["coarse_celltype", "age", "group", "subcluster_name"])
    table_dir = phase_table_dir(cfg, "07e_subcluster_composition")
    csv_path = table_dir / f"07e_composition_{tissue}.csv"
    out.to_csv(csv_path, index=False)
    print(f"\n  Source CSV (both plots reproducible from this): {csv_path}")
    print(f"  rows={len(out)}  cell types={out['coarse_celltype'].nunique()}")


if __name__ == "__main__":
    main()
