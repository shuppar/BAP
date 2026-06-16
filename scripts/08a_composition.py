#!/usr/bin/env python
"""
08a_composition.py — Phase 8a: cell-type composition analysis (propeller / speckle).

Tests whether stress changes cell-type PROPORTIONS, across a (sex stratum x
contrast x age x level x granularity) grid, every test run per-donor via
scripts/run_propeller.R (speckle + limma; reused unchanged).

Statistical unit is the ANIMAL (donor_id): composition is per-donor cell-type
counts (one row per pup), the compositional analog of pseudobulk. propeller's
limma empirical-Bayes moderation borrows variance across cell types, which suits
the small n here. Pairwise contrasts -> moderated t-test; 3-group -> ANOVA F.
Confounders (sex, pool) enter as extra design columns AFTER the group columns;
the R worker auto-drops any confounder that is constant in a slice, so single-sex
strata and single-pool slices never produce a rank-deficient design.

GRID
  sex stratum : combined (sex is a covariate) | each sex value (subset; sex drops
                out of the covariate list automatically).
  contrast    : Early-vs-Relaxed, Late-vs-Relaxed (from contrasts.yaml, kind="de"),
                3-group omnibus (CSV only), and — brain only — Early-vs-Late
                (synthesized here; not added to the shared YAML so Phase 8b is
                untouched). Placenta keeps its original set (no Early-vs-Late):
                the two stress arms never share an age (E12.5=Early+Relaxed,
                E18.5=Late+Relaxed), so it is skipped with an announcement.
  age         : per age (group_by from the contrast spec).
  level       : whole tissue; and (brain only, if celltypist_region is present)
                each region. Placenta = whole only.
  granularity : brain  -> broad (celltypist_broad) AND class (celltypist_class)
                          AND subtype (mixed: focal coarse types exploded to
                          subcluster_name, everything else at broad).
                placenta -> broad (celltype_majority) AND subtype.

CONTAMINATION = DROP. Cells flagged Contamination_*/unresolved at the 07d
subcluster stage are removed from BOTH numerator and denominator, once, up front,
so every granularity sees the same cleaned cells. Caveats (documented on outputs):
purity correction applies only to the subclustered types; contaminants are removed
not reassigned, so their true types are mildly undercounted for absolute-baseline
reading (irrelevant to the stress contrast, which is what propeller tests).

Composition needs only obs, never .X — so the whole pipeline runs on one obs
DataFrame (fast; no repeated AnnData slicing).

Usage:
  uv run python scripts/08a_composition.py --config config/brain.yaml
  uv run python scripts/08a_composition.py --config config/placenta.yaml --min-donors 2
  uv run python scripts/08a_composition.py --config config/brain.yaml --rscript /usr/bin/Rscript

Inputs:
  {results_dir}/h5ad/08_annotated/all_samples.h5ad         (labels + denominators)
  {results_dir}/h5ad/08c_subclustered/{slug}.h5ad          (subcluster_name, per focal type)

Outputs:
  {results_dir}/plots/08a_composition/heatmaps/{granularity}/{sex}/{contrast}_{age}.png
  {results_dir}/plots/08a_composition/{sex}/{whole|region/<region>}/{all_cells|<celltype>}/makeup.png
  {results_dir}/tables/08a_composition/08a_composition_results.csv   (master)
"""

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from math import ceil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata as ad

from _utils import (load_config, load_contrasts, phase_table_dir, iter_strata,
                    parallel_map, unassigned_mask)

GROUP_ORDER = ["Relaxed", "Early_Stress", "Late_Stress"]
AGE_ORDER = ["P1", "4W", "3mo", "E12.5", "E18.5"]
CONTAM_PREFIX = "Contamination"

# Per-tissue label tiers. Focal coarse types are the ones subclustered in 07b;
# their value strings must match the base label tier exactly (we guard if not).
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


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def slugify(name: str) -> str:
    return re.sub(r"[^0-9a-zA-Z]+", "_", str(name)).strip("_").lower()


def ordered(values, order):
    vals = list(dict.fromkeys(values))
    return [v for v in order if v in vals] + sorted(v for v in vals if v not in order)


def is_contam(name) -> bool:
    s = str(name)
    return s.startswith(CONTAM_PREFIX) or s == "unresolved"


def read_obs(path, cols=None):
    """Read obs only (backed); return a DataFrame restricted to `cols` if given."""
    a = ad.read_h5ad(path, backed="r")
    obs = a.obs if cols is None else a.obs[[c for c in cols if c in a.obs.columns]]
    df = obs.copy()
    a.file.close()
    return df


def per_donor_counts(df, label_col, covariates):
    """Per-donor cell-type count matrix + sample-level covariates.

    Returns (out_df indexed by donor_id with covariate cols then one integer
    count column per category, list_of_count_columns).
    """
    counts = pd.crosstab(df["donor_id"], df[label_col])
    counts.columns = [str(c) for c in counts.columns]
    cov = (df[["donor_id"] + covariates]
           .drop_duplicates("donor_id").set_index("donor_id"))
    dup = cov.index[cov.index.duplicated()]
    if len(dup):
        raise ValueError(f"covariate(s) {covariates} vary within donor(s) {list(dup)}")
    out = cov.join(counts)
    out[counts.columns] = out[counts.columns].fillna(0).astype(int)
    return out, list(counts.columns)


def run_propeller(cmat, ct_cols, covariates, test_factor, levels, rscript, transform="logit"):
    """Write CSV, call run_propeller.R, read results. levels=None -> omnibus ANOVA."""
    with tempfile.TemporaryDirectory(prefix="propeller_") as td:
        td = Path(td)
        in_csv, out_csv = td / "counts.csv", td / "res.csv"
        cmat.to_csv(in_csv)  # donor_id is the index
        cmd = [
            rscript, "scripts/run_propeller.R",
            "--counts", str(in_csv),
            "--celltypes", ",".join(ct_cols),
            "--covariates", ",".join(covariates),
            "--test", test_factor,
            "--levels", (",".join(levels) if levels else ""),
            "--transform", transform,
            "--out", str(out_csv),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError("propeller R subprocess failed:\n"
                               f"  stdout: {proc.stdout.strip()}\n"
                               f"  stderr: {proc.stderr.strip()[-900:]}")
        return pd.read_csv(out_csv)


def col_ci(df, *cands):
    """First column in df matching any candidate (case-insensitive)."""
    low = {c.lower(): c for c in df.columns}
    for c in cands:
        if c.lower() in low:
            return low[c.lower()]
    return None


def aliased_with(df, factor, cov):
    """True if `cov` is perfectly collinear with `factor` in this slice (either
    direction), which makes a `~ ... + cov + factor` design rank-deficient
    (propeller: 'coefficients not estimable'). Classic case: at P1, Late_Stress
    is Pool3-only, so `pool` is functionally determined by `group`. We drop such
    covariates and flag the group effect as confounded with them, rather than
    crash. Evaluated per donor (covariates are donor-level constants)."""
    dd = df[["donor_id", factor, cov]].drop_duplicates()
    if dd[factor].nunique() < 2 or dd[cov].nunique() < 2:
        return False
    f_determines_c = dd.groupby(factor, observed=True)[cov].nunique().max() <= 1
    c_determines_f = dd.groupby(cov, observed=True)[factor].nunique().max() <= 1
    return bool(f_determines_c or c_determines_f)


# ---------------------------------------------------------------------------
# label construction
# ---------------------------------------------------------------------------

def build_label(meta, granularity, tiers):
    """Return a Series of labels for the requested granularity.

    broad / class -> the corresponding base tier column.
    subtype       -> mixed: subtype_base for everything, focal cells replaced by
                     subcluster_name (contaminants already dropped upstream).
    """
    if granularity in tiers["granularities"]:
        return meta[tiers["granularities"][granularity]].astype(str)
    if granularity == "subtype":
        lab = meta[tiers["subtype_base"]].astype(str).copy()
        foc = meta["subcluster_name"].notna()
        lab[foc] = meta.loc[foc, "subcluster_name"].astype(str)
        return lab
    raise ValueError(f"unknown granularity {granularity}")


def granularities_for(tissue, tiers):
    g = list(tiers["granularities"].keys())  # broad (+ class for brain)
    g.append("subtype")
    return g


# ---------------------------------------------------------------------------
# makeup (descriptive) bars — pooled, faceted by age
# ---------------------------------------------------------------------------

def plot_makeup(meta_slice, label_col, title, footnote, out):
    d = meta_slice
    if d.empty or d[label_col].nunique() == 0:
        return
    ages = ordered(d["age"].astype(str).unique(), AGE_ORDER)
    cats = sorted(d[label_col].astype(str).unique())
    cmap = plt.get_cmap("tab20")
    colors = {c: cmap(i % 20) for i, c in enumerate(cats)}
    fig, axes = plt.subplots(1, len(ages), figsize=(max(4.0, 3.2 * len(ages)), 5),
                             squeeze=False)
    for ax, age in zip(axes[0], ages):
        sub = d[d["age"].astype(str) == age]
        groups = ordered(sub["group"].astype(str).unique(), GROUP_ORDER)
        ct = (pd.crosstab(sub["group"].astype(str), sub[label_col].astype(str))
                .reindex(index=groups, columns=cats, fill_value=0))
        frac = ct.div(ct.sum(axis=1), axis=0).fillna(0)
        bottom = np.zeros(len(groups)); x = np.arange(len(groups))
        for c in cats:
            ax.bar(x, frac[c].values, bottom=bottom, color=colors[c], label=c,
                   width=0.7, edgecolor="white", linewidth=0.3)
            bottom += frac[c].values
        for i, g in enumerate(groups):
            ax.text(i, 1.01, f"n={int(ct.loc[g].sum()):,}", ha="center",
                    va="bottom", fontsize=7, color="0.3")
        ax.set_xticks(x); ax.set_xticklabels([g.replace("_Stress", "") for g in groups])
        ax.set_title(age, fontsize=11); ax.set_ylim(0, 1)
        ax.set_ylabel("fraction of cells" if age == ages[0] else "")
        ax.spines[["top", "right"]].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[c]) for c in cats]
    fig.legend(handles, cats, loc="center left", bbox_to_anchor=(1.0, 0.5),
               fontsize=7, frameon=False, ncol=1 if len(cats) <= 16 else 2)
    fig.suptitle(title, y=1.03, fontsize=12)
    fig.text(0.5, -0.06, footnote, ha="center", fontsize=7, style="italic")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------------------
# change heatmap — rows=category, cols=levels, color=log2(prop_ratio), *=FDR<.05
# ---------------------------------------------------------------------------

def plot_heatmap(df_cell, title, out, vlim=2.0):
    """df_cell: rows for one (sex, granularity, pairwise contrast, age)."""
    if df_cell.empty:
        return
    levels = ordered(df_cell["level"].unique(), ["whole"])
    piv = df_cell.pivot_table(index="category", columns="level",
                              values="log2_prop_ratio", aggfunc="first")
    fdrp = df_cell.pivot_table(index="category", columns="level",
                               values="fdr", aggfunc="first")
    levels = [l for l in levels if l in piv.columns]
    piv = piv.reindex(columns=levels)
    fdrp = fdrp.reindex(index=piv.index, columns=levels)
    piv = piv.sort_index()
    fdrp = fdrp.reindex(index=piv.index)

    nrow, ncol = piv.shape
    fig, ax = plt.subplots(figsize=(max(4.0, 1.1 * ncol + 2.5),
                                    max(2.5, 0.34 * nrow + 1.2)))
    data = piv.values.astype(float)
    im = ax.imshow(np.ma.masked_invalid(data), cmap="RdBu_r",
                   vmin=-vlim, vmax=vlim, aspect="auto")
    ax.set_xticks(range(ncol)); ax.set_xticklabels(levels, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(nrow)); ax.set_yticklabels(piv.index, fontsize=7)
    for i in range(nrow):
        for j in range(ncol):
            v = data[i, j]
            if np.isnan(v):
                ax.text(j, i, "·", ha="center", va="center", fontsize=8, color="0.6")
                continue
            f = fdrp.values[i, j]
            sig = bool(pd.notna(f) and f < 0.05)
            ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6.5,
                    fontweight="bold" if sig else "normal",
                    color="white" if abs(v) > vlim * 0.6 else "black")
            if sig:   # prominent outline = FDR<0.05
                ax.add_patch(plt.Rectangle((j - 0.46, i - 0.46), 0.92, 0.92, fill=False,
                                           edgecolor="black", lw=2.2, zorder=5))
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("log2(prop ratio, stress/ref)", fontsize=8)
    ax.set_title(title, fontsize=10)
    fig.text(0.5, -0.04,
             "Black-outlined cells: propeller FDR<0.05.  · = stratum too thin / category absent.  "
             "Contaminants/unassigned dropped; purity applies to subclustered types only.",
             ha="center", fontsize=6.5, style="italic")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Phase 8a: composition (propeller)")
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--min-donors", type=int, default=None,
                    help="Min donors/group to RUN a stratum. CLI > YAML > 2.")
    ap.add_argument("--reliable-donors", type=int, default=None,
                    help="Donors/group at/above which a stratum is 'ok' (else low_n). "
                         "CLI > YAML > 3.")
    ap.add_argument("--rscript", default=None, help="Path to Rscript (default: PATH)")
    ap.add_argument("--n-jobs", type=int, default=8,
                    help="Concurrent propeller R subprocesses (default 8). Each is a "
                         "lightweight R process (the cost is R startup + speckle load), "
                         "so this parallelizes the slow part. Bump to ~16-24 on the WS.")
    args = ap.parse_args()

    print("\n=== Phase 8a: composition analysis (propeller / speckle) ===")
    cfg = load_config(args.config)
    tissue = cfg["tissue"]
    if tissue not in TISSUE_TIERS:
        sys.exit(f"ERROR: unknown tissue '{tissue}' (expected one of {list(TISSUE_TIERS)}).")
    tiers = TISSUE_TIERS[tissue]

    comp_cfg = cfg.get("composition", {})
    min_donors = (args.min_donors if args.min_donors is not None
                  else int(comp_cfg.get("min_donors", 2)))
    reliable = (args.reliable_donors if args.reliable_donors is not None
                else int(comp_cfg.get("reliable_donors", 3)))
    print(f"  tissue={tissue}  min_donors/group={min_donors}  reliable>={reliable}")

    rscript = args.rscript or shutil.which("Rscript")
    if not rscript:
        sys.exit("ERROR: Rscript not found on PATH. Install R + speckle, or pass --rscript.\n"
                 "  R deps:  BiocManager::install(c('speckle','limma'))  + optparse")
    n_jobs = max(1, args.n_jobs)
    print(f"  Rscript: {rscript}  | propeller workers: {n_jobs}")

    h5 = Path(cfg["results_dir"]) / "h5ad"
    annotated = h5 / "08_annotated" / "all_samples.h5ad"
    if not annotated.is_file():
        sys.exit(f"ERROR: annotated input not found: {annotated}")
    print(f"  Input: {annotated}")

    # ---- assemble the working obs DataFrame --------------------------------
    base_keys = list(tiers["granularities"].values()) + [tiers["subtype_base"]]
    base_keys = list(dict.fromkeys(base_keys))
    region_key = tiers["region_key"]
    want = ["donor_id", "group", "age", "sex", "pool"] + base_keys + (
        [region_key] if region_key else [])
    meta = read_obs(annotated, cols=want)
    missing_core = [c for c in ("donor_id", "group", "age") if c not in meta.columns]
    if missing_core:
        sys.exit(f"ERROR: 08_annotated obs missing required columns: {missing_core}")
    for c in ("donor_id", "group", "age", "sex", "pool"):
        if c in meta.columns:
            meta[c] = meta[c].astype(str)

    has_sex = "sex" in meta.columns and meta["sex"].nunique() > 1
    has_region = bool(region_key) and region_key in meta.columns \
        and meta[region_key].notna().any()
    if region_key and not has_region:
        print(f"  [note] region key '{region_key}' absent/empty -> region levels skipped.")

    # ---- join subcluster_name from each focal type's 08c object ------------
    sub_base = h5 / "08c_subclustered"
    subname = pd.Series(index=meta.index, dtype="object")
    focal_present = []
    base_for_focal = tiers["subtype_base"]
    for fl in tiers["focal"]:
        if fl not in set(meta[base_for_focal].astype(str)):
            print(f"  [note] focal '{fl}' not a value of {base_for_focal}; skipping.")
            continue
        p = sub_base / f"{slugify(fl)}.h5ad"
        if not p.is_file():
            print(f"  [note] no 08c object for focal '{fl}' ({p.name}); subtype view will "
                  f"omit it.")
            continue
        sobs = read_obs(p, cols=["subcluster_name"])
        if "subcluster_name" not in sobs.columns:
            print(f"  [note] {p.name} lacks subcluster_name (run 07d); skipping focal '{fl}'.")
            continue
        s = sobs["subcluster_name"].reindex(meta.index)
        subname = subname.where(s.isna(), s)
        focal_present.append(fl)
    meta["subcluster_name"] = subname
    n_focal_cells = int(meta["subcluster_name"].notna().sum())
    print(f"  focal types with subclusters: {focal_present or '(none)'}  "
          f"({n_focal_cells:,} cells carry subcluster_name)")

    # ---- Record, then DROP non-cell-type cells (contaminants + unassigned) --
    # Both are excluded from the TESTED composition (numerator AND denominator):
    # contaminants are misassigned cells; 'unassigned*' are cells the Phase-7
    # gate couldn't annotate. Neither is a real cell type, so testing them would
    # pollute FDR and distort proportions. But they are NOT lost silently — the
    # per-donor counts/fractions are written to a diagnostic table (see below),
    # so you can see how much mass was dropped and whether it shifts between
    # groups. Masks computed on the full object, before dropping.
    ua_keys = list(dict.fromkeys(list(tiers["granularities"].values())
                                 + [tiers["subtype_base"]]))
    contam_mask = meta["subcluster_name"].notna() & meta["subcluster_name"].map(is_contam)
    ua_mask = unassigned_mask(meta, ua_keys)
    grp_keys = ["donor_id", "age", "group"] + (["sex"] if "sex" in meta.columns else [])
    dropped_diag = (meta.assign(_contam=contam_mask, _unassigned=ua_mask)
                    .groupby(grp_keys, observed=True)
                    .agg(n_total=("_contam", "size"),
                         n_contaminant=("_contam", "sum"),
                         n_unassigned=("_unassigned", "sum")).reset_index())
    dropped_diag["frac_contaminant"] = dropped_diag.n_contaminant / dropped_diag.n_total
    dropped_diag["frac_unassigned"] = dropped_diag.n_unassigned / dropped_diag.n_total
    meta = meta.loc[~(contam_mask | ua_mask)].copy()
    print(f"  dropped {int(contam_mask.sum()):,} contaminant + {int(ua_mask.sum()):,} "
          f"unassigned cells (recorded in diagnostic) -> {len(meta):,} cells remain")

    # ---- contrasts (declarative; early_vs_late_per_age is already in the YAML,
    #      so it inherits here automatically — no synthesis) -------------------
    contrasts = load_contrasts(cfg, kind="de")

    # ---- sex strata (declarative; one definition shared across all 8x stages) --
    sex_strata = iter_strata(cfg, axis="sex")
    if "sex" not in meta.columns:
        sex_strata = [("combined", None)]
        print("  [note] no 'sex' column -> only the combined stratum.")
    print(f"  sex strata: {[s for s, _ in sex_strata]}")

    plot_root = Path(cfg["results_dir"]) / "plots" / "08a_composition"
    table_dir = phase_table_dir(cfg, "08a_composition")
    table_dir.mkdir(parents=True, exist_ok=True)

    # diagnostic: contaminant + unassigned cells that were dropped (per donor)
    dropped_diag.to_csv(table_dir / "08a_dropped_cells_per_donor.csv", index=False)
    print(f"  dropped-cell diagnostic -> {table_dir / '08a_dropped_cells_per_donor.csv'}")

    grans = granularities_for(tissue, tiers)
    rows = []
    jobs = []   # collected first (cheap pandas), then run concurrently below

    for cname, spec in contrasts.items():
        test_raw = spec.get("test", "")
        gb_raw = spec.get("group_by")
        # 8a tests group proportions per age. Skip anything that isn't that:
        # interactions (group:age), across-age (test='age'), and the old list
        # group_by sex-stratified contrast — sex is handled by strata now.
        if (":" in str(test_raw) or test_raw not in ("group", "group_omnibus")
                or isinstance(gb_raw, list)):
            print(f"\n  [skip] {cname}: not a per-age group/omnibus contrast "
                  f"(interaction/across-age/handled by sex strata).")
            continue

        is_omnibus = (test_raw == "group_omnibus")
        test_factor = "group" if is_omnibus else test_raw   # omnibus tests the 'group' column
        levels = spec.get("levels")
        group_by = gb_raw or "age"
        design = spec.get("design", f"~ {test_factor}")
        confound = spec.get("confound_warnings", {}) or {}
        cflag = spec.get("flag")

        # Early-vs-Late is brain-only; placenta has no such contrast in its YAML,
        # but guard anyway in case one is added later.
        if levels and set(map(str, levels)) == {"Early_Stress", "Late_Stress"} \
                and tissue != "brain":
            print(f"\n  [skip] {cname}: Early-vs-Late impossible for {tissue} "
                  f"(stress arms never share an age).")
            continue

        cov_terms = [t.strip() for t in design.replace("~", "").split("+")]
        cov_terms = [t for t in cov_terms if t and t != test_factor and "*" not in t]
        ages = ordered(meta[group_by].astype(str).unique(), AGE_ORDER)
        print(f"\n  Contrast {cname}: test={test_factor} levels={levels or 'omnibus'} "
              f"flag={cflag}")

        for sex_label, sex_val in sex_strata:
            m_sex = meta if sex_val is None else meta[meta["sex"] == sex_val]
            if m_sex.empty:
                continue
            for age in ages:
                m_age = m_sex[m_sex[group_by].astype(str) == age]
                if not is_omnibus:
                    m_age = m_age[m_age[test_factor].astype(str).isin(list(map(str, levels)))]
                if m_age.empty:
                    continue
                # build level list (whole + regions for brain)
                level_specs = [("whole", None)]
                if has_region:
                    for r in ordered(m_age[region_key].dropna().astype(str).unique(), []):
                        level_specs.append((str(r), str(r)))

                for level_name, region_val in level_specs:
                    m_lvl = (m_age if region_val is None
                             else m_age[m_age[region_key].astype(str) == region_val])
                    if m_lvl.empty:
                        continue
                    # per-group donor gate
                    gd = m_lvl.groupby(test_factor, observed=True)["donor_id"].nunique()
                    need_groups = 3 if is_omnibus else 2
                    if len(gd) < need_groups or (gd < min_donors).any():
                        continue
                    reliability = "ok" if gd.min() >= reliable else "low_n"
                    note = confound.get(age, "")
                    if sex_val is not None:
                        note = (note + "; " if note else "") + "sex-specific (low power)"

                    covs = [c for c in cov_terms
                            if c in m_lvl.columns and m_lvl[c].nunique() > 1]
                    # Drop covariates perfectly collinear with the test factor
                    # (e.g. P1: Late_Stress is Pool3-only → pool == group → design
                    # rank-deficient, propeller "coefficients not estimable").
                    # Drop + flag rather than crash; the group effect is then
                    # confounded with the dropped covariate.
                    aliased = [c for c in covs if aliased_with(m_lvl, test_factor, c)]
                    if aliased:
                        covs = [c for c in covs if c not in aliased]
                        amsg = (f"{'/'.join(aliased)} aliased with {test_factor} → "
                                f"dropped from design (confounded_with_pool)")
                        note = (note + "; " if note else "") + amsg
                    covariates = [test_factor] + covs

                    for gran in grans:
                        ml = m_lvl.copy()
                        ml["_label"] = build_label(ml, gran, tiers)
                        if ml["_label"].nunique() < 2:
                            continue
                        cmat, ct_cols = per_donor_counts(ml, "_label", covariates)
                        jobs.append(dict(
                            cname=cname, flag=cflag, sex=sex_label, age=age,
                            level=level_name, granularity=gran,
                            reliability=reliability, note=note,
                            test_factor=test_factor,
                            levels=(None if is_omnibus else list(map(str, levels))),
                            cmat=cmat, ct_cols=ct_cols, covariates=covariates))
                    print(f"    {cname}|{sex_label}|{age}|{level_name}: "
                          f"donors/group {dict(gd)} [{reliability}]")

    # ---- run all propeller jobs concurrently -------------------------------
    # parallel_map (see _utils) handles workers, error capture, and progress.
    # Each job is an independent R subprocess (cost = R startup + speckle load),
    # so threads overlap the slow part. The standard parallel pattern for phases.
    def _run_job(job):
        return run_propeller(job["cmat"], job["ct_cols"], job["covariates"],
                             job["test_factor"], job["levels"], rscript)

    def _meta(job, **extra):
        return dict(tissue=tissue, sex=job["sex"], contrast=job["cname"],
                    flag=job["flag"], age=job["age"], level=job["level"],
                    granularity=job["granularity"], reliability=job["reliability"],
                    note=job["note"], **extra)

    print(f"\n  Running {len(jobs)} propeller jobs across {n_jobs} workers...")
    for job, res, err in parallel_map(_run_job, jobs, n_jobs=n_jobs, desc="propeller"):
        if err:
            last = err.strip().splitlines()[-1] if err.strip() else err
            print(f"    [warn] {job['cname']}|{job['sex']}|{job['age']}|"
                  f"{job['level']}|{job['granularity']}: {last[:160]}")
            r0 = _meta(job, category=None, test_type=None, prop_ratio=None,
                       log2_prop_ratio=None, statistic=None, pvalue=None, fdr=None)
            r0["note"] = (r0["note"] + "; propeller failed").strip("; ")
            rows.append(r0)
            continue
        ratio_col = col_ci(res, "PropRatio")
        fdr_col = col_ci(res, "FDR")
        p_col = col_ci(res, "P.Value", "pvalue", "p_value")
        stat_col = col_ci(res, "Tstatistic", "Fstatistic", "statistic", "t", "f")
        for _, r in res.iterrows():
            ratio = r.get(ratio_col) if ratio_col else None
            log2r = (float(np.log2(ratio))
                     if (ratio is not None and pd.notna(ratio) and ratio > 0)
                     else np.nan)
            rows.append(_meta(
                job, category=r.get("celltype"), test_type=r.get("test_type"),
                prop_ratio=(ratio if ratio_col else None), log2_prop_ratio=log2r,
                statistic=r.get(stat_col) if stat_col else None,
                pvalue=r.get(p_col) if p_col else None,
                fdr=r.get(fdr_col) if fdr_col else None))

    # ---- master CSV --------------------------------------------------------
    res_df = pd.DataFrame(rows)
    out_csv = table_dir / "08a_composition_results.csv"
    res_df.to_csv(out_csv, index=False)
    n_sig = int((res_df["fdr"] < 0.05).sum()) if ("fdr" in res_df and len(res_df)) else 0
    print(f"\n  Master table: {out_csv}  ({len(res_df)} rows, {n_sig} at FDR<0.05)")

    # ---- change heatmaps (pairwise contrasts only) -------------------------
    if len(res_df):
        pair = res_df[res_df["log2_prop_ratio"].notna()]
        for (sex_label, gran, cname, age), grp in pair.groupby(
                ["sex", "granularity", "contrast", "age"], observed=True):
            out = (plot_root / "heatmaps" / gran / sex_label / f"{cname}_{age}.png")
            title = f"{tissue} | {cname} | {age} | {gran} | sex={sex_label}"
            plot_heatmap(grp, title, out)
        print(f"  Heatmaps: {plot_root / 'heatmaps'}")

    # ---- makeup bars (descriptive) -----------------------------------------
    #   informative scopes only: whole x {all_cells + each focal subtype view},
    #   and (brain) region x all_cells (broad). region x subtype is too sparse.
    foot = ("Pooled cells, descriptive — propeller (per-donor) does the test. "
            "Contaminants dropped.")
    broad_key = tiers["subtype_base"]
    for sex_label, sex_val in sex_strata:
        m_sex = meta if sex_val is None else meta[meta["sex"] == sex_val]
        if m_sex.empty:
            continue
        # whole / all_cells (broad makeup)
        plot_makeup(m_sex.assign(_l=m_sex[broad_key].astype(str)), "_l",
                    f"{tissue} whole — all cells (broad), sex={sex_label}", foot,
                    plot_root / sex_label / "whole" / "all_cells" / "makeup.png")
        # whole / <focal type> (subtype makeup within that type)
        for fl in focal_present:
            sub = m_sex[(m_sex[broad_key].astype(str) == fl) & m_sex["subcluster_name"].notna()]
            if sub.empty:
                continue
            plot_makeup(sub.assign(_l=sub["subcluster_name"].astype(str)), "_l",
                        f"{tissue} whole — {fl} subtypes, sex={sex_label}", foot,
                        plot_root / sex_label / "whole" / slugify(fl) / "makeup.png")
        # region / all_cells (broad) — brain only
        if has_region:
            for r in ordered(m_sex[region_key].dropna().astype(str).unique(), []):
                mr = m_sex[m_sex[region_key].astype(str) == r]
                plot_makeup(mr.assign(_l=mr[broad_key].astype(str)), "_l",
                            f"{tissue} {r} — all cells (broad), sex={sex_label}", foot,
                            plot_root / sex_label / "region" / slugify(r) / "all_cells" / "makeup.png")
    print(f"  Makeup bars: {plot_root}")

    print("\n✓ Phase 8a complete.")
    print("  Read 'fdr' with 'reliability' (low_n = underpowered, e.g. sex-specific "
          "or thin regions) and 'note'. Heatmap color = log2(prop ratio), * = FDR<0.05.")
    print("  Caveat: contaminants were dropped (not reassigned); purity correction "
          "applies to subclustered types only — absolute baselines of their true types "
          "are mildly undercounted (irrelevant to the stress contrast).\n")


if __name__ == "__main__":
    main()
