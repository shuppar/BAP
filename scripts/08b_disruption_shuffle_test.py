#!/usr/bin/env python
"""
08b_disruption_shuffle_test.py — k-preserving label-permutation null for the
LOST and GAINED counts observed in the developmental-disruption analysis.

The observed asymmetry (more "LOST" than "GAINED" age-trajectory genes under
stress) needs a real null. A naive null that samples sig sets independently
for each group OVERESTIMATES both LOST and GAINED by ignoring gene-level
overlap (most age-DE genes are sig in multiple groups simultaneously —
"universal" developmental signal). So that null is uninformative.

This script uses a k-preserving null instead. For each (sex × level × celltype):

  1. Compute observed sig sets in Relaxed / Early / Late and the observed
     LOST (R-only) and GAINED (E∩L only) counts.
  2. Per gene, count k_i = number of groups in which it is sig (0, 1, 2, 3).
  3. Generate n_perm null shuffles. In each shuffle, each gene KEEPS its k_i
     but RANDOMIZES which of the 3 groups it is sig in. Universal (k=3) and
     never-sig (k=0) genes are unchanged.
  4. Recompute LOST and GAINED in each null. Report observed vs null mean,
     5-95% range, one-sided p-value, z-score.

Interpretation:
  - obs_lost ≫ null_lost  -> R has a biased share of exclusive-sig genes
                              ("disruption" is real, beyond marginal-rate maths).
  - obs_gained ≫ null_gained -> E and L converge on the same genes more often
                                than chance (convergent stress response).
  - Either or both can hold independently.

Outputs:
  tables/08b_de/08b_disruption_shuffle_test.csv          — one row per slice
  plots/08b_de/summary/shuffle_test/{sex}/{level}.png    — MIRROR bar:
        LEFT=obs_lost (red bar) with null 5-95% whisker overlay,
        RIGHT=obs_gained (blue bar) with null 5-95% whisker overlay.

Significance thresholds (LOCKED, match 08b_de.py): padj<0.05 AND |log2FC|>1.

Usage:
  uv run python scripts/08b_disruption_shuffle_test.py --config config/brain.yaml
  uv run python scripts/08b_disruption_shuffle_test.py --config config/brain.yaml \\
      --n-perm 1000 --n-jobs 16
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _utils import load_config, phase_table_dir, parallel_map


PADJ_THR = 0.05
LFC_THR  = 1.0
GROUPS   = ["Relaxed", "Early_Stress", "Late_Stress"]
NEEDED_COLS = ["contrast", "test_method", "sex", "group_level", "pair",
               "level", "celltype", "gene", "log2FC", "padj"]

# Match the disruption-plot threshold (skip tiny slices to avoid noise)
MIN_GENES = 50           # minimum gene universe per slice to run a shuffle
MIN_OBS_TOTAL = 5        # require obs_lost + obs_gained >= this; else skip


def slugify(s):
    return re.sub(r"[^0-9a-zA-Z]+", "_", str(s)).strip("_").lower()


def safe_fig(fig, out, dpi=140):
    out.parent.mkdir(parents=True, exist_ok=True)
    using_constrained = getattr(fig, "get_constrained_layout",
                                lambda: False)()
    if using_constrained:
        fig.savefig(out, dpi=dpi)
    else:
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def collapse_age_pairs(w):
    """One row per (slice × group × gene), most-significant pair kept."""
    return (w.sort_values("padj")
             .drop_duplicates(
                 ["sex", "level", "celltype", "group_level", "gene"],
                 keep="first"))


# ---------------------------------------------------------------------------
# Worker: shuffle test for one slice
# ---------------------------------------------------------------------------

def shuffle_test_slice(job):
    """Worker. job is a dict packaging one slice's data + parameters.

    Returns a dict of result statistics or None if the slice is too small.

    Null: per-gene k-preserving (each gene keeps its #sig groups but
    randomizes WHICH groups; vectorized via per-row argsort of random scores).
    """
    sex, level, ct = job["sex"], job["level"], job["celltype"]
    rel, early, late = job["rel"], job["early"], job["late"]
    n_perm, seed = job["n_perm"], job["seed"]

    # Build universe and (n, 3) significance matrix
    genes = rel.index.union(early.index).union(late.index)
    n = len(genes)
    if n < MIN_GENES:
        return None
    gene_to_idx = {g: i for i, g in enumerate(genes)}

    sig = np.zeros((n, 3), dtype=bool)
    for j, frame in enumerate([rel, early, late]):
        ix = frame.index.intersection(genes)
        if len(ix) == 0:
            continue
        positions = np.fromiter((gene_to_idx[g] for g in ix), dtype=np.int64,
                                count=len(ix))
        padj_vals = frame.loc[ix, "padj"].values
        lfc_vals  = frame.loc[ix, "log2FC"].values
        is_sig = ((~np.isnan(padj_vals)) & (padj_vals < PADJ_THR)
                  & (~np.isnan(lfc_vals)) & (np.abs(lfc_vals) > LFC_THR))
        sig[positions[is_sig], j] = True

    # Marginal sig counts (informational — preserved exactly under the
    # k-preserving null in expectation, exactly per-gene by construction)
    marg = sig.sum(axis=0)
    r_count, e_count, l_count = int(marg[0]), int(marg[1]), int(marg[2])

    # Per-gene total sig count k_i (0/1/2/3). The null preserves this exactly.
    k_per_gene = sig.sum(axis=1).astype(np.int8)   # (n,)
    n_k0 = int((k_per_gene == 0).sum())
    n_k1 = int((k_per_gene == 1).sum())
    n_k2 = int((k_per_gene == 2).sum())
    n_k3 = int((k_per_gene == 3).sum())

    # ------- All 6 disjoint single-pattern observed counts -------
    # k=1 trio:
    obs_r_only = int((sig[:, 0] & ~sig[:, 1] & ~sig[:, 2]).sum())  # = LOST
    obs_e_only = int((~sig[:, 0] & sig[:, 1] & ~sig[:, 2]).sum())
    obs_l_only = int((~sig[:, 0] & ~sig[:, 1] & sig[:, 2]).sum())
    # k=2 trio:
    obs_re_only = int((sig[:, 0] & sig[:, 1] & ~sig[:, 2]).sum())
    obs_rl_only = int((sig[:, 0] & ~sig[:, 1] & sig[:, 2]).sum())
    obs_el_only = int((~sig[:, 0] & sig[:, 1] & sig[:, 2]).sum())  # = GAINED
    # Universal:
    obs_universal = int((sig[:, 0] & sig[:, 1] & sig[:, 2]).sum())

    # Observed LOST / GAINED (aliases for clarity)
    obs_lost   = obs_r_only
    obs_gained = obs_el_only
    obs_diff   = obs_lost - obs_gained
    obs_total  = obs_lost + obs_gained

    if obs_total < MIN_OBS_TOTAL:
        # Slice too sparse — record but skip nulls
        nan = float("nan")
        return dict(
            sex=sex, level=level, celltype=ct, n_genes=n,
            r_count=r_count, e_count=e_count, l_count=l_count,
            n_k0=n_k0, n_k1=n_k1, n_k2=n_k2, n_k3=n_k3,
            obs_r_only=obs_r_only, obs_e_only=obs_e_only,
            obs_l_only=obs_l_only,
            obs_re_only=obs_re_only, obs_rl_only=obs_rl_only,
            obs_el_only=obs_el_only,
            obs_universal=obs_universal,
            obs_lost=obs_lost, obs_gained=obs_gained, obs_diff=obs_diff,
            obs_ratio=nan,
            null_lost_mean=nan, null_lost_p5=nan, null_lost_p95=nan,
            null_gained_mean=nan, null_gained_p5=nan, null_gained_p95=nan,
            null_diff_mean=nan, null_diff_p5=nan, null_diff_p95=nan,
            z_lost=nan, z_gained=nan, z_diff=nan,
            p_lost=nan, p_lost_dep=nan,
            p_gained=nan, p_gained_dep=nan,
            p_diff=nan,
            p_r_only_enr=nan, p_r_only_dep=nan,
            p_e_only_enr=nan, p_e_only_dep=nan,
            p_l_only_enr=nan, p_l_only_dep=nan,
            p_re_only_enr=nan, p_re_only_dep=nan,
            p_rl_only_enr=nan, p_rl_only_dep=nan,
            p_el_only_enr=nan, p_el_only_dep=nan,
            chi2_k1=nan, p_chi2_k1=nan,
            chi2_k2=nan, p_chi2_k2=nan,
            n_perm=0, reliability="low_signal",
            _null_lost=np.array([]), _null_gained=np.array([]),
        )

    # k-preserving null (vectorized):
    #   For each gene with k_i in {1, 2}, randomly choose which k_i of the 3
    #   group slots are sig. Genes with k_i == 0 or 3 are unchanged.
    #   Implementation: argsort random scores per row, then mark positions
    #   whose rank-in-row is < k_i as sig.
    rng = np.random.default_rng(seed)
    null_lost   = np.empty(n_perm, dtype=np.int64)
    null_gained = np.empty(n_perm, dtype=np.int64)

    k_col = k_per_gene[:, None].astype(np.int32)   # (n, 1) broadcast target
    for k in range(n_perm):
        scores = rng.random((n, 3))
        order  = np.argsort(scores, axis=1)
        ranks  = np.argsort(order, axis=1)        # per-row rank (0..2)
        null_sig = ranks < k_col                  # (n, 3) bool; row has k_i Trues
        null_lost[k]   = int((null_sig[:, 0] & ~null_sig[:, 1] & ~null_sig[:, 2]).sum())
        null_gained[k] = int((~null_sig[:, 0] & null_sig[:, 1] & null_sig[:, 2]).sum())

    null_diff = null_lost - null_gained

    def summarize(obs, null_arr):
        mean = float(null_arr.mean())
        std  = float(null_arr.std(ddof=1)) if n_perm > 1 else float("nan")
        p95  = float(np.percentile(null_arr, 95))
        p5   = float(np.percentile(null_arr, 5))
        z    = ((obs - mean) / std
                if (std and not np.isnan(std) and std > 0) else float("nan"))
        # Two one-sided p-values (Mid-p style with +1 in num/denom)
        p_enr = (1 + int((null_arr >= obs).sum())) / (1 + n_perm)
        p_dep = (1 + int((null_arr <= obs).sum())) / (1 + n_perm)
        return mean, p5, p95, z, p_enr, p_dep

    lost_stats   = summarize(obs_lost,   null_lost)
    gained_stats = summarize(obs_gained, null_gained)
    diff_stats   = summarize(obs_diff,   null_diff)

    # ------- Within-stratum tests (the "is Relaxed special?" question) -------
    # Under k-preserving null, the 3 categories within each stratum should
    # each contain n_k/3 genes. Per-category binomial tests use the analytic
    # null Binom(n_k, 1/3). Chi-square tests goodness-of-fit to uniform.
    from scipy import stats as sps
    def _binom_two_sided(obs, n, p=1.0/3):
        # Returns (p_enrichment, p_depletion); BH later
        if n <= 0:
            return float("nan"), float("nan")
        # Enrichment: P(X >= obs)
        p_enr = float(sps.binom.sf(obs - 1, n, p))
        # Depletion: P(X <= obs)
        p_dep = float(sps.binom.cdf(obs, n, p))
        return p_enr, p_dep

    p_r_only_enr,  p_r_only_dep  = _binom_two_sided(obs_r_only,  n_k1)
    p_e_only_enr,  p_e_only_dep  = _binom_two_sided(obs_e_only,  n_k1)
    p_l_only_enr,  p_l_only_dep  = _binom_two_sided(obs_l_only,  n_k1)
    p_re_only_enr, p_re_only_dep = _binom_two_sided(obs_re_only, n_k2)
    p_rl_only_enr, p_rl_only_dep = _binom_two_sided(obs_rl_only, n_k2)
    p_el_only_enr, p_el_only_dep = _binom_two_sided(obs_el_only, n_k2)

    # Chi-square goodness-of-fit within each stratum (H0: 3 cats equal)
    def _chi2_uniform(counts):
        counts = np.asarray(counts, dtype=float)
        n_total = counts.sum()
        if n_total <= 0:
            return float("nan"), float("nan")
        expected = np.full_like(counts, n_total / len(counts))
        chi2 = float(((counts - expected) ** 2 / expected).sum())
        p = float(sps.chi2.sf(chi2, df=len(counts) - 1))
        return chi2, p

    chi2_k1, p_chi2_k1 = _chi2_uniform([obs_r_only, obs_e_only, obs_l_only])
    chi2_k2, p_chi2_k2 = _chi2_uniform(
        [obs_re_only, obs_rl_only, obs_el_only])

    return dict(
        sex=sex, level=level, celltype=ct, n_genes=n,
        r_count=r_count, e_count=e_count, l_count=l_count,
        n_k0=n_k0, n_k1=n_k1, n_k2=n_k2, n_k3=n_k3,
        # disjoint category counts
        obs_r_only=obs_r_only, obs_e_only=obs_e_only, obs_l_only=obs_l_only,
        obs_re_only=obs_re_only, obs_rl_only=obs_rl_only,
        obs_el_only=obs_el_only,
        obs_universal=obs_universal,
        # primary aliases
        obs_lost=obs_lost, obs_gained=obs_gained, obs_diff=obs_diff,
        obs_ratio=(obs_lost / obs_gained) if obs_gained > 0 else float("inf"),
        # permutation null stats
        null_lost_mean=lost_stats[0], null_lost_p5=lost_stats[1],
        null_lost_p95=lost_stats[2],
        null_gained_mean=gained_stats[0], null_gained_p5=gained_stats[1],
        null_gained_p95=gained_stats[2],
        null_diff_mean=diff_stats[0],
        null_diff_p5=diff_stats[1], null_diff_p95=diff_stats[2],
        z_lost=lost_stats[3], z_gained=gained_stats[3], z_diff=diff_stats[3],
        # permutation p-values (LOST/GAINED only; same as binomial in expectation)
        p_lost=lost_stats[4], p_lost_dep=lost_stats[5],
        p_gained=gained_stats[4], p_gained_dep=gained_stats[5],
        p_diff=diff_stats[4],
        # per-category analytic binomial p-values (6 categories x 2 directions)
        p_r_only_enr=p_r_only_enr,   p_r_only_dep=p_r_only_dep,
        p_e_only_enr=p_e_only_enr,   p_e_only_dep=p_e_only_dep,
        p_l_only_enr=p_l_only_enr,   p_l_only_dep=p_l_only_dep,
        p_re_only_enr=p_re_only_enr, p_re_only_dep=p_re_only_dep,
        p_rl_only_enr=p_rl_only_enr, p_rl_only_dep=p_rl_only_dep,
        p_el_only_enr=p_el_only_enr, p_el_only_dep=p_el_only_dep,
        # within-stratum chi-square goodness-of-fit
        chi2_k1=chi2_k1, p_chi2_k1=p_chi2_k1,
        chi2_k2=chi2_k2, p_chi2_k2=p_chi2_k2,
        n_perm=n_perm,
        reliability="ok",
        _null_lost=null_lost, _null_gained=null_gained,
    )


# ---------------------------------------------------------------------------
# Plotting: MIRROR bar — obs LOST left + obs GAINED right, null whiskers overlaid
# ---------------------------------------------------------------------------

COL_LOST   = "#c0392b"   # red — LOST observed (matches disruption plot)
COL_GAINED = "#2980b9"   # blue — GAINED observed
COL_NULL   = "#7f8c8d"   # gray — null 5-95% range
COL_NULLM  = "black"     # null mean tick


def _sig_marker(p_bh):
    """Return a short significance string for annotation."""
    if not np.isfinite(p_bh):
        return ""
    if p_bh < 0.001:  return "***"
    if p_bh < 0.01:   return "**"
    if p_bh < 0.05:   return "*"
    return "ns"


def _direction_marker(p_enr_bh, p_dep_bh):
    """Return (text, color) for a direction-aware significance label."""
    def stars(p):
        if not np.isfinite(p): return ""
        if p < 0.001: return "***"
        if p < 0.01:  return "**"
        if p < 0.05:  return "*"
        return ""
    s_enr = stars(p_enr_bh)
    s_dep = stars(p_dep_bh)
    if s_enr:
        return f"↑ {s_enr}", COL_LOST
    if s_dep:
        return f"↓ {s_dep}", "#2c3e50"
    return "ns", "0.55"


# Stratum colours: k=1 trio in red/orange/yellow,
# k=2 trio in purple/teal/blue (E∩L is the canonical blue = GAINED).
COL_K1 = {
    "r_only": COL_LOST,
    "e_only": "#e67e22",
    "l_only": "#f1c40f",
}
COL_K2 = {
    "re_only": "#9b59b6",
    "rl_only": "#16a085",
    "el_only": COL_GAINED,
}
COL_LOST_FADED   = "#e6cccc"   # pale red — for depleted LOST (Δ<0)
COL_GAINED_FADED = "#ccd9e6"   # pale blue — for depleted GAINED (Δ<0)
LABEL_MAP = {
    "r_only": "R-only", "e_only": "E-only", "l_only": "L-only",
    "re_only": "R∩E", "rl_only": "R∩L", "el_only": "E∩L",
}


def plot_shuffle(results, summary, out_dir, tissue, min_obs_total=5):
    """Two-panel figure per (sex × level), modelled on the disruption plot.

    Panel A (LEFT, narrower) — mirror bar:
        LEFT half (red):  |Δ_LOST|   going LEFT  (Δ = obs_lost  − null_lost_mean)
        RIGHT half (blue): |Δ_GAINED| going RIGHT (Δ = obs_gained − null_gained_mean)
        Bars are SOLID coloured when Δ > 0 (the direction of disruption-related
        enrichment for LOST, of stress-convergence enrichment for GAINED) and
        FADED coloured when Δ < 0 (depletion in that direction).
        In-bar annotation shows the signed Δ + direction marker (↑/↓/ns).

    Panel B (RIGHT, wider) — within-stratum 6-bar breakdown per cell type:
        For each cell type, 6 mini horizontal bars stacked vertically:
            R-only / E-only / L-only   (k=1 trio, red shades)
            R∩E    / R∩L    / E∩L      (k=2 trio, blue shades)
        Dashed vertical line at n_k1/3 (k=1 expected) and n_k2/3 (k=2 expected).
        Per-bar significance from BH-corrected binomial test (obs ~ Binom(n_k, 1/3)).
    """
    if not results:
        return 0
    sm = summary.set_index(["sex", "level", "celltype"])

    by_slice = {}
    for r in results:
        if r is None:
            continue
        by_slice.setdefault((r["sex"], r["level"]), []).append(r)

    drawn = 0
    for (sex_label, level), rows in by_slice.items():
        rows = [r for r in rows
                if r["reliability"] == "ok"
                and (r["obs_lost"] + r["obs_gained"]) >= min_obs_total]
        if not rows:
            continue

        # Sort by descending |Δ_LOST| (the strongest disruption signal at top).
        def _delta_lost(r):
            d = r["obs_lost"] - r["null_lost_mean"]
            return abs(d) if np.isfinite(d) else 0.0
        rows = sorted(rows, key=lambda r: -_delta_lost(r))
        celltypes = [r["celltype"] for r in rows]
        nrow = len(rows)

        # Disruption-plot-style figure size (matches the proportions of
        # 08b_followup_plots.py:plot_disruption — Panel A : Panel B width
        # 1 : 1.6, height scales with cell-type count).
        fig_w = 16
        fig_h = max(5.0, 0.95 * nrow + 2.4)
        fig, (axA, axB) = plt.subplots(
            1, 2, figsize=(fig_w, fig_h),
            gridspec_kw=dict(width_ratios=[1, 1.6], wspace=0.25),
            constrained_layout=True,
        )

        # ===========================================================
        # Panel A: Mirror bar of |Δ| (LOST left, GAINED right)
        # ===========================================================
        y_pos = np.arange(nrow)
        bar_h = 0.7

        delta_lost   = np.array([r["obs_lost"]   - r["null_lost_mean"]
                                  for r in rows])
        delta_gained = np.array([r["obs_gained"] - r["null_gained_mean"]
                                  for r in rows])

        # Bars (drawn at |Δ|; sign is encoded by color + annotation arrow)
        for i, r in enumerate(rows):
            d_l = delta_lost[i]
            d_g = delta_gained[i]
            c_l = COL_LOST   if d_l >= 0 else COL_LOST_FADED
            c_g = COL_GAINED if d_g >= 0 else COL_GAINED_FADED
            axA.barh(i, -abs(d_l), bar_h, color=c_l,
                     edgecolor="black", lw=0.4)
            axA.barh(i,  abs(d_g), bar_h, color=c_g,
                     edgecolor="black", lw=0.4)

        # x-extent: symmetric, with extra headroom for labels placed
        # OUTSIDE the bars. Labels go to the LEFT of the LOST bar tip and
        # to the RIGHT of the GAINED bar tip. Extra 0.50 * max_n on each
        # side gives room for "Δ=±NNN  ↑ ***" without clipping.
        max_n = float(max(np.abs(delta_lost).max(),
                          np.abs(delta_gained).max(), 1.0))
        label_offset = max_n * 0.02

        for i, r in enumerate(rows):
            d_l = delta_lost[i]
            d_g = delta_gained[i]
            # Significance markers from BH p-values
            try:
                bh_l_enr = sm.loc[(sex_label, level, r["celltype"]), "p_lost_BH"]
                bh_l_dep = sm.loc[(sex_label, level, r["celltype"]), "p_lost_dep_BH"]
                bh_g_enr = sm.loc[(sex_label, level, r["celltype"]), "p_gained_BH"]
                bh_g_dep = sm.loc[(sex_label, level, r["celltype"]), "p_gained_dep_BH"]
            except Exception:
                bh_l_enr = bh_l_dep = bh_g_enr = bh_g_dep = np.nan
            star_l, _ = _direction_marker(bh_l_enr, bh_l_dep)
            star_g, _ = _direction_marker(bh_g_enr, bh_g_dep)

            # LOST label — ALWAYS outside (to the LEFT of the left-bar tip).
            label_l = f"Δ={int(round(d_l)):+d}  {star_l}"
            axA.text(-abs(d_l) - label_offset, i, label_l,
                     ha="right", va="center",
                     fontsize=9, color="black", fontweight="bold")

            # GAINED label — ALWAYS outside (to the RIGHT of the right-bar tip).
            label_g = f"Δ={int(round(d_g)):+d}  {star_g}"
            axA.text(abs(d_g) + label_offset, i, label_g,
                     ha="left", va="center",
                     fontsize=9, color="black", fontweight="bold")

        axA.axvline(0, color="k", lw=0.7)
        axA.set_yticks(y_pos)
        axA.set_yticklabels(celltypes, fontsize=10)
        axA.invert_yaxis()
        # Bump xlim to make room for outside labels (was 1.30 with inside).
        axA.set_xlim(-max_n * 1.65, max_n * 1.65)
        # Symmetric absolute-value tick labels
        xt = axA.get_xticks()
        axA.set_xticklabels([f"{abs(int(x))}" for x in xt])
        axA.set_xlabel("|Δ| = |observed − null mean|  (# age-DE genes)",
                       fontsize=10)
        axA.set_title(
            "A. Asymmetry vs k-preserving null\n"
            "← LOST (R-only)     GAINED (E∩L only) →",
            fontsize=11)
        axA.legend(handles=[
            mpatches.Patch(color=COL_LOST,
                           label="LOST enriched   (Δ>0)"),
            mpatches.Patch(color=COL_LOST_FADED,
                           label="LOST depleted   (Δ<0)"),
            mpatches.Patch(color=COL_GAINED,
                           label="GAINED enriched (Δ>0)"),
            mpatches.Patch(color=COL_GAINED_FADED,
                           label="GAINED depleted (Δ<0)"),
        ], fontsize=7, loc="upper left", bbox_to_anchor=(0.0, -0.10),
           frameon=False, ncol=2)
        axA.spines[["top", "right"]].set_visible(False)

        # ===========================================================
        # Panel B: Within-stratum 6-bar breakdown per cell type
        # ===========================================================
        # Per cell type, 6 horizontal mini-bars stacked vertically inside a
        # band of height ~0.85; 3 red shades on top (k=1 trio), 3 blue shades
        # on bottom (k=2 trio), with a thin grey line separating them.
        cat_order_k1 = ["r_only", "e_only", "l_only"]
        cat_order_k2 = ["re_only", "rl_only", "el_only"]
        cat_color = {**COL_K1, **COL_K2}

        max_obs_B = 0
        for r in rows:
            for c in cat_order_k1 + cat_order_k2:
                v = r.get(f"obs_{c}", 0)
                if v > max_obs_B:
                    max_obs_B = v
        max_obs_B = max(max_obs_B, 1)
        text_pad = max_obs_B * 0.012

        band = 0.85
        n_per_ct = 6
        sub_h = band / n_per_ct

        for i, r in enumerate(rows):
            y0 = i - band / 2
            cats = cat_order_k1 + cat_order_k2
            for j, cat in enumerate(cats):
                sub_y = y0 + (j + 0.5) * sub_h
                obs = r.get(f"obs_{cat}", 0)
                axB.barh(sub_y, obs, sub_h * 0.82, color=cat_color[cat],
                         edgecolor="black", lw=0.25, zorder=3)
                # Per-bar significance
                try:
                    bh_enr = sm.loc[(sex_label, level, r["celltype"]),
                                    f"p_{cat}_enr_BH"]
                    bh_dep = sm.loc[(sex_label, level, r["celltype"]),
                                    f"p_{cat}_dep_BH"]
                except Exception:
                    bh_enr = bh_dep = np.nan
                text, _ = _direction_marker(bh_enr, bh_dep)
                lab = f"{LABEL_MAP[cat]}: {obs}  {text}".rstrip()
                axB.text(obs + text_pad, sub_y, lab,
                         ha="left", va="center", fontsize=7.5,
                         color="0.20", zorder=6)
            # Reference lines (n_k1/3 within k=1 band, n_k2/3 within k=2 band)
            ref_k1 = r["n_k1"] / 3.0
            ref_k2 = r["n_k2"] / 3.0
            y_k1_top = y0 + 0 * sub_h
            y_k1_bot = y0 + 3 * sub_h
            y_k2_top = y0 + 3 * sub_h
            y_k2_bot = y0 + 6 * sub_h
            axB.plot([ref_k1, ref_k1], [y_k1_top, y_k1_bot],
                     color="black", lw=1.1, ls="--", alpha=0.7, zorder=4)
            axB.plot([ref_k2, ref_k2], [y_k2_top, y_k2_bot],
                     color="black", lw=1.1, ls="--", alpha=0.7, zorder=4)
            # Thin separator between k=1 and k=2 within this cell type's band
            axB.axhline(y0 + 3 * sub_h, color="0.85", lw=0.4, zorder=2)
            # Cell-type-level separator (above each cell type, except the first)
            if i > 0:
                axB.axhline(y0, color="0.65", lw=0.6, zorder=2)

        axB.set_yticks(y_pos)
        axB.set_yticklabels([""] * nrow)        # cell-type names already on A
        axB.tick_params(axis="y", length=0)
        axB.invert_yaxis()
        axB.set_xlim(0, max_obs_B * 1.35)
        axB.set_xlabel("# age-DE genes in category", fontsize=10)
        axB.set_title(
            "B. Within-stratum breakdown\n"
            "(top 3 bars per cell type = k=1 trio  •  bottom 3 = k=2 trio  "
            "•  dashed line = n_k/3 expected under null)",
            fontsize=11)
        axB.legend(handles=[
            mpatches.Patch(color=COL_K1["r_only"], label="R-only (=LOST)"),
            mpatches.Patch(color=COL_K1["e_only"], label="E-only"),
            mpatches.Patch(color=COL_K1["l_only"], label="L-only"),
            mpatches.Patch(color=COL_K2["re_only"], label="R∩E"),
            mpatches.Patch(color=COL_K2["rl_only"], label="R∩L"),
            mpatches.Patch(color=COL_K2["el_only"], label="E∩L (=GAINED)"),
        ], fontsize=7, loc="upper left", bbox_to_anchor=(0.0, -0.10),
           frameon=False, ncol=6)
        axB.spines[["top", "right"]].set_visible(False)
        axB.grid(axis="x", color="0.92", lw=0.4, zorder=0)
        axB.set_axisbelow(True)

        # Figure-level title
        fig.suptitle(
            f"{tissue} | sex={sex_label} | level={level} | "
            f"within_group_across_age "
            f"(padj<{PADJ_THR} & |log2FC|>{LFC_THR})\n"
            f"Disruption shuffle test — k-preserving null   "
            f"(↑ = obs > null  ↓ = obs < null   "
            f"* p_BH<0.05  ** <0.01  *** <0.001)",
            fontsize=11)
        # Pool-confound footnote (consistent with disruption plot)
        fig.text(0.5, 0.005,
                 "Pool-confounded contrast — interpret with care.",
                 ha="center", fontsize=7, style="italic", color="0.4")

        out = out_dir / sex_label / f"{slugify(level)}.png"
        safe_fig(fig, out)
        drawn += 1
    return drawn


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--subcluster", default=None)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    print("\n=== 08b disruption shuffle test ===")
    cfg = load_config(args.config)
    tissue = cfg.get("tissue")
    print(f"  Tissue: {tissue}; n_perm={args.n_perm}; n_jobs={args.n_jobs}")
    if args.subcluster:
        print(f"  Subcluster mode: {args.subcluster}")

    table_dir = phase_table_dir(cfg, "08b_de")
    suffix = f"_subcluster_{args.subcluster}" if args.subcluster else ""
    csv_path = table_dir / f"08b_de_results{suffix}.csv"
    if not csv_path.is_file():
        sys.exit(f"ERROR: master CSV not found: {csv_path}")
    print(f"Reading {csv_path} ({csv_path.stat().st_size / 1e6:.1f} MB)...")
    df = pd.read_csv(csv_path, usecols=lambda c: c in NEEDED_COLS,
                     low_memory=False)
    print(f"  {len(df):,} rows loaded.")

    w = df[(df["test_method"] == "Wald")
           & (df["contrast"] == "within_group_across_age")]
    if w.empty:
        sys.exit("No within_group_across_age Wald rows — nothing to test.")
    print(f"  within_group_across_age Wald rows: {len(w):,}")
    w = collapse_age_pairs(w)
    print(f"  After collapsing age-pairs: {len(w):,} rows")

    # Build jobs
    jobs = []
    rng_seed = args.seed
    n_total_slices = 0
    n_skipped = 0
    for (sex_label, level, ct), g in w.groupby(
            ["sex", "level", "celltype"], observed=True):
        n_total_slices += 1
        per_grp = {grp: gg.set_index("gene")[["log2FC", "padj"]]
                   for grp, gg in g.groupby("group_level", observed=True)}
        if not all(k in per_grp for k in GROUPS):
            n_skipped += 1
            continue
        jobs.append(dict(
            sex=sex_label, level=level, celltype=ct,
            rel=per_grp["Relaxed"],
            early=per_grp["Early_Stress"],
            late=per_grp["Late_Stress"],
            n_perm=args.n_perm,
            seed=rng_seed,
        ))
        rng_seed += 1  # distinct seed per slice for reproducibility
    print(f"  Slices to test: {len(jobs)} "
          f"(skipped {n_skipped}/{n_total_slices} missing one of the 3 groups)")

    if not jobs:
        sys.exit("No slices to test.")

    # Run in parallel (CPU-bound -> processes via use_threads=False)
    print(f"  Running {len(jobs)} shuffle tests with {args.n_jobs} workers...")
    results = []
    for job, result, err in parallel_map(
            shuffle_test_slice, jobs,
            n_jobs=args.n_jobs, use_threads=False, desc="shuffle"):
        if err:
            print(f"  [warn] {job['sex']}|{job['level']}|{job['celltype']}: "
                  f"{err[:200]}")
            continue
        if result is None:
            continue
        results.append(result)
    print(f"  Completed: {len(results)} slices.")

    if not results:
        sys.exit("No usable results.")

    # Save CSV (strip the raw null arrays before writing)
    summary_rows = []
    for r in results:
        row = {k: v for k, v in r.items() if not k.startswith("_")}
        row["tissue"] = tissue
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    # Multiple-testing correction (BH within each p-value column)
    try:
        from statsmodels.stats.multitest import multipletests
        bh_pairs = [
            ("p_lost", "p_lost_BH"),
            ("p_lost_dep", "p_lost_dep_BH"),
            ("p_gained", "p_gained_BH"),
            ("p_gained_dep", "p_gained_dep_BH"),
            ("p_diff", "p_diff_BH"),
            # Within-stratum chi-square
            ("p_chi2_k1", "p_chi2_k1_BH"),
            ("p_chi2_k2", "p_chi2_k2_BH"),
        ]
        # Per-category binomial p-values
        for cat in ("r_only", "e_only", "l_only",
                    "re_only", "rl_only", "el_only"):
            bh_pairs.append((f"p_{cat}_enr", f"p_{cat}_enr_BH"))
            bh_pairs.append((f"p_{cat}_dep", f"p_{cat}_dep_BH"))
        for src, dst in bh_pairs:
            if src not in summary.columns:
                continue
            ok = summary[src].notna()
            p_adj = np.full(len(summary), np.nan)
            if ok.any():
                p_adj[ok.values] = multipletests(
                    summary.loc[ok, src].values, method="fdr_bh")[1]
            summary[dst] = p_adj
    except ImportError:
        # Best-effort fallback if statsmodels missing
        for _, dst in [("p_lost", "p_lost_BH"), ("p_gained", "p_gained_BH")]:
            summary[dst] = np.nan

    # Order columns sensibly
    cats = ["r_only", "e_only", "l_only", "re_only", "rl_only", "el_only"]
    cat_cols = []
    for c in cats:
        cat_cols.extend([f"obs_{c}",
                         f"p_{c}_enr", f"p_{c}_enr_BH",
                         f"p_{c}_dep", f"p_{c}_dep_BH"])
    col_order = (
        ["tissue", "sex", "level", "celltype", "n_genes",
         "r_count", "e_count", "l_count",
         "n_k0", "n_k1", "n_k2", "n_k3",
         "obs_lost", "obs_gained", "obs_universal",
         "obs_diff", "obs_ratio"]
        + cat_cols
        + ["chi2_k1", "p_chi2_k1", "p_chi2_k1_BH",
           "chi2_k2", "p_chi2_k2", "p_chi2_k2_BH",
           "null_lost_mean", "null_lost_p5", "null_lost_p95",
           "null_gained_mean", "null_gained_p5", "null_gained_p95",
           "null_diff_mean", "null_diff_p5", "null_diff_p95",
           "z_lost", "z_gained", "z_diff",
           "p_lost", "p_lost_BH",
           "p_lost_dep", "p_lost_dep_BH",
           "p_gained", "p_gained_BH",
           "p_gained_dep", "p_gained_dep_BH",
           "p_diff", "p_diff_BH",
           "n_perm", "reliability"]
    )
    summary = summary[[c for c in col_order if c in summary.columns]]
    summary = summary.sort_values(["sex", "level", "celltype"]).reset_index(drop=True)

    out_csv = table_dir / f"08b_disruption_shuffle_test{suffix}.csv"
    summary.to_csv(out_csv, index=False)
    print(f"\n  Wrote {len(summary):,} rows -> {out_csv}")

    # Plots
    plot_root = (Path(cfg["results_dir"]) / "plots"
                 / ("08b_de" + suffix) / "summary" / "shuffle_test")
    n_drawn = plot_shuffle(results, summary, plot_root, tissue)
    print(f"  Drew {n_drawn} plots -> {plot_root}")

    # Headline (sex=combined, level=whole)
    print("\nHeadline (sex=combined, level=whole):")
    h = summary[(summary["sex"] == "combined")
                & (summary["level"] == "whole")
                & (summary["reliability"] == "ok")]
    if h.empty:
        print("  (no combined×whole slices with sufficient signal)")
    else:
        show = h[["celltype",
                  "obs_lost", "null_lost_mean", "null_lost_p95",
                  "z_lost", "p_lost_BH",
                  "obs_gained", "null_gained_mean", "null_gained_p95",
                  "z_gained", "p_gained_BH"]].copy()
        for c in ("null_lost_mean", "null_lost_p95",
                   "null_gained_mean", "null_gained_p95",
                   "z_lost", "z_gained"):
            if c in show.columns:
                show[c] = show[c].round(1)
        for c in ("p_lost_BH", "p_gained_BH"):
            if c in show.columns:
                show[c] = show[c].apply(
                    lambda x: f"{x:.3g}" if pd.notna(x) else "nan")
        print(show.to_string(index=False))


if __name__ == "__main__":
    main()
