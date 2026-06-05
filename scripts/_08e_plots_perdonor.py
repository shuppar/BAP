"""
_08e_plots_perdonor.py — per-donor quantification plot functions for Phase 8e.

Imported by 08e_communication.py. Not a standalone entry point.

Input: per_donor_df (08e_lr_per_donor.csv) + quant_df (08e_lr_quantified.csv)
Both cover all three group comparisons: ES-v-Rel, LS-v-Rel, ES-v-LS.

Functions:
  plot_per_donor_quantification — two-panel per contrast:
    1. Ranked Δ-activity bar chart (effect size + FDR annotation)
    2. Donor-level stripplot for significant LR pairs (distribution behind stats)
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path


def _slug(s: str) -> str:
    return s.replace(" ", "_").replace("/", "-").replace(".", "")


def plot_per_donor_quantification(per_donor_df, quant_df, ref_group, top_n, pdir):
    """Per-donor quantification plots for all group comparisons in quant_df.

    Called with data already filtered to one age (done in main).
    Generates two plots per contrast (test_group × ctrl_group combination):
      1. Ranked bar of Δ activity (top_n pairs by |Δ|) — effect size figure
      2. Donor-level stripplot for top significant pairs — shows the distribution
         the statistics are based on (key for reviewers)

    All three comparisons are covered (ES-v-Rel, LS-v-Rel, ES-v-LS) because
    quant_df contains all pairwise tests.
    Reproducible from 08e_lr_per_donor.csv + 08e_lr_quantified.csv.
    """
    if per_donor_df.empty or quant_df.empty:
        return

    # Ensure derived columns exist
    pdf = per_donor_df.copy()
    if "activity" not in pdf.columns:
        pdf["activity"] = 1 - pdf["magnitude_rank"]
    if "lr_pair" not in pdf.columns:
        pdf["lr_pair"] = pdf["ligand_complex"] + "→" + pdf["receptor_complex"]
    if "ct_pair" not in pdf.columns:
        pdf["ct_pair"] = pdf["source"] + "→" + pdf["target"]

    groups_present = sorted(pdf["group"].unique())
    group_colors = {g: plt.cm.Set1(i / max(len(groups_present) - 1, 1))
                    for i, g in enumerate(groups_present)}

    # Iterate over all contrast pairs present in quant_df
    contrasts = quant_df[["test_group", "ctrl_group", "contrast"]].drop_duplicates()

    for _, row in contrasts.iterrows():
        test_grp = row["test_group"]
        ctrl_grp = row["ctrl_group"]
        contrast_label = row["contrast"]

        slice_df = quant_df[quant_df["contrast"] == contrast_label].copy()
        if slice_df.empty:
            continue

        slice_df["abs_delta"] = slice_df["delta_activity"].abs()
        top = slice_df.nlargest(top_n, "abs_delta").sort_values("delta_activity")

        if top.empty:
            continue

        # ---- Plot 1: Ranked Δ-activity bar chart ----
        colors = ["#d73027" if d > 0 else "#4575b4" for d in top["delta_activity"]]
        fig, ax = plt.subplots(figsize=(8, max(5, len(top) * 0.35 + 1)))
        ax.barh(range(len(top)), top["delta_activity"], color=colors, alpha=0.85)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(
            [f"{r['lr_pair']}  [{r['ct_pair']}]" for _, r in top.iterrows()],
            fontsize=7
        )
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel(f"Δ activity ({test_grp} − {ctrl_grp};\n1−magnitude_rank)")
        ax.set_title(f"Top {top_n} differential LR pairs\n"
                     f"{test_grp} vs {ctrl_grp}\n"
                     f"(red={test_grp} stronger, blue={ctrl_grp} stronger; * FDR<0.05)")

        for i, (_, r) in enumerate(top.iterrows()):
            if r.get("fdr", 1.0) < 0.05:
                x = r["delta_activity"]
                ax.text(x + (0.01 if x >= 0 else -0.01), i, "*",
                        va="center", ha="left" if x >= 0 else "right",
                        fontsize=10, fontweight="bold")

        fig.tight_layout()
        out = pdir / f"quantified_delta_bar_{_slug(contrast_label)}.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Plot: {out.name}")

        # ---- Plot 2: Donor-level stripplot for significant pairs ----
        sig = slice_df[slice_df.get("significant", False)].nlargest(
            min(12, top_n), "abs_delta")
        if sig.empty:
            # No significant pairs — show top 6 anyway, unlabeled
            sig = slice_df.nlargest(min(6, top_n), "abs_delta")
            labeled_sig = False
        else:
            labeled_sig = True

        n_pairs = len(sig)
        if n_pairs == 0:
            continue

        ncols = min(3, n_pairs)
        nrows = int(np.ceil(n_pairs / ncols))
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(4.5 * ncols, 3.5 * nrows),
                                 constrained_layout=True)
        axes_flat = np.array(axes).flatten() if n_pairs > 1 else [axes]

        for ax, (_, r) in zip(axes_flat, sig.iterrows()):
            lr, ct = r["lr_pair"], r["ct_pair"]
            pair_sub = pdf[(pdf["lr_pair"] == lr) & (pdf["ct_pair"] == ct)]

            for i, grp in enumerate(groups_present):
                scores = pair_sub.loc[pair_sub["group"] == grp, "activity"].values
                if len(scores) == 0:
                    continue
                jx = np.random.default_rng(42).uniform(-0.15, 0.15, len(scores)) + i
                ax.scatter(jx, scores, color=group_colors[grp],
                           s=40, alpha=0.8, zorder=3)
                ax.plot([i - 0.2, i + 0.2],
                        [scores.mean(), scores.mean()],
                        color=group_colors[grp], lw=2.5, zorder=4)

            ax.set_xticks(range(len(groups_present)))
            ax.set_xticklabels([g.replace("_", "\n") for g in groups_present],
                               fontsize=7)
            ax.set_ylabel("activity (1−magnitude_rank)", fontsize=7)
            fdr_val = r.get("fdr", np.nan)
            fdr_str = f"  FDR={fdr_val:.2e}" if labeled_sig and fdr_val < 0.05 else ""
            ax.set_title(f"{lr}\n[{ct}]{fdr_str}", fontsize=7)
            ax.set_ylim(0, 1)

        for ax in axes_flat[n_pairs:]:
            ax.set_visible(False)

        sig_label = "significant pairs" if labeled_sig else "top pairs (none sig.)"
        fig.suptitle(f"Donor-level LR activity: {test_grp} vs {ctrl_grp}\n"
                     f"(line = group mean; {sig_label})", fontsize=9)
        out = pdir / f"quantified_stripplot_{_slug(contrast_label)}.png"
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"  Plot: {out.name}")
