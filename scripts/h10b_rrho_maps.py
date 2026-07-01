#!/usr/bin/env python
"""h10b_rrho_maps.py -- rank-rank (RRHO) heatmaps for the brain arm, FULL grid.

The map is the GROUND TRUTH the summary label compresses away: it shows whether a high peak
is a clean single-corner hotspot (directional concordance) or both tails lighting up
(strong-but-bidirectional, which the argmax label can't capture).

Reads ONLY saved outputs (rankings parquet + rrho summary) -- never recomputes DE. The cheap
RRHO matrix is recomputed from the saved ranking vectors via the verbatim rrho_matrix.

Layout: one figure per (contrast x level); within it rows = mouse age (P1/4W/3mo),
cols = celltype pair (mouse_ct -> human_ct). Each panel titled with peak / empirical_p /
robust_class. Corner guide (top-left = shared most-UP genes, bottom-right = shared most-DOWN)
printed once per figure.

Usage (Mac or WS, from project root):
  uv run python scripts/h10b_rrho_maps.py --dataset velmeshev
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from h09e_cross_species_rrho import rrho_matrix  # noqa: E402

BRAIN = Path("data/human_validation/brain")
SUB = {"velmeshev": "velmeshev_2019_autism", "maitra": "maitra_2023_GSE213982",
       "nagy": "nagy_2020_GSE144136", "macnair": "macnair_2025_MS"}
AGES = ["P1", "4W", "3mo"]


def load_rankings(rk):
    mouse, human = {}, {}
    for (c, a, l, ct), g in rk[rk.side == "mouse"].groupby(
            ["contrast", "age", "level", "celltype"]):
        mouse[(c, a, l, ct)] = g.set_index("gene")["stat"]
    for ct, g in rk[rk.side == "human"].groupby("celltype"):
        human[ct] = g.set_index("gene")["stat"]
    return mouse, human


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=sorted(SUB))
    args = ap.parse_args()
    tab = BRAIN / SUB[args.dataset] / "tables"
    plot = BRAIN / SUB[args.dataset] / "plots" / "h10b_rrho_maps"
    plot.mkdir(parents=True, exist_ok=True)

    rk = pd.read_parquet(tab / f"h10b_{args.dataset}_rankings.parquet")
    summ = pd.read_csv(tab / f"h10b_{args.dataset}_rrho_summary.csv")
    mouse, human = load_rankings(rk)

    # global color scale = max peak across all cells, so panels are comparable
    vmax = float(summ["concordance_peak"].max())

    for (contrast, level), sub in summ.groupby(["contrast", "level"]):
        cts = sorted(sub[["mouse_ct", "human_ct"]].drop_duplicates()
                     .itertuples(index=False, name=None),
                     key=lambda p: (p[0], p[1]))
        ages = [a for a in AGES if a in sub["mouse_age"].unique()]
        nr, nc = len(ages), max(len(cts), 1)
        fig, axes = plt.subplots(nr, nc, figsize=(2.5 * nc, 2.7 * nr), squeeze=False)

        for i, age in enumerate(ages):
            for j, (mct, hct) in enumerate(cts):
                ax = axes[i][j]
                row = sub[(sub.mouse_age == age) & (sub.mouse_ct == mct)
                          & (sub.human_ct == hct)]
                mk = (contrast, age, level, mct)
                if row.empty or mk not in mouse or hct not in human:
                    ax.axis("off"); continue
                r = row.iloc[0]
                mat, _ = rrho_matrix(mouse[mk], human[hct])
                if mat is None:
                    ax.axis("off"); continue
                ax.imshow(mat, cmap="viridis", origin="upper", aspect="auto",
                          vmin=0, vmax=vmax)
                # title carries the biology-relevant numbers
                lab = (f"{mct}\u2192{hct}\npeak {r['concordance_peak']:.0f}  "
                       f"p={r['empirical_p']:.0e}\n{r['robust_class']}")
                star = " *" if r["age_matched"] else ""
                ax.set_title(lab + star, fontsize=7.5,
                             color=("firebrick" if r["age_matched"] else "black"))
                ax.set_xticks([]); ax.set_yticks([])
                if j == 0:
                    ax.set_ylabel(f"{age}\nmouse rank\nup\u2192down", fontsize=7)
                if i == nr - 1:
                    ax.set_xlabel("human rank\nup\u2192down", fontsize=7)

        # hide any unused trailing axes
        for i in range(nr):
            for j in range(len(cts), nc):
                axes[i][j].axis("off")

        sm = plt.cm.ScalarMappable(cmap="viridis",
                                   norm=plt.Normalize(vmin=0, vmax=vmax))
        cbar = fig.colorbar(sm, ax=axes, fraction=0.015, pad=0.01)
        cbar.set_label("-log10 hypergeometric p", fontsize=8)
        fig.suptitle(
            f"{args.dataset} RRHO maps \u2014 mouse {contrast} \u00d7 {level} "
            f"(top-left=shared UP genes, bottom-right=shared DOWN; * = age-matched)",
            fontsize=10, y=0.99)
        out = plot / f"rrho_grid_{contrast}_{level}.png"
        fig.savefig(out, dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"[h10b_maps] {contrast} x {level}: {len(cts)} celltype cols x "
              f"{len(ages)} ages -> {out}")

    print(f"[h10b_maps] all grids -> {plot}")


if __name__ == "__main__":
    main()
