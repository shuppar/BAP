#!/usr/bin/env python
"""h09k_rrho_maps.py -- the actual RRHO rank-rank heatmaps for the 2x2 (were missing).

Reads h09k_rankings.parquet (saved by the patched h09k) and renders, per compartment,
a 2x2 grid of RRHO maps (rows=mouse age, cols=human onset). No DE recompute.

Usage:  uv run python scripts/h09k_rrho_maps.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from h09e_cross_species_rrho import rrho_matrix, OUT_DIR  # noqa: E402

TAB = OUT_DIR / "tables"
PLOT = OUT_DIR / "plots" / "h09k_admati_2x2" / "rrho_maps"
GA_MATCHED = {("E12.5", "eoPE"), ("E18.5", "loPE")}
COMPS = ["trophoblast", "decidua_stromal", "vascular", "immune"]
MOUSE_AGES = ["E12.5", "E18.5"]
HUMAN_ONSETS = ["eoPE", "loPE"]


def main():
    rk = pd.read_parquet(TAB / "h09k_rankings.parquet")
    # {(side, arm, comp): Series(stat indexed by gene)}
    R = {}
    for (side, arm, comp), g in rk.groupby(["side", "arm", "compartment"]):
        R[(side, arm, comp)] = g.set_index("gene")["stat"]

    PLOT.mkdir(parents=True, exist_ok=True)
    summ = pd.read_csv(TAB / "h09k_rrho_2x2_summary.csv")
    cls = {(r.mouse_age, r.human_onset, r.compartment): (r.rrho_class, r.concordance_peak,
                                                         r.empirical_p, r.spearman_r)
           for r in summ.itertuples()}

    for comp in COMPS:
        fig, axes = plt.subplots(len(MOUSE_AGES), len(HUMAN_ONSETS), figsize=(9, 8))
        any_plotted = False
        for i, age in enumerate(MOUSE_AGES):
            for j, onset in enumerate(HUMAN_ONSETS):
                ax = axes[i, j]
                m = R.get(("mouse", age, comp)); h = R.get(("human", onset, comp))
                if m is None or h is None:
                    ax.axis("off"); continue
                mat, cut = rrho_matrix(m, h)
                if mat is None:
                    ax.axis("off"); continue
                any_plotted = True
                im = ax.imshow(mat, origin="lower", aspect="auto", cmap="inferno")
                info = cls.get((age, onset, comp), ("", np.nan, np.nan, np.nan))
                star = " *" if (age, onset) in GA_MATCHED else ""
                ax.set_title(f"mouse {age} x {onset}{star}\n{info[0]} peak={info[1]:.1f} "
                             f"p={info[2]:.0e} r={info[3]:.2f}", fontsize=8)
                ax.set_xlabel(f"{onset} rank", fontsize=7)
                ax.set_ylabel(f"mouse {age} rank", fontsize=7)
                fig.colorbar(im, ax=ax, fraction=0.046, label="-log10 p")
        if not any_plotted:
            plt.close(fig); continue
        fig.suptitle(f"{comp}: cross-species RRHO maps (* = GA-matched diagonal)", fontsize=11)
        fig.tight_layout()
        fig.savefig(PLOT / f"rrho_grid_{comp}.png", dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote rrho_grid_{comp}.png")

    print(f"\n[h09k maps] -> {PLOT}")


if __name__ == "__main__":
    main()
