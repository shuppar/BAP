#!/usr/bin/env python
"""h10b_diagnostics.py -- interrogate the brain-RRHO anomaly BEFORE interpretation.

Reads ONLY the saved h10b outputs (rankings parquet + rrho summary) -- never recomputes DE.
Recomputes the cheap RRHO matrix from saved rankings to decompose quadrants.

Three questions (mirrors the h09k_diagnostics discipline):
  1. PEAK-vs-VECTOR-SIZE: is concordance_peak mechanically driven by the SIZE/COHERENCE of
     the mouse DE (P1 has the biggest perinatal DE)? Correlate peak vs mouse-vector strength.
  2. QUADRANT DECOMPOSITION: for the high-peak `discordant` cells, are all four quadrant
     maxima comparable (=> tail artifact, one corner flips the label) or is one quadrant
     genuinely dominant?
  3. P1 vs AGE-MATCHED: is P1 systematically higher across ALL celltypes (mechanical) or
     only in specific cells (biological)?

Usage (Mac or WS, from project root):
  uv run python scripts/h10b_diagnostics.py --dataset velmeshev
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from h09e_cross_species_rrho import rrho_matrix  # noqa: E402

BRAIN = Path("data/human_validation/brain")
SUB = {"velmeshev": "velmeshev_2019_autism", "maitra": "maitra_2023_GSE213982",
       "nagy": "nagy_2020_GSE144136", "macnair": "macnair_2025_MS"}
MOUSE2HUMAN_BROAD = {"Ast": ["Ast"], "ExN": ["ExN"], "InN": ["InN"], "Mic": ["Mic"],
                     "Endo": ["Endo"], "Oli_OPC": ["Oli", "OPC"]}


def vector_strength(s: pd.Series):
    """Descriptors of a ranking vector's DE strength (size + tail coherence)."""
    a = s.to_numpy()
    return {
        "n_genes": len(a),
        "frac_abs_gt2": float(np.mean(np.abs(a) > 2)),
        "mean_top250_abs": float(np.mean(np.sort(np.abs(a))[-250:])),
        "sd_stat": float(np.std(a)),
    }


def quad_decompose(mat):
    """All four quadrant maxima from the RRHO matrix (rows=mouse, cols=human, top=up)."""
    k = mat.shape[0]; h = k // 2
    return {
        "q_concordant_up": float(mat[:h, :h].max()),
        "q_concordant_down": float(mat[h:, h:].max()),
        "q_discordant_tr": float(mat[:h, h:].max()),
        "q_discordant_bl": float(mat[h:, :h].max()),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=sorted(SUB))
    args = ap.parse_args()
    tab = BRAIN / SUB[args.dataset] / "tables"

    rk = pd.read_parquet(tab / f"h10b_{args.dataset}_rankings.parquet")
    summ = pd.read_csv(tab / f"h10b_{args.dataset}_rrho_summary.csv")

    # reconstruct ranking Series
    mouse, human = {}, {}
    for (c, a, l, ct), g in rk[rk.side == "mouse"].groupby(["contrast", "age", "level", "celltype"]):
        mouse[(c, a, l, ct)] = g.set_index("gene")["stat"]
    for ct, g in rk[rk.side == "human"].groupby("celltype"):
        human[ct] = g.set_index("gene")["stat"]

    # ---- 1. peak vs mouse-vector strength ----------------------------------
    print("=" * 70)
    print("1. PEAK vs MOUSE-VECTOR STRENGTH (is peak mechanically size-driven?)")
    print("=" * 70)
    vrows = []
    for (c, a, l, ct), s in mouse.items():
        vs = vector_strength(s)
        vs.update({"contrast": c, "mouse_age": a, "level": l, "mouse_ct": ct})
        vrows.append(vs)
    vdf = pd.DataFrame(vrows)
    print("\nmouse-vector strength by age (mean across cells):")
    print(vdf.groupby("mouse_age")[["n_genes", "frac_abs_gt2", "mean_top250_abs", "sd_stat"]]
          .mean().to_string())

    merged = summ.merge(vdf, on=["contrast", "mouse_age", "level", "mouse_ct"], how="left")
    print("\ncorrelation (across all 50 cells) concordance_peak vs:")
    for col in ["n_genes", "frac_abs_gt2", "mean_top250_abs", "sd_stat"]:
        r = merged[["concordance_peak", col]].corr().iloc[0, 1]
        print(f"  {col:18s}  r = {r:+.3f}")
    print("\n  -> if peak strongly tracks mean_top250_abs/frac_abs_gt2, the high P1 peaks are "
          "mechanical (bigger perinatal DE), not human concordance.")

    # ---- 2. quadrant decomposition of the top cells ------------------------
    print("\n" + "=" * 70)
    print("2. QUADRANT DECOMPOSITION of top-10 peak cells (tail artifact check)")
    print("=" * 70)
    qrows = []
    for _, row in summ.sort_values("concordance_peak", ascending=False).head(10).iterrows():
        mk = (row["contrast"], row["mouse_age"], row["level"], row["mouse_ct"])
        if mk not in mouse or row["human_ct"] not in human:
            continue
        mat, _ = rrho_matrix(mouse[mk], human[row["human_ct"]])
        if mat is None:
            continue
        q = quad_decompose(mat)
        # "label margin": how much the winning quadrant beats the runner-up.
        vals = sorted(q.values(), reverse=True)
        q.update({"cell": f'{row["contrast"]}/{row["mouse_age"]}/{row["level"]}/'
                          f'{row["mouse_ct"]}->{row["human_ct"]}',
                  "class": row["rrho_class"], "peak": row["concordance_peak"],
                  "spearman_r": row["spearman_r"], "margin_top_vs_2nd": round(vals[0] - vals[1], 2)})
        qrows.append(q)
    qdf = pd.DataFrame(qrows)[["cell", "class", "peak", "spearman_r",
                               "q_concordant_up", "q_concordant_down",
                               "q_discordant_tr", "q_discordant_bl", "margin_top_vs_2nd"]]
    pd.set_option("display.width", 200, "display.max_columns", 20)
    print(qdf.to_string(index=False))
    print("\n  -> small margin_top_vs_2nd + flat spearman_r => tail artifact (one corner flips "
          "the label). Large margin + matching spearman sign => real directional structure.")

    # ---- 3. P1 vs others, per celltype -------------------------------------
    print("\n" + "=" * 70)
    print("3. P1 vs 4W vs 3mo, per celltype (mechanical-across-the-board vs cell-specific)")
    print("=" * 70)
    piv = (summ.pivot_table(index=["mouse_ct", "human_ct"], columns="mouse_age",
                            values="concordance_peak", aggfunc="max")
           .reindex(columns=["P1", "4W", "3mo"]))
    print(piv.to_string())
    print("\nmean peak by age (all cells):")
    print(summ.groupby("mouse_age")["concordance_peak"].agg(["mean", "max"]).to_string())
    print("\n  -> P1 highest in EVERY celltype row => mechanical. P1 high only in a few "
          "specific rows => potentially biological (perinatal-specific).")

    # save the merged table for plotting/inspection
    out = tab / f"h10b_{args.dataset}_diagnostics.csv"
    merged.to_csv(out, index=False)
    qdf.to_csv(tab / f"h10b_{args.dataset}_quadrant_decomp.csv", index=False)
    print(f"\n[h10b_diag] merged peak/strength -> {out}")
    print(f"[h10b_diag] quadrant decomp -> {tab / f'h10b_{args.dataset}_quadrant_decomp.csv'}")


if __name__ == "__main__":
    main()
