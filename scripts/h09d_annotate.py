#!/usr/bin/env python
"""h09d_annotate.py -- annotate the integrated human placenta clusters.

Primary track: marker-majority on the Gunter-Rahman paper scheme (config/
human_placenta_markers.yaml). Mirrors mouse 07d: sc.tl.score_genes on lognorm
-> per-Leiden-cluster mean -> argmax = subtype; subtype -> compartment via the config.

Corroboration track: SingleR (h_run_singler.R) vs Vento-Tormo, collapsed to the 5
compartments. Honest ceiling = compartment-level only (ref is first-trimester, no term
SCT substates / little erythroid). Writes a per-cluster agreement table; does NOT
override marker labels -- it's a check.

Outputs (in place on h09c_integrated.h5ad):
  obs['subtype']            paper-exact label (marker argmax per cluster)
  obs['compartment']        5-way compartment (marker-derived)
  obs['singler_compartment'] SingleR per-cell compartment (corroboration)
  tables: h09d_cluster_annotation.csv (per-cluster subtype/compartment + marker scores
          + SingleR majority compartment + agreement flag)

Usage (from project root):
  uv run python scripts/h09d_annotate.py
  uv run python scripts/h09d_annotate.py --no-singler   # marker-only
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp
import yaml

GSE_DIR = Path("data/human_validation/placenta/gunter_rahman_2025_GSE271976")
MARKERS_YAML = Path("config/human_placenta_markers.yaml")
REF_H5AD = Path("refs/vento_tormo_2018/vento_tormo_2018.h5ad")
RWORKER = "scripts/h_run_singler.R"

# Vento-Tormo author_cell_type / cell_type -> our 5 compartments (compartment-level check)
VT_COMPARTMENT = {
    # trophoblast
    "VCT": "trophoblast", "SCT": "trophoblast", "EVT": "trophoblast",
    "placental villous trophoblast": "trophoblast", "syncytiotrophoblast cell": "trophoblast",
    "extravillous trophoblast": "trophoblast",
    # decidua / stromal
    "dS1": "decidua_stromal", "dS2": "decidua_stromal", "dS3": "decidua_stromal",
    "fFB1": "decidua_stromal", "fFB2": "decidua_stromal",
    "decidual cell": "decidua_stromal", "fibroblast": "decidua_stromal",
    "glandular secretory epithelial cell": "decidua_stromal",
    # vascular
    "endothelial cell": "vascular", "pericyte": "vascular",
    # immune
    "dNK1": "immune", "dNK2": "immune", "dNK3": "immune", "dNK p": "immune",
    "NK CD16+": "immune", "NK CD16-": "immune", "Tcells": "immune", "ILC3": "immune",
    "dM1": "immune", "dM2": "immune", "dM3": "immune", "MO": "immune", "HB": "immune",
    "DC1": "immune", "DC2": "immune", "Plasma": "immune", "Granulocytes": "immune",
    "decidual natural killer cell, human": "immune", "T cell": "immune",
    "macrophage": "immune", "Hofbauer cell": "immune", "monocyte": "immune",
    "plasma cell": "immune", "conventional dendritic cell": "immune",
    "granulocyte": "immune", "innate lymphoid cell": "immune",
    "CD16-positive, CD56-dim natural killer cell, human": "immune",
    "CD16-negative, CD56-bright natural killer cell, human": "immune",
}


def marker_majority(adata, cfg):
    """score_genes per subtype on lognorm layer -> per-cluster mean -> argmax."""
    subs = cfg["subtypes"]
    scores = {}
    for name, spec in subs.items():
        present = [g for g in spec["markers"] if g in adata.var_names]
        missing = [g for g in spec["markers"] if g not in adata.var_names]
        if missing:
            print(f"    [{name}] missing {missing}")
        if len(present) < 2:
            print(f"    [{name}] <2 markers present -- skipped")
            continue
        key = f"_ms_{name}"
        sc.tl.score_genes(adata, gene_list=present, score_name=key,
                          use_raw=False, layer="lognorm")
        scores[name] = adata.obs.pop(key).values
    if not scores:
        sys.exit("no marker sets scored -- check gene symbols")
    score_df = pd.DataFrame(scores, index=adata.obs_names)
    clusters = adata.obs["leiden"].astype(str)
    cluster_mean = score_df.assign(_c=clusters.values).groupby("_c")[list(scores)].mean()
    sub2comp = {s: subs[s]["compartment"] for s in subs}
    cl_subtype = cluster_mean.idxmax(axis=1)
    cl_comp = cl_subtype.map(sub2comp)
    return cl_subtype, cl_comp, cluster_mean


def run_singler(adata, cfg, n_jobs=16):
    """Build query (our lognorm) + Vento-Tormo ref (compartment labels), call SingleR."""
    if not REF_H5AD.is_file():
        print(f"  [warn] {REF_H5AD} not found -- skipping SingleR")
        return None
    ref = sc.read_h5ad(REF_H5AD)
    # map ref var_names (Ensembl) -> symbols
    sym = ref.var["gene_symbols"] if "gene_symbols" in ref.var else ref.var["feature_name"]
    ref.var_names = sym.astype(str).values
    ref.var_names_make_unique()
    # compartment label per ref cell from author_cell_type (fallback cell_type)
    lab = ref.obs["author_cell_type"].astype(str).map(VT_COMPARTMENT)
    lab2 = ref.obs["cell_type"].astype(str).map(VT_COMPARTMENT)
    comp = lab.fillna(lab2)
    keep = comp.notna().values
    ref = ref[keep].copy()
    comp = comp[keep]
    print(f"  SingleR ref: {ref.n_obs} cells across {comp.nunique()} compartments")

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        # query lognorm (genes x cells)
        sio.mmwrite(td / "q.mtx", sp.csr_matrix(adata.layers["lognorm"]).T)
        (td / "q_genes.tsv").write_text("\n".join(adata.var_names) + "\n")
        (td / "q_bc.tsv").write_text("\n".join(adata.obs_names) + "\n")
        # ref lognorm (.X is lognorm in the CELLxGENE object)
        sio.mmwrite(td / "r.mtx", sp.csr_matrix(ref.X).T)
        (td / "r_genes.tsv").write_text("\n".join(ref.var_names) + "\n")
        (td / "r_lab.tsv").write_text("\n".join(comp.astype(str)) + "\n")
        out = td / "singler.tsv"
        cmd = ["Rscript", RWORKER,
               "--query-mtx", str(td / "q.mtx"), "--query-genes", str(td / "q_genes.tsv"),
               "--query-barcodes", str(td / "q_bc.tsv"),
               "--ref-mtx", str(td / "r.mtx"), "--ref-genes", str(td / "r_genes.tsv"),
               "--ref-labels", str(td / "r_lab.tsv"), "--output", str(out),
               "--n-jobs", str(n_jobs)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        print(r.stdout.strip())
        if r.returncode != 0:
            print(f"  [warn] SingleR failed (corroboration skipped):\n{r.stderr[-1500:]}")
            return None
        d = pd.read_csv(out, sep="\t").set_index("barcode")
    return d.reindex(adata.obs_names)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-singler", action="store_true")
    ap.add_argument("--n-jobs", type=int, default=16)
    args = ap.parse_args()

    in_path = GSE_DIR / "h5ad" / "h09c_integrated.h5ad"
    if not in_path.is_file():
        sys.exit(f"missing {in_path} (run h09c first)")
    cfg = yaml.safe_load(MARKERS_YAML.read_text())

    adata = sc.read_h5ad(in_path)
    # h09c left lognorm in .X (post normalize_total+log1p) and raw in layers['counts'].
    # Standardize: lognorm in a named layer so score_genes(layer="lognorm") matches mouse.
    if "lognorm" not in adata.layers:
        adata.layers["lognorm"] = adata.X.copy()
    print(f"[h09d] {adata.n_obs:,} cells, {adata.obs['leiden'].nunique()} clusters")

    print("[h09d] marker-majority (primary)")
    cl_sub, cl_comp, cl_scores = marker_majority(adata, cfg)
    clusters = adata.obs["leiden"].astype(str)
    adata.obs["subtype"] = clusters.map(cl_sub.to_dict()).astype("category")
    adata.obs["compartment"] = clusters.map(cl_comp.to_dict()).astype("category")

    sgl = None
    if not args.no_singler:
        print("[h09d] SingleR corroboration (compartment-level)")
        sgl = run_singler(adata, cfg, n_jobs=args.n_jobs)
        if sgl is not None:
            adata.obs["singler_compartment"] = sgl["singler_label"].values

    # per-cluster annotation table + agreement
    rows = []
    for cl in sorted(clusters.unique(), key=int):
        m = clusters == cl
        row = {"leiden": cl, "n_cells": int(m.sum()),
               "subtype": cl_sub[cl], "compartment": cl_comp[cl]}
        if sgl is not None:
            sm = adata.obs.loc[m.values, "singler_compartment"].value_counts()
            row["singler_compartment"] = sm.index[0] if len(sm) else "NA"
            row["singler_frac"] = round(sm.iloc[0] / m.sum(), 3) if len(sm) else 0.0
            row["agree"] = row["singler_compartment"] == row["compartment"]
        rows.append(row)
    tab = pd.DataFrame(rows)
    tab_dir = GSE_DIR / "tables"
    tab_dir.mkdir(parents=True, exist_ok=True)
    tab.to_csv(tab_dir / "h09d_cluster_annotation.csv", index=False)

    # UMAP by subtype + compartment
    plot_dir = GSE_DIR / "plots" / "h09d"
    plot_dir.mkdir(parents=True, exist_ok=True)
    for col in ["compartment", "subtype"]:
        fig, ax = plt.subplots(figsize=(8, 7))
        sc.pl.umap(adata, color=col, ax=ax, show=False, frameon=False,
                   legend_loc="on data", legend_fontsize=6, size=5)
        for c in ax.collections:
            c.set_rasterized(True)
        fig.tight_layout()
        fig.savefig(plot_dir / f"umap_{col}.png", dpi=140, bbox_inches="tight")
        plt.close(fig)

    if "lognorm" in adata.layers:
        del adata.layers["lognorm"]
    adata.write_h5ad(in_path)

    print(f"\n[h09d] wrote annotation -> {tab_dir/'h09d_cluster_annotation.csv'}")
    print(tab.to_string(index=False))
    if sgl is not None:
        n_agree = tab["agree"].sum()
        print(f"\n[h09d] SingleR compartment agreement: {n_agree}/{len(tab)} clusters")
        print("  (disagreements are clusters to inspect; marker label is primary)")
    print(f"[h09d] compartment counts:\n{adata.obs['compartment'].value_counts().to_string()}")


if __name__ == "__main__":
    main()
