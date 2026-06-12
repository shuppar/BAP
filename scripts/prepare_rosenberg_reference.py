#!/usr/bin/env python
"""
prepare_rosenberg_reference.py — build a labeled P1 reference from Rosenberg 2018.

Rosenberg et al. 2018 (Science, GSE110823) SPLiT-seq: P2 + P11 mouse brain +
spinal cord, 156,049 nuclei, whole CNS. We use the P2 BRAIN subset as the
reference for scANVI label transfer onto our P1 whole-brain samples. Di Bella
(cortex-only) mislabeled ~42% of P1 cells as "erythrocyte" because it has no
classes for cerebellum / midbrain / hindbrain / thalamus etc.; Rosenberg P2
brain covers all those regions.

DESIGN (important): we keep ONLY the ~50 PUBLISHED Rosenberg fine labels as the
reference annotation (`rosenberg_fine`). We do NOT invent coarser labels here.
scANVI transfers `rosenberg_fine` → query `subclass`; the coarser tiers
(class / region / broad) are DERIVED downstream as deterministic groupings of
that transferred label (see config CSVs written by this script). Nothing is
ever trained on a label the authors did not publish.

Input: GSM3017261_150000_CNS_nuclei.mat  (fields: DGE [156049 x 26894 float],
  genes [26894], sample_type [156049], cluster_assignment [156049],
  spinal_cluster_assignment, barcodes)

Steps:
  1. Filter to sample_type == 'p2_brain'           (P2 ≈ P1; drop P11 + spine)
  2. Drop '53 Unresolved' / '54 Unresolved Kcng1'  (junk buckets)
  3. Strip trailing whitespace from genes + labels
  4. Write labeled h5ad: obs['rosenberg_fine'] (raw counts in .X, int)
  5. Emit grouping tables used downstream to derive class/region/broad:
       config/rosenberg_subclass_to_class.csv
       config/rosenberg_subclass_to_region.csv
       config/rosenberg_class_to_broad.csv

Output:
  refs/rosenberg_p2brain_reference.h5ad
  config/rosenberg_subclass_to_class.csv
  config/rosenberg_subclass_to_region.csv
  config/rosenberg_class_to_broad.csv

Usage:
  uv run python scripts/prepare_rosenberg_reference.py \
      --mat refs/GSM3017261_150000_CNS_nuclei.mat \
      --out refs/rosenberg_p2brain_reference.h5ad \
      --config-dir config
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad
import pandas as pd


# ---------------------------------------------------------------------------
# DERIVATION TABLES (fine -> class, fine -> region, class -> broad)
# Written to config/*.csv so they're auditable + editable; downstream reads
# the CSVs, not these dicts.
#
# 'class' is region-tagged (parallels ABC's "19 MB Glut"). 'broad' uses the
# SAME vocabulary as refs/abc_class_to_broad.csv so P1 and 4W/3mo align at the
# broad level. 'region' uses ABC anatomical-division spirit (CTX, CB, TH, HPF,
# OLF, STR, MB) with 'non-regional' for spatially-distributed types
# (glia/immune/vascular/migrating-interneurons).
#
# EXHAUSTIVE over the 72 P2-brain labels. _DROP_ = excluded.
# ---------------------------------------------------------------------------

# fine label -> (class [region-tagged], region)
FINE_TO_CLASS_REGION = {
    # ---- DROP ----
    "53 Unresolved":               ("_DROP_", "_DROP_"),
    "54 Unresolved Kcng1":         ("_DROP_", "_DROP_"),

    # ---- Cortex glutamatergic (CTX) ----
    "5 CTX PyrL2/L3 Pappa2":       ("CTX Glut", "CTX"),
    "6 CTX PyrL2/L3/L4 Ntf3":      ("CTX Glut", "CTX"),
    "7 CTX PyrL2/L3 Met":          ("CTX Glut", "CTX"),
    "8 CTX PyrL4 Wnt5b":           ("CTX Glut", "CTX"),
    "9 CTX PyrL2/L3/L4 Mef2c":     ("CTX Glut", "CTX"),
    "10 CTX PyrL4 Rorb":           ("CTX Glut", "CTX"),
    "11 CTX PyrL4/L5":             ("CTX Glut", "CTX"),
    "12 CTX PyrL5 Itgb3":          ("CTX Glut", "CTX"),
    "13 CTX PyrL5 Fezf2":          ("CTX Glut", "CTX"),
    "14 CTX PyrL6a":               ("CTX Glut", "CTX"),
    "15 CTX PyrL5/L6 Sulf1":       ("CTX Glut", "CTX"),
    "16 CTX PyrL5/L6 Npr3":        ("CTX Glut", "CTX"),
    "17 CTX PyrL6":                ("CTX Glut", "CTX"),
    "18 CLAU Pyr":                 ("CTX Glut", "CTX"),
    "34 SUB Pyr":                  ("HPF Glut", "HPF"),

    # ---- Hippocampal formation (HPF) ----
    "33 HIPP Pyr Cr2":             ("HPF Glut", "HPF"),
    "35 HIPP Pyr Crym":            ("HPF Glut", "HPF"),
    "37 HIPP Pyr Precursor":       ("HPF Glut", "HPF"),
    "38 HIPP Pyr Grik4":           ("HPF Glut", "HPF"),
    "41 HIPP Pyr Npy2r":           ("HPF Glut", "HPF"),
    "36 HIPP Granule Mki67":       ("HPF Glut", "HPF"),
    "39 HIPP Granule Nrp2":        ("HPF Glut", "HPF"),
    "40 HIPP Granule/PyrCA3":      ("HPF Glut", "HPF"),

    # ---- Olfactory (OLF) ----
    "1 OB Mitral/Tufted Eomes":    ("OLF Glut", "OLF"),
    "2 OB Mitral/Tufted Ms4a15":   ("OLF Glut", "OLF"),
    "3 OB Mitral/Tufted Svil":     ("OLF Glut", "OLF"),

    # ---- Thalamus (TH) ----
    "20 THAL Glut":                ("TH Glut", "TH"),
    "21 THAL Int Six3":            ("TH GABA", "TH"),

    # ---- Midbrain / tegmentum / SC (MB) ----
    "19 MTt Glut":                 ("MB Glut", "MB"),
    "42 SC Glut Hmga2":            ("MB Glut", "MB"),
    "30 MD Glyc Int":              ("MB GABA", "MB"),
    "31 MD Int Rxfp2":             ("MB GABA", "MB"),
    "32 Nigral Dopaminergic":      ("MB Dopa", "MB"),

    # ---- Cerebellum (CB) ----
    "25 CB Granule Precursor":     ("CB Glut", "CB"),
    "28 CB Granule":               ("CB Glut", "CB"),
    "22 Purkinje Early":           ("CB GABA", "CB"),
    "23 Purkinje Late":            ("CB GABA", "CB"),
    "24 CB Int Progenitor":        ("CB GABA", "CB"),
    "26 CB Int Stellate/Basket":   ("CB GABA", "CB"),
    "27 CB Int Golgi/Stellate/Basket": ("CB GABA", "CB"),
    "29 CB Int Precursor":         ("CB GABA", "CB"),

    # ---- Striatum (STR) ----
    "4 Medium Spiny Neurons":      ("STR GABA", "STR"),

    # ---- Migrating interneurons (distributed) ----
    "44 Migrating Int Lhx6":       ("Migrating Int", "non-regional"),
    "45 Migrating Int Trdn":       ("Migrating Int", "non-regional"),
    "46 Migrating Int Cpa6":       ("Migrating Int", "non-regional"),
    "47 Migrating Int Foxp2":      ("Migrating Int", "non-regional"),
    "48 Migrating Int Pbx3":       ("Migrating Int", "non-regional"),
    "49 Migrating Int Lgr6":       ("Migrating Int", "non-regional"),
    "50 Migrating Int Adarb2":     ("Migrating Int", "non-regional"),

    # ---- Cajal-Retzius (distributed) ----
    "52 Cajal-Retzius":            ("Cajal-Retzius", "non-regional"),

    # ---- Astrocytes / Ependymal (distributed) ----
    "68 Astro Slc7a10":            ("Astro-Epen", "non-regional"),
    "69 Astro Prdm16":             ("Astro-Epen", "non-regional"),
    "70 Astro Gfap":               ("Astro-Epen", "non-regional"),
    "71 Bergmann Glia":            ("Astro-Epen", "non-regional"),
    "72 Ependyma":                 ("Astro-Epen", "non-regional"),
    "51 SVZ Stem":                 ("Astro-Epen", "non-regional"),

    # ---- OPC / Oligo (distributed) ----
    "61 OPC":                      ("OPC-Oligo", "non-regional"),
    "59 Oligo COP1":               ("OPC-Oligo", "non-regional"),
    "60 Oligo COP2":               ("OPC-Oligo", "non-regional"),
    "58 Oligo NFOL1":              ("OPC-Oligo", "non-regional"),
    "56 Oligo MFOL1":              ("OPC-Oligo", "non-regional"),
    "55 Oligo MFOL2":              ("OPC-Oligo", "non-regional"),
    "57 Oligo MOL":                ("OPC-Oligo", "non-regional"),

    # ---- Immune (distributed) ----
    "63 Microglia":                ("Immune", "non-regional"),
    "62 Macrophage":               ("Immune", "non-regional"),

    # ---- Vascular (distributed) ----
    "64 Endothelia":               ("Vascular", "non-regional"),
    "65 SMC":                      ("Vascular", "non-regional"),
    "66 VLMC Slc6a13":             ("Vascular", "non-regional"),
    "67 VLMC Slc47a1":             ("Vascular", "non-regional"),

    # ---- OEC ----
    "73 OEC":                      ("OEC", "OLF"),
}

# class (region-tagged) -> broad  (broad matches refs/abc_class_to_broad.csv)
CLASS_TO_BROAD = {
    "CTX Glut":        "Excitatory neurons",
    "HPF Glut":        "Excitatory neurons",
    "OLF Glut":        "Excitatory neurons",
    "TH Glut":         "Excitatory neurons",
    "MB Glut":         "Excitatory neurons",
    "CB Glut":         "Excitatory neurons",
    "TH GABA":         "Inhibitory neurons",
    "MB GABA":         "Inhibitory neurons",
    "CB GABA":         "Inhibitory neurons",
    "STR GABA":        "Inhibitory neurons",
    "Migrating Int":   "Inhibitory neurons",
    "Cajal-Retzius":   "Excitatory neurons",
    "MB Dopa":         "Dopaminergic neurons",
    "Astro-Epen":      "Astrocytes/Ependymal",
    "OPC-Oligo":       "OPC/Oligodendrocytes",
    "Immune":          "Immune",
    "Vascular":        "Vascular",
    "OEC":             "Olfactory ensheathing cells",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mat", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--config-dir", type=Path, default=Path("config"))
    ap.add_argument("--min-cells-per-label", type=int, default=10,
                    help="drop fine labels with fewer than this many P2-brain cells")
    args = ap.parse_args()

    if not args.mat.is_file():
        sys.exit(f"ERROR: {args.mat} not found")

    print(f"[rosenberg] loading {args.mat} ...")
    m = sio.loadmat(str(args.mat), verify_compressed_data_integrity=False)

    DGE   = m["DGE"]
    genes = np.array([g.strip() for g in m["genes"].ravel()])
    stype = np.array([s.strip() for s in m["sample_type"].ravel()])
    clst  = np.array([c.strip() for c in m["cluster_assignment"].ravel()])
    print(f"  DGE: {DGE.shape}  genes: {genes.shape}  cells: {stype.shape}")

    # 1. P2 brain only
    keep = stype == "p2_brain"
    print(f"  P2 brain cells: {int(keep.sum())} / {keep.size}")
    X = DGE[keep]
    labels_fine = clst[keep]

    # 2. validate every observed label is mapped (fail loud)
    observed = set(np.unique(labels_fine))
    missing = sorted(l for l in observed if l not in FINE_TO_CLASS_REGION)
    if missing:
        sys.exit("ERROR: P2-brain labels not in FINE_TO_CLASS_REGION:\n  "
                 + "\n  ".join(missing))

    is_drop = np.array([FINE_TO_CLASS_REGION[l][0] == "_DROP_" for l in labels_fine])
    print(f"  dropping {int(is_drop.sum())} cells in _DROP_ buckets (Unresolved)")
    X = X[~is_drop]
    labels_fine = labels_fine[~is_drop]

    # 3. sparse + int
    X = sp.csr_matrix(X)
    if X.dtype.kind == "f":
        X.data = np.rint(X.data)
        X = X.astype(np.int32)

    # 4. build AnnData with ONLY the real fine label
    var = pd.DataFrame(index=pd.Index(genes))
    obs = pd.DataFrame({"rosenberg_fine": pd.Categorical(labels_fine)})
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names_make_unique()

    vc = adata.obs["rosenberg_fine"].value_counts()
    small = vc[vc < args.min_cells_per_label].index.tolist()
    if small:
        print(f"  dropping {len(small)} fine label(s) < {args.min_cells_per_label} "
              f"cells: {small}")
        adata = adata[~adata.obs["rosenberg_fine"].isin(small)].copy()
    adata.obs["rosenberg_fine"] = adata.obs["rosenberg_fine"].cat.remove_unused_categories()

    print(f"\n  Final reference: {adata.n_obs} cells x {adata.n_vars} genes, "
          f"{adata.obs['rosenberg_fine'].nunique()} fine labels")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    print(f"  wrote {args.out}")

    # 5. emit derivation CSVs (only for labels that survived)
    args.config_dir.mkdir(parents=True, exist_ok=True)
    kept = sorted(adata.obs["rosenberg_fine"].cat.categories.tolist())

    sc_rows = [{"rosenberg_fine": f,
                "class": FINE_TO_CLASS_REGION[f][0],
                "region": FINE_TO_CLASS_REGION[f][1]} for f in kept]
    sc_df = pd.DataFrame(sc_rows)
    sc_df[["rosenberg_fine", "class"]].to_csv(
        args.config_dir / "rosenberg_subclass_to_class.csv", index=False)
    sc_df[["rosenberg_fine", "region"]].to_csv(
        args.config_dir / "rosenberg_subclass_to_region.csv", index=False)

    classes = sorted(sc_df["class"].unique())
    cb_missing = [c for c in classes if c not in CLASS_TO_BROAD]
    if cb_missing:
        sys.exit(f"ERROR: classes missing from CLASS_TO_BROAD: {cb_missing}")
    cb_df = pd.DataFrame([{"class": c, "broad": CLASS_TO_BROAD[c]} for c in classes])
    cb_df.to_csv(args.config_dir / "rosenberg_class_to_broad.csv", index=False)

    print(f"\n  wrote derivation tables to {args.config_dir}/:")
    print(f"    rosenberg_subclass_to_class.csv   ({len(sc_df)} fine labels)")
    print(f"    rosenberg_subclass_to_region.csv")
    print(f"    rosenberg_class_to_broad.csv      ({len(cb_df)} classes)")
    print(f"\n  class vocabulary ({len(classes)}): {classes}")
    print(f"  region vocabulary: {sorted(sc_df['region'].unique())}")
    print(f"  broad vocabulary: {sorted(cb_df['broad'].unique())}")


if __name__ == "__main__":
    main()
