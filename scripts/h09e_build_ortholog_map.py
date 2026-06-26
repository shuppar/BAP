#!/usr/bin/env python
"""h09e_build_ortholog_map.py -- mouse<->human 1:1 ortholog table for cross-species RRHO.

Pulls orthologs from Ensembl biomart (mouse mart, human homolog attributes), keeps
high-confidence one-to-one pairs, writes refs/mouse_human_orthologs.tsv. Built once;
reused by h09e (placenta) and the brain RRHO. Cached -- re-run only to refresh.

Usage (from project root, needs network):
  uv run python scripts/h09e_build_ortholog_map.py
"""
import sys
from pathlib import Path

import pandas as pd
from pybiomart import Server

OUT = Path("refs/mouse_human_orthologs.tsv")


def main():
    print("[ortho] connecting to Ensembl biomart...")
    server = Server(host="http://www.ensembl.org")
    mouse = server.marts["ENSEMBL_MART_ENSEMBL"].datasets["mmusculus_gene_ensembl"]

    attrs = [
        "external_gene_name",                       # mouse symbol
        "ensembl_gene_id",                          # mouse ensembl
        "hsapiens_homolog_associated_gene_name",    # human symbol
        "hsapiens_homolog_ensembl_gene",            # human ensembl
        "hsapiens_homolog_orthology_type",          # ortholog_one2one / one2many / many2many
        "hsapiens_homolog_orthology_confidence",    # 1 = high confidence
    ]
    print("[ortho] querying mouse->human homologs (this can take ~1-2 min)...")
    df = mouse.query(attributes=attrs)
    df.columns = ["mouse_symbol", "mouse_ensembl", "human_symbol",
                  "human_ensembl", "orthology_type", "confidence"]

    n_raw = len(df)
    # keep high-confidence one-to-one with both symbols present
    df = df[
        (df["orthology_type"] == "ortholog_one2one")
        & (df["confidence"] == 1)
        & df["mouse_symbol"].notna() & (df["mouse_symbol"] != "")
        & df["human_symbol"].notna() & (df["human_symbol"] != "")
    ].copy()

    # enforce strict 1:1 on symbols (drop any symbol appearing >once on either side)
    for col in ("mouse_symbol", "human_symbol"):
        dup = df[col].value_counts()
        df = df[df[col].isin(dup[dup == 1].index)]

    df = df[["mouse_symbol", "human_symbol", "mouse_ensembl", "human_ensembl"]].reset_index(drop=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False)

    print(f"[ortho] {n_raw} raw homolog rows -> {len(df)} strict 1:1 ortholog pairs")
    print(f"[ortho] wrote {OUT}")
    # spot-check a few canonical placenta/stress genes
    for ms in ["Cga", "Hbb-bs", "Pgr", "Flt1", "Gcm1", "Nr3c1"]:
        hit = df.loc[df["mouse_symbol"] == ms, "human_symbol"]
        print(f"  {ms:10s} -> {hit.iloc[0] if len(hit) else '(no 1:1 ortholog)'}")


if __name__ == "__main__":
    main()
