#!/usr/bin/env python
"""
build_go_relevance.py

Read refs/msigdb_mouse.tsv (GO:BP / GO:MF / GO:CC subcollections only),
apply per-cell-type keyword filter from config/celltype_go_keywords.yaml,
write config/celltype_go_relevance.yaml.

The output is consumed by scripts/08c_pathways_summary.py at label-time:
in M5 panels, only GO terms in the cell type's relevance list are labeled.
All other panels (MH, M2, M8) and all data points are unaffected — this is
purely a label-prioritization layer.

Matching is case-insensitive substring against the full gs_name (e.g.,
keyword 'INFLAMMATORY' matches 'GOBP_REGULATION_OF_INFLAMMATORY_RESPONSE').

Usage:
    uv run python scripts/build_go_relevance.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parent.parent
KEYWORDS_PATH = REPO / "config" / "celltype_go_keywords.yaml"
RELEVANCE_PATH = REPO / "config" / "celltype_go_relevance.yaml"
MSIGDB_PATH = REPO / "refs" / "msigdb_mouse.tsv"

GO_SUBCOLLECTIONS = {"GO:BP", "GO:MF", "GO:CC"}


def make_patterns(keywords: list[str]) -> list[re.Pattern]:
    """Compile keywords into word-boundary regex patterns.

    A keyword matches only if it appears delimited by underscores or term edges,
    so 'T_CELL' matches 'GOBP_T_CELL_ACTIVATION' but NOT the 'crest_cell' chunk
    of 'GOBP_NEURAL_CREST_CELL_FATE_COMMITMENT'. Also 'BMP' matches '_BMP_*'
    but not 'BMPR2', 'WNT' matches '_WNT_*' but not 'WNT5A'.
    """
    return [
        re.compile(rf"(?:^|_){re.escape(kw)}(?:_|$)", re.IGNORECASE)
        for kw in keywords
    ]


def matches_any(term: str, patterns: list[re.Pattern]) -> bool:
    return any(p.search(term) for p in patterns)


def main() -> None:
    if not KEYWORDS_PATH.is_file():
        sys.exit(f"Missing keywords config: {KEYWORDS_PATH}")
    if not MSIGDB_PATH.is_file():
        sys.exit(f"Missing MSigDB file: {MSIGDB_PATH}")

    print(f"Loading {MSIGDB_PATH} ...")
    msigdb = pd.read_csv(
        MSIGDB_PATH,
        sep="\t",
        usecols=["subcollection", "gs_name"],
        dtype=str,
    )
    mask = msigdb["subcollection"].isin(GO_SUBCOLLECTIONS)
    go_terms = sorted(msigdb.loc[mask, "gs_name"].dropna().drop_duplicates().tolist())
    print(f"  {len(go_terms):,} unique GO terms in {sorted(GO_SUBCOLLECTIONS)}")

    by_sub = msigdb.loc[mask].drop_duplicates(["subcollection", "gs_name"])
    counts = by_sub.groupby("subcollection")["gs_name"].nunique().to_dict()
    for sub, n in counts.items():
        print(f"    {sub:8s}: {n:,}")

    with open(KEYWORDS_PATH) as fh:
        cfg = yaml.safe_load(fh)

    output: dict[str, dict[str, dict]] = {}
    for tissue, celltypes in cfg.items():
        output[tissue] = {}
        shared_keywords: list[str] = []
        if "_shared" in celltypes:
            shared_keywords = list(celltypes["_shared"].get("keywords", []))

        print(f"\n[{tissue}]")
        for ct, spec in celltypes.items():
            if ct == "_shared":
                continue
            ct_keywords = list(spec.get("keywords", []))
            kw = ct_keywords + list(shared_keywords)
            patterns = make_patterns(kw)
            matched = [t for t in go_terms if matches_any(t, patterns)]
            output[tissue][ct] = {
                "n_matched": len(matched),
                "keywords_used": kw,
                "go_terms": matched,
            }
            print(f"  {ct:25s}  {len(matched):4d} GO terms matched  "
                  f"({len(ct_keywords)} kw + {len(shared_keywords)} shared)")

    RELEVANCE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(RELEVANCE_PATH, "w") as fh:
        yaml.safe_dump(
            output,
            fh,
            default_flow_style=False,
            sort_keys=False,
            width=2000,
            allow_unicode=True,
        )
    print(f"\nWrote {RELEVANCE_PATH}")
    print("\nReview the output YAML. Trim/add GO terms manually as needed.")
    print("Re-run this script after editing keywords; or hand-edit the relevance "
          "YAML directly (build script won't overwrite manual edits if you skip re-running).")


if __name__ == "__main__":
    main()
