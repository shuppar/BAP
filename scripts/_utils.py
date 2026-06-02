"""
_utils.py — shared helpers used by phase scripts (01_validate.py … 0N_*.py)
and by notebooks. NOT a phase entry point; the leading underscore marks that.

Import pattern (inside scripts/):
    from _utils import load_config, add_lognorm, phase_paths, select_accelerator

Import pattern (from notebooks at repo root):
    import sys; sys.path.insert(0, "scripts")
    from _utils import add_lognorm

This file is deliberately small. The project conventions (INSTRUCTIONS.md) say
"plain Python scripts, not a package" — so helpers live here as flat functions,
not as a `src/` Python package with subpackages and __init__ trees.
"""

from __future__ import annotations

import sys
from pathlib import Path

import yaml


# ============================================================================
# Config loading
# ============================================================================

def load_config(path: Path) -> dict:
    """Load YAML config and apply the dev.yaml indirection conventions.

    Supports:
      - `samples_from: <other.yaml>`  — pulls the samples list from another file.
      - `subset.enabled: true` + `subset.sample_ids: [...]`  — filters samples
        to that ID list (order preserved from the source manifest).

    Also resolves any relative `h5` paths in samples to absolute paths,
    rooted at the current working directory (= repo root if you run scripts
    from there, as convention).

    Returns a plain dict — no schema validation. Failures surface at use site.
    """
    path = Path(path)
    with path.open() as f:
        cfg = yaml.safe_load(f)

    # --- Indirection: pull records from another YAML if requested ---
    # `samples_from` brings in the sample list AND any tissue-level config blocks
    # the dev config doesn't redefine (annotation, reference). This keeps
    # dev.yaml minimal: it inherits brain.yaml's reference/annotation config
    # rather than having to duplicate it. Locally-defined keys win over inherited.
    if "samples_from" in cfg:
        src_path = Path(cfg["samples_from"])
        with src_path.open() as f:
            src_cfg = yaml.safe_load(f)
        if "samples" not in src_cfg:
            sys.exit(f"ERROR: samples_from={src_path} has no 'samples' key")
        cfg["samples"] = src_cfg["samples"]
        # Inherit tissue-level blocks the dev config didn't override.
        for key in ("annotation", "reference", "contrasts", "stress_focused_cell_types",
                    "composition", "pathways"):
            if key not in cfg and key in src_cfg:
                cfg[key] = src_cfg[key]

    if "samples" not in cfg:
        sys.exit(f"ERROR: config {path} has no 'samples' (and no samples_from)")

    # --- Subset by explicit ID list (dev mode) ---
    subset = cfg.get("subset", {})
    if subset.get("enabled", False):
        ids = set(subset.get("sample_ids", []))
        if not ids:
            sys.exit("ERROR: subset.enabled=true but subset.sample_ids is empty")
        before = len(cfg["samples"])
        cfg["samples"] = [s for s in cfg["samples"] if s["id"] in ids]
        missing = ids - {s["id"] for s in cfg["samples"]}
        if missing:
            sys.exit(f"ERROR: subset.sample_ids not in manifest: {sorted(missing)}")
        print(f"  Subset: {len(cfg['samples'])}/{before} samples")

    # --- Resolve relative h5 paths (and raw_h5 if present) ---
    cwd = Path.cwd()
    for s in cfg["samples"]:
        for key in ("h5", "raw_h5"):
            if key in s and s[key]:
                p = Path(s[key])
                if not p.is_absolute():
                    s[key] = str((cwd / p).resolve())

    # --- Resolve results_dir to absolute too ---
    # Phase scripts build input/output paths from results_dir. Leaving it
    # relative makes every script silently depend on the current working
    # directory — "file not found" when the file actually exists, just under a
    # different CWD. Anchoring it here removes that footgun: scripts work from
    # anywhere, not only from the repo root.
    if "results_dir" in cfg:
        rp = Path(cfg["results_dir"])
        if not rp.is_absolute():
            cfg["results_dir"] = str((cwd / rp).resolve())

    return cfg


# ============================================================================
# Lognorm layer (recomputed on demand from raw counts)
# ============================================================================

def add_lognorm(adata, layer_name: str = "lognorm", target_sum: float = 1e4) -> None:
    """Add a log-normalized layer to `adata` IN PLACE.

    Idempotent: if the layer already exists, it's overwritten.

    Pipeline:  layer = log1p( normalize_total(.X, target_sum) )

    Why this exists: the production pipeline stores raw counts in .X (scVI
    needs them) and DROPS the lognorm layer after Phase 5 to save disk
    (project doc §3 — don't carry redundant layers at 600K-cell scale).
    Notebooks and downstream analyses recompute it on demand via this function.

    Args:
        adata: AnnData with raw counts in .X
        layer_name: where to store the result (default "lognorm")
        target_sum: cell sum target before log1p (default 1e4, Scanpy default)
    """
    # Local imports so importing _utils doesn't drag scanpy in for callers that
    # don't need it (e.g. 01_validate.py).
    import scanpy as sc

    # Operate on a temp copy so we don't disturb .X. scanpy's layer= kwarg
    # to normalize_total/log1p works in modern scanpy (>=1.9), but going via
    # a copy is bulletproof across versions.
    tmp = adata.copy()
    sc.pp.normalize_total(tmp, target_sum=target_sum)
    sc.pp.log1p(tmp)
    adata.layers[layer_name] = tmp.X
    del tmp


# ============================================================================
# Per-phase output paths
# ============================================================================

# Canonical phase → directory mapping. Mirrors project doc §6 layout.
# Note: Phase 0 (validation) writes to results/{tissue}/validation/ as a single
# flat dir — not split into plots/h5ad/tables — and is handled directly by
# 01_validate.py without going through phase_paths().
_PHASE_H5AD_SUBDIRS = {
    "qc":              "03_qc_filtered",
    "doublets":        "04_doublets_removed",
    "integration_prep":"05_integration_ready",
    "integration":     "06_integrated",
}

_PHASE_PLOT_SUBDIRS = {
    "qc":               "02_qc",
    "doublets":         "03_doublets",
    "integration_prep": "04_integration_prep",
    "integration":      "05_integration",
}


def phase_paths(cfg: dict, phase: str) -> dict:
    """Return a dict of standard per-phase output paths, with parents created.

    Keys returned:
      - results_dir : Path  — the tissue's root (results/brain/, results/dev/, etc.)
      - h5ad        : Path  — h5ad checkpoint dir for this phase
      - plots       : Path  — plot output dir for this phase
      - tables      : Path  — shared tables dir (across phases)

    Phases recognized: qc, doublets, integration_prep, integration
    (Phase 0 validation uses its own path convention — see 01_validate.py.)
    """
    if phase not in _PHASE_PLOT_SUBDIRS:
        raise ValueError(f"unknown phase: {phase!r}. "
                         f"Known: {sorted(_PHASE_PLOT_SUBDIRS)}")

    results_dir = Path(cfg["results_dir"])
    plots = results_dir / "plots" / _PHASE_PLOT_SUBDIRS[phase]
    tables = results_dir / "tables"
    h5ad = results_dir / "h5ad" / _PHASE_H5AD_SUBDIRS[phase]
    h5ad.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)

    return {"results_dir": results_dir, "h5ad": h5ad,
            "plots": plots, "tables": tables}


# ============================================================================
# Compute selection (GPU when present, CPU fallback)
# ============================================================================

def select_accelerator(force_cpu: bool = False) -> tuple[str, str]:
    """Return (accelerator, precision) for scvi-tools / pytorch-lightning kwargs.

    Returns ("cpu", "32-true") on Mac / no-GPU systems.
    Returns ("gpu", "bf16-mixed") when CUDA is available (workstation Ada).

    Per project doc §3: BF16 on Ada gives ~1.3-1.5× speedup vs FP32 with no
    accuracy loss; the 4500 Ada has strong BF16 throughput.
    """
    if force_cpu:
        return "cpu", "32-true"
    try:
        import torch
        if torch.cuda.is_available():
            return "gpu", "bf16-mixed"
    except Exception:
        pass
    return "cpu", "32-true"


# ============================================================================
# Declarative contrast spec (Phase 8)
# ============================================================================

# Recognized flags — anything else is a typo we want to catch early.
_VALID_CONTRAST_FLAGS = {
    "primary", "secondary", "confounded_with_pool",
    "underpowered_exploratory", "derived",
}


def load_contrasts(cfg: dict, kind: str = "de") -> dict:
    """Return the validated `contrasts:` block from a tissue config.

    kind: "de" returns only DE-style contrasts (those with a `test`); "derived"
          returns only post-hoc set-operation contrasts (those with
          `source_contrast`); "all" returns everything. Phase 8a/8b/8c/8e want
          "de"; 8g wants "derived".

    Validates each contrast so malformed entries fail HERE (config load), not
    deep inside a multi-hour DE run. Hard-fails on:
      - missing `contrasts` block
      - unknown flag
      - pairwise contrast missing `levels` (unless omnibus or derived)
      - derived contrast missing `source_contrast`
    """
    contrasts = cfg.get("contrasts")
    if not contrasts:
        import sys
        sys.exit(
            "ERROR: no 'contrasts:' block in config. Phase 8 needs one.\n"
            "  Re-run build_yaml.py (which now emits it) for this tissue."
        )

    import sys
    validated = {}
    for name, spec in contrasts.items():
        flag = spec.get("flag")
        if flag not in _VALID_CONTRAST_FLAGS:
            sys.exit(f"ERROR: contrast '{name}' has unknown flag {flag!r}. "
                     f"Valid: {sorted(_VALID_CONTRAST_FLAGS)}")

        is_derived = "source_contrast" in spec
        is_omnibus = spec.get("test") == "group_omnibus"
        is_interaction = spec.get("test", "").find(":") >= 0
        is_pairwise_age = "pairwise" in spec

        if is_derived:
            # post-hoc set op: must point at a real source contrast
            src_name = spec["source_contrast"]
            if src_name not in contrasts:
                sys.exit(f"ERROR: derived contrast '{name}' references unknown "
                         f"source_contrast '{src_name}'.")
        else:
            # DE-style: needs a design + test
            if "design" not in spec:
                sys.exit(f"ERROR: contrast '{name}' missing 'design'.")
            if "test" not in spec:
                sys.exit(f"ERROR: contrast '{name}' missing 'test'.")
            # pairwise level contrasts need exactly [test_level, reference_level]
            if not (is_omnibus or is_interaction or is_pairwise_age):
                lv = spec.get("levels")
                if not (isinstance(lv, list) and len(lv) == 2):
                    sys.exit(f"ERROR: contrast '{name}' needs levels=[test, reference] "
                             f"(got {lv!r}).")

        # Filter by requested kind
        if kind == "de" and is_derived:
            continue
        if kind == "derived" and not is_derived:
            continue
        validated[name] = spec

    if not validated:
        sys.exit(f"ERROR: no contrasts of kind '{kind}' found in config.")
    return validated
