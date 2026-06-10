#!/usr/bin/env bash
# setup-remote.sh -- bootstrap the snRNA-seq pipeline environment on the workstation.
#
# Idempotent: safe to re-run. Skips steps already done.
#
# Prerequisites on the workstation:
#   - Linux x86_64
#   - sudo access (for installing R if not present)
#   - Network access to PyPI, CRAN, Bioconductor, GitHub, AWS S3 (no conda needed)
#
# What this script does:
#   1.   Installs uv (single static binary in ~/.local/bin) if not present
#   2.   Creates ./.venv via `uv sync`
#   2.5  Patches CellTypist's train.py to remove `multi_class='ovr'` argument
#        (sklearn 1.7+ removed it; CellTypist still passes it -> hard crash on
#        any LogisticRegression-based training. See ADR comment in step 2.5).
#   3.   Installs abc_atlas_access from GitHub (Phase 7 brain reference download)
#   4.   Installs R via apt if not present
#   5.   Installs R packages (scDblFinder, edgeR, CellChat, speckle, msigdbr, ...)
#   6.   Sets up the CellBender sidecar venv (.venv-cellbender)
#        NOTE: CellBender is currently NOT used in the pipeline (skipped 2026-06-05
#        due to an unresolved checkpoint pickle bug; see INSTRUCTIONS.md). Venv
#        kept bootstrapped in case we revisit via Docker image or SoupX.
#   7.   Fetches MSigDB gene sets (refs/msigdb_mouse.tsv) if not present
#   8.   Builds the brain reference -- two sub-steps:
#          8a. refs/abc_brain_ref.h5ad        (validated reference for Phase 7c)
#          8b. Three CellTypist models + two mapping CSVs from
#              scripts/train_celltypist_brain.py:
#                refs/celltypist_brain_adult_class.pkl     (34 ABC classes)
#                refs/celltypist_brain_adult_subclass.pkl  (~334 ABC subclasses)
#                refs/celltypist_brain_adult_region.pkl    (~12 ABC regions)
#                refs/abc_class_to_broad.csv
#                refs/abc_subclass_to_region.csv
#   9.   Prints next-step instructions
#
# Usage:
#   ./setup-remote.sh                    # full setup
#   ./setup-remote.sh --skip-r           # env-only (skip R + Bioc + msigdb + references)
#   ./setup-remote.sh --skip-cellbender  # skip CellBender sidecar venv
#   ./setup-remote.sh --skip-references  # skip steps 7+8 (long-running atlas download)

set -euo pipefail

# -- Flags --
SKIP_R=0
SKIP_CELLBENDER=0
SKIP_REFERENCES=0
for arg in "$@"; do
  case $arg in
    --skip-r)          SKIP_R=1 ;;
    --skip-cellbender) SKIP_CELLBENDER=1 ;;
    --skip-references) SKIP_REFERENCES=1 ;;
    -h|--help)
      grep '^#' "$0" | head -50
      exit 0 ;;
    *) echo "Unknown flag: $arg" >&2; exit 1 ;;
  esac
done

# -- Pretty logging --
log()  { printf "\033[1;34m[setup]\033[0m %s\n" "$*"; }
warn() { printf "\033[1;33m[warn ]\033[0m %s\n" "$*" >&2; }
err()  { printf "\033[1;31m[error]\033[0m %s\n" "$*" >&2; }

# -- Sanity check --
if [[ "$(uname -s)" != "Linux" ]]; then
  err "This script targets Linux. Detected: $(uname -s)"
  exit 1
fi
if [[ "$(uname -m)" != "x86_64" ]]; then
  warn "Architecture $(uname -m) is not x86_64; PyTorch wheel selection may fail."
fi

# ============================================================================
# 1. Install uv
# ============================================================================
log "Checking for uv..."
if ! command -v uv >/dev/null 2>&1; then
  log "Installing uv (single binary, no admin required)..."
  curl -LsSf https://astral.sh/uv/install.sh | sh
  export PATH="$HOME/.local/bin:$PATH"
  log "uv installed. Add ~/.local/bin to your PATH permanently in ~/.bashrc:"
  log "    echo 'export PATH=\"\$HOME/.local/bin:\$PATH\"' >> ~/.bashrc"
else
  log "uv already installed: $(uv --version)"
fi

# ============================================================================
# 2. Create Python venv and install dependencies
# ============================================================================
log "Creating Python venv with uv sync..."
uv sync
log "Python venv created at $(pwd)/.venv"

# ============================================================================
# 2.5. Patch CellTypist for sklearn 1.7+ compatibility
# ============================================================================
# WHY: CellTypist's train.py hardcodes `multi_class='ovr'` in two
# LogisticRegression() calls (lines 126 and 146 in celltypist 1.6.x).
# scikit-learn 1.7 REMOVED the multi_class parameter entirely, so any call
# to celltypist.train() with use_SGD=False crashes with:
#   TypeError: LogisticRegression.__init__() got an unexpected keyword argument 'multi_class'
#
# Dropping the argument is also semantically correct: sklearn 1.7+ uses
# multinomial (softmax) by default, which gives properly calibrated
# probabilities. The previous OvR sigmoid behaviour was what made
# `celltypist_class_conf` come out bimodal at 0/1 on our first round of
# brain training -- so this patch also fixes that downstream pain.
#
# This patch must be re-applied after every `uv sync` because the file lives
# in .venv/ and uv has no concept of "vendor patches". When CellTypist ships
# a fix upstream (or we pin sklearn<1.7), this step can be removed.
log "Patching CellTypist for sklearn 1.7+ (remove multi_class='ovr')..."
CELLTYPIST_TRAIN_PY=$(find ./.venv -path '*/celltypist/train.py' -type f 2>/dev/null | head -1 || true)
if [[ -z "${CELLTYPIST_TRAIN_PY:-}" ]]; then
  warn "  could not locate celltypist/train.py inside .venv -- skipping patch."
  warn "  if CellTypist is meant to be installed, check uv.lock + 'uv sync' output."
elif grep -q "multi_class = 'ovr'" "$CELLTYPIST_TRAIN_PY"; then
  log "  found multi_class='ovr' in $CELLTYPIST_TRAIN_PY -- patching..."
  # Two patterns: one for "multi_class = 'ovr', " mid-arglist (line 126),
  # and one for "multi_class = 'ovr'" standalone (line 146).
  sed -i "s/multi_class = 'ovr', //" "$CELLTYPIST_TRAIN_PY"
  sed -i "s/multi_class = 'ovr'//" "$CELLTYPIST_TRAIN_PY"
  if grep -q "multi_class = 'ovr'" "$CELLTYPIST_TRAIN_PY"; then
    err "  patch incomplete -- multi_class still present. Inspect manually:"
    err "    grep -n multi_class $CELLTYPIST_TRAIN_PY"
    exit 1
  fi
  log "  patched OK."
else
  log "  CellTypist already patched (no multi_class refs)."
fi

# ============================================================================
# 3. Install abc_atlas_access from GitHub (Phase 7 reference build)
# ============================================================================
log "Checking abc_atlas_access (Allen Brain Cell Atlas API)..."
if uv run python -c "import abc_atlas_access" 2>/dev/null; then
  log "abc_atlas_access already installed."
else
  log "Installing abc_atlas_access from GitHub..."
  uv pip install "git+https://github.com/AllenInstitute/abc_atlas_access.git"
fi

# ============================================================================
# 4. Install R
# ============================================================================
if [[ $SKIP_R -eq 0 ]]; then
  log "Checking for R..."
  if ! command -v R >/dev/null 2>&1; then
    log "R not found -- installing via apt (requires sudo)..."
    if ! command -v sudo >/dev/null 2>&1; then
      err "sudo not available and R not installed. Install R manually, then re-run with --skip-r."
      exit 1
    fi
    sudo apt-get update
    sudo apt-get install -y --no-install-recommends \
      r-base r-base-dev \
      libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
      libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
      libtiff5-dev libjpeg-dev libhdf5-dev
    log "R installed: $(R --version | head -1)"
  else
    log "R already installed: $(R --version | head -1)"
  fi

  # ==========================================================================
  # 5. Install R packages (scDblFinder, edgeR, CellChat, speckle, msigdbr, ...)
  # ==========================================================================
  log "Installing R packages..."
  log "This step can take 20-40 minutes on first run (compiling from source)."
  Rscript scripts/install-r-packages.R
  log "R packages installed."
else
  log "Skipping R setup (--skip-r given)"
fi

# ============================================================================
# 6. Set up isolated CellBender venv
# ============================================================================
# NOTE: CellBender is currently NOT in the pipeline (see INSTRUCTIONS.md
# "Isolate fragile dependency stacks"). We still bootstrap the venv here so
# the path stays valid if we revisit via the Docker image or SoupX alternative.
if [[ $SKIP_CELLBENDER -eq 0 ]]; then
  log "Setting up isolated CellBender venv..."
  if [[ ! -d ".venv-cellbender" ]]; then
    uv venv .venv-cellbender --python 3.10
    log "Created .venv-cellbender (Python 3.10)"
  fi
  log "Installing CellBender..."
  uv pip install --python .venv-cellbender/bin/python \
    cellbender \
    "torch>=2.0,<2.5"
  log "CellBender venv ready (currently unused by pipeline -- see INSTRUCTIONS.md)."
else
  log "Skipping CellBender venv (--skip-cellbender given)"
fi

# ============================================================================
# 7. Fetch MSigDB gene sets (one-time, ~1 min) -- Phase 8c GSEA + leading-edge
# ============================================================================
if [[ $SKIP_REFERENCES -eq 0 && $SKIP_R -eq 0 ]]; then
  if [[ -f "refs/msigdb_mouse.tsv" ]]; then
    log "refs/msigdb_mouse.tsv already present -- skipping fetch."
  else
    log "Fetching MSigDB mouse gene sets -> refs/msigdb_mouse.tsv ..."
    mkdir -p refs
    Rscript scripts/fetch_genesets.R --out refs/msigdb_mouse.tsv
    log "MSigDB gene sets ready."
  fi
fi

# ============================================================================
# 8. Build brain reference (ABC atlas h5ad + three CellTypist models)
# ============================================================================
# 8a. Validated reference h5ad with obs columns: class, subclass,
#     anatomical_division_label. Produced by scripts/prepare_brain_reference.py
#     (downloads ABC WMB-10Xv3 + joins per-cell metadata + writes
#     refs/abc_brain_ref.h5ad). LONG step: ~tens of GB downloaded from AWS S3.
# 8b. Three CellTypist models + two mapping CSVs from
#     scripts/train_celltypist_brain.py:
#       - class    (34 ABC classes)         ~20-40 min
#       - subclass (~334 ABC subclasses)    ~45-90 min
#       - region   (~12 ABC regions)        ~10-20 min
#       - abc_class_to_broad.csv            (instant)
#       - abc_subclass_to_region.csv        (instant)
#     Idempotent: skips any output that already exists. Use --force on the
#     script to rebuild from scratch.
if [[ $SKIP_REFERENCES -eq 0 ]]; then
  log "=== Step 8a: brain reference h5ad ==="
  if [[ -f "refs/abc_brain_ref.h5ad" ]]; then
    log "refs/abc_brain_ref.h5ad already present -- skipping ABC download."
  else
    log "Building refs/abc_brain_ref.h5ad from ABC atlas (LONG step, run in tmux)..."
    log "  This downloads ~tens of GB from AWS S3."
    uv run python scripts/prepare_brain_reference.py
  fi

  log "=== Step 8b: CellTypist brain models (class/subclass/region) ==="
  brain_outputs=(
    refs/celltypist_brain_adult_class.pkl
    refs/celltypist_brain_adult_subclass.pkl
    refs/celltypist_brain_adult_region.pkl
    refs/abc_class_to_broad.csv
    refs/abc_subclass_to_region.csv
  )
  missing=0
  for f in "${brain_outputs[@]}"; do
    if [[ ! -f "$f" ]]; then
      missing=1
      log "  missing: $f"
    fi
  done
  if [[ $missing -eq 1 ]]; then
    if [[ ! -f "refs/abc_brain_ref.h5ad" ]]; then
      warn "refs/abc_brain_ref.h5ad is missing -- step 8b skipped."
      warn "  Re-run setup-remote.sh after the reference download completes."
    else
      log "Training CellTypist brain models (~1.5-2.5 hours total on WS)..."
      log "  Run in tmux if you intend to disconnect:"
      log "    tmux new -s train 'uv run python -u scripts/train_celltypist_brain.py --config config/brain.yaml 2>&1 | tee logs/train_celltypist_brain.log'"
      log "  Running inline now (will block this script until done)..."
      uv run python -u scripts/train_celltypist_brain.py --config config/brain.yaml
    fi
  else
    log "All CellTypist brain outputs present -- skipping training."
    log "  To rebuild: uv run python scripts/train_celltypist_brain.py --config config/brain.yaml --force"
  fi
else
  log "Skipping reference build (--skip-references given)"
fi

# ============================================================================
# 9. Summary & next steps
# ============================================================================
echo
log "=========================================="
log "Environment setup complete."
log "=========================================="
log ""
log "What got created:"
log "  ./.venv/                                       -- main Python env (scanpy, scvi-tools, celltypist, ...)"
log "  ./.venv-cellbender/                            -- isolated CellBender env (currently unused)"
log "  ./.renv-cache/                                 -- R package library (project-local)"
log "  ./renv.lock                                    -- R lock file (COMMIT TO GIT)"
log "  ./uv.lock                                      -- Python lock file (COMMIT TO GIT)"
log "  ./refs/msigdb_mouse.tsv                        -- Phase 8c pathway gene sets"
log "  ./refs/abc_brain_ref.h5ad                      -- Phase 7c scANVI reference (brain)"
log "  ./refs/celltypist_brain_adult_class.pkl        -- Phase 7 class model (34 ABC classes)"
log "  ./refs/celltypist_brain_adult_subclass.pkl     -- Phase 7 subclass model (~334)"
log "  ./refs/celltypist_brain_adult_region.pkl       -- Phase 7 region model (~12)"
log "  ./refs/abc_class_to_broad.csv                  -- deterministic class -> broad map"
log "  ./refs/abc_subclass_to_region.csv              -- subclass -> region majority (informational)"
log ""
log "Tech debt / known issues:"
log "  - The sklearn 1.7+ CellTypist patch (step 2.5) lives in .venv/ and is"
log "    re-applied on every setup-remote.sh run. If you ever 'uv sync' WITHOUT"
log "    re-running setup-remote.sh, the patch will be silently reverted and"
log "    LogisticRegression-based CellTypist training will crash again."
log "    Long-term fix: pin sklearn<1.7 in pyproject.toml or wait for upstream."
log ""
log "Next steps:"
log "  1. Add ~/.local/bin to PATH in ~/.bashrc (if uv was installed fresh)"
log "  2. Re-generate brain.yaml / placenta.yaml so they pick up the new reference paths:"
log "       uv run python scripts/build_yaml.py"
log "  3. Verify GPU is accessible:"
log "       uv run python -c 'import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))'"
log "  4. Verify main env imports OK:"
log "       uv run python -c 'import scanpy, scvi, anndata, celltypist, pydeseq2, decoupler, omnipath, liana; print(\"deps OK\")'"
log "  5. Begin pipeline at Phase 0:"
log "       uv run python scripts/01_validate.py --config config/brain.yaml"
log ""
