#!/usr/bin/env bash
# setup-remote.sh — bootstrap the snRNA-seq pipeline environment on the workstation.
#
# Idempotent: safe to re-run. Skips steps already done.
#
# Prerequisites on the workstation:
#   - Linux x86_64
#   - sudo access (for installing R if not present)
#   - Network access to PyPI, CRAN, Bioconductor, GitHub, AWS S3 (no conda needed)
#
# What this script does:
#   1. Installs uv (single static binary in ~/.local/bin) if not present
#   2. Creates ./.venv via `uv sync`
#   3. Installs abc_atlas_access from GitHub (Phase 7 brain reference)
#   4. Installs R via apt if not present
#   5. Installs R packages (scDblFinder, edgeR, CellChat, speckle, msigdbr, ...)
#   6. Sets up the CellBender sidecar venv (.venv-cellbender)
#   7. Fetches MSigDB gene sets (refs/msigdb_mouse.tsv) if not present
#   8. Builds the brain reference (ABC atlas + CellTypist .pkl) if not present
#   9. Prints next-step instructions
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
      grep '^#' "$0" | head -30
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
    log "R not found — installing via apt (requires sudo)..."
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
  log "CellBender venv ready."
else
  log "Skipping CellBender venv (--skip-cellbender given)"
fi

# ============================================================================
# 7. Fetch MSigDB gene sets (one-time, ~1 min) — Phase 8c GSEA + leading-edge
# ============================================================================
if [[ $SKIP_REFERENCES -eq 0 && $SKIP_R -eq 0 ]]; then
  if [[ -f "refs/msigdb_mouse.tsv" ]]; then
    log "refs/msigdb_mouse.tsv already present — skipping fetch."
  else
    log "Fetching MSigDB mouse gene sets → refs/msigdb_mouse.tsv ..."
    mkdir -p refs
    Rscript scripts/fetch_genesets.R --out refs/msigdb_mouse.tsv
    log "MSigDB gene sets ready."
  fi
fi

# ============================================================================
# 8. Build brain reference (ABC atlas + CellTypist .pkl) — Phase 7 + 7c
# ============================================================================
# This is the long step: ~tens of GB downloaded from AWS S3, then ~30 min
# CellTypist training. Idempotent — skipped if outputs exist. Run in tmux.
if [[ $SKIP_REFERENCES -eq 0 ]]; then
  if [[ -f "refs/abc_brain_ref.h5ad" && -f "refs/celltypist_brain_adult.pkl" ]]; then
    log "Brain reference outputs already exist — skipping build."
  else
    log "Building brain reference from ABC atlas (LONG step, run in tmux)..."
    log "  This downloads ~tens of GB and trains a CellTypist model."
    log "  Outputs: refs/abc_brain_ref.h5ad + refs/celltypist_brain_adult.pkl"
    uv run python scripts/prepare_brain_reference.py
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
log "  ./.venv/                            — main Python env (scanpy, scvi-tools, ...)"
log "  ./.venv-cellbender/                 — isolated CellBender env"
log "  ./.renv-cache/                      — R package library (project-local)"
log "  ./renv.lock                         — R lock file (COMMIT TO GIT)"
log "  ./uv.lock                           — Python lock file (COMMIT TO GIT)"
log "  ./refs/msigdb_mouse.tsv             — Phase 8c pathway gene sets"
log "  ./refs/abc_brain_ref.h5ad           — Phase 7c scANVI reference (brain)"
log "  ./refs/celltypist_brain_adult.pkl   — Phase 7 CellTypist model (4W + 3mo brain)"
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
