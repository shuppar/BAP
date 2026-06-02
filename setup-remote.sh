#!/usr/bin/env bash
# setup-remote.sh — bootstrap the snRNA-seq pipeline environment on the workstation.
#
# Idempotent: safe to re-run. Skips steps already done.
#
# Prerequisites on the workstation:
#   - Linux x86_64
#   - sudo access (for installing R if not present)
#   - Network access to PyPI, CRAN, Bioconductor (verified — conda channels NOT required)
#
# What this script does:
#   1. Installs uv (single static binary in ~/.local/bin) if not present
#   2. Creates a project-local Python venv at ./.venv via `uv sync`
#   3. Installs R via apt if not present
#   4. Installs the R packages needed for scDblFinder / edgeR / CellChat / speckle
#   5. Sets up a separate venv for CellBender (different PyTorch pin)
#   6. Prints next-step instructions
#
# Usage:
#   ./setup-remote.sh                  # full setup
#   ./setup-remote.sh --skip-r         # skip R steps (if R already configured)
#   ./setup-remote.sh --skip-cellbender # skip CellBender env (set up later)

set -euo pipefail

# -- Flags --
SKIP_R=0
SKIP_CELLBENDER=0
for arg in "$@"; do
  case $arg in
    --skip-r) SKIP_R=1 ;;
    --skip-cellbender) SKIP_CELLBENDER=1 ;;
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
  # Add to PATH for this session; user should source their shell rc afterward
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
# uv reads pyproject.toml + .python-version, downloads Python 3.11 if needed,
# creates .venv/, installs all pinned deps, writes uv.lock.
uv sync

log "Python venv created at $(pwd)/.venv"
log "  To activate manually:  source .venv/bin/activate"
log "  Or use uv run:         uv run python run.py ..."

# ============================================================================
# 3. Install R
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
  # 4. Install R packages (scDblFinder, edgeR, CellChat, etc.)
  # ==========================================================================
  log "Installing R packages (scDblFinder, edgeR, CellChat, ...)"
  log "This step can take 20-40 minutes on first run (compiling from source)."
  Rscript scripts/install-r-packages.R
  log "R packages installed and snapshotted to renv.lock"
else
  log "Skipping R setup (--skip-r given)"
fi

# ============================================================================
# 5. Set up isolated CellBender venv
# ============================================================================
if [[ $SKIP_CELLBENDER -eq 0 ]]; then
  log "Setting up isolated CellBender venv..."

  # CellBender pins PyTorch hard and tends to conflict with scvi-tools.
  # We give it its own venv. Python pipeline calls it via subprocess.
  if [[ ! -d ".venv-cellbender" ]]; then
    uv venv .venv-cellbender --python 3.10
    log "Created .venv-cellbender (Python 3.10)"
  fi

  log "Installing CellBender..."
  # --python flag points uv pip at the alternate venv
  uv pip install --python .venv-cellbender/bin/python \
    cellbender \
    "torch>=2.0,<2.5"   # CellBender's tested PyTorch range as of mid-2026

  log "CellBender venv ready."
  log "  To run:  .venv-cellbender/bin/cellbender remove-background ..."
else
  log "Skipping CellBender venv (--skip-cellbender given)"
fi

# ============================================================================
# 5b. (removed) scCODA sidecar — composition analysis now uses propeller
#     (speckle + limma) via R subprocess. speckle/limma install with the other
#     Bioconductor packages in scripts/install-r-packages.R (step 4 above), so
#     there is no separate venv to build here. scDblFinder/edgeR/CellChat/speckle
#     all share the project-local renv library.
# ============================================================================

# ============================================================================
# 6. Summary & next steps
# ============================================================================
echo
log "=========================================="
log "Environment setup complete."
log "=========================================="
log ""
log "What got created:"
log "  ./.venv/                   — main Python env (scanpy, scvi-tools, ...)"
log "  ./.venv-cellbender/        — isolated CellBender env"
log "  ./.renv-cache/             — R package library (project-local)"
log "  ./renv.lock                — R lock file (COMMIT TO GIT)"
log "  ./uv.lock                  — Python lock file (COMMIT TO GIT)"
log ""
log "Next steps:"
log "  1. Add ~/.local/bin to PATH in ~/.bashrc (if uv was installed fresh)"
log "  2. Verify GPU is accessible:"
log "     uv run python -c 'import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))'"
log "  3. Run Phase 0 validation:"
log "     uv run python run.py --config config/brain.yaml --step validate"
log ""
