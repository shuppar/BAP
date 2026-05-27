# Environment Setup

## Why this layout (and not conda)

The remote workstation has **conda channels blocked at the firewall**, but PyPI, CRAN,
and Bioconductor are all reachable. So:

- **Python:** managed by [uv](https://github.com/astral-sh/uv) — fast, single-binary,
  PyPI-only, no admin needed
- **R:** managed by [renv](https://rstudio.github.io/renv/) — packages from CRAN +
  Bioconductor, project-local library
- **Python ↔ R bridge:** Python calls R as subprocesses (see `scripts/run-scdblfinder.R`
  as an example). No rpy2 to keep things robust.

Reproducibility: both `uv.lock` and `renv.lock` are committed to git.

## One-time bootstrap on the workstation

```bash
git clone <repo> snrna-project
cd snrna-project
./setup-remote.sh
```

This installs uv (if missing), creates `.venv/` with all Python deps, installs R
(if missing), installs all R packages, and creates an isolated `.venv-cellbender/`
for CellBender. Takes 20–40 minutes on first run.

## Day-to-day usage

```bash
# From the Mac, push code changes to the remote
make sync

# On the remote (via SSH or VSCode Remote-SSH)
make validate         # Phase 0 — mandatory first run
make qc               # Phase 2
make integration      # Phase 5 (GPU-heavy)
# ...etc — see `make help` for full list

# Verify GPU is accessible
make gpu-check

# Verify R packages load
make r-check
```

## Reproducing the exact environment elsewhere

Anyone with the repo + a Linux x86_64 machine + network access to PyPI/CRAN/Bioc
can run `./setup-remote.sh` and get the same environment (modulo system-library
differences). The lock files (`uv.lock`, `renv.lock`) pin every transitive
dependency to exact versions.

## Updating dependencies

```bash
# Add a Python package
uv add <package>          # updates pyproject.toml + uv.lock

# Add an R package (interactively in R)
R
> renv::install("<package>")
> renv::snapshot()         # updates renv.lock

# Then commit the lock files
git add uv.lock renv.lock
git commit -m "deps: add <package>"
```

## Why CellBender has its own venv

CellBender pins PyTorch tightly and often conflicts with scvi-tools' PyTorch
requirements. We give it its own `.venv-cellbender/` and call it as a subprocess
from `src/ambient.py`. Standard pattern in single-cell pipelines.

## File map

| File | Purpose |
|---|---|
| `pyproject.toml` | Python deps spec (uv reads this) |
| `uv.lock` | Python deps locked to exact versions (commit to git) |
| `.python-version` | Pins Python to 3.11 |
| `renv.lock` | R deps locked to exact versions (commit to git) |
| `scripts/install-r-packages.R` | R-side installer, run by `setup-remote.sh` |
| `scripts/run-scdblfinder.R` | Example R subprocess driver |
| `setup-remote.sh` | One-shot bootstrap for the workstation |
| `Makefile` | Convenience commands (`make sync`, `make validate`, etc.) |
