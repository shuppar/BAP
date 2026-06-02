#!/usr/bin/env Rscript
# Install R packages needed for the snRNA-seq pipeline.
# Called from setup-remote.sh. Idempotent — safe to re-run.
#
# Packages installed:
#   - renv               : R env management (lock file = renv.lock)
#   - BiocManager        : Bioconductor installer
#   - scDblFinder        : doublet detection (primary, called as subprocess from Python)
#   - edgeR              : DE cross-check (secondary)
#   - CellChat           : cell-cell communication (called as subprocess from Python)
#   - jsonlite           : JSON I/O for Python <-> R data exchange
#   - Matrix             : sparse matrices (required by scDblFinder)
#   - SingleCellExperiment : SCE object format (scDblFinder input)
#   - DropletUtils       : 10x-format I/O on the R side
#   - speckle / limma    : propeller composition analysis (Phase 8a)
#
# All packages installed into a project-local library managed by renv,
# so the system R installation isn't touched.

# -- 0. Fail loudly on any error --
options(error = function() {
  traceback(2)
  if (!interactive()) quit(status = 1)
})

# -- 1. Tell renv where to put things --
# Use a project-local library. Avoids polluting system R.
RENV_PATHS_ROOT <- file.path(getwd(), ".renv-cache")
Sys.setenv(RENV_PATHS_ROOT = RENV_PATHS_ROOT)

# -- 2. Bootstrap renv itself --
if (!requireNamespace("renv", quietly = TRUE)) {
  message("[setup] Installing renv from CRAN...")
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# -- 3. Initialize / restore renv project --
# If renv.lock exists, restore from it (reproducible install).
# Otherwise, initialize a new renv project and install packages fresh.
if (file.exists("renv.lock")) {
  message("[setup] renv.lock found — restoring exact package versions...")
  renv::restore(prompt = FALSE)
} else {
  message("[setup] No renv.lock — initializing new renv project...")
  renv::init(bare = TRUE, force = TRUE)

  # -- 4. Install BiocManager (gateway to Bioconductor) --
  message("[setup] Installing BiocManager...")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

  # Pin Bioconductor version. 3.19 is current stable for R 4.4 as of mid-2026.
  # Update this when bumping R.
  BIOC_VERSION <- "3.19"
  BiocManager::install(version = BIOC_VERSION, ask = FALSE, update = FALSE)

  # -- 5. Install CRAN packages --
  cran_packages <- c(
    "jsonlite",       # JSON I/O for Python <-> R bridge
    "Matrix",         # sparse matrix support
    "data.table",     # fast tables
    "optparse",       # CLI parsing in R scripts
    "remotes"         # for installing from GitHub (CellChat)
  )
  message("[setup] Installing CRAN packages: ", paste(cran_packages, collapse = ", "))
  install.packages(cran_packages, repos = "https://cloud.r-project.org")

  # -- 6. Install Bioconductor packages --
  bioc_packages <- c(
    "scDblFinder",
    "edgeR",
    "SingleCellExperiment",
    "DropletUtils",
    "scran",          # used by scDblFinder for normalization
    "scuttle",        # used by scDblFinder for QC
    "speckle",        # propeller — Phase 8a composition (replaces scCODA)
    "limma"           # required by propeller (empirical-Bayes moderation)
  )
  message("[setup] Installing Bioconductor packages: ", paste(bioc_packages, collapse = ", "))
  BiocManager::install(bioc_packages, ask = FALSE, update = FALSE)

  # -- 7. Install CellChat from GitHub --
  # CellChat is not on CRAN/Bioconductor — install from the official repo.
  message("[setup] Installing CellChat from GitHub...")
  remotes::install_github("jinworks/CellChat", upgrade = "never")

  # -- 8. Snapshot for reproducibility --
  # Writes renv.lock with exact versions of every installed package.
  # COMMIT renv.lock TO GIT.
  message("[setup] Writing renv.lock...")
  renv::snapshot(prompt = FALSE)
}

# -- 9. Sanity check: load every package once to confirm it works --
message("[setup] Verifying packages load correctly...")
required <- c("scDblFinder", "edgeR", "CellChat", "jsonlite",
              "SingleCellExperiment", "DropletUtils", "speckle", "limma")
for (pkg in required) {
  ok <- suppressMessages(requireNamespace(pkg, quietly = TRUE))
  if (!ok) {
    stop(sprintf("[setup] FAILED to load %s after install", pkg))
  }
  message(sprintf("  ✓ %s", pkg))
}

message("\n[setup] R environment ready.")
message("[setup] To use from R: setwd(<project>) and source this script's directory.")
message("[setup] To regenerate lock file after adding packages: renv::snapshot()")
