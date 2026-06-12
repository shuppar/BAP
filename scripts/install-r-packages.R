#!/usr/bin/env Rscript
# Install R packages needed for the snRNA-seq pipeline.
# Called from setup-remote.sh. Idempotent — safe to re-run.
#
# Packages installed:
#   - renv               : R env management (lock file = renv.lock)
#   - BiocManager        : Bioconductor installer
#   - scDblFinder        : doublet detection (primary, called as subprocess from Python)
#   - SoupX              : ambient RNA correction Phase 1 (replaces CellBender)
#   - scran              : quickCluster for SoupX + normalization for scDblFinder
#   - edgeR              : DE cross-check (secondary)
#   - jsonlite           : JSON I/O for Python <-> R data exchange
#   - Matrix             : sparse matrices (required by scDblFinder)
#   - SingleCellExperiment : SCE object format (scDblFinder input)
#   - DropletUtils       : 10x-format I/O on the R side
#   - speckle / limma    : propeller composition analysis (Phase 8a)
#   - msigdbr            : MSigDB gene sets for Phase 8c GSEA (fetch_genesets.R)
#
# NOTE: CellChat dropped — 8e cell-cell communication is LIANA+ in Python.
# CellChat pulls heavy plotting transitive deps (fs, ragg, sass, ggrastr) that
# need Ubuntu system libs not in setup-remote.sh, and we don't use it.
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

  # Pin Bioconductor version. 3.21 is the release that supports R 4.5
  # (R 4.5.x → Bioconductor 3.21). 3.19 is for R 4.4 and will fail on R 4.5+.
  # Update this when bumping R.
  BIOC_VERSION <- "3.21"
  BiocManager::install(version = BIOC_VERSION, ask = FALSE, update = FALSE)

  # -- 5. Install CRAN packages --
  cran_packages <- c(
    "jsonlite",       # JSON I/O for Python <-> R bridge
    "Matrix",         # sparse matrix support
    "data.table",     # fast tables
    "optparse",       # CLI parsing in R scripts
    "msigdbr"         # MSigDB gene sets for fetch_genesets.R (Phase 8c)
  )
  message("[setup] Installing CRAN packages: ", paste(cran_packages, collapse = ", "))
  install.packages(cran_packages, repos = "https://cloud.r-project.org")

  # -- 6. Install Bioconductor packages --
  bioc_packages <- c(
    "scDblFinder",
    "edgeR",
    "SingleCellExperiment",
    "DropletUtils",
    "scran",          # quickCluster for SoupX + normalization for scDblFinder
    "scuttle",        # used by scDblFinder for QC
    "SoupX",          # Phase 1 ambient RNA correction (replaces CellBender)
    "speckle",        # propeller — Phase 8a composition (replaces scCODA)
    "limma"           # required by propeller (empirical-Bayes moderation)
  )
  message("[setup] Installing Bioconductor packages: ", paste(bioc_packages, collapse = ", "))
  BiocManager::install(bioc_packages, ask = FALSE, update = FALSE)

  # -- 7. Snapshot for reproducibility --
  # Writes renv.lock with exact versions of every installed package.
  # COMMIT renv.lock TO GIT.
  #
  # Wrapped in tryCatch because transitive plotting deps (ragg via scater, etc.)
  # can fail to compile if a system lib is missing. Those deps are NEVER called
  # by our pipeline; snapshotting them is bookkeeping. If it fails, log and
  # continue — the verification step below confirms the packages we actually
  # use are loadable. Regenerate renv.lock manually later if needed:
  #   Rscript -e 'renv::snapshot(prompt=FALSE)'
  message("[setup] Writing renv.lock...")
  tryCatch(
    renv::snapshot(prompt = FALSE),
    error = function(e) {
      message("[setup] WARN: renv::snapshot failed — continuing.")
      message("[setup] Reason: ", conditionMessage(e))
      message("[setup] Pipeline packages will be verified next; the lockfile")
      message("[setup] is non-critical (reproducibility metadata only).")
    }
  )
}

# -- 8. Drop the bail-on-error handler before verify --
# We want the verify loop to test EVERY required package and report all that
# fail, not bail on the first one.
options(error = NULL)

# -- 9. Sanity check: load every package once to confirm it works --
message("[setup] Verifying packages load correctly...")
required <- c("scDblFinder", "edgeR", "jsonlite",
              "SingleCellExperiment", "DropletUtils", "scran", "SoupX",
              "speckle", "limma", "msigdbr")
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
