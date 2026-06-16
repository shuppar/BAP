#!/usr/bin/env Rscript
# Install R packages needed for the snRNA-seq pipeline.
# Called from setup-remote.sh. Idempotent — safe to re-run.
#
# Packages installed:
#   - renv               : R env management (lock file = renv.lock)
#   - BiocManager        : Bioconductor installer
#   - scDblFinder        : doublet detection (primary, called as subprocess from Python)
#   - SoupX              : ambient RNA correction Phase 1 (replaces CellBender)
#   - scran / scuttle    : quickCluster for SoupX + normalization/QC for scDblFinder
#   - edgeR              : DE cross-check (secondary)
#   - jsonlite           : JSON I/O for Python <-> R data exchange
#   - Matrix             : sparse matrices (required by scDblFinder)
#   - data.table         : fast tables
#   - SingleCellExperiment : SCE object format (scDblFinder input)
#   - DropletUtils       : 10x-format I/O on the R side (SoupX)
#   - optparse           : CLI parsing in run_*.R subprocess scripts
#   - speckle / limma    : propeller composition analysis (Phase 8a)
#   - msigdbr            : MSigDB gene sets for Phase 8c GSEA (fetch_genesets.R)
#
# SELF-HEALING (added 2026-06-14): the declared lists below are the single source
# of truth. On EVERY run — fresh init OR restore-from-lock — the script installs
# any declared package that isn't actually present, then re-snapshots. This fixes
# the failure mode where a package added to this list never got installed because
# renv.lock predated it and `renv::restore()` only installs what the lock names.
# So: add a package here, re-run setup-remote.sh, done.
#
# NOTE: CellChat dropped — 8e cell-cell communication is LIANA+ in Python.
# CellChat pulls heavy plotting transitive deps (fs, ragg, sass, ggrastr) that
# need Ubuntu system libs not in setup-remote.sh, and we don't use it.
#
# All packages install into a project-local library managed by renv, so the
# system R installation isn't touched.

# -- 0. Fail loudly on any error (dropped before the verify loop) --
options(error = function() {
  traceback(2)
  if (!interactive()) quit(status = 1)
})

# -- 1. Declared package lists — SINGLE SOURCE OF TRUTH --
# To add a package: put it in the right list and re-run. The top-up step (§6)
# installs it even if renv.lock already exists.
CRAN_PACKAGES <- c(
  "jsonlite",       # JSON I/O for Python <-> R bridge
  "Matrix",         # sparse matrix support
  "data.table",     # fast tables
  "optparse",       # CLI parsing in run_*.R (scDblFinder/SoupX/propeller)
  "msigdbr"         # MSigDB gene sets for fetch_genesets.R (Phase 8c)
)
BIOC_PACKAGES <- c(
  "scDblFinder",          # Phase 3 doublets
  "edgeR",                # DE cross-check
  "SingleCellExperiment", # SCE container
  "DropletUtils",         # 10x I/O (SoupX)
  "scran",                # quickCluster (SoupX) + normalization (scDblFinder)
  "scuttle",              # QC helpers used by scDblFinder
  "SoupX",                # Phase 1 ambient RNA correction
  "speckle",              # propeller — Phase 8a composition (replaces scCODA)
  "limma"                 # empirical-Bayes moderation (propeller dependency)
)
# Everything the pipeline actually loads — used by the top-up and verify steps.
REQUIRED <- c(CRAN_PACKAGES, BIOC_PACKAGES)

# Bioconductor release pinned to the R version. 3.21 supports R 4.5.x.
# Update when bumping R (R 4.4 -> 3.19, R 4.5 -> 3.21).
BIOC_VERSION <- "3.21"

# -- 2. renv library location + Suggests-skip + bundled libuv --
# Project-local library, so system R is untouched.
RENV_PATHS_ROOT <- file.path(getwd(), ".renv-cache")
Sys.setenv(RENV_PATHS_ROOT = RENV_PATHS_ROOT)
# Belt-and-suspenders for the fs -> libuv-dev system-lib chain (project doc §29).
Sys.setenv(USE_BUNDLED_LIBUV = "1")

# -- 3. Bootstrap renv itself --
if (!requireNamespace("renv", quietly = TRUE)) {
  message("[setup] Installing renv from CRAN...")
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# -- 4. Initialize / restore the renv project --
if (file.exists("renv.lock")) {
  message("[setup] renv.lock found — restoring exact package versions...")
  renv::restore(prompt = FALSE)
} else {
  message("[setup] No renv.lock — initializing new renv project...")
  renv::init(bare = TRUE, force = TRUE)
}

# Skip Suggests for ALL installs in this project (Depends/Imports/LinkingTo only).
# Prevents speckle/SoupX/etc. pulling Seurat -> shiny -> bslib -> fs -> libuv-dev
# (project doc §29). persist=TRUE writes it to renv/settings.dcf.
renv::settings$package.dependency.fields(
  c("Depends", "Imports", "LinkingTo"), persist = TRUE
)

# -- 5. Ensure BiocManager is present + the Bioc version is pinned --
# Needed by the top-up below on BOTH branches (restore path may lack it).
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("[setup] Installing BiocManager...")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(version = BIOC_VERSION, ask = FALSE, update = FALSE)

# -- 6. TOP-UP: install any declared package that isn't actually present --
# On a fresh init this installs everything; on a restore it fills only the gaps
# (e.g. a package added to the lists above after renv.lock was written).
missing <- Filter(function(p) !requireNamespace(p, quietly = TRUE), REQUIRED)
if (length(missing) > 0) {
  miss_cran <- intersect(missing, CRAN_PACKAGES)
  miss_bioc <- intersect(missing, BIOC_PACKAGES)
  if (length(miss_cran)) {
    message("[setup] Installing missing CRAN: ", paste(miss_cran, collapse = ", "))
    install.packages(miss_cran, repos = "https://cloud.r-project.org")
  }
  if (length(miss_bioc)) {
    message("[setup] Installing missing Bioconductor: ", paste(miss_bioc, collapse = ", "))
    BiocManager::install(miss_bioc, ask = FALSE, update = FALSE)
  }

  # -- 7. Snapshot so renv.lock reflects the now-complete library --
  # COMMIT renv.lock TO GIT. Wrapped in tryCatch because transitive plotting
  # deps can fail to compile if a system lib is missing; those are never called
  # by our pipeline, so a snapshot failure is non-critical (the verify loop
  # below confirms the packages we actually use are loadable).
  message("[setup] Packages changed — writing renv.lock...")
  tryCatch(
    renv::snapshot(prompt = FALSE),
    error = function(e) {
      message("[setup] WARN: renv::snapshot failed — continuing.")
      message("[setup] Reason: ", conditionMessage(e))
      message("[setup] Regenerate later with: Rscript -e 'renv::snapshot(prompt=FALSE)'")
    }
  )
} else {
  message("[setup] All declared packages already present — nothing to install.")
}

# -- 8. Drop the bail-on-error handler before verify --
# Test EVERY required package and report all failures, not just the first.
options(error = NULL)

# -- 9. Sanity check: load every required package once --
message("[setup] Verifying packages load correctly...")
failed <- character(0)
for (pkg in REQUIRED) {
  ok <- suppressMessages(requireNamespace(pkg, quietly = TRUE))
  if (ok) {
    message(sprintf("  OK %s", pkg))
  } else {
    failed <- c(failed, pkg)
    message(sprintf("  FAIL %s", pkg))
  }
}
if (length(failed) > 0) {
  stop("[setup] FAILED to load after install: ", paste(failed, collapse = ", "))
}

message("\n[setup] R environment ready.")
message("[setup] Regenerate lock file after adding packages: renv::snapshot()")
