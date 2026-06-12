#!/usr/bin/env Rscript
# run_soupx.R -- per-sample SoupX ambient RNA correction.
#
# Invoked by scripts/02_soupx.py once per sample. The contract is:
#   IN:  cellranger filtered counts (h5 or MTX dir)
#        cellranger raw counts      (h5 or MTX dir)
#   OUT: <output_dir>/matrix.mtx          -- corrected counts (genes x cells)
#        <output_dir>/barcodes.tsv        -- cell barcodes (one per col)
#        <output_dir>/features.tsv        -- gene metadata (ensembl, symbol, type)
#        <output_dir>/soupx_summary.json  -- rho, pct_removed, timings
#
# rho estimation: scran::quickCluster -> SoupX::setClusters -> autoEstCont.
# Auto, data-driven. Manual rho override available via --rho.
#
# Dependencies (Bioconductor + CRAN): SoupX DropletUtils scran Matrix
#                                     optparse jsonlite

suppressPackageStartupMessages({
  library(optparse)
  library(SoupX)
  library(DropletUtils)
  library(Matrix)
  library(scran)
  library(jsonlite)
})

option_list <- list(
  make_option("--filtered",   type = "character", help = "Cellranger filtered counts: .h5 OR MTX directory"),
  make_option("--raw",        type = "character", help = "Cellranger raw counts: .h5 OR MTX directory"),
  make_option("--output_dir", type = "character", help = "Where to write outputs"),
  make_option("--sample_id",  type = "character", help = "Sample ID (for logging only)"),
  make_option("--rho",        type = "double",    default = NA, help = "Manual contamination fraction (default: autoEst via scran clusters)"),
  make_option("--min_cluster_size", type = "integer", default = 20, help = "scran::quickCluster min.size [default 20]")
)
args <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(args$filtered), !is.null(args$raw),
          !is.null(args$output_dir), !is.null(args$sample_id))

t0 <- Sys.time()
cat(sprintf("[soupx] %s: starting (R %s)\n", args$sample_id, getRversion()))

# ------------------------------------------------------------------
# Load matrices. read10xCounts handles both .h5 and MTX directories.
# ------------------------------------------------------------------
cat(sprintf("[soupx] %s: loading filtered from %s\n", args$sample_id, args$filtered))
filtered <- read10xCounts(args$filtered, col.names = TRUE, type = "auto")
toc <- counts(filtered)  # table of counts (filtered cells), genes x cells

cat(sprintf("[soupx] %s: loading raw from %s\n", args$sample_id, args$raw))
raw <- read10xCounts(args$raw, col.names = TRUE, type = "auto")
tod <- counts(raw)       # table of droplets (raw), genes x droplets

cat(sprintf("[soupx] %s: filtered=%d cells, raw=%d droplets, %d genes\n",
            args$sample_id, ncol(toc), ncol(tod), nrow(toc)))

# Align features via rowData$ID (Ensembl). read10xCounts doesn't always set
# matrix rownames identically for the H5 vs MTX paths, so don't rely on
# `identical(rownames(...))` -- use the curated rowData$ID column.
if (!"ID" %in% colnames(rowData(filtered)) || !"ID" %in% colnames(rowData(raw))) {
  stop(sprintf("[soupx] %s: missing rowData$ID column on filtered or raw",
               args$sample_id))
}
toc_ids <- rowData(filtered)$ID
tod_ids <- rowData(raw)$ID

if (!identical(toc_ids, tod_ids)) {
  common <- intersect(toc_ids, tod_ids)
  if (length(common) < 0.9 * length(toc_ids)) {
    stop(sprintf(
      "[soupx] %s: only %d common features (<90%% of filtered's %d). Cannot align.",
      args$sample_id, length(common), length(toc_ids)
    ))
  }
  filtered <- filtered[match(common, toc_ids), ]
  raw      <- raw[match(common, tod_ids), ]
  toc      <- counts(filtered)
  tod      <- counts(raw)
  cat(sprintf(
    "[soupx] %s: aligned to %d common features (filtered=%d, raw=%d)\n",
    args$sample_id, length(common), length(toc_ids), length(tod_ids)
  ))
}

# Force consistent rownames using Ensembl IDs (SoupX checks rownames match).
rownames(toc) <- rowData(filtered)$ID
rownames(tod) <- rowData(raw)$ID

# Final sanity check
if (!identical(rownames(toc), rownames(tod))) {
  stop(sprintf("[soupx] %s: rownames still differ after alignment (bug)",
               args$sample_id))
}

# Materialize to in-memory dgCMatrix. read10xCounts returns a TENxMatrix
# (HDF5-backed DelayedArray) for H5 inputs, which SoupX::quickMarkers can't
# coerce to "dMatrix" downstream. Force in-memory CsparseMatrix/dgCMatrix.
# Memory cost: filtered ~50-150 MB; raw ~1-3 GB at 707K droplets x 19K
# genes. Comfortable under the workstation's 258 GB budget.
cat(sprintf("[soupx] %s: materializing matrices to dgCMatrix (in-memory sparse)\n",
            args$sample_id))
toc <- as(toc, "CsparseMatrix")
if (!inherits(toc, "dgCMatrix")) toc <- as(toc, "dgCMatrix")
tod <- as(tod, "CsparseMatrix")
if (!inherits(tod, "dgCMatrix")) tod <- as(tod, "dgCMatrix")

# ------------------------------------------------------------------
# Build SoupChannel + cluster + estimate rho + adjust counts.
# ------------------------------------------------------------------
sc <- SoupChannel(tod, toc)

if (is.na(args$rho)) {
  cat(sprintf("[soupx] %s: scran::quickCluster (min.size=%d)\n",
              args$sample_id, args$min_cluster_size))
  clust <- as.character(quickCluster(toc, min.size = args$min_cluster_size))
  sc <- setClusters(sc, clust)

  cat(sprintf("[soupx] %s: autoEstCont (%d clusters)\n",
              args$sample_id, length(unique(clust))))
  sc <- autoEstCont(sc, doPlot = FALSE, verbose = TRUE)
} else {
  cat(sprintf("[soupx] %s: setContaminationFraction(%.4f) [manual]\n",
              args$sample_id, args$rho))
  sc <- setContaminationFraction(sc, args$rho)
  clust <- rep("manual_rho", ncol(toc))
}

rho_per_cell <- sc$metaData$rho
rho_mean <- mean(rho_per_cell)
cat(sprintf("[soupx] %s: rho mean=%.4f min=%.4f max=%.4f\n",
            args$sample_id, rho_mean, min(rho_per_cell), max(rho_per_cell)))

cat(sprintf("[soupx] %s: adjustCounts\n", args$sample_id))
corrected <- adjustCounts(sc, roundToInt = TRUE)

# ------------------------------------------------------------------
# Write outputs.
# ------------------------------------------------------------------
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

# Matrix Market sparse format (Python's scipy.io.mmread reads this).
writeMM(corrected, file.path(args$output_dir, "matrix.mtx"))

# Barcodes (one per cell column).
writeLines(colnames(corrected), file.path(args$output_dir, "barcodes.tsv"))

# Features (gene metadata). rowData typically has columns: ID, Symbol, Type.
feat_df <- as.data.frame(rowData(filtered))
# Make sure ID is the first column for downstream parsing.
if ("ID" %in% colnames(feat_df)) {
  feat_df <- feat_df[, c("ID", setdiff(colnames(feat_df), "ID"))]
}
write.table(feat_df, file.path(args$output_dir, "features.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Summary JSON.
n_before <- sum(toc)
n_after  <- sum(corrected)
pct_removed <- (n_before - n_after) / n_before * 100
elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

summary <- list(
  sample_id      = args$sample_id,
  rho_mean       = rho_mean,
  rho_min        = min(rho_per_cell),
  rho_max        = max(rho_per_cell),
  n_cells        = ncol(corrected),
  n_genes        = nrow(corrected),
  n_total_before = n_before,
  n_total_after  = n_after,
  pct_removed    = pct_removed,
  n_clusters     = length(unique(clust)),
  elapsed_sec    = elapsed_sec,
  mode           = if (is.na(args$rho)) "autoEst" else "manual"
)
write_json(summary, file.path(args$output_dir, "soupx_summary.json"),
           auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("[soupx] %s: DONE rho=%.4f removed=%.2f%% (%d/%d counts) time=%.1fs\n",
            args$sample_id, rho_mean, pct_removed,
            n_before - n_after, n_before, elapsed_sec))
