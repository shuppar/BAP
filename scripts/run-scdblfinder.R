#!/usr/bin/env Rscript
# scDblFinder driver — called as a subprocess from src/doublets.py
#
# Usage (from Python):
#   subprocess.run([
#       "Rscript", "scripts/run-scdblfinder.R",
#       "--input", "/path/to/counts.mtx",
#       "--barcodes", "/path/to/barcodes.tsv",
#       "--genes", "/path/to/genes.tsv",
#       "--output", "/path/to/doublet_scores.tsv",
#       "--library", "lib1"
#   ], check=True)
#
# Inputs:
#   --input      : path to .mtx counts matrix (CellRanger style) OR a .rds SCE object
#   --barcodes   : path to barcodes.tsv (one barcode per line, matches matrix columns)
#   --genes      : path to genes.tsv (matches matrix rows)
#   --output     : path to write TSV with columns [barcode, doublet_score, doublet_class]
#   --library    : library/capture identifier (scDblFinder runs per library)
#
# This script does NOT modify the AnnData object. Python reads the output TSV
# and joins it back onto adata.obs by barcode.
#
# Rationale for subprocess pattern (vs rpy2):
#   - Process isolation: R crash doesn't take down Python
#   - Easier debugging: can run the R script standalone
#   - No rpy2 build-against-system-R fragility
#   - JSON/TSV-based contract is language-agnostic

suppressPackageStartupMessages({
  library(optparse)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(Matrix)
  library(data.table)
})

# -- CLI --
option_list <- list(
  make_option("--input",    type = "character", help = "Path to counts matrix (.mtx or .rds)"),
  make_option("--barcodes", type = "character", default = NULL, help = "barcodes.tsv (if .mtx input)"),
  make_option("--genes",    type = "character", default = NULL, help = "genes.tsv (if .mtx input)"),
  make_option("--output",   type = "character", help = "Output TSV path"),
  make_option("--library",  type = "character", default = "lib1", help = "Library ID for logging"),
  make_option("--seed",     type = "integer",   default = 42L,    help = "Random seed")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$input), !is.null(opt$output))

# -- Load counts --
message(sprintf("[scDblFinder] Loading %s...", opt$input))
if (grepl("\\.rds$", opt$input)) {
  sce <- readRDS(opt$input)
} else if (grepl("\\.mtx", opt$input)) {
  stopifnot(!is.null(opt$barcodes), !is.null(opt$genes))
  counts <- as(Matrix::readMM(opt$input), "CsparseMatrix")
  barcodes <- readLines(opt$barcodes)
  genes <- readLines(opt$genes)
  rownames(counts) <- genes
  colnames(counts) <- barcodes
  sce <- SingleCellExperiment(assays = list(counts = counts))
} else {
  stop("Unsupported input format. Use .mtx or .rds.")
}

message(sprintf("[scDblFinder] library=%s, n_cells=%d, n_genes=%d",
                opt$library, ncol(sce), nrow(sce)))

# -- Run scDblFinder --
set.seed(opt$seed)
sce <- scDblFinder(sce, verbose = FALSE)

# -- Write results --
results <- data.table(
  barcode        = colnames(sce),
  doublet_score  = sce$scDblFinder.score,
  doublet_class  = as.character(sce$scDblFinder.class),
  library        = opt$library
)
fwrite(results, opt$output, sep = "\t")

n_dbl <- sum(results$doublet_class == "doublet")
message(sprintf("[scDblFinder] %d / %d cells flagged as doublets (%.1f%%)",
                n_dbl, nrow(results), 100 * n_dbl / nrow(results)))
message(sprintf("[scDblFinder] Wrote %s", opt$output))
