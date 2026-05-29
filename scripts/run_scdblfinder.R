#!/usr/bin/env Rscript
# run_scdblfinder.R — run scDblFinder on a per-pool combined matrix.
#
# Called by scripts/03_doublets.py as a subprocess. Standalone-runnable for
# debugging:
#
#   Rscript scripts/run_scdblfinder.R \
#     --matrix /tmp/pool1_counts.mtx \
#     --barcodes /tmp/pool1_barcodes.tsv \
#     --features /tmp/pool1_features.tsv \
#     --samples /tmp/pool1_samples.tsv \
#     --output /tmp/pool1_doublets.tsv
#
# Inputs (all required):
#   --matrix     MatrixMarket (.mtx) genes x cells (Scanpy writes cells x genes,
#                Python side handles the transpose before writing)
#   --barcodes   one barcode per line (cell IDs in matrix column order)
#   --features   one feature/gene per line (matrix row order)
#   --samples    one sample_id per line, same order as barcodes — passed to
#                scDblFinder as `samples=` so simulated doublets respect
#                sample boundaries within the pool
#   --output     output TSV: barcode, sample_id, doublet_score, doublet_class
#
# Optional:
#   --seed       RNG seed (default 42)

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(BiocParallel)
})

opt_list <- list(
  make_option("--matrix",   type="character"),
  make_option("--barcodes", type="character"),
  make_option("--features", type="character"),
  make_option("--samples",  type="character"),
  make_option("--output",   type="character"),
  make_option("--seed",     type="integer", default=42L)
)
opt <- parse_args(OptionParser(option_list=opt_list))

stopifnot(file.exists(opt$matrix), file.exists(opt$barcodes),
          file.exists(opt$features), file.exists(opt$samples))

set.seed(opt$seed)

cat(sprintf("[R] loading matrix: %s\n", opt$matrix))
m <- readMM(opt$matrix)                   # genes x cells
barcodes <- readLines(opt$barcodes)
features <- readLines(opt$features)
samples  <- readLines(opt$samples)

stopifnot(ncol(m) == length(barcodes),
          nrow(m) == length(features),
          length(samples) == length(barcodes))

rownames(m) <- features
colnames(m) <- barcodes
sce <- SingleCellExperiment(assays = list(counts = m))

cat(sprintf("[R] running scDblFinder on %d cells, %d samples...\n",
            ncol(sce), length(unique(samples))))
# samples= ensures simulated doublets respect sample boundaries (per scDblFinder docs);
# critical for pooled/multiplexed data where cells from different samples shouldn't
# be combined into artificial doublets.
sce <- scDblFinder(sce, samples = samples, BPPARAM = SerialParam())

out <- data.frame(
  barcode        = colnames(sce),
  sample_id      = samples,
  doublet_score  = sce$scDblFinder.score,
  doublet_class  = as.character(sce$scDblFinder.class),
  stringsAsFactors = FALSE
)
write.table(out, file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE)

n_dbl <- sum(out$doublet_class == "doublet")
cat(sprintf("[R] done. doublets: %d / %d (%.2f%%)\n",
            n_dbl, nrow(out), 100 * n_dbl / nrow(out)))
cat(sprintf("[R] wrote %s\n", opt$output))
