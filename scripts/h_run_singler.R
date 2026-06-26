#!/usr/bin/env Rscript
# h_run_singler.R -- SingleR label transfer (Spearman correlation), generic worker.
# The STAMP-analog for the human side: per-cell correlation of query against
# reference cell-type profiles. Compartment collapse is done Python-side (h09d);
# this worker just transfers whatever labels it's given.
#
# Inputs (all MatrixMarket genes x cells, symbols as feature rows, shared gene space
# is intersected internally):
#   --query-mtx / --query-genes / --query-barcodes   query (our placenta, lognorm)
#   --ref-mtx   / --ref-genes   / --ref-labels        reference (Vento-Tormo, lognorm)
#   --output    TSV: barcode, singler_label, delta_next (score margin), pruned_label
#
# Both matrices must be log-normalized (SingleR correlates on lognorm). de.method
# 'classic' (correlation) matches the STAMP Spearman approach.

suppressPackageStartupMessages({
  library(optparse); library(Matrix); library(SingleR); library(SummarizedExperiment)
  library(BiocParallel)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--query-mtx", type = "character"),
  make_option("--query-genes", type = "character"),
  make_option("--query-barcodes", type = "character"),
  make_option("--ref-mtx", type = "character"),
  make_option("--ref-genes", type = "character"),
  make_option("--ref-labels", type = "character"),
  make_option("--n-jobs", type = "integer", default = 8),
  make_option("--output", type = "character")
)))

read_mtx <- function(mtx, genes, barcodes = NULL) {
  m <- readMM(mtx)
  rownames(m) <- readLines(genes)
  if (!is.null(barcodes)) colnames(m) <- readLines(barcodes)
  as(m, "CsparseMatrix")
}

q <- read_mtx(opt$`query-mtx`, opt$`query-genes`, opt$`query-barcodes`)
r <- read_mtx(opt$`ref-mtx`, opt$`ref-genes`)
ref_labels <- readLines(opt$`ref-labels`)
stopifnot(length(ref_labels) == ncol(r))

# shared gene space (symbols); collapse duplicate symbols by first occurrence
q <- q[!duplicated(rownames(q)), , drop = FALSE]
r <- r[!duplicated(rownames(r)), , drop = FALSE]
common <- intersect(rownames(q), rownames(r))
if (length(common) < 200) stop(sprintf("only %d shared genes -- check symbol mapping", length(common)))
cat(sprintf("[singler] %d shared genes; query %d cells, ref %d cells\n",
            length(common), ncol(q), ncol(r)))
q <- q[common, , drop = FALSE]
r <- r[common, , drop = FALSE]

bp <- if (opt$`n-jobs` > 1) MulticoreParam(workers = opt$`n-jobs`) else SerialParam()
cat(sprintf("[singler] running with %d worker(s)\n", opt$`n-jobs`))
pred <- SingleR(test = q, ref = r, labels = ref_labels,
                de.method = "classic", BPPARAM = bp)   # correlation-based (Spearman), STAMP-like

out <- data.frame(
  barcode      = colnames(q),
  singler_label = pred$labels,
  pruned_label  = pred$pruned.labels,
  stringsAsFactors = FALSE
)
# score margin: top minus second-best tuning score (confidence proxy)
sc <- pred$scores
ord <- t(apply(sc, 1, function(x) sort(x, decreasing = TRUE)[1:2]))
out$delta_next <- ord[, 1] - ord[, 2]

write.table(out, opt$output, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[singler] wrote %d labels -> %s\n", nrow(out), opt$output))
