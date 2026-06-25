#!/usr/bin/env Rscript
# run_soupx_from_raw.R -- knee cell-call + SoupX on a RAW-only 10x h5.
# For human placenta (Gunter-Rahman GSE271976): GEO shipped raw_feature_bc_matrix only,
# so we cell-call first (mouse run_soupx.R expects cellranger filtered+raw; we have raw only),
# then run SoupX with the same scran/autoEstCont logic.
#
# Per sample: read raw h5 -> barcodeRanks knee/inflection -> call cells
#   -> SoupChannel(tod=all droplets, toc=called cells) -> quickCluster -> autoEstCont
#   -> adjustCounts -> write gzipped MTX trio (cellranger-style) + summary.json.
#
# Usage:
#   Rscript scripts/run_soupx_from_raw.R --h5 <raw.h5> --out-dir <dir> --sample-id <id> \
#           [--cutoff knee|inflection] [--rho <manual fraction>]

suppressPackageStartupMessages({
  library(DropletUtils); library(SoupX); library(scran); library(Matrix); library(jsonlite)
  library(DelayedArray)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
h5      <- getarg("--h5")
out_dir <- getarg("--out-dir")
sid     <- getarg("--sample-id")
cutoff  <- getarg("--cutoff", "knee")
rho_man <- getarg("--rho", NA)

if (is.null(h5) || is.null(out_dir) || is.null(sid)) {
  stop("need --h5 --out-dir --sample-id")
}
if (!file.exists(h5)) stop(sprintf("[%s] h5 not found: %s", sid, h5))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- read raw droplets ---------------------------------------------------------
sce       <- read10xCounts(h5, col.names = TRUE)
m         <- counts(sce)
m         <- as(m, "CsparseMatrix")            # realize HDF5/Delayed -> in-memory dgCMatrix (SoupX needs plain sparse)
gene_ids  <- as.character(rowData(sce)$ID)
gene_syms <- make.unique(as.character(rowData(sce)$Symbol))
rownames(m) <- gene_syms                       # symbols as row names (match Python var_names)

tot <- Matrix::colSums(m)
keep <- tot > 0
m    <- m[, keep]; tot <- tot[keep]            # drop all-zero droplets early
n_droplets <- ncol(m)

# --- knee cell-call ------------------------------------------------------------
br   <- barcodeRanks(m)
knee <- metadata(br)$knee
infl <- metadata(br)$inflection
thr  <- if (identical(cutoff, "inflection")) infl else knee
cells <- colnames(m)[tot >= thr]
if (length(cells) < 100) {
  stop(sprintf("[%s] only %d cells at %s=%.0f -- cutoff too aggressive (knee=%.0f infl=%.0f)",
               sid, length(cells), cutoff, thr, knee, infl))
}

toc <- m[, cells]                              # called cells (table of counts)
tod <- m                                       # all droplets (table of droplets / background)

# --- SoupX ---------------------------------------------------------------------
sc <- SoupChannel(tod, toc)
set.seed(1)
clus <- tryCatch(quickCluster(toc, min.size = 50), error = function(e) NULL)
if (!is.null(clus)) sc <- setClusters(sc, as.character(clus))

if (!is.na(rho_man)) {
  sc  <- setContaminationFraction(sc, as.numeric(rho_man))
  rho <- as.numeric(rho_man)
  rho_method <- "manual"
} else {
  est <- tryCatch(autoEstCont(sc, doPlot = FALSE), error = function(e) NULL)
  if (is.null(est)) {                          # fallback: fixed 10% if autoEst fails
    sc  <- setContaminationFraction(sc, 0.10)
    rho <- 0.10; rho_method <- "fallback_0.10"
  } else {
    sc  <- est; rho <- mean(sc$metaData$rho); rho_method <- "autoEstCont"
  }
}

out <- adjustCounts(sc, roundToInt = TRUE)     # rows == rownames(m) order, preserved
pct_removed <- 100 * (1 - sum(out) / sum(toc))

# --- write gzipped MTX trio + summary -----------------------------------------
mtx_path  <- file.path(out_dir, "matrix.mtx")
bc_path   <- file.path(out_dir, "barcodes.tsv")
feat_path <- file.path(out_dir, "features.tsv")

writeMM(out, mtx_path)
writeLines(colnames(out), bc_path)
# features.tsv: gene_id \t symbol \t feature_type (cellranger order, aligned to out rows)
write.table(
  data.frame(gene_ids, rownames(out), "Gene Expression"),
  feat_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
for (p in c(mtx_path, bc_path, feat_path)) system2("gzip", c("-f", p))

summary <- list(
  sample_id = sid, cutoff = cutoff, knee = knee, inflection = infl,
  n_droplets = n_droplets, n_cells = ncol(toc),
  rho_mean = rho, rho_method = rho_method, pct_removed = pct_removed
)
write_json(summary, file.path(out_dir, "summary.json"), auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("[%s] cells=%d (cutoff=%s thr=%.0f; knee=%.0f infl=%.0f) rho=%.3f (%s) removed=%.1f%%\n",
            sid, ncol(toc), cutoff, thr, knee, infl, rho, rho_method, pct_removed))
