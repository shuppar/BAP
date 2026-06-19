#!/usr/bin/env Rscript
# run_deseq2_lrt.R — DESeq2 likelihood-ratio test (Phase 8b subprocess worker)
#
# Called from scripts/08b_de.py. Two contrast types route here:
#   omnibus_3group_per_age : full=~sex+pool+group  reduced=~sex+pool        (df=2)
#   group_x_age_interaction: full=~sex+pool+group*age  reduced=~sex+pool+group+age (df=4)
#
# Why R: PyDESeq2 only exposes Wald in its public API. A joint LRT (testing
# several coefficients at once against a chi-squared distribution) needs the
# per-gene log-likelihoods of the full and reduced fits — internal arrays in
# PyDESeq2 whose layout shifted between 0.4 and 0.5. R's DESeq2 has done this
# in one line since 2014 (Love, Huber, Anders 2014, Genome Biology):
#     DESeq(dds, test="LRT", reduced=~...)
# So we delegate. Matches the project's existing R-subprocess pattern
# (run_propeller.R, run_soupx.R, run_scdblfinder.R).
#
# Reference level: group=Relaxed (positive log2FC = upregulated in stress) —
# LRT itself is invariant to the reference level, but we set it for consistency
# with the Wald (PyDESeq2) side. log2FC + lfcSE are NOT returned: in an LRT
# they refer to one coefficient among several being tested jointly, so a single
# LFC would mislead the reader. Python sets log2FC=NaN on every LRT row.
#
# I/O contract with 08b_de.py:
#   counts CSV : rows = gene (header = "gene" or unnamed), cols = donor
#   meta   CSV : rows = donor (first column = donor_id), cols = covariate
#   out    CSV : columns = gene, baseMean, stat, pvalue, padj  (Python reads
#                this with pd.read_csv(out, index_col=0))
#
# Usage:
#   Rscript scripts/run_deseq2_lrt.R \
#     --counts  counts.csv \
#     --meta    meta.csv \
#     --full    "~sex+pool+group" \
#     --reduced "~sex+pool" \
#     --out     results.csv \
#     --label   "omnibus_3group"      # informational; goes to stderr only
#
# Install (one-time, on the workstation):
#   Rscript -e 'BiocManager::install("DESeq2", ask=FALSE, update=FALSE)'

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
})

# ---- CLI ------------------------------------------------------------------

option_list <- list(
  make_option("--counts",  type="character", help="counts CSV (rows=gene, cols=donor)"),
  make_option("--meta",    type="character", help="metadata CSV (rows=donor, cols=covariate)"),
  make_option("--full",    type="character", help="full design formula, e.g. '~sex+pool+group*age'"),
  make_option("--reduced", type="character", help="reduced design formula"),
  make_option("--out",     type="character", help="output CSV path"),
  make_option("--label",   type="character", default="LRT",
              help="informational label written to stderr (default: LRT)")
)
opt <- parse_args(OptionParser(option_list=option_list))

req <- c("counts", "meta", "full", "reduced", "out")
missing_args <- req[sapply(req, function(x) is.null(opt[[x]]))]
if (length(missing_args)) {
  message("ERROR: missing required arg(s): ", paste(missing_args, collapse = ", "))
  quit(status = 2)
}

# ---- read inputs ----------------------------------------------------------

counts <- tryCatch(
  as.matrix(read.csv(opt$counts, row.names = 1, check.names = FALSE)),
  error = function(e) { message("ERROR reading counts: ", conditionMessage(e)); quit(status = 1) }
)
storage.mode(counts) <- "integer"          # DESeq2 requires integer counts

meta <- tryCatch(
  read.csv(opt$meta, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE),
  error = function(e) { message("ERROR reading meta: ", conditionMessage(e)); quit(status = 1) }
)

if (!identical(colnames(counts), rownames(meta))) {
  # try to reorder meta to match counts columns (Python writes both from the
  # same sorted donor list, but be defensive)
  if (setequal(colnames(counts), rownames(meta))) {
    meta <- meta[colnames(counts), , drop = FALSE]
  } else {
    message("ERROR: counts columns and meta rows don't match.")
    message("  counts: ", paste(head(colnames(counts), 5), collapse = ","), "...")
    message("  meta  : ", paste(head(rownames(meta), 5), collapse = ","), "...")
    quit(status = 1)
  }
}

# ---- coerce factor columns + set reference levels -------------------------
# LRT is invariant to reference levels (it tests the model as a whole), but
# we set them for consistency with the Wald (PyDESeq2) path elsewhere in 8b.

factor_cols <- c("group", "age", "sex", "pool")
for (cn in intersect(factor_cols, colnames(meta))) {
  meta[[cn]] <- factor(meta[[cn]])
}
if ("group" %in% colnames(meta) && "Relaxed" %in% levels(meta$group)) {
  meta$group <- relevel(meta$group, ref = "Relaxed")
}

full_f    <- as.formula(opt$full)
reduced_f <- as.formula(opt$reduced)

# guard: identical designs -> LRT is undefined
if (identical(opt$full, opt$reduced)) {
  message("ERROR: full and reduced formulas are identical -> LRT undefined.")
  quit(status = 1)
}

message(sprintf("[%s] full=%s  reduced=%s  n_genes=%d  n_donors=%d",
                opt$label, opt$full, opt$reduced, nrow(counts), ncol(counts)))

# ---- build + fit ----------------------------------------------------------

dds <- tryCatch(
  DESeqDataSetFromMatrix(countData = counts, colData = meta, design = full_f),
  error = function(e) {
    message("ERROR: DESeqDataSetFromMatrix failed: ", conditionMessage(e))
    quit(status = 1)
  }
)

# Drop genes with zero counts across all donors (DESeq2 handles them but they
# add noise to mean-dispersion estimation at this n).
keep <- rowSums(counts(dds)) > 0
dds  <- dds[keep, ]

dds <- tryCatch(
  DESeq(dds, test = "LRT", reduced = reduced_f, quiet = TRUE),
  error = function(e) {
    message("ERROR: DESeq() LRT failed: ", conditionMessage(e))
    quit(status = 1)
  }
)

res <- results(dds)
df  <- as.data.frame(res)

# ---- write results --------------------------------------------------------
# Columns Python reads: gene (becomes the index via pd.read_csv(...,
# index_col=0)), baseMean, stat (LRT chi-squared), pvalue, padj. log2FC +
# lfcSE are NOT exported — see header comment.

out <- data.frame(
  gene     = rownames(df),
  baseMean = df$baseMean,
  stat     = df$stat,        # LRT chi-squared statistic
  pvalue   = df$pvalue,
  padj     = df$padj,
  stringsAsFactors = FALSE
)

write.csv(out, opt$out, row.names = FALSE)
n_sig <- sum(out$padj < 0.05, na.rm = TRUE)
message(sprintf("[%s] wrote %d genes (%d at padj<0.05) -> %s",
                opt$label, nrow(out), n_sig, opt$out))
