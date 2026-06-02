#!/usr/bin/env Rscript
# run_propeller.R — composition analysis via speckle::propeller. Called by 08a
# via subprocess (same pattern as run_scdblfinder.R). Replaces scCODA.
#
# Why propeller over scCODA: scCODA's TF/TFP/arviz/matplotlib/numpy/setuptools
# stack is a dependency tar pit. propeller (speckle + limma) is a clean
# Bioconductor install, and its limma empirical-Bayes variance moderation is
# well-suited to the small n here (borrows information across cell types).
#
# Contract (CSV in, CSV out):
#   --counts  : CSV, rows = donors (index col 'donor_id'), columns = covariates
#               + one integer-count column per cell type.
#   --celltypes : comma-separated cell type column names (the count columns)
#   --covariates: comma-separated covariate column names (e.g. "group,sex,pool")
#   --test      : the factor being tested (e.g. "group")
#   --levels    : comma-separated "test,ref" for pairwise (e.g. "Early_Stress,Relaxed");
#                 empty => omnibus ANOVA across all levels of --test
#   --transform : "logit" (default) or "asin"
#   --out       : output CSV: [celltype, baseline_prop, log_fc_or_stat,
#                 statistic, pvalue, fdr, test_type]
#
# Covariate adjustment: confounders (sex, pool) are added as extra columns in the
# limma design matrix, AFTER the group columns (as speckle requires).

suppressPackageStartupMessages({
  library(optparse); library(speckle); library(limma)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts", type = "character"),
  make_option("--celltypes", type = "character"),
  make_option("--covariates", type = "character"),
  make_option("--test", type = "character"),
  make_option("--levels", type = "character", default = ""),
  make_option("--transform", type = "character", default = "logit"),
  make_option("--out", type = "character")
)))

df <- read.csv(opt$counts, row.names = 1, check.names = FALSE)
celltypes  <- strsplit(opt$celltypes, ",")[[1]]
covariates <- strsplit(opt$covariates, ",")[[1]]
test_fac   <- opt$test
levs       <- if (nchar(opt$levels) > 0) strsplit(opt$levels, ",")[[1]] else NULL

# Rebuild cell-level cluster/sample/covariate vectors from the per-donor count
# matrix (speckle works on cell-level factors). Each donor contributes its counts.
donors <- rownames(df)
counts_mat <- as.matrix(df[, celltypes, drop = FALSE])
storage.mode(counts_mat) <- "integer"

clusters <- character(0); sample <- character(0)
covdata  <- setNames(vector("list", length(covariates)), covariates)
for (cv in covariates) covdata[[cv]] <- character(0)

for (i in seq_along(donors)) {
  d <- donors[i]
  for (ct in celltypes) {
    n <- counts_mat[i, ct]
    if (n > 0) {
      clusters <- c(clusters, rep(ct, n))
      sample   <- c(sample,   rep(d,  n))
      for (cv in covariates) covdata[[cv]] <- c(covdata[[cv]], rep(as.character(df[i, cv]), n))
    }
  }
}
clusters <- factor(clusters); sample <- factor(sample)

# Transformed proportions (per sample)
props <- getTransformedProps(clusters = clusters, sample = sample, transform = opt$transform)

# Per-sample covariate table (one row per sample, aligned to props columns)
samp_levels <- colnames(props$TransformedProps)
covtab <- data.frame(sample = samp_levels, stringsAsFactors = FALSE)
# map each sample -> its (constant) covariate values
cov_cell <- data.frame(sample = sample, lapply(covdata, factor), check.names = FALSE)
cov_uniq <- cov_cell[!duplicated(cov_cell$sample), ]
rownames(cov_uniq) <- as.character(cov_uniq$sample)
cov_uniq <- cov_uniq[samp_levels, , drop = FALSE]

# Build design matrix: test factor first (no intercept), then confounders.
test_vec <- factor(cov_uniq[[test_fac]])
confounders <- setdiff(covariates, test_fac)
# drop confounders that don't vary (limma errors on constant columns).
# vapply (not sapply) guarantees a logical vector even when confounders is
# empty or length 1 — sapply can return a list and break the subscript.
if (length(confounders) > 0) {
  varies <- vapply(confounders,
                   function(cc) nlevels(factor(cov_uniq[[cc]])) > 1,
                   logical(1))
  confounders <- confounders[varies]
}

design_formula <- paste0("~ 0 + test_vec",
                         if (length(confounders)) paste0(" + ",
                           paste0("factor(cov_uniq[['", confounders, "']])", collapse = " + ")) else "")
design <- model.matrix(as.formula(design_formula))
# name the group columns cleanly
grp_cols <- grep("^test_vec", colnames(design))
colnames(design)[grp_cols] <- levels(test_vec)

if (is.null(levs)) {
  # Omnibus ANOVA across all levels of the test factor
  res <- propeller.anova(props, design = design, coef = grp_cols,
                         robust = TRUE, trend = FALSE, sort = FALSE)
  res$test_type <- "anova_omnibus"
  res$celltype <- rownames(res)
} else {
  # Pairwise t-test: contrast = test level - reference level
  ct_vec <- rep(0, ncol(design))
  ct_vec[which(colnames(design) == levs[1])] <-  1
  ct_vec[which(colnames(design) == levs[2])] <- -1
  res <- propeller.ttest(props, design = design, contrasts = ct_vec,
                         robust = TRUE, trend = FALSE, sort = FALSE)
  res$test_type <- paste0("ttest_", levs[1], "_vs_", levs[2])
  res$celltype <- rownames(res)
}

write.csv(res, opt$out, row.names = FALSE)
cat(sprintf("propeller OK: %d cell types, test=%s\n", nrow(res), res$test_type[1]))
