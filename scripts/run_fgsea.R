#!/usr/bin/env Rscript
# scripts/run_fgsea.R — fgsea-multilevel worker for 08c_pathways.py
#
# Replaces decoupler.mt.gsea() per (slice × collection). Called as subprocess.
# fgsea-multilevel (Korotkevich et al. 2021, bioRxiv 060012) estimates
# arbitrarily small p-values via adaptive multilevel splitting — fixes the
# 1/nperm permutation floor that produces bimodal FDRs in decoupler.
#
# Usage:
#   Rscript scripts/run_fgsea.R <ranks.tsv> <pathways.tsv> <out.tsv> \
#     [min_size=15] [max_size=500] [seed=42]
#
# Inputs:
#   ranks.tsv:    header `gene\tstat`. stat = Wald statistic (continuous, signed).
#                 Duplicate genes are collapsed by max |stat|.
#                 Non-finite stats are dropped.
#   pathways.tsv: header `pathway_name\tgene`. Long-form gene-set membership.
#                 Genes that don't appear in ranks are silently dropped by fgsea.
#
# Output (tab-separated, header):
#   pathway       — pathway name (matches input `pathway_name`)
#   ES            — enrichment score (raw running-sum max-deviation)
#   NES           — normalized enrichment score
#   pval          — continuous p-value from fgseaMultilevel (eps=0)
#   padj          — BH-adjusted p across the pathways passed in this call
#                   (= per-collection BH when caller invokes per-collection)
#   size          — number of pathway genes present in ranks
#   leading_edge  — "|"-separated gene symbols driving the enrichment

suppressPackageStartupMessages({
  library(fgsea)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript run_fgsea.R <ranks.tsv> <pathways.tsv> <out.tsv> ",
       "[min_size] [max_size] [seed]")
}
ranks_file    <- args[1]
pathways_file <- args[2]
out_file      <- args[3]
min_size <- if (length(args) >= 4L) as.integer(args[4]) else 15L
max_size <- if (length(args) >= 5L) as.integer(args[5]) else 500L
seed     <- if (length(args) >= 6L) as.integer(args[6]) else 42L

set.seed(seed)

write_empty <- function(path) {
  empty <- data.table(pathway = character(0), ES = numeric(0), NES = numeric(0),
                      pval = numeric(0), padj = numeric(0), size = integer(0),
                      leading_edge = character(0))
  fwrite(empty, path, sep = "\t")
  quit(save = "no", status = 0)
}

# --- ranks ---
ranks_df <- fread(ranks_file, header = TRUE)
if (!all(c("gene", "stat") %in% colnames(ranks_df))) {
  stop("ranks file must have header columns 'gene' and 'stat'")
}
ranks <- as.numeric(ranks_df$stat)
names(ranks) <- as.character(ranks_df$gene)
ranks <- ranks[is.finite(ranks)]
# Collapse duplicates by max |stat|
if (anyDuplicated(names(ranks))) {
  d <- data.table(gene = names(ranks), stat = ranks, abs_stat = abs(ranks))
  setorder(d, -abs_stat)
  d <- d[!duplicated(gene)]
  ranks <- d$stat
  names(ranks) <- d$gene
}
ranks <- sort(ranks, decreasing = TRUE)
if (length(ranks) < min_size) write_empty(out_file)

# --- pathways ---
path_df <- fread(pathways_file, header = TRUE)
if (!all(c("pathway_name", "gene") %in% colnames(path_df))) {
  stop("pathways file must have header columns 'pathway_name' and 'gene'")
}
pathways <- split(as.character(path_df$gene), as.character(path_df$pathway_name))
sz <- vapply(pathways, length, integer(1L))
pathways <- pathways[sz >= min_size & sz <= max_size]
if (length(pathways) == 0L) write_empty(out_file)

# --- run fgsea-multilevel (eps=0 for maximum p-value precision) ---
res <- tryCatch({
  fgsea(pathways    = pathways,
        stats       = ranks,
        minSize     = min_size,
        maxSize     = max_size,
        eps         = 0)
}, error = function(e) {
  stop(sprintf("fgsea failed (n_ranks=%d, n_pathways=%d): %s",
               length(ranks), length(pathways), conditionMessage(e)))
})

# Pipe-separate leading edge for portable TSV write
res$leading_edge <- vapply(res$leadingEdge,
                           function(x) paste(x, collapse = "|"),
                           FUN.VALUE = character(1L))
res$leadingEdge  <- NULL

# fgsea uses 'pathway'; keep canonical column order
setcolorder(res, c("pathway", "ES", "NES", "pval", "padj", "size", "leading_edge"))

fwrite(res, out_file, sep = "\t")
