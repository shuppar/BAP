#!/usr/bin/env Rscript
# h_fetch_genesets.R -- export HUMAN MSigDB gene sets to a long TSV for Phase 9 GSEA (h09g).
# Human analog of fetch_genesets.R. Same output schema so load_genesets_tsv reads it unchanged.
#
# Collections (human codes; mouse->human equivalents):
#   H               -- hallmark              (mouse MH)
#   C2 / CP:REACTOME -- Reactome canonical    (mouse M2)
#   C5 / GO:BP      -- GO Biological Process  (mouse M5)
# C8 (cell-type signatures, = mouse M8) is OMITTED -- 8c drops M8 from plots as noise.
#
# Usage:  Rscript scripts/h_fetch_genesets.R --out refs/msigdb_human.tsv
# Output columns (tab-separated): collection, subcollection, gs_name, gene_symbol

suppressPackageStartupMessages({
  library(optparse); library(msigdbr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--out", type = "character", default = "refs/msigdb_human.tsv")
)))

# (collection, subcollection) pairs. NA subcollection = whole collection.
specs <- list(c("H", NA), c("C2", "CP:REACTOME"), c("C5", "GO:BP"))

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)

parts <- list()
for (s in specs) {
  coll <- s[1]; subcoll <- s[2]
  msg <- if (is.na(subcoll)) coll else paste0(coll, " / ", subcoll)
  message("[genesets] fetching ", msg, " ...")
  df <- if (is.na(subcoll)) {
    msigdbr(db_species = "HS", species = "Homo sapiens", collection = coll)
  } else {
    msigdbr(db_species = "HS", species = "Homo sapiens",
            collection = coll, subcollection = subcoll)
  }
  out <- data.frame(
    collection    = coll,
    subcollection = if (is.na(subcoll)) "" else subcoll,
    gs_name       = df$gs_name,
    gene_symbol   = df$gene_symbol,
    stringsAsFactors = FALSE
  )
  out <- unique(out)
  message("    ", length(unique(out$gs_name)), " sets, ", nrow(out), " gene rows")
  parts[[msg]] <- out
}

final <- unique(do.call(rbind, parts))
write.table(final, file = opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
message("\n[genesets] wrote ", nrow(final), " rows, ",
        length(unique(final$gs_name)), " gene sets -> ", opt$out)
message("[genesets] collections: ", paste(unique(final$collection), collapse = ", "))
