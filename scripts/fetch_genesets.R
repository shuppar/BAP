#!/usr/bin/env Rscript
# fetch_genesets.R — export mouse MSigDB gene sets to a long TSV for Phase 8c.
#
# Uses msigdbr (mouse-native gene sets — no ortholog mapping; decoupler's Python
# get_resource('MSigDB', organism='mouse') is broken). Base R only besides
# msigdbr/optparse — no tidyverse needed.
#
# Collections (mouse codes):
#   MH               — hallmark (inflammation, OXPHOS, UPR/ER-stress, hypoxia...)
#   M2 / CP:REACTOME — Reactome canonical pathways
#   M5 / GO:BP       — GO Biological Process (glucocorticoid response, synaptic,
#                      neuroinflammation, mitochondrial respiration...)
#   M8               — cell type signatures from mouse scRNA-seq (esp. 7b subclusters)
#
# msigdbr API verified: db_species="MM", collection=, subcollection=; data ships
# in the msigdbdf package (a msigdbr dependency).
#
# Usage:  Rscript scripts/fetch_genesets.R --out refs/msigdb_mouse.tsv
# Output columns (tab-separated): collection, subcollection, gs_name, gene_symbol

suppressPackageStartupMessages({
  library(optparse); library(msigdbr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--out", type = "character", default = "refs/msigdb_mouse.tsv")
)))

# (collection, subcollection) pairs. NA subcollection = whole collection.
specs <- list(c("MH", NA), c("M2", "CP:REACTOME"), c("M5", "GO:BP"), c("M8", NA))

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)

parts <- list()
for (s in specs) {
  coll <- s[1]; subcoll <- s[2]
  msg <- if (is.na(subcoll)) coll else paste0(coll, " / ", subcoll)
  message("[genesets] fetching ", msg, " ...")
  df <- if (is.na(subcoll)) {
    msigdbr(db_species = "MM", species = "Mus musculus", collection = coll)
  } else {
    msigdbr(db_species = "MM", species = "Mus musculus",
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
