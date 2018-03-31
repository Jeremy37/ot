#!/usr/bin/env Rscript
library(tidyverse)
root = "/Users/jeremys/work/opentargets"

# Map gene symbols to Ensembl gene IDs
hgnc.df = readr::read_tsv(file.path(root, "reference/hgnc_complete_set.txt"))
hgnc.df = hgnc.df %>% dplyr::select(symbol, alias_symbol, prev_symbol, ensembl_gene_id)

# Make a column for each potential alias of a gene
for (i in 1:20) {
  hgnc.df[,paste0("alias", i)] = sapply(hgnc.df$alias_symbol, function(x) {strsplit(x, split="|", fixed=T)[[1]][i]})
}
# Do the same for previous symbols of a gene
for (i in 1:20) {
  hgnc.df[,paste0("prev_symbol", i)] = sapply(hgnc.df$prev_symbol, function(x) {strsplit(x, split="|", fixed=T)[[1]][i]})
}

hgnc.gather = hgnc.df %>% tidyr::gather(symbol, alias1:alias20, prev_symbol1:prev_symbol20, key="symbol_type", value="symbol") %>%
  dplyr::select(ensembl_gene_id, symbol_type, symbol) %>%
  na.omit()

# Reorder to put the canonical symbol first, so that when joining this
# table to others you get those first (which is best if you then remove dups)
symbol_types = c("symbol", "alias1", "prev_symbol1", unique(as.character(hgnc.gather$symbol_type)))
symbol_types = symbol_types[!duplicated(symbol_types)]
hgnc.gather$symbol_type = factor(as.character(hgnc.gather$symbol_type), levels = symbol_types)

hgnc.gather = hgnc.gather %>% arrange(symbol_type, symbol)

readr::write_tsv(hgnc.gather, path=file.path(root, "reference/hgnc.ensembl.map.txt"), col_names = T)
