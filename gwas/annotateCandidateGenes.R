#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
genesFile = args[1]
rpkmsFile = args[2]

trimGeneID = function(geneID) {
  gsub("\\..*", "", geneID)
}

genes.df = readr::read_tsv(genesFile) %>%
  dplyr::mutate(geneID = trimGeneID(geneID))

rpkm.df = readr::read_tsv(rpkmsFile)

genes.df = genes.df %>% dplyr::left_join(rpkm.df, by=c("geneID" = "gene_id"))

write.table(genes.df, file="", sep="\t", quote=F, na="", row.names=F, col.names=T)
