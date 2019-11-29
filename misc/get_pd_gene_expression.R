#!/usr/bin/env Rscript
library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("--genes"), type="character", default=NULL,
              help="File of PD genes and locations from bedtools intersect"),
  make_option(c("--expr"), type="character", default=NULL,
              help="Table of gene expression")
)
opt <- parse_args(OptionParser(option_list=option_list))
# opt = list(genes = "/Users/jeremys/work/opentargets/PD/Nalls_2019_genes.1Mb_window.gene_overlaps.tsv",
#            expr = "/Users/jeremys/work/opentargets/reference/tissueExpr/tissues.selected.tpm.tsv.gz")

trimGeneID = function(geneID) {
  gsub("\\..*", "", geneID)
}

genes.df = readr::read_tsv(opt$genes, col_names=T) %>%
  dplyr::mutate(gene_id = trimGeneID(gene_id))

expr.df = readr::read_tsv(opt$expr)

# Round values to 3 decimal places
roundToString = function(x) {
  if (is.na(x)) { return("") }
  sprintf("%.3g", round(x, digits = 3))
}
expr.df.short = cbind(expr.df[,1], apply(expr.df[,-1:-2], MARGIN=c(1, 2), FUN=roundToString))


getDist = function(pos, region_start, region_end) {
  if (pos < region_start) {
    pos - region_start
  } else if (pos > region_end) {
    pos - region_end
  } else {
    0
  }
}

genedist.df = genes.df %>%
  group_by(locus_number, gene_id) %>%
  mutate(dist = min(getDist(as.integer((windowStart+windowEnd)/2), gene_start, gene_end))) %>%
  ungroup() %>%
  mutate(chrInt = as.integer(gsub("chr", "", chr))) %>%
  arrange(chrInt, locus_number, dist) %>%
  select(-chrInt)

genedist.expr.df = genedist.df %>% dplyr::left_join(expr.df.short, by=c("gene_id" = "ensgene"))

write.table(genedist.expr.df, file="", sep="\t", quote=F, na="", row.names=F, col.names=T)

