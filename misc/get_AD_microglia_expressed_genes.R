library(tidyverse)

genes.df = readr::read_tsv("genes/AD.loci.1Mb_window.gene_overlaps.rpkms.tsv") %>%
  select(leadSNP, locus, geneID, symbol, chr, gene_start, gene_end, microglia_gaffney_rpkm = microglia_gaffney)

trimGeneID = function(geneID) {
  gsub("\\..*", "", geneID)
}

genedist.df = readr::read_tsv("genes/AD.loci.1Mb_window.genedist.tsv") %>%
  select(locus, geneID, lead_pos, region_start, region_end, dist) %>%
  dplyr::mutate(geneID = trimGeneID(geneID))

genes.ann.df = genes.df %>%
  left_join(genedist.df, by=c("locus", "geneID")) %>%
  arrange(locus, desc(microglia_gaffney_rpkm)) %>%
  group_by(locus) %>%
  mutate(locus_rank = row_number())


write.table(genes.ann.df, file="genes/AD.loci.genes.microglia_expression.tsv",
            sep="\t", quote=F, na="", row.names=F, col.names=T)

