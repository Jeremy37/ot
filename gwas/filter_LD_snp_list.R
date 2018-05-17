#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
gwas.fname = args[1]
ldsnps.fname = args[2]

gwas.df = readr::read_tsv(gwas.fname)
ldsnps.df = readr::read_tsv(ldsnps.fname)

gwas.lead.df = gwas.df %>% dplyr::filter(snp %in% ldsnps.df$variation1)
ldsnps.df = ldsnps.df %>% dplyr::left_join(gwas.lead.df %>% dplyr::select(locus, snp), by=c("variation1"="snp"))

leadSnps.df = ldsnps.df %>% dplyr::filter(!duplicated(variation1)) %>% as.data.frame()
leadSnps.df$variation2 = leadSnps.df$variation1
leadSnps.df$population_r2 = 1
rownames(leadSnps.df) = leadSnps.df$variation1
leadSnps.df$locusRank = 1:nrow(leadSnps.df)
ldsnps.df$locusRank = leadSnps.df[ldsnps.df$variation1,]$locusRank

for (curlocus in unique(gwas.lead.df$locus)) {
  ldsnps.df = ldsnps.df %>% dplyr::filter(locus != curlocus |
                                          (!variation2 %in% (gwas.df %>% filter(locus == curlocus) %>% .$snp)))
}

ldsnps.df = rbind(leadSnps.df, ldsnps.df) %>%
  dplyr::arrange(locusRank) %>%
  dplyr::select(-locusRank)


write.table(ldsnps.df %>% dplyr::select(locus, leadsnp=variation1, ld_snp=variation2, population_r2),
            file="", sep="\t", quote=F, na="", row.names=F, col.names=T)
