#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
missing.fname = args[1]
stats.fname = args[2]
AD=F

missing.df = readr::read_tsv(missing.fname)
if (AD) {
  stats.df = readr::read_tsv(stats.fname) %>%
    dplyr::select(SNP, CHR, BP, A1, A2, GWAS_BETA, GWAS_P, GWAX_P, META_P)
} else {
  stats.df = readr::read_tsv(stats.fname) %>%
    dplyr::select(SNP, CHR, BP, ALLELE1, ALLELE2, EFFECT, P.META)
}

df = missing.df %>% dplyr::left_join(stats.df, by=c("ld_snp" = "SNP"))
write.table(df, file="", row.names=F, col.names=T, sep="\t", quote=F, na="")
