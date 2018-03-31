#!/usr/bin/env Rscript
library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
inFname = args[1]
mafFname = args[2]
outFname = args[3]

snp.df = readr::read_tsv(inFname, col_names = c("chr", "pos", "snp"))
maf.df = readr::read_tsv(mafFname, col_names = c("snp", "maf"))

snp.df = snp.df %>% dplyr::left_join(maf.df, by="snp") %>%
  dplyr::mutate(ref = "-") %>%
  dplyr::mutate(alt = "-") %>%
  dplyr::select(chr, pos, snp, ref, alt, maf)

gzf = gzfile(outFname, open = "w")
readr::write_tsv(snp.df, gzf, col_names = F, na = "")
close(gzf)
