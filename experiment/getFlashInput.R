#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
regionsFile = args[1]
fastqMetaFile = args[2]

regions.df = readr::read_tsv(regionsFile)

meta.df = readr::read_tsv(fastqMetaFile)

getMinOverlap = function(amplicon_size, read_len) {
  minOverlap = 10
  if (2*read_len - amplicon_size < 10) {
    minOverlap = 2*read_len - amplicon_size
    if (minOverlap < 0) {
      minOverlap = 4
    }
  }
  minOverlap
}
df = meta.df %>% dplyr::left_join(regions.df, by="name") %>%
  dplyr::rename(amplicon_size = end) %>%
  dplyr::mutate(read_len = 150, amplicon_sd = 20,
                expected_overlap = 2*read_len - amplicon_size) %>%
  dplyr::mutate(output_file = paste0(replicate, ".extendedFrags.fastq.gz"))

df$min_overlap = sapply(1:nrow(df), FUN = function(i) getMinOverlap(df$amplicon_size[i], df$read_len[i]))
df$max_mismatch_dens = 0.05 + 1 / df$expected_overlap
df$max_mismatch_dens = sprintf("%.3g", df$max_mismatch_dens)

df = df %>% dplyr::select(replicate, Read_1_file, Read_2_file, read_len, amplicon_size, amplicon_sd, min_overlap, max_mismatch_dens, output_file)

write.table(df, file="", sep="\t", quote=F, row.names=F, col.names=T)
