#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
file1 = args[1]
file2 = args[2]

df1 = readr::read_tsv(file1) %>% dplyr::mutate(population_r2 = sprintf("%s_r2_%.2f", population_name, r2))
df2 = readr::read_tsv(file2) %>% dplyr::mutate(population_r2 = sprintf("%s_r2_%.2f", population_name, r2))

df = rbind(df1, df2) %>% dplyr::group_by(variation1, variation2) %>% dplyr::summarise(population_r2=paste0(population_r2, collapse = ","))
write.table(df, file="", row.names=F, col.names=T, sep="\t", quote=F, na="")
