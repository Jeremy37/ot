#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
regionsFile = args[1]
fastqMetaFile = args[2]
readLenParam = args[3]
sdParam = args[4]

if (is.null(readLenParam) | is.na(readLenParam)) {
  readLenParam = 150
} else {
  readLenParam = strtoi(readLenParam)
}
if (is.null(sdParam) | is.na(sdParam)) {
  sdParam = 20
} else {
  sdParam = strtoi(sdParam)
}

regions.df = readr::read_tsv(regionsFile)

meta.df = readr::read_tsv(fastqMetaFile)

getMinOverlap = function(amplicon_size, read_len) {
  minOverlap = 10
  if (2*read_len - amplicon_size < 10) {
    minOverlap = 2*read_len - amplicon_size - 1
    if (minOverlap < 4) {
      minOverlap = 4
    }
  }
  minOverlap
}
if (any(!meta.df$name %in% regions.df$name)) {
  missingNames = meta.df$name[!meta.df$name %in% regions.df$name]
  stop(sprintf("Some region names in metadata file %s are not found in regions file %s:\n%s",
               fastqMetaFile, regionsFile, paste(missingNames, collapse=", ")))
}
df = meta.df %>% dplyr::left_join(regions.df, by="name") %>%
  dplyr::mutate(amplicon_size = end - start + 1) %>%
  group_by(name, replicate_full) %>%
  dplyr::mutate(read_len = readLenParam, amplicon_sd = sdParam,
                expected_overlap = max(5, 2*read_len - amplicon_size)) %>%
  dplyr::mutate(output_file = paste0(replicate_full, ".extendedFrags.fastq.gz"),
                min_overlap = getMinOverlap(amplicon_size, read_len),
                max_mismatch_dens = min(0.15, 0.1 + 1 / expected_overlap)) %>%
  ungroup()

df$max_mismatch_dens = sprintf("%.3g", df$max_mismatch_dens)

df = df %>% dplyr::select(replicate_full, Read_1_file, Read_2_file, read_len, amplicon_size, amplicon_sd, min_overlap, max_mismatch_dens, output_file)

write.table(df, file="", sep="\t", quote=F, row.names=F, col.names=T)
