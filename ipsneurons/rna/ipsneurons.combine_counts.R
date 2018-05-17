#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
samples.fname = args[1]
samples.dir = args[2]
output.fname = args[3]
counts.suffix = args[4]
if (is.na(counts.suffix)) {
  counts.suffix = ".counts.txt"
}

# Load featureCounts output
loadCounts <- function(sample_dir, sample_names, counts_suffix = ".counts.txt", sub_dir = TRUE) {
  matrix = c()
  for (i in c(1:length(sample_names))){
    if (sub_dir == TRUE){
      path = file.path(sample_dir, sample_names[i], paste(sample_names[i], counts_suffix, sep = ""))
    } else {
      path = file.path(sample_dir, paste(sample_names[i], counts_suffix, sep = ""))      
    }
    print(sample_names[i])
    table = readr::read_tsv(path, skip = 1, col_types = "cccccii")
    print(head(table))
    if (i == 1){
      matrix = table[,c(1,6,7)]
    }
    else{
      matrix = cbind(matrix, table[,7])
    }
  }
  colnames(matrix) = c("gene_id", "length", sample_names)
  return(matrix)
}

sample_names = read.table(samples.fname, sep ="\t", comment.char = "", stringsAsFactors = FALSE)[,1]

counts = loadCounts(samples.dir, sample_names, sub_dir = TRUE, counts_suffix = counts.suffix)

if (grepl(".gz$", output.fname)) {
  outfile = gzfile(output.fname, "w")
} else {
  outfile = file(output.fname, "w")
}
write.table(counts, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
close(outfile)

