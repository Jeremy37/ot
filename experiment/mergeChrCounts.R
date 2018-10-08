#!/usr/bin/env Rscript
library(tidyverse)
options(stringsAsFactors=F)

args <- commandArgs(trailingOnly = TRUE)
fileList = args[1]
chrCountFiles = read.delim(fileList, header = F)[[1]]

counts.df = data.frame(file = chrCountFiles[1], chr1 = 0)
curFileIndex = 1

for (f in chrCountFiles) {
  df = read.delim(file = f, header = F, colClasses = c("numeric", "character"), col.names = c("count", "chr"))
  counts.df[curFileIndex,]$file = f
  for (i in 1:nrow(df)) {
    chrName = df[i,]$chr
    if (!chrName %in% colnames(counts.df)) {
      counts.df[,chrName] = NA
    }
    counts.df[curFileIndex, chrName] = df[i,]$count
  }
  curFileIndex = curFileIndex + 1
}

counts.df[is.na(counts.df)] = 0
write.table(counts.df[, order(colnames(counts.df))], file="", col.names=T, row.names=F, quote=F, sep="\t")
