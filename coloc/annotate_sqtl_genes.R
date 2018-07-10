library(data.table)
library(dplyr)
library(rtracklayer)
options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly = TRUE)
sqtl_introns_fname = args[1]
sqtl_col = args[2]
exons_fname = args[3]
genes_fname = args[4]

sqtl.introns <- readr::read_tsv(sqtl_introns_fname) %>%
  dplyr::rename(clusterID = feature)

#Inefficient way of getting the chr, start, end, but gets the job done
sqtl.introns$chr = sapply(sqtl.introns$clusterID, function(x) {strsplit(x, split=":", fixed=T)[[1]][1]})
sqtl.introns$start = as.integer(sapply(sqtl.introns$clusterID, function(x) {strsplit(x, split=":", fixed=T)[[1]][2]}))
sqtl.introns$end = as.integer(sapply(sqtl.introns$clusterID, function(x) {strsplit(x, split=":", fixed=T)[[1]][3]}))
sqtl.introns$cluster = sapply(sqtl.introns$clusterID, function(x) {strsplit(x, split=":", fixed=T)[[1]][4]})

sqtl.coords = sqtl.introns %>% dplyr::select(chr, start, end, cluster) %>%
  dplyr::mutate(start = start - 50, end = end + 50) %>%
  dplyr::arrange(chr, start)

write.table(sqtl.coords, file="sqtl.clusters.expanded.bed", quote=F, row.names=F, col.names=F, sep="\t")

overlapFname = paste0(sqtl_introns_fname, ".exons.txt")
cmd = sprintf("bedtools closest -d -a sqtl.clusters.expanded.bed -b %s | cut -f 1-4,8 | sort | uniq > %s", exons_fname, overlapFname)
output = system(cmd, intern = T)

library(data.table)
sqtl.overlaps = data.table(readr::read_tsv(overlapFname, col_name=c("chr", "start", "end", "cluster", "geneid"), col_types="ciicc"))
sqtl.overlaps$start = sqtl.overlaps$start + 50
sqtl.overlaps$end = sqtl.overlaps$end - 50

# Annotated the ensembl IDs with gene symbols
gene.df = readr::read_tsv(genes_fname)
gene.df$ensemblID = gsub("\\.[\\d]+", "", gene.df$ensemblID, perl=T)
sqtl.overlaps = sqtl.overlaps %>% dplyr::left_join(gene.df, by=c("geneid"="ensemblID")) %>% as.data.table()

sqtl.overlaps = sqtl.overlaps[, .(paste(c(geneid), collapse=","), paste(c(geneSymbol), collapse=",")), by = .(chr, start, end, cluster)] %>%
  dplyr::rename(geneid = V1, geneSymbol = V2)
sqtl.overlaps$clusterKey = paste(sqtl.overlaps$chr, sqtl.overlaps$start, sqtl.overlaps$end, sqtl.overlaps$cluster, sep=":")

sqtl.introns.new = sqtl.introns %>% dplyr::left_join(sqtl.overlaps %>% dplyr::select(clusterKey, geneid, geneSymbol), by=c("clusterID" = "clusterKey"))
write.table(sqtl.introns.new %>% dplyr::select(-start, -end, -cluster),
            file=paste0(sqtl_introns_fname, ".ann.txt"),
            quote=F, row.names=F, col.names=T, sep="\t")

