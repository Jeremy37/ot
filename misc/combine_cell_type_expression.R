#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
root = "/Users/jeremys/work/opentargets"

###########################################################################
# Barres lab primary brain sorted cell types
barres.meta = readr::read_tsv(file.path(root, "reference/Barres_RNAseq/barres_rnaseq.human.metadata.txt"))
barres.df = readr::read_tsv(file.path(root, "reference/Barres_RNAseq/barres_rnaseq.human.txt"))

barres.celltypeAvg.df = data.frame(symbol = barres.df$Gene)
for (celltype in unique(barres.meta$celltype)) {
  sampleIDs = barres.meta$sampleID[barres.meta$celltype == celltype]
  celltype.df = barres.df %>% dplyr::select(one_of(sampleIDs))
  barres.celltypeAvg.df = cbind(barres.celltypeAvg.df, rowMeans(celltype.df))
}
colnames(barres.celltypeAvg.df)[-1] = paste0("barres.", unique(barres.meta$celltype))


###########################################################################
# iPSC-derived sensory neurons
sn.rpkm.df = readr::read_tsv(file.path(root, "sensoryneurons/all_basic_counts.v5.rpkm.txt.gz"))
sn.rpkmavg.df = data.frame(geneID = sn.rpkm.df$gene_id, sensory_neuron=rowMeans(sn.rpkm.df[,-1:-2]))


###########################################################################
# iPSC-derived macrophages
ipsMacrophage.counts.df = read.delim(file.path(root, "macrophage/RNA_count_matrix.txt.gz"))
ipsMacrophage.samplemeta = read.delim(file.path(root, "macrophage/RNA_sample_metadata.txt.gz"))
ipsMacrophage.genemeta = read.delim(file.path(root, "macrophage/RNA_gene_metadata.txt.gz"))
ipsMacrophage.rpkm.df = apply(ipsMacrophage.counts.df, MARGIN=2, FUN=function(x) x / sum(x)) * 1e6 * 1e3
ipsMacrophage.rpkm.df = as.data.frame(apply(ipsMacrophage.rpkm.df, MARGIN=2, FUN=function(x) x / ipsMacrophage.genemeta$length))
ipsMacrophage.rpkm.df$geneID = ipsMacrophage.genemeta$gene_id

macrophage.conditions.df = data.frame(geneID = ipsMacrophage.rpkm.df$geneID)
for (condition in unique(ipsMacrophage.samplemeta$condition_name)) {
  sampleIDs = ipsMacrophage.samplemeta$sample_id[ipsMacrophage.samplemeta$condition_name == condition]
  condition.df = ipsMacrophage.rpkm.df %>% dplyr::select(one_of(sampleIDs))
  macrophage.conditions.df = cbind(macrophage.conditions.df, rowMeans(condition.df))
}
colnames(macrophage.conditions.df)[-1] = paste0("ipsMacrophage.", unique(ipsMacrophage.samplemeta$condition_name))


###########################################################################
# iPSC-derived neurons
ipsneurons.rpkm.df = readr::read_tsv(file.path(root, "ipsneurons/GRCh38/RNA/analysis/ipsneurons.fpkm.txt.gz"))
colnames(ipsneurons.rpkm.df)[-1] = toupper(colnames(ipsneurons.rpkm.df)[-1])
ipsneurons.meta = readr::read_tsv(file.path(root, "ipsneurons/GRCh38/RNA/meta/ipsneurons.rna.metadata.txt"))
ipsneurons.meta = ipsneurons.meta %>% filter(sampleID %in% colnames(ipsneurons.rpkm.df))

ipsneurons.avg.df = data.frame(geneID = ipsneurons.rpkm.df$gene_id)
for (celltype in unique(ipsneurons.meta$celltype2)) {
  sampleIDs = ipsneurons.meta$sampleID[ipsneurons.meta$celltype2 == celltype]
  celltype.df = ipsneurons.rpkm.df %>% dplyr::select(one_of(sampleIDs))
  ipsneurons.avg.df = cbind(ipsneurons.avg.df, rowMeans(celltype.df))
}
colnames(ipsneurons.avg.df)[-1] = paste0("ipsneuron.", unique(ipsneurons.meta$celltype2))


write.table(rpkm.combined.df, file=gzfile(file.path(root, "reference/celltype.rpkm.txt.gz")),
            sep="\t", col.name=T, row.names=F, quote=F)

