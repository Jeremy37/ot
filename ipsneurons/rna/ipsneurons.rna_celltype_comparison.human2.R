#!/usr/bin/env Rscript
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(DESeq2)
options(stringsAsFactors = F)

#root = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys"
root = "/Users/jeremys/work/opentargets"

#source("~/js29/src/R/pca.projection.R")
source(file.path(root, "src/R/pca.projection.R"))
source(file.path(root, "src/R/countUtils.R"))

outputPath = file.path(root, "ipsneurons/GRCh38/RNA/analysis")
gtex.rpkmfile = file.path("/Users/jeremys/work/sensoryneurons/snqtl/results/MDS", "rpkms.gtex.ips.drg.sn.rds")
gtex.metafile = file.path("/Users/jeremys/work/sensoryneurons/snqtl/results/MDS", "rpkms.sample_meta.txt")
counts.fname = file.path(outputPath, "ipsneurons.counts.txt.gz")
txRefPath = file.path(root, "reference/GRCh38/Homo_sapiens.GRCh38.89.transcript_data.rds")

sample.meta.full = readtsv(file.path(root, "ipsneurons/GRCh38/RNA/meta/ipsneurons.rna.metadata.txt"))
rownames(sample.meta.full) = sample.meta.full$sampleID

counts = readCounts(counts.fname)
counts.sampleIDs = colnames(counts)[3:15]

sample.meta = sample.meta.full[counts.sampleIDs,]
sample.meta = data.frame(ID=sample.meta$Sample_name, group=paste0("IPS_", sample.meta$celltype))

# Change the column names to friendly sample names (e.g. BOB_iN_D11 rather than 4860STDY7028458)
colnames(counts)[3:15] = sample.meta.full[counts.sampleIDs,]$Sample_name

sample.mat = countsToRPKM(counts)
#rpkm.df = cbind(data.frame(gene_id = rownames(rpkm)), rpkm)

# Let's add in some of our other iPSC-derived samples
inhouse.df = read.delim(file.path(root, "reference/IPS/inhouse.combined.rpkm.txt.gz"))
inhouse.meta = read.delim(file.path(root, "reference/IPS/inhouse.combined.meta.txt"))
toselect = c("HPSI0114i-bezi_1", "HPSI0114i-kolf_2", "nukw_1a", "auim_2", "T31996", "T32042", "xugn_A", "vorx_B")
inhouse.sel = inhouse.df %>% dplyr::select(gene_id, which(colnames(inhouse.df) %in% toselect))
#inhouse.sel.meta = inhouse.meta %>% dplyr::filter(sampleID %in% toselect) %>%
#  dplyr::rename(ID=sampleID)
inhouse.sel.meta = data.frame(ID=c("IPSC-bezi_1", "IPSC-kolf_2", "SN_nukw_1a", "SN_auim_2", "DRG_T31996", "DRG_T32042", "Macro_xugn_naive", "Macro_vorx_IFNg"),
                              group=c("IPSC", "IPSC", "Sensory neuron", "Sensory neuron", "DRG", "DRG", "Macrophage", "Macrophage"))
colnames(inhouse.sel)[2:9] = inhouse.sel.meta$ID

rpkm.df = as.data.frame(sample.mat)
rpkm.df$gene_id = rownames(rpkm.df)
sample.combined.mat = rpkm.df %>% dplyr::inner_join(inhouse.sel, by="gene_id")
rownames(sample.combined.mat) = sample.combined.mat$gene_id
sample.combined.mat = sample.combined.mat %>% dplyr::select(-gene_id) %>% as.matrix()
sample.combined.meta = rbind(sample.meta, inhouse.sel.meta)

barres.df = readtsv(file.path(root, "reference/Barres_RNAseq/barres_rnaseq.human.txt"))
# Remove duplicate gene IDs. First order by expression level so that if one version of the
# gene is expressed we capture that one
barres.df = barres.df[order(-rowMeans(barres.df[,-1])),]
barres.df = barres.df[!duplicated(barres.df$Gene),]
barres.mat = as.matrix(barres.df[,-1])
rownames(barres.mat) = barres.df$Gene
barres.mat = HGNCToEnsemblRownames(barres.mat, txRefPath)

barres.meta.full = readtsv(file.path(root, "reference/Barres_RNAseq/barres_rnaseq.human.metadata.txt"))
rownames(barres.meta.full) = barres.meta.full$sampleID

# Make the metadata DFs expected by the plotting functions (ID & group columns)
# Make sure the metadata is in the same order as the cols of the matrix
barres.meta.full = barres.meta.full[colnames(barres.mat),]
barres.meta = data.frame(ID=barres.meta.full$sampleID, group=paste0("BAR_", barres.meta.full$celltype))

pdf(file.path(outputPath, "ipsneurons.rna_comparison.barres.allgenes.pdf"), width=10, height=8)
heatmapWithReference(barres.mat, barres.meta, sample.mat, sample.meta, averageGroups = T, title="Barres lab brain - gene expression correlations")

# Show all samples, but remove a few of the numerous mature astrocytes
barres.meta.reduced = barres.meta %>% dplyr::filter(!grepl('mature_astrocyte_([4-9]|10|11|12)|tumor_astrocyte|fetal_astrocyte_([4-6])|hippocampi_astrocyte_4', ID))
barres.mat.reduced = barres.mat[, barres.meta.reduced$ID]
heatmapWithReference(barres.mat.reduced, barres.meta.reduced, sample.mat, sample.meta, averageGroups = F, title="Barres lab brain - gene expression correlations\nindividual samples")

heatmapWithReference(barres.mat, barres.meta, sample.combined.mat, sample.combined.meta, averageGroups = T, title="Barres lab brain - gene expression correlations\nincluding iPSC, DRG, SN, Macrophage")

doPcaProjection(barres.mat, barres.meta, sample.combined.mat, sample.combined.meta, maxPC=3)
dev.off()

# Do the heatmap and PCA plots based on a subset of genes only. First we need to
# get the ensembl IDs for the genes.
neuronal.gene.df = readtsv(file.path(outputPath, "neuronal_genelist.txt"))

txdb = readRDS(txRefPath)
genedb = txdb %>% dplyr::select(ensembl_gene_id, external_gene_name, gene_biotype) %>% 
  unique() %>% dplyr::rename("gene_id" = "ensembl_gene_id", "gene_name" = "external_gene_name")

neuronal.gene.df %<>% dplyr::inner_join(genedb %>% dplyr::select(gene_id, gene_name), by="gene_name")

barres.mat.neuronal = barres.mat[neuronal.gene.df$gene_id[neuronal.gene.df$gene_id %in% rownames(barres.mat)],]
sample.mat.neuronal = sample.mat[neuronal.gene.df$gene_id[neuronal.gene.df$gene_id %in% rownames(sample.mat)],]

pdf(file.path(outputPath, "ipsneurons.rna_comparison.barres.neuronal_genes.pdf"), width=10, height=8)
heatmapWithReference(barres.mat.neuronal, barres.meta, sample.mat.neuronal, sample.meta, averageGroups = T, title="Barres lab brain correlations - ~160 neuronal genes")

barres.mat.neuronal.reduced = barres.mat.neuronal[, barres.meta.reduced$ID]
heatmapWithReference(barres.mat.neuronal.reduced, barres.meta.reduced, sample.mat.neuronal, sample.meta, averageGroups = F, title="Barres lab brain - ~160 neuronal genes\nindividual samples")

sample.combined.mat.neuronal = sample.combined.mat[neuronal.gene.df$gene_id[neuronal.gene.df$gene_id %in% rownames(sample.mat)],]
heatmapWithReference(barres.mat.neuronal, barres.meta, sample.combined.mat.neuronal, sample.combined.meta, averageGroups = T, title="Barres lab brain - ~160 neuronal genes\nincluding iPSC, DRG, SN, Macrophage")

doPcaProjection(barres.mat.neuronal, barres.meta, sample.combined.mat.neuronal, sample.combined.meta, maxPC=3)
dev.off()

