#!/usr/bin/env Rscript
library(dplyr)
source("../R/countUtils.R")
root = "/Users/jeremys/work/opentargets"

###########################################################################
# In-house samples

# The GTEx RPKM file I produced before included HIPSCI iPSCs,
# sensory neurons, and DRG
#gtex.rpkmfile = file.path("/lustre/scratch117/cellgen/team170/js29/sensoryneurons/snqtl/results/MDS", "rpkms.gtex.ips.drg.sn.rds")
gtex.rpkmfile = file.path("/Users/jeremys/work/sensoryneurons/snqtl/results/MDS", "rpkms.gtex.ips.drg.sn.rds")
gtex.metafile = file.path("/Users/jeremys/work/sensoryneurons/snqtl/results/MDS", "rpkms.sample_meta.txt")

gtex.df = readRDS(gtex.rpkmfile) %>% dplyr::rename(gene_id = Name)
gtex.meta.all = read.delim(gtex.metafile, stringsAsFactors = F)
gtex.meta = gtex.meta.all %>%
  dplyr::filter(SAMPID %in% colnames(gtex.df)) %>%
  dplyr::rename(sampleID = SAMPID)

combined.df = gtex.df[, !grepl("GTEX", colnames(gtex.df))]
sample.meta = gtex.meta %>%
  dplyr::filter(sampleID %in% colnames(combined.df)) %>%
  dplyr::select(sampleID, SMTSD) %>%
  dplyr::rename(group = SMTSD)
sample.meta$group[grepl("IPSDSN", sample.meta$group)] = "IPSDSN"

# While we're at it, get RPKM averages for GTEx tissues
gtex.smtsAvg.df = data.frame(gene_id = gtex.df$gene_id)
gtex.smtsd = unique(gtex.meta$SMTSD[gtex.meta$sampleID %in% colnames(gtex.df)[grepl("GTEX", colnames(gtex.df))]])
for (tissue in gtex.smtsd) {
  sampleIDs = gtex.meta$sampleID[gtex.meta$SMTSD == tissue]
  tissue.df = gtex.df %>% dplyr::select(one_of(sampleIDs))
  gtex.smtsAvg.df = cbind(gtex.smtsAvg.df, rowMeans(tissue.df))
}
colnames(gtex.smtsAvg.df)[-1] = paste0("gtex.", gtex.smtsd)

gzf = gzfile(file.path(root, "reference/tissueRPKM/gtex.rpkm_averages.txt.gz"), "w")
write.table(gtex.smtsAvg.df, gzf, row.names=F, col.names=T, quote=F, sep="\t")
close(gzf)

###########################################################################
# Kaur's macrophage RPKMs
macrophage.counts = read.delim("/Users/jeremys/work/opentargets/macrophage/RNA_count_matrix.txt.gz")
macrophage.genes = read.delim("/Users/jeremys/work/opentargets/macrophage/RNA_gene_metadata.txt.gz")
macrophage.samples = read.delim("/Users/jeremys/work/opentargets/macrophage/RNA_sample_metadata.txt.gz")

source("/Users/jeremys/work/opentargets/src/R/countUtils.R")
macrophage.rpkm = as.data.frame(countsToRPKM(macrophage.counts, macrophage.genes, useDeseq=T))
macrophage.rpkm$gene_id = macrophage.genes$gene_id
macrophage.meta = data.frame(sampleID = macrophage.samples$sample_id,
                             group = paste0("macrophage_", macrophage.samples$condition_name))

# Add the macrophage RPKMs to our other samples
combined.meta = rbind(sample.meta, macrophage.meta)
combined.df = combined.df %>%
  dplyr::full_join(macrophage.rpkm, by="gene_id")
# 32,431 genes overlap

gzf = gzfile(file.path(root, "reference/tissueRPKM/inhouse.combined.rpkm.txt.gz"), "w")
write.table(combined.df, file=gzf, row.names=F, col.names=T, quote=F, sep="\t")
close(gzf)

write.table(combined.meta, file=file.path(root, "reference/inhouse.combined.meta.txt"), row.names=F, col.names=T, quote=F, sep="\t")


###########################################################################
# Barres lab primary brain sorted cell types
barres.meta = readr::read_tsv(file.path(root, "reference/Barres_RNAseq/barres_rnaseq.human.metadata.txt")) %>%
  dplyr::rename(group = celltype)
barres.df = readr::read_tsv(file.path(root, "reference/Barres_RNAseq/barres_rnaseq.human.txt")) %>%
  dplyr::rename(symbol = Gene) %>%
  filter(!duplicated(symbol))
barres.df$symbol = toupper(barres.df$symbol)

# Map Ensembl gene IDs to gene symbols
hgnc.map.df = readr::read_tsv(file.path(root, "reference/hgnc.ensembl.map.txt")) %>%
  dplyr::rename(gene_id = ensembl_gene_id)
hgnc.map.df$symbol = toupper(hgnc.map.df$symbol)

# Join the barres data with ensembl gene IDs
barres.joined.df = barres.df %>% dplyr::left_join(hgnc.map.df %>% dplyr::select(symbol, symbol_type, gene_id), by="symbol") %>%
  dplyr::filter(!duplicated(symbol)) %>% dplyr::select(symbol, symbol_type, gene_id, everything())
barres.joined.df = barres.joined.df %>%
  dplyr::arrange(-rowMeans(barres.joined.df %>% dplyr::select(-symbol, -symbol_type, -gene_id))) %>%
  dplyr::filter(!duplicated(gene_id))
# Now we have no duplicated symbols or ensembl IDs

barres.meta$group = paste0("Barres.", barres.meta$group)
combined.meta = rbind(combined.meta, barres.meta %>% dplyr::select(sampleID, group))

combined.df = combined.df %>%
  dplyr::full_join(barres.joined.df, by="gene_id") %>%
  dplyr::rename(barres_symbol = symbol, barres_symbol_type = symbol_type)
#View(combined.df[,seq(1,760,10)])


###########################################################################
# ipsNeurons iPSC-derived cells
ipsneurons.rpkm.df = readr::read_tsv(file.path(root, "ipsneurons/GRCh38/RNA/analysis/ipsneurons.fpkm.txt.gz"))
colnames(ipsneurons.rpkm.df)[-1] = toupper(colnames(ipsneurons.rpkm.df)[-1])
ipsneurons.meta = readr::read_tsv(file.path(root, "ipsneurons/GRCh38/RNA/meta/ipsneurons.rna.metadata.txt"))
ipsneurons.meta = ipsneurons.meta %>% filter(sampleID %in% colnames(ipsneurons.rpkm.df)) %>%
  dplyr::rename(group = celltype2)

combined.df = combined.df %>%
  dplyr::full_join(ipsneurons.rpkm.df, by="gene_id")
combined.meta = rbind(combined.meta, ipsneurons.meta %>% dplyr::select(sampleID, group))


###########################################################################
# Fiona / Gaffney microglia
microglia.counts.df = readr::read_tsv(file.path(root, "microglia", "microglia_counts_matrix.gz"))
transcript_data.df = readr::read_tsv(file.path(root, "reference", "GRCh38", "GRCh38.p12.geneid_hgnc_length.txt.gz"))

gene_data.df = transcript_data.df %>% group_by(gene_id) %>%
  arrange(gene_id, -transcript_length) %>%
  dplyr::summarise(length = max(transcript_length), transcript_id = first(transcript_id), hgnc_symbol = first(hgnc_symbol))
microglia.genes.df = microglia.counts.df %>% dplyr::select(gene_id = Geneid) %>%
  dplyr::left_join(gene_data.df, by="gene_id")

microglia.counts.mat = microglia.counts.df %>% dplyr::select(-Geneid) %>% as.matrix()
rownames(microglia.counts.mat) = microglia.counts.df$Geneid
microglia.rpkm.df = as.data.frame(countsToRPKM(microglia.counts.mat, microglia.genes.df, useDeseq=T))

microglia.meta = data.frame(sampleID = colnames(microglia.rpkm.df),
                            group = paste0("microglia_gaffney"))
microglia.rpkm.df$gene_id = microglia.genes.df$gene_id

combined.df = combined.df %>%
  dplyr::full_join(microglia.rpkm.df, by="gene_id")
combined.meta = rbind(combined.meta, microglia.meta)


###########################################################################
# Blueprint



###########################################################################

# Now annotate the full combined table with HGNC gene symbols (and their previous symbols)
hgnc.df = readr::read_tsv(file.path(root, "reference/hgnc_complete_set.txt")) %>%
  dplyr::select(ensembl_gene_id, symbol, alias_symbol, prev_symbol) %>%
  dplyr::rename(gene_id = ensembl_gene_id) %>%
  dplyr::filter(!is.na(gene_id), !duplicated(gene_id))

# We want 1 column containing symbol aliases and previous symbols
mergeSymbols = function(x) {
  if (is.na(x[1]) & is.na(x[2])) { "" }
  else if (is.na(x[1])) { x[2] }
  else if (is.na(x[2])) { x[1] }
  else paste(x[1], x[2], sep = ";")
}
hgnc.df$other_symbols = apply(hgnc.df %>% dplyr::select(alias_symbol, prev_symbol),
                              MARGIN = 1, FUN = mergeSymbols)
#View(hgnc.df %>% dplyr::select(gene_id, symbol, alias_symbol, prev_symbol, other_symbols))
sum(duplicated(hgnc.df$gene_id))

combined.df = combined.df %>% dplyr::left_join(hgnc.df %>% dplyr::select(symbol, other_symbols, gene_id), by="gene_id") %>%
  dplyr::select(gene_id, symbol, other_symbols, barres_symbol, barres_symbol_type, everything())

# Write out the full RPKM table of all samples
gzf = gzfile(file.path(root, "reference/tissueRPKM/tissues.combined.rpkm.txt.gz"), "w")
write.table(combined.df, file=gzf, row.names=F, col.names=T, quote=F, sep="\t")
close(gzf)

write.table(combined.meta, file=file.path(root, "reference/tissueRPKM/tissues.combined.meta.txt"),
            row.names=F, col.names=T, quote=F, sep="\t")


###########################################################################
# Compute averages for different tissue/cell type groups

groups = unique(combined.meta$group)
group.avg.df = combined.df %>% dplyr::select(gene_id, symbol, other_symbols, barres_symbol, barres_symbol_type)
group.sd.df = combined.df %>% dplyr::select(gene_id, symbol, other_symbols, barres_symbol, barres_symbol_type)
for (grp in groups) {
  sampleIDs = combined.meta$sampleID[combined.meta$group == grp]
  group.df = combined.df %>% dplyr::select(one_of(sampleIDs))
  group.avg.df = cbind(group.avg.df, rowMeans(as.matrix(group.df)))
  group.sd.df = cbind(group.sd.df, rowSds(as.matrix(group.df)))
}
colnames(group.avg.df)[-1:-5] = unique(combined.meta$group)
colnames(group.sd.df)[-1:-5] = unique(combined.meta$group)
# Get the coefficient of variation for each tissue/gene: sd / mean
group.cv.df = group.sd.df
group.cv.df[,-1:-5] = group.sd.df[,-1:-5] / group.avg.df[,-1:-5]

combined.avg.df = group.avg.df %>% dplyr::full_join(gtex.smtsAvg.df, by="gene_id")

gzf = gzfile(file.path(root, "reference/tissueRPKM/tissues.combined.rpkm_average.txt.gz"), "w")
write.table(combined.avg.df %>% dplyr::select(-barres_symbol_type),
            file=gzf, row.names=F, col.names=T, quote=F, sep="\t")
close(gzf)

gzf = gzfile(file.path(root, "reference/tissueRPKM/tissues.combined.rpkm_sd.txt.gz"), "w")
write.table(group.sd.df %>% dplyr::select(-barres_symbol_type),
            file=gzf, row.names=F, col.names=T, quote=F, sep="\t")
close(gzf)

gzf = gzfile(file.path(root, "reference/tissueRPKM/tissues.combined.rpkm_cv.txt.gz"), "w")
write.table(group.cv.df %>% dplyr::select(-barres_symbol_type),
            file=gzf, row.names=F, col.names=T, quote=F, sep="\t")
close(gzf)

# Save the same files, but with values rounded to 3 decimal places
roundToString = function(x) {
  if (is.na(x)) { return("") }
  sprintf("%.3g", round(x, digits = 3))
}
combined.avg.df.short = cbind(combined.avg.df[,1:4], apply(combined.avg.df[,-1:-5], MARGIN=c(1, 2), FUN=roundToString))
group.sd.df.short = cbind(group.sd.df[,1:4], apply(group.sd.df[,-1:-5], MARGIN=c(1, 2), FUN=roundToString))
group.cv.df.short = cbind(group.cv.df[,1:4], apply(group.cv.df[,-1:-5], MARGIN=c(1, 2), FUN=roundToString))

write.table(combined.avg.df.short,
            file=file.path(root, "reference/tissueRPKM/tissues.combined.rpkm_avg.rounded.txt"),
            row.names=F, col.names=T, quote=F, sep="\t")

write.table(group.sd.df.short,
            file=file.path(root, "reference/tissueRPKM/tissues.combined.rpkm_sd.rounded.txt"),
            row.names=F, col.names=T, quote=F, sep="\t")

write.table(group.cv.df.short,
            file=file.path(root, "reference/tissueRPKM/tissues.combined.rpkm_cv.rounded.txt"),
            row.names=F, col.names=T, quote=F, sep="\t")


###########################################################################
# Write a smaller subset of cells/tissues of greatest interest

gzf = gzfile(file.path(root, "reference/tissueRPKM/tissues.selected.rpkm_average.txt.gz"), "w")
combined.avg.df.selected = combined.avg.df.short %>%
  dplyr::select(gene_id, IPSC, iNeuron_d9, iNeuron_d11, NPC, neuron, IPSDSN,
                macrophage_naive, macrophage_IFNg, macrophage_SL1344, macrophage_IFNg_SL1344, microglia_gaffney,
                `gtex.Brain - Hippocampus`, `gtex.Brain - Frontal Cortex (BA9)`)
write.table(combined.avg.df.selected, file=gzf, row.names=F, col.names=T, quote=F, sep="\t")
close(gzf)

gzf = gzfile(file.path(root, "reference/tissueRPKM/tissues.selected.expr_ranks.txt.gz"), "w")
combined.avg.df.selected.ranks = cbind(combined.avg.df[,1:4],
                                       apply(-combined.avg.df[,-1:-5], MARGIN=2,
                                             FUN=function(x) rank(x, na.last = "keep", ties.method = "average")))
write.table(combined.avg.df.selected.ranks, file=gzf, row.names=F, col.names=T, quote=F, sep="\t", na="")
close(gzf)



