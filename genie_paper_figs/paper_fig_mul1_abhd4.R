library(tidyverse)
library(wiggleplotr)
library(cowplot)
# library(GenomicRanges)
library("GenomicFeatures")
library(biomaRt)

setwd("/Users/jeremys/work/opentargets/")

library("EnsDb.Hsapiens.v86")
tx <- transcripts(EnsDb.Hsapiens.v86, filter = GeneNameFilter("MUL1"))
plotTranscriptsFromEnsembldb(EnsDb.Hsapiens.v86, gene_names = "MUL1", 
                             transcript_ids = c("ENST00000264198"))

ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", host = "sep2019.archive.ensembl.org")
ensembl_dataset = useDataset("hsapiens_gene_ensembl", mart=ensembl_mart)
selected_attributes = c("ensembl_transcript_id", "ensembl_gene_id", 
                        "external_gene_name", "strand", 
                        "gene_biotype", "transcript_biotype")
data = getBM(attributes = selected_attributes, mart = ensembl_dataset)
head(data)
transcript_metadata = dplyr::rename(data, 
                                    transcript_id = ensembl_transcript_id, 
                                    gene_id = ensembl_gene_id, 
                                    gene_name = external_gene_name)


txdb_file = "/Users/jeremys/work/opentargets/reference/hsapiens_gene_ensembl.txdb.rds"
if (file.exists(txdb_file)) {
  txdb = loadDb(txdb_file)
} else {
  txdb = makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL",
                             dataset = "hsapiens_gene_ensembl",
                             host="sep2019.archive.ensembl.org")
  saveDb(txdb, txdb_file)
}

exons = exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)

############## MUL1
bwDir = "/Users/jeremys/work/opentargets/ipsneurons/GRCh38/ATAC/bw"
track_data = data.frame(sample_id = c("Kolf2_iPS_ATAC_1", "Kolf2_iPS_ATAC_2", "Kolf2_iPS_ATAC_3"),
                        bigWig = c(file.path(bwDir, "4859STDY7028454.bw"),
                                   file.path(bwDir, "4859STDY7028455.bw"),
                                   file.path(bwDir, "4859STDY7028456.bw")),
                        track_id = c("iPSC", "iPSC", "iPSC"),
                        colour_group = c("iPSC", "iPSC", "iPSC"),
                        scaling_factor = c(1, 258/147, 258/100))

gene_tx = transcript_metadata %>%
  dplyr::filter(gene_name == "MUL1", transcript_biotype == "protein_coding")
gene_txid = gene_tx$transcript_id
#plotTranscripts(exons[gene_txid], cdss[gene_txid], transcript_metadata, rescale_introns = TRUE)

gene_exons = exons[gene_txid]
gene_cdss = cdss[gene_txid]
seqlevelsStyle(gene_exons) <- "UCSC"
seqlevelsStyle(gene_cdss) <- "UCSC"

mul1.plots = plotCoverage(gene_exons, gene_cdss,
                      transcript_metadata, track_data,
                      heights = c(1,1), fill_palette = getGenotypePalette(),
                      transcript_label = T, rescale_introns = F, return_subplots_list = T)
mul1.plot.full = cowplot::plot_grid(mul1.plots$coverage_plot + geom_vline(xintercept = (20508117), col = "blue"),
                   mul1.plots$tx_structure + geom_vline(xintercept = (20508117), col = "blue"),
                   ncol = 1, align = "v", rel_heights = c(1,1))

mul1.plots = plotCoverage(gene_exons, gene_cdss,
                          transcript_metadata, track_data,
                          heights = c(1,1), fill_palette = getGenotypePalette(),
                          transcript_label = F, rescale_introns = F, region_coords = c(20506500, 20508500), return_subplots_list = T)
mul1.plot.zoom = cowplot::plot_grid(mul1.plots$coverage_plot + geom_vline(xintercept = (20508117), col = "blue"),
                                    mul1.plots$tx_structure + geom_vline(xintercept = (20508117), col = "blue"),
                                    ncol = 1, align = "v", rel_heights = c(1,1))


############## ABHD4
gene_tx = transcript_metadata %>%
  dplyr::filter(gene_name == "ABHD4", transcript_biotype == "protein_coding")
gene_txid = "ENST00000428304"
gene_exons = exons[gene_txid]
gene_cdss = cdss[gene_txid]
seqlevelsStyle(gene_exons) <- "UCSC"
seqlevelsStyle(gene_cdss) <- "UCSC"

abhd4.plots = plotCoverage(gene_exons, gene_cdss,
                      transcript_metadata, track_data,
                      heights = c(1,1), fill_palette = getGenotypePalette(),
                      transcript_label = T, rescale_introns = F, return_subplots_list = T)
abhd4.plot.full = cowplot::plot_grid(abhd4.plots$coverage_plot + geom_vline(xintercept = (22599022), col = "blue"),
                                     abhd4.plots$tx_structure + geom_vline(xintercept = (22599022), col = "blue"),
                                     ncol = 1, align = "v", rel_heights = c(1,1))

abhd4.plots = plotCoverage(gene_exons, gene_cdss,
                            transcript_metadata, track_data,
                            heights = c(1,1), fill_palette = getGenotypePalette(),
                            transcript_label = F, rescale_introns = F, region_coords = c(22598000, 22599500), return_subplots_list = T)
abhd4.plot.zoom = cowplot::plot_grid(abhd4.plots$coverage_plot + geom_vline(xintercept = (22599022), col = "blue"),
                                     abhd4.plots$tx_structure + geom_vline(xintercept = (22599022), col = "blue"),
                                     ncol = 1, align = "v", rel_heights = c(1,1))


pdf("plots/mul1_abhd4.pdf", width = 6, height = 1.8)
mul1.plot.full
mul1.plot.zoom
abhd4.plot.full
abhd4.plot.zoom
dev.off()


# There's no reason to plot ATAC coverage for these slicing QTLs, and indeed there's nothing interesting when we do
############## TAF1C
gene_tx = transcript_metadata %>%
  dplyr::filter(gene_name == "TAF1C", transcript_biotype == "protein_coding")
gene_txid = gene_tx$transcript_id
gene_exons = exons[gene_tx$transcript_id]
gene_cdss = cdss[gene_tx$transcript_id]
seqlevelsStyle(gene_exons) <- "UCSC"
seqlevelsStyle(gene_cdss) <- "UCSC"

p.taf1c.full = plotCoverage(gene_exons, gene_cdss,
                            transcript_metadata, track_data,
                            heights = c(1,2), fill_palette = getGenotypePalette(),
                            transcript_label = T, rescale_introns = F)

gene_txid = "ENST00000567759"
gene_exons = exons[gene_txid]
gene_cdss = cdss[gene_txid]
seqlevelsStyle(gene_exons) <- "UCSC"
seqlevelsStyle(gene_cdss) <- "UCSC"
p.taf1c.zoom = plotCoverage(gene_exons, gene_cdss,
                             transcript_metadata, track_data,
                             heights = c(2,1), fill_palette = getGenotypePalette(),
                             transcript_label = F, rescale_introns = F, region_coords = c(84184500, 84185500))


############## DRAM2
gene_tx = transcript_metadata %>%
  dplyr::filter(gene_name == "DRAM2", transcript_biotype == "protein_coding")
gene_txid = gene_tx$transcript_id
gene_exons = exons[gene_tx$transcript_id]
gene_cdss = cdss[gene_tx$transcript_id]
seqlevelsStyle(gene_exons) <- "UCSC"
seqlevelsStyle(gene_cdss) <- "UCSC"

p.dram2.full = plotCoverage(gene_exons, gene_cdss,
                            transcript_metadata, track_data,
                            heights = c(1,1), fill_palette = getGenotypePalette(),
                            transcript_label = T, rescale_introns = F)

gene_txid = "ENST00000286692"
gene_exons = exons[gene_txid]
gene_cdss = cdss[gene_txid]
seqlevelsStyle(gene_exons) <- "UCSC"
seqlevelsStyle(gene_cdss) <- "UCSC"
p.dram2.zoom = plotCoverage(gene_exons, gene_cdss,
                            transcript_metadata, track_data,
                            heights = c(2,1), fill_palette = getGenotypePalette(),
                            transcript_label = F, rescale_introns = F, region_coords = c(111137000, 111140200))


############## SDF4
gene_tx = transcript_metadata %>%
  dplyr::filter(gene_name == "SDF4", transcript_biotype == "protein_coding")
gene_txid = gene_tx$transcript_id
gene_exons = exons[gene_tx$transcript_id]
gene_cdss = cdss[gene_tx$transcript_id]
seqlevelsStyle(gene_exons) <- "UCSC"
seqlevelsStyle(gene_cdss) <- "UCSC"

p.sdf4.full = plotCoverage(gene_exons, gene_cdss,
                            transcript_metadata, track_data,
                            heights = c(1,1), fill_palette = getGenotypePalette(),
                            transcript_label = T, rescale_introns = F)

gene_txid = "ENST00000660930"
gene_exons = exons[gene_txid]
gene_cdss = cdss[gene_txid]
seqlevelsStyle(gene_exons) <- "UCSC"
seqlevelsStyle(gene_cdss) <- "UCSC"
p.sdf4.zoom = plotCoverage(gene_exons, gene_cdss,
                            transcript_metadata, track_data,
                            heights = c(2,1), fill_palette = getGenotypePalette(),
                            transcript_label = F, rescale_introns = F, region_coords = c(1228400, 1232500))







