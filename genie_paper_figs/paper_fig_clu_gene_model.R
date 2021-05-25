library(tidyverse)
library(wiggleplotr)
library(cowplot)
# library(GenomicRanges)
# library("GenomicFeatures")
library(biomaRt)

setwd("/Users/jeremys/work/opentargets/")

library("EnsDb.Hsapiens.v86")

clu_main_txid = "ENST00000405140"
tx <- transcripts(EnsDb.Hsapiens.v86, filter = GeneNameFilter("CLU"))
plotTranscriptsFromEnsembldb(EnsDb.Hsapiens.v86, gene_names = "CLU", rescale_introns = FALSE)
plotTranscriptsFromEnsembldb(EnsDb.Hsapiens.v86, gene_names = "CLU", rescale_introns = FALSE, transcript_ids = clu_main_txid)

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

bwDir = "/Users/jeremys/work/opentargets/ipsneurons/GRCh38/ATAC/bw"
track_data = data.frame(sample_id = c("Kolf2_iPS_ATAC_1", "Kolf2_iPS_ATAC_2", "Kolf2_iPS_ATAC_3",
                                      "Kolf2_iPSneuron_ATAC_1", "Kolf2_iPSneuron_ATAC_2", "Kolf2_iPSneuron_ATAC_3",
                                      "Kolf2_iNeuron_ATAC_1", "Kolf2_iNeuron_ATAC_2", "Kolf2_iNeuron_ATAC_3", "Kolf2_iNeuron_ATAC_4", "Kolf2_iNeuron_ATAC_5", "Kolf2_iNeuron_ATAC_6"),
                        bigWig = c(file.path(bwDir, "4859STDY7028455.bw"),
                                   file.path(bwDir, "4859STDY7028456.bw"),
                                   file.path(bwDir, "4859STDY7028457.bw"),
                                   file.path(bwDir, "4859STDY7079824.bw"),
                                   file.path(bwDir, "4859STDY7079825.bw"),
                                   file.path(bwDir, "4859STDY7079826.bw"),
                                   file.path(bwDir, "4859STDY7028449.bw"),
                                   file.path(bwDir, "4859STDY7028450.bw"),
                                   file.path(bwDir, "4859STDY7028451.bw"),
                                   file.path(bwDir, "4859STDY7028452.bw"),
                                   file.path(bwDir, "4859STDY7028453.bw"),
                                   file.path(bwDir, "4859STDY7028454.bw")),
                        track_id = c("iPSC", "iPSC", "iPSC", "neuron", "neuron", "neuron", "iNeuron", "iNeuron", "iNeuron", "iNeuron", "iNeuron", "iNeuron"),
                        colour_group = c("iPSC", "iPSC", "iPSC", "neuron", "neuron", "neuron", "iNeuron", "iNeuron", "iNeuron", "iNeuron", "iNeuron", "iNeuron"),
                        scaling_factor = c(1, 258/147, 258/100, 1, 516/488, 400/488, 1, 308/382, 406/382, 406/382, 459/382, 258/382))

gene_tx = transcript_metadata %>%
  dplyr::filter(gene_name == "CLU", transcript_biotype == "protein_coding")
#plotTranscripts(exons[clu_main_txid], cdss[clu_main_txid], transcript_metadata, rescale_introns = FALSE)

gene_exons = exons[clu_main_txid]
gene_cdss = cdss[clu_main_txid]
seqlevelsStyle(gene_exons) <- "UCSC"
seqlevelsStyle(gene_cdss) <- "UCSC"

clu.plots = plotCoverage(gene_exons, gene_cdss,
                         gene_tx, track_data,
                      heights = c(1,1), fill_palette = getGenotypePalette(),
                      transcript_label = T, rescale_introns = F, return_subplots_list = T)

clu.plot.full = cowplot::plot_grid(clu.plots$tx_structure + theme(text = element_text(size=8)),
                                   clu.plots$coverage_plot + theme(text = element_text(size=8)) + facet_grid(track_id~., scales="free"),
                                   ncol = 1, align = "v", rel_heights = c(2, 1.5))

pdf("plots/clu_genemodel_atac.pdf", width = 6.7, height = 2)
clu.plot.full
dev.off()

clu.plot.full.vlines = cowplot::plot_grid(clu.plots$tx_structure + theme(text = element_text(size=8)) +
                                            geom_vline(xintercept = (27598736), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27604964), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27607002), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27607412), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27607795), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27608640), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27608664), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27608798), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27610169), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27610304), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27610986), col = "red", size=0.3, alpha=0.7),
                                          clu.plots$coverage_plot + theme(text = element_text(size=8)) + facet_grid(track_id~., scales="free") +
                                            geom_vline(xintercept = (27598736), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27604964), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27607002), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27607412), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27607795), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27608640), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27608664), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27608798), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27610169), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27610304), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27610986), col = "blue", size=0.3, alpha=0.7),
                                          ncol = 1, align = "v", rel_heights = c(2, 1.5))

pdf("plots/clu_genemodel_atac.lines.pdf", width = 6.7, height = 2.1)
clu.plot.full.vlines
dev.off()


###############################################################################
# Add a coverage track for microglia... but these are in GRCh37 coords
#tbls = supportedUCSCtables(genome="hg19")
txdb_hg19 = makeTxDbFromUCSC( genome="hg19",
                  tablename="knownGene",
                  transcript_ids=clu_main_txid,
                  circ_seqs=DEFAULT_CIRC_SEQS,
                  url="http://genome.ucsc.edu/cgi-bin/",
                  goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
                  taxonomyId=NA,
                  miRBaseBuild=NA)

txdb_hg19_file = "/Users/jeremys/work/opentargets/reference/hsapiens_gene_ensembl.clu.hg19.txdb.rds"
if (file.exists(txdb_hg19_file)) {
  txdb_hg19 = loadDb(txdb_hg19_file)
} else {
  txdb_hg19 = makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL",
                                  transcript_ids=clu_main_txid,
                                  dataset = "hsapiens_gene_ensembl",
                                  host="http://grch37.ensembl.org/")
  saveDb(txdb_hg19, txdb_hg19_file)
}

exons_hg19 = exonsBy(txdb_hg19, by = "tx", use.names = TRUE)
cdss_hg19 = cdsBy(txdb_hg19, by = "tx", use.names = TRUE)
gene_exons = exons_hg19[clu_main_txid]
gene_cdss = cdss_hg19[clu_main_txid]

microglia_bwdir = "/Users/jeremys/work/opentargets/datasets/glass_microglia/ATAC/bw/"
track_data_mg = data.frame(sample_id = c("Kolf2_iPS_ATAC_1", "Kolf2_iPS_ATAC_2", "Kolf2_iPS_ATAC_3"),
                        bigWig = c(file.path(microglia_bwdir, "SRR5955079.bw"),
                                   file.path(microglia_bwdir, "SRR5955080.bw"),
                                   file.path(microglia_bwdir, "SRR5955081.bw"),
                                   file.path(microglia_bwdir, "SRR5955084.bw"),
                                   file.path(microglia_bwdir, "SRR5955085.bw"),
                                   file.path(microglia_bwdir, "SRR5955087.bw"),
                                   file.path(microglia_bwdir, "SRR5955090.bw"),
                                   file.path(microglia_bwdir, "SRR5955091.bw"),
                                   file.path(microglia_bwdir, "SRR5955092.bw")),
                        track_id = c("mg", "mg", "mg", "mg", "mg", "mg", "mg", "mg", "mg"),
                        colour_group = c("mg", "mg", "mg", "mg", "mg", "mg", "mg", "mg", "mg"),
                        scaling_factor = c(1, 58/81, 93/81, 138/81, 119/81, 133/81, 130/81, 47/81, 54/81))

clu.plots.mg = plotCoverage(gene_exons, gene_cdss,
                         gene_tx, track_data_mg,
                         heights = c(1,1), fill_palette = getGenotypePalette()[3],
                         transcript_label = T, rescale_introns = F, return_subplots_list = T)

clu.plot.full.vlines = cowplot::plot_grid(clu.plots.mg$tx_structure + theme(text = element_text(size=8)) +
                                            geom_vline(xintercept = (27456253), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27462481), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27464519), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27464929), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27465312), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466157), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466181), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466315), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27467686), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27467821), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27468503), col = "red", size=0.3, alpha=0.7),
                                          clu.plots.mg$coverage_plot + theme(text = element_text(size=8)) +
                                            geom_vline(xintercept = (27456253), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27462481), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27464519), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27464929), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27465312), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466157), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466181), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466315), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27467686), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27467821), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27468503), col = "blue", size=0.3, alpha=0.7),
                                          ncol = 1, align = "v", rel_heights = c(2,1))

pdf("plots/clu_genemodel_atac.microglia.lines.pdf", width = 6.7, height = 0.9)
clu.plot.full.vlines
dev.off()

#### All atac tracks together
clu.plot.full.vlines = cowplot::plot_grid(clu.plots$tx_structure + theme(text = element_text(size=8)) +
                                            geom_vline(xintercept = (27456253), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27462481), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27464519), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27464929), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27465312), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466157), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466181), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466315), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27467686), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27467821), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27468503), col = "red", size=0.3, alpha=0.7),
                                          clu.plots$coverage_plot + theme(text = element_text(size=8)) + facet_grid(track_id~., scales="free") +
                                            geom_vline(xintercept = (27598736), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27604964), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27607002), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27607412), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27607795), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27608640), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27608664), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27608798), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27610169), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27610304), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27610986), col = "blue", size=0.3, alpha=0.7),
                                          clu.plots.mg$coverage_plot + theme(text = element_text(size=8)) +
                                            geom_vline(xintercept = (27456253), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27462481), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27464519), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27464929), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27465312), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466157), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466181), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466315), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27467686), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27467821), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27468503), col = "blue", size=0.3, alpha=0.7),
                                          ncol = 1, align = "v", rel_heights = c(2,2,0.65))

pdf("plots/clu_genemodel_atac.all.lines.pdf", width = 6.7, height = 2.5)
clu.plot.full.vlines
dev.off()


#### All atac tracks together, except iNeuron
track_data = data.frame(sample_id = c("Kolf2_iPS_ATAC_1", "Kolf2_iPS_ATAC_2", "Kolf2_iPS_ATAC_3",
                                      "Kolf2_iPSneuron_ATAC_1", "Kolf2_iPSneuron_ATAC_2", "Kolf2_iPSneuron_ATAC_3"),
                        bigWig = c(file.path(bwDir, "4859STDY7028455.bw"),
                                   file.path(bwDir, "4859STDY7028456.bw"),
                                   file.path(bwDir, "4859STDY7028457.bw"),
                                   file.path(bwDir, "4859STDY7079824.bw"),
                                   file.path(bwDir, "4859STDY7079825.bw"),
                                   file.path(bwDir, "4859STDY7079826.bw")),
                        track_id = c("iPSC", "iPSC", "iPSC", "neuron", "neuron", "neuron"),
                        colour_group = c("iPSC", "iPSC", "iPSC", "neuron", "neuron", "neuron"),
                        scaling_factor = c(1, 258/147, 258/100, 1, 516/488, 400/488))

gene_exons = exons[clu_main_txid]
gene_cdss = cdss[clu_main_txid]
seqlevelsStyle(gene_exons) <- "UCSC"
seqlevelsStyle(gene_cdss) <- "UCSC"

clu.plots = plotCoverage(gene_exons, gene_cdss,
                         gene_tx, track_data,
                         heights = c(1,1), fill_palette = getGenotypePalette(),
                         transcript_label = T, rescale_introns = F, return_subplots_list = T)

clu.plot.full.vlines = cowplot::plot_grid(clu.plots$tx_structure + theme(text = element_text(size=8)) +
                                            geom_vline(xintercept = (27456253), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27462481), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27464519), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27464929), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27465312), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466157), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466181), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466315), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27467686), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27467821), col = "red", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27468503), col = "red", size=0.3, alpha=0.7),
                                          clu.plots$coverage_plot + theme(text = element_text(size=8)) + facet_grid(track_id~., scales="free") +
                                            geom_vline(xintercept = (27598736), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27604964), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27607002), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27607412), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27607795), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27608640), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27608664), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27608798), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27610169), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27610304), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27610986), col = "blue", size=0.3, alpha=0.7),
                                          clu.plots.mg$coverage_plot + theme(text = element_text(size=8)) +
                                            geom_vline(xintercept = (27456253), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27462481), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27464519), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27464929), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27465312), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466157), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466181), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27466315), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27467686), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27467821), col = "blue", size=0.3, alpha=0.7) +
                                            geom_vline(xintercept = (27468503), col = "blue", size=0.3, alpha=0.7),
                                          ncol = 1, align = "v", rel_heights = c(3,2,1))

pdf("plots/clu_genemodel_atac.ipsc_neuron_mg.lines.pdf", width = 6.7, height = 2)
clu.plot.full.vlines
dev.off()


clu.plot.full.vlines = cowplot::plot_grid(clu.plots$tx_structure + theme(text = element_text(size=8)),
                                          clu.plots$coverage_plot + theme(text = element_text(size=8)) + facet_grid(track_id~., scales="free"),
                                          clu.plots.mg$coverage_plot + theme(text = element_text(size=8)),
                                          ncol = 1, align = "v", rel_heights = c(3,2,1))

pdf("plots/clu_genemodel_atac.ipsc_neuron_mg.pdf", width = 6.7, height = 2)
clu.plot.full.vlines
dev.off()
