#!/usr/bin/env Rscript
library(tidyverse)
library(biomaRt)
library(GenomicFeatures)

root = "/Users/jeremys/work/opentargets"
root = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys"

txdb = makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL", 
                           dataset = "hsapiens_gene_ensembl", 
                           host="dec2017.archive.ensembl.org")

txdb_file = file.path(root, "reference/GRCh38/ensembldb.91.rds")
saveDb(txdb, file = txdb_file)


#txdb = loadDb(txdb_file)
#cdss = cdsBy(txdb, by = "tx", use.names = TRUE)
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
rtracklayer::export.bed(reduce(unlist(exons)), paste0(root, "reference/GRCh38/ensembldb.91.exons.bed"))


ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", host = "dec2017.archive.ensembl.org")
ensembl_dataset = useDataset("hsapiens_gene_ensembl",mart=ensembl_mart)

selected_attributes = c("ensembl_transcript_id", "ensembl_gene_id", 
                        "external_gene_name", "strand", 
                        "gene_biotype", "transcript_biotype")
txmeta = getBM(attributes = selected_attributes, mart = ensembl_dataset)
txmeta = dplyr::rename(txmeta,
                       transcript_id = ensembl_transcript_id, 
                       gene_id = ensembl_gene_id, 
                       gene_name = external_gene_name)
head(txmeta)

txmeta_file = file.path(root, "reference/GRCh38/ensembldb.91.transcripts.rds")
saveRDS(txmeta, txmeta_file)



selected_transcripts = txmeta %>%
  dplyr::filter(gene_name == "NCOA7", transcript_biotype == "protein_coding")
tx_ids = selected_transcripts$transcript_id


###########################################################################
# GRCh37
txdb = makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL", 
                           dataset = "hsapiens_gene_ensembl", 
                           host="grch37.ensembl.org")

txdb_file = file.path(root, "reference/GRCh37/ensembldb.91.grch37.rds")
saveDb(txdb, file = txdb_file)
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
rtracklayer::export.bed(reduce(unlist(exons)), file.path(root, "reference/GRCh37/ensembldb.91.grch37.exons.bed"))

ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org")
ensembl_dataset = useDataset("hsapiens_gene_ensembl",mart=ensembl_mart)
selected_attributes = c("ensembl_transcript_id", "ensembl_gene_id", 
                        "external_gene_name", "strand", 
                        "gene_biotype", "transcript_biotype")
txmeta = getBM(attributes = selected_attributes, mart = ensembl_dataset)
txmeta = dplyr::rename(txmeta,
                       transcript_id = ensembl_transcript_id, 
                       gene_id = ensembl_gene_id, 
                       gene_name = external_gene_name)
txmeta_file = file.path(root, "reference/GRCh37/ensembldb.91.grch37.transcripts.rds")
saveRDS(txmeta, txmeta_file)

