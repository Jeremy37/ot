library("GenomicFeatures")
library("biomaRt")
library("dplyr")

ensemblHost = "may2017.archive.ensembl.org"
#Make TranscriptDb object from biomart
txdb89 = makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", 
                             host=ensemblHost)
saveDb(txdb89, "../../../reference/GRCh38/TranscriptDb_GRCh38.89.db")
saveDb(txdb89, "TranscriptDb_GRCh38.89.db")

#Downlaod transcript metadata from Ensembl
ensembl = useMart("ENSEMBL_MART_ENSEMBL", host = ensemblHost)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes = listAttributes(ensembl)

#Define attributes to be downloaded from biomart
biomart_attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype","status", 
                       "chromosome_name","strand", "transcript_start","transcript_end","ensembl_transcript_id", 
                       "transcript_status", "transcript_tsl","transcript_version", "transcript_appris",
                       "transcript_gencode_basic","external_transcript_name", 
                       "transcript_length", "transcript_biotype", "ccds")
refseq_attributes = c("ensembl_gene_id","ensembl_transcript_id","refseq_mrna")

# Download data for sample genes
genes = c("ENSG00000111912", "ENSG00000266094","ENSG00000068028")
data = getBM(attributes = biomart_attributes, 
             filters = c("ensembl_gene_id"), 
             values = genes, 
             mart = ensembl)
data

data1 = getBM(attributes = refseq_attributes, 
              filters = c("ensembl_gene_id"), 
              values = genes, 
              mart = ensembl)
data1

# Download data for all transcripts
transcript_data = getBM(attributes = biomart_attributes, mart = ensembl)
saveRDS(transcript_data, "../../../reference/GRCh38/Homo_sapiens.GRCh38.89.transcript_data.rds")
saveRDS(transcript_data, "Homo_sapiens.GRCh38.89.transcript_data.rds")
# refseq_data = getBM(attributes = refseq_attributes, mart = ensembl)
# saveRDS(refseq_data, "../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.89.refseq_data.rds")

# Download data for gene exons
#exons.all = exonsBy(txdb89, by = "tx", use.names = TRUE)
exon_attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype", "exon_chrom_start", "exon_chrom_end")
exon_data = getBM(attributes = exon_attributes, mart = ensembl)
exon_data = exon_data %>% dplyr::left_join(transcript_data %>% dplyr::select(ensembl_gene_id, chromosome_name) %>% unique(), by="ensembl_gene_id")
write.table(exon_data, "../../../reference/GRCh38/Homo_sapiens.GRCh38.89.exons.all.txt", sep = "\t", row.names=F, quote=F)

exon_data.pc = exon_data %>%
  dplyr::filter(gene_biotype == "protein_coding") %>%
  dplyr::rename(gene_id = ensembl_gene_id, gene_name = external_gene_name)
write.table(exon_data.pc, "../../../reference/GRCh38/Homo_sapiens.GRCh38.89.exons.protein_coding.txt", sep = "\t", row.names=F, quote=F)

exon_data.bed = exon_data %>% dplyr::select(chromosome_name, exon_chrom_start, exon_chrom_end, ensembl_gene_id) %>%
  dplyr::rename(chr = chromosome_name, start = exon_chrom_start, end = exon_chrom_end, gene_id = ensembl_gene_id)
write.table(exon_data.bed, "../../../reference/GRCh38/Homo_sapiens.GRCh38.89.exons.bed", sep = "\t", row.names=F, col.names=F, quote=F)

exon_data.pc.bed = exon_data.pc %>% dplyr::select(chromosome_name, exon_chrom_start, exon_chrom_end, gene_id) %>%
  dplyr::rename(chr = chromosome_name, start = exon_chrom_start, end = exon_chrom_end)
write.table(exon_data.pc.bed, "../../../reference/GRCh38/Homo_sapiens.GRCh38.89.exons.protein_coding.bed", sep = "\t", row.names=F, col.names=F, quote=F)


#Extract GENCODE basic transcript ids
# gencode_basic = dplyr::tbl_df(transcript_data) %>% 
#   dplyr::filter(transcript_gencode_basic == "GENCODE basic") %>% 
#   dplyr::select(ensembl_gene_id, ensembl_transcript_id)
# write.table(gencode_basic, "../../annotations/GRCh38/genes/Ensembl_87/Homo_sapiens.GRCh38.87.gencode_basic.txt", sep = "\t", row.names = FALSE, quote = FALSE)



