#!/usr/bin/env Rscript
library(tidyverse)
options(stringsAsFactors = F)

root = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys"
#root = "/Users/jeremys/work/opentargets"

df.runtable = read.delim(file.path(root, "datasets/microglia/SraRunTable.txt"))

phenotype_file = file.path(root, "datasets/microglia/PhenoGenotypeFiles/RootStudyConsentSet_phs001373.HumanMicroglia.v1.p1.c1.HMB/PhenotypeFiles",
                           "phs001373.v1.pht006703.v1.p1.c1.Human_Microglia_Subject_Phenotypes.HMB.txt")
phenotype.df = read.delim(phenotype_file, blank.lines.skip = T, comment.char = "#")

df.meta = df.runtable %>% dplyr::select(-DATASTORE_filetype, -DATASTORE_provider, -InsertSize, -sex, SUBJECT_ID = submitted_subject_id) %>%
  dplyr::left_join(phenotype.df, by="SUBJECT_ID") %>%
  dplyr::select(SampleID = Run, SUBJECT_ID, everything())

write.table(df.meta, file.path(root, "datasets/microglia/glass_microglia.metadata.txt"), row.names=F, col.names=T, quote=F, sep="\t")

write.table(df.meta %>% dplyr::filter(Assay_Type == "ATAC-seq"),
            file.path(root, "datasets/microglia/glass_microglia.metadata.atac.txt"), row.names=F, col.names=T, quote=F, sep="\t")

df.fastq = df.meta %>% dplyr::select(SampleID, Library_Name)
df.fastq$fastq = paste0(root, "/datasets/microglia/fastq/", df.fastq$SampleID, ".fastq.gz")
df.fastq = df.fastq %>% dplyr::select(SampleID, fastq, Library_Name)
write.table(df.fastq, file.path(root, "datasets/microglia/glass_microglia.sample.fastq.txt"), row.names=F, col.names=T, quote=F, sep="\t")

df.fastq = df.meta %>% dplyr::select(SampleID, Library_Name)
df.fastq$fastq = paste0(root, "/datasets/microglia/fastq/", df.fastq$SampleID, ".fastq.gz")
write.table(df.fastq %>% dplyr::filter(grepl("ATAC", Library_Name)) %>% dplyr::select(-Library_Name),
            file.path(root, "datasets/microglia/glass_microglia.atac.fastq.txt"),
            row.names=F, col.names=T, quote=F, sep="\t")
