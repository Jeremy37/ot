#!/usr/bin/env Rscript
library(tidyverse)

root = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys"
#root = "/Users/jeremys/work/opentargets"
#basename = "AD.finemap"
basename = "AD.finemap.BIN1"

# Change the RPKM columns to only use 3 positions after the decimal
tissue.rpkm.df = readr::read_tsv(file.path(root, "reference/tissueRPKM/tissues.selected.rpkm_average.txt.gz"))
for (col in 2:11) {
  tissue.rpkm.df[,col] = sapply(tissue.rpkm.df[,col], function(x) sprintf("%.3f", x))
}
colnames(tissue.rpkm.df)[-1] = paste0("rpkm.", colnames(tissue.rpkm.df)[-1])

# Read in AD file of VEP annotated transcribed SNPs, but retain only
# distinct SNPs (ignoring different transcripts affected)
df = readr::read_tsv(file.path(root, "gwas/AD/finemap/", paste0(basename, ".vep.transcribed.pcausal_gt_0.locusnames.txt"))) %>%
  dplyr::rename(snp_id = `#Uploaded_variation`, gene_id = Gene, pos_hg38 = Location) %>%
  dplyr::distinct(snp_id, .keep_all = T)
df$pos_hg38 = sapply(df$pos_hg38, FUN = function(x)strsplit(x, split=":|-", perl = T)[[1]][2])

# Use the summary stats file to add in the chr:pos fields in hg19 coords
snpChrPos.df = readr::read_tsv(file.path(root, "gwas/AD/AD.IGAP1_GWAX.meta.chrpos.txt.gz")) %>%
  dplyr::rename(snp_id = SNP, pos_hg19 = BP, Chr = CHR)
df.withpos = dplyr::left_join(df, snpChrPos.df) %>%
  dplyr::select(snp_id, Chr, pos_hg19, pos_hg38, everything())

df.ann = dplyr::left_join(df.withpos, tissue.rpkm.df, by="gene_id")
write.table(df.ann, file.path(root, "gwas/AD/finemap/", paste0(basename, ".vep.transcribed.pcausal_gt_0.locusnames.rpkms.txt")), row.names=F, col.names=T, quote=F, sep="\t")


###############################################################################
# PD
df = readr::read_tsv(file.path(root, "gwas/PD/finemap/PD.finemap.vep.transcribed.pcausal_gt_0.txt")) %>%
  dplyr::rename(snp_id = `#Uploaded_variation`, gene_id = Gene, pos_hg38 = Location) %>%
  dplyr::distinct(snp_id, .keep_all = T)
df$pos_hg38 = sapply(df$pos_hg38, FUN = function(x)strsplit(x, split=":|-", perl = T)[[1]][2])

# Use the credible sets file to add in the chr:pos fields in hg19 coords
#credsetFile = file.path(root, "gwas/PD/finemap/PD.credible_sets.set")
credsetFile = file.path(root, "gwas/PD/finemap/pd_gwas_gwax.combined.set")
snpChrPos.df = readr::read_tsv(credsetFile) %>%
  dplyr::select(SNP, Chr, pos) %>%
  dplyr::rename(snp_id = SNP, pos_hg19 = pos)
df.withpos = dplyr::left_join(df, snpChrPos.df) %>%
  dplyr::select(snp_id, Chr, pos_hg19, pos_hg38, everything())

df.ann = dplyr::left_join(df.withpos, tissue.rpkm.df, by="gene_id") %>%
  dplyr::select(locus, snp_id, Chr, pos_hg19, pos_hg38, Allele, p, prob, everything(), -beta, -se, -pos, -Feature_type)
write.table(df.ann, file.path(root, "gwas/PD/finemap/PD.finemap.vep.transcribed.pcausal_gt_0.rpkms.txt"), row.names=F, col.names=T, quote=F, sep="\t")
