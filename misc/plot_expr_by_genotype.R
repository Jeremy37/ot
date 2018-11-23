#!/usr/bin/env Rscript
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
counts_file = args[1]
vcf_file = args[2]
sample_metadata_file = args[3] # Maps from cols in counts_file to cols in vcf_file
gene_metadata_file = args[4]
gene_id = args[5]
snp_id = args[6]
log_axis_param = args[7]

log_y_axis = F
if (!is.na(log_axis_param)) {
  log_y_axis = T
}
# counts_file = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/macrophage/RNA_count_matrix.txt.gz"
# vcf_file = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/macrophage/genotypes/imputed.86_samples.sorted.filtered.named.vcf.gz"
# sample_metadata_file = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/macrophage/RNA_sample_metadata.txt.gz"
# gene_metadata_file = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/macrophage/RNA_gene_metadata.txt.gz"

# gene_id_in = "ENSG00000120899" # PTK2B
# gene_id_in = "ENSG00000120885" # CLU
# gene_id_in = "ENSG00000168077" # SCARA3
# snp_id = "rs28834970"

genes.df = readr::read_tsv(gene_metadata_file)

counts = read.delim(counts_file, sep="\t")
rpm = apply(counts, MARGIN=2, FUN=function(x) x / sum(x)) * 1e6
rpkm = rpm * 1e3
rpkm = apply(rpkm, MARGIN=2, FUN=function(x) x / genes.df$length)

# Get table mapping from count table IDs to VCF IDs
sample.meta.df = readr::read_tsv(sample_metadata_file) %>%
  dplyr::select(sample_id, donor, condition, condition_name, genotype_id)
sum(colnames(counts) %in% sample.meta.df$sample_id)

plotGene = function(gene_id_in, snp_id, log_y_axis = F) {
  # Just extract a 1-Mb region around the gene using tabix to try to find the SNP
  gene_row = genes.df %>% dplyr::filter(gene_id == gene_id_in)
  gene_chr = gene_row$chr[1]
  gene_start = gene_row$start[1]
  gene_end = gene_row$end[1]
  vcf_cmd = sprintf("gunzip -c %s | head -n 5000 | grep -e '#CHROM'; tabix %s %s:%d-%d | grep -e '%s\t'", vcf_file, vcf_file, gene_chr, gene_start-1000000, gene_end+1000000, snp_id)
  
  snp.df = readr::read_tsv(file = pipe(vcf_cmd))
  genotype.df = snp.df[, 10:ncol(snp.df)] %>% tidyr::gather(key = "genotype_id", value = "gt")
  genotype.df$snp_count = sapply(genotype.df$gt, FUN = function(s) sum(strtoi(strsplit(s, "[|/:]")[[1]][1:2])))
  
  gene.rpkm.df = as.data.frame(rpkm[gene_id_in, , drop=F]) %>% tidyr::gather(key = "sample_id", value = "rpkm")
  
  conditions = unique(sample.meta.df$condition)
  condition_names = unique(sample.meta.df$condition_name)
  
  gene.rpkm.df = gene.rpkm.df %>%
    dplyr::left_join(sample.meta.df, by = "sample_id") %>%
    dplyr::left_join(genotype.df, by = "genotype_id")
  gene.rpkm.df$snp_count = factor(as.character(gene.rpkm.df$snp_count), levels=c("0", "1", "2"))
  p = ggplot(gene.rpkm.df, aes(x=snp_count, y=rpkm)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.15, alpha=0.6) +
    facet_wrap(~condition_name) +
    theme_bw() + ggtitle(gene_row$gene_name)
  if (log_y_axis) {
    p = p + scale_y_log10()
  }
  p
}

print(plotGene(gene_id_in, log_y_axis))

# log_y_axis = F
# pdf(sprintf("%s.gene_expr.pdf", snp_id), width=7, height=7)
# print(plotGene("ENSG00000120899", snp_id = snp_id, log_y_axis))
# print(plotGene("ENSG00000120885", snp_id = snp_id, log_y_axis))
# print(plotGene("ENSG00000168077", snp_id = snp_id, log_y_axis))
# dev.off()
# 
# log_y_axis = T
# pdf(sprintf("%s.gene_expr.log10.pdf", snp_id), width=7, height=7)
# print(plotGene("ENSG00000120899", snp_id = snp_id, log_y_axis))
# print(plotGene("ENSG00000120885", snp_id = snp_id, log_y_axis))
# print(plotGene("ENSG00000168077", snp_id = snp_id, log_y_axis))
# dev.off()

