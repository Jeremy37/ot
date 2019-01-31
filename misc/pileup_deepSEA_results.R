#!/usr/bin/env Rscript
library(tidyverse)
options(stringsAsFactors=F)

args = commandArgs(trailingOnly = TRUE)
deepSEA.root = args[1] # Root name for DeepSEA files
celltype_filter = args[2]

# There should be 7 files:
# root.diff
# root.funsig
# root.evalue
# root.logfoldchange
# root.snpclass
# root.alt
# root.ref

# celltype_filter should a regex, e.g. be of the format below separated with bar |:
# celltype_filter="GM12878|BE2_C|Gliobla|H1-hESC|H9ES|HBMEC|iPS|Medullo|K562|Monocytes-CD14+RO01746|SH-SY5Y|SK-N-MC|SK-N-SH"

#setwd("/Users/jeremys/Downloads/")
#deepSEA.root = "DeepSEA.CCDC6_rs1171832_13bp_satmut/infile.vcf.out"

diff.df = readr::read_csv(paste0(deepSEA.root, ".diff")) %>%
  dplyr::mutate(colname = "Diff")

logfc.df = readr::read_csv(paste0(deepSEA.root, ".logfoldchange")) %>%
  dplyr::mutate(colname = "Log2FC")

evalue.df = readr::read_csv(paste0(deepSEA.root, ".evalue")) %>%
  dplyr::mutate(colname = "E_value")

ref.df = readr::read_csv(paste0(deepSEA.root, ".ref")) %>%
  dplyr::mutate(colname = "P_reference")

alt.df = readr::read_csv(paste0(deepSEA.root, ".alt")) %>%
  dplyr::mutate(colname = "P_alt")

df = bind_rows(diff.df, logfc.df, evalue.df, ref.df, alt.df) %>%
  tidyr::gather(assay, value, -(X1:alt), -colname) %>%
  tidyr::spread(colname, value) %>%
  dplyr::rename(index = X1)

# Split the cell type / assay / treatment column into separate cols
df$cell_type = sapply(df$assay, FUN=function(s) strsplit(s, "|", fixed=T)[[1]][1])
df$feature = sapply(df$assay, FUN=function(s) strsplit(s, "|", fixed=T)[[1]][2])
df$treatment = sapply(df$assay, FUN=function(s) strsplit(s, "|", fixed=T)[[1]][3])

# Filter cell types now
if (!is.na(celltype_filter)) {
  df = df %>% dplyr::filter(grepl(pattern=celltype_filter, x=cell_type, perl=T))
}

funsig.df = readr::read_csv(paste0(deepSEA.root, ".funsig")) %>%
  dplyr::rename(funsig_score = `DeepSEA score`, index = X1)

summary.df = df %>% dplyr::group_by(index, cell_type) %>%
  dplyr::summarise(max_abs_log2fc = max(abs(Log2FC)),
                   min_evalue = min(E_value),
                   max_p_ref = max(P_reference),
                   max_p_alt = max(P_alt))

out.df = df %>%
  dplyr::select(-assay) %>%
  dplyr::left_join(funsig.df %>% dplyr::select(index, funsig_score), by="index") %>%
  dplyr::left_join(summary.df, by=c("index", "cell_type")) %>%
  dplyr::arrange(funsig_score, index, -abs(Log2FC))

write.table(out.df, "", quote=F, sep="\t", na="", row.names=F, col.names=T)
#write.table(out.df, "DeepSEA.CCDC6_rs1171832_13bp_satmut/merged.tsv", quote=F, sep="\t", na="", row.names=F, col.names=T)
