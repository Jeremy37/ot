#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(magrittr)
options(stringsAsFactors = F)

#root = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys"
root = "/Users/jeremys/work/opentargets"

coloc.macrophage.eqtl.df = readr::read_tsv(file.path(root, "coloc", "ipsmacrophage", "coloc.AD.meta.featureCounts.1e+05.txt"))
coloc.macrophage.sqtl.df = readr::read_tsv(file.path(root, "coloc", "ipsmacrophage", "coloc.AD.meta.leafcutter.1e+05.txt.ann.txt"))
coloc.macrophage.reviseAnn.df = readr::read_tsv(file.path(root, "coloc", "ipsmacrophage", "coloc.AD.meta.reviseAnnotations.1e+05.txt"))

coloc.mono.eqtl.df = readr::read_tsv(file.path(root, "coloc", "blueprint", "coloc.AD.meta.mono_gene_eQTL.1e+05.txt"))

getTopColocs = function(coloc.df, useCondition=F) {
  # Remove any duplicate entries with the same gwas_lead and qtl_lead SNP.
  coloc.df %<>% dplyr::arrange(qtl_pval) %>%
    dplyr::mutate(sortkey = paste(gwas_lead, qtl_lead, condition, ensemblID, sep=".")) %>%
    dplyr::filter(!duplicated(sortkey)) %>%
    dplyr::select(-sortkey) %>%
    dplyr::arrange(-PP.H4) %>%
    dplyr::mutate(H4_rel_H3H4 = PP.H4 / (PP.H3 + PP.H4)) %>%
    dplyr::rename(H4 = PP.H4) %>%
    dplyr::filter(gwas_pval < 1e-7)
  if (useCondition) {
    coloc.df$colocStr = sprintf("%s,%.2g,%.2g,%s", coloc.df$geneSymbol, coloc.df$H4, coloc.df$H4_rel_H3H4, coloc.df$condition)
  } else {
    coloc.df$colocStr = sprintf("%s,%.2g,%.2g", coloc.df$geneSymbol, coloc.df$H4, coloc.df$H4_rel_H3H4)
  }
  
  coloc.df %<>% dplyr::filter(qtl_pval < 1e-3)
  topColocs = data.frame(lead = unique(coloc.df$gwas_lead), colocStr="")
  for (i in 1:nrow(topColocs)) {
    coloc.locus.df = coloc.df %>% filter(gwas_lead == topColocs$lead[i])
    topColocs[i,]$colocStr = paste(coloc.locus.df$colocStr, collapse=" / ")
  }
  topColocs %>% dplyr::filter()
}


AD.annotated.df = readr::read_tsv(file.path(root, "gwas", "AD", "toby.jimmy.finemap.annotated.roadmapEnhPromLinks.txt"))

# Make a table mapping from each SNP to the lead SNP at a locus
AD.annotated.df$gwas_lead = NA
for (locus in unique(AD.annotated.df$signal)) {
  sig.df = AD.annotated.df %>% filter(signal == locus) %>%
    dplyr::arrange(finemap.p)
  AD.annotated.df[AD.annotated.df$signal == locus,]$gwas_lead = sig.df[1,]$snp
}
#########################
topColocs.eqtl = getTopColocs(coloc.macrophage.eqtl.df, useCondition=T)
topColocs.eqtl %<>% dplyr::left_join(AD.annotated.df %>% dplyr::select(snp, gwas_lead),
                                     by=c("lead" = "snp")) %>%
  dplyr::filter(!is.na(gwas_lead) & !duplicated(gwas_lead)) %>%
  dplyr::select(gwas_lead, colocStr)
AD.annotated.df %<>% dplyr::left_join(topColocs.eqtl, by=c("snp" = "gwas_lead")) %>%
  dplyr::rename(mac_eQTL_coloc = colocStr)
#########################
topColocs.sqtl = getTopColocs(coloc.macrophage.sqtl.df, useCondition=T)
topColocs.sqtl %<>% dplyr::left_join(AD.annotated.df %>% dplyr::select(snp, gwas_lead),
                                     by=c("lead" = "snp")) %>%
  dplyr::filter(!is.na(gwas_lead) & !duplicated(gwas_lead)) %>%
  dplyr::select(gwas_lead, colocStr)
AD.annotated.df %<>% dplyr::left_join(topColocs.sqtl, by=c("snp" = "gwas_lead")) %>%
  dplyr::rename(mac_sQTL_coloc = colocStr)
#########################
topColocs.reviseAnn = getTopColocs(coloc.macrophage.reviseAnn.df, useCondition=T)
topColocs.reviseAnn %<>% dplyr::left_join(AD.annotated.df %>% dplyr::select(snp, gwas_lead),
                                     by=c("lead" = "snp")) %>%
  dplyr::filter(!is.na(gwas_lead) & !duplicated(gwas_lead)) %>%
  dplyr::select(gwas_lead, colocStr)
AD.annotated.df %<>% dplyr::left_join(topColocs.reviseAnn, by=c("snp" = "gwas_lead")) %>%
  dplyr::rename(mac_reviseAnn_coloc = colocStr)
#########################
blueprintColocs.eqtl = getTopColocs(coloc.mono.eqtl.df, useCondition=F)
blueprintColocs.eqtl %<>% dplyr::left_join(AD.annotated.df %>% dplyr::select(snp, gwas_lead),
                                           by=c("lead" = "snp")) %>%
  dplyr::filter(!is.na(gwas_lead) & !duplicated(gwas_lead)) %>%
  dplyr::select(gwas_lead, colocStr)
AD.annotated.df %<>% dplyr::left_join(blueprintColocs.eqtl, by=c("snp" = "gwas_lead")) %>%
  dplyr::rename(mono_eQTL_coloc = colocStr,
                enhPromScores = geneScores, enhPromScoreDetails = geneScoreDetails) %>%
  dplyr::select(-gwas_lead)


write.table(AD.annotated.df, file=file.path(root, "gwas", "toby.jimmy.finemap.annotated.EnhProm.coloc.txt"),
            sep="\t", row.names=F, col.names=T, quote=F, na="")
