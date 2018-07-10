#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(magrittr)
options(stringsAsFactors = F)

root = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys"
#root = "/Users/jeremys/work/opentargets"

prefix = "coloc.AD.meta"
outputroot = file.path(root, "gwas", "AD", "AD.finemap.annotated")
annotated_gwas = file.path(root, "gwas", "AD", "JEME", "AD.finemap.annotated.roadmapEnhPromLinks.txt")
output_file = paste0(outputroot, ".colocs.txt")

# prefix = "coloc.PD.meta"
# outputroot = file.path(root, "gwas", "PD", "PD.finemap.annotated")
# annotated_gwas = file.path(root, "gwas", "PD", "JEME", "PD.finemap.annotated.roadmapEnhPromLinks.txt")
# output_file = paste0(outputroot, ".colocs.txt")

###############################################################################
# Macrophage
mac_conditions = c("naive", "IFNg", "SL1344", "IFNg_SL1344")
coloc.macrophage.eqtl.df = data.frame()
coloc.macrophage.sqtl.df = data.frame()
coloc.macrophage.txrevise_promoters.df = data.frame()
coloc.macrophage.txrevise_ends.df = data.frame()
for (condition in mac_conditions) {
  eqtl.df = readr::read_tsv(file.path(root, "coloc", "new", sprintf("%s.ips_macrophage.eQTL_%s.5e+05.txt", prefix, condition))) %>%
    dplyr::mutate(dataset = paste0("Macrophage eQTL ", condition), dataset_short = paste0("mac_eqtl_", condition), ensembl_id = feature, geneSymbol = NA)
  coloc.macrophage.eqtl.df = bind_rows(coloc.macrophage.eqtl.df, eqtl.df)
  
  sqtl.df = readr::read_tsv(file.path(root, "coloc", "new", sprintf("%s.ips_macrophage.leafcutter_QTL_%s.5e+05.txt.ann.txt", prefix, condition))) %>%
    dplyr::mutate(dataset = paste0("Macrophage sQTL ", condition), dataset_short = paste0("mac_sqtl_", condition)) %>%
    dplyr::rename(feature = clusterID, ensembl_id = geneid) %>%
    dplyr::select(-ensembl_id, -geneSymbol, everything(), ensembl_id, geneSymbol)
  coloc.macrophage.sqtl.df = bind_rows(coloc.macrophage.sqtl.df, sqtl.df)
  
  txrevise_promoters.df = readr::read_tsv(file.path(root, "coloc", "new", sprintf("%s.ips_macrophage.txrevise_promoters_%s.5e+05.txt", prefix, condition))) %>%
    dplyr::mutate(dataset = paste0("Macrophage TxRevise promoters ", condition), dataset_short = paste0("mac_txrevise_prom_", condition), ensembl_id = NA, geneSymbol = NA)
  coloc.macrophage.txrevise_promoters.df = bind_rows(coloc.macrophage.txrevise_promoters.df, txrevise_promoters.df)
  
  txrevise_ends.df = readr::read_tsv(file.path(root, "coloc", "new", sprintf("%s.ips_macrophage.txrevise_ends_%s.5e+05.txt", prefix, condition))) %>%
    dplyr::mutate(dataset = paste0("Macrophage TxRevise ends ", condition), dataset_short = paste0("mac_txrevise_ends_", condition), ensembl_id = NA, geneSymbol = NA)
  coloc.macrophage.txrevise_ends.df = bind_rows(coloc.macrophage.txrevise_ends.df, txrevise_ends.df)
}

# Update the ensembl ID for each of these
coloc.macrophage.txrevise_promoters.df$ensembl_id = sapply(coloc.macrophage.txrevise_promoters.df$feature, FUN = function(x) {strsplit(x, "\\.")[[1]][1]})
coloc.macrophage.txrevise_ends.df$ensembl_id = sapply(coloc.macrophage.txrevise_ends.df$feature, FUN = function(x) {strsplit(x, "\\.")[[1]][1]})


###############################################################################
# Blueprint
coloc.mono.eqtl.df = readr::read_tsv(file.path(root, "coloc", "new", paste0(prefix, ".blueprint.mono_gene_eQTL.5e+05.txt"))) %>%
  dplyr::mutate(dataset = "Blueprint monocyte eQTL", dataset_short = "mono_eQTL",
                ensembl_id = feature, geneSymbol = NA)

coloc.mono.h3k27ac.df = readr::read_tsv(file.path(root, "coloc", "new", paste0(prefix, ".blueprint.mono_K27AC.5e+05.txt"))) %>%
  dplyr::mutate(dataset = "Blueprint monocyte H3K27ac", dataset_short = "mono_h3k27ac",
                ensembl_id = NA, geneSymbol = paste0("peak_", feature))

coloc.mono.h3k4me1.df = readr::read_tsv(file.path(root, "coloc", "new", paste0(prefix, ".blueprint.mono_K4ME1.5e+05.txt"))) %>%
  dplyr::mutate(dataset = "Blueprint monocyte H3K4me1", dataset_short = "mono_h3k4me1",
                ensembl_id = NA, geneSymbol = paste0("peak_", feature))


###############################################################################
# xQTL
# xQTL colocs have the gene symbol in the EnsemblID column
coloc.xqtl.eqtl.df = readr::read_tsv(file.path(root, "coloc", "new", paste0(prefix, ".xQTL_eQTL.5e+05.txt"))) %>%
  dplyr::mutate(dataset = "Brain ROSMAP eQTL", dataset_short = "xQTL_eQTL",
                ensembl_id = NA, geneSymbol = feature)

# When getting haQTL and mQTLs, we want to provide some indication as to which
# peak / CpG a QTL refers to. We need to extract this from another file, in this
# case the *.minp file for each peak.
coloc.xqtl.haqtl.df = readr::read_tsv(file.path(root, "coloc", "new", paste0(prefix, ".xQTL_haQTL.5e+05.txt"))) %>%
  dplyr::mutate(dataset = "Brain ROSMAP H3K27ac", dataset_short = "xQTL_h3k27ac", ensembl_id = NA)
haqtl.colnames = c("featureName", "SNPchromosome", "SNPpos", "SNPid", "featureChromosome", "featurePositionStart", "SpearmanRho", "pValue")
xqtl.haqtl.minp.df =  readr::read_tsv(file.path(root, "reference", "xQTL", "xQTL_haQTL.gene_minp.txt.gz"), col_names = haqtl.colnames)
xqtl.haqtl.minp.df = xqtl.haqtl.minp.df %>% dplyr::mutate(geneSymbol = paste(featureName, featureChromosome, featurePositionStart, sep = "_"))
coloc.xqtl.haqtl.df = coloc.xqtl.haqtl.df %>%
  dplyr::left_join(xqtl.haqtl.minp.df %>% dplyr::select(featureName, geneSymbol), by=c("feature" = "featureName"))

coloc.xqtl.mqtl.df = readr::read_tsv(file.path(root, "coloc", "new", paste0(prefix, ".xQTL_mQTL.5e+05.txt"))) %>%
  dplyr::mutate(dataset = "Brain ROSMAP DNA methylation", dataset_short = "xQTL_meth", ensembl_id = NA)
mqtl.colnames = c("featureName", "SNPchromosome", "SNPpos", "SNPid", "featureChromosome", "featurePositionStart", "SpearmanRho", "pValue")
xqtl.mqtl.minp.df =  readr::read_tsv(file.path(root, "reference", "xQTL", "xQTL_mQTL.gene_minp.txt.gz"), col_names = mqtl.colnames)
xqtl.mqtl.minp.df = xqtl.mqtl.minp.df %>% dplyr::mutate(geneSymbol = paste(featureName, featureChromosome, featurePositionStart, sep = "_"))
coloc.xqtl.mqtl.df = coloc.xqtl.mqtl.df %>%
  dplyr::left_join(xqtl.mqtl.minp.df %>% dplyr::select(featureName, geneSymbol), by=c("feature" = "featureName"))

rm(xqtl.haqtl.minp.df, xqtl.mqtl.minp.df)


###############################################################################
# Sensory neurons
coloc.sensneur.eqtl.df = readr::read_tsv(file.path(root, "coloc", "new", paste0(prefix, ".sens_neur.eqtl.500k.5e+05.txt"))) %>%
  dplyr::mutate(dataset = "Sensory neuron eQTL", dataset_short = "sensneur_eqtl",
                ensembl_id = feature, geneSymbol = NA)

coloc.sensneur.sqtl.df = readr::read_tsv(file.path(root, "coloc", "new", paste0(prefix, ".sens_neur.sqtl.5e+05.txt.ann.txt"))) %>%
  dplyr::mutate(dataset = paste0("Sensory neuron sQTL ", condition), dataset_short = "sensneur_sqtl") %>%
  dplyr::rename(feature = clusterID, ensembl_id = geneid) %>%
  dplyr::select(-ensembl_id, -geneSymbol, everything(), ensembl_id, geneSymbol)


# # Use Ensembl to HGNC map to add in gene symbol
# coloc.sensneur.eqtl.df %<>% dplyr::left_join(ensemblMap, by=c("phenotype_id" = "ensembl_gene_id")) %>%
#   dplyr::filter(!duplicated(phenotype_id, symbol)) %>%
#   dplyr::select(-symbol_type) %>%
#   dplyr::rename(geneSymbol = symbol) %>%
#   dplyr::mutate(dataset = "Sensory neuron eQTL", dataset_short = "sensneur_eqtl")
# 
# # For sensory neuron sQTLs, get "gene symbol" to use from the original lead SNP file
# sensneurLeadSQTL.df = readr::read_tsv(file.path(root, "sensoryneurons", "GRCh38", "sqtl", "fastqtl.sqtl.permutations.10k.PCs.5.fdr0.1.txt"))
# geneSummary = function(i) {
#   if (grepl("[A-Za-z]", sensneurLeadSQTL.df[i,]$symbols)) {
#     return(sensneurLeadSQTL.df[i,]$symbols)
#   }
#   return(sensneurLeadSQTL.df[i,]$ensembl_ids)
# }
# sensneurLeadSQTL.df$geneSymbol = sapply(1:nrow(sensneurLeadSQTL.df), FUN = geneSummary)
# coloc.sensneur.sqtl.df %<>% dplyr::left_join(sensneurLeadSQTL.df %>% dplyr::select(phenotype_id, geneSymbol),
#                                              by="phenotype_id") %>%
#   dplyr::mutate(dataset = "Sensory neuron sQTL", dataset_short = "sensneur_sqtl")



###############################################################################
# Make a big table of all the colocs

# Write out a smaller version of the table, including each region-coloc as one row.
colocs.df = rbind(coloc.macrophage.eqtl.df,
                  coloc.macrophage.sqtl.df,
                  coloc.macrophage.txrevise_promoters.df,
                  coloc.macrophage.txrevise_ends.df,
                  coloc.mono.eqtl.df,
                  coloc.mono.h3k27ac.df,
                  coloc.mono.h3k4me1.df,
                  coloc.xqtl.eqtl.df,
                  coloc.xqtl.haqtl.df,
                  coloc.xqtl.mqtl.df,
                  coloc.sensneur.eqtl.df,
                  coloc.sensneur.sqtl.df
) %>% dplyr::rename(gwas_snp = gwas_lead)


###############################################################################
# Use Ensembl to HGNC map to add in gene symbol where it isn't already present

ensemblMap = readr::read_tsv(file.path(root, "reference", "hgnc.ensembl.map.txt")) %>%
  dplyr::filter(symbol_type == "symbol")
colocs.df %<>% dplyr::left_join(ensemblMap, by=c("ensembl_id" = "ensembl_gene_id")) %>%
  dplyr::select(-symbol_type)
colocs.df$geneSymbol[is.na(colocs.df$geneSymbol)] = colocs.df$symbol[is.na(colocs.df$geneSymbol)]
colocs.df = colocs.df %>% dplyr::select(-symbol)



###############################################################################
# Identify the GWAS signal that each coloc correspond with
# First read in the GWAS fine mapping table
annotated.df = readr::read_tsv(annotated_gwas) %>%
  dplyr::rename(enhPromScores = geneScores, enhPromScoreDetails = geneScoreDetails)

# UGLY CODE - harmonize column names between AD and PD gwas annotated files
if ("signal" %in% colnames(annotated.df)) {
  annotated.df %<>% dplyr::rename(locus = signal)
}
if (!("finemap.p" %in% colnames(annotated.df))) {
  annotated.df %<>% dplyr::rename(finemap.p = p, snp = SNP)
}

annotated.df$gwas_lead = NA
for (locusname in unique(annotated.df$locus)) {
  sig.df = annotated.df %>% filter(locus == locusname) %>%
    dplyr::arrange(finemap.p)
  annotated.df[annotated.df$locus == locusname,]$gwas_lead = sig.df[1,]$snp
}
signals.df = annotated.df %>% 
  dplyr::mutate(chrpos = paste0(Chr, ".", pos)) %>%
  dplyr::select(locus, snp, chrpos, gwas_lead)

# First replace the coloc GWAS SNPs with the rsID from the annotated table,
# matching on chr:pos. This is necessary because the coloc code changes the
# GWAS lead SNP IDs to those from the QTL variant info, and if the QTL variant
# info used a different ID (e.g. chr17_40051621 rather than rs...) then the
# IDs won't match.
colocs.df$chrpos = paste0(colocs.df$chr, ".", colocs.df$gwas_pos)
colocs.df %<>% dplyr::left_join(signals.df, by="chrpos")
colocs.df %<>% dplyr::select(-chrpos)

# Remove any duplicate entries with the same gwas_lead and qtl_lead SNP,
# and arrange in order of decreasing H4 probability
colocs.df %<>% dplyr::arrange(qtl_pval) %>%
  dplyr::mutate(sortkey = paste(gwas_lead, qtl_lead, feature, sep=".")) %>%
  dplyr::filter(!duplicated(sortkey)) %>%
  dplyr::select(-sortkey) %>%
  dplyr::mutate(H4_rel_H3H4 = PP.H4 / (PP.H3 + PP.H4)) %>%
  dplyr::rename(H4 = PP.H4) %>%
  dplyr::arrange(-H4) %>%
  dplyr::filter(gwas_pval < 1e-7) %>%
  dplyr::filter(qtl_pval < 1e-3)

colocs.top.df  = colocs.df %>%
  dplyr::filter(H4 > 0.5 | (H4 > 0.2 & H4_rel_H3H4 > 0.9))
# AD.colocs.top.df %<>% dplyr::filter(dataset_short != "xQTL_meth")

# Summarize colocs for a locus across all QTL datasets
getColocsPerLocus = function(coloc.df, includeDataset = F) {
  coloc.df %<>% dplyr::arrange(-H4, -H4_rel_H3H4)
  
  if (includeDataset) {
    coloc.df$colocStr = sprintf("%s,%s,%.2g,%.2g", coloc.df$geneSymbol, coloc.df$dataset_short, coloc.df$H4, coloc.df$H4_rel_H3H4)
  } else {
    coloc.df$colocStr = sprintf("%s,%.2g,%.2g", coloc.df$geneSymbol, coloc.df$H4, coloc.df$H4_rel_H3H4)
  }
  topColocs = data.frame(gwas_lead = unique(coloc.df$gwas_lead), colocStr="")
  for (i in 1:nrow(topColocs)) {
    coloc.locus.df = coloc.df %>% filter(gwas_lead == topColocs$gwas_lead[i])
    topColocs[i,]$colocStr = paste(coloc.locus.df$colocStr, collapse=" / ")
  }
  topColocs
}

topColocs.df = getColocsPerLocus(colocs.top.df, includeDataset = T)

annotated.df %<>% dplyr::left_join(topColocs.df, by=c("snp" = "gwas_lead")) %>%
  dplyr::rename(topColocs = colocStr)

# Write out the full table, with the top colocs annotated only in the rows of GWAS
# lead SNPs
write.table(annotated.df %>% dplyr::select(-gwas_lead) %>% dplyr::arrange(Chr, locus, finemap.p),
            file=output_file,
            sep="\t", row.names=F, col.names=T, quote=F, na="")



# Write out a more detailed table including all colocs
coloc.summary.df = data.frame()
for (curDataset in unique(colocs.df$dataset)) {
  topColocs.df = getColocsPerLocus(colocs.df %>% dplyr::filter(dataset == curDataset))
  topColocs.df$dataset = curDataset
  coloc.summary.df = rbind(coloc.summary.df, topColocs.df)
}

signals.df %<>% dplyr::distinct(locus, gwas_lead) %>% na.omit()
coloc.summary.df = coloc.summary.df %>% dplyr::left_join(signals.df, by="gwas_lead") %>%
  na.omit() %>%
  dplyr::left_join(annotated.df %>% dplyr::select(snp, Chr, pos),
                   by=c("gwas_lead" = "snp")) %>%
  dplyr::arrange(Chr, pos, dataset) %>%
  dplyr::rename(colocDetails = colocStr) %>%
  dplyr::select(locus, gwas_lead, dataset, colocDetails)

write.table(coloc.summary.df, file=paste0(outputroot, ".coloc_details.txt"),
            sep="\t", row.names=F, col.names=T, quote=F, na="")

