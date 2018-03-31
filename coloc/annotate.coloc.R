#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(magrittr)
options(stringsAsFactors = F)

#root = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys"
root = "/Users/jeremys/work/opentargets"

# prefix = "coloc.AD.meta"
# outputroot = file.path(root, "gwas", "AD", "toby.jimmy.finemap.annotated")
# annotated_gwas = file.path(root, "gwas", "AD", "toby.jimmy.finemap.annotated.roadmapEnhPromLinks.txt")
# output_file = paste0(outputroot, ".colocs.txt")

prefix = "coloc.PD.meta"
outputroot = file.path(root, "gwas", "PD", "PD.finemap.annotated")
annotated_gwas = file.path(root, "gwas", "PD", "JEME", "PD.finemap.annotated.roadmapEnhPromLinks.txt")
output_file = paste0(outputroot, ".colocs.txt")

coloc.macrophage.eqtl.df = readr::read_tsv(file.path(root, "coloc", "ipsmacrophage", paste0(prefix, ".featureCounts.1e+05.txt"))) %>%
  dplyr::rename(phenotype_id = ensemblID) %>%
  dplyr::mutate(dataset = paste0("Macrophage eQTL ", condition), dataset_short = paste0("mac_eqtl_", condition)) %>%
  dplyr::select(-condition)
coloc.macrophage.sqtl.df = readr::read_tsv(file.path(root, "coloc", "ipsmacrophage", paste0(prefix, ".leafcutter.1e+05.txt.ann.txt"))) %>%
  dplyr::rename(phenotype_id = clusterID) %>%
  dplyr::mutate(dataset = paste0("Macrophage sQTL ", condition), dataset_short = paste0("mac_eqtl_", condition)) %>%
  dplyr::select(-condition, -geneid)
coloc.macrophage.reviseAnn.df = readr::read_tsv(file.path(root, "coloc", "ipsmacrophage", paste0(prefix, ".reviseAnnotations.1e+05.txt"))) %>%
  dplyr::rename(phenotype_id = ensemblID) %>%
  dplyr::mutate(dataset = paste0("Macrophage reviseAnnotations ", condition), dataset_short = paste0("mac_eqtl_", condition)) %>%
  dplyr::select(-condition, -reviseAnn)

coloc.macrophage.sqtl.df$geneSymbol[is.na(coloc.macrophage.sqtl.df$geneSymbol)] = coloc.macrophage.sqtl.df$phenotype_id[is.na(coloc.macrophage.sqtl.df$geneSymbol)]


coloc.mono.eqtl.df = readr::read_tsv(file.path(root, "coloc", "blueprint", paste0(prefix, ".mono_gene_eQTL.1e+05.txt")))
# Use Ensembl to HGNC map to add in gene symbol
ensemblMap = readr::read_tsv(file.path(root, "reference", "hgnc.ensembl.map.txt"))
coloc.mono.eqtl.df %<>% dplyr::left_join(ensemblMap, by=c("phenotype_id" = "ensembl_gene_id")) %>%
  dplyr::filter(!duplicated(phenotype_id, symbol)) %>%
  dplyr::select(-symbol_type) %>%
  dplyr::rename(geneSymbol = symbol) %>%
  dplyr::mutate(dataset = "Blueprint monocyte eQTL", dataset_short = "mono_eqtl")
  
coloc.mono.h3k27ac.df = readr::read_tsv(file.path(root, "coloc", "blueprint", paste0(prefix, ".mono_K27AC.1e+05.txt"))) %>%
  dplyr::mutate(dataset = "Blueprint monocyte H3K27ac", dataset_short = "mono_h3k27ac",
                geneSymbol = paste0("peak_", phenotype_id))

coloc.mono.h3k4me1.df = readr::read_tsv(file.path(root, "coloc", "blueprint", paste0(prefix, ".mono_K4ME1.1e+05.txt"))) %>%
  dplyr::mutate(dataset = "Blueprint monocyte H3K4me1", dataset_short = "mono_h3k4me1",
                geneSymbol = paste0("peak_", phenotype_id))

# xQTL colocs have the gene symbol in the EnsemblID column
coloc.xqtl.eqtl.df = readr::read_tsv(file.path(root, "coloc", "xQTL", paste0(prefix, ".xQTL_eQTL.1e+05.txt"))) %>%
  dplyr::mutate(dataset = "Brain ROSMAP eQTL", dataset_short = "xQTL_eQTL",
                geneSymbol = phenotype_id)

# When getting haQTL and mQTLs, we want to provide some indication as to which
# peak / CpG a QTL refers to. We need to extract this from another file, in this
# case the *.minp file for each peak.
coloc.xqtl.haqtl.df = readr::read_tsv(file.path(root, "coloc", "xQTL", paste0(prefix, ".xQTL_haQTL.1e+05.txt"))) %>%
  dplyr::mutate(dataset = "Brain ROSMAP H3K27ac", dataset_short = "xQTL_h3k27ac")
haqtl.colnames = c("featureName", "SNPchromosome", "SNPpos", "SNPid", "featureChromosome", "featurePositionStart", "SpearmanRho", "pValue")
xqtl.haqtl.minp.df =  readr::read_tsv(file.path(root, "reference", "xQTL", "xQTL_haQTL.gene_minp.txt.gz"), col_names = haqtl.colnames)
xqtl.haqtl.minp.df = xqtl.haqtl.minp.df %>% dplyr::mutate(geneSymbol = paste(featureName, featureChromosome, featurePositionStart, sep = "_")) %>%
  dplyr::rename(ensemblID = featureName)
coloc.xqtl.haqtl.df = coloc.xqtl.haqtl.df %>%
  dplyr::left_join(xqtl.haqtl.minp.df %>% dplyr::select(ensemblID, geneSymbol), by=c("phenotype_id" = "ensemblID"))

coloc.xqtl.mqtl.df = readr::read_tsv(file.path(root, "coloc", "xQTL", paste0(prefix, ".xQTL_mQTL.1e+05.txt"))) %>%
  dplyr::mutate(dataset = "Brain ROSMAP DNA methylation", dataset_short = "xQTL_meth")
mqtl.colnames = c("featureName", "SNPchromosome", "SNPpos", "SNPid", "featureChromosome", "featurePositionStart", "SpearmanRho", "pValue")
xqtl.mqtl.minp.df =  readr::read_tsv(file.path(root, "reference", "xQTL", "xQTL_mQTL.gene_minp.txt.gz"), col_names = mqtl.colnames)
xqtl.mqtl.minp.df = xqtl.mqtl.minp.df %>% dplyr::mutate(geneSymbol = paste(featureName, featureChromosome, featurePositionStart, sep = "_")) %>%
  dplyr::rename(ensemblID = featureName)
coloc.xqtl.mqtl.df = coloc.xqtl.mqtl.df %>%
  dplyr::left_join(xqtl.mqtl.minp.df %>% dplyr::select(ensemblID, geneSymbol), by=c("phenotype_id" = "ensemblID"))

rm(xqtl.haqtl.minp.df, xqtl.mqtl.minp.df)

coloc.sensneur.eqtl.df = readr::read_tsv(file.path(root, "coloc", "sensoryneuron", paste0(prefix, ".sens_neur.500k.1e+05.txt")))
# Use Ensembl to HGNC map to add in gene symbol
coloc.sensneur.eqtl.df %<>% dplyr::left_join(ensemblMap, by=c("phenotype_id" = "ensembl_gene_id")) %>%
  dplyr::filter(!duplicated(phenotype_id, symbol)) %>%
  dplyr::select(-symbol_type) %>%
  dplyr::rename(geneSymbol = symbol) %>%
  dplyr::mutate(dataset = "Sensory neuron eQTL", dataset_short = "sensneur_eqtl")

# For sensory neuron sQTLs, get "gene symbol" to use from the original lead SNP file
coloc.sensneur.sqtl.df = readr::read_tsv(file.path(root, "coloc", "sensoryneuron", paste0(prefix, ".sens_neur.sqtl.1e+05.txt")))
sensneurLeadSQTL.df = readr::read_tsv(file.path(root, "sensoryneurons", "GRCh38", "fastqtl.sqtl.permutations.10k.PCs.5.fdr0.1.txt"))
geneSummary = function(i) {
  if (grepl("[A-Za-z]", sensneurLeadSQTL.df[i,]$symbols)) {
    return(sensneurLeadSQTL.df[i,]$symbols)
  }
  return(sensneurLeadSQTL.df[i,]$ensembl_ids)
}
sensneurLeadSQTL.df$geneSymbol = sapply(1:nrow(sensneurLeadSQTL.df), FUN = geneSummary)
coloc.sensneur.sqtl.df %<>% dplyr::left_join(sensneurLeadSQTL.df %>% dplyr::select(phenotype_id, geneSymbol),
                                             by="phenotype_id") %>%
  dplyr::mutate(dataset = "Sensory neuron sQTL", dataset_short = "sensneur_sqtl")


###############################################################################
# Make a big table of all the colocs

# Write out a smaller version of the table, including each region-coloc as one row.
colocs.df = rbind(coloc.macrophage.eqtl.df,
                     coloc.macrophage.sqtl.df,
                     coloc.macrophage.reviseAnn.df,
                     coloc.mono.eqtl.df,
                     coloc.mono.h3k27ac.df,
                     coloc.mono.h3k4me1.df,
                     coloc.xqtl.eqtl.df,
                     coloc.xqtl.haqtl.df,
                     coloc.xqtl.mqtl.df,
                     coloc.sensneur.eqtl.df,
                     coloc.sensneur.sqtl.df
) %>% dplyr::rename(gwas_snp = gwas_lead)
colocs.df = rbind(coloc.mono.eqtl.df,
                  coloc.mono.h3k27ac.df,
                  coloc.mono.h3k4me1.df,
                  coloc.xqtl.eqtl.df,
                  coloc.xqtl.haqtl.df,
                  coloc.sensneur.eqtl.df,
                  coloc.sensneur.sqtl.df
) %>% dplyr::rename(gwas_snp = gwas_lead)

colocs.df$geneSymbol[is.na(colocs.df$geneSymbol)] = colocs.df$phenotype_id[is.na(colocs.df$geneSymbol)]

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
colocs.df$chrpos = paste0(colocs.df$gwas_chr, ".", colocs.df$gwas_pos)
colocs.df %<>% dplyr::left_join(signals.df, by="chrpos")
colocs.df %<>% dplyr::select(-chrpos)

# Remove any duplicate entries with the same gwas_lead and qtl_lead SNP,
# and arrange in order of decreasing H4 probability
colocs.df %<>% dplyr::arrange(qtl_pval) %>%
  dplyr::mutate(sortkey = paste(gwas_lead, qtl_lead, phenotype_id, sep=".")) %>%
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
write.table(annotated.df %>% dplyr::select(-gwas_lead) %>%
              dplyr::arrange(Chr, locus, finemap.p), file=output_file,
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

