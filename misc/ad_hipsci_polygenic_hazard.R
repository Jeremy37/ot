#!/usr/bin/env Rscript
library(tidyverse)
options(stringsAsFactors=F)

args = commandArgs(trailingOnly = TRUE)
hipsciVCF.file = args[1]
phsSNPDetails.file = args[2]
outputPath = args[3]

# setwd("/Users/jeremys/work/opentargets/polygenic_hazard/")
# hipsciVCF.file = "hipsci.merged.phs.snps.vcf"
# phsSNPDetails.file = "phs.snp.details.txt"
# outputPath = "hipsci"

hipsci.df = readr::read_tsv(hipsciVCF.file)
phs.snp.df = readr::read_tsv(phsSNPDetails.file)

getGTString = function(gtStr) {
  sapply(strsplit(gtStr, ":", fixed=T), function(l) l[[1]])
}

hipsci.gather.df = hipsci.df %>%
  dplyr::select(ID, starts_with("HPSI")) %>%
  tidyr::gather(sampleID, gtCol, -ID) %>%
  dplyr::mutate(donor = sapply(strsplit(sampleID, "[-_]"), function(l) l[[2]])) %>%
  dplyr::mutate(gtStr = getGTString(gtCol)) %>%
  dplyr::mutate(allele1 = strtoi(sapply(strsplit(gtStr, "|", fixed=T), function(l) l[[1]]))) %>%
  dplyr::mutate(allele2 = strtoi(sapply(strsplit(gtStr, "|", fixed=T), function(l) l[[2]]))) %>%
  dplyr::mutate(gtCount = allele1 + allele2)

getDonorFromSampleID = function(sampleID) {
  sapply(strsplit(sampleID, "[-_]"), function(l) l[[2]])
}

# Get APOE genotypes
getApoeAllele = function(ID, alleles) {
  if (ID[1] == "rs7412") {
    ID = rev(ID)
    allele = rev(alleles)
  }
  # alleles[1] == 0 --> rs429358-T
  # alleles[1] == 1 --> rs429358-C
  # 
  # alleles[2] == 0 --> rs7412-C
  # alleles[2] == 1 --> rs7412-T
  if (alleles[1] == 0 & alleles[2] == 1) {
    "e2"
  } else if (alleles[1] == 0 & alleles[2] == 0) {
    "e3"
  } else if (alleles[1] == 1 & alleles[2] == 0) {
    "e4"
  } else {
    "E5"
    #stop("Unexpected APOE allele! rs7412-T, rs429358-C")
  }
}


hipsci.apoe.df = hipsci.gather.df %>%
  dplyr::group_by(sampleID) %>%
  dplyr::filter(ID %in% c("rs7412", "rs429358")) %>%
  dplyr::summarise(apoeAllele1 = getApoeAllele(ID, allele1),
                   apoeAllele2 = getApoeAllele(ID, allele2))

table(hipsci.apoe.df$apoeAllele1)
table(hipsci.apoe.df$apoeAllele2)

# Check allele frequency of APOE SNPs for sanity
rs7412_af = hipsci.gather.df %>%
  dplyr::filter(ID %in% c("rs7412")) %>%
  .$gtCount %>% sum() / nrow(hipsci.apoe.df) / 2
print(sprintf("rs7412 AF = %.3g", rs7412_af)) # Should be about 8%

rs429358_af = hipsci.gather.df %>%
  dplyr::filter(ID %in% c("rs429358")) %>%
  .$gtCount %>% sum() / nrow(hipsci.apoe.df) / 2
print(sprintf("rs429358 AF = %.3g", rs429358_af)) # Should be about 15%

# Compute hazard ratios. First do it for APOE, then for other SNPs,
# and finally combine these together.

apoe.snp.df = data.frame(ID = c("e2", "e3", "e4"), hr = c(exp(-0.47), 1, exp(1.03)))
hipsci.apoe.hr = hipsci.apoe.df %>%
  dplyr::left_join(apoe.snp.df %>% dplyr::rename(allele1HR = hr), by=c("apoeAllele1" = "ID")) %>%
  dplyr::left_join(apoe.snp.df %>% dplyr::rename(allele2HR = hr), by=c("apoeAllele2" = "ID")) %>%
  dplyr::mutate(apoeHR = allele1HR * allele2HR)

hipsci.all.hr = hipsci.gather.df %>%
  dplyr::left_join(phs.snp.df %>% dplyr::select(ID, log_hr), by="ID") %>%
  na.omit() %>%
  dplyr::group_by(sampleID) %>%
  dplyr::summarise(polygenicHR = exp(sum(log_hr * gtCount))) %>%
  dplyr::left_join(hipsci.apoe.hr %>% dplyr::select(sampleID, apoeAllele1, apoeAllele2, apoeHR), by="sampleID") %>%
  dplyr::mutate(overallHR = polygenicHR * apoeHR) %>%
  dplyr::mutate(donor = getDonorFromSampleID(sampleID)) %>%
  dplyr::arrange(-overallHR)


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

donor.hr = hipsci.all.hr %>%
  dplyr::group_by(donor) %>%
  dplyr::summarise(polygenicHR = mean(polygenicHR),
                   apoeAllele1 = first(apoeAllele1),
                   apoeAllele2 = first(apoeAllele2),
                   apoeHR = getmode(apoeHR),
                   overallHR = mean(overallHR)) %>%
  dplyr::arrange(-overallHR)

donor.hr$polygenicHR_quantile = ecdf(donor.hr$polygenicHR)(donor.hr$polygenicHR)
donor.hr$overallHR_quantile = ecdf(donor.hr$overallHR)(donor.hr$overallHR)
hipsci.all.hr$polygenicHR_quantile = ecdf(donor.hr$polygenicHR)(hipsci.all.hr$polygenicHR)
hipsci.all.hr$overallHR_quantile = ecdf(donor.hr$overallHR)(hipsci.all.hr$overallHR)

#ggplot(donor.hr, aes(x=polygenicHR)) + geom_histogram(binwidth = 0.05)
#ggplot(donor.hr, aes(x=overallHR)) + geom_histogram(binwidth = 0.05)


write.table(hipsci.all.hr, file=paste0(outputPath, ".AD_polygenic_hazard.cell_lines.txt"),
            sep="\t", quote=F, col.names=T, row.names=F)

write.table(donor.hr, file=paste0(outputPath, ".AD_polygenic_hazard.donors.txt"),
            sep="\t", quote=F, col.names=T, row.names=F)

