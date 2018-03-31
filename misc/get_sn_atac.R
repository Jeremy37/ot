#!/usr/bin/env Rscript
# This script is to get candidate causal ATAC QTL SNPs, which can then be
# used for Sarah's CRISPR assays for ATAC-altering variants
library(tidyverse)
#setwd("/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/ATAC/")

args <- commandArgs(trailingOnly = TRUE)
leadSnpFile = args[1]  #rasqual.1k.leadSNPs.fdr0.1.ann.txt
allSnpFile = args[2]   #rasqual.1k.pthresh0.01.ppathresh.0.0001.txt
snPeakFile = args[3]   #sensoryneuron_consensus_atac_peaks.GRCh38.bed
ipscPeakFile = args[4] #atac_ipsc_peaks.narrowPeak
leadSnpFile = "rasqual.1k.leadSNPs.fdr0.1.ann.txt"
allSnpFile = "rasqual.1k.pthresh0.01.ppathresh.0.0001.txt"
snPeakFile = "sensoryneuron_consensus_atac_peaks.GRCh38.txt"
ipscPeakFile = "atac_ipsc_peaks.narrowPeak"

leadsnp.df = readr::read_tsv(leadSnpFile, col_types="ccciccddddddddiiidddcc")
#allsnp.df = readr::read_tsv(allSnpFile) # Doesn't work... not sure why, "Error in make.names(x) : invalid multibyte string 1"
allsnp.df = read.delim(allSnpFile) %>% dplyr::rename(peak = gene)
snPeaks.df = readr::read_tsv(snPeakFile, col_types="cccii") %>%
  dplyr::rename(peak=gene_id, chr=chromosome_name, sn_peak_start=exon_starts, sn_peak_end=exon_ends)
ipscPeaks.df = readr::read_tsv(ipscPeakFile, col_types="ciicicdddi", col_names=c("chr", "start", "end", "name", "score", "unk", "fc", "log10pval", "log10qval", "summitpos"))
ipscPeaks.df$chr = gsub("^chr", "", ipscPeaks.df$chr)

snp.df = leadsnp.df %>% dplyr::select(peak, FDR, geneid, symbol, pvalue) %>%
  dplyr::rename(leadSnpPval = pvalue) %>%
  dplyr::left_join(allsnp.df)

# Identify which iPSC peaks have at least one candidate causal SNP in our list,
# and then subset our SNP table to include all SNPs from those peaks where there
# was an overlap with an iPSC peak
candidate.snp.df = snp.df %>% dplyr::filter(PPA > 0.25)
candidate.snp.gr = GRanges(seqnames=candidate.snp.df$chr, 
                           ranges=IRanges(start=candidate.snp.df$pos, end=candidate.snp.df$pos),
                           strand = NA,
                           candidate.snp.df)
ipsc.peaks.gr = GRanges(seqnames=ipscPeaks.df$chr, IRanges(start=ipscPeaks.df$start, ipscPeaks.df$end), strand=NA,
                        ipscPeaks.df %>% dplyr::select(-chr, -start, -end))
overlaps = as.data.frame(findOverlaps(candidate.snp.gr, ipsc.peaks.gr)) %>% dplyr::filter(!duplicated(queryHits))
peak.overlaps.df = data.frame(sn_peak = candidate.snp.df[overlaps$queryHits,]$peak,
                              ipsc_peak = ipscPeaks.df[overlaps$subjectHits,]$name)

#overlapping.snp.df = snp.df %>% dplyr::filter(peak %in% peak.overlaps) %>%
#  dplyr::arrange(leadSnpPval)

overlapping.snp.df = peak.overlaps.df %>%
  dplyr::left_join(ipscPeaks.df %>% dplyr::select(ipsc_peak=name, ipsc_peak_start=start, ipsc_peak_end=end)) %>%
  dplyr::left_join(snPeaks.df %>% dplyr::select(sn_peak=peak, sn_peak_start, sn_peak_end)) %>%
  dplyr::left_join(snp.df, by=c("sn_peak" = "peak")) %>%
  dplyr::arrange(leadSnpPval, -PPA) %>%
  dplyr::filter(FDR < 0.01, !duplicated(snpid)) %>%
  dplyr::select(sn_peak, ipsc_peak, chr, ipsc_peak_start, ipsc_peak_end, ipsc_peak_end, sn_peak_start, sn_peak_end,
                geneid, symbol, leadSnpPval, FDR, snpid, chr, pos, af, imputation_quality, effect_size, fsnps, rsnps, fsnp_genotype_corr, rsnp_genotype_corr, pvalue, PPA)

readr::write_tsv(overlapping.snp.df, "sensoryneurons.caqtls.in_ipsc_peaks.txt", col_names = T)
