#!/usr/bin/env Rscript
library(tidyverse)

source("/Users/jeremys/work/opentargets/src/R/ppaFunctions.R")
source("/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/src/R/ppaFunctions.R")

args <- commandArgs(trailingOnly = TRUE)
regionsFile = args[1]
finemap_v3meta_ncausal1_file = args[2]
finemap_v3meta_ncausal2_file = args[3]
finemap_v3meta_ncausal3_file = args[4]
finemap_v3meta_ncausal4_file = args[5]
gcta_cond_file = args[6]
out = args[7]

window = 1000000

setwd("/Users/jeremys/work/opentargets/AD_finemap/")
regions_file = "/Users/jeremys/work/opentargets/AD_finemap/AD.compare_loci.finemap.tsv"
finemap_v3meta_ncausal1_file = "/Users/jeremys/work/opentargets/AD_finemap/compare_versions/v3.meta.finemap.ncausal_1.snp"
finemap_v3meta_ncausal2_file = "/Users/jeremys/work/opentargets/AD_finemap/compare_versions/v3.meta.finemap.ncausal_2.snp"
finemap_v3meta_ncausal3_file = "/Users/jeremys/work/opentargets/AD_finemap/compare_versions/v3.meta.finemap.ncausal_3.snp"
finemap_v3meta_ncausal4_file = "/Users/jeremys/work/opentargets/AD_finemap/compare_versions/v3.meta.finemap.ncausal_4.snp"
gcta_cond_file = "/Users/jeremys/work/opentargets/AD_finemap/gcta/output_1e-5/cond/merged_loci.cond.out.flt.tsv"
out = "/Users/jeremys/work/opentargets/AD_finemap/compare_versions/finemap.ncausal_comparison"

regions.df = readr::read_tsv(regions_file)

ViewDups = function(df, col) {
  df = as.data.frame(df)
  dupVals = df[duplicated(df[,col]), col]
  dups = df[df[,col] %in% dupVals,]
  View(dups[order(dups[,col]),])
}

# Add locus labels to the SNPs
addLocusToDF = function(df, regions.df) {
  df$locus = ""
  df$locus_ncausal = 0
  for (i in 1:nrow(regions.df)) {
    locus = regions.df$locus_name[i]
    chr = regions.df$Chr[i]
    start = regions.df$start[i] - window
    end = regions.df$stop[i] + window
    locus_ncausal = regions.df$n_snps[i]
    inlocus = (df$chr == chr & start <= df$pos & df$pos <= end)
    print(sprintf("Num SNPs in locus = %d", sum(inlocus)))
    if (sum(inlocus) > 0) {
      df[inlocus,]$locus = locus
      df[inlocus,]$locus_ncausal = locus_ncausal
    }
  }
  df
}


###############################################################################
finemap.ncausal_1.df = readr::read_tsv(finemap_v3meta_ncausal1_file) %>%
  select(snp=rsid, chr=chromosome, pos = position, prob1 = prob)
finemap.ncausal_2.df = readr::read_tsv(finemap_v3meta_ncausal2_file) %>%
  select(snp=rsid, prob2 = prob)
#  select(snp=rsid, chr=chromosome, pos = position, prob2 = prob)

finemap.df = finemap.ncausal_1.df %>%
  full_join(finemap.ncausal_2.df, by="snp")

finemap.df = addLocusToDF(finemap.df, regions.df) %>%
  arrange(locus, desc(prob2)) %>%
  filter(locus_ncausal >= 2)

# Plot the loci individually
plotCompareFinemapLocus = function(df, xlabel, ylabel) {
  topSnps1 = df %>% arrange(-prob1) %>% filter(prob1 >= 0.01 & prob2 >= 0.01) %>% .[1:5,] %>% na.omit()
  topSnps2 = df %>% arrange(-prob2) %>% filter(prob1 >= 0.01 & prob2 >= 0.01) %>% .[1:5,] %>% na.omit()
  topSnpsBoth = bind_rows(topSnps1, topSnps2 %>% filter(!snp %in% topSnps1$snp))
  topSnpsBoth$snp[is.na(topSnpsBoth$snp)] = topSnpsBoth$snp[is.na(topSnpsBoth$snp)]
  
  topSnps1 = df %>% arrange(-prob1) %>% filter(prob2 < 0.01) %>% .[1:5,] %>% na.omit()
  topSnps2 = df %>% arrange(-prob2) %>% filter(prob1 < 0.01) %>% .[1:5,] %>% na.omit()
  topSnpsEither = bind_rows(topSnps1, topSnps2 %>% filter(!snp %in% topSnps1$snp))
  topSnpsEither$snp = topSnpsEither$snp
  topSnpsEither$snp[is.na(topSnpsEither$snp)] = topSnpsEither$snp[is.na(topSnpsEither$snp)]
  xmax = max(topSnpsEither$prob1, na.rm=T)
  ymax = max(topSnpsEither$prob2, na.rm=T)
  
  locusName = df$locus[!is.na(df$locus)][1]
  p = ggplot(df, aes(x=prob1, y=prob2)) +
    geom_point(alpha=0.4) +
    theme_bw() +
    ggtitle(locusName) + xlab(xlabel) + ylab(ylabel)
  if (nrow(topSnpsBoth) > 0) {
    p = p + geom_text(mapping=aes(label=snp, hjust="left", vjust="center"), data=topSnpsBoth, size=2.5, nudge_x=0.007*xmax)
  }
  if (nrow(topSnpsEither) > 0) {
    p = p + geom_text(mapping=aes(label=snp, hjust="left", vjust="center"), data=topSnpsEither, size=2.5, nudge_x=0.007*xmax, angle=40)
  }
  p
}

finemap.cmp.df = finemap.df
finemap.cmp.df$prob1[is.na(finemap.cmp.df$prob1)] = 0
finemap.cmp.df$prob2[is.na(finemap.cmp.df$prob2)] = 0
finemap.cmp.df = finemap.cmp.df %>%
  filter(is.na(prob1) | is.na(prob2) | prob1 > 0.0001 | prob2 > 0.0001)

plot_list = by(finemap.cmp.df, finemap.cmp.df$locus, FUN=function(df) plotCompareFinemapLocus(df, xlabel="ncausal=1 SNP prob", ylabel="ncausal=2 SNP prob"))

pdf(file = paste0(out, ".2_vs_1.loci.pdf"), width=10, height=8)
plot_list
dev.off()


# Compare more causals at the TREM2 locus
finemap.ncausal_3.df = readr::read_tsv(finemap_v3meta_ncausal3_file) %>%
  select(snp=rsid, prob3 = prob)

finemap.ncausal_4.df = readr::read_tsv(finemap_v3meta_ncausal4_file) %>%
  select(snp=rsid, prob4 = prob)

finemap.df = finemap.ncausal_1.df %>%
  full_join(finemap.ncausal_2.df, by="snp") %>%
  full_join(finemap.ncausal_3.df, by="snp") %>%
  full_join(finemap.ncausal_4.df, by="snp")

finemap.df = addLocusToDF(finemap.df, regions.df) %>%
  arrange(locus, desc(prob2))

finemap.cmp.df = finemap.df %>%
  filter(locus_ncausal >= 3)
finemap.cmp.df$prob1[is.na(finemap.cmp.df$prob1)] = 0
finemap.cmp.df$prob2[is.na(finemap.cmp.df$prob2)] = 0
finemap.cmp.df$prob3[is.na(finemap.cmp.df$prob3)] = 0
finemap.cmp.df$prob4[is.na(finemap.cmp.df$prob4)] = 0
finemap.cmp.df = finemap.cmp.df %>%
  filter(is.na(prob1) | is.na(prob2) | is.na(prob3) | is.na(prob3) | prob1 > 0.0001 | prob2 > 0.0001 | prob3 > 0.0001 | prob4 > 0.0001)


# Ncausal 2 vs 1
plot1 = by(finemap.cmp.df, finemap.cmp.df$locus, FUN=function(df) plotCompareFinemapLocus(df, xlabel="ncausal=1 SNP prob", ylabel="ncausal=2 SNP prob"))

# Ncausal 3 vs 2
finemap.cmp.df.tmp = finemap.cmp.df %>%
  select(-prob1) %>%
  rename(prob1 = prob2, prob2 = prob3)
plot2 = by(finemap.cmp.df.tmp, finemap.cmp.df.tmp$locus, FUN=function(df) plotCompareFinemapLocus(df, xlabel="ncausal=2 SNP prob", ylabel="ncausal=3 SNP prob"))

# Ncausal 3 vs 1
finemap.cmp.df.tmp = finemap.cmp.df %>%
  select(-prob2) %>%
  rename(prob2 = prob3)
plot3 = by(finemap.cmp.df.tmp, finemap.cmp.df.tmp$locus, FUN=function(df) plotCompareFinemapLocus(df, xlabel="ncausal=1 SNP prob", ylabel="ncausal=3 SNP prob"))

# Ncausal 4 vs 3
finemap.cmp.df.tmp = finemap.cmp.df %>%
  select(-prob1, -prob2) %>%
  rename(prob1 = prob3, prob2 = prob4)
plot4 = by(finemap.cmp.df.tmp, finemap.cmp.df.tmp$locus, FUN=function(df) plotCompareFinemapLocus(df, xlabel="ncausal=3 SNP prob", ylabel="ncausal=4 SNP prob"))

# Ncausal 4 vs 2
finemap.cmp.df.tmp = finemap.cmp.df %>%
  select(-prob1) %>%
  rename(prob1 = prob2, prob2 = prob4)
plot5 = by(finemap.cmp.df.tmp, finemap.cmp.df.tmp$locus, FUN=function(df) plotCompareFinemapLocus(df, xlabel="ncausal=2 SNP prob", ylabel="ncausal=4 SNP prob"))


pdf(file = paste0(out, ".ncausal_3.loci.pdf"), width=10, height=8)
plot1
plot2
plot3
plot4
plot5
dev.off()


###############################################################################
# Compare FINEMAP results with GCTA conditional analyses

gcta.df = readr::read_tsv(gcta_cond_file) %>%
  rename(cond_lead_snp = cond_snp, pCond = pC, snp=SNP) %>%
  left_join(regions.df %>% select(locus=locus_name, n_snps), by = "locus")

gcta.df$signal = paste(gcta.df$locus, gcta.df$cond_lead_snp)
gcta.df$gcta_condProb = getPPAs(gcta.df %>% rename(F = freq) %>% arrange(signal),
                                segmentCol = "signal", pCol = "pCond", defaultN = 10000)

# For each locus, collapse the conditional probabilities for each SNP to
# select that with the max probability
gctaCond.simple.df = gcta.df %>%
  group_by(locus, snp) %>%
  summarise(gcta_maxCondProb = max(gcta_condProb),
            gcta_pCond = pCond[which.max(gcta_condProb)],
            gcta_bCond = bC[which.max(gcta_condProb)],
            gcta_bCondSE = bC_se[which.max(gcta_condProb)])

# Make a column that has the Finemap probabilities corresponding to the
# number of SNPs selected by GCTA for the locus
getFinemapProb = function(df) {
  nCausal = df$locus_ncausal[1]
  if (nCausal == 1) {
    df$prob1
  } else if (nCausal == 2) {
    df$prob2
  } else if (nCausal == 3) {
    df$prob3
  } else if (nCausal == 4) {
    df$prob4
  } else {
    rep(NA, nrow(df))
  }
}
finemap.df$finemap_prob_nc = unlist(by(finemap.df, finemap.df$locus, getFinemapProb))

finemap.cmp.df = finemap.df %>%
  left_join(gctaCond.simple.df, by=c("locus", "snp")) %>%
  select(-prob1, -prob2, -prob3, -prob4) %>%
  rename(prob1 = gcta_maxCondProb, prob2 = finemap_prob_nc) %>%
  filter(prob1 > 0.005 | prob2 > 0.005)

plot_finemap_vs_gcta = by(finemap.cmp.df, finemap.cmp.df$locus,
           FUN=function(df) plotCompareFinemapLocus(df, xlabel="GCTA conditional SNP prob", ylabel="Finemap SNP prob"))

pdf(file = paste0(out, ".finemap_vs_gcta.loci.all.pdf"), width=10, height=8)
plot_finemap_vs_gcta
dev.off()


finemap.nochange.df = finemap.cmp.df %>%
  filter(locus %in% c("ACE", "ADAM10", "APP-ADAMTS1", "BIN1", "EPHA1", "NCK2", "PTK2B-CLU"))
plot_finemap_vs_gcta_nochange = ggplot(finemap.nochange.df, aes(x=prob1, y=prob2)) +
  geom_point(alpha=0.4) +
  theme_bw(14) + facet_wrap(~locus, scales = "free") +
  xlab("GCTA conditional SNP prob") + ylab("Finemap SNP pro")

pdf(file = paste0(out, ".finemap_vs_gcta.loci.nochange.pdf"), width=10, height=8)
plot_finemap_vs_gcta_nochange
dev.off()


finemap.changed.df = finemap.cmp.df %>%
  filter(locus %in% c("ABCA7", "MS4A4A", "PLCG2", "TREM2"))
plot_finemap_vs_gcta_changed = ggplot(finemap.changed.df, aes(x=prob1, y=prob2)) +
  geom_point(alpha=0.4) +
  theme_bw(14) + facet_wrap(~locus, scales = "free") +
  xlab("GCTA conditional SNP prob") + ylab("Finemap SNP pro")

pdf(file = paste0(out, ".finemap_vs_gcta.loci.changed.pdf"), width=10, height=8)
plot_finemap_vs_gcta_changed
dev.off()

