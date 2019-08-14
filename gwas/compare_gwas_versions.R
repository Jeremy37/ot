#!/usr/bin/env Rscript
library(tidyverse)

source("/Users/jeremys/work/opentargets/src/R/ppaFunctions.R")
source("/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/src/R/ppaFunctions.R")

  
args <- commandArgs(trailingOnly = TRUE)
regionsFile = args[1]
gwas_v2meta_file = args[2]
gwas_v3meta_file = args[3]
gwas_kunkle_file = args[4]
gwas_lambert_file = args[5]
gwas_gwax_file = args[6]
finemap_v2meta_file = args[7]
finemap_v3meta_file = args[8]
finemap_v3meta_nc_file = args[9]
out = args[10]

window = 1000000

setwd("/Users/jeremys/work/opentargets/AD_finemap/")
regions_file = "/Users/jeremys/work/opentargets/AD_finemap/AD.compare_loci.tsv"
gwas_v2meta_file = "/Users/jeremys/work/opentargets/AD_finemap/compare_versions/v2.meta.loci.gz"
gwas_v3meta_file = "/Users/jeremys/work/opentargets/AD_finemap/compare_versions/v3.meta.loci.gz"
gwas_kunkle_file = "/Users/jeremys/work/opentargets/AD_finemap/compare_versions/Kunkle.loci.gz"
gwas_lambert_file = "/Users/jeremys/work/opentargets/AD_finemap/compare_versions/Lambert.loci.gz"
gwas_gwax_file = "/Users/jeremys/work/opentargets/AD_finemap/compare_versions/gwax.loci.gz"
finemap_v2meta_file = "/Users/jeremys/work/opentargets/AD_finemap/compare_versions/v2.meta.finemap.snp"
finemap_v3meta_file = "/Users/jeremys/work/opentargets/AD_finemap/compare_versions/v3.meta.finemap.snp"
finemap_v3meta_nc_file = "/Users/jeremys/work/opentargets/AD_finemap/compare_versions/v3.meta.finemap_nc.tsv"
out = "/Users/jeremys/work/opentargets/AD_finemap/compare_versions/result"

regions.df = readr::read_tsv(regions_file)

gwas.v2meta.df = readr::read_tsv(gwas_v2meta_file) %>%
  select(chr=CHR, pos=BP, meta_v2_snp=SNP, meta_v2_a1=A1, meta_v2_a2=A2, meta_v2_p=META_P, meta_v2_beta=META_BETA, gwax_british_p1=GWAX_BRITISH_P) %>%
  arrange(meta_v2_p) %>% filter(!duplicated(paste(chr, pos)))

gwas.v3meta.df = readr::read_tsv(gwas_v3meta_file) %>%
  select(chr=CHR, pos=BP, meta_v3_snp=SNP, meta_v3_a1=A1, meta_v3_a2=A2, meta_v3_p=META_P, meta_v3_beta=META_BETA, gwax_british_p2=GWAX_BRITISH_P)

gwas.gwax.df = readr::read_tsv(gwas_gwax_file) %>%
  select(chr=CHR, pos=BP, gwax_snp=SNP, gwax_a1=ALLELE1, gwax_a2=ALLELE0, gwax_p=P_BOLT_LMM, gwax_beta=BETA, AF=A1FREQ, everything()) %>%
  mutate(MAF=sapply(AF, FUN=function(x) min(x, 1-x))) %>%
  arrange(gwax_p) %>% filter(!duplicated(paste(chr, pos)), MAF > 0)

gwas.kunkle.df = readr::read_tsv(gwas_kunkle_file) %>%
  select(chr=Chromosome, pos=Position, kunkle_snp=MarkerName, kunkle_a1=Effect_allele, kunkle_a2=Non_Effect_allele, kunkle_p=Pvalue, kunkle_beta=Beta, everything()) %>%
  arrange(kunkle_p) %>% filter(!duplicated(paste(chr, pos)))

gwas.lambert.df = readr::read_tsv(gwas_lambert_file) %>%
  select(chr=Chr, pos=Pos, lambert_snp=RSid, lambert_a1=Eff_allele, lambert_p=pval, lambert_beta=beta) %>%
  mutate(lambert_a2 = NA) %>%
  arrange(lambert_p) %>% filter(!duplicated(paste(chr, pos)))

ViewDups = function(df, col) {
  df = as.data.frame(df)
  dupVals = df[duplicated(df[,col]), col]
  dups = df[df[,col] %in% dupVals,]
  View(dups[order(dups[,col]),])
}

sameLength = function(a, b) { nchar(a) == nchar(b) }

gwas.df = gwas.v2meta.df %>%
  full_join(gwas.v3meta.df, by=c("chr", "pos")) %>%
  full_join(gwas.gwax.df, by=c("chr", "pos")) %>%
  full_join(gwas.kunkle.df, by=c("chr", "pos")) %>%
  full_join(gwas.lambert.df, by=c("chr", "pos")) %>%
  filter(is.na(meta_v2_p) | is.na(meta_v3_p) | is.na(kunkle_p) | is.na(lambert_p) |
           meta_v2_p < 0.001 | meta_v3_p < 0.001 | gwax_p < 0.001 | kunkle_p < 0.001 | lambert_p < 0.001,
         !(is.na(meta_v2_p) & is.na(meta_v3_p) & is.na(kunkle_p) & is.na(lambert_p)))

isIndel = function(i) {
  !sameLength(gwas.df$gwax_a1[i], gwas.df$gwax_a2[i]) |
  !sameLength(gwas.df$kunkle_a1[i], gwas.df$kunkle_a2[i]) |
  !sameLength(gwas.df$lambert_a1[i], gwas.df$lambert_a2[i])
}
gwas.df$is_indel = sapply(1:nrow(gwas.df), FUN = isIndel)
sum(gwas.df$is_indel, na.rm=T)


# Check how many SNPs have missing data in different GWAS datasets
sum(is.na(gwas.df$MAF))
sum(is.na(gwas.df$gwax_p))
sum(is.na(gwas.df$meta_v2_p))
sum(is.na(gwas.df$meta_v3_p))
sum(is.na(gwas.df$kunkle_p))
sum(is.na(gwas.df$lambert_p))

# SNPs in the old meta, but not the new one
sum(!is.na(gwas.df$meta_v2_p) & (is.na(gwas.df$kunkle_p) & is.na(gwas.df$meta_v3_p)))
# SNPs in Lambert but not in the updated GWASes
sum(!is.na(gwas.df$lambert_p) & (is.na(gwas.df$kunkle_p) & is.na(gwas.df$meta_v3_p)))
# SNPs in our new Meta, but not in the previous GWASes
sum(!is.na(gwas.df$meta_v3_p) & (is.na(gwas.df$meta_v2_p) & is.na(gwas.df$lambert_p)))
# SNPs in Kunkle but not Lambert or our new meta-analysis. Why are so many
# Kunkle SNPs absent from our meta-analysis?
sum(!is.na(gwas.df$kunkle_p) & (is.na(gwas.df$lambert_p) & is.na(gwas.df$meta_v3_p)))
View(gwas.df[!is.na(gwas.df$kunkle_p) & (is.na(gwas.df$lambert_p) & is.na(gwas.df$meta_v3_p)),])
sum(!is.na(gwas.df$kunkle_p) & (is.na(gwas.df$lambert_p) & is.na(gwas.df$gwax_p)))
View(gwas.df[!is.na(gwas.df$kunkle_p) & (is.na(gwas.df$lambert_p) & is.na(gwas.df$gwax_p)),])
View(gwas.df[!is.na(gwas.df$kunkle_p) & (is.na(gwas.df$gwax_p)),])

# SNPs in Kunkle but not in our meta - despite being in the GWAX! I expected this to be zero
sum(!is.na(gwas.df$kunkle_p) & (is.na(gwas.df$meta_v3_p) & !is.na(gwas.df$gwax_p)))
View(gwas.df[!is.na(gwas.df$lambert_p) & (is.na(gwas.df$kunkle_p) & is.na(gwas.df$meta_v3_p)),])


# Add locus labels to the SNPs
addLocusToDF = function(df, regions.df) {
  df$locus = ""
  for (i in 1:nrow(regions.df)) {
    locus = regions.df$locus_name[i]
    chr = regions.df$Chr[i]
    start = regions.df$start[i] - window
    end = regions.df$stop[i] + window
    inlocus = (df$chr == chr & start <= df$pos & df$pos <= end)
    print(sprintf("Num SNPs in locus = %d", sum(inlocus)))
    df[inlocus,]$locus = locus
  }
  df
}

gwas.df = addLocusToDF(gwas.df, regions.df) %>%
  arrange(locus)

# Compute WTCCC-style posterior probabilities (assuming 1 causal variant)
# for each locus in each GWAS version.
gwas.df$meta_v2_ppa = getPPAs(gwas.df, "locus", "meta_v2_p")
gwas.df$meta_v3_ppa = getPPAs(gwas.df, "locus", "meta_v3_p")
gwas.df$kunkle_ppa = getPPAs(gwas.df, "locus", "kunkle_p")
gwas.df$lambert_ppa = getPPAs(gwas.df, "locus", "lambert_p")
gwas.df$gwax_ppa = getPPAs(gwas.df, "locus", "gwax_p")

# Add a column indicating whether the SNP has another at the same chr:pos.
gwas.df$chr_pos = paste(gwas.df$chr, gwas.df$pos, sep=":")
dupVals = gwas.df$chr_pos[duplicated(gwas.df$chr_pos)]
gwas.df$dup_chrpos = gwas.df$chr_pos %in% dupVals

# For each locus, get a table of the SNPs that have a reasonably high
# PPA in one of the studies, but are absent or have low PPA in our
# v3 meta-analysis.
gwas.df$missing = is.na(gwas.df$meta_v3_p) & (!is.na(gwas.df$lambert_p) | !is.na(gwas.df$meta_v2_p) | !is.na(gwas.df$gwax_p))

PPA_THRESHOLD = 0.05
gwas.df$missing_of_interest = is.na(gwas.df$meta_v3_p) &
  (  (!is.na(gwas.df$kunkle_ppa) & gwas.df$kunkle_ppa > PPA_THRESHOLD) |
     (!is.na(gwas.df$lambert_ppa) & gwas.df$lambert_ppa > PPA_THRESHOLD) |
     (!is.na(gwas.df$meta_v2_ppa) & gwas.df$meta_v2_ppa > PPA_THRESHOLD))

gwas.df$added = !is.na(gwas.df$meta_v3_ppa) & is.na(gwas.df$meta_v2_p)
gwas.df$added_of_interest = gwas.df$added & gwas.df$meta_v3_ppa > PPA_THRESHOLD

sum(gwas.df %>% filter(!is.na(meta_v3_ppa), is.na(meta_v2_ppa)) %>% .$meta_v3_ppa)
# Total ppa of new variants is 5.6

sum(gwas.df %>% filter(!is.na(meta_v3_ppa), is.na(lambert_ppa)) %>% .$meta_v3_ppa)
sum(gwas.df %>% filter(!is.na(meta_v3_ppa), is.na(meta_v2_ppa)) %>% .$meta_v3_ppa > 0.05)
sum(gwas.df$added_of_interest)
sum(gwas.df$missing_of_interest)

# How much of the new variant PPA is indels? Only 0.7
sum(gwas.df %>% filter(is_indel, !is.na(meta_v3_ppa), is.na(meta_v2_ppa)) %>% .$meta_v3_ppa)
sum(gwas.df %>% filter(is_indel) %>% .$meta_v3_ppa, na.rm=T)
View(gwas.df %>% select(-meta_v2_a1, -meta_v2_a2, -meta_v2_beta) %>% filter(added, meta_v3_ppa > 0.01))

# For plotting, any SNP which isn't in a study is given a P value of 1,
# so that it is still plotted and we can see that it was absent.
gwas.df_nona = gwas.df
gwas.df_nona$meta_v2_p[is.na(gwas.df$meta_v2_p)] = 1
gwas.df_nona$meta_v3_p[is.na(gwas.df$meta_v3_p)] = 1
gwas.df_nona$kunkle_p[is.na(gwas.df$kunkle_p)] = 1
gwas.df_nona$lambert_p[is.na(gwas.df$lambert_p)] = 1
gwas.df_nona$gwax_p[is.na(gwas.df$gwax_p)] = 1

pdf(file = paste0(out, ".meta_v3_vs_v2.overview.pdf"), width=15, height=12)
ggplot(gwas.df_nona, aes(x=-log10(meta_v2_p), y=-log10(meta_v3_p))) +
  geom_point(alpha=0.4) +
  theme_bw() + facet_wrap(~locus, scales = "free") +
  xlab("V2 AD meta -log10(p)") + ylab("V3 AD meta -log10(p)")
dev.off()

# Plot the loci individually
plotCompareGwasLocus = function(df, xlabel, ylabel) {
  topSnps1 = df %>% arrange(p1) %>% filter(p2 < 1) %>% .[1:5,]
  topSnps2 = df %>% arrange(p2) %>% filter(p1 < 1) %>% .[1:5,]
  topSnpsBoth = bind_rows(topSnps1, topSnps2 %>% filter(!snp2 %in% topSnps1$snp2))
  topSnpsBoth$snp = topSnpsBoth$snp1
  topSnpsBoth$snp[is.na(topSnpsBoth$snp)] = topSnpsBoth$snp2[is.na(topSnpsBoth$snp)]
  xmax = max(-log10(topSnpsBoth$p1), na.rm=T)
  ymax = max(-log10(topSnpsBoth$p2), na.rm=T)
  
  topSnps1 = df %>% arrange(p1) %>% filter(p2 == 1) %>% .[1:5,]
  topSnps2 = df %>% arrange(p2) %>% filter(p1 == 1) %>% .[1:5,]
  topSnpsEither = bind_rows(topSnps1, topSnps2 %>% filter(!snp2 %in% topSnps1$snp2))
  topSnpsEither$snp = topSnpsEither$snp1
  topSnpsEither$snp[is.na(topSnpsEither$snp)] = topSnpsEither$snp2[is.na(topSnpsEither$snp)]
  
  minxySide = min(xmax, ymax)
  ggplot(df, aes(x=-log10(p1), y=-log10(p2))) +
    geom_point(alpha=0.4) +
    theme_bw() +
    ggtitle(df$locus[1]) + xlab(xlabel) + ylab(ylabel) +
    geom_text(mapping=aes(label=snp, hjust="left", vjust="center"), data=topSnpsBoth, size=2.5, nudge_x=0.007*xmax) +
    geom_text(mapping=aes(label=snp, hjust="left", vjust="center"), data=topSnpsEither, size=2.5, nudge_x=0.007*xmax, angle=40) +
    annotate("segment", x = 0, xend = minxySide, y = 0, yend = minxySide, colour = "blue", alpha = 0.4) +
    annotate("segment", x = 0, xend = minxySide, y = -log10(5e-8), yend = -log10(5e-8), colour = "red", alpha = 0.3, linetype="dashed")
}


pdf(file = paste0(out, ".meta_v3_vs_v2.loci.pdf"), width=10, height=8)
plot_list = by(gwas.df_nona %>% rename(p1=meta_v2_p, p2=meta_v3_p, snp1=meta_v2_snp, snp2=meta_v3_snp),
               gwas.df_nona$locus, FUN=function(df) plotCompareGwasLocus(df, xlabel="-log10(p) meta_v2", ylabel="-log10(p) meta_v3"))
print(plot_list)
dev.off()

pdf(file = paste0(out, ".meta_v3_vs_v2.loci.little_change.pdf"), width=10, height=8)
gwas.df.tmp = gwas.df_nona %>% rename(p1=meta_v2_p, p2=meta_v3_p, snp1=meta_v2_snp, snp2=meta_v3_snp) %>%
  filter(locus %in% c("ADAM10", "BIN1", "CASS4", "PILRA"))
ggplot(gwas.df.tmp, aes(x=-log10(p1), y=-log10(p2))) +
  geom_point(alpha=0.4) +
  theme_bw(14) + facet_wrap(~locus, scales = "free") +
  xlab("V2 AD meta -log10(p)") + ylab("V3 AD meta -log10(p)")
dev.off()

pdf(file = paste0(out, ".meta_v3_vs_v2.loci.new_variants.pdf"), width=10, height=8)
gwas.df.tmp = gwas.df_nona %>% rename(p1=meta_v2_p, p2=meta_v3_p, snp1=meta_v2_snp, snp2=meta_v3_snp) %>%
  filter(locus %in% c("PTK2B-CLU", "SCIMP", "TMEM163", "MS4A4A"))
ggplot(gwas.df.tmp, aes(x=-log10(p1), y=-log10(p2))) +
  geom_point(alpha=0.4) +
  theme_bw(14) + facet_wrap(~locus, scales = "free") +
  xlab("V2 AD meta -log10(p)") + ylab("V3 AD meta -log10(p)")
dev.off()

pdf(file = paste0(out, ".meta_v3_vs_v2.loci.new_loci.pdf"), width=10, height=4.5)
gwas.df.tmp = gwas.df_nona %>% rename(p1=meta_v2_p, p2=meta_v3_p, snp1=meta_v2_snp, snp2=meta_v3_snp) %>%
  filter(locus %in% c("NCK2", "TSPAN14"))
ggplot(gwas.df.tmp, aes(x=-log10(p1), y=-log10(p2))) +
  geom_point(alpha=0.4) +
  theme_bw(14) + facet_wrap(~locus, scales = "free") +
  xlab("V2 AD meta -log10(p)") + ylab("V3 AD meta -log10(p)")
dev.off()



pdf(file = paste0(out, ".kunkle_vs_lambert.loci.pdf"), width=10, height=8)
plot_list = by(gwas.df_nona %>% rename(p1=lambert_p, p2=kunkle_p, snp1=lambert_snp, snp2=kunkle_snp),
               gwas.df_nona$locus, FUN=function(df) plotCompareGwasLocus(df, xlabel="-log10(p) Lambert", ylabel="-log10(p) Kunkle"))
print(plot_list)
dev.off()


pdf(file = paste0(out, ".meta_v3_vs_kunkle.loci.pdf"), width=10, height=8)
plot_list = by(gwas.df_nona %>% rename(p1=kunkle_p, p2=meta_v3_p, snp1=kunkle_snp, snp2=meta_v3_snp),
               gwas.df_nona$locus, FUN=function(df) plotCompareGwasLocus(df, xlabel="-log10(p) Kunkle", ylabel="-log10(p) meta_v3"))
print(plot_list)
dev.off()


###############################################################################
finemap1.df = readr::read_tsv(finemap_v2meta_file) %>%
  rename(allele1=noneff_allele, allele2=eff_allele, prob1=finemap_prob, snp=rsID)
finemap2.df = readr::read_tsv(finemap_v3meta_file) %>%
  rename(prob2 = prob, snp=rsid, chr=chromosome, pos = position)

finemap.df = finemap1.df %>%
  full_join(finemap2.df, by=c("chr", "pos"))

finemap.df = addLocusToDF(finemap.df, regions.df) %>%
  arrange(locus, desc(prob2)) %>%
  filter(locus != "")

# Plot the loci individually
plotCompareFinemapLocus = function(df, xlabel, ylabel) {
  topSnps1 = df %>% arrange(-prob1) %>% filter(prob2 > 0) %>% .[1:5,] %>% na.omit()
  topSnps2 = df %>% arrange(-prob2) %>% filter(prob1 > 0) %>% .[1:5,] %>% na.omit()
  topSnpsBoth = bind_rows(topSnps1, topSnps2 %>% filter(!snp.y %in% topSnps1$snp.y))
  topSnpsBoth$snp = topSnpsBoth$snp.x
  topSnpsBoth$snp[is.na(topSnpsBoth$snp)] = topSnpsBoth$snp.y[is.na(topSnpsBoth$snp)]

  topSnps1 = df %>% arrange(-prob1) %>% filter(prob2 == 0) %>% .[1:5,]
  topSnps2 = df %>% arrange(-prob2) %>% filter(prob1 == 0) %>% .[1:5,]
  topSnpsEither = bind_rows(topSnps1, topSnps2 %>% filter(!snp.y %in% topSnps1$snp.y))
  topSnpsEither$snp = topSnpsEither$snp.x
  topSnpsEither$snp[is.na(topSnpsEither$snp)] = topSnpsEither$snp.y[is.na(topSnpsEither$snp)]
  xmax = max(c(topSnpsEither$prob1, topSnpsBoth$prob1), na.rm=T)
  ymax = max(c(topSnpsEither$prob2, topSnpsBoth$prob2), na.rm=T)
  
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

finemap.cmp.df = finemap.df %>%
  filter(!is.na(prob1) | !is.na(prob2))
finemap.cmp.df$prob1[is.na(finemap.cmp.df$prob1)] = 0
finemap.cmp.df$prob2[is.na(finemap.cmp.df$prob2)] = 0
finemap.cmp.df = finemap.cmp.df %>%
  filter(is.na(prob1) | is.na(prob2) | prob1 > 0.0001 | prob2 > 0.0001)

plot_list = by(finemap.cmp.df, finemap.cmp.df$locus, FUN=function(df) plotCompareFinemapLocus(df, xlabel="V2 meta SNP prob", ylabel="V3 meta SNP prob"))

pdf(file = paste0(out, ".finemap.loci.pdf"), width=10, height=8)
plot_list
dev.off()

# Make a plot of credible set sizes at each locus
finemap.cmp.df = finemap.cmp.df %>% arrange(locus, desc(prob1))
finemap.cmp.df$v2_cum_ppa = unlist(by(finemap.cmp.df, finemap.cmp.df$locus, FUN=function(df) {cumsum(df$prob1)}))
finemap.cmp.df = finemap.cmp.df %>% arrange(locus, desc(prob2))
finemap.cmp.df$v3_cum_ppa = unlist(by(finemap.cmp.df, finemap.cmp.df$locus, FUN=function(df) {cumsum(df$prob2)}))
locus_max_ppa = finemap.cmp.df %>% group_by(locus) %>%
  summarise(v2_prob_sum = max(v2_cum_ppa, na.rm=T), v3_prob_sum = max(v3_cum_ppa, na.rm=T))

finemap.onlyv2.df = finemap.cmp.df %>%
  filter(!is.na(prob1) & is.na(prob2))
finemap.onlyv3.df = finemap.cmp.df %>%
  filter(is.na(prob1) & !is.na(prob2))
finemap.inv2.df = finemap.cmp.df %>%
  filter(!is.na(prob1))
finemap.inv3.df = finemap.cmp.df %>%
  filter(!is.na(prob2))



finemap.df = finemap.df %>% arrange(desc(prob2), desc(prob1))
finemap.df = finemap.df %>% filter(!duplicated(paste(chr, pos)))
# Join the finemap details with the gwas summary table
gwas.df = gwas.df %>%
  left_join(finemap.df %>% select(chr, pos, meta_v2_finemap=prob1, meta_v3_finemap=prob2),
            by=c("chr", "pos"))

gwas.df$finemap_diff = gwas.df$meta_v3_finemap - gwas.df$meta_v2_finemap
gwas.df$divergent_finemap_down = !is.na(gwas.df$meta_v3_finemap) &
  !is.na(gwas.df$meta_v2_finemap) & gwas.df$meta_v2_finemap > 0.05 &
  (log10(gwas.df$meta_v3_finemap) - log10(gwas.df$meta_v2_finemap)) < -1
gwas.df$divergent_finemap_up = !is.na(gwas.df$meta_v3_finemap) &
  !is.na(gwas.df$meta_v2_finemap) & gwas.df$meta_v3_finemap > 0.05 &
  (log10(gwas.df$meta_v3_finemap) - log10(gwas.df$meta_v2_finemap)) > 1

# Add in the "finemap_nc" column from our annotated spreadsheet, representing
# the finemap probs from the expected number of causal variants at each locus.
finemap.nc.df = readr::read_tsv(finemap_v3meta_nc_file) %>%
  rename(meta_v3_snp = snp, meta_v3_finemap_nc = finemap_prob_nc)

gwas.df = gwas.df %>% 
  left_join(finemap.nc.df, by="meta_v3_snp") %>%
  select(locus, chr, pos, meta_v3_snp, meta_v3_a1, meta_v3_a2,
         meta_v3_finemap_nc, meta_v3_finemap, meta_v2_finemap,
         finemap_diff, divergent_finemap_down, divergent_finemap_up,
         meta_v3_p, meta_v2_p, kunkle_p, lambert_p, gwax_p,
         meta_v3_ppa, meta_v2_ppa, kunkle_ppa, lambert_ppa, gwax_ppa,
         meta_v2_snp, kunkle_snp, lambert_snp, gwax_snp,
         meta_v2_a1, meta_v2_a2, kunkle_a1, kunkle_a2, lambert_a1, lambert_a2, gwax_a1, gwax_a2,
         missing_of_interest, added_of_interest)

write.table(gwas.df,
            file = paste0(out, ".meta_v3_merged.tsv"),
            quote=F, sep="\t", row.names=F, col.names=T, na="")

write.table(gwas.df %>% filter(divergent_finemap_down | divergent_finemap_up | missing_of_interest | added_of_interest ),
            file = paste0(out, ".meta_v3_merged.snps_of_interest.tsv"),
            quote=F, sep="\t", row.names=F, col.names=T, na="")
