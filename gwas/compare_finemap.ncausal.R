#!/usr/bin/env Rscript
library(tidyverse)

# Here we compare the causal probabilities across loci estimated by FINEMAP when
# it's allowed to select 2 causal variants at all loci, rather than the number
# of causals (1 or 2) determined based on GCTA analyses.

root = "/Users/jeremys/work/opentargets"
finemapdir = file.path(root, "gwas/AD/finemap")

credset.df = readr::read_tsv(file.path(finemapdir, "AD.credible_sets.set")) %>%
  dplyr::select(locus, SNP, everything(), -log10bf) %>%
  dplyr::rename(wtccc_prob=prob)
finemap.df = readr::read_tsv(file.path(finemapdir, "AD.finemap.snp")) %>%
  dplyr::rename(SNP=snp, finemap_prob=snp_prob) %>%
  dplyr::select(locus, SNP, finemap_prob)

loci.df = credset.df %>% inner_join(finemap.df %>% dplyr::select(-locus), by="SNP") %>%
  arrange(Chr, locus, -wtccc_prob)
assertthat::assert_that(nrow(loci.df) == nrow(credset.df))

pdf(file=file.path(finemapdir, "finemap_cmp_wtccc.pdf"), width = 6, height=4.5)
for (locusName in unique(loci.df$locus)) {
  locus.df = loci.df %>% filter(locus == locusName)
  titleStr = paste0(locusName, "    prob sum: ", sum(locus.df$finemap_prob))
  print(titleStr)

  p = ggplot(locus.df, aes(x=wtccc_prob, y=finemap_prob)) +
    geom_jitter(size=1.2, alpha=0.5, color="blue", width=max(locus.df$wtccc_prob) / 100) +
    geom_abline(color="red", alpha=0.5) +
    ggtitle(titleStr)
  print(p)
}
dev.off()

# Do the same with the FINEMAP runs using max_causal=2
finemap.df = readr::read_tsv(file.path(finemapdir, "AD.finemap.max_causal_2.snp")) %>%
  dplyr::rename(SNP=snp, finemap_prob=snp_prob) %>%
  dplyr::select(locus, SNP, finemap_prob)

loci.df = credset.df %>% inner_join(finemap.df %>% dplyr::select(-locus), by="SNP")
#assertthat::assert_that(nrow(loci.df) == nrow(credset.df))
lostSNP.df = credset.df %>% filter(!(SNP %in% loci.df$SNP))
max(lostSNP.df$wtccc_prob)

pdf(file=file.path(finemapdir, "finemap_cmp_wtccc.max_causal_2.pdf"), width = 6, height=4.5)
for (locusName in unique(loci.df$locus)) {
  locus.df = loci.df %>% filter(locus == locusName)
  titleStr = paste0(locusName, "    prob sum: ", sum(locus.df$finemap_prob))
  print(titleStr)
  
  p = ggplot(locus.df, aes(x=wtccc_prob, y=finemap_prob)) +
    geom_jitter(size=1.2, alpha=0.5, color="blue", width=max(locus.df$wtccc_prob) / 100) +
    geom_abline(color="red", alpha=0.5) +
    ggtitle(titleStr)
  print(p)
}
dev.off()


###############################################################################
# Now for Parkinson's
finemapdir = file.path(root, "gwas/PD/finemap")

finemap.df = readr::read_tsv(file.path(finemapdir, "PD.finemap.snp")) %>%
  dplyr::rename(SNP=snp, finemap_prob=snp_prob) %>%
  dplyr::select(locus, SNP, finemap_prob)

loci.df = credset.df %>% inner_join(finemap.df %>% dplyr::select(-locus), by="SNP")
#assertthat::assert_that(nrow(loci.df) == nrow(credset.df))
lostSNP.df = credset.df %>% filter(!(SNP %in% loci.df$SNP))
max(lostSNP.df$wtccc_prob)

pdf(file=file.path(finemapdir, "finemap_cmp_wtccc.pdf"), width = 6, height=4.5)
for (locusName in unique(loci.df$locus)) {
  locus.df = loci.df %>% filter(locus == locusName)
  titleStr = paste0(locusName, "    prob sum: ", sum(locus.df$finemap_prob))
  print(titleStr)
  
  p = ggplot(locus.df, aes(x=wtccc_prob, y=finemap_prob)) +
    geom_jitter(size=1.2, alpha=0.5, color="blue", width=max(locus.df$wtccc_prob) / 100) +
    geom_abline(color="red", alpha=0.5) +
    ggtitle(titleStr)
  print(p)
}
dev.off()



