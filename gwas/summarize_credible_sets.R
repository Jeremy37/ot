#!/usr/bin/env Rscript
library(tidyverse)

# Here we make a summary of the number of credible set SNPs across loci,
# for both AD and PD
root = "/Users/jeremys/work/opentargets/gwas/"
ad.df = readr::read_tsv(file.path(root, "AD/finemap", "AD.credible_sets.set"))
pd.df = readr::read_tsv(file.path(root, "PD/finemap", "pd_gwas_gwax.combined.set"))

makeSummaryPlots = function(df, trait, credsetName) {
  summary.df = df %>% dplyr::filter(cs > 0) %>% dplyr::group_by(locus) %>% dplyr::summarise(num_variants = n())
  p = ggplot(summary.df, aes(x=num_variants)) + geom_histogram(binwidth = 1) + theme_bw() +
    ggtitle(paste(trait, "Number of loci with each credible set size:", credsetName))
  print(p)
  p = ggplot(summary.df, aes(x=num_variants)) + stat_ecdf(geom = "step") + theme_bw() +
    ggtitle(paste(trait, "CDF of number of credible set variants per locus:", credsetName))
  print(p)
}


makeSummaryPlots(ad.df %>% dplyr::mutate(cs = cs95), "AD", "95% credible set")
makeSummaryPlots(ad.df %>% dplyr::mutate(cs = cs99), "AD", "99% credible set")

makeSummaryPlots(pd.df %>% dplyr::mutate(cs = cs95), "PD", "95% credible set")
makeSummaryPlots(pd.df %>% dplyr::mutate(cs = cs99), "PD", "99% credible set")
