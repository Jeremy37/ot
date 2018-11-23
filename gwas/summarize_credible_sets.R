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


# UPDATE: Do it again for the v2 AD summary stats
root = "/Users/jeremys/work/opentargets"
ad.loci = readr::read_tsv(file.path(root, "AD_finemap", "AD.locusnames.txt"))
ad.df = readr::read_tsv(file.path(root, "AD_finemap", "annotated", "AD.credible_sets.set.tsv")) %>%
  dplyr::left_join(ad.loci, by="locus")

ad.df$is_ukbb_variant = is.na(ad.df$prob)
ad.flt = ad.df %>% dplyr::group_by(locus_name) %>%
  dplyr::filter(cs95 > 0)
summary.df = ad.flt %>% dplyr::summarise(meta_variants = sum(cs95 & !is_ukbb_variant),
                                         ukbb_variants = sum(is_ukbb_variant)) %>%
  dplyr::arrange(meta_variants)
sortedLoci = as.character(summary.df$locus_name)

summary.df = summary.df %>% tidyr::gather(key="variant_type", value="num_variants", -locus_name)
summary.df$locus_name = factor(as.character(summary.df$locus_name), levels=sortedLoci)
summary.df$variant_type = factor(as.character(summary.df$variant_type), levels=c("ukbb_variants", "meta_variants"))

pdf(file = file.path(root, "AD_finemap", "AD.credset_sizes.all.pdf"), width=8, height=3.4)
ggplot(summary.df %>% dplyr::filter(locus_name != "NME8"), aes(x=locus_name, y=num_variants, fill=variant_type)) +
  geom_bar(stat="identity") +
  theme_bw(15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  ylim(c(0, 180)) + scale_fill_manual(values = c("ukbb_variants"="coral2", "meta_variants"="dodgerblue"), guide=F) +
  xlab("Locus") + ylab("Variants in\n95% credible set")
dev.off()

pdf(file = file.path(root, "AD_finemap", "AD.credset_sizes.meta.pdf"), width=8, height=3.4)
ggplot(summary.df %>% dplyr::filter(locus_name != "NME8", variant_type == "meta_variants"), aes(x=locus_name, y=num_variants, fill=variant_type)) +
  geom_bar(stat="identity") +
  theme_bw(15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  ylim(c(0, 180)) + scale_fill_manual(values = c("ukbb_variants"="coral2", "meta_variants"="dodgerblue"), guide=F) +
  xlab("Locus") + ylab("Variants in\n95% credible set")
dev.off()

