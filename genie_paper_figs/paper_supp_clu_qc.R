library(tidyverse)
library(cowplot)

setwd("/Users/jeremys/work/opentargets")


# First save a supplementary table with details of the stats at CLU
# clu.regions.grep.df = readr::read_tsv("experiment/transcribed/clu_reanalysis/analysis/clu_reanalysis.20bp_edit_window.grep_analysis.region_stats.tsv") %>%
#   mutate(name = gsub("_g1", "", name,),
#          name = gsub("_g2", "", name))
# clu.regions.dels.df = readr::read_tsv("experiment/transcribed/clu_reanalysis/analysis/clu_reanalysis.20bp_edit_window.region_stats.tsv") %>%
#   mutate(name = gsub("_g1", "", name),
#          name = gsub("_g2", "", name))

clu.regions.grep.df = readr::read_tsv("experiment/transcribed/clu_reanalysis/analysis/clu_reanalysis.20bp_edit_window.subsampled.grep_analysis.region_stats.tsv") %>%
   mutate(name = gsub("_[^_]+$", "", name))
clu.regions.dels.df = readr::read_tsv("experiment/transcribed/clu_reanalysis/analysis/clu_reanalysis.20bp_edit_window.subsampled.region_stats.tsv") %>%
  mutate(name = gsub("_[^_]+$", "", name))

clu.stats.df = clu.regions.grep.df %>%
  left_join(clu.regions.dels.df %>% select(name, del_fraction = del.rate, del.pval, del.effect, del.effect_confint_lo, del.effect_confint_hi), by="name") %>%
  select(name, hdr_fraction = hdr.rate, everything(), -hdr.effect_sd)

readr::write_tsv(clu.stats.df, "experiment/transcribed/clu_reanalysis/analysis/clu_reanalysis.20bp_edit_window.supp_table_stats.tsv")

# Supplementary Figure: make barplots of HDR:WT ratio and Del:WT ratio for all CLU SNPs
barplot_theme = theme_bw(8) +
  theme(axis.line = element_line(colour = "black", size=0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

#clu.replicates.grep.df = readr::read_tsv("experiment/transcribed/clu_reanalysis/analysis/clu_reanalysis.20bp_edit_window.grep_analysis.replicate_stats.tsv") %>%
clu.replicates.grep.df = readr::read_tsv("experiment/transcribed/clu_reanalysis/analysis/clu_reanalysis.20bp_edit_window.subsampled.grep_analysis.replicate_stats.tsv") %>%
  mutate(name = gsub("_[^_]+$", "", name)) %>%
  mutate(name = factor(as.character(name), levels = paste0("CLU_", seq(1, 11))))


color_values = c(`cDNA`="firebrick1", `cDNA outlier`="orange1", `gDNA`="dodgerblue3", `gDNA outlier`="turquoise1")

p.replicate_hdr = ggplot(clu.replicates.grep.df, aes(x=type, y=HDR_WT_ratio, col=type)) +
  geom_boxplot(alpha = 0.8, position="dodge", outlier.shape=NA) +
  geom_jitter(alpha = 0.7, width = 0.25) +
  theme_bw(10) + theme(legend.title=element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(values = color_values, guide=F) +
  facet_wrap(~name, nrow=3, scales = "free")

#clu.replicates.del.df = readr::read_tsv("experiment/transcribed/clu_reanalysis/analysis/clu_reanalysis.20bp_edit_window.replicate_stats.tsv") %>%
clu.replicates.del.df = readr::read_tsv("experiment/transcribed/clu_reanalysis/analysis/clu_reanalysis.20bp_edit_window.subsampled.replicate_stats.tsv") %>%
  mutate(name = gsub("_[^_]+$", "", name)) %>%
  mutate(name = factor(as.character(name), levels = paste0("CLU_", seq(1, 11))))

p.replicate_dels = ggplot(clu.replicates.del.df, aes(x=type, y=DEL_WT_ratio, col=type)) +
  geom_boxplot(alpha = 0.8, position="dodge", outlier.shape=NA) +
  geom_jitter(alpha = 0.7, width = 0.25) +
  theme_bw(10) + theme(legend.title=element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(values = color_values, guide=F) +
  facet_wrap(~name, nrow=3, scales = "free")
p.replicate_dels

pdf("plots/clu_qc_hdr.pdf", width=8, height=3.5)
p.replicate_hdr
dev.off()

pdf("plots/clu_qc_dels.pdf", width=8, height=3.5)
p.replicate_dels
dev.off()
