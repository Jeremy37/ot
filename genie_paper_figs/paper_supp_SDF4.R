library(tidyverse)
library(cowplot)

setwd("/Users/jeremys/work/opentargets")

# Supplementary Figure showing HDR:WT effects of SDF4 edits
barplot_theme = theme_bw(8) +
  theme(axis.line = element_line(colour = "black", size=0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

replicates.grep.df = readr::read_tsv("experiment/transcribed/batch3/analysis/batch3.SDF4_exonic.grep.grep_analysis.replicate_stats.tsv") %>%
  filter(grepl("_vs_TG", name))

color_values = c(`cDNA`="firebrick1", `cDNA outlier`="orange1", `gDNA`="dodgerblue3", `gDNA outlier`="turquoise1")

p.CC = ggplot(replicates.grep.df %>% filter(name == "SDF4_CC_vs_TG"),
              aes(x=type, y=HDR_WT_ratio, col=type)) +
  geom_boxplot(alpha = 0.8, position="dodge") +
  geom_jitter(alpha = 0.7, width = 0.25) +
  theme_bw(10) + theme(legend.title=element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(values = color_values, guide=F) +
  ggtitle("Allele CC")

p.CG = ggplot(replicates.grep.df %>% filter(name == "SDF4_CG_haplotype_vs_TG"),
              aes(x=type, y=HDR_WT_ratio, col=type)) +
  geom_boxplot(alpha = 0.8, position="dodge") +
  geom_jitter(alpha = 0.7, width = 0.25) +
  theme_bw(10) + theme(legend.title=element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(values = color_values, guide=F) +
  ggtitle("Allele CG")

p.TC = ggplot(replicates.grep.df %>% filter(name == "SDF4_TC_vs_TG"),
              aes(x=type, y=HDR_WT_ratio, col=type)) +
  geom_boxplot(alpha = 0.8, position="dodge") +
  geom_jitter(alpha = 0.7, width = 0.25) +
  theme_bw(10) + theme(legend.title=element_blank(), axis.title.x = element_blank()) +
  scale_color_manual(values = color_values, guide=F) +
  coord_cartesian(ylim = c(0, 0.4)) +
  ggtitle("Allele TC")

pdf("plots/sdf4_qc_hdr.pdf", width=7, height=2.2)
cowplot::plot_grid(p.CC, p.CG, p.TC, nrow = 1)
dev.off()

