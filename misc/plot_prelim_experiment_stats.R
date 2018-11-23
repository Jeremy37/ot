#!/usr/bin/env Rscript
library(tidyverse)
root="/Users/jeremys/work/opentargets/experiment/transcribed"
fname = file.path(root, "batch1.region_stats.edited.txt")
df = readr::read_tsv(fname) %>% dplyr::arrange(index) %>%
  dplyr::filter(index != 25)
df$short_name = factor(df$short_name, levels=df$short_name)
df$nReps = 3
df$nReps[df$batch_redo == T] = 9
df$nReps = factor(df$nReps)
p = ggplot(df, aes(x=short_name, y=effect, fill=nReps)) + geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=effect-1.96*effect_sd, ymax=effect+1.96*effect_sd), width=.35) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim=c(0, 1.6)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  scale_fill_manual(values=c("3"="grey", "9"="dodgerblue")) +
  ylab("HDR effect relative to WT") + xlab("Variant / gene")
  
p

pdf(file.path(root, "batch1.region_stats.plots.pdf", width=8, height=5))
p
dev.off()
