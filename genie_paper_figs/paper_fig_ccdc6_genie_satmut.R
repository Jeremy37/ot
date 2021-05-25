library(tidyverse)
library(cowplot)

setwd("/Users/jeremys/work/opentargets")

getSignificanceStr = function(pval) {
  if (is.na(pval) | pval >= 0.01) { "p >= 0.01" }
  else if (pval < 0.001) { "p < 0.001" }
  else { "p < 0.01" }
}

barplot_theme = theme_bw(12) +
  theme(axis.line.x = element_line(colour = "black", size=0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

satmut.grep.df = readr::read_tsv("experiment/transcribed/CCDC6_satmut_genie/analysis/CCDC6_satmut_genie.HDR.by_batch.region_stats.tsv")
satmut.dels.df = readr::read_tsv("experiment/transcribed/CCDC6_satmut_genie/analysis/CCDC6_satmut_genie.HDR.by_batch.region_stats.tsv")

edit_groups.df = readr::read_tsv("experiment/transcribed/CCDC6_satmut_genie/ccdc6.satmut.edit_groups.tsv")
satmut.grep.df = satmut.grep.df %>%
  left_join(edit_groups.df, by=c("name"="edit"))

satmut.grep.df$hdr_significance = factor(sapply(satmut.grep.df$hdr.pval, FUN = getSignificanceStr), levels=c("p >= 0.01", "p < 0.01", "p < 0.001"))

p.barplot = ggplot(satmut.grep.df, aes(x=display_name, y=hdr.effect, fill = hdr_significance)) +
  geom_bar(aes(fill = hdr_significance), stat = "identity", position = "dodge", width=0.8) + 
  geom_errorbar(aes(ymin = hdr.effect_confint_lo, ymax = hdr.effect_confint_hi),
                position = "dodge", width = 0.5, col = "grey20") +
  geom_hline(yintercept = 1, col = "red", size=0.1) +
  scale_fill_manual(values=c(`p >= 0.01`="grey70", `p < 0.01`="cornflowerblue", `p < 0.001`="red"), guide=F) +
  barplot_theme + ylab("Effect size") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
p.barplot

pdf("plots/ccdc6_satmut_barplot.pdf", width = 6.5, height = 2.6)
p.barplot
dev.off()

satmut.grep.df$barwidth = c(1,2,rep(3, 31))/3
fill_vals = c(rep("cornflowerblue", 33), "grey70", "cornflowerblue", "red")
names(fill_vals) = c(satmut.grep.df$name, "p >= 0.01", "p < 0.01", "p < 0.001")
p.barplot_grouped = ggplot(satmut.grep.df, aes(x=edit_group, y=hdr.effect, fill = name)) +
  geom_bar(aes(group = name, fill = hdr_significance), stat = "identity", position = position_dodge(width = 1), width=0.8) +
  geom_errorbar(aes(ymin = hdr.effect_confint_lo, ymax = hdr.effect_confint_hi),
                position = position_dodge(width = 1), width = 0.8, col = "grey20") +
  geom_hline(yintercept = 1, col = "red", size=0.1) +
  scale_fill_manual(values=fill_vals, guide=F) +
  barplot_theme + ylab("Effect size") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
#scale_fill_manual(values=c(`p >= 0.01`="grey70", `p < 0.01`="cornflowerblue", `p < 0.001`="red"), guide=F) +
  
pdf("plots/ccdc6_satmut_barplot_grouped.pdf", width = 6.5, height = 2.1)
p.barplot_grouped
dev.off()

satmut.grep.df = satmut.grep.df %>% arrange(display_name)
satmut.grep.df$barpos = c(0, 1.2, 2.4, 3.1, 3.8, 5, 5.7, 6.4, 7.6, 8.3, 9, 10.2, 10.9, 11.6, 12.8, 13.5, 14.2, 15.4, 16.6, 17.3, 18, 19.2, 19.9, 20.6, 21.8, 23, 23.7, 24.4, 25.6, 26.3, 27, 28.2, 29.4)
p.barplot_grouped2 = ggplot(satmut.grep.df, aes(x=barpos, y=hdr.effect)) +
  geom_bar(aes(group = name, fill = hdr_significance), stat = "identity", position = position_dodge(width = 1), width=0.6) +
  geom_errorbar(aes(ymin = hdr.effect_confint_lo, ymax = hdr.effect_confint_hi),
                position = position_dodge(width = 1), width = 0.5, col = "grey20") +
  geom_hline(yintercept = 1, col = "red", size=0.1) +
  scale_fill_manual(values=fill_vals, guide=F) +
  barplot_theme + ylab("Effect size") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) +
  scale_x_reverse()
p.barplot_grouped2
pdf("plots/ccdc6_satmut_barplot_grouped2.pdf", width = 6.5, height = 1.75)
p.barplot_grouped2
dev.off()


p.barplot_grouped3 = ggplot(satmut.grep.df, aes(x=barpos, y=hdr.effect)) +
  geom_hline(yintercept = 1, col = "red", size=0.1) +
  geom_errorbar(aes(ymin = hdr.effect_confint_lo, ymax = hdr.effect_confint_hi),
                position = position_dodge(width = 1), width = 0, size = 1, col = "dodgerblue3") +
  geom_point(aes(group = name, fill = hdr_significance), stat = "identity", position = position_dodge(width = 1), col = "dodgerblue3", size = 2) +
  scale_fill_manual(values=fill_vals, guide=F) +
  barplot_theme + ylab("Effect size") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
p.barplot_grouped3
pdf("plots/ccdc6_satmut_barplot_grouped3.pdf", width = 6.5, height = 1.75)
p.barplot_grouped3
dev.off()

