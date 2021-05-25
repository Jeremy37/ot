library(tidyverse)
library(cowplot)

setwd("/Users/jeremys/work/opentargets")

barplot_theme = theme_bw(10) +
  theme(axis.line = element_line(colour = "black", size=0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 10))

# First we plot barplots of the effect size per locus, estimated from all replicates
# Because cDNA and gDNA samples aren't paired, we can only look at the cDNA:gDNA ratio
# across all replicates.
#clu_reanalysis.grep.df = readr::read_tsv("experiment/transcribed/clu_reanalysis/analysis/clu_reanalysis.20bp_edit_window.subsampled.grep_analysis.region_stats.tsv")
#clu_reanalysis.dels.df = readr::read_tsv("experiment/transcribed/clu_reanalysis/analysis/clu_reanalysis.20bp_edit_window.subsampled.region_stats.tsv")
clu_reanalysis.grep.df = readr::read_tsv("experiment/transcribed/clu_reanalysis/analysis/clu_reanalysis.subsampled.20bp_edit_window.grep_analysis.region_stats.tsv")
clu_reanalysis.dels.df = readr::read_tsv("experiment/transcribed/clu_reanalysis/analysis/clu_reanalysis.subsampled.20bp_edit_window.region_stats.tsv")

clu_reanalysis.grep.df$hdr.effect[clu_reanalysis.grep.df$name == "CLU_4_b1"] = 0
clu_reanalysis.grep.df$hdr.effect_confint_lo[clu_reanalysis.grep.df$name == "CLU_4_b1"] = NA
clu_reanalysis.grep.df$hdr.effect_confint_hi[clu_reanalysis.grep.df$name == "CLU_4_b1"] = NA
clu_reanalysis.grep.df$name = as.character(seq(1, 11, 1))

clu_reanalysis.dels.df$del.effect[clu_reanalysis.dels.df$name == "CLU_4_g1"] = 0
clu_reanalysis.dels.df$del.effect_confint_lo[clu_reanalysis.dels.df$name == "CLU_4_b1"] = NA
clu_reanalysis.dels.df$del.effect_confint_hi[clu_reanalysis.dels.df$name == "CLU_4_b1"] = NA
clu_reanalysis.dels.df$name = as.character(seq(1, 11, 1))

###############################################################################
library(egg)

statsSummaryPlot = function(grep.summary.df, dels.summary.df, titles = F, legends = T) {
  getSignificanceStr = function(pval) {
    if (is.na(pval) | pval >= 0.01) { "p >= 0.01" }
    else if (pval < 0.001) { "p < 0.001" }
    else { "p < 0.01" }
  }

  p.grep.effect = NULL
  p.stats.effect = NULL
  effect_size_theme = barplot_theme + theme(axis.text.x = element_blank(),
                                           legend.title = element_blank(),
                                           axis.title.x = element_blank(),
                                           plot.margin = unit(c(0.1, 0, 0.1 ,1), "cm"))
  exp_names = grep.summary.df$name
  grep.summary.df$name = factor(as.character(exp_names), levels=exp_names)
  grep.summary.df$hdr_significance = factor(sapply(grep.summary.df$hdr.pval, FUN = getSignificanceStr), levels=c("p >= 0.01", "p < 0.01", "p < 0.001"))

  p.grep.effect = ggplot(grep.summary.df, aes(x=name, y=hdr.effect, fill=hdr_significance)) +
    geom_bar(stat = "identity", width=0.5) +
    geom_hline(yintercept = 1, col = "red") +
    geom_errorbar(aes(ymin = hdr.effect_confint_lo, ymax = hdr.effect_confint_hi),
                  width = 0.2, col = "grey30") +
    effect_size_theme +
    scale_y_continuous(breaks = c(0, 0.4, 0.8, 1.2)) +
    ylab("effect size") +
    scale_fill_manual(values=c(`p >= 0.01`="grey70", `p < 0.01`="cornflowerblue", `p < 0.001`="red3")) +
    coord_cartesian(ylim = c(0, max(1, max(grep.summary.df$hdr.effect * 1.25, na.rm = T))))
  if (!legends) {
    p.grep.effect = p.grep.effect + guides(fill = FALSE)
  }

  exp_names = dels.summary.df$name
  dels.summary.df$name = factor(as.character(exp_names), levels=exp_names)
  dels.summary.df$del_significance = factor(sapply(dels.summary.df$del.pval, FUN = getSignificanceStr), levels=c("p >= 0.01", "p < 0.01", "p < 0.001"))

  p.stats.del = ggplot(dels.summary.df, aes(x=name, y=del.effect, fill=del_significance)) +
    geom_bar(stat = "identity", width=0.5) +
    geom_hline(yintercept = 1, col = "red") +
    geom_errorbar(aes(ymin = del.effect_confint_lo, ymax = del.effect_confint_hi),
                  width = 0.2, col = "grey30") +
    effect_size_theme +
    ylab("effect size") +
    scale_fill_manual(values=c(`p >= 0.01`="grey70", `p < 0.01`="cornflowerblue", `p < 0.001`="red3"), guide = F) +
    coord_cartesian(ylim = c(0, min(3, max(dels.summary.df$del.effect * 1.2, na.rm = T))))

  plot.df = dels.summary.df %>% dplyr::select(name, HDR=hdr.rate, NHEJ=del.rate) %>%
    tidyr::gather(key = "type", value = "value", -name)
  plot.df$type = factor(as.character(plot.df$type), levels = c("NHEJ", "HDR"))
  p.stats.editing = ggplot(plot.df, aes(x=name, y=value*100, fill=type)) +
    geom_bar(stat = "identity", position = position_stack(), width=0.5) +
    barplot_theme + theme(legend.title = element_blank(),
                          axis.title.x = element_blank(),
                          plot.margin = unit(c(0.1, 0, 0.1 ,1), "cm")) +
    ylab("% editing") +
    scale_fill_manual(values=c(HDR="darkorange", NHEJ="cornflowerblue")) +
    coord_cartesian(ylim = c(0, 100))
  if (!legends) {
    p.stats.editing = p.stats.editing + guides(fill = FALSE)
  }

  if (titles == T) {
    p.grep.effect = p.grep.effect + ggtitle("HDR effect size")
    p.stats.del = p.stats.del + ggtitle("Deletion effect size")
    p.stats.editing = p.stats.editing + ggtitle("Editing rates")
  }
  # p.res = cowplot::plot_grid(plotlist = list(p.grep.effect, p.stats.del, p.stats.editing),
  #                            ncol=1, rel_heights=c(1,1,1.1), align="v", axis = "right")
  p.res = egg::ggarrange(plots = list(p.grep.effect, p.stats.del, p.stats.editing),
                         ncol=1, heights=c(1,1,1.1), draw=F)
  p.res
}

p1 = statsSummaryPlot(clu_reanalysis.grep.df, clu_reanalysis.dels.df, titles = T)

pdf("plots/genie_clu_reananalysis.pdf", width=6, height=2.6)
p1
dev.off()

p1 = statsSummaryPlot(clu_reanalysis.grep.df, clu_reanalysis.dels.df, titles = F)

pdf("plots/genie_clu_reananalysis.no_titles.pdf", width=6, height=2.6)
p1
dev.off()



###############################################################################
# BOB_CLU
# bob_clu_ccdc6.grep.df = readr::read_tsv("experiment/transcribed/BOB_CLU_CCDC6/analysis/BOB_CLU_CCDC6.new2.grep_analysis.region_stats.tsv")
# bob_clu_ccdc6.dels.df = readr::read_tsv("experiment/transcribed/BOB_CLU_CCDC6/analysis/BOB_CLU_CCDC6.new2.region_stats.tsv")
bob_clu_ccdc6.grep.df = readr::read_tsv("experiment/transcribed/BOB_CLU_CCDC6/analysis/BOB_CLU_CCDC6.grep_analysis.region_stats.tsv")
bob_clu_ccdc6.dels.df = readr::read_tsv("experiment/transcribed/BOB_CLU_CCDC6/analysis/BOB_CLU_CCDC6.region_stats.tsv")

clu.grep.df = bob_clu_ccdc6.grep.df %>% filter(grepl("CLU", name)) %>%
  mutate(name = gsub("CLU_12", "6", name),
         name = gsub("CLU_14", "7", name),
         name = gsub("CLU_18", "8", name),
         name = factor(as.character(name), levels=c("6", "7", "8")))
clu.dels.df = bob_clu_ccdc6.dels.df %>% filter(grepl("CLU", name)) %>%
  mutate(name = gsub("CLU_12", "6", name),
         name = gsub("CLU_14", "7", name),
         name = gsub("CLU_18", "8", name),
         name = factor(as.character(name), levels=c("6", "7", "8")))

p.bob = statsSummaryPlot(clu.grep.df, clu.dels.df, titles = F, legends = F)

#p2 = statsSummaryPlot(clu.dels.df %>% rename(effect = del.effect, effect_lo = del.effect_confint_lo, effect_hi = del.effect_confint_hi, pval = del.pval), colour = T) +
#  ylab("Del effect")

pdf("plots/genie_clu_bob_plots.pdf", width=2, height=2.4)
p.bob
dev.off()
