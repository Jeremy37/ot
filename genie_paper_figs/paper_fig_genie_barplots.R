library(tidyverse)
library(cowplot)

setwd("/Users/jeremys/work/opentargets")

barplot_theme = theme_bw(8) +
  theme(axis.line = element_line(colour = "black", size=0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# First we plot barplots of the effect size per locus, estimated from all replicates
# Because cDNA and gDNA samples aren't paired, we can only look at the cDNA:gDNA ratio
# across all replicates.
batch3.grep.df = readr::read_tsv("experiment/transcribed/batch3/analysis/batch3.20bp_edit_window.grep_analysis.region_stats.tsv")
batch3.dels.df = readr::read_tsv("experiment/transcribed/batch3/analysis/batch3.20bp_edit_window.region_stats.tsv")

# batch3.grep.df$name = factor(as.character(batch3.grep.df$name), levels=unique(batch3.grep.df$name))
# batch3.dels.df$name = factor(as.character(batch3.dels.df$name), levels=unique(batch3.dels.df$name))

genieBarplot = function(stats.df, colour = F) {
  stats.df$sig = "p >= 0.01"
  stats.df$sig[stats.df$pval < 0.01] = "p < 0.01"
  stats.df$sig[stats.df$pval < 0.001] = "p < 0.001"
  colorVals = c(`p >= 0.01`="grey70", `p < 0.01`="cornflowerblue", `p < 0.001`="red3")
  if (!colour) {
    colorVals = c(`p >= 0.01`="cornflowerblue", `p < 0.01`="cornflowerblue", `p < 0.001`="cornflowerblue")
  }
  ggplot(stats.df, aes(x=name, y=effect, fill = sig)) +
    geom_bar(stat = "identity", width=0.6) +
    geom_hline(yintercept = 1, col = "red") +
    geom_errorbar(aes(ymin = effect_lo, ymax = effect_hi),
                  width = 0.25, size = 0.3, col = "grey30") +
    scale_fill_manual(values=colorVals, guide=F) +
    barplot_theme +
    theme(axis.title.x = element_blank())
#  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
}

genieEditingPieChart = function(dels.df) {
  editing.plot.df = dels.df %>% dplyr::select(name, HDR=hdr.rate, NHEJ=del.rate) %>%
    mutate(WT = 1 - HDR - NHEJ) %>%
    tidyr::gather(key = "type", value = "value", -name)
  editing.plot.df$type = factor(as.character(editing.plot.df$type), levels = c("WT", "NHEJ", "HDR"))
  editing.plot.df = editing.plot.df %>%
    arrange(desc(type)) %>%
    mutate(prop = value / sum(value) * 100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  ggpie = ggplot(editing.plot.df, aes(x="", y=prop, fill=type)) +
    geom_bar(stat="identity", width=0.1) +
    coord_polar("y", start=0) +
    #    geom_text(aes(x=1.2, y = ypos, label = type), size=3) +
    scale_fill_manual(values = c(HDR="darkorange", NHEJ="cornflowerblue", WT="grey85"), guide=F) +
    labs(x = editing.plot.df$name[1], y = NULL, fill = NULL) +
    theme_classic(9) + theme(axis.line = element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank())
  ggpie
}


###############################################################################
# MUL1 - ABHD4

mul1_abhd4.grep.df = batch3.grep.df %>% filter(grepl("MUL1|ABHD4", name)) %>%
  mutate(name = gsub("MUL1_rs6700034", "MUL1", name),
         name = gsub("ABHD4_rs8011143", "ABHD4", name),
         name = factor(as.character(name), levels=c("MUL1", "ABHD4"))) %>%
  select(name, effect = hdr.effect, effect_lo = hdr.effect_confint_lo, effect_hi = hdr.effect_confint_hi, pval = hdr.pval)
mul1_abhd4.dels.df = batch3.dels.df %>% filter(grepl("MUL1|ABHD4", name)) %>%
  mutate(name = gsub("MUL1_rs6700034", "MUL1", name),
         name = gsub("ABHD4_rs8011143", "ABHD4", name),
         name = factor(as.character(name), levels=c("MUL1", "ABHD4"))) %>%
  select(name, effect = del.effect, effect_lo = del.effect_confint_lo, effect_hi = del.effect_confint_hi, pval = del.pval)

mul1.bars.df = bind_rows(mul1_abhd4.grep.df %>% filter(name == "MUL1") %>% mutate(name = gsub("MUL1", "HDR", name)),
                         mul1_abhd4.dels.df %>% filter(name == "MUL1") %>% mutate(name = gsub("MUL1", "Dels", name)))
mul1.bars.df$name = factor(as.character(mul1.bars.df$name), levels = c("HDR", "Dels"))

p1 = genieBarplot(mul1.bars.df) + ylab("Effect")
pdf("plots/genie_mul1_plot.pdf", width=1, height=1)
p1
dev.off()

abhd4.bars.df = bind_rows(mul1_abhd4.grep.df %>% filter(name == "ABHD4") %>% mutate(name = gsub("ABHD4", "HDR", name)),
                          mul1_abhd4.dels.df %>% filter(name == "ABHD4") %>% mutate(name = gsub("ABHD4", "Dels", name)))
abhd4.bars.df$name = factor(as.character(abhd4.bars.df$name), levels = c("HDR", "Dels"))

p2 = genieBarplot(abhd4.bars.df) + ylab("Effect")
pdf("plots/genie_abhd4_plot.pdf", width=1, height=1)
p2
dev.off()

pdf("plots/genie_abhd4_mul1_merged_bars.pdf", width=1.2, height=2)
cowplot::plot_grid(p1, p2, ncol = 1, align = "v")
dev.off()


pdf("plots/genie_mul1_edit_pie.pdf", width=1, height=1)
genieEditingPieChart(batch3.dels.df %>% filter(grepl("MUL1", name)))
dev.off()

pdf("plots/genie_abhd4_edit_pie.pdf", width=1, height=1)
genieEditingPieChart(batch3.dels.df %>% filter(grepl("ABHD4", name)))
dev.off()

###############################################################################
# MUL1 guide 2

batch4.grep.df = readr::read_tsv("experiment/transcribed/batch4_MUL1_CD33/analysis/batch4_MUL1.20bp_edit_window.grep_analysis.region_stats.tsv")
batch4.dels.df = readr::read_tsv("experiment/transcribed/batch4_MUL1_CD33/analysis/batch4_MUL1.20bp_edit_window.region_stats.tsv")

mul1_g2.grep.df = batch4.grep.df %>% filter(grepl("MUL1", name)) %>%
  mutate(name = gsub("MUL1_rs6700034", "MUL1", name),
         name = factor(as.character(name), levels=c("MUL1"))) %>%
  select(name, effect = hdr.effect, effect_lo = hdr.effect_confint_lo, effect_hi = hdr.effect_confint_hi, pval = hdr.pval)
mul1_g2.dels.df = batch4.dels.df %>% filter(grepl("MUL1", name)) %>%
  mutate(name = gsub("MUL1_rs6700034", "MUL1", name),
         name = factor(as.character(name), levels=c("MUL1"))) %>%
  select(name, effect = del.effect, effect_lo = del.effect_confint_lo, effect_hi = del.effect_confint_hi, pval = del.pval)

mul1.bars.df = bind_rows(mul1_g2.grep.df %>% filter(name == "MUL1") %>% mutate(name = gsub("MUL1", "HDR", name)),
                         mul1_g2.dels.df %>% filter(name == "MUL1") %>% mutate(name = gsub("MUL1", "Dels", name)))
mul1.bars.df$name = factor(as.character(mul1.bars.df$name), levels = c("HDR", "Dels"))

p1 = genieBarplot(mul1.bars.df) + ylab("Effect")
pdf("plots/genie_mul1_g2_plot.pdf", width=1, height=1)
p1
dev.off()

pdf("plots/genie_mul1_g2_edit_pie.pdf", width=1, height=1)
genieEditingPieChart(batch4.dels.df %>% filter(grepl("MUL1", name)))
dev.off()


####################### OLD ##############################
# p1 = genieBarplot(mul1_abhd4.grep.df %>% rename(effect = hdr.effect, effect_lo = hdr.effect_confint_lo, effect_hi = hdr.effect_confint_hi, pval = hdr.pval)) +
#   ylab("HDR effect")
#
# p2 = genieBarplot(mul1_abhd4.dels.df %>% rename(effect = del.effect, effect_lo = del.effect_confint_lo, effect_hi = del.effect_confint_hi, pval = hdr.pval)) +
#   ylab("Del effect")
#
# genieEditingPlot = function(dels.df) {
#   editing.plot.df = dels.df %>% dplyr::select(name, HDR=hdr.rate, NHEJ=del.rate) %>%
#     tidyr::gather(key = "type", value = "value", -name)
#   editing.plot.df$type = factor(as.character(editing.plot.df$type), levels = c("NHEJ", "HDR"))
#   ggplot(editing.plot.df, aes(x=name, y=value*100, fill=type)) +
#     geom_bar(stat = "identity", position = position_stack(), width=0.6) +
#     ylab("% editing") +
#     scale_fill_manual(values=c(HDR="darkorange", NHEJ="cornflowerblue")) +
#     barplot_theme +
#     theme(legend.title = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
# }
#
# p1
# p2
# p.editing = genieEditingPlot(mul1_abhd4.dels.df)
# p.editing
# p.editing = p.editing + scale_fill_manual(values=c(HDR="darkorange", NHEJ="cornflowerblue"), guide=F)
# p.editing
# dev.off()
#
# pdf("plots/genie_mul1_abhd4_plots.merged.pdf", width=1, height=2.5)
# cowplot::plot_grid(plotlist = list(p1 + theme(axis.text.x = element_blank()),
#                                    p2 + theme(axis.text.x = element_blank()),
#                                    p.editing), ncol = 1, rel_heights = c(0.3, 0.3, 0.4), align = "v")
# dev.off()
####################### end OLD ##############################


# batch3.repl.grep.df = readr::read_tsv("experiment/transcribed/batch3/analysis/batch3.20bp_edit_window.grep_analysis.replicate_stats.tsv")
# batch3.repl.dels.df = readr::read_tsv("experiment/transcribed/batch3/analysis/batch3.20bp_edit_window.replicate_stats.tsv")
# batch3.repl.grep.df$name = factor(as.character(batch3.repl.grep.df$name), levels=unique(batch3.repl.grep.df$name))
# batch3.repl.dels.df$name = factor(as.character(batch3.repl.dels.df$name), levels=unique(batch3.repl.dels.df$name))
# mul1_abhd4.repl.grep.df = batch3.repl.grep.df %>% filter(grepl("MUL1|ABHD4", name))
# mul1_abhd4.repl.dels.df = batch3.repl.dels.df %>% filter(grepl("MUL1|ABHD4", name))


###############################################################################
# TAF1C - DRAM2

# taf1c_dram2.grep.df = batch3.grep.df %>% filter(grepl("TAF1C|DRAM2", name)) %>%
#   mutate(name = gsub("TAF1C_rs4150126", "TAF1C", name),
#          name = gsub("DRAM2_rs3762374_TT", "DRAM2-TT", name),
#          name = gsub("DRAM2_rs3762374_TC", "DRAM2-TC", name),
#          name = gsub("DRAM2_rs3762374_AT_haplotype", "DRAM2-AT", name),
#          name = factor(as.character(name), levels=c("TAF1C", "DRAM2-TT", "DRAM2-TC", "DRAM2-AT")))
taf1c.grep.df = batch3.grep.df %>% filter(grepl("TAF1C", name))
taf1c.dels.df = batch3.dels.df %>% filter(grepl("TAF1C", name))

taf1c.bars.df = bind_rows(taf1c.grep.df %>%
                            select(name, effect = hdr.effect, effect_lo = hdr.effect_confint_lo, effect_hi = hdr.effect_confint_hi, pval = hdr.pval) %>%
                            mutate(name = gsub("TAF1C_rs4150126", "HDR", name)),
                          taf1c.dels.df %>%
                            select(name, effect = del.effect, effect_lo = del.effect_confint_lo, effect_hi = del.effect_confint_hi, pval = del.pval) %>%
                            mutate(name = gsub("TAF1C_rs4150126", "Dels", name)))
taf1c.bars.df$name = factor(as.character(taf1c.bars.df$name), levels = c("HDR", "Dels"))

p1 = genieBarplot(taf1c.bars.df) + ylab("Effect")
pdf("plots/genie_taf1c.barplot.pdf", width=1, height=1)
p1
dev.off()

sdf4.grep.df = readr::read_tsv("experiment/transcribed/batch3/analysis/batch3.SDF4_exonic.grep.grep_analysis.region_stats.tsv") %>%
  filter(grepl("_vs_TG", name)) %>%
  mutate(name = gsub("SDF4_CC_vs_TG", "CC", name),
         name = gsub("SDF4_TC_vs_TG", "TC", name),
         name = gsub("SDF4_CG_haplotype_vs_TG", "CG", name))

# Get deletion DF from intronic analysis for SDF4 only to get editing rate
# Set deletion effect to null, since it's not reliable
sdf4.dels.df = readr::read_tsv("experiment/transcribed/batch3/analysis/batch3.SDF4_intronic.grep.region_stats.tsv") %>%
  filter(grepl("SDF4_CC_vs_TG", name)) %>%
  mutate(del.effect = NA, del.effect_confint_lo = NA, del.effect_confint_hi = NA)

sdf4.bars.df = sdf4.grep.df %>%
  select(name, effect = hdr.effect, effect_lo = hdr.effect_confint_lo, effect_hi = hdr.effect_confint_hi, pval = hdr.pval) %>%
  mutate(name = factor(as.character(name), levels=c("CC", "CG", "TC")))

p2 = genieBarplot(sdf4.bars.df) + ylab("HDR Effect")
pdf("plots/genie_sdf4.barplot.pdf", width=1, height=1)
p2
dev.off()

pdf("plots/genie_taf1c_sdf4_plots.all.pdf", width=1, height=1)
p1
p2

p.taf1c.pie = genieEditingPieChart(taf1c.dels.df)
p.taf1c.pie
p.sdf4.pie = genieEditingPieChart(sdf4.dels.df)
p.sdf4.pie

dev.off()

####################### OLD ##############################

# taf1c_sdf4.grep.df = bind_rows(taf1c.grep.df, sdf4.grep.df) %>%
#   mutate(name = gsub("TAF1C_rs4150126", "TAF1C", name),
#          name = gsub("SDF4_CC_vs_TG", "CC", name),
#          name = gsub("SDF4_TC_vs_TG", "TC", name),
#          name = gsub("SDF4_CG_haplotype_vs_TG", "CG", name),
#          name = factor(as.character(name), levels=c("TAF1C", "CC", "TC", "CG")))
#
# taf1c_sdf4.dels.df = bind_rows(taf1c.dels.df, sdf4.dels.df) %>%
#   mutate(name = gsub("TAF1C_rs4150126", "TAF1C", name),
#          name = gsub("SDF4_CC_vs_TG", "SDF4", name),
#          name = factor(as.character(name), levels=c("TAF1C", "SDF4")))
#
# pdf("plots/genie_taf1c_sdf4_plots.all.pdf", width=1, height=1)
#
# p1 = genieBarplot(taf1c_sdf4.grep.df %>% rename(effect = hdr.effect, effect_lo = hdr.effect_confint_lo, effect_hi = hdr.effect_confint_hi, pval = hdr.pval)) +
#   ylab("HDR effect")
#
# p2 = genieBarplot(taf1c_sdf4.dels.df %>% rename(effect = del.effect, effect_lo = del.effect_confint_lo, effect_hi = del.effect_confint_hi, pval = del.pval)) +
#   ylab("Del effect")
#
# p1
# p2
#
# p.editing = genieEditingPlot(taf1c_sdf4.dels.df)
# p.editing
# p.editing = p.editing + scale_fill_manual(values=c(HDR="darkorange", NHEJ="cornflowerblue"), guide=F)
# p.editing
#
# dev.off()
#
# pdf("plots/genie_taf1c_sdf4_plots.hdr.pdf", width=1, height=0.85)
# p1
# dev.off()
#
# pdf("plots/genie_taf1c_sdf4_plots.merged.pdf", width=1, height=1.5)
# cowplot::plot_grid(plotlist = list(p2 + theme(axis.text.x = element_blank()),
#                                    p.editing), ncol = 1, rel_heights = c(0.4, 0.55), align = "v")
# dev.off()
####################### end OLD ##############################


###############################################################################
# BOB_CLU
bob_clu_ccdc6.grep.df = readr::read_tsv("experiment/transcribed/BOB_CLU_CCDC6/analysis/BOB_CLU_CCDC6.new2.grep_analysis.region_stats.tsv")
bob_clu_ccdc6.dels.df = readr::read_tsv("experiment/transcribed/BOB_CLU_CCDC6/analysis/BOB_CLU_CCDC6.new2.region_stats.tsv")

clu.grep.df = bob_clu_ccdc6.grep.df %>% filter(grepl("CLU", name)) %>%
  mutate(name = gsub("CLU_12", "CLU-6", name),
         name = gsub("CLU_14", "CLU-7", name),
         name = gsub("CLU_18", "CLU-8", name),
         name = factor(as.character(name), levels=c("CLU-6", "CLU-7", "CLU-8")))
clu.dels.df = bob_clu_ccdc6.dels.df %>% filter(grepl("CLU", name)) %>%
  mutate(name = gsub("CLU_12", "CLU-6", name),
         name = gsub("CLU_14", "CLU-7", name),
         name = gsub("CLU_18", "CLU-8", name),
         name = factor(as.character(name), levels=c("CLU-6", "CLU-7", "CLU-8")))

barplot_theme = theme_bw(10) +
  theme(axis.line = element_line(colour = "black", size=0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

p1 = genieBarplot(clu.grep.df %>% rename(effect = hdr.effect, effect_lo = hdr.effect_confint_lo, effect_hi = hdr.effect_confint_hi, pval = hdr.pval), colour = T) +
  ylab("HDR effect")

p2 = genieBarplot(clu.dels.df %>% rename(effect = del.effect, effect_lo = del.effect_confint_lo, effect_hi = del.effect_confint_hi, pval = del.pval), colour = T) +
  ylab("Del effect")

pdf("plots/genie_clu_bob_plots.2.pdf", width=2, height=1.3)
p1
p2

p.editing = genieEditingPlot(clu.dels.df)
p.editing
p.editing = p.editing + scale_fill_manual(values=c(HDR="darkorange", NHEJ="cornflowerblue"), guide=F)
p.editing

dev.off()

pdf("plots/genie_clu_bob_plots.merged2.pdf", width=2, height=3)
cowplot::plot_grid(plotlist = list(p1 + theme(axis.text.x = element_blank()),
                                   p2 + theme(axis.text.x = element_blank()),
                                   p.editing), ncol = 1, rel_heights = c(0.3, 0.3, 0.4), align = "v")
dev.off()

