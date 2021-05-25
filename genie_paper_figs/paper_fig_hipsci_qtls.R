library(tidyverse)
library(egg) # theme_article()
library(ggbeeswarm)
setwd("/Users/jeremys/work/opentargets")

hipsci.df = read_tsv("plots/hipsci_qtl/hipsci_qtl_data.tsv")
hipsci.df = hipsci.df %>% mutate(index = row_number())

qtlViolinPlot = function(df, ptSize = 0.9) {
  ggplot(df, aes(x=genotype, y=value, fill="blue")) +
    geom_beeswarm(col="dodgerblue3", size=ptSize, cex=0.7) +
    geom_violin(alpha=0.2) +
    geom_boxplot(outlier.shape = NA, width=0.12, fill="white") +
    scale_fill_manual(values = c("grey50"), guide=F) +
    theme_article(24) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  #geom_jitter(col="dodgerblue3", size=1) +
  #geom_dotplot(binaxis = "y", stackdir = "center", dotsize=0.1, position=position_dodge2()) +
  #geom_jitter()
}

p.dram2 = qtlViolinPlot(hipsci.df %>% filter(locus == "chr1:111680208:111682123:clu_25929"))
p.dram2

p.taf1c = qtlViolinPlot(hipsci.df %>% filter(locus == "chr16:84218595:84220507:clu_8787"))

p.mul1 = qtlViolinPlot(hipsci.df %>% filter(locus == "ENSG00000090432"))

p.abhd4 = qtlViolinPlot(hipsci.df %>% filter(locus == "ENSG00000100439"))

#p.sdf4 = qtlViolinPlot(hipsci.df %>% filter(locus == "chr1:1164326:1166887:clu_24394"))
#p.sdf4
p.sdf4 = qtlViolinPlot(hipsci.df %>% filter(locus == "chr1:1164326:1167272:clu_24394"), ptSize = 0.7)
p.sdf4


pdf("plots/hipsci_qtl/hipsci_qtls.pdf", width=3.5, height=3.2)
p.sdf4
p.taf1c
p.mul1
p.abhd4
dev.off()
