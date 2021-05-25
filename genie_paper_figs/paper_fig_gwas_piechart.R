library(tidyverse)

setwd("/Users/jeremys/work/opentargets")

#gc.df = readr::read_tsv("/Users/jeremys/work/opentargets/reference/gwas_catalog/gwas_catalog_v1.0.2-associations_e96_r2019-10-14.tsv.gz")
#gc.df = gc.df  %>% filter(!duplicated(PUBMEDID, SNPS))
#write_tsv(gc.df, "reference/gwas_catalog/gwas_catalog_v1.0.2-associations_e96_r2019-10-14.dedup.tsv")
gc.df = read_tsv("reference/gwas_catalog/gwas_catalog_v1.0.2-associations_e96_r2019-10-14.dedup.tsv")

freqs = table(gc.df$CONTEXT)
freqs.plot = freqs[freqs > 5]

rename_df = data.frame(annot = c("3_prime_UTR_variant", "5_prime_UTR_variant", "intergenic_variant", "intron_variant", "missense_variant", "non_coding_transcript_exon_variant", "regulatory_region_variant", "splice_acceptor_variant", "splice_region_variant", "stop_gained", "synonymous_variant", "TF_binding_site_variant"),
                       annot_name = c("UTR", "UTR", "intergenic", "intron", "missense", "non-coding transcript", "regulatory region", "splice / NMD", "splice / NMD", "splice / NMD", "synonymous", "TF binding site"),
                       annot_group = c("UTR", "UTR", "intergenic", "intron", "missense", "non-coding transcript", "intergenic", "splice / NMD", "splice / NMD", "splice / NMD", "synonymous", "intergenic"))

freqs.plot = data.frame(annot = names(freqs.plot), count_det = as.numeric(freqs.plot)) %>%
  left_join(rename_df, by="annot")

freqs.plot.df = freqs.plot %>%
  group_by(annot_name) %>% summarise(count = sum(count_det)) %>%
  mutate(pct = count / sum(freqs.plot$count_det) * 100,
         pctlabel = paste0(round(pct), "%"))

freqs.group.plot.df = freqs.plot %>%
  group_by(annot_group) %>% summarise(count = sum(count_det)) %>%
  mutate(pct = count / sum(freqs.plot$count_det) * 100,
         pctlabel = paste0(round(pct), "%"))

pdf("plots/gwas_catalog_pie.pdf", width = 8, height = 6)
################################################
# Standard R pie chart
library(RColorBrewer)
myPalette <- brewer.pal(nrow(freqs.plot.df), "Set2") 
pie(freqs.plot.df$pct, labels = freqs.plot.df$annot_name, col=rainbow(nrow(freqs.plot.df)), cex = 0.7)

myPalette <- brewer.pal(nrow(freqs.group.plot.df), "Set2") 
pie(freqs.group.plot.df$pct, labels = sprintf("%s (%s)", freqs.group.plot.df$annot_group, freqs.group.plot.df$pctlabel), col=myPalette, radius = 0.7, cex = 0.7)

################################################
# 3d Pie
library(plotrix)
pie3D(freqs.group.plot.df$pct,
      labels = sprintf("%s (%s)", freqs.group.plot.df$annot_group, freqs.group.plot.df$pctlabel),
      explode = 0, radius = 1, labelcex = 0.7,
      start = 0.5,
      col = rainbow(nrow(freqs.group.plot.df)))


################################################
# ggplot pie chart using coord_polar
palette1 = c("#33658A", # intergenic
             "#F6AE2D", # intron
             "#F26419", # missense
             "#55DDE0", # non-coding
             "#2F4858", # splice / NMD
             "#999999", # synonymous
             "#992999") # UTR

ggpie = ggplot(freqs.group.plot.df, aes(x="", y=pct, fill=annot_group)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=1.25) +
  scale_fill_manual(values = palette1) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic(9) + theme(axis.line = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks = element_blank())

ggpie

ggpie + geom_text(aes(label = pctlabel), position = position_stack(vjust = 0.5), size = 3)

# Make a palette which shows editable types in light colours and non-editable dark
palette2 = c("#33658A", # intergenic
             "#F6AE2D", # intron
             "#F26419", # missense
             "#ffe342", # non-coding
             "#ff4242", # splice / NMD
             "#bdff52", # synonymous
             "#ff9e61") # UTR
ggpie + scale_fill_manual(values = palette2)


# Order categories to place noncoding transcribed next to each other
freqs.group.plot.df = freqs.group.plot.df %>%
  mutate(annot_group = factor(as.character(annot_group), levels = c("intergenic", "intron", "UTR", "non-coding transcript", "splice / NMD", "missense", "synonymous")))

ggpie = ggplot(freqs.group.plot.df, aes(x="", y=pct, fill=annot_group)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=1.25) +
  scale_fill_manual(values = palette1) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic(9) + theme(axis.line = element_blank(),
                           axis.text = element_blank(),
                           axis.ticks = element_blank())
# display.brewer.pal(4, "YlOrRd")
# display.brewer.pal(9, "YlOrRd")
# display.brewer.pal(9, "Blues")
# display.brewer.pal(9, "Greens")
# display.brewer.pal(9, "YlOrRd")

palette3 = c(brewer.pal(9, "Blues")[7], # intergenic
             brewer.pal(9, "YlOrRd")[4], # intron
             brewer.pal(9, "YlOrRd")[5], # UTR
             brewer.pal(9, "YlOrRd")[3], # non-coding transcript
             brewer.pal(9, "YlOrRd")[6], # splice / NMD
             brewer.pal(9, "Greens")[7], # missense
             brewer.pal(9, "Greens")[5]) # synonymous
ggpie + scale_fill_manual(values = palette3)


dev.off()

