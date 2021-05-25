library(tidyverse)
library(cowplot)

setwd("/Users/jeremys/work/opentargets")

getUDPDistributionPlot = function(udp.df, sites) {
  # We allow there to be up to 2 deletions in the UDP
  udp.df = udp.df %>% filter(has_crispr_deletion, !is.na(deletion_start))
  show_udp_sharing = !is.null(udp.df$udp_sharing)
  if (is.null(udp.df$udp_sharing)) {
    udp.df$udp_sharing = "both"
  }
  udp.dels.df = udp.df %>% dplyr::select(deletion_start, deletion_end, deletion2_start, deletion2_end, sharing=udp_sharing)
  udp.dels.df = udp.dels.df %>% arrange(deletion_start, deletion2_start)
  udp.dels.df$y = 1:nrow(udp.dels.df)
  udp.plot.df = bind_rows(udp.dels.df %>% dplyr::select(y, deletion_start, deletion_end, sharing),
                          udp.dels.df %>% dplyr::select(y, deletion_start=deletion2_start, deletion_end=deletion2_end, sharing))
  udp.plot.df = udp.plot.df %>% dplyr::filter(!is.na(deletion_start))
  
  fake_udp.df = data.frame(y = c(0,0,0),
                           deletion_start = c(0,0,0),
                           deletion_end = c(0,0,0),
                           sharing = c("cDNA only", "gDNA only", "both"))
  udp.plot.df = bind_rows(udp.plot.df, fake_udp.df)
  xmax = nchar(udp.df$udp[1])
  segment_size = 0.5
  if (nrow(udp.dels.df) > 200) {
    segment_size = 0.3
  }
  simple_theme = theme_bw(8) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  sharing_colors = c("cDNA only"="red", "gDNA only"="orange", "both"="dodgerblue3", "unclear"="grey")
  p.udp_dist = ggplot(udp.plot.df) +
    geom_segment(aes(x = deletion_start, xend = deletion_end, y = y, yend = y, color = sharing), size = segment_size) +
    simple_theme + coord_cartesian(xlim=c(sites$start, sites$end)) + scale_y_reverse() +
    xlab("Cut site") + ylab("Unique deletions")
  if (show_udp_sharing) {
    p.udp_dist = p.udp_dist + scale_color_manual(values=sharing_colors) +
      theme(legend.position = c(0.3, 0.25), legend.spacing.y = unit(0.1, "cm"), legend.key.height = unit(0.45, "cm"))
  } else {
    p.udp_dist = p.udp_dist + scale_color_manual(values=sharing_colors, guide=F)
  }
  return(p.udp_dist)
}


#udp.df = readr::read_tsv("experiment/transcribed/batch3/analysis/batch3.20bp_edit_window.merged_udps.tsv")
udp.df = readr::read_tsv("experiment/transcribed/batch3/analysis/batch3.20bp_edit_window.MUL1_ABDH4_TAF1C_DRAM2.merged_udps.tsv")

mul1_cut_site = 133
mul1.df = udp.df %>% filter(name == "MUL1_rs6700034")
mul1.df = mul1.df %>%
  mutate(deletion_start = deletion_start - mul1_cut_site + 1,
         deletion_end = deletion_end - mul1_cut_site + 1,
         deletion2_start = deletion2_start - mul1_cut_site + 1,
         deletion2_end = deletion2_end - mul1_cut_site + 1)

sites = list(start = -50, end = 30)

cdna_plot = getUDPDistributionPlot(mul1.df %>% filter(type == "cDNA"), sites)
gdna_plot = getUDPDistributionPlot(mul1.df %>% filter(type == "gDNA"), sites)

cdna_plot = cdna_plot +
  guides(color = FALSE) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())

gdna_plot = gdna_plot +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), legend.background = element_blank(), axis.title.x = element_blank())

paired_udp_plot = cowplot::plot_grid(plotlist = list(gdna_plot, cdna_plot),
                                     ncol = 2)

pdf(file = "plots/genie_udp_plot.pdf", width = 3, height = 2.3)
paired_udp_plot
dev.off()

