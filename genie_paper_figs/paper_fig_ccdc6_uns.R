library(tidyverse)
source("/Users/jeremys/work/opentargets/src/experiment/genie.functions.R")

setwd("/Users/jeremys/work/opentargets")

plot_theme = theme_bw(11) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

getUNSPlot_CCDC6 = function(replicate.udp.df, replicates.df, sites, reverse_amplicon = F, min_gDNA_count = 10, min_cDNA_count = 0, max_udps = 20) {
  uns_res = getUNSData(replicate.udp.df, replicates.df, sites, region_name=NULL, min_gDNA_count, min_cDNA_count, max_udps = max_udps)
  udp.dels.df = uns_res$udp.dels.df %>% dplyr::filter(!is_wt_allele)
  udp.del_replicates.df = uns_res$replicate.dels.df %>% dplyr::filter(!is_wt_allele)
  numUDPs = nrow(udp.dels.df)
  
  if (!is.na(max_udps) & numUDPs > max_udps) {
    udp.dels.df = udp.dels.df %>% .[1:max_udps,]
    udp.del_replicates.df = udp.del_replicates.df %>% dplyr::filter(udp %in% udp.dels.df$udp)
    numUDPs = max_udps
  }
  
  udp.dels.df$uns_gt_3 = udp.dels.df$uns > 3
  #udp.dels.df = udp.dels.df %>% arrange(desc(uns)) %>% mutate(y = 1:nrow(udp.dels.df))
  udp.dels.df = udp.dels.df %>% arrange(desc(uns_gt_3), desc(deletion_start)) %>% mutate(y = 1:nrow(udp.dels.df))
  
  # matrix containing the deletion binary code
  M = as.matrix(do.call(rbind, lapply(as.list(udp.dels.df$udp), udpToBinary)))
  M_chars = as.matrix(do.call(rbind, lapply(as.list(udp.dels.df$udp), strToChars)))
  
  model <- hclust(dist(M))
  dhc <- as.dendrogram(model)
  ddata <- dendro_data(dhc, type = "rectangle")
  dendro_span = max(ddata$segments$y, ddata$segments$yend) - min(ddata$segments$y, ddata$segments$yend)
  
  #udp.order.map = data.frame(y = label(ddata)$x, udp = udp.dels.df[model$order,]$udp)
  #udp.dels.df = udp.dels.df[model$order,] %>% dplyr::left_join(udp.order.map, by="udp")
  
  plot.df = as.data.frame(M_chars)
  # plot.df = as.data.frame(M_chars[model$order,])
  pos = 1:ncol(plot.df)
  if (reverse_amplicon) {
    s = ncol(plot.df) + 1
    pos = ncol(plot.df):1
    sites$start = s - sites$start
    sites$end = s - sites$end
    sites$cut_site = s - sites$cut_site
    sites$highlight_site = s - sites$highlight_site
  }
  colnames(plot.df) = as.character(pos)
  profileSpan = ncol(plot.df)
  plot.df$id = c(1:nrow(plot.df))
  plot.df[,"HDR allele"] = udp.dels.df$is_hdr_allele
  
  plot.gather.df = plot.df %>% tidyr::gather(key = "position", value = "udpchar", 1:profileSpan)
  plot.gather.df$pos = as.numeric(plot.gather.df$position)
  
  p.dendro = ggplot() + 
    geom_segment(aes(x = segment(ddata)$x, y = segment(ddata)$y, xend = segment(ddata)$xend, yend = segment(ddata)$yend)) +
    coord_cartesian(xlim=c(min(plot.gather.df$id), max(plot.gather.df$id)), ylim=c(sites$start, sites$end)) +
    ylab("") + scale_y_continuous(expand = c(0, 0), trans = "reverse") +
    coord_flip() +
    plot_theme + theme(axis.title.y = element_blank(),
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       axis.line = element_blank(),
                       panel.border = element_blank(),
                       panel.grid = element_blank(),
                       plot.margin = unit(c(0,0,0.05,0), "cm"))
  
  numBases = sites$end - sites$start + 1
  dash_size = 2
  dot_size = 0.5
  if (numBases > 60) {
    dash_size = 1.2
    dot_size = 0.4
  }
  if (numBases > 101) {
    dash_size = 0.8
    dot_size = 0.25
  }
  if (any(plot.gather.df$`HDR allele`)) {
    p.udp = ggplot() + 
      geom_rect(aes(xmin=min(pos), xmax=max(pos), ymin=(id-0.5), ymax=(id+0.5)), fill = "palegreen", data = plot.gather.df[plot.gather.df$`HDR allele`,]) +
      geom_point(aes(x = pos, y = id), size = dot_size, shape = 19, alpha = 0.2, data = plot.gather.df[plot.gather.df$udpchar == '-' & plot.gather.df$`HDR allele`,]) +
      geom_point(aes(x = pos, y = id), size = dash_size, shape = '-', alpha = 0.8, data = plot.gather.df[plot.gather.df$udpchar == '*' & plot.gather.df$`HDR allele`,])
  } else {
    p.udp = ggplot()
  }
  p.udp = p.udp + 
    geom_point(aes(x = pos, y = id), size = dot_size, shape = 19, alpha = 0.5, data = plot.gather.df[plot.gather.df$udpchar == '-',]) +
    geom_point(aes(x = pos, y = id), size = dash_size, shape = '-', alpha = 0.8, data = plot.gather.df[plot.gather.df$udpchar == '*',]) +
    scale_color_manual(values = c("grey20", "green4")) +
    geom_point(aes(x = pos, y = id, shape = udpchar), color = "black", size = 2.5, data = plot.gather.df[!(plot.gather.df$udpchar == '*' | plot.gather.df$udpchar == '-'),]) + scale_shape_identity() +
    coord_cartesian(ylim=c(min(plot.gather.df$id), max(plot.gather.df$id)), xlim=c(min(c(sites$start, sites$end)), max(c(sites$start, sites$end)))) +
    xlab("Amplicon position") + scale_x_continuous(expand = c(0.01, 0)) + 
    plot_theme + theme(legend.position = "none",
                       axis.title.y = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       axis.line = element_blank(),
                       panel.border = element_blank(),
                       plot.margin = unit(c(0,0,0.05,0), "cm"))
  # Used to colour the HDR allele text, but now we show a rectangle around HDR alleles
  #geom_point(aes(x = pos, y = id, colour = `HDR allele`), size = dot_size, shape = 19, alpha = 0.7, data = plot.gather.df[plot.gather.df$udpchar == '-',]) +
  
  if (!is.na(sites$cut_site)) {
    p.udp = p.udp + geom_vline(xintercept = sites$cut_site, color="grey20", linetype = "longdash", alpha=0.5)
  }
  
  # Left side of the plot, showing UDP-normalised scores (UNS)
  # uns.plot.df = data.frame(y = label(ddata)$x,
  #                          x = 0,
  #                          uns = udp.dels.df$uns,
  #                          cDNA_lower = udp.dels.df$uns < 1,
  #                          gDNA = udp.dels.df$gDNA_total_count,
  #                          uns_conf_hi = udp.dels.df$uns_confint_lo,
  #                          uns_conf_lo = udp.dels.df$uns_confint_hi)
  uns.plot.df = data.frame(y = 1:nrow(udp.dels.df),
                           x = 0,
                           uns = udp.dels.df$uns,
                           cDNA_lower = udp.dels.df$uns < 1,
                           gDNA = udp.dels.df$gDNA_total_count,
                           uns_conf_hi = udp.dels.df$uns_confint_lo,
                           uns_conf_lo = udp.dels.df$uns_confint_hi)
  
  # uns.replicates.plot.df = udp.del_replicates.df %>% 
  #   dplyr::select(udp, uns, gDNA_count) %>%
  #   dplyr::left_join(udp.order.map, by="udp")
  
  maxDotSize = 4.5
  if (numUDPs > 60) {
    maxDotSize = 3.5
  } else if (numUDPs > 30) {
    maxDotSize = 4.5
  }
  
  min_uns_threshold = 0.25
  uns.plot.df$uns = sapply(uns.plot.df$uns, FUN = function(x) max(x, min_uns_threshold))
  uns.plot.df$uns_conf_lo = sapply(uns.plot.df$uns_conf_lo, FUN = function(x) max(x, 0.01))
  max_uns_display = max(1.5, ceiling(max(uns.plot.df$uns, na.rm = T)))
  min_uns_display = min(uns.plot.df$uns, na.rm = T)
  
  # Code testing different visualisation for representing confidence in UNS value
  #log2_uns_range = log2(max(uns.plot.df$uns)) - log2(min(uns.plot.df$uns))
  #max_log_gDNA = max(log10(uns.plot.df$gDNA))
  #uns.plot.df$confidence = sapply(1:nrow(uns.plot.df), FUN=function(i) min(log10(max(uns.plot.df[i,]$gDNA-10, 10)) / max_log_gDNA, 1) * abs(log2(uns.plot.df[i,]$uns)) / log2_uns_range)
  #uns.plot.df$confidence = (uns.plot.df$alpha - min(uns.plot.df$alpha)) / (max(uns.plot.df$alpha) - min(uns.plot.df$alpha))
  uns.plot.df$expr = "higher"
  uns.plot.df$expr[uns.plot.df$cDNA_lower] = "lower"
  uns.plot.df$spliceDisrupting = uns.plot.df$uns > 3
  p.uns = ggplot(uns.plot.df) +
    #annotate("segment", x = uns.plot.df$x, xend = uns.plot.df$uns, y = uns.plot.df$y, yend = uns.plot.df$y, colour = 'gray') +
    geom_point(aes(x = uns, y = y, colour = spliceDisrupting, size = log10(gDNA)), alpha = 0.7) +
    geom_errorbarh(aes(y = y, xmin = uns_conf_lo, xmax = uns_conf_hi, height = 0.4), colour = "grey20") +
    scale_size(range = c(1, maxDotSize)) +
    scale_color_manual(values=c("TRUE"="red", "FALSE"="blue"), guide=F) +
    scale_x_continuous(trans="log2", breaks = c(0.5, 1, 2, 4, 8)) +
    #scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, limits = c(0.1, 4)) +
    #geom_point(data = uns.replicates.plot.df, mapping = aes(y = y, x = uns, size = 2, alpha = 0.6)) +
    scale_y_continuous(position = "right") +
    coord_cartesian(xlim = c(min(0.8, min_uns_display), max_uns_display),
                    ylim=c(min(plot.gather.df$id), max(plot.gather.df$id))) +
    ylab("Unique deletion profile") + xlab("Relative expr") +
    geom_vline(xintercept = 1.0, alpha = 0.25, linetype = "longdash") +
    plot_theme + theme(legend.position = "right",
                       #axis.line = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       panel.border = element_blank(),
                       plot.margin = unit(c(0,0,0.05,0), "cm"))
#  scale_color_manual(values=c(higher="red", lower="blue")) +
    #  geom_point(aes(x = uns, y = y, colour = expr, size = log10(gDNA), alpha = log10(gDNA))) +
    #  scale_alpha_continuous(range = c(0.2, 1)) + 
    
  p.udp_profile = egg::ggarrange(plots = list(p.udp, p.uns), nrow=1, ncol=2, widths=c(4,1), draw = F)
  #p.udp_profile = plot_grid(p.uns, p.udp, p.dendro, nrow=1, ncol=3, rel_widths = c(2,4,1))
  return(p.udp_profile)
}


###############################################################################

default_opts = list(
  grep_analysis = T,
  grep_window = 10,
  grep_window_left = NULL,
  grep_window_right = NULL,
  deletion_analysis = T,
  minMapQ = 0,
  subsample = NULL,
  max_mismatch_frac = 0.05,
  viewing_window = 40,
  editing_window = 20, # window around the cut site for counting deletions as edits
  min_window_overlap = 30, # minimum number of read bases correctly aligned within the region of interest
  exclude_multiple_deletions = F,
  exclude_nonspanning_reads = T,
  exclude_nonspanning_deletions = T,
  ratio_to_total_reads = F,
  use_cdna_dels_only = F,
  qc_plot_max_udps = 10,
  qc_plot_min_udp_fraction = 0.002,
  qc_plot_exclude_wt = T,
  del_span_start = NULL,
  del_span_end = NULL,
  uns_plot_min_gDNA = 10,
  uns_plot_min_cDNA = 0,
  uns_plot_max_udps = 40,
  allele_profile = F,
  no_site_profile = F,
  no_udp_profile = F,
  no_uns_profile = F,
  no_stats = F,
  replicate_qc_plots = T,
  no_replicate_plots = T,
  variance_analysis = F,
  power_analysis = F,
  plot_width = 8,
  plot_height = 7,
  variance_analysis_min_count = 100,
  variance_analysis_min_fraction = 0.001)

opts <<- default_opts
opts$regions = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1_updated/batch1_updated.regions.ref.tsv"
opts$replicates = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1_updated/batch1_updated.meta.tsv"
opts$read_data = "analysis/batch1_updated.read_data.rds"
opt = opts
#system.time(runGenIE(opts))

replicates.df = readr::read_tsv(opts$replicates) %>% filter(name == "26_CCDC6, rs1171830")

replicate.udp.df = readr::read_tsv("experiment/transcribed/batch1_updated/analysis/batch1_updated.edit_dels.20bp_window.new.replicate_udps.tsv") %>%
  filter(name == "26_CCDC6, rs1171830")
replicate.udp.df.single_dels = replicate.udp.df %>% filter(!grepl("\\*+[^\\*]+\\*+", udp))

sites = list(start = 27, end = 132, cut_site = 82, highlight_site = 82)

pdf("plots/genie_ccdc6_uns.v2.pdf", width = 6.5, height = 2.8)
getUNSPlot_CCDC6(replicate.udp.df.single_dels, replicates.df, sites, reverse_amplicon = T, min_gDNA_count = 10, min_cDNA_count = 0, max_udps = 16)
dev.off()

