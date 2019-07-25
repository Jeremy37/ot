#!/usr/bin/env Rscript

# This script reads in a VCF file and an "ASECounts" file. The VCF file has the genotypes
# for a single individual in the region of 1 gene. The ASECounts file has allele-specific
# read counts for possibly multiple samples that come from the individual and which cover
# the gene. The script performs statistical tests for allele-specific expression, and
# produces some plots.
suppressMessages(library(tidyverse))
suppressMessages(library(stringr))
suppressMessages(library(egg)) # for ggarrange
suppressMessages(library(cowplot)) # could probably use just one of these two plotting libraries
suppressMessages(library(aod)) # for betabin
suppressMessages(library(ggdendro))
suppressMessages(library(variancePartition))
suppressMessages(library(doParallel))
suppressMessages(library(isofor))
suppressMessages(library(FNN))
suppressMessages(library(broom))
suppressMessages(library(lmerTest))

cl <- makeCluster(4)
registerDoParallel(cl)
#library(profvis)
options(stringsAsFactors = F)

opt = NULL

runGenIE = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  if (!is.null(opt$del_span_start) & is.null(opt$del_span_end) | is.null(opt$del_span_start) & !is.null(opt$del_span_end)) {
    stop("Both options 'del_span_start' and 'del_span_end' must be specified together. Only one was given.")
  }
  opt$custom_del_span <<- F
  if (!is.null(opt$del_span_start) & !is.null(opt$del_span_end)) {
    opt$custom_del_span <<- T
  }
  
  if (!opt$deletion_analysis & !opt$grep_analysis) {
    stop("One of the options 'deletion_analysis' and 'grep_analysis' must be specified.")
  }
  
  if (is.null(opt$grep_window_left)) { opt$grep_window_left <<- opt$grep_window }
  if (is.null(opt$grep_window_right)) { opt$grep_window_right <<- opt$grep_window }
  
  regions.df = readr::read_tsv(opt$regions, col_types="cciiiiccc")
  if (any(duplicated(regions.df$name))) {
    duplicated_region = regions.df$name[ which(duplicated(regions.df$name))[1] ]
    stop(sprintf("All region names must be unique. Found duplicated region name: %s.", duplicated_region))
  }
  region_names = unique(regions.df$name)
  
  replicates.df = readr::read_tsv(opt$replicates)
  replicates.df$replicate = as.character(replicates.df$replicate)
  # Check that each replicate has a distinct name; otherwise it will
  # be impossible to distinguish replicates in the output files
  for (region_name in region_names) {
    df = replicates.df %>% dplyr::filter(name == region_name)
    if (any(duplicated(df$replicate))) {
      stop(sprintf("Region %s: not all replicates have distinct names. It will be impossible to distinguish results in output files that relate to specific replicates.", region_name))
    }
  }
  
  read_data = NULL
  opt$save_read_data <<- F
  opt$save_read_data <- F
  if (!is.null(opt$read_data)) {
    opt$save_read_data <<- !file.exists(opt$read_data)
    opt$save_read_data <- !file.exists(opt$read_data)
    if (file.exists(opt$read_data)) {
      cat(sprintf("Loading read data from file: %s", opt$read_data))
      read_data = readRDS(opt$read_data)
    }
  }
  
  all_region_plots = list()
  del_results = list()
  grep_results = list()
  region_index = 1
  #profvis({
  for (region_name in region_names) {
    df = regions.df %>% dplyr::filter(name == region_name)
    if (nrow(df) > 1) {
      stop("Error: multiple regions with the same region name (%s). Each input line for the --regions file should have a unique name.")
    } else if (nrow(df) == 0) {
      stop("Internal error.")
    }
    cur_replicates.df = replicates.df %>% dplyr::filter(name == region_name)
    if (nrow(cur_replicates.df) == 0) {
      warning(sprintf("No replicates found for region %s", region_name))
      next
    }
    cur_region = as.list(df[1,])
    
    result = doFullRegionAnalysis(locus_name = region_name, region = cur_region, replicates.df = cur_replicates.df, read_data)
    result$region = region_name
    result$replicate.df = cur_replicates.df
    
    grep_results[[region_index]] = result$grep_res
    del_results[[region_index]] = result$del_res
    all_region_plots[[region_index]] = list()
    if (!is.null(result$grep_res)) {
      all_region_plots[[region_index]] = list(grep_stats = result$grep_res$p.stats)
    }
    if (!is.null(result$grep_res)) {
      all_region_plots[[region_index]] = c(all_region_plots[[region_index]], result$del_res$plot_list)
    }
    region_index = region_index + 1
  }
  #})
  
  # Also write out a summary of stats per region, rather than per replicate
  #hdr_prop_error_pval = lapply(del_results, function(res) res$stats.df)
  getListItemField = function(l, fieldList) {
    for (field in fieldList) {
      l = l[[field]]
      if (is.null(l)) {
        break
      }
    }
    if (is.null(l)) {NA} else {l}
  }
  
  grep.summary.df = NULL
  if (opt$grep_analysis) {
    grep_stats.df = bind_rows(lapply(grep_results, function(res) res$stats.df))
    fname = sprintf("%s.grep_analysis.replicate_stats.tsv", opt$out)
    write.table(grep_stats.df, fname, quote=F, row.names=F, col.names=T, na="", sep="\t")
    
    grep.summary.df = data.frame(
      name = region_names,
      hdr.rate          = sapply(grep_results, FUN = function(res) getListItemField(res, c("stats.summary", "mean_hdr_rate"))),
      hdr.pval          = sapply(grep_results, FUN = function(res) getListItemField(res, c("stats.summary", "hdr_wt.res", "pval"))),
      hdr.effect        = sapply(grep_results, FUN = function(res) getListItemField(res, c("stats.summary", "hdr_wt.res", "effect"))),
      hdr.effect_sd     = sapply(grep_results, FUN = function(res) getListItemField(res, c("stats.summary", "hdr_wt.res", "effect_sd"))),
      hdr.effect_confint_lo = sapply(grep_results, FUN = function(res) getListItemField(res, c("stats.summary", "hdr_wt.res", "effect_confint_lo"))),
      hdr.effect_confint_hi = sapply(grep_results, FUN = function(res) getListItemField(res, c("stats.summary", "hdr_wt.res", "effect_confint_hi")))
    )
    fname = sprintf("%s.grep_analysis.region_stats.tsv", opt$out)
    write.table(grep.summary.df, fname, quote=F, row.names=F, col.names=T, sep="\t")
  }
  
  if (opt$deletion_analysis) {
    # Save read_data object if the option is set
    if (opt$save_read_data) {
      read_data_list = unlist(lapply(del_results, FUN = function(x) x$read_data),
                              recursive = F)
      datadir = base::dirname(opt$read_data)
      if (!dir.exists(datadir)) {
        dir.create(datadir)
      }
      saveRDS(read_data_list, file=opt$read_data)
    }
    
    if (opt$allele_profile) {
      df = bind_rows(lapply(del_results, function(res) res$wt_hdr.df))
      fname = sprintf("%s.mismatch_profile.tsv", opt$out)
      write.table(df, fname, quote=F, row.names=F, col.names=T, na="", sep="\t")
    }
    if (!opt$no_udp_profile) {
      df = bind_rows(lapply(del_results, function(res) res$replicate.udp.df))
      fname = sprintf("%s.replicate_udps.tsv", opt$out)
      write.table(df, fname, quote=F, row.names=F, col.names=T, na="", sep="\t")
      
      df = bind_rows(lapply(del_results, function(res) res$merged.udp.df))
      fname = sprintf("%s.merged_udps.tsv", opt$out)
      write.table(df, fname, quote=F, row.names=F, col.names=T, na="", sep="\t")
    }
    if (!opt$no_uns_profile) {
      df = bind_rows(lapply(del_results, function(res) res$uns.df))
      fname = sprintf("%s.uns.tsv", opt$out)
      write.table(df, fname, quote=F, row.names=F, col.names=T, na="", sep="\t")
    }
    if (!opt$no_site_profile) {
      df = bind_rows(lapply(del_results, function(res) res$site.profiles.df))
      fname = sprintf("%s.site_profiles.tsv", opt$out)
      write.table(df, fname, quote=F, row.names=F, col.names=T, na="", sep="\t")
    }
  }
  
  stats_plot_del = NULL
  del.summary.df = NULL
  if (!opt$no_stats & opt$deletion_analysis) {
    stats.df = bind_rows(lapply(del_results, function(res) res$stats.df))
    fname = sprintf("%s.replicate_stats.tsv", opt$out)
    write.table(stats.df, fname, quote=F, row.names=F, col.names=T, na="", sep="\t")
    
    del.summary.df = data.frame(
      name = region_names,
      hdr.rate          = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "mean_hdr_rate"))),
      del.rate          = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "mean_del_rate"))),
      editing.rate      = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "mean_edit_rate"))),
      hdr.pval          = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "hdr_wt.res", "pval"))),
      hdr.effect        = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "hdr_wt.res", "effect"))),
      hdr.effect_sd     = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "hdr_wt.res", "effect_sd"))),
      hdr.effect_confint_lo = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "hdr_wt.res", "effect_confint_lo"))),
      hdr.effect_confint_hi = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "hdr_wt.res", "effect_confint_hi"))),
      del.pval          = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.res", "pval"))),
      del.effect        = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.res", "effect"))),
      del.effect_sd     = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.res", "effect_sd"))),
      del.effect_confint_lo = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.res", "effect_confint_lo"))),
      del.effect_confint_hi = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.res", "effect_confint_hi"))),
      del.2bp.pval      = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.2bp.res", "pval"))),
      del.2bp.effect    = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.2bp.res", "effect"))),
      del.2bp.effect_sd = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.2bp.res", "effect_sd"))),
      del.2bp.effect_confint_lo = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.2bp.res", "effect_confint_lo"))),
      del.2bp.effect_confint_hi = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.2bp.res", "effect_confint_hi"))),
      del.10bp.pval      = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.10bp.res", "pval"))),
      del.10bp.effect    = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.10bp.res", "effect"))),
      del.10bp.effect_sd = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.10bp.res", "effect_sd"))),
      del.10bp.effect_confint_lo = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.10bp.res", "effect_confint_lo"))),
      del.10bp.effect_confint_hi = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.10bp.res", "effect_confint_hi"))),
      del.20bp.pval      = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.20bp.res", "pval"))),
      del.20bp.effect    = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.20bp.res", "effect"))),
      del.20bp.effect_sd = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.20bp.res", "effect_sd"))),
      del.20bp.effect_confint_lo = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.20bp.res", "effect_confint_lo"))),
      del.20bp.effect_confint_hi = sapply(del_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.20bp.res", "effect_confint_hi")))
    )
    fname = sprintf("%s.region_stats.tsv", opt$out)
    write.table(del.summary.df, fname, quote=F, row.names=F, col.names=T, sep="\t")
  }

  stats_summary_plot = statsSummaryPlot(grep.summary.df, del.summary.df)
  
  #"regions", "replicates", "out", "read_data", 
  settings.df = data.frame(setting=names(opt_in), value=as.character(opt_in)) %>%
    filter(setting %in% c("grep_analysis", "grep_window", "grep_window_left", "grep_window_right",
                          "deletion_analysis", "minMapQ", "max_mismatch_frac", "viewing_window", "editing_window", "min_window_overlap",
                          "exclude_multiple_deletions", "exclude_nonspanning_reads", "exclude_nonspanning_deletions",
                          "ratio_to_total_reads", "qc_plot_max_udps", "qc_plot_min_udp_fraction", "qc_plot_exclude_wt",
                          "del_span_start", "del_span_end", "uns_plot_min_gDNA", "uns_plot_min_cDNA", "uns_plot_max_udps"))
  ggThemeBlank = theme_bw() + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.border = element_blank())
  myTableTheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.70)),
                                            colhead = list(fg_params=list(cex = 0.8)),
                                            rowhead = list(fg_params=list(cex = 0.8)))
  settingsPlot = ggplot(data.frame(x=1:10, y=1:10), aes(x,y)) + geom_blank() + ggThemeBlank +
    annotation_custom(tableGrob(settings.df, theme = myTableTheme, rows = NULL), xmin=1, xmax=9, ymin=1, ymax=10) +
    ggtitle("CRISPR editing analysis settings")
  
  varcomp_plots = NULL
  if (opt$variance_analysis & opt$deletion_analysis) {
    varcomp_plots = list()
    # Do variance components analysis across all regions together
    all_replicates.udp.df = bind_rows(lapply(del_results, FUN = function(res) getListItemField(res, c("replicate.udp.df"))))
    #if (length(unique(all_replicates.udp.df$name)) > 1) {
      varcomp = getVarianceComponents(replicate.udp.df = all_replicates.udp.df,
                                      replicates.df = replicates.df %>% dplyr::filter(name %in% regions.df$name),
                                      method = "read_fraction",
                                      min_udp_total_count = opt$variance_analysis_min_count,
                                      min_udp_fraction = opt$variance_analysis_min_fraction)
      varcomp_plots[[1]] = getVarianceComponentsPlots(varcomp$vp.cDNA, varcomp$vp.gDNA,
                                                      min_udp_total_count = opt$variance_analysis_min_count,
                                                      min_udp_fraction = opt$variance_analysis_min_fraction,
                                                      plot_title = "Variance components across regions")
      fname = sprintf("%s.variance_components.tsv", opt$out)
      vp.df = rbind(varcomp$vp.cDNA %>% mutate(type = "cDNA"), varcomp$vp.gDNA %>% mutate(type = "gDNA")) %>%
        select(name, udp, type, everything()) %>%
        arrange(name, udp, -frac)
      write.table(vp.df, fname, quote=F, row.names=F, col.names=T, sep="\t")
      
      # Get variance components for subsets of UDPs defined by their read fraction
      varcomp_plots[[2]] = getVarianceComponentsPlots(varcomp$vp.cDNA %>% filter(frac < 0.005),
                                                      varcomp$vp.gDNA %>% filter(frac < 0.005),
                                                      min_udp_total_count = opt$variance_analysis_min_count,
                                                      min_udp_fraction = opt$variance_analysis_min_fraction,
                                                      plot_title = "Variance components across regions: fraction < 0.5%")
      varcomp_plots[[3]] = getVarianceComponentsPlots(varcomp$vp.cDNA %>% filter(frac >= 0.005, frac < 0.02),
                                                      varcomp$vp.gDNA %>% filter(frac >= 0.005, frac < 0.02),
                                                      min_udp_total_count = opt$variance_analysis_min_count,
                                                      min_udp_fraction = opt$variance_analysis_min_fraction,
                                                      plot_title = "Variance components across regions: fraction 0.5% - 2%")
      varcomp_plots[[4]] = getVarianceComponentsPlots(varcomp$vp.cDNA %>% filter(frac >= 0.02),
                                                      varcomp$vp.gDNA %>% filter(frac >= 0.02),
                                                      min_udp_total_count = opt$variance_analysis_min_count,
                                                      min_udp_fraction = opt$variance_analysis_min_fraction,
                                                      plot_title = "Variance components across regions: fraction > 2%")
    #}
  }
  
  
  fname = sprintf("%s.plots.pdf", opt$out)
  pdf(file = fname, width = opt$plot_width, height = opt$plot_height)
  print(settingsPlot)
  if (!is.null(stats_summary_plot)) { print(stats_summary_plot) }
  print(all_region_plots)
  if (!is.null(varcomp_plots)) { print(varcomp_plots) }
  dev.off()
}


statsSummaryPlot = function(grep.summary.df, stats.summary.df) {
  getSignificanceStr = function(pval) {
    if (is.na(pval) | pval >= 0.01) { "p >= 0.01" }
    else if (pval < 0.001) { "p < 0.001" }
    else { "p < 0.01" }
  }
  
  p.grep.effect = NULL
  p.stats.effect = NULL
  effect_size_theme = theme_bw(10) + theme(axis.text.x = element_blank(),
                                           legend.title = element_blank(),
                                           axis.title.x = element_blank(),
                                           plot.margin = unit(c(0.1, 0, 0.1 ,1), "cm"))
  if (!is.null(grep.summary.df)) {
    exp_names = sapply(grep.summary.df$name, FUN = function(s) strsplit(s, ",", T)[[1]][1])
    if (any(duplicated(exp_names))) {
      exp_names = grep.summary.df$name
    }
    grep.summary.df$name = factor(as.character(exp_names), levels=exp_names)
    grep.summary.df$hdr_significance = factor(sapply(grep.summary.df$hdr.pval, FUN = getSignificanceStr), levels=c("p >= 0.01", "p < 0.01", "p < 0.001"))
    
    p.grep.effect = ggplot(grep.summary.df, aes(x=name, y=hdr.effect, fill=hdr_significance)) +
      geom_bar(stat = "identity", width=0.5) + 
      geom_errorbar(aes(ymin = hdr.effect_confint_lo, ymax = hdr.effect_confint_hi),
                    width = 0.2, col = "grey30") +
      geom_hline(yintercept = 1, col = "red") +
      effect_size_theme +
      ylab("HDR effect size") + ggtitle("HDR effect size - grep analysis") +
      scale_fill_manual(values=c(`p >= 0.01`="grey70", `p < 0.01`="cornflowerblue", `p < 0.001`="red3")) +
      coord_cartesian(ylim = c(0, max(1, max(grep.summary.df$hdr.effect * 1.05, na.rm = T))))
    
    hdr.df = grep.summary.df
  }
  
  if (!is.null(stats.summary.df)) {
    exp_names = sapply(stats.summary.df$name, FUN = function(s) strsplit(s, ",", T)[[1]][1])
    if (any(duplicated(exp_names))) {
      exp_names = stats.summary.df$name
    }
    stats.summary.df$name = factor(as.character(exp_names), levels=exp_names)
    stats.summary.df$hdr_significance = factor(sapply(stats.summary.df$hdr.pval, FUN = getSignificanceStr), levels=c("p >= 0.01", "p < 0.01", "p < 0.001"))
    
    p.stats.effect = ggplot(stats.summary.df, aes(x=name, y=hdr.effect, fill=hdr_significance)) +
      geom_bar(stat = "identity", width=0.5) + 
      geom_errorbar(aes(ymin = hdr.effect_confint_lo, ymax = hdr.effect_confint_hi),
                    width = 0.2, col = "grey30") +
      geom_hline(yintercept = 1, col = "red") +
      effect_size_theme +
      ylab("HDR effect size") + ggtitle("HDR effect size - alignment deletion analysis") +
      scale_fill_manual(values=c(`p >= 0.01`="grey70", `p < 0.01`="cornflowerblue", `p < 0.001`="red3")) +
      coord_cartesian(ylim = c(0, max(1, max(stats.summary.df$hdr.effect * 1.05, na.rm = T))))
    
    stats.summary.df$del_significance = factor(sapply(stats.summary.df$del.pval, FUN = getSignificanceStr), levels=c("p >= 0.01", "p < 0.01", "p < 0.001"))
    p.stats.del = ggplot(stats.summary.df, aes(x=name, y=del.effect, fill=del_significance)) +
      geom_bar(stat = "identity", width=0.5) + 
      geom_errorbar(aes(ymin = del.effect_confint_lo, ymax = del.effect_confint_hi),
                    width = 0.2, col = "grey30") +
      geom_hline(yintercept = 1, col = "red") +
      effect_size_theme +
      ylab("Del effect size") + ggtitle("Deletion effect size") +
      scale_fill_manual(values=c(`p >= 0.01`="grey70", `p < 0.01`="cornflowerblue", `p < 0.001`="red3")) +
      coord_cartesian(ylim = c(0, min(3, max(stats.summary.df$del.effect * 1.2, na.rm = T))))
    
    plot.df = stats.summary.df %>% dplyr::select(name, HDR=hdr.rate, NHEJ=del.rate) %>%
      tidyr::gather(key = "type", value = "value", -name)
    plot.df$type = factor(as.character(plot.df$type), levels = c("NHEJ", "HDR"))
    p.stats.editing = ggplot(plot.df, aes(x=name, y=value*100, fill=type)) +
      geom_bar(stat = "identity", position = position_stack(), width=0.5) + 
      theme_bw(10) + theme(axis.text.x = element_blank(),
                           legend.title = element_blank(),
                           axis.title.x = element_blank(),
                           plot.margin = unit(c(0.1, 0, 0.1 ,1), "cm")) +
      ylab("% editing") + ggtitle("Editing rates") +
      scale_fill_manual(values=c(HDR="darkorange", NHEJ="cornflowerblue"))
    
    hdr.df = stats.summary.df
  }
  
  hdr.df$type = "HDR"
  p.stats.hdr = ggplot(hdr.df, aes(x=name, y=hdr.rate * 100, fill=type)) +
    geom_bar(stat = "identity", position = position_stack(), width=0.5) + 
    theme_bw(10) + theme(axis.text.x = element_text(angle = 37, hjust = 1, size=7),
                         legend.title = element_blank(),
                         axis.title.x = element_blank(),
                         plot.margin = unit(c(0.1, 0, 0.1 ,1), "cm")) +
    ylab("% HDR") + ggtitle("HDR rates") +
    scale_fill_manual(values=c(HDR="darkorange", NHEJ="cornflowerblue"))
  
  p.title = ggdraw() + draw_label("Experiment summary", fontface='bold')
  if (!is.null(p.grep.effect) & !is.null(p.stats.effect)) {
    p.res = egg::ggarrange(p.title, p.grep.effect, p.stats.effect, p.stats.del, p.stats.editing, p.stats.hdr, ncol=1, heights=c(1.2,2.4,2.4,2.4,2.4,2.4), draw = F)
  } else if (!is.null(p.grep.effect)) {
    p.res = egg::ggarrange(p.title, p.grep.effect, p.stats.del, p.stats.editing, p.stats.hdr, ncol=1, heights=c(1.2,3,3,3,3), draw = F)
  } else {
    p.res = egg::ggarrange(p.title, p.stats.effect, p.stats.del, p.stats.editing, p.stats.hdr, ncol=1, heights=c(1.2,3,3,3,3), draw = F)
  }
  p.res
}


doFullRegionAnalysis = function(locus_name, region, replicates.df, read_data = NULL)
{
  if (nrow(replicates.df) < 1) {
    stop(sprintf("doFullRegionAnalysis: No replicates specified for region %s.", region$name))
  }
  # stats_list = list()
  sites = list(start = region$start,
               end = region$end,
               highlight_site = region$highlight_site,
               cut_site = region$cut_site)
  rel_sites = getRelativeSites(sites, opt$viewing_window)
  
  replicate_grep_analyses = list()
  replicate_del_analyses = list()
  for (i in 1:nrow(replicates.df)) {
    bam_file = replicates.df[i,]$bam
    replicate = replicates.df[i,]$replicate
    type = replicates.df[i,]$type
    hdr_profile = str_to_upper(region$hdr_allele_profile)
    wt_profile = str_to_upper(region$wt_allele_profile)
    ref_sequence = str_to_upper(region$ref_sequence)
    
    cat(sprintf("\n\nAnalysing region %s, replicate %s, file %s, %s:%d-%d, highlight site %d, cut site %d\n",
                locus_name, replicate, bam_file, region$sequence_name, sites$start, sites$end, sites$highlight_site, sites$cut_site))
    cat(sprintf("HDR allele: %s\n", hdr_profile))
    cat(sprintf("WT allele:  %s\n", wt_profile))
    cat(sprintf("REF sequence: %s\n", ref_sequence))
    
    sam_reads = NULL
    if (opt$grep_analysis | (opt$deletion_analysis & !opt$grep_analysis & !is.null(read_data))) {
      sam_reads = getRegionReadsFromBam(locus_name, replicate, bam_file, region$sequence_name, sites$start, sites$end, ref_sequence)
      if (length(sam_reads) <= 0) {
        cat(sprintf("\nERROR: No reads retrieved at the specified region for replicate %s from file %s\n", replicate, bam_file))
        stop()
      }
    }
    
    if (opt$grep_analysis) {
      result = doReplicateGrepAnalysis(sam_reads = sam_reads,
                                       hdr_profile = hdr_profile,
                                       wt_profile = wt_profile,
                                       ref_sequence = ref_sequence)
      result$name = locus_name
      result$replicate = replicate
      result$type = type
      result$num_reads = length(sam_reads)
      replicate_grep_analyses[[i]] = result
    }
    
    if (opt$deletion_analysis) {
      result = doReplicateDeletionAnalysis(name = locus_name, replicate = replicate, type = type,
                                           sam_reads = sam_reads, sites = sites,
                                           hdr_profile = hdr_profile,
                                           wt_profile = wt_profile,
                                           ref_sequence = ref_sequence,
                                           read_data_all = read_data)
      result$name = locus_name
      result$replicate = replicate
      result$type = type
      replicate_del_analyses[[i]] = result
    }
  }
  
  grep_res = del_res = NULL
  if (opt$grep_analysis) {
    counts.df = bind_rows(replicate_grep_analyses) %>%
      dplyr::select(name, replicate, type, num_reads, everything())
    grep_res = doRegionGrepAnalysis(counts.df, replicates.df)
  }
  if (opt$deletion_analysis) {
    del_res = doRegionDeletionAnalysis(replicate_del_analyses, replicates.df, rel_sites, locus_name)
  }
  return(list(grep_res = grep_res, del_res = del_res))
}


doRegionGrepAnalysis = function(counts.df, replicates.df) {
  stats.df = counts.df %>% rename(num_hdr_reads = hdr_read_count, num_wt_reads = wt_read_count)
  stats.df$HDR_WT_ratio = stats.df$num_hdr_reads / stats.df$num_wt_reads
  stats.df$HDR_rate = stats.df$num_hdr_reads / stats.df$num_reads
  stats.df$WT_rate = stats.df$num_wt_reads / stats.df$num_reads
  
  # Summary of stats from the propagation of errors method
  hdr_ratio_res = NULL
  method = ""
  stats.gDNA = stats.df %>% filter(type == "gDNA")
  stats.cDNA = stats.df %>% filter(type == "cDNA")

  if (nrow(stats.gDNA) < 2 | nrow(stats.cDNA) < 2) {
    summary.left = sprintf("Unable to calculate stats with < 2 replicates")
    summary.right = ""
  } else {
    denom = "num_wt_reads"
    ratioTo = "WT"
    if (opts$ratio_to_total_reads) {
      denom = "num_reads"
      ratioTo = "N_tot"
    }
    
    confIntervalString = function(ratioRes) {
      sprintf("95%% CI: (%.3g, %.3g)", ratioRes$effect_confint_lo, ratioRes$effect_confint_hi)
    }
    hdr_ratio_res = getUDPRatioEstimate(stats.df, replicates.df, numerator = "num_hdr_reads", denominator = denom,
                                        batchCol = opt$batch_col, randomEffectsCols = opt$random_effects_cols)
    hdr.conf.interval.str = confIntervalString(hdr_ratio_res)
    hdr.summary = sprintf("cDNA:gDNA ratio (HDR/%s): %.3g\n%s,    p = %.3g",
                          ratioTo, hdr_ratio_res$effect, hdr.conf.interval.str, hdr_ratio_res$pval)
    method = sprintf("Method: %s", hdr_ratio_res$method)
    
    hdr.rate.str = sprintf("Mean HDR rate gDNA: %.2g%%,  cDNA: %.2g%%", mean(stats.gDNA$HDR_rate) * 100, mean(stats.cDNA$HDR_rate) * 100)
    wt.rate.str = sprintf("Mean WT rate gDNA: %.2g%%,  cDNA: %.2g%%", mean(stats.gDNA$WT_rate) * 100, mean(stats.cDNA$WT_rate) * 100)
    summary.left = paste(hdr.rate.str, wt.rate.str, sep = "\n")
    summary.right = paste(method, hdr.summary, sep = "\n")
  }
  
  # Convert to strings for nice printing (with 3 significant digits)
  stats.plot.df = stats.df %>%
    dplyr::select(replicate, type, num_reads, "HDR reads" = num_hdr_reads, "WT reads" = num_wt_reads, HDR_WT_ratio, HDR_rate, WT_rate)
  stats.plot.df$HDR_WT_ratio = sprintf("%.3g", stats.plot.df$HDR_WT_ratio)
  stats.plot.df$HDR_rate = sprintf("%.2f%%", 100 * stats.plot.df$HDR_rate)
  stats.plot.df$WT_rate = sprintf("%.2f%%", 100 * stats.plot.df$WT_rate)
  ggThemeBlank = theme_bw() + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.border = element_blank())
  
  mycex = min(0.70, 0.75 * 11 / nrow(stats.plot.df))
  myTableTheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = mycex)),
                                            colhead = list(fg_params=list(cex = mycex)),
                                            rowhead = list(fg_params=list(cex = mycex)),
                                            padding = unit(c(2, 4), "mm"))
  # p.stats = ggplot(data.frame(x=1:10, y=1:10), aes(x,y)) + geom_blank() + ggThemeBlank +
  #   annotation_custom(tableGrob(t(stats.plot.df), theme = myTableTheme), xmin=1.2, xmax=10, ymin=1, ymax=8) +
  #   annotate("text", x=5, y=10, label = sprintf("%s grep summary", stats.df$name[1]), vjust = 1, fontface = 2, size = 5) +
  #   annotate("text", x=1, y=9.3, label = summary.left.prop, vjust = 1, hjust = 0, size = 3.1)
  plot_title = sprintf("%s grep summary", stats.df$name[1])
  p.stats = ggplot(data.frame(x=1:10, y=1:10), aes(x,y)) + geom_blank() + ggThemeBlank +
    annotation_custom(tableGrob(t(stats.plot.df), theme = myTableTheme), xmin=1.2, xmax=10, ymin=1, ymax=8.5) +
    annotate("text", x=1, y=10, label = summary.left, vjust = 1, hjust = 0, size = 2.9) +
    annotate("text", x=9.5, y=10, label = summary.right, vjust = 1, hjust = 1, size = 2.9)
  
  color_values = c(`cDNA`="firebrick1", `cDNA outlier`="orange1", `gDNA`="dodgerblue3", `gDNA outlier`="turquoise1")
  p1 = ggplot(stats.df, aes(x=replicate, y=num_reads, fill=type)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%d", num_reads)), size = 2.4, position = position_stack(vjust = 0.5)) +
    theme_bw(10) + theme(axis.text.x=element_blank(), legend.title=element_blank(), axis.title.x = element_blank()) +
    scale_fill_manual(values = color_values) +
    ylab("Number of reads") +
    ggtitle("Number of reads")
  
  p2 = ggplot(stats.df, aes(x=replicate, y=HDR_WT_ratio, fill=type)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.3g", HDR_WT_ratio)), size = 2.4, position = position_stack(vjust = 0.5)) +
    theme_bw(10) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.title=element_blank()) +
    scale_fill_manual(values = color_values, guide=F) +
    ylab("HDR:WT ratio") +
    ggtitle("HDR:WT ratio")
  
  p.title = ggdraw() + draw_label(plot_title, fontface='bold')
  #p.replicate_qc = plot_grid(p.udp_fractions, p.udp_avg_deviation, ncol=1)
  #p.res = plot_grid(p.title, p.replicate_qc, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
  p.res = egg::ggarrange(p.title, p.stats, p1, p2, ncol=1, heights=c(0.5, 4.5, 1.5, 1.5), draw = F)
  
  stats.summary = list(mean_hdr_rate = mean(stats.gDNA$HDR_rate),
                       mean_wt_rate = mean(stats.gDNA$WT_rate),
                       hdr_wt.res = hdr_ratio_res,
                       error_prop.summary = summary.left.prop)
  
  result_list = list(p.stats = p.res,
                     stats.df = stats.df,
                     stats.summary = stats.summary)
}


doReplicateGrepAnalysis = function(sam_reads, hdr_profile, wt_profile, ref_sequence) {
  # Get sequences from full sam_read entries
  sam_seqs = vector(mode = "character", length = length(sam_reads)) # Preallocate a vector of discarded reads
  for (i in 1:length(sam_reads)) {
    line_split = strsplit(sam_reads[i], "\t", fixed = T)[[1]]
    sam_seqs[i] = line_split[10]
  }
  hdr_profile_chars = strsplit(hdr_profile, "")[[1]]
  wt_profile_chars = strsplit(wt_profile, "")[[1]]
  ref_seq_chars = strsplit(ref_sequence, "")[[1]]
  
  getGrepSequence = function(profile_chars, ref_chars) {
    profile_is_letter = isDNALetter(profile_chars)
    char_positions = which(profile_is_letter)
    minpos = max(1, min(char_positions) - opt$grep_window_left)
    maxpos = min(length(ref_chars), max(char_positions) + opt$grep_window_right)
    grep_chars = ref_chars
    grep_chars[profile_is_letter] = profile_chars[profile_is_letter]
    paste(grep_chars[minpos:maxpos], collapse = "")
  }
  hdr_seq = getGrepSequence(hdr_profile_chars, ref_seq_chars)
  wt_seq = getGrepSequence(wt_profile_chars, ref_seq_chars)
  
  hdr_read_count = sum(sapply(sam_seqs, FUN = function(s) grepl(hdr_seq, s, fixed=T)))
  wt_read_count = sum(sapply(sam_seqs, FUN = function(s) grepl(wt_seq, s, fixed=T)))
  return(list(hdr_read_count = hdr_read_count, wt_read_count = wt_read_count))
}


doRegionDeletionAnalysis = function(replicate_del_analyses, replicates.df, rel_sites, locus_name) {
  site.profiles.df = NULL
  wt_hdr.df = NULL
  
  if (opt$allele_profile) {
    replicate_wt_hdr = lapply(replicate_del_analyses, FUN = function(res) res$wt_hdr.df)
    wt_hdr.df = bind_rows(replicate_wt_hdr)
  }
  
  if (!opt$no_site_profile) {
    # Make a table of the "site profiles" - combinations of alleles at sites of interest
    replicate_site_profiles = lapply(replicate_del_analyses, FUN = function(res) res$sites.profile.df)
    site.profiles.df = bind_rows(replicate_site_profiles)
    site.profiles.df$replicate = as.character(site.profiles.df$replicate)
    # Combine replicates and add merged counts as a replicate named "all"
    site.profile.merged = site.profiles.df %>% group_by(sites_profile, type) %>%
      dplyr::summarise(name = first(name),
                       replicate = "all",
                       count = sum(count)) %>%
      dplyr::select(name, replicate, type, sites_profile, count)
    site.profiles.df = bind_rows(site.profiles.df, site.profile.merged) %>%
      dplyr::arrange(replicate, type, sites_profile)
  }
  
  stats_res = getFullReplicateStats(replicate_del_analyses, rel_sites, replicates.df)
  
  # Merge together results from all replicates
  replicate_udps = lapply(replicate_del_analyses, FUN = function(res) res$udp.df)
  replicate.udp.df = bind_rows(replicate_udps)
  
  replicate_plots = list()
  if (!opt$no_replicate_plots) {
    for (i in 1:nrow(replicates.df)) {
      cur_replicate = replicates.df[i,]$replicate
      cur_type = replicates.df[i,]$type
      plot_title = sprintf("%s replicate %s, %s", locus_name, cur_replicate, cur_type)
      udp.df = replicate.udp.df %>% dplyr::filter(name == locus_name, replicate == cur_replicate, type == cur_type)
      replicate_plots[[i]] = getUDPPlot(udp.df, plot_title, rel_sites)
    }
  }
  
  # Determine which UDPs are shared between cDNA and gDNA
  getSharing = function(types) {
    if ("cDNA" %in% types) {
      if ("gDNA" %in% types) {
        "both"
      } else {
        "cDNA only"
      }
    } else {
      if (!"gDNA" %in% types) { stop("Unexpected: neither cDNA nor gDNA among replicate types.") }
      "gDNA only"
    }
  }

  gDNACounts.df = replicate.udp.df %>% filter(type == "gDNA") %>% group_by(udp) %>%
    summarise(udpcount_gDNA = sum(num_reads))
  cDNACounts.df = replicate.udp.df %>% filter(type == "cDNA") %>% group_by(udp) %>%
    summarise(udpcount_cDNA = sum(num_reads))
  replicate.udp.df = replicate.udp.df %>%
    left_join(gDNACounts.df, by=c("udp")) %>%
    left_join(cDNACounts.df, by=c("udp"))
  threshold = 10
  replicate.udp.df$udp_sharing = "unclear"
  replicate.udp.df$udp_sharing[replicate.udp.df$udpcount_gDNA >= threshold & replicate.udp.df$udpcount_cDNA >= threshold] = "both"
  replicate.udp.df$udp_sharing[replicate.udp.df$udpcount_gDNA > 0 & is.na(replicate.udp.df$udpcount_cDNA)] = "gDNA only"
  replicate.udp.df$udp_sharing[replicate.udp.df$udpcount_cDNA > 0 & is.na(replicate.udp.df$udpcount_gDNA)] = "cDNA only"

  merged.udp.df = summarise(replicate.udp.df %>% group_by(name, type, udp),
                            num_reads = sum(num_reads),
                            is_hdr_allele = first(is_hdr_allele),
                            is_wt_allele = first(is_wt_allele),
                            has_any_deletion = first(has_any_deletion),
                            has_custom_deletion = first(has_custom_deletion),
                            has_crispr_deletion = first(has_crispr_deletion),
                            deletion_start = first(deletion_start),
                            deletion_end = first(deletion_end),
                            deletion2_start = first(deletion2_start),
                            deletion2_end = first(deletion2_end),
                            avg_seq_length = mean(avg_seq_length),
                            avg_mismatch_count = mean(avg_mismatch_count),
                            udp_sharing = first(udp_sharing)) %>%
    arrange(-num_reads) %>% ungroup()
  
  # Merged UDPs plot
  region_title = sprintf("%s deletion alleles", locus_name)
  p.udp_gDNA = getUDPPlot(merged.udp.df %>% dplyr::filter(type == "gDNA"), plot_title="gDNA", sites=rel_sites)
  p.udp_cDNA = getUDPPlot(merged.udp.df %>% dplyr::filter(type == "cDNA"), plot_title="cDNA", sites=rel_sites)
  p.title = ggdraw() + draw_label(region_title, fontface='bold')
  p.cDNA_gDNA = plot_grid(p.udp_gDNA, p.udp_cDNA, nrow=1)
  p.merged_UDP = plot_grid(p.title, p.cDNA_gDNA, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
  
  # Deletion profiles plot - averaged + separate replicates
  delprofile.udp.df = replicate.udp.df
  delprofile.udp.df$count_udp = 1
  if (opt$custom_del_span) {
    delprofile.udp.df = replicate.udp.df %>% dplyr::mutate(count_udp = ifelse(has_custom_deletion, 1, 0))
  }
  if (opt$use_cdna_dels_only) {
    delprofile.udp.df = delprofile.udp.df %>% dplyr::mutate(count_udp = count_udp & ifelse(udp_sharing == "both", 1, 0))
  }
  p.del_profile_pct = getDeletionProfilePlot(delprofile.udp.df,
                                                rel_sites,
                                                plot_title = "Relative to all reads",
                                                show_average = T, show_replicates = T, ratioToWT = F)
  p.del_profile_pct = p.del_profile_pct + theme(axis.title.x=element_blank())
  p.del_profile_wtratio = getDeletionProfilePlot(delprofile.udp.df,
                                                   rel_sites,
                                                   plot_title = "Relative to WT",
                                                   show_average = T, show_replicates = T, ratioToWT = T)
  del_profile_title = sprintf("%s deletion profile %s", locus_name, ifelse(opt$custom_del_span, "\ncustom dels only", ""))
  p.del_profile = egg::ggarrange(ggdraw() + draw_label(del_profile_title, fontface='bold'),
                                 p.del_profile_pct, p.del_profile_wtratio, 
                                 ncol=1, heights=c(0.1, 0.5, 0.5), draw = F)
  
  replicate_qc_plots = list(p1 = NULL, p2 = NULL)
  if (opt$replicate_qc_plots) {
    qc_metrics = getReplicateQCMetrics(stats.df = stats_res$stats.df, replicate.udp.df = replicate.udp.df,
                                       max_udps = opt$qc_plot_max_udps, min_avg_udp_fraction = opt$qc_plot_min_udp_fraction,
                                       exclude_wt = opt$qc_plot_exclude_wt)
    stats_res$stats.df = qc_metrics$stats.df
    replicate_qc_plots = getReplicateQCPlots(stats.df = stats_res$stats.df,
                                             replicate.udp.fractions.df = qc_metrics$replicate.udp.fractions.df,
                                             min_avg_udp_fraction = opt$qc_plot_min_udp_fraction,
                                             region_title = locus_name)
  }
  
  # UNS plot
  p.merged_UNS = getUNSPlot(replicate.udp.df, replicates.df = replicates.df,
                            sites = rel_sites, plot_title = region_title,
                            min_gDNA_count = opt$uns_plot_min_gDNA, min_cDNA_count = opt$uns_plot_min_cDNA,
                            max_udps = opt$uns_plot_max_udps)
  
  uns_res = getUNSData(replicate.udp.df, replicates.df = replicates.df,
                       sites = rel_sites, region_name = region_title,
                       min_gDNA_count = opt$uns_plot_min_gDNA, min_cDNA_count = opt$uns_plot_min_cDNA)
  uns.df = NULL
  uns.replicates.df = NULL
  if (!is.null(uns_res)) {
    uns.df = uns_res$udp.dels.df %>% dplyr::mutate(name = locus_name) %>% dplyr::select(name, everything())
    uns.replicates.df = uns_res$replicate.dels.df %>% dplyr::mutate(name = locus_name) %>% dplyr::select(name, everything())
  }
  
  p.variance_components = NULL
  if (opt$variance_analysis) {
    varcomp = getVarianceComponents(replicate.udp.df, replicates.df, method = "read_fraction",
                                    min_udp_total_count = opt$variance_analysis_min_count,
                                    min_udp_fraction = opt$variance_analysis_min_fraction)
    p.variance_components = getVarianceComponentsPlots(varcomp$vp.cDNA, varcomp$vp.gDNA,
                                                       min_udp_total_count = opt$variance_analysis_min_count,
                                                       min_udp_fraction = opt$variance_analysis_min_fraction,
                                                       plot_title = sprintf("%s variance components", replicates.df$name[1]))
  }
  
  p.variance_fit = NULL
  p.power = NULL
  p.replicate_allocation = NULL
  if (opt$power_analysis) {
    power_plots = getPowerPlots(replicate.udp.df, replicates.df, titlestr = locus_name, min_udp_total_count = opt$variance_analysis_min_count)
    p.variance_fit = power_plots$cv_plots
    p.power = power_plots$power
    p.replicate_allocation = power_plots$replicate_allocation
  }
  
  plot_list = list(stats = stats_res$p.stats,
                   merged_udp = p.merged_UDP,
                   merged_del_profile = p.del_profile,
                   replicate_qc_1 = replicate_qc_plots$p1,
                   replicate_qc_2 = replicate_qc_plots$p2,
                   merged_UNS = p.merged_UNS,
                   variance_components = p.variance_components,
                   variance_fit = p.variance_fit,
                   power = p.power,
                   replicate_allocation = p.replicate_allocation,
                   replicate_plots = replicate_plots)
  
  read_data = lapply(replicate_del_analyses, FUN = function(x) x$read_data)
  names(read_data) = paste(replicates.df$name, replicates.df$replicate, replicates.df$type, sep="_")
  
  result_list = list(replicate_list = replicate_del_analyses,
                     plot_list = plot_list,
                     replicate.udp.df = replicate.udp.df,
                     merged.udp.df = merged.udp.df,
                     uns.df = uns.df,
                     uns.replicates.df = uns.replicates.df,
                     site.profiles.df = site.profiles.df,
                     wt_hdr.df = wt_hdr.df,
                     stats.df = stats_res$stats.df,
                     stats.summary = stats_res$stats.summary,
                     read_data = read_data)
  return(result_list)
}

######################################################################################

getFullReplicateStats = function(replicate_del_analyses, rel_sites, replicates.df) {
  stats.df = bind_rows(lapply(replicate_del_analyses, function(res) as.data.frame(res$stats)))
  stats.df$HDR_WT_ratio = (stats.df$num_hdr_reads / stats.df$num_wt_reads)
  stats.df$DEL_WT_ratio = (stats.df$num_deletion_reads / stats.df$num_wt_reads)
  stats.df$HDR_rate = stats.df$num_hdr_reads / stats.df$num_kept_reads
  stats.df$DEL_rate = stats.df$num_deletion_reads / stats.df$num_kept_reads
  stats.df$editing_rate = stats.df$num_edit_reads / stats.df$num_kept_reads
  stats.df$WT_rate = stats.df$num_wt_reads / stats.df$num_kept_reads
  
  # Summary of stats from the propagation of errors method
  hdr_ratio_res = NULL
  del_ratio_res = NULL
  del_ratio_res_2bp = NULL
  del_ratio_res_10bp = NULL
  del_ratio_res_20bp = NULL
  del_ratio_res_custom = NULL
  method = ""
  stats.gDNA = stats.df %>% filter(type == "gDNA")
  stats.cDNA = stats.df %>% filter(type == "cDNA")
  hdr.rate = sprintf("Mean HDR rate gDNA: %.2g%%,  cDNA: %.2g%%", mean(stats.gDNA$HDR_rate) * 100, mean(stats.cDNA$HDR_rate) * 100)
  del.rate = sprintf("Mean DEL rate gDNA: %.2g%%,  cDNA: %.2g%%", mean(stats.gDNA$DEL_rate) * 100, mean(stats.cDNA$DEL_rate) * 100)
  wt.rate = sprintf("Mean WT rate gDNA: %.2g%%,  cDNA: %.2g%%", mean(stats.gDNA$WT_rate) * 100, mean(stats.cDNA$WT_rate) * 100)

  if (nrow(stats.gDNA) < 2 | nrow(stats.cDNA) < 2) {
    summary.left.prop = sprintf("Unable to calculate stats with < 2 replicates")
    summary.right.prop = ""
  } else {
    denom = "num_wt_reads"
    ratioTo = "WT"
    if (opts$ratio_to_total_reads) {
      denom = "num_kept_reads"
      ratioTo = "N_tot"
    }
    
    confIntervalString = function(ratioRes) {
      sprintf("95%% CI: (%.3g, %.3g)", ratioRes$effect_confint_lo, ratioRes$effect_confint_hi)
    }
    hdr_ratio_res = getUDPRatioEstimate(stats.df, replicates.df, numerator = "num_hdr_reads", denominator = denom, batchCol = opt$batch_col, randomEffectsCols = opt$random_effects_cols)
    hdr.conf.interval.str = confIntervalString(hdr_ratio_res)
    hdr.summary.prop = sprintf("cDNA:gDNA ratio (HDR/%s): %.3g\n%s,    p = %.3g",
                               ratioTo, hdr_ratio_res$effect, hdr.conf.interval.str, hdr_ratio_res$pval)
    method = sprintf("\nMethod: %s", hdr_ratio_res$method)

    summary.left.prop = paste(hdr.rate, del.rate, wt.rate, method, hdr.summary.prop, sep = "\n")
    
    del_ratio_res = getUDPRatioEstimate(stats.df, replicates.df, numerator = "num_deletion_reads", denominator = denom, batchCol = opt$batch_col, randomEffectsCols = opt$random_effects_cols)
    del.conf.interval.str = confIntervalString(del_ratio_res)
    if (opt$custom_del_span) {
      del.summary.prop = sprintf("Custom del site %d, span %d-%d\ncDNA:gDNA ratio (DEL/WT): %.3g\n%s,    p = %.3g",
                                 ratioTo, rel_sites$highlight_site, opt$del_span_start, opt$del_span_end,
                                 del_ratio_res$effect, del.conf.interval.str, del_ratio_res$pval)
    } else {
      del.summary.prop = sprintf("cDNA:gDNA ratio (DEL/%s): %.3g\n%s,    p = %.3g",
                                 ratioTo, del_ratio_res$effect, del.conf.interval.str, del_ratio_res$pval)
    }
    
    del_ratio_res_2bp = getUDPRatioEstimate(stats.df, replicates.df, numerator = "num_deletions_2bp_window", denominator = denom, batchCol = opt$batch_col, randomEffectsCols = opt$random_effects_cols)
    del.conf.interval.str = confIntervalString(del_ratio_res_2bp)
    del.summary.prop_2bp = sprintf("cDNA:gDNA ratio (DEL/%s) - 2 bp span: %.3g\n%s,    p = %.3g",
                                   ratioTo, del_ratio_res_2bp$effect, del.conf.interval.str, del_ratio_res_2bp$pval)
    
    del_ratio_res_10bp = getUDPRatioEstimate(stats.df, replicates.df, numerator = "num_deletions_10bp_window", denominator = denom, batchCol = opt$batch_col, randomEffectsCols = opt$random_effects_cols)
    del.conf.interval.str = confIntervalString(del_ratio_res_10bp)
    del.summary.prop_10bp = sprintf("cDNA:gDNA ratio (DEL/%s) - 10 bp span: %.3g\n%s,    p = %.3g",
                                    ratioTo, del_ratio_res_10bp$effect, del.conf.interval.str, del_ratio_res_10bp$pval)
    
    del_ratio_res_20bp = getUDPRatioEstimate(stats.df, replicates.df, numerator = "num_deletions_20bp_window", denominator = denom, batchCol = opt$batch_col, randomEffectsCols = opt$random_effects_cols)
    del.conf.interval.str = confIntervalString(del_ratio_res_20bp)
    del.summary.prop_20bp = sprintf("cDNA:gDNA ratio (DEL/%s) - 20 bp span: %.3g\n%s,    p = %.3g",
                                    ratioTo, del_ratio_res_20bp$effect, del.conf.interval.str, del_ratio_res_20bp$pval)
    
    summary.right.prop = paste(del.summary.prop, del.summary.prop_2bp, del.summary.prop_10bp, sep = "\n")
  }
  
  # Convert to strings for nice printing (with 3 significant digits)
  stats.plot.df = stats.df %>% dplyr::select(replicate, type, num_udps, HDR_WT_ratio, DEL_WT_ratio,
                                             HDR_rate, DEL_rate, editing_rate, WT_rate, num_reads, num_hdr_reads, num_wt_reads,
                                             num_deletion_reads, num_insertion, reads_excluded_for_minoverlap, reads_excluded_for_mismatches,
                                             reads_excluded_nonspanning, reads_excluded_for_multiple_deletions)
  stats.plot.df = stats.plot.df %>% dplyr::rename("HDR reads" = num_hdr_reads,
                                                  "WT reads" = num_wt_reads,
                                                  "Deletion reads" = num_deletion_reads,
                                                  "excluded-insertion" = num_insertion,
                                                  "excluded-minoverlap" = reads_excluded_for_minoverlap,
                                                  "excluded-mismatches" = reads_excluded_for_mismatches,
                                                  "excluded-nonspanning" = reads_excluded_nonspanning,
                                                  "excluded-mult.deletions" = reads_excluded_for_multiple_deletions)
  stats.plot.df$HDR_WT_ratio = sprintf("%.3g", stats.plot.df$HDR_WT_ratio)
  stats.plot.df$DEL_WT_ratio = sprintf("%.3g", stats.plot.df$DEL_WT_ratio)
  stats.plot.df$HDR_rate = sprintf("%.2f%%", 100 * stats.plot.df$HDR_rate)
  stats.plot.df$DEL_rate = sprintf("%.2f%%", 100 * stats.plot.df$DEL_rate)
  stats.plot.df$editing_rate = sprintf("%.2f%%", 100 * stats.plot.df$editing_rate)
  stats.plot.df$WT_rate = sprintf("%.2f%%", 100 * stats.plot.df$WT_rate)
  ggThemeBlank = theme_bw() + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.border = element_blank())
  
  mycex = min(0.70, 0.75 * 11 / nrow(stats.plot.df))
  myTableTheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = mycex)),
                                            colhead = list(fg_params=list(cex = mycex)),
                                            rowhead = list(fg_params=list(cex = mycex)),
                                            padding = unit(c(2, 3), "mm"))
  p.stats = ggplot(data.frame(x=1:10, y=1:10), aes(x,y)) + geom_blank() + ggThemeBlank +
    annotation_custom(tableGrob(t(stats.plot.df), theme = myTableTheme), xmin=1.2, xmax=10, ymin=1, ymax=8) +
    annotate("text", x=5, y=10, label = sprintf("%s analysis summary", stats.df$name[1]), vjust = 1, fontface = 2, size = 5) +
    annotate("text", x=1, y=9.4, label = summary.left.prop, vjust = 1, hjust = 0, size = 2.9) +
    annotate("text", x=9.5, y=9.4, label = summary.right.prop, vjust = 1, hjust = 1, size = 2.9)
  if (!opt$no_stats) {
    p.stats = p.stats + annotate("text", x=1, y=1, label = "Data saved in *.replicate_stats.tsv", vjust = 1, hjust = 0, size = 3)
  }
  
  stats.summary = list(mean_hdr_rate = mean(stats.gDNA$HDR_rate),
                       mean_del_rate = mean(stats.gDNA$DEL_rate),
                       mean_edit_rate = mean(stats.gDNA$editing_rate),
                       hdr_wt.res = hdr_ratio_res,
                       del_wt.res = del_ratio_res,
                       del_wt.2bp.res = del_ratio_res_2bp,
                       del_wt.10bp.res = del_ratio_res_10bp,
                       del_wt.20bp.res = del_ratio_res_20bp,
                       del_wt.custom.res = del_ratio_res_custom,
                       error_prop.summary = paste(summary.left.prop, summary.right.prop, sep = "\n"))
  
  result_list = list(p.stats = p.stats,
                     stats.df = stats.df,
                     stats.summary = stats.summary)
  return(result_list)
}

######################################################################################

doReplicateDeletionAnalysis = function(name, replicate, type, sam_reads, sites, hdr_profile, wt_profile, ref_sequence, read_data_all = NULL)
{
  stats = list(name = name, replicate = replicate, type = type)
  
  region_length = (sites$end - sites$start + 1)
  if (nchar(wt_profile) != region_length) {
    stop(sprintf("The WT allele profile given has length %d, but the region size (end - start + 1) is %d", nchar(wt_profile), region_length))
  }
  if (is.na(hdr_profile)) {
    hdr_profile = ""
  }
  
  replicate_name = paste(name, replicate, type, sep="_")
  if (is.null(read_data_all)) {
    read_data = getAlignedReads(sam_reads, ref_sequence, sites$start, sites$end)
  } else {
    read_data = read_data_all[[replicate_name]]
    if (is.null(read_data)) {
      cat(sprintf("\nERROR: Failed to get read data for replicate %s from file %s\n", replicate, opt$read_data))
      stop()
    }
  }
  reads.df = read_data$alignedReads
  
  cat(sprintf("%d starting cDNA reads\n", sum(reads.df$count)))
  cat(sprintf("%d cDNA reads of %d total (%.1f%%) were soft-clipped\n", read_data$num_softclipped, read_data$num_reads, 100.0 * read_data$num_softclipped / read_data$num_reads))
  cat(sprintf("%d cDNA reads of %d total (%.1f%%) were hard-clipped\n", read_data$num_hardclipped, read_data$num_reads, 100.0 * read_data$num_hardclipped / read_data$num_reads))
  cat(sprintf("%d cDNA reads of %d total (%.1f%%) had insertions\n", read_data$num_insertion, read_data$num_reads, 100.0 * read_data$num_insertion / read_data$num_reads))
  cat(sprintf("%d cDNA reads of %d total (%.1f%%) were likely primer dimers\n", read_data$num_primerdimer, read_data$num_reads, 100.0 * read_data$num_primerdimer / read_data$num_reads))
  stats$num_softclipped = read_data$num_softclipped
  stats$num_hardclipped = read_data$num_hardclipped
  stats$num_insertion = read_data$num_insertion
  stats$num_primerdimer = read_data$num_primerdimer
  stats$num_reads = read_data$num_reads
  
  hdr_profile_chars = strsplit(hdr_profile, "")[[1]]
  wt_profile_chars = strsplit(wt_profile, "")[[1]]
  ref_seq_chars = strsplit(ref_sequence, "")[[1]]
  
  reads.df$read_chars = sapply(reads.df$region_read, FUN=function(s) strsplit(s, ""))
  reads.df$seq_length = sapply(reads.df$read_chars, FUN=function(s) sum(isDNALetter(s)))
  
  #reads.df$seq_length = sapply(reads.df$read_chars, FUN=function(s) sum(s %in% c("A", "C", "G", "T")))  # SLOWER
  #reads.df$seq_length = sapply(reads.df$region_read, FUN=function(s) str_count(s, "[ACGT]"))  # SLOWER
  #reads.df$seq_length = sapply(reads.df$region_read, FUN=function(s) nchar(gsub("[-*]", "", s)))  # SLOWER
  
  reads.df$mismatch_count = sapply(reads.df$read_chars, FUN=getMismatchCharsCount, ref_seq_chars)
  
  #reads.df$mismatch_count = sapply(reads.df$read_chars, getMismatchCharsCount, ref_seq_chars)  # SLOWER
  #reads.df$mismatch_count = sapply(reads.df$region_read, getMismatchCount, ref_sequence)  # SLOWER
  
  reads.df$spanning_read = T
  span_site = sites$cut_site - sites$start + 1
  if (!is.na(sites$highlight_site)) {
    span_site = sites$highlight_site - sites$start + 1
  }
  reads.df$spanning_read = sapply(reads.df$region_read, FUN=function(s) (substring(s, span_site, span_site) != '-'))
  #reads.df$spanning_read = sapply(reads.df$read_chars, FUN=function(s) (s[span_site] != '-'))
  if (opt$exclude_nonspanning_reads) {
    stats$reads_excluded_nonspanning = sum((!reads.df$spanning_read) * reads.df$count)
    cat(sprintf("%d of %d reads (%.2f%%) excluded due to not spanning the site of interest (position %d)\n",
                stats$reads_excluded_nonspanning, stats$num_reads, 100.0 * stats$reads_excluded_nonspanning / stats$num_reads,
                span_site))
    reads.df = reads.df[reads.df$spanning_read, ]
  }

  # Exclude reads that don't cover enough of the region of interest
  exclude_for_overlap = (reads.df$seq_length < opt$min_window_overlap)
  stats$reads_excluded_for_minoverlap = sum(exclude_for_overlap * reads.df$count)
  cat(sprintf("%d of %d reads (%.2f%%) excluded due to aligning to less than %d bp in the region of interest\n",
              stats$reads_excluded_for_minoverlap, stats$num_reads, 100.0 * stats$reads_excluded_for_minoverlap / stats$num_reads,
              opt$min_window_overlap))
  reads.df = reads.df[!exclude_for_overlap, ]

  # Identify unique deletion profile (UDP) for each read
  #reads.df$udp = getReadUDPs(reads.df$region_read, wt_profile_chars)
  reads.df$udp = sapply(reads.df$read_chars, FUN=getReadCharsUDP, wt_profile_chars)
  #reads.df$udp = sapply(reads.df$region_read, FUN=getReadUDP, wt_profile_chars) # SLOWER
  
  reads.df$is_wt_allele = (reads.df$udp == wt_profile)
  # make sure we only count as WT those reads which actually cover the HDR site
  reads.df$is_wt_allele[!reads.df$spanning_read] = NA
  
  reads.df$has_any_deletion = sapply(reads.df$udp, FUN=function(s) grepl("[*]", s))
  rel_cut_site = sites$cut_site - sites$start + 1
  if (opt$exclude_nonspanning_deletions) {
    reads.df$has_crispr_deletion = (!reads.df$is_wt_allele & grepl("[*]", substr(reads.df$udp, rel_cut_site - opt$editing_window, rel_cut_site + opt$editing_window)) )
    
    # We change the UDPs to "zero out" any deletions that don't overlap
    # with the "editing region" defined by the cut site and editing_window
    updateUDPEdits = function(udp, cut_site, editing_window) {
      dels = str_locate_all(udp, "[*]+")[[1]]
      if (nrow(dels) == 0)
        return(udp)
      for (i in 1:nrow(dels)) {
        start = dels[i,1]
        end = dels[i,2]
        if (start > cut_site + editing_window | end < cut_site - editing_window) {
          substr(udp, start, end) <- strrep("-", end - start + 1)
        }
      }
      return(udp)
    }
    if (sum(reads.df$has_any_deletion) > 0) {
      reads.df$udp[reads.df$has_any_deletion & !reads.df$is_wt_allele] = sapply(reads.df$udp[reads.df$has_any_deletion & !reads.df$is_wt_allele], updateUDPEdits, rel_cut_site, opt$editing_window)
      # If we have zeroed out a deletion, then it's important that we also zero out
      # the flag for whether it has any deletion.
      reads.df$has_any_deletion[reads.df$has_any_deletion] = sapply(reads.df$udp[reads.df$has_any_deletion], FUN=function(s) grepl("[*]", s))
    }
  } else {
    # Accept deletions anywhere in the read
    reads.df$has_crispr_deletion = reads.df$has_any_deletion
  }
  # Include deletions anywhere in the read
  reads.df$has_multiple_deletions = reads.df$has_any_deletion
  if (sum(reads.df$has_multiple_deletions) > 0) {
    reads.df$has_multiple_deletions[reads.df$has_multiple_deletions] = 
      sapply(reads.df$udp[reads.df$has_multiple_deletions], FUN=function(s) grepl("[*]+[^*]+[*]+", s))
  }
  
  # Do this again, to count as WT any reads which had non-CRISPR deletions zeroed out
  reads.df$is_wt_allele = (reads.df$udp == wt_profile)
  # make sure we only count as WT those reads which actually cover the HDR site
  reads.df$is_wt_allele[!reads.df$spanning_read] = NA
  # If the WT profile has a deletion, then these should not be counted as CRISPR deletions
  reads.df$has_crispr_deletion[reads.df$is_wt_allele] = F
  
  n_wt_vars = sum(isDNALetter(wt_profile_chars) & wt_profile_chars != ref_seq_chars)
  if (n_wt_vars > 0) {
    reads.df$mismatch_count[reads.df$is_wt_allele] = reads.df$mismatch_count[reads.df$is_wt_allele] - n_wt_vars
  }
  
  # Identify which reads are HDR.
  reads.df$is_hdr_allele = F
  if (hdr_profile != "") {
    if (nchar(hdr_profile) != region_length) {
      stop(sprintf("The HDR allele profile given has length %d, but the region size (end - start + 1) is %d", nchar(hdr_profile), region_length))
    }
    # Check that the HDR profile and WT profile have DNA characters at
    # the same positions
    # 
    #if (any((hdr_profile_chars == "-") != (wt_profile_chars == "-"))) {
    #  stop("Error: HDR profile and WT profile should both indicate the expected sequence letters at the same positions")
    #}
    if (grepl("B|D|H|V", hdr_profile)) {
      hdr_profile_set = getHDRProfiles(hdr_profile)
      reads.df$is_hdr_allele = sapply(reads.df$udp, function(udp) any(udp == hdr_profile_set))
    } else {
      reads.df$is_hdr_allele = (reads.df$udp == hdr_profile)
    }
    
    if (!any(hdr_profile_chars == "*")) {
      # If the HDR profile doesn't have a deletion itself, then we don't expect any HDR alleles to have
      # deletions near the CRISPR cut site
      reads.df$is_hdr_allele = reads.df$is_hdr_allele & !reads.df$has_crispr_deletion
    }
    # We don't want to count the HDR site itself as a mismatch
    n_hdr_vars = sum(isDNALetter(hdr_profile_chars) & hdr_profile_chars != ref_seq_chars)
    if (n_hdr_vars > 0) {
      reads.df$mismatch_count[reads.df$is_hdr_allele] = reads.df$mismatch_count[reads.df$is_hdr_allele] - n_hdr_vars
    }
  }
  if (any(reads.df$mismatch_count < 0)) {
    stop("ERROR: something went wrong - got a negative mismatch count.")
  }
  
  # Exclude reads with too many mismatches
  exclude_for_mismatches = (reads.df$mismatch_count / reads.df$seq_length) > opt$max_mismatch_frac
  stats$reads_excluded_for_mismatches = sum(exclude_for_mismatches * reads.df$count)
  cat(sprintf("%d of %d reads (%.2f%%) excluded due to having more than %.2f%% of the read being mismatches\n",
              stats$reads_excluded_for_mismatches, stats$num_reads, 100.0 * stats$reads_excluded_for_mismatches / stats$num_reads,
              opt$max_mismatch_frac * 100))
  reads.df = reads.df[!exclude_for_mismatches, ]
  
  n_wt_reads = sum(reads.df$is_wt_allele * reads.df$count, na.rm = T)
  if (n_wt_reads < 1) {
    warning("ERROR: no wild-type reads found. Check the WT allele profile. You may also want to check that your reads span the edit site, given your read length and amplicon coords.")
  } else if (n_wt_reads < 100) {
    warning(sprintf("Warning: only %d wild-type reads found in experiment.", n_wt_reads))
  }
  
  stats$num_wt_reads = n_wt_reads
  stats$num_hdr_reads = sum(reads.df$is_hdr_allele * reads.df$count, na.rm = T)
  stats$num_deletion_reads = sum((reads.df$has_crispr_deletion & !reads.df$is_hdr_allele) * reads.df$count, na.rm = T)
  stats$num_kept_reads = sum(reads.df$count, na.rm = T)
  
  # Count deletion reads with different windows around the site of interest
  has_deletion_2bp_window = grepl("[*]", substr(reads.df$udp, span_site - 1, span_site + 1))
  has_deletion_10bp_window = grepl("[*]", substr(reads.df$udp, span_site - 5, span_site + 5))
  has_deletion_20bp_window = grepl("[*]", substr(reads.df$udp, span_site - 10, span_site + 10))
  stats$num_deletions_2bp_window = sum((has_deletion_2bp_window & !reads.df$is_hdr_allele) * reads.df$count, na.rm = T)
  stats$num_deletions_10bp_window = sum((has_deletion_10bp_window & !reads.df$is_hdr_allele) * reads.df$count, na.rm = T)
  stats$num_deletions_20bp_window = sum((has_deletion_20bp_window & !reads.df$is_hdr_allele) * reads.df$count, na.rm = T)
  
  reads.df$has_custom_deletion = F
  if (opt$custom_del_span) {
    reads.df$has_custom_deletion = ( grepl("[*]", substr(reads.df$udp, span_site, span_site)) &
                                       !grepl("[*]", substr(reads.df$udp, 1, opt$del_span_start)) &
                                       !grepl("[*]", substr(reads.df$udp, opt$del_span_end, str_length(reads.df$udp))) )
    stats$num_deletion_reads = sum((reads.df$has_custom_deletion & !reads.df$is_hdr_allele) * reads.df$count, na.rm = T)
  }
  
  # Aggregate reads according to their UDP
  udp.df = summarise(reads.df %>% group_by(udp),
                     num_reads = sum(count),
                     is_hdr_allele = first(is_hdr_allele),
                     is_wt_allele = first(is_wt_allele),
                     spanning_read = first(spanning_read),
                     has_any_deletion = first(has_any_deletion),
                     has_crispr_deletion = first(has_crispr_deletion),
                     has_custom_deletion = first(has_custom_deletion),
                     has_multiple_deletions = first(has_multiple_deletions),
                     avg_seq_length = mean(seq_length),
                     avg_mismatch_count = mean(mismatch_count)) %>%
    mutate(name = name, replicate = replicate, type = type) %>%
    select(name, replicate, type, everything()) %>%
    arrange(-num_reads)
  stats$num_udps = nrow(udp.df)
  
  udp.df$has_edit = udp.df$has_crispr_deletion | udp.df$is_hdr_allele
  stats$num_edit_reads = sum(udp.df$num_reads[udp.df$has_edit], na.rm = T)
  
  stats$reads_excluded_for_multiple_deletions = 0
  if (opt$exclude_multiple_deletions) {
    stats$reads_excluded_for_multiple_deletions = sum(udp.df$num_reads[udp.df$has_multiple_deletions])
    stats$udps_excluded_for_multiple_deletions = sum(udp.df$has_multiple_deletions)
    cat(sprintf("%d of %d UDPs (%d of %d reads, %.2f%%) excluded due to having more than one separate deletion\n",
                stats$udps_excluded_for_multiple_deletions, stats$num_udps,
                stats$reads_excluded_for_multiple_deletions, stats$num_reads, 100.0 * stats$reads_excluded_for_multiple_deletions / stats$num_reads))
    udp.df = udp.df[!udp.df$has_multiple_deletions, ]
  }
  
  # Order UDPs by the deletion start position and plot
  getDeletionLoc = function(i) {
    if (udp.df$has_crispr_deletion[i]) {
      return(sapply(udp.df$udp[i], FUN=function(udp) str_locate_all(udp, "[\\*]+")))
    }
    return(NA)
  }
  udp_dels = sapply(1:nrow(udp.df), FUN=getDeletionLoc)
  udp.df$deletion_start = sapply(udp_dels, FUN=function(x) {if (length(x) == 1) return(NA); x[1,1]})
  udp.df$deletion_end = sapply(udp_dels, FUN=function(x) {if (length(x) == 1) return(NA); x[1,2] + 1})
  udp.df$deletion2_start = sapply(udp_dels, FUN=function(x) {if (length(x) == 1) return(NA); if(nrow(x) > 1) {x[2,1]} else {NA}})
  udp.df$deletion2_end = sapply(udp_dels, FUN=function(x) {if (length(x) == 1) return(NA); if(nrow(x) > 1) {x[2,2]+1} else {NA}})
  #udp.df$deletion_length = sapply(udp_dels, FUN=function(x) x[2]-x[1]+1)
  
  # Make a table which has just the different versions of the WT and HDR alleles
  wt_hdr.df = NULL
  if (opt$allele_profile) {
    wt_hdr.reads.df = reads.df %>% dplyr::filter(is_wt_allele | is_hdr_allele)
    #wt_hdr.reads.df$mismatch_profile = getReadMismatchProfiles(wt_hdr.reads.df$region_read, ref_seq_chars)
    wt_hdr.reads.df$mismatch_profile = sapply(wt_hdr.reads.df$read_chars, FUN=getReadCharsMismatchProfile, ref_seq_chars)
    
    wt_hdr.df = summarise(wt_hdr.reads.df %>% group_by(mismatch_profile),
                          num_reads = sum(count),
                          mismatch_count = first(mismatch_count),
                          spanning_read = first(spanning_read),
                          is_wt_allele = first(is_wt_allele),
                          is_hdr_allele = first(is_hdr_allele)) %>%
      mutate(name = locus_name, replicate = replicate, type = type) %>%
      dplyr::select(name, replicate, type, everything()) %>%
      dplyr::arrange(-num_reads)
  }
  
  # Make a table which has all variations of read sequences at the "profile" sites,
  # i.e. those sites used to identify WT and HDR alleles
  sites.profile.df = NULL
  if (!opt$no_site_profile) {
    #reads.df$sites_profile = getReadSiteProfiles(reads.df$region_read, wt_profile_chars)
    profile_positions = which(isDNALetter(wt_profile_chars))
    reads.df$sites_profile = sapply(reads.df$region_read, FUN=getReadSiteProfile, profile_positions)
    sites.profile.df = reads.df %>% group_by(sites_profile) %>%
      summarise(count = sum(count))  %>%
      mutate(name = locus_name, replicate = replicate, type = type) %>%
      dplyr::select(name, replicate, type, everything())
  }
  return(list(udp.df = udp.df, wt_hdr.df = wt_hdr.df, sites.profile.df = sites.profile.df, stats = stats, read_data = read_data))
}

getRelativeSites = function(sites, viewing_window = NULL) {
  rel_sites = list()
  region_length = (sites$end - sites$start + 1)
  rel_sites$cut_site = sites$cut_site - sites$start + 1
  rel_sites$highlight_site = sites$highlight_site - sites$start + 1
  rel_sites$start = 1
  rel_sites$end = region_length
  if (!is.null(viewing_window)) {
    # Subset everything to a window around the cut site - including HDR profile, WT profile, ref_sequence, region read
    rel_sites$start = rel_sites$cut_site - viewing_window
    if (rel_sites$start < 1) {
      rel_sites$start = 1
    }
    rel_sites$end = rel_sites$cut_site + viewing_window - 1
    if (rel_sites$end > region_length) {
      rel_sites$end = region_length
    }
    rel_sites$window_cut_site = rel_sites$cut_site - rel_sites$start + 1
    rel_sites$window_highlight_site = rel_sites$highlight_site - rel_sites$start + 1
  }
  return(rel_sites)
}

getRegionReadsFromBam = function(name, replicate, bam_file, chr, start, end, ref_sequence) {
  alignmentFlag = "-F 0x904 " # only include primary alignment for each read, and must be mapped
  mapqFlag = sprintf("-q %d ", opt$minMapQ)
  subsampleStr = ""
  if (!is.null(opt$subsample)) {
    if (opt$subsample < 1) {
      subsampleStr = sprintf("-s %f ", opt$subsample)
    }
  }
  regionStr = sprintf("%s:%d-%d", chr, start, end)
  cmd = paste0("samtools view ", alignmentFlag, mapqFlag, subsampleStr, bam_file, " ", regionStr)
  #cmd = "samtools view " + alignmentFlag + mapqFlag + subsampleStr + bam_file + " " + regionStr + " | head -n 50"
  cat(paste0("Calling:\n", cmd, "\n"))
  sam_reads = system(cmd, intern=TRUE)
  return(sam_reads)
}


# this function takes a .sam read and registers to the reference genome, relative
# to the region of interest in the reference
# '-' not sequenced nucleotides
# '*' deleted nucleotides
# if the cigars contains H or I the variable 'read_ok' will be False
# the returned read should be the same length as the reference genome
registerRead = function(relative_pos, sam_cigar, sam_read, region_length) {
  sam_cigar_parsed = parseCigar(sam_cigar)
  read_ok = TRUE
  read_status = ""
  expanded_cigar = ""
  cigar_bits = ""
  cursor_read = 1
  match_length = 0
  del_length = 0
  read_span_to = relative_pos + 1
  left_softclip_len = 0
  right_softclip_pos = -1
  right_softclip_len = 0
  
  for (i in 1:length(sam_cigar_parsed$letters)) {
    type = sam_cigar_parsed$letters[i]
    
    if (type == 'M') {
      match_length = match_length + sam_cigar_parsed$numbers[i]
      cigar_bits[i] = substr(sam_read, cursor_read, cursor_read + sam_cigar_parsed$numbers[i] - 1)
      #expanded_cigar = paste0(expanded_cigar, substr(sam_read, cursor_read, cursor_read + sam_cigar_parsed$numbers[i] - 1))
      cursor_read = cursor_read + sam_cigar_parsed$numbers[i]
      read_span_to = read_span_to + sam_cigar_parsed$numbers[i]
    }
    else if (type == 'D') {
      cigar_bits[i] = strrep('*', sam_cigar_parsed$numbers[i])
      #expanded_cigar = paste0(expanded_cigar, strrep('*', sam_cigar_parsed$numbers[i]))
      del_length = del_length + sam_cigar_parsed$numbers[i]
      read_span_to = read_span_to + sam_cigar_parsed$numbers[i]
    } else if (type == 'N') {
      read_status = "spliced"
      cigar_bits[i] = strrep('-', sam_cigar_parsed$numbers[i])
    } else if (type == 'S') {
      if (left_softclip_len > 0) {
        right_softclip_pos = read_span_to
        right_softclip_len = sam_cigar_parsed$numbers[i]
      } else {
        left_softclip_len = sam_cigar_parsed$numbers[i]
      }
      cursor_read = cursor_read + sam_cigar_parsed$numbers[i]
      read_status = "softclipped"
    } else if (type == 'H') {
      read_status = "hardclipped"
      read_ok = F
      break
    } else if (type == 'I') {
      read_status = "insertion"
      read_ok = F
      break
    } else {
      read_ok = F
      break
    }
  }
  expanded_cigar = paste0(cigar_bits, collapse="")
  registered_read = ""
  if (read_ok) {
    if (relative_pos < 0) {
      if ((nchar(expanded_cigar) - relative_pos) > 0) {
        registered_read = substr(expanded_cigar, (-relative_pos)+1, nchar(expanded_cigar))
      }
    } else if (relative_pos == 0) {
      registered_read = expanded_cigar
    } else {
      registered_read = paste0(strrep('-', relative_pos), expanded_cigar)
    }
    
    len = nchar(registered_read)
    if (len > region_length) {
      registered_read = substr(registered_read, 1, region_length)
    } else if (len < region_length){
      registered_read = paste0(registered_read, strrep('-', region_length - len))
    }
  }
  
  registered_left_softclip = vector(mode = "integer", length = region_length)
  registered_right_softclip = vector(mode = "integer", length = region_length)
  if (left_softclip_len > 0 & relative_pos > 0) {
    # relative_pos is the read start position relative to the ref sequence, but soft-clipping
    # occurs to the left of the read start. E.g. relative_pos == 0 means the non-clipped read
    # sequence starts at the first base of the ref sequence.
    startpos = relative_pos + 1 - left_softclip_len
    if (startpos < 1) {
      left_softclip_len = left_softclip_len - (1 - startpos)
      startpos = 1
    }
    if (startpos + left_softclip_len - 1 > region_length) {
      left_softclip_len = region_length - startpos + 1
    }
    registered_left_softclip[startpos:(startpos + left_softclip_len - 1)] = 1
  }
  if (right_softclip_len > 0 & relative_pos < region_length) {
    startpos = right_softclip_pos
    if (startpos < 1) {
      right_softclip_len = right_softclip_len - (1 - startpos)
      startpos = 1
    }
    if (startpos + right_softclip_len - 1 > region_length) {
      right_softclip_len = region_length - startpos + 1
    }
    registered_right_softclip[startpos:(startpos + right_softclip_len - 1)] = 1
  }
  
  #print registered_read
  return(list("read" = registered_read, "left_softclip" = registered_left_softclip, "right_softclip" = registered_right_softclip,
              "read_ok" = read_ok, "read_status" = read_status, "match_length" = match_length, "deletion_length" = del_length))
}

# this function takes a string cigar "10M20D7M" and converts it to arrays of letters ['M','D','M'] and numbers [10,20,7]
parseCigar = function(input_cigar) {
  #print(input_cigar)
  X = strsplit(input_cigar, split="")[[1]] # splits the string into array of char
  numbers = numeric()
  letters = character()
  idx_last_letter = 0
  idx_current = 1
  
  for (x in X) {
    if (is.na(strtoi(x))) {
      letters = c(letters, x)
      idx_1 = idx_last_letter + 1
      idx_2 = idx_current - 1
      #numbers = c(numbers, as.integer(paste0(X[idx_1:idx_2], collapse="")))
      numbers = c(numbers, strtoi(substr(input_cigar, idx_1, idx_2)))
      idx_last_letter = idx_current
    }
    idx_current = idx_current + 1
  }
  return(list("letters" = letters, "numbers" = numbers))
}


# Goes through each SAM file read and gets its sequence with respect to the reference
# zone of interest, where "-" indicates that the read does not cover the position, and
# * indicates a deletion.
getAlignedReads = function(sam_reads, ref_sequence, start, end) {
  # reads_input: array of lines from a SAM file
  # ref_sequence: sequence of the amplicon
  # start, end: coordinates of the region of interest in the reference
  ref_seq_length = nchar(ref_sequence) # size of the ref genome
  if (ref_seq_length != (end - start + 1)) {
    stop(sprintf("Error: reference sequence length (%d) should be the same as the region size (%d = end - start + 1)",
                 ref_seq_length, end - start + 1))
  }
  
  sam_positions = vector(mode = "integer", length = length(sam_reads)) # Preallocate a vector of aligned reads
  sam_cigars = vector(mode = "character", length = length(sam_reads)) # Preallocate a vector of aligned reads
  discard = vector(mode = "logical", length = length(sam_reads)) # Preallocate a vector of discarded reads
  seq = vector(mode = "character", length = length(sam_reads)) # Preallocate a vector of discarded reads
  
  num_outsidewindow = 0
  for (i in 1:length(sam_reads)) {
    line_split = strsplit(sam_reads[i], "\t", fixed = T)[[1]]
    sam_positions[i] = as.integer(line_split[4])
    sam_cigars[i] = line_split[6]
    seq[i] = line_split[10]
    
    if (sam_cigars[i] == "*") {
      discard[i] = T
      next # Ignore this read, as it doesn't have a valid CIGAR
    }
    relative_pos = sam_positions[i] - start
    if (abs(relative_pos) > 1000) {
      num_outsidewindow = num_outsidewindow + 1
      stop("Likely ERROR: read position and reference coordinates differ by more than 1000 bp")
    }
  }
  reads.df = data.frame(sam_read = sam_reads,
                        seq = seq,
                        sam_position = sam_positions,
                        sam_cigar = sam_cigars,
                        discard = discard)
  reads.summary.df = reads.df %>% group_by(sam_position, sam_cigar, seq) %>%
    summarise(count = n(), discard = first(discard)) %>%
    dplyr::mutate(relative_pos = sam_position - start)
  
  discardedReads.df = reads.summary.df %>% dplyr::filter(discard == TRUE)
  keptReads.df = reads.summary.df %>% dplyr::filter(!discard)
  
  alignedReads = vector(mode = "character", length = nrow(keptReads.df)) # Preallocate a vector of aligned reads
  
  num_reads = length(sam_reads)
  num_softclipped = 0
  num_hardclipped = 0
  num_insertion = 0
  
  # As we go through the aligned reads we check for possible primer dimers
  num_primerdimer = 0
  primer_dimer_start = substr(ref_sequence, 1, 20)
  primer_dimer_end = substr(ref_sequence, nchar(ref_sequence) - 19, nchar(ref_sequence))
  isPrimerDimer = function(reg_read, match_length) {
    (match_length < 40 & (substr(reg_read, 1, 20) == primer_dimer_start | substr(reg_read, nchar(reg_read) - 19, nchar(reg_read)) == primer_dimer_end))
  }
  
  for (i in 1:nrow(keptReads.df)) {
    # computes the read sequence within the coords of interest
    output = registerRead(keptReads.df$relative_pos[i], keptReads.df$sam_cigar[i], keptReads.df$seq[i], ref_seq_length)
    if (output$read_ok) {
      alignedReads[i] = output$read
    }
    if (output$read_status == "hardclipped") {
      num_hardclipped = num_hardclipped + keptReads.df$count[i]
      if (isPrimerDimer(output$read, output$match_length)) {
        num_primerdimer = num_primerdimer + keptReads.df$count[i]
      }
    } else if (output$read_status == "insertion") {
      num_insertion = num_insertion + keptReads.df$count[i]
    } else if (output$read_status == "softclipped") {
      num_softclipped = num_softclipped + keptReads.df$count[i]
      if (isPrimerDimer(output$read, output$match_length)) {
        num_primerdimer = num_primerdimer + keptReads.df$count[i]
      }
    }
  }
  keptReads.df$region_read = alignedReads
  keptReads.df = keptReads.df %>% dplyr::filter(region_read != "") %>%
    dplyr::select(region_read, count, sam_position, sam_cigar, seq)
  # Some reads may have a distinct sequence yet have the same registered
  # read sequence, because the soft-clipped sequence could differ.
  
  return(list("alignedReads" = keptReads.df, "discardedReads" = discardedReads.df, "num_reads" = num_reads,
              "num_hardclipped" = num_hardclipped, "num_insertion" = num_insertion,
              "num_softclipped" = num_softclipped, "num_primerdimer" = num_primerdimer,
              "num_outsidewindow" = num_outsidewindow))
}


getUDPPlot = function(udp.df, plot_title, sites) {
  plot_list = list()
  # We allow there to be up to 2 deletions in the UDP
  udp.df = udp.df %>% filter(has_crispr_deletion, !is.na(deletion_start))
  if (nrow(udp.df) == 0) {
    return(egg::ggarrange(textPlot("No UDPs to plot."), top=plot_title, draw = F))
  }
  udp.dels.df = udp.df %>% dplyr::select(deletion_start, deletion_end, deletion2_start, deletion2_end, has_custom_deletion, sharing=udp_sharing)
  udp.dels.df = udp.dels.df %>% arrange(deletion_start, deletion2_start)
  udp.dels.df$y = 1:nrow(udp.dels.df)
  udp.plot.df = bind_rows(udp.dels.df %>% dplyr::select(y, deletion_start, deletion_end, has_custom_deletion, sharing),
                          udp.dels.df %>% dplyr::select(y, deletion_start=deletion2_start, deletion_end=deletion2_end, has_custom_deletion, sharing))
  udp.plot.df = udp.plot.df %>% dplyr::filter(!is.na(deletion_start))
  
  xmax = nchar(udp.df$udp[1])
  segment_size = 0.5
  if (nrow(udp.dels.df) > 200) {
    segment_size = 0.3
  }
  simple_theme = theme_bw() + theme(panel.grid.minor = element_line(colour="grey90", size=0.1), panel.grid.major = element_line(colour="grey90", size=0.1))
  if (any(udp.plot.df$has_custom_deletion)) {
    p.udp_dist = ggplot(udp.plot.df) +
      geom_segment(aes(x = deletion_start, xend = deletion_end, y = y, yend = y, color = has_custom_deletion), size = segment_size) +
      simple_theme + coord_cartesian(xlim=c(sites$start, sites$end)) + scale_y_reverse() +
      scale_color_manual(values=c("TRUE"="red", "FALSE"="dodgerblue3")) +
      theme(legend.position = c(0.25, 0.2), legend.spacing.y = unit(0.1, "cm"), legend.key.height = unit(0.45, "cm")) +
      xlab("Nucleotide position") + ylab("Cumulative UDP count")
  } else {
    p.udp_dist = ggplot(udp.plot.df) +
      geom_segment(aes(x = deletion_start, xend = deletion_end, y = y, yend = y, color = sharing), size = segment_size) +
      simple_theme + coord_cartesian(xlim=c(sites$start, sites$end)) + scale_y_reverse() +
      scale_color_manual(values=c("cDNA only"="red", "gDNA only"="orange", "both"="dodgerblue3", "unclear"="grey")) +
      theme(legend.position = c(0.25, 0.2), legend.spacing.y = unit(0.1, "cm"), legend.key.height = unit(0.45, "cm")) +
      xlab("Nucleotide position") + ylab("Cumulative UDP count")
  }
  
  count.plot.df = data.frame(x=1:xmax)
  udp.char.matrix = str_split_fixed(udp.df$udp, "", n = nchar(udp.df$udp[1]))
  # isPositionDel = function(i) {
  #   sapply(udp.df$udp, FUN=function(x) substr(x, i, i) == '*')
  # }
  isPositionDel = function(i) {
    (udp.char.matrix[,i] == '*')
  }
  getDelCount = function(i) {
    sum(udp.char.matrix[,i] == '*')
  }
  getDelReadCount = function(i, udp.df) {
    sum(isPositionDel(i) * udp.df$num_reads)
  }
  
  count.plot.df$udp_count = sapply(count.plot.df$x, FUN=getDelCount)
  count.plot.df$read_count = sapply(count.plot.df$x, FUN=getDelReadCount, udp.df)
  
  count.plot.df$udp_count_custom_del = 0
  count.plot.df$read_count_custom_del = 0
  if (any(udp.df$has_custom_deletion)) {
    udp.custom_del.df = udp.df %>% dplyr::filter(has_custom_deletion)
    udp.char.matrix = str_split_fixed(udp.custom_del.df$udp, "", n = nchar(udp.custom_del.df$udp[1]))
    count.plot.df$udp_count_custom_del = sapply(count.plot.df$x, FUN=getDelCount)
    count.plot.df$read_count_custom_del = sapply(count.plot.df$x, FUN=getDelReadCount, udp.custom_del.df)
  }
  
  p.udp_count = ggplot(count.plot.df, aes(x=x, y=udp_count)) +
    geom_bar(stat="identity", fill="dodgerblue2") +
    geom_bar(mapping = aes(x=x, y=udp_count_custom_del), stat="identity", fill="red") +
    simple_theme + theme(axis.title.x=element_blank()) + ylab("UDP count") +
    coord_cartesian(xlim=c(sites$start, sites$end))
  p.read_count = ggplot(count.plot.df, aes(x=x, y=read_count)) +
    geom_bar(stat="identity", fill="dodgerblue2") +
    geom_bar(mapping = aes(x=x, y=read_count_custom_del), stat="identity", fill="red") +
    simple_theme + theme(axis.title.x=element_blank()) + ylab("Read count") +
    coord_cartesian(xlim=c(sites$start, sites$end))
  
  if (!is.na(sites$highlight_site)) {
    p.udp_dist = p.udp_dist + geom_vline(xintercept = sites$highlight_site, color="darkgreen", alpha=0.5)
    p.udp_count = p.udp_count + geom_vline(xintercept = sites$highlight_site, color="darkgreen", alpha=0.5)
    p.read_count = p.read_count + geom_vline(xintercept = sites$highlight_site, color="darkgreen", alpha=0.5)
  }
  if (!is.na(sites$cut_site)) {
    p.udp_dist = p.udp_dist + geom_vline(xintercept = sites$cut_site, color="grey20", linetype = "longdash", alpha=0.5)
    p.udp_count = p.udp_count + geom_vline(xintercept = sites$cut_site, color="grey20", linetype = "longdash", alpha=0.5)
    p.read_count = p.read_count + geom_vline(xintercept = sites$cut_site, color="grey20", linetype = "longdash", alpha=0.5)
  }
  if (opt$custom_del_span) {
    p.udp_dist = p.udp_dist + geom_vline(xintercept = opt$del_span_start, color="red", size=0.5, alpha=0.3) + geom_vline(xintercept = opt$del_span_end, color="red", size=0.5, alpha=0.3)
    p.udp_count = p.udp_count + geom_vline(xintercept = opt$del_span_start, color="red", size=0.5, alpha=0.3) + geom_vline(xintercept = opt$del_span_end, color="red", size=0.5, alpha=0.3)
    p.read_count = p.read_count + geom_vline(xintercept = opt$del_span_start, color="red", size=0.5, alpha=0.3) + geom_vline(xintercept = opt$del_span_end, color="red", size=0.5, alpha=0.3)
  }
  p.full = egg::ggarrange(p.read_count, p.udp_count, p.udp_dist, ncol=1, heights=c(1,1,4),
                          top = plot_title, draw = F)
  return(p.full)
}

getDeletionProfilePlot = function(replicate.udp.df, sites, plot_title = NA, show_average = T, show_replicates = T, ratioToWT = F) {
  if (!show_average & !show_replicates) {
    stop("getDeletionProfilePlot: One of show_average and show_replicates should be true.")
  }
  if (nrow(replicate.udp.df) == 0) {
    return(egg::ggarrange(textPlot("No UDPs to plot."), top=plot_title, draw = F))
  }
  xmax = nchar(replicate.udp.df$udp[1])
  #udp.char.matrix = NULL
  #cur.df = NULL
  # Get the deletion profile for each replicate separately
  isPositionDel = function(udp.char.matrix, i) {
    (udp.char.matrix[,i] == '*')
  }
  getDelReadCount = function(cur.df, udp.char.matrix, i, filterUdps) {
    if (filterUdps) {
      # only count deletions in UDPs where count_udp == 1
      sum(isPositionDel(udp.char.matrix, i) * cur.df$num_reads * cur.df$count_udp)
    } else {
      sum(isPositionDel(udp.char.matrix, i) * cur.df$num_reads)
    }
  }
  getDelpctDataframe = function(unique_reps.df, filterUdps = T, ratioToWT = F) {
    counts.list = list()
    for (i in 1:nrow(unique_reps.df)) {
      curType = unique_reps.df[i,]$type
      curReplicate = unique_reps.df[i,]$replicate
      count.df = data.frame(x=1:xmax, type = curType, replicate = curReplicate)
      cur.df = replicate.udp.df %>% dplyr::filter(replicate == curReplicate, type == curType)
      udp.char.matrix = str_split_fixed(cur.df$udp, "", n = nchar(cur.df$udp[1]))
      if (ratioToWT) {
        wtCount = sum(cur.df$is_wt_allele * cur.df$num_reads)
        count.df$del_pct = sapply(count.df$x, FUN=function(i) {getDelReadCount(cur.df, udp.char.matrix, i, filterUdps)}) / wtCount
      } else {
        count.df$del_pct = 100 * sapply(count.df$x, FUN=function(i) {getDelReadCount(cur.df, udp.char.matrix, i, filterUdps)}) / sum(cur.df$num_reads)
      }
      counts.list = c(counts.list, list(count.df))
    }
    # Combine deletion profiles for replicates
    bind_rows(counts.list)
  }
  
  unique_reps.df = unique(replicate.udp.df %>% dplyr::select(type, replicate))
  delpct.plot.df = getDelpctDataframe(unique_reps.df, filterUdps=T, ratioToWT=ratioToWT)
  delpct.plot.df$replicate = as.character(delpct.plot.df$replicate)
  delpct.plot.df$type = paste(delpct.plot.df$type, "replicate")
  
  if (show_average) {
    # Merge gDNA replicates together (and similarly for cDNA), and
    # determine deletion profiles
    # merged.udp.df = summarise(replicate.udp.df %>% group_by(udp, type),
    #                           num_reads = sum(num_reads)) %>%
    #   arrange(-num_reads)
    merged.udp.df = summarise(replicate.udp.df %>% group_by(type, udp),
                              num_reads = sum(num_reads), count_udp = first(count_udp), is_wt_allele = first(is_wt_allele)) %>%
      arrange(-num_reads)
    unique_reps.df = unique(merged.udp.df %>% ungroup() %>% dplyr::select(type))
    counts.list = list()
    for (i in 1:nrow(unique_reps.df)) {
      curType = unique_reps.df[i,]$type
      count.df = data.frame(x=1:xmax, type = curType)
      cur.df = merged.udp.df %>% dplyr::filter(type == curType)
      udp.char.matrix = str_split_fixed(cur.df$udp, "", n = nchar(cur.df$udp[1]))
      if (ratioToWT) {
        wtCount = sum(cur.df$is_wt_allele * cur.df$num_reads)
        count.df$del_pct = sapply(count.df$x, FUN=function(i) {getDelReadCount(cur.df, udp.char.matrix, i, filterUdps=T)}) / wtCount
      } else {
        count.df$del_pct = 100 * sapply(count.df$x, FUN=function(i) {getDelReadCount(cur.df, udp.char.matrix, i, filterUdps=T)}) / sum(cur.df$num_reads)
      }
      counts.list = c(counts.list, list(count.df))
    }
    merged.delpct.plot.df = bind_rows(counts.list)
    merged.delpct.plot.df$type[merged.delpct.plot.df$type == "gDNA"] = "gDNA average"
    merged.delpct.plot.df$type[merged.delpct.plot.df$type == "cDNA"] = "cDNA average"
    merged.delpct.plot.df$replicate = merged.delpct.plot.df$type
  }
  
  if (any(!replicate.udp.df$count_udp)) {
    # Make a deletion profile for unfiltered UDPs
    unique_reps.df = unique(replicate.udp.df %>% dplyr::select(type, replicate))
    delpct.plot.unfiltered.df = getDelpctDataframe(unique_reps.df, filterUdps=F)
    delpct.plot.unfiltered.df$replicate = as.character(delpct.plot.unfiltered.df$replicate)
    delpct.plot.unfiltered.df$type = paste(delpct.plot.unfiltered.df$type, "rep unflt")
    delpct.plot.df = bind_rows(delpct.plot.df, delpct.plot.unfiltered.df)
  }
  
  alpha_values = c(`cDNA average`=0.4, `gDNA average`=0.4, `cDNA replicate`=0.9, `gDNA replicate`=0.9, `cDNA rep unflt`=0.9, `gDNA rep unflt`=0.9)
  color_values = c(`cDNA average`="darkorange", `gDNA average`="cyan3", `cDNA replicate`="red", `gDNA replicate`="blue", `cDNA rep unflt`="orange", `gDNA rep unflt`="turquoise4")
  size_values = c(`cDNA average`=1.8, `gDNA average`=1.8, `cDNA replicate`=0.3, `gDNA replicate`=0.3, `cDNA rep unflt`=0.3, `gDNA rep unflt`=0.3)
  if (show_average & show_replicates) {
    plot.df = bind_rows(delpct.plot.df, merged.delpct.plot.df)
  } else if (show_average) {
    plot.df = merged.delpct.plot.df
    size_values["cDNA average"] = size_values["gDNA average"] = 1.8
    alpha_values["cDNA average"] = alpha_values["gDNA average"] = 0.8
  } else {
    plot.df = delpct.plot.df
    size_values["cDNA replicate"] = size_values["gDNA replicate"] = 0.4
  }
  plot.df$type = factor(plot.df$type, levels = c("cDNA average", "gDNA average", "cDNA replicate", "gDNA replicate", "cDNA rep unflt", "gDNA rep unflt"))
  plot.df$lty_dash = plot.df$type %in% c("cDNA rep unflt", "gDNA rep unflt")
  simple_theme = theme_bw() + theme(panel.grid.minor = element_line(colour="grey90", size=0.1), panel.grid.major = element_line(colour="grey90", size=0.1))
  p.gDNA_cDNA = ggplot() +
    geom_line(aes(x=x, y=del_pct, color=type, size=type, alpha=type, linetype=lty_dash, group=paste(type, replicate)), data = plot.df) +
    simple_theme + xlab("Nucleotide position") + ylab("Deletion frequency (%)") +
    coord_cartesian(xlim=c(sites$start, sites$end)) +
    scale_color_manual(values = color_values) +
    scale_size_manual(values = size_values) +
    scale_alpha_manual(values = alpha_values) +
    scale_linetype_manual(values = c("TRUE"="dashed", "FALSE"="solid"), guide=F)
    
  if (!is.na(plot_title)) {
    p.gDNA_cDNA = p.gDNA_cDNA + ggtitle(plot_title)
    
    num_udps = length(unique(replicate.udp.df$udp))
    included_udps = length(unique(replicate.udp.df %>% filter(count_udp != 0) %>% .$udp))
    if (included_udps < num_udps) {
      p.gDNA_cDNA = p.gDNA_cDNA + labs(subtitle=sprintf("%d UDPs included out of %d", included_udps, num_udps))
    }
  }
  if (!is.na(sites$highlight_site)) {
    p.gDNA_cDNA = p.gDNA_cDNA + geom_vline(xintercept = sites$highlight_site, color="darkgreen", alpha=0.5)
  }
  if (!is.na(sites$cut_site)) {
    p.gDNA_cDNA = p.gDNA_cDNA + geom_vline(xintercept = sites$cut_site, color="grey20", linetype = "longdash", alpha=0.5)
  }
  if (ratioToWT) {
    p.gDNA_cDNA = p.gDNA_cDNA + ylab("Del:WT ratio")
  }
  return(p.gDNA_cDNA)
}


getReplicateStatsPlot = function(stats.df, region_title) {
  color_values = c(`cDNA`="firebrick1", `cDNA outlier`="orange1", `gDNA`="dodgerblue3", `gDNA outlier`="turquoise1")
  stats.df$replicate = fct_reorder(stats.df$replicate, as.integer(factor(stats.df$type)))
  stats.df = stats.df %>% group_by(type) %>%
    mutate(type2 = ifelse(!is.na(is_outlier) & is_outlier, paste(type, "outlier"), type)) %>%
    arrange(replicate)
  p1 = ggplot(stats.df, aes(x=replicate, y=num_reads, fill=type2)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%d", num_reads)), size = 2.4, position = position_stack(vjust = 0.5)) +
    theme_bw(10) + theme(axis.text.x=element_blank(), legend.title=element_blank(), axis.title.x = element_blank()) +
    scale_fill_manual(values = color_values) +
    ylab("Number of reads") +
    ggtitle("Number of reads")
  
  p2 = ggplot(stats.df, aes(x=replicate, y=num_udps, fill=type2)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%d", num_udps)), size = 2.4, position = position_stack(vjust = 0.5)) +
    theme_bw(10) + theme(axis.text.x=element_blank(), legend.title=element_blank(), axis.title.x = element_blank()) +
    scale_fill_manual(values = color_values, guide=F) +
    ylab("Number of UDPs") +
    ggtitle("Number of UDPs")
  
  p3 = ggplot(stats.df, aes(x=replicate, y=HDR_WT_ratio, fill=type2)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.3g", HDR_WT_ratio)), size = 2.4, position = position_stack(vjust = 0.5)) +
    theme_bw(10) + theme(axis.text.x=element_blank(), legend.title=element_blank(), axis.title.x = element_blank()) +
    scale_fill_manual(values = color_values, guide=F) +
    ylab("HDR:WT ratio") +
    ggtitle("HDR:WT ratio")
  
  p4 = ggplot(stats.df, aes(x=replicate, y=DEL_WT_ratio, fill=type2)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.3g", DEL_WT_ratio)), size = 2.4, position = position_stack(vjust = 0.5)) +
    theme_bw(10) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.title=element_blank()) +
    scale_fill_manual(values = color_values, guide=F) +
    ylab("Del:WT ratio") +
    ggtitle("Deletion:WT ratio")
  
  p.title = ggdraw() + draw_label(sprintf("%s replicate QC", region_title), fontface='bold')
  #p.replicate_qc = plot_grid(p.udp_fractions, p.udp_avg_deviation, ncol=1)
  #p.res = plot_grid(p.title, p.replicate_qc, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
  p.res = egg::ggarrange(p.title, p1, p2, p3, p4, ncol=1, heights=c(1,3,3,3,3), draw = F)
  p.res
}

getReplicateQCMetrics = function(stats.df, replicate.udp.df, max_udps = 20, min_avg_udp_fraction = 0.005, exclude_wt = FALSE) {
  # Here we get various metrics for each replicate that enable visual or automatic
  # QC. One of the key inputs is, for each replicate, the UDP fraction for the
  # top N UDPs. When plotted, this should enable visually identifying outlier
  # samples based on the variability of the UDP fractions.
  
  # We need to have a value for each fo the top N UDPs. To do this,
  # spread the replicates out into columns (cDNA_1, cDNA_2, etc.), with
  # fill = 0 for when a replicate is missing, and then gather back.
  getType = function(type_rep) { sapply(type_rep, FUN=function(x) strsplit(x, "^" , fixed=T)[[1]][1]) }
  getReplicate = function(type_rep) { sapply(type_rep, FUN=function(x) strsplit(x, "^" , fixed=T)[[1]][2]) }
  replicate.udp.filled.df = replicate.udp.df %>%
    dplyr::mutate(type_replicate = paste(type, replicate, sep="^")) %>%
    dplyr::select(udp, type_replicate, is_wt_allele, num_reads) %>%
    tidyr::spread(type_replicate, num_reads, fill = 0) %>%
    tidyr::gather(key="type_replicate", value="num_reads", -udp, -is_wt_allele) %>%
    dplyr::mutate(type = getType(type_replicate),
                  replicate = getReplicate(type_replicate)) %>%
    dplyr::select(udp, type, is_wt_allele, replicate, num_reads)
  
  replicate.totalreads.df = replicate.udp.df %>%
    dplyr::group_by(replicate) %>%
    dplyr::summarise(replicate_total_reads = sum(num_reads))
  udp.totalreads.df = replicate.udp.df %>%
    dplyr::group_by(udp) %>%
    dplyr::summarise(udp_total_reads = sum(num_reads),
                     udp_reads_gDNA = sum(num_reads[type == "gDNA"]),
                     udp_reads_cDNA = sum(num_reads[type == "cDNA"]),
                     is_wt_allele = first(is_wt_allele)) %>%
    dplyr::arrange(-udp_total_reads)
  
  if (exclude_wt) {
    replicate.udp.filled.df = replicate.udp.filled.df %>% dplyr::filter(!is_wt_allele)
    udp.totalreads.df = udp.totalreads.df %>% dplyr::filter(!is_wt_allele)
  }
  
  if (max_udps < nrow(udp.totalreads.df)) {
    udp.totalreads.df = udp.totalreads.df[1:max_udps,]
  }
  
  replicate.udp.filled.df = replicate.udp.filled.df %>%
    dplyr::inner_join(udp.totalreads.df, by = "udp") %>%
    dplyr::left_join(replicate.totalreads.df, by="replicate") %>%
    dplyr::mutate(udp_fraction = num_reads / replicate_total_reads) %>%
    dplyr::arrange(-udp_total_reads)
  udp.avg.fractions = replicate.udp.filled.df %>%
    dplyr::group_by(udp) %>%
    dplyr::summarise(avg_udp_fraction = mean(udp_fraction),
                     avg_gdna_udp_fraction = mean(udp_fraction[type == "gDNA"]))
  
  replicate.udp.filled.df$replicate = fct_reorder(replicate.udp.filled.df$replicate, as.integer(factor(replicate.udp.filled.df$type)))
  replicate.udp.filled.df = replicate.udp.filled.df %>%
    dplyr::left_join(udp.avg.fractions, by = "udp") %>%
    dplyr::filter(avg_udp_fraction >= min_avg_udp_fraction)
  
  num_udps = length(unique(replicate.udp.filled.df$udp))
  replicate.udp.filled.df$udp = factor(replicate.udp.filled.df$udp, levels=unique(replicate.udp.filled.df$udp))
  replicate.udp.filled.df$udp_id = factor(as.integer(replicate.udp.filled.df$udp))
  
  # Compute a statistic which, for each replicate, is the mean deviation of the
  # replicate's UDP fractions from the mean UDP fractions across replicates.
  replicate.udp_avg_deviation.df = replicate.udp.filled.df %>%
    dplyr::group_by(replicate) %>%
    dplyr::summarise(avg_udp_deviation = mean(abs(udp_fraction - avg_udp_fraction)),
                     type = first(type))
  
  # Compute the deviation relative to the mean of gDNA only
  replicate.udp_avg_deviation_from_gDNA.df = replicate.udp.filled.df %>%
    dplyr::group_by(replicate) %>%
    dplyr::summarise(avg_udp_deviation_from_gDNA = mean(abs(udp_fraction - avg_gdna_udp_fraction)),
                     type = first(type))
  
  # Compute the average RELATIVE deviation (CV) across UDPs for each replicate.
  replicate.udp_rel_deviation_from_gDNA.df = replicate.udp.filled.df %>%
    dplyr::group_by(replicate) %>%
    dplyr::summarise(avg_udp_rel_deviation_from_gDNA = mean(abs(udp_fraction / avg_gdna_udp_fraction - 1)),
                     type = first(type))
  
  stats.new.df = stats.df %>%
    dplyr::left_join(replicate.udp_avg_deviation.df, by=c("replicate", "type")) %>%
    dplyr::left_join(replicate.udp_avg_deviation_from_gDNA.df, by=c("replicate", "type")) %>%
    dplyr::left_join(replicate.udp_rel_deviation_from_gDNA.df, by=c("replicate", "type"))
  
  # Use multiple metrics across replicates to calculate outlier scores using KNN
  # and isolation forest, separately for cDNA and gDNA
  n_cDNA = sum(stats.df$type == "cDNA")
  n_gDNA = sum(stats.df$type == "gDNA")
  stats.new.df$outlier_score_knn = NA
  stats.new.df$outlier_score_iso = NA
  stats.new.df$is_outlier = NA
  if (n_cDNA >= 4 & !any(is.na(stats.new.df$avg_udp_deviation))) {
    stats.cDNA.df = stats.new.df %>%
      filter(type == "cDNA") %>%
      select(num_udps, HDR_WT_ratio, DEL_WT_ratio, avg_udp_deviation, avg_udp_deviation_from_gDNA)
    # Compute KNN outlier metric. First scale the data, and for UDP deviation set low values to
    # zero since these are good and shouldn't be considered outliers.
    stats_scaled <- scale(stats.cDNA.df)
    stats_scaled[is.na(stats_scaled)] = 0
    stats_scaled[, "avg_udp_deviation"] = sapply(stats_scaled[, "avg_udp_deviation"], function(x) max(0, x))
    stats_scaled[, "avg_udp_deviation_from_gDNA"] = sapply(stats_scaled[, "avg_udp_deviation_from_gDNA"], function(x) max(0, x))
    stats_nn <- get.knn(stats_scaled, k = max(2, floor(n_cDNA / 2)))
    stats_nnd <- rowMeans(stats_nn$nn.dist)
    stats.new.df[stats.new.df$type == "cDNA",]$outlier_score_knn = stats_nnd
    # Compute isolation forest outlier metric
    stats_forest <- iForest(stats.cDNA.df, nt = 100, phi = n_cDNA)
    forest_score <- predict(stats_forest, newdata = stats.cDNA.df)
    stats.new.df[stats.new.df$type == "cDNA",]$outlier_score_iso <- forest_score
    
    stats.new.df[stats.new.df$type == "cDNA",]$is_outlier = (stats.new.df[stats.new.df$type == "cDNA",]$outlier_score_knn > 3)
  }
  if (n_gDNA >= 4 & !any(is.na(stats.new.df$avg_udp_deviation))) {
    stats.gDNA.df = stats.new.df %>%
      filter(type == "gDNA") %>%
      select(num_udps, HDR_WT_ratio, DEL_WT_ratio, avg_udp_deviation, avg_udp_deviation_from_gDNA)
    stats_scaled <- scale(stats.gDNA.df)
    stats_scaled[is.na(stats_scaled)] = 0
    stats_scaled[, "avg_udp_deviation"] = sapply(stats_scaled[, "avg_udp_deviation"], function(x) max(0, x))
    stats_scaled[, "avg_udp_deviation_from_gDNA"] = sapply(stats_scaled[, "avg_udp_deviation_from_gDNA"], function(x) max(0, x))
    stats_nn <- get.knn(stats_scaled, k = max(2, floor(n_gDNA / 2)))
    stats_nnd <- rowMeans(stats_nn$nn.dist)
    stats.new.df[stats.new.df$type == "gDNA",]$outlier_score_knn = stats_nnd
    
    stats_forest <- iForest(stats.gDNA.df, nt = 100, phi = n_gDNA)
    forest_score <- predict(stats_forest, newdata = stats.gDNA.df)
    stats.new.df[stats.new.df$type == "gDNA",]$outlier_score_iso <- forest_score
    
    stats.new.df[stats.new.df$type == "gDNA",]$is_outlier = (stats.new.df[stats.new.df$type == "gDNA",]$outlier_score_knn > 3)
  }
  
  return( list(stats.df = stats.new.df,
               replicate.udp.fractions.df = replicate.udp.filled.df) )
}

getReplicateQCPlots = function(stats.df, replicate.udp.fractions.df, min_avg_udp_fraction, region_title) {
  color_values = c(`cDNA`="firebrick1", `cDNA outlier`="orange1", `gDNA`="dodgerblue3", `gDNA outlier`="turquoise1")
  shapes = c(16,17,15,3,12,8,6,5,0,1,11,10,18,7,9,2,3,4,13,14)
  num_gDNA_reps = length(unique(replicate.udp.fractions.df %>% dplyr::filter(type == "gDNA") %>% .$replicate))
  num_cDNA_reps = length(unique(replicate.udp.fractions.df %>% dplyr::filter(type == "cDNA") %>% .$replicate))
  shape_values = c(shapes[1:min(20, num_cDNA_reps)], shapes[1:min(20, num_gDNA_reps)])
  if (length(unique(replicate.udp.fractions.df$replicate)) <= 20) {
    p.udp_fractions = ggplot(replicate.udp.fractions.df, aes(x=udp_id, y=udp_fraction, group=replicate, col=type, shape=replicate))
  } else {
    p.udp_fractions = ggplot(replicate.udp.fractions.df, aes(x=udp_id, y=udp_fraction, group=replicate, col=type))
  }
  num_udps = length(unique(replicate.udp.fractions.df$udp))
  p.udp_fractions = p.udp_fractions +
    geom_line() + geom_point(size=3, alpha=0.5) +
    scale_color_manual(values = color_values, guide = F) +
    scale_shape_manual(values = shape_values) +
    theme_bw(10) + xlab("UDP") + ylab("UDP fraction") +
    scale_y_log10() +
    ggtitle(sprintf("UDP fractions (min %.2g%%, %d UDPs)", min_avg_udp_fraction*100, num_udps))
  
  stats.df = stats.df %>% group_by(type) %>%
    mutate(type2 = ifelse(!is.na(is_outlier) & is_outlier, paste(type, "outlier"), type)) %>%
    arrange(replicate)
  stats.df$replicate = fct_reorder(stats.df$replicate, as.integer(factor(stats.df$type)))
  
  p.udp_avg_deviation = ggplot(stats.df %>% filter(!is.na(avg_udp_deviation)),
                               aes(x=replicate, y=avg_udp_deviation, fill=type2)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.2g%%", avg_udp_deviation*100)), size = 2.4, position = position_stack(vjust = 0.5)) +
    theme_bw(10) + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.title=element_blank(), legend.position="none") +
    scale_fill_manual(values = color_values) +
    ylab("UDP frac deviation") +
    ggtitle("Mean UDP fraction deviation")
  
  p.udp_avg_deviation_gDNA = ggplot(stats.df %>% filter(!is.na(avg_udp_deviation_from_gDNA)),
                                    aes(x=replicate, y=avg_udp_deviation_from_gDNA, fill=type2)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.2g%%", avg_udp_deviation_from_gDNA*100)), size = 2.4, position = position_stack(vjust = 0.5)) +
    theme_bw(10) + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.title=element_blank()) +
    scale_fill_manual(values = color_values) +
    ylab("UDP frac deviation") +
    ggtitle("Mean UDP fraction deviation (compared to gDNA)")
  
  # p.udp_avg_deviation_relative = ggplot(stats.df, aes(x=replicate, y=avg_udp_rel_deviation_from_gDNA, fill=type2)) +
  #   geom_bar(stat = "identity", alpha = 0.8) +
  #   geom_text(aes(label = sprintf("%.0f%%", avg_udp_rel_deviation_from_gDNA*100)), size = 2.4, position = position_stack(vjust = 0.5)) +
  #   theme_bw(10) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.title=element_blank()) +
  #   scale_fill_manual(values = color_values) +
  #   ylab("UDP rel. deviation") +
  #   ggtitle("Mean UDP fraction relative deviation (rel. to gDNA UDP mean)")
  
  outlier.df = stats.df
  if (all(is.na(outlier.df$outlier_score_knn))) {
    p.outlier_scores = ggplot(outlier.df, aes(x=replicate, y=outlier_score_knn)) +
      geom_blank() +
      theme_bw(10) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.title=element_blank(), legend.position="none") +
      xlab("Replicate") + ylab("Outlier score") +
      ggtitle("KNN outlier score") +
      annotate("text", x=outlier.df$replicate[1], y=NA, label="Outlier scores are only calculated with sufficient UDPs and >= 4 replicates", hjust = 0)
  } else {
    p.outlier_scores = ggplot(outlier.df, aes(x=replicate, y=outlier_score_knn, fill=type2)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      geom_text(aes(label = sprintf("%.3g", outlier_score_knn)), size = 2.4, position = position_stack(vjust = 0.5)) +
      theme_bw(10) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.title=element_blank(), legend.position="none") +
      scale_fill_manual(values = color_values) +
      xlab("Replicate") + ylab("Outlier score") +
      ggtitle("KNN outlier score")
  }
    
  p.title = ggdraw() + draw_label(sprintf("%s replicate QC", region_title), fontface='bold')
  #p.replicate_qc = plot_grid(p.udp_fractions, p.udp_avg_deviation, ncol=1)
  #p.res = plot_grid(p.title, p.replicate_qc, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
  p.res = egg::ggarrange(p.title, p.udp_fractions, p.udp_avg_deviation, p.udp_avg_deviation_gDNA, p.outlier_scores, ncol=1, heights=c(1.5,4,4,4,4), draw = F)

  return( list(p1 = getReplicateStatsPlot(stats.df, region_title), p2 = p.res) )
}


getUNSDataOld2 = function(replicate.udp.df, replicates.df, sites, region_name, min_gDNA_count = 10, min_cDNA_count = 0) {
  udp_types = unique(replicate.udp.df$type)
  if (!("cDNA" %in% udp_types & "gDNA" %in% udp_types)) {
    return(NULL)
  }
  if (min_gDNA_count < 1) {
    warning(sprintf("getUNSData: min_gDNA_count was set to %d, but must be an integer greater than zero. Setting it to 1.", min_gDNA_count))
    min_gDNA_count = 1
  }
  
  # We need to have the same number of replicates for every UDP. We do this
  # by spreading the replicates out into columns (cDNA_1, cDNA_2, etc.), with
  # fill = 0 for when a replicate is missing, and then gather back into separate
  # replicates
  getType = function(type_rep) { as.character(sapply(type_rep, FUN=function(x) strsplit(x, "^" , fixed=T)[[1]][1])) }
  getReplicate = function(type_rep) { as.character(sapply(type_rep, FUN=function(x) strsplit(x, "^" , fixed=T)[[1]][2])) }
  replicate.udp.filled.df = replicate.udp.df %>%
    dplyr::mutate(type_replicate = paste(type, replicate, sep="^")) %>%
    dplyr::select(-name, -replicate, -type, -avg_seq_length, -avg_mismatch_count) %>%
    tidyr::spread(type_replicate, num_reads, fill = 0) %>%
    tidyr::gather(key="type_replicate", value="num_reads", -(udp:deletion2_end)) %>%
    dplyr::mutate(type = getType(type_replicate),
                  replicate = getReplicate(type_replicate)) %>%
    dplyr::select(-type_replicate)
  
  # We use the same function to calculate stats for all UDPs as we do for the HDR allele.
  # For this we need a column with num_wt_reads.
  replicate.wt_reads.df = replicate.udp.filled.df %>%
    dplyr::filter(is_wt_allele, !is_hdr_allele) %>%
    dplyr::select(replicate, type, num_wt_reads = num_reads) %>%
    dplyr::group_by(replicate, type) %>%
    dplyr::summarise(num_wt_reads = sum(num_wt_reads)) %>%
    ungroup()
  
  # Results not likely to be stable if the WT count is too low
  all_wt_gDNA_count = sum(replicate.wt_reads.df %>% dplyr::filter(type == "gDNA") %>% .$num_wt_reads)
  all_wt_cDNA_count = sum(replicate.wt_reads.df %>% dplyr::filter(type == "cDNA") %>% .$num_wt_reads)
  if (all_wt_gDNA_count < 100) {
    warning(sprintf("%s wild-type gDNA count (%d) is low and may lead to unstable estimates.", region_name, all_wt_gDNA_count))
    if (all_wt_gDNA_count < 1) {
      warning(sprintf("%s wild-type gDNA count is zero!", region_name))
      return(NULL)
    }
    if (all_wt_gDNA_count < min_gDNA_count) {
      warning(sprintf("%s wild-type gDNA count (%d) is below minimum (%d).", region_name, all_wt_gDNA_count, min_gDNA_count))
      return(NULL)
    }
  }
  if (all_wt_cDNA_count < 100) {
    warning(sprintf("%s wild-type cDNA count (%d) is low and may lead to unstable estimates.", region_name, all_wt_cDNA_count))
    if (all_wt_cDNA_count < 1) {
      warning(sprintf("%s wild-type cDNA count is zero!", region_name))
      return(NULL)
    }
  }

  # Spread num_reads out so that we have a column for gDNA and cDNA. For a given
  # replicate/UDP, one of these will be zero in each row.
  replicate.udp.filled.spread.df = replicate.udp.filled.df %>%
    dplyr::mutate(type2 = type) %>%
    spread(type2, num_reads, fill = 0) %>%
    dplyr::rename(gDNA_count = gDNA, cDNA_count = cDNA)

  # Get total read counts by UDP
  udp.total_counts.df = replicate.udp.filled.spread.df %>% 
    group_by(udp) %>%
    dplyr::summarise(gDNA_total_count = sum(gDNA_count),
                     cDNA_total_count = sum(cDNA_count),
                     total_count = gDNA_total_count + cDNA_total_count)
  
  # Following this, we have a row per replicate per UDP, with counts for cDNA
  # and gDNA, which are zero if the UDP wasn't observed in the replicate.
  replicate.dels.df = replicate.udp.filled.spread.df %>%
    dplyr::left_join(udp.total_counts.df, by="udp") %>%
    dplyr::left_join(replicate.wt_reads.df %>% dplyr::select(-type), by="replicate") %>%
    dplyr::filter(is_hdr_allele | is_wt_allele | (has_crispr_deletion & gDNA_total_count >= min_gDNA_count & cDNA_total_count >= min_cDNA_count)) %>%
    dplyr::arrange(!is_wt_allele, !is_hdr_allele, -total_count)
  
  N_cDNA = sum(replicate.wt_reads.df$type == "cDNA")
  N_gDNA = sum(replicate.wt_reads.df$type == "gDNA")
  
  udp.dels.df = replicate.dels.df %>% group_by(udp) %>%
    dplyr::filter(!is_wt_allele) %>%
    summarise(gDNA_total_count = first(gDNA_total_count),
              cDNA_total_count = first(cDNA_total_count),
              total_count = first(total_count),
              mean_gDNA_count = mean(gDNA_count[type == "gDNA"], na.rm = T),
              mean_cDNA_count = mean(cDNA_count[type == "cDNA"], na.rm = T),
              cDNA_ratio = mean(cDNA_count[type == "cDNA"] / num_wt_reads[type == "cDNA"], na.rm = T),
              sd_cDNA_ratio = sd(cDNA_count[type == "cDNA"] / num_wt_reads[type == "cDNA"], na.rm = T),
              se_cDNA_ratio = sd_cDNA_ratio / sqrt(n()),
              gDNA_ratio = mean(gDNA_count[type == "gDNA"] / num_wt_reads[type == "gDNA"], na.rm = T),
              sd_gDNA_ratio = sd(gDNA_count[type == "gDNA"] / num_wt_reads[type == "gDNA"], na.rm = T),
              se_gDNA_ratio = sd_gDNA_ratio / sqrt(n()),
              is_hdr_allele = first(is_hdr_allele),
              is_wt_allele = first(is_wt_allele),
              has_crispr_deletion = first(has_crispr_deletion),
              deletion_start = first(deletion_start),
              deletion_end = first(deletion_end),
              ratio_est = list(getUDPRatioEstimate(cDNA_UDP_counts = cDNA_count[type == "cDNA"], cDNA_wt_counts = num_wt_reads[type == "cDNA"],
                                                   gDNA_UDP_counts = gDNA_count[type == "gDNA"], gDNA_wt_counts = num_wt_reads[type == "gDNA"])),
              uns = ratio_est[[1]]$effect,
              uns_se = ratio_est[[1]]$effect_sd,
              uns_df_est = welch_df_estimate(sd_cDNA_ratio, N_cDNA, sd_gDNA_ratio, N_gDNA),
              uns_confint = uns_se * qt(0.975, df = uns_df_est, lower.tail=T)) %>%
    dplyr::select(-ratio_est) %>%
    dplyr::arrange(!is_wt_allele, !is_hdr_allele, -total_count)

  return(list(udp.dels.df = udp.dels.df, replicate.dels.df = replicate.dels.df))
}

getUNSData = function(replicate.udp.df, replicates.df, sites, region_name, min_gDNA_count = 10, min_cDNA_count = 0, max_udps = 10000) {
  udp_types = unique(replicate.udp.df$type)
  if (!("cDNA" %in% udp_types & "gDNA" %in% udp_types)) {
    return(NULL)
  }
  if (min_gDNA_count < 1) {
    warning(sprintf("getUNSData: min_gDNA_count was set to %d, but must be an integer greater than zero. Setting it to 1.", min_gDNA_count))
    min_gDNA_count = 1
  }
  
  # We need to have the same number of replicates for every UDP. We do this
  # by spreading the replicates out into columns (cDNA_1, cDNA_2, etc.), with
  # fill = 0 for when a replicate is missing, and then gather back into separate
  # replicates
  getType = function(type_rep) { as.character(sapply(type_rep, FUN=function(x) strsplit(x, "^" , fixed=T)[[1]][1])) }
  getReplicate = function(type_rep) { as.character(sapply(type_rep, FUN=function(x) strsplit(x, "^" , fixed=T)[[1]][2])) }
  replicate.udp.filled.df = replicate.udp.df %>%
    dplyr::mutate(type_replicate = paste(type, replicate, sep="^")) %>%
    dplyr::select(udp, num_reads, is_hdr_allele, is_wt_allele, has_crispr_deletion, deletion_start, deletion_end, type_replicate) %>%
    tidyr::spread(type_replicate, num_reads, fill = 0) %>%
    tidyr::gather(key="type_replicate", value="num_reads", -(udp:deletion_end)) %>%
    dplyr::mutate(type = getType(type_replicate),
                  replicate = getReplicate(type_replicate)) %>%
    dplyr::select(-type_replicate)
  
  # We use the same function to calculate stats for all UDPs as we do for the HDR allele.
  # For this we need a column with num_wt_reads.
  replicate.read_counts.df = replicate.udp.filled.df %>%
    dplyr::select(replicate, type, is_wt_allele, num_reads) %>%
    dplyr::group_by(replicate, type) %>%
    dplyr::summarise(replicate_num_reads = sum(num_reads),
                     num_wt_reads = sum(num_reads[is_wt_allele])) %>%
    ungroup()
  
  # Results not likely to be stable if the WT count is too low
  all_wt_gDNA_count = sum(replicate.read_counts.df %>% dplyr::filter(type == "gDNA") %>% .$num_wt_reads)
  all_wt_cDNA_count = sum(replicate.read_counts.df %>% dplyr::filter(type == "cDNA") %>% .$num_wt_reads)
  if (all_wt_gDNA_count < 100) {
    warning(sprintf("%s wild-type gDNA count (%d) is low and may lead to unstable estimates.", region_name, all_wt_gDNA_count))
    if (all_wt_gDNA_count < 1) {
      warning(sprintf("%s wild-type gDNA count is zero!", region_name))
      return(NULL)
    }
    if (all_wt_gDNA_count < min_gDNA_count) {
      warning(sprintf("%s wild-type gDNA count (%d) is below minimum (%d).", region_name, all_wt_gDNA_count, min_gDNA_count))
      return(NULL)
    }
  }
  if (all_wt_cDNA_count < 100) {
    warning(sprintf("%s wild-type cDNA count (%d) is low and may lead to unstable estimates.", region_name, all_wt_cDNA_count))
    if (all_wt_cDNA_count < 1) {
      warning(sprintf("%s wild-type cDNA count is zero!", region_name))
      return(NULL)
    }
  }
  
  # Get total read counts by UDP, for gDNA and cDNA
  udp.total_counts.df = replicate.udp.filled.df %>% 
    group_by(udp) %>%
    dplyr::summarise(gDNA_total_count = sum(num_reads[type == "gDNA"]),
                     cDNA_total_count = sum(num_reads[type == "cDNA"]),
                     total_count = gDNA_total_count + cDNA_total_count)
  
  # Following this, we have a row per replicate per UDP, with counts for cDNA
  # and gDNA, which are zero if the UDP wasn't observed in the replicate.
  replicate.dels.df = replicate.udp.filled.df %>%
    dplyr::left_join(udp.total_counts.df, by="udp") %>%
    dplyr::left_join(replicate.read_counts.df %>% dplyr::select(-type), by="replicate") %>%
    dplyr::filter(gDNA_total_count >= min_gDNA_count & cDNA_total_count >= min_cDNA_count) %>%
    dplyr::arrange(!is_wt_allele, !is_hdr_allele, desc(total_count))
  
  N_cDNA = sum(replicate.read_counts.df$type == "cDNA")
  N_gDNA = sum(replicate.read_counts.df$type == "gDNA")
  
  # Summarize details per uDP
  udp.dels.df = replicate.dels.df %>% group_by(udp) %>%
    dplyr::filter(!is_wt_allele) %>%
    summarise(gDNA_total_count = first(gDNA_total_count),
              cDNA_total_count = first(cDNA_total_count),
              total_count = first(total_count),
              mean_cDNA_count = mean(num_reads[type == "cDNA"], na.rm = T),
              mean_gDNA_count = mean(num_reads[type == "gDNA"], na.rm = T),
              cDNA_ratio = mean(num_reads[type == "cDNA"] / num_wt_reads[type == "cDNA"], na.rm = T),
              sd_cDNA_ratio = sd(num_reads[type == "cDNA"] / num_wt_reads[type == "cDNA"], na.rm = T),
              se_cDNA_ratio = sd_cDNA_ratio / sqrt(n()),
              gDNA_ratio = mean(num_reads[type == "gDNA"] / num_wt_reads[type == "gDNA"], na.rm = T),
              sd_gDNA_ratio = sd(num_reads[type == "gDNA"] / num_wt_reads[type == "gDNA"], na.rm = T),
              se_gDNA_ratio = sd_gDNA_ratio / sqrt(n()),
              is_hdr_allele = first(is_hdr_allele),
              is_wt_allele = first(is_wt_allele),
              has_crispr_deletion = first(has_crispr_deletion),
              deletion_start = first(deletion_start),
              deletion_end = first(deletion_end)) %>%
    dplyr::arrange(!is_wt_allele, !is_hdr_allele, desc(total_count))
  
  # Get UNS and associated stats per UDP
  denom = "num_wt_reads"
  if (opts$ratio_to_total_reads) {
    denom = "replicate_num_reads"
  }
  udp.uns_stats.df = replicate.dels.df %>%
    dplyr::filter(!is_wt_allele, udp %in% udp.dels.df$udp[1:max_udps]) %>%
    group_by(udp) %>%
    do( getUDPRatioEstimate(., replicates.df, numerator = "num_reads", denominator = denom, batchCol = opt$batch_col, randomEffectsCols = opt$random_effects_cols) )
  
  udp.dels.df = udp.dels.df %>%
    dplyr::left_join(udp.uns_stats.df %>% rename(uns = effect, uns_se = effect_sd, uns_confint_lo = effect_confint_lo,
                                                 uns_confint_hi = effect_confint_hi, uns_df_est = df_estimated),
                     by="udp") %>%
    dplyr::arrange(!is_wt_allele, !is_hdr_allele, -total_count)
  
  return(list(udp.dels.df = udp.dels.df, replicate.dels.df = replicate.dels.df))
}

getUNSPlot = function(replicate.udp.df, replicates.df, sites, plot_title, min_gDNA_count = 10, min_cDNA_count = 0, max_udps = 40) {
  udp_types = unique(replicate.udp.df$type)
  if (!("cDNA" %in% udp_types & "gDNA" %in% udp_types)) {
    return(egg::ggarrange(textPlot("Cannot make UNS plot without both cDNA and gDNA replicates"), top=plot_title, draw = F))
  }
  
  uns_res = getUNSData(replicate.udp.df, replicates.df, sites, region_name=plot_title, min_gDNA_count, min_cDNA_count, max_udps = max_udps)
  udp.dels.df = uns_res$udp.dels.df %>% dplyr::filter(!is_wt_allele)
  udp.del_replicates.df = uns_res$replicate.dels.df %>% dplyr::filter(!is_wt_allele)
  
  numUDPs = nrow(udp.dels.df)
  if (numUDPs < 2) {
    return(egg::ggarrange(textPlot("No UDPs pass the thresholds for min read counts."), top=plot_title, draw = F))
  }
  if (!is.na(max_udps) & numUDPs > max_udps) {
    udp.dels.df = udp.dels.df %>% .[1:max_udps,]
    udp.del_replicates.df = udp.del_replicates.df %>% dplyr::filter(udp %in% udp.dels.df$udp)
    numUDPs = max_udps
  }
  
  # matrix containing the deletion binary code
  M = as.matrix(do.call(rbind, lapply(as.list(udp.dels.df$udp), udpToBinary)))
  M_chars = as.matrix(do.call(rbind, lapply(as.list(udp.dels.df$udp), strToChars)))
  
  model <- hclust(dist(M))
  dhc <- as.dendrogram(model)
  ddata <- dendro_data(dhc, type = "rectangle")
  dendro_span = max(ddata$segments$y, ddata$segments$yend) - min(ddata$segments$y, ddata$segments$yend)
  
  udp.order.map = data.frame(y = label(ddata)$x, udp = udp.dels.df[model$order,]$udp)
  udp.dels.df = udp.dels.df[model$order,] %>% dplyr::left_join(udp.order.map, by="udp")
  
  plot.df = as.data.frame(M_chars[model$order,])
  colnames(plot.df) = as.character(c(1:ncol(plot.df)))
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
    theme_bw() + theme(axis.title.y = element_blank(),
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
      geom_point(aes(x = pos, y = id), size = dot_size, shape = 19, alpha = 0.7, data = plot.gather.df[plot.gather.df$udpchar == '-' & plot.gather.df$`HDR allele`,]) +
      geom_point(aes(x = pos, y = id), size = dash_size, shape = '-', alpha = 0.8, data = plot.gather.df[plot.gather.df$udpchar == '*' & plot.gather.df$`HDR allele`,])
  } else {
    p.udp = ggplot()
  }
  p.udp = p.udp + 
    geom_point(aes(x = pos, y = id), size = dot_size, shape = 19, alpha = 0.7, data = plot.gather.df[plot.gather.df$udpchar == '-',]) +
    geom_point(aes(x = pos, y = id), size = dash_size, shape = '-', alpha = 0.8, data = plot.gather.df[plot.gather.df$udpchar == '*',]) +
    scale_color_manual(values = c("grey20", "green4")) +
    geom_point(aes(x = pos, y = id, shape = udpchar), color = "black", size = 2.5, data = plot.gather.df[!(plot.gather.df$udpchar == '*' | plot.gather.df$udpchar == '-'),]) + scale_shape_identity() +
    coord_cartesian(ylim=c(min(plot.gather.df$id), max(plot.gather.df$id)), xlim=c(sites$start, sites$end)) +
    xlab("Position") + scale_x_continuous(expand = c(0.01, 0)) + 
    theme_bw() + theme(legend.position = "none",
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
  uns.plot.df = data.frame(y = label(ddata)$x,
                           x = 0,
                           uns = udp.dels.df$uns,
                           cDNA_lower = udp.dels.df$uns < 1,
                           gDNA = udp.dels.df$gDNA_total_count,
                           uns_conf_hi = udp.dels.df$uns_confint_lo,
                           uns_conf_lo = udp.dels.df$uns_confint_hi)
  
  # uns.replicates.plot.df = udp.del_replicates.df %>% 
  #   dplyr::select(udp, uns, gDNA_count) %>%
  #   dplyr::left_join(udp.order.map, by="udp")
  
  maxDotSize = 5.5
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
  uns.plot.df$`cDNA expr` = "higher"
  uns.plot.df$`cDNA expr`[uns.plot.df$cDNA_lower] = "lower"
  
  p.uns = ggplot(uns.plot.df) +
    #annotate("segment", x = uns.plot.df$x, xend = uns.plot.df$uns, y = uns.plot.df$y, yend = uns.plot.df$y, colour = 'gray') +
    geom_errorbarh(aes(y = y, x = uns, xmin = uns_conf_lo, xmax = uns_conf_hi, height = 0.4), colour = "grey70") +
    geom_point(aes(x = uns, y = y, colour = `cDNA expr`, size = log10(gDNA), alpha = log10(gDNA))) +
    scale_size(range = c(1, maxDotSize)) +
    scale_color_manual(values=c(higher="red", lower="blue")) +
    scale_alpha_continuous(range = c(0.2, 1)) + 
    scale_x_continuous(trans="log2", breaks = c(0.5, 1, 2, 4, 8)) +
    #scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, limits = c(0.1, 4)) +
    #geom_point(data = uns.replicates.plot.df, mapping = aes(y = y, x = uns, size = 2, alpha = 0.6)) +
    scale_y_continuous(position = "right") +
    coord_cartesian(xlim = c(min_uns_display, max_uns_display),
                    ylim=c(min(plot.gather.df$id), max(plot.gather.df$id))) +
    ylab("Unique deletion profile") + xlab("UNS") +
    geom_vline(xintercept = 1.0, alpha = 0.25, linetype = "longdash") +
    theme_bw() + theme(legend.position = "right",
                       #axis.line = element_blank(),
                       panel.border = element_blank(),
                       plot.margin = unit(c(0,0,0.05,0), "cm"))
  
  p.udp_profile = egg::ggarrange(p.dendro, p.udp, p.uns, nrow=1, ncol=3, widths=c(0.6,4,1), top = plot_title, draw = F)
  #p.udp_profile = plot_grid(p.uns, p.udp, p.dendro, nrow=1, ncol=3, rel_widths = c(2,4,1))
  return(p.udp_profile)
}


fitVarianceComponents = function(replicate.udp.spread.df, replicates.type.df) {
  replicate.udp.spread.df = replicate.udp.spread.df %>% as.data.frame()
  rownames(replicate.udp.spread.df) = replicate.udp.spread.df$udp
  
  # Get all columns of the replicate dataframe which begin with "replicate",
  # but ensure that they have more than one factor level
  replicate_cols = colnames(replicates.type.df)[grepl("^replicate_", colnames(replicates.type.df))]
  remove_cols = c()
  for (col in replicate_cols) {
    replicates.type.df[, col] = as.factor(replicates.type.df[, col, drop=T])
    if (length(unique(replicates.type.df[, col, drop=T])) <= 1) {
      warning(sprintf("Column %s cannot be used in variance components analysis since it has only one level for region %s.", col, replicates.type.df$name[1]))
      remove_cols = c(col, remove_cols)
    }
  }
  replicate_cols = setdiff(replicate_cols, remove_cols)
  if (length(replicate_cols) == 0) {
    warning(sprintf("Region %s: no replicate variables are available to use in variance components analysis.", replicates.type.df$name[1]))
    return(NULL)
  }
  
  formulaStr = sprintf("~ (1|%s)", replicate_cols[1])
  for (col in replicate_cols[-1]) {
    formulaStr = sprintf("%s + (1|%s)", formulaStr, col)
  }
  
  # IMPORTANT: Make sure that the columns of the expression input to
  # fitExtractVarPartMoel are in the same order of samples as the rows
  # of the "data" input!
  replicate.udp.expr = replicate.udp.spread.df[, replicates.type.df$replicate]
  
  vp <- fitExtractVarPartModel(exprObj = replicate.udp.expr,
                               formula = as.formula(formulaStr),
                               data = replicates.type.df,
                               useWeights = F)
  #View(data.frame(vp))
  vp
}

doVariancePartitionPlot = function(vp.df, residuals=T, pointColor = NA) {
  vp.df = vp.df %>%
    gather(variable, variance, -udp, -name, -frac) %>%
    mutate(pctvariance = variance*100)
  
  if (!residuals) {
    vp.df = vp.df %>% dplyr::filter(variable != "Residuals")
  }
  
  medians = vp.df %>% group_by(variable) %>% summarise(median=median(pctvariance)) %>%
    arrange(-median) %>% filter(variable != "Residuals")
  vp.df$variable = factor(as.character(vp.df$variable),
                          levels=c(medians$variable, "Residuals"))
  numVariables = length(levels(vp.df$variable))
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  getFractionCategory = function(x) {
    if (x < 0.005) { "< 0.5%" }
    else if (x < 0.02) { "< 2%" }
    else { "> 2%" }
  }
  vp.df$fraction = sapply(vp.df$frac, getFractionCategory)

  plotaes = aes(color = "black")
  if (!is.na(pointColor)) {
    if (pointColor == "fraction") {
      plotaes = aes(color = fraction)
    } else {
      plotaes = aes(color = pointColor)
    }
  }
  p = ggplot(vp.df, aes(x=variable, y=pctvariance, fill="grey")) +
    geom_violin(scale="width") +
    geom_boxplot(width=0.2, fill="grey80", outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.75, mapping = aes(color = fraction)) +
    theme_bw() + guides(fill=FALSE) +
    theme(legend.position="none") + 
    theme(axis.text.x=element_text(angle=20, hjust = 1), axis.title.x=element_blank(), panel.grid=element_blank()) +
    ylab("Variance explained (%)") +
    ylim(c(0,100)) +
    scale_fill_manual(values = c("gray95"))
  #  scale_fill_manual(values = c(gg_color_hue(numVariables-1), "gray85"))
  if (!is.na(pointColor) & pointColor == "fraction") {
    #p = p + scale_color_manual(values = c(rep(color, numVariables)))
    p = p + scale_color_manual(values = c("< 0.5%"="chartreuse", "< 2%"="blue", "> 2%"="red2"))
  }
  return(p)
}

getVarianceComponents = function(replicate.udp.df, replicates.df, method = "read_fraction", min_udp_total_count = 100, min_udp_fraction = 0.001) {
  replicate_cols = colnames(replicates.df)[grepl("^replicate_", colnames(replicates.df))]
  if (length(replicate_cols) < 1) {
    stop("--variance_analysis option given but no columns beginning with 'replicate_' found in the --replicates file")
  }
  # Check that each replicate column has more than 1 distinct value,
  # and fewer distinct values than the number of replicates
  isValidCol = function(col) {
    n = length(unique(replicates.df[,col,drop=T]))
    return(n > 1 & n < nrow(replicates.df))
  }
  validCols = sapply(replicate_cols, isValidCol)
  if (any(!validCols)) {
    warning(sprintf("Not including the following replicate column(s) in variance analysis: %s",
                    paste(replicate_cols[!validCols], collapse = ", ")))
    replicate_cols = replicate_cols[validCols]
  }
  # Ensure that replicate cols are treated as factors
  for (col in replicate_cols) {
    replicates.df[, col] = as.factor(replicates.df[, col, drop=T])
  }
  
  # Subset replicate info dataframe to those replicates in replicate.udp.df
  replicates.df = replicates.df %>% dplyr::filter(replicate %in% unique(replicate.udp.df$replicate))
  
  filterUDPs = function(udp.df, min_udp_total_count, min_udp_fraction) {
    # Spread the replicates out into columns (cDNA_1, cDNA_2, etc.), with
    # fill = 0 for when a replicate has no value for a given UDP, then
    # gather again.
    udp.filled.df = udp.df %>%
      dplyr::select(udp, type, replicate, num_reads) %>%
      tidyr::spread(replicate, num_reads, fill = 0) %>%
      tidyr::gather(key="replicate", value="num_reads", -udp, -type)
    totalReads = sum(udp.filled.df$num_reads)
    udp.summary = udp.filled.df %>% group_by(udp, type) %>%
      dplyr::summarise(sum = sum(num_reads), frac = sum / totalReads)
    
    # Remove UDPs with too few total reads
    udp.filled.df %>% dplyr::inner_join(udp.summary %>% dplyr::filter(sum >= min_udp_total_count, frac > min_udp_fraction), by=c("udp", "type"))
  }
  
  # This function returns a value based on the "num_reads" for each replicate,
  # but which is what we use to determine variance components. Raw number of
  # reads is not really suitable. The recommended value is "read_fraction",
  # where this represents the fraction of reads that a UDP represents for a
  # given replicate.
  getValueForVariance = function(replicate.udp.df, method) {
    if (method == "read_fraction") {
      replicate.numreads = replicate.udp.df %>%
        dplyr::group_by(replicate) %>%
        dplyr::summarise(replicate_num_reads = sum(num_reads))
      replicate.udp.df = replicate.udp.df %>%
        dplyr::left_join(replicate.numreads, by="replicate") %>%
        dplyr::mutate(udp_fraction = num_reads / replicate_num_reads)
      value = replicate.udp.df$udp_fraction
    }
    else if (method == "reads") {
      value = replicate.udp.df$num_reads
    }
    else if (method == "log_reads") {
      value = log(replicate.udp.df$num_reads + 1)
    }
    else if (method == "deviation_from_poisson") {
      mean_reads_per_udp.df = replicate.udp.df %>% dplyr::group_by(udp) %>%
        dplyr::summarise(mean_reads = mean(num_reads))
      replicate.udp.df = replicate.udp.df %>% dplyr::left_join(mean_reads_per_udp.df, by="udp")
      value = abs(replicate.udp.df$num_reads - replicate.udp.df$mean_reads) / sqrt(replicate.udp.df$mean_reads)
    } else {
      stop(sprintf("Unrecognized variance plots method: '%s'", method))
    }
    value
  }
  
  # Function to get variance components result for cDNA and gDNA for a single region
  getRegionVarianceComponents = function(replicate.udp.df, replicates.type.df, dnaType, min_udp_total_count, min_udp_fraction) {
    region = replicate.udp.df$name[1]
    replicates.type.df %>% dplyr::filter(name == region, type == dnaType)
    replicate.udp.type.df = filterUDPs(replicate.udp.df %>% dplyr::filter(type == dnaType),
                                       min_udp_total_count = min_udp_total_count, min_udp_fraction = min_udp_fraction)
    
    replicate.udp.type.df$value = getValueForVariance(replicate.udp.type.df, method)
    
    # Spread the replicates out into columns (cDNA_1, cDNA_2, etc.), with
    # fill = 0 for when a replicate has no value for a given UDP
    replicate.udp.type.spread.df = replicate.udp.type.df %>%
      dplyr::select(udp, frac, replicate, value) %>%
      tidyr::spread(replicate, value, fill = 0)
    udpFraction = rowMeans(replicate.udp.type.spread.df %>% dplyr::select(-udp))
    vp = fitVarianceComponents(replicate.udp.spread.df = replicate.udp.type.spread.df %>% filter(frac > 0.001) %>% select(-frac),
                               replicates.type.df = replicates.type.df)
    vp.df = NULL
    if (!is.null(vp)) {
      vp.df = data.frame(vp) %>% 
        mutate(udp = rownames(vp)) %>%
        left_join(replicate.udp.type.spread.df %>% select(udp, frac), by="udp")
    }
    vp.df
  }
  
  #########################################################################
  dnaType = "cDNA"
  # Get variance components separately for each region (if there's more than one)
  # by.vp.dfs = by(replicate.udp.df, replicate.udp.df$name, FUN = getRegionVarianceComponents,
  #                replicates.type.df = replicates.df %>% dplyr::filter(type == dnaType),
  #                dnaType = dnaType, min_udp_total_count, min_udp_fraction)
  # # Merge each region's variance components result together
  # for (i in 1:length(by.vp.dfs)) {
  #   by.vp.dfs[[i]]$name = names(by.vp.dfs)[i]
  # }
  # vp.cDNA.df = bind_rows(lapply(by.vp.dfs, function(x) x))
  
  regions = unique(replicate.udp.df$name)
  dflist = list()
  for (i in 1:length(regions)) {
    dflist[[i]] = getRegionVarianceComponents(replicate.udp.df %>% filter(name == regions[i]),
                                           replicates.df %>% dplyr::filter(name == regions[i], type == dnaType),
                                           dnaType, min_udp_total_count, min_udp_fraction)
    dflist[[i]]$name = regions[i]
  }
  vp.cDNA.df = bind_rows(dflist)
  
  # Do the same for gDNA
  dnaType = "gDNA"
  # by.vp.dfs = by(replicate.udp.df, replicate.udp.df$name, FUN = getRegionVarianceComponents,
  #                replicates.type.df = replicates.df %>% dplyr::filter(type == dnaType),
  #                dnaType = dnaType, min_udp_total_count, min_udp_fraction)
  # for (i in 1:length(by.vp.dfs)) {
  #   by.vp.dfs[[i]]$name = names(by.vp.dfs)[i]
  # }
  # vp.gDNA.df = bind_rows(lapply(by.vp.dfs, function(x) x))
  regions = unique(replicate.udp.df$name)
  dflist = list()
  for (i in 1:length(regions)) {
    dflist[[i]] = getRegionVarianceComponents(replicate.udp.df %>% filter(name == regions[i]),
                                              replicates.df %>% dplyr::filter(name == regions[i], type == dnaType),
                                              dnaType, min_udp_total_count, min_udp_fraction)
    dflist[[i]]$name = regions[i]
  }
  vp.gDNA.df = bind_rows(dflist)
  
  return(list(vp.cDNA = vp.cDNA.df, vp.gDNA = vp.gDNA.df))
}

getVarianceComponentsPlots = function(vp.cDNA.df, vp.gDNA.df, min_udp_total_count, min_udp_fraction, method = "read_fraction", plot_title = NULL) {
  if (nrow(vp.cDNA.df) <= 0) {
    p.cDNA = textPlot("Cannot make variance components plot.")
    cDNA.summary.df = data.frame(x="empty")
    cDNA_note = ""
  } else {
    p.cDNA = doVariancePartitionPlot(vp.cDNA.df, pointColor="fraction") + ggtitle("cDNA")
    cDNA.summary.df = variancePartitionSummaryTable(vp.cDNA.df %>% dplyr::select(-name, -udp, -frac))
    cDNA_note = sprintf("Based on %d UDPs with >%d reads, >%.1g%% udp fraction", nrow(vp.cDNA.df), min_udp_total_count, min_udp_fraction*100)
  }
  
  if (nrow(vp.gDNA.df) <= 0) {
    p.gDNA = textPlot("Cannot make variance components plot.")
    gDNA.summary.df = data.frame(x="empty")
    gDNA_note = ""
  } else {
    p.gDNA = doVariancePartitionPlot(vp.gDNA.df, pointColor="fraction") + ggtitle("gDNA") + theme(legend.position="left")
    gDNA.summary.df = variancePartitionSummaryTable(vp.gDNA.df %>% dplyr::select(-name, -udp, -frac))
    gDNA_note = sprintf("Based on %d UDPs with >%d reads, >%.1g%% udp fraction", nrow(vp.gDNA.df), min_udp_total_count, min_udp_fraction*100)
  }
  
  ggThemeBlank = theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                                    axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
                                    panel.border = element_blank())
  myTableTheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.75)),
                                            colhead = list(fg_params=list(cex = 0.75)),
                                            rowhead = list(fg_params=list(cex = 0.75)))
  
  if (is.null(plot_title)) {
    plot_title = "Variance components"
  }
  p.stats = ggplot(data.frame(x=1:10, y=1:10), aes(x,y)) + geom_blank() + ggThemeBlank +
    annotation_custom(tableGrob(cDNA.summary.df, theme = myTableTheme), xmin=1, xmax=5, ymin=2, ymax=8.5) +
    annotation_custom(tableGrob(gDNA.summary.df, theme = myTableTheme), xmin=6, xmax=10, ymin=2, ymax=8.5) +
    # annotate("text", x=3.5, y=8, label = "cDNA", hjust = 0.5, fontface = 2, size = 5) +
    # annotate("text", x=7.5, y=8, label = "gDNA", hjust = 0.5, fontface = 2, size = 5) +
    annotate("text", x=5, y=10, label = plot_title, vjust = 1, fontface = 2, size = 5) +
    annotate("text", x=5, y=9, label = sprintf("Method = %s", method), vjust = 1, size = 3) +
    annotate("text", x=1, y=1, label = cDNA_note, hjust = 0, vjust = 0, size = 3) +
    annotate("text", x=6, y=1, label = gDNA_note, hjust = 0, vjust = 0, size = 3)
  
  #egg::ggarrange(p.stats, grid.arrange(p.cDNA, p.gDNA, ncol=2), nrow = 1)
  return( plot_grid(p.stats, plot_grid(p.cDNA, p.gDNA, nrow = 1, rel_widths=c(0.46, 0.54)), nrow = 2, rel_heights = c(0.4, 0.6)) )
}

variancePartitionSummaryTable = function(vp.df) {
  # Get some summary stats for each column of the variance partition data.frame
  col = colnames(vp.df)[1]
  summary.df = summary(vp.df[, col])
  for (col in colnames(vp.df)[-1]) {
    summary.df = rbind(summary.df, summary(vp.df[, col]))
  }
  summary.df = as.data.frame(summary.df)
  rownames(summary.df) = colnames(vp.df)
  for (col in colnames(summary.df)) {
    summary.df[,col] = sapply(summary.df[,col], FUN = function(x) sprintf("%.3g", x))
  }
  summary.df[, c("1st Qu.", "Median", "3rd Qu.")]
}


getUDPStats = function(replicate.udp.df, dnaType, min_udp_total_count) {
  # Spread the replicates out into columns (cDNA_1, cDNA_2, etc.), with
  # fill = 0 for when a replicate has no value for a given UDP, then
  # gather again
  replicate.udp.filled.df = replicate.udp.df %>%
    dplyr::filter(type == dnaType) %>%
    dplyr::select(udp, replicate, num_reads) %>%
    tidyr::spread(replicate, num_reads, fill = 0)
  
  replicate.udp.filled.df$total_count = apply(replicate.udp.filled.df[,-1], MARGIN = 1, FUN = sum)
  replicate.udp.filled.df = replicate.udp.filled.df %>%
    dplyr::filter(total_count > min_udp_total_count) %>%
    dplyr::select(-total_count)
  replicate.udp.filled.df = replicate.udp.filled.df %>%
    tidyr::gather(key="replicate", value="num_reads", -udp)
  
  replicate_reads = replicate.udp.filled.df %>% dplyr::group_by(replicate) %>% dplyr::summarise(replicate_reads = sum(num_reads))
  replicate.udp.filled.df = replicate.udp.filled.df %>%
    dplyr::left_join(replicate_reads, by="replicate")
  mean_replicate_reads = mean(replicate_reads$replicate_reads)
  stats.df = replicate.udp.filled.df %>% dplyr::group_by(udp) %>%
    dplyr::summarise(type = dnaType,
                     udp_mean = mean(num_reads),
                     udp_sd = sd(num_reads),
                     udp_cv = udp_sd / udp_mean,
                     udp_frac_mean = mean(num_reads) / mean_replicate_reads,
                     udp_frac_sd = sd(num_reads / replicate_reads),
                     udp_frac_cv = udp_frac_sd / udp_frac_mean)
  stats.df
}


getPowerPlots = function(replicate.udp.df, replicates.df, titlestr = "", min_udp_total_count = 100) {
  # We plot power against %editing, i.e. the % of all reads represented by a
  # given edit. The %editing can be considered equivalent to mean UDP read
  # count across replicates.
  # To determine power, we need the standard deviation expected across cDNA
  # replicates for a given %editing. We can get this by fitting a curve to
  # the empirical curve of coefficient of variation (CV) vs. UDP read count.
  powerPlotTitle = sprintf("%s power estimate", titlestr)
  
  # An error may occur in NLS, so we need to catch it
  res = tryCatch({
    # First fit a curve to a table of CV vs. mean read count.
    ## TODO: consider changing the modeling form to a free exponent, e.g.
    ## udp_frac_cv ~ a + b * udp_frac_mean^c
    ## along with constraints such as a > 0, b > 0, -1 <= c <= 0
    cDNA.stats.df = getUDPStats(replicate.udp.df, "cDNA", min_udp_total_count)
    fit.cDNA <- nls(udp_frac_cv ~ a + b * (1/sqrt(udp_frac_mean)), data = cDNA.stats.df,
                    weights = log(cDNA.stats.df$udp_mean), start = list(a = 0, b = 1))
    beta_cDNA <- coefficients(fit.cDNA)
    
    gDNA.stats.df = getUDPStats(replicate.udp.df, "gDNA", min_udp_total_count)
    fit.gDNA <- nls(udp_frac_cv ~ a + b * (1/sqrt(udp_frac_mean)), data = gDNA.stats.df,
                    weights = log(gDNA.stats.df$udp_mean), start = list(a = 0, b = 1))
    beta_gDNA <- coefficients(fit.gDNA)
  }
  , error = function(e) e, warning = function(w) w)
  if (is(res, "error")) {
    warning(sprintf("%s: Error in NLS fitting to UDP CV vs. UDP fraction.\nPower plots cannot be generated.", titlestr))
    plot_list = list(cv_plots = NULL,
                     power = ggarrange(textPlot("Error in NLS fitting to UDP CV vs. UDP fraction.\nPower plots cannot be generated."),
                                       top=powerPlotTitle, draw = F),
                     replicate_allocation = NULL)
    return(plot_list)
  }
  
  fitcDNAUDP_CV = function(udp_frac_mean) { beta_cDNA["a"] + beta_cDNA["b"]/sqrt(udp_frac_mean) }
  fitgDNAUDP_CV = function(udp_frac_mean) { beta_gDNA["a"] + beta_gDNA["b"]/sqrt(udp_frac_mean) }
  
  # Make some plots that show the relationship between UDP read fraction and CV
  # Below I plot the UDP read % (fraction of total reads for each UDP) against CV
  # I've also plotted UDP read count against CV, but read fraction makes more sense
  # since it controls for the number of reads in each replicate. (For very low UDP
  # read counts it might make sense to use read count though.)
  stats.df = rbind(cDNA.stats.df, gDNA.stats.df)
  p.udpfrac_cv = ggplot(stats.df, aes(x=udp_frac_mean * 100, y=udp_cv, col = type)) + geom_point(alpha = 0.7) +
    scale_x_log10(breaks = c(0.1, 1, 10)) +
    theme_bw() + xlab("Mean UDP read %") + ylab("UDP coeff of variation") +
    stat_function(fun = function(x) 1/sqrt(x / 100 * sum(gDNA.stats.df$udp_mean)), color = "grey50", size = 1) + # poisson rate of CV
    stat_function(fun = function(x) {fitcDNAUDP_CV(x / 100)}, color = "red", alpha = 0.5, size = 1.3) +
    stat_function(fun = function(x) {fitgDNAUDP_CV(x / 100)}, color = "blue", alpha = 0.5, size = 1.3) +
    scale_color_manual(values = c(cDNA = "red", gDNA = "blue"))
  
  # Add a label showing the equations of the fit
  fit_label = sprintf("cDNA: %.3g + %.3g/sqrt(x)\ngDNA: %.3g + %.3g/sqrt(x)",
                      beta_cDNA["a"], beta_cDNA["b"], beta_gDNA["a"], beta_gDNA["b"])
  p.udpfrac_cv = p.udpfrac_cv + annotate("text", Inf, Inf, label = fit_label, hjust = 1, vjust = 1, size = 2.8, color = "grey30")
  
  p.cv_density = ggplot(stats.df, aes(x = udp_frac_cv, color = type)) + geom_density() +
    scale_color_manual(values = c(cDNA = "red", gDNA = "blue")) +
    theme_bw() + xlab("UDP frac CV") + ylab("Density")
  p.cv_plots = egg::ggarrange(p.udpfrac_cv, p.cv_density, ncol=1, heights = c(2,1),
                              top = sprintf("%s variance summary", titlestr), draw = F)
  
  # The general formula for propagation of uncertainty is (where f = A/B):
  #     sd(f) = f * sqrt( (sd(A) / A)^2 + (sd(B) / B)^2 - 2*cov(A,B) )
  # In some cases we may be able to assume that cov(A,B) is zero, e.g.
  # where A and B are independent replicates. Note that sd(A)/A is the
  # coefficient of variation, which is why we use CV below. Frac is the
  # fraction of all reads that come from one UDP, so it is basically A/B,
  # where A is number for one UDP and B is the total number of reads.
  sd_cDNA_udp_frac = function(frac) { frac * sqrt(fitcDNAUDP_CV(frac)^2 + fitcDNAUDP_CV(1-frac)^2) }
  sd_gDNA_udp_frac = function(frac) { frac * sqrt(fitgDNAUDP_CV(frac)^2 + fitgDNAUDP_CV(1-frac)^2) }
  
  n_cDNA_rep = length(replicates.df %>% dplyr::filter(type == "cDNA") %>% .$replicate)
  n_gDNA_rep = length(replicates.df %>% dplyr::filter(type == "gDNA") %>% .$replicate)
  # Here we are computing the SD of the cDNA / gDNA ratio estimate, where this
  # is based on all replicates done for each. Hence, we use the standard error
  # for each of these estimates (cDNA:WT ratio and gDNA:WT ratio). The SE of
  # these estimates is the SD / sqrt(n).
  sd_ratio_cDNA_gDNA = function(cDNA_frac, gDNA_frac, n_cDNA, n_gDNA) {
    (cDNA_frac / gDNA_frac) * sqrt( (sd_cDNA_udp_frac(cDNA_frac) / cDNA_frac / sqrt(n_cDNA))^2 +
                                      (sd_gDNA_udp_frac(gDNA_frac) / gDNA_frac / sqrt(n_gDNA))^2 )
  }
  
  df = data.frame(x=seq(0.0005, 0.9995, by=0.0005))
  # ggplot(df, aes(x = x*100)) +
  #   scale_x_log10(breaks = c(0.1, 1, 10)) +
  #   stat_function(fun = function(x) sd_ratio_cDNA_gDNA(x/100, x/100, n_cDNA_rep, n_gDNA_rep), color = "blue") +
  #   theme_bw() + theme(panel.grid.major = element_line(size = 1, colour = "grey95")) +
  #   xlab("UDP read %") + ylab("SE of cDNA:gDNA ratio")
  
  power = function(effect, cDNA_frac, gDNA_frac, n_cDNA, n_gDNA) {
    effect_se = sd_ratio_cDNA_gDNA(cDNA_frac, gDNA_frac, n_cDNA, n_gDNA)
    Zscore = (effect - 1) / effect_se
    pvalue = 2 * pnorm(-abs(Zscore))
    1 - pvalue
  }
  p.power = ggplot(df, aes(x = x*100)) + 
    scale_x_log10(breaks = c(0.1, 1, 10)) +
    stat_function(fun = function(x) power(1.05, x/100, x/100, n_cDNA_rep, n_gDNA_rep), mapping = aes(color = "1.05x"), size=1) +
    stat_function(fun = function(x) power(1.1, x/100, x/100, n_cDNA_rep, n_gDNA_rep), mapping = aes(color = "1.1x"), size=1) +
    stat_function(fun = function(x) power(1.2, x/100, x/100, n_cDNA_rep, n_gDNA_rep), mapping = aes(color = "1.2x"), size=1) +
    stat_function(fun = function(x) power(1.5, x/100, x/100, n_cDNA_rep, n_gDNA_rep), mapping = aes(color = "1.5x"), size=1) +
    stat_function(fun = function(x) power(2, x/100, x/100, n_cDNA_rep, n_gDNA_rep), mapping = aes(color = "2.0x"), size=1) +
    theme_bw() + xlab("UDP read %") + ylab("Power") +
    scale_color_discrete(name = "Effect size") +
    ggtitle(label = powerPlotTitle,
            subtitle = sprintf("Based on %d cDNA and %d gDNA replicates", n_cDNA_rep, n_gDNA_rep)) +
    theme(panel.grid.major = element_line(size = 1, colour = "grey95"),
          plot.title = element_text(size=14, hjust=0.5, face="bold"),
          plot.subtitle = element_text(size=12, hjust=0.5))
  
  # We also make a plot that shows the optimal allocation of replicates to cDNA
  # vs. gDNA.
  replicateAllocationPlot = function(udpFraction, showTitle=T) {
    getPowerDF = function(nRep, udpFraction) {
      df1.1 = data.frame(x = seq(1:(nRep-1)), effectSize = 1.1, udpFraction = udpFraction)
      df1.1$power = power(df1.1$effectSize, df1.1$udpFraction, df1.1$udpFraction, df1.1$x, nRep - df1.1$x)
      df1.2 = data.frame(x = seq(1:(nRep-1)), effectSize = 1.2, udpFraction = udpFraction)
      df1.2$power = power(df1.2$effectSize, df1.2$udpFraction, df1.2$udpFraction, df1.2$x, nRep - df1.2$x)
      df = rbind(df1.1, df1.2)
      df$effectSize = sprintf("%.1fx", df$effectSize)
      df
    }
    df.fraction = data.frame(x=seq(0.01, 0.99, by=0.01))
    df6 = getPowerDF(6, udpFraction)
    df10 = getPowerDF(10, udpFraction)
    df15 = getPowerDF(15, udpFraction)
    df25 = getPowerDF(25, udpFraction)
    nReps = factor(as.character(1:25), levels = as.character(1:25))
    #p = ggplot(df, aes(x = x)) +
    p = ggplot(df.fraction, aes(x = x)) +
      geom_point(data = df6,  mapping = aes(x = x / 6, y = power, color = effectSize, shape = "6"), size = 2.5) +
      geom_point(data = df10, mapping = aes(x = x / 10, y = power, color = effectSize, shape = "10"), size = 2.5) +
      geom_point(data = df15, mapping = aes(x = x / 15, y = power, color = effectSize, shape = "15"), size = 2.5) +
      geom_point(data = df25, mapping = aes(x = x / 25, y = power, color = effectSize, shape = "25"), size = 2.5) +
      stat_function(fun = function(x) power(1.1, udpFraction, udpFraction, 6*x, 6*(1-x)),   mapping = aes(color = "1.1x"), size=0.8, alpha = 0.8) +
      stat_function(fun = function(x) power(1.1, udpFraction, udpFraction, 10*x, 10*(1-x)), mapping = aes(color = "1.1x"), size=0.8, alpha = 0.8) +
      stat_function(fun = function(x) power(1.1, udpFraction, udpFraction, 15*x, 15*(1-x)), mapping = aes(color = "1.1x"), size=0.8, alpha = 0.8) +
      stat_function(fun = function(x) power(1.1, udpFraction, udpFraction, 25*x, 25*(1-x)), mapping = aes(color = "1.1x"), size=0.8, alpha = 0.8) +
      stat_function(fun = function(x) power(1.2, udpFraction, udpFraction, 6*x, 6*(1-x)),   mapping = aes(color = "1.2x"), size=0.8, alpha = 0.8) +
      stat_function(fun = function(x) power(1.2, udpFraction, udpFraction, 10*x, 10*(1-x)), mapping = aes(color = "1.2x"), size=0.8, alpha = 0.8) +
      stat_function(fun = function(x) power(1.2, udpFraction, udpFraction, 15*x, 15*(1-x)), mapping = aes(color = "1.2x"), size=0.8, alpha = 0.8) +
      stat_function(fun = function(x) power(1.2, udpFraction, udpFraction, 25*x, 25*(1-x)), mapping = aes(color = "1.2x"), size=0.8, alpha = 0.8) +
      theme_bw() + xlab("Fraction of replicates cDNA") + ylab("Power") +
      scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
      scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
      scale_color_discrete(name = "Effect size") +
      scale_shape_manual(name = "N Replicates", values = c("6"=16, "10"=17, "15"=15, "25"=8), breaks=c("6", "10", "15", "25")) +
      theme(panel.grid.major = element_line(size = 1, colour = "grey95"),
            plot.title = element_text(size=14, hjust=0.5, face="bold"),
            plot.subtitle = element_text(size=12, hjust=0.5))
    subtitleStr = sprintf("UDP fraction %d%%", as.integer(udpFraction*100))
    if (showTitle) {
      p = p + ggtitle(label = powerPlotTitle, subtitle = subtitleStr)
    } else {
      p = p + ggtitle(label = NULL, subtitle = subtitleStr)
    }
    p
  }
  # replicateAllocationPlot = function(showTitle=T) {
  #   getPowerDF = function(nRep, udpFraction) {
  #     df1.1 = data.frame(x = seq(1:(nRep-1)), effectSize = 1.1, udpFraction = udpFraction)
  #     df1.1$power = power(df1.1$effectSize, df1.1$udpFraction, df1.1$udpFraction, df1.1$x, nRep - df1.1$x)
  #     df1.2 = data.frame(x = seq(1:(nRep-1)), effectSize = 1.2, udpFraction = udpFraction)
  #     df1.2$power = power(df1.2$effectSize, df1.2$udpFraction, df1.2$udpFraction, df1.2$x, nRep - df1.2$x)
  #     df = rbind(df1.1, df1.2)
  #     #df$effectSize = sprintf("%.1fx", df$effectSize)
  #     df$effectSize_fraction = sprintf("%.1fx, %g%%", df$effectSize, df$udpFraction*100)
  #     df
  #   }
  #   df.fraction = data.frame(x=seq(0.01, 0.99, by=0.01))
  #   df6 = rbind(getPowerDF(6, 0.01), getPowerDF(6, 0.1))
  #   df10 = rbind(getPowerDF(10, 0.01), getPowerDF(10, 0.1))
  #   df15 = rbind(getPowerDF(15, 0.01), getPowerDF(15, 0.1))
  #   df25 = rbind(getPowerDF(25, 0.01), getPowerDF(25, 0.1))
  #   p = ggplot(df.fraction, aes(x = x, alpha = "val")) +
  #     geom_point(data = df6,  mapping = aes(x = x / 6, y = power, color = effectSize_fraction, shape = "6"), size = 2.5) +
  #     geom_point(data = df10, mapping = aes(x = x / 10, y = power, color = effectSize_fraction, shape = "10"), size = 2.5) +
  #     geom_point(data = df15, mapping = aes(x = x / 15, y = power, color = effectSize_fraction, shape = "15"), size = 2.5) +
  #     geom_point(data = df25, mapping = aes(x = x / 25, y = power, color = effectSize_fraction, shape = "25"), size = 2.5) +
  #     stat_function(fun = function(x) power(1.1, 0.01, 0.01, 6*x, 6*(1-x)),   mapping = aes(color = "1.1x, 1%"), size=0.8) +
  #     stat_function(fun = function(x) power(1.1, 0.01, 0.01, 10*x, 10*(1-x)), mapping = aes(color = "1.1x, 1%"), size=0.8) +
  #     stat_function(fun = function(x) power(1.1, 0.01, 0.01, 15*x, 15*(1-x)), mapping = aes(color = "1.1x, 1%"), size=0.8) +
  #     stat_function(fun = function(x) power(1.1, 0.01, 0.01, 25*x, 25*(1-x)), mapping = aes(color = "1.1x, 1%"), size=0.8) +
  #     stat_function(fun = function(x) power(1.2, 0.01, 0.01, 6*x, 6*(1-x)),   mapping = aes(color = "1.2x, 1%"), size=0.8) +
  #     stat_function(fun = function(x) power(1.2, 0.01, 0.01, 10*x, 10*(1-x)), mapping = aes(color = "1.2x, 1%"), size=0.8) +
  #     stat_function(fun = function(x) power(1.2, 0.01, 0.01, 15*x, 15*(1-x)), mapping = aes(color = "1.2x, 1%"), size=0.8) +
  #     stat_function(fun = function(x) power(1.2, 0.01, 0.01, 25*x, 25*(1-x)), mapping = aes(color = "1.2x, 1%"), size=0.8) +
  #     stat_function(fun = function(x) power(1.1, 0.1, 0.1, 6*x, 6*(1-x)),   mapping = aes(color = "1.1x, 10%"), size=0.8) +
  #     stat_function(fun = function(x) power(1.1, 0.1, 0.1, 10*x, 10*(1-x)), mapping = aes(color = "1.1x, 10%"), size=0.8) +
  #     stat_function(fun = function(x) power(1.1, 0.1, 0.1, 15*x, 15*(1-x)), mapping = aes(color = "1.1x, 10%"), size=0.8) +
  #     stat_function(fun = function(x) power(1.1, 0.1, 0.1, 25*x, 25*(1-x)), mapping = aes(color = "1.1x, 10%"), size=0.8) +
  #     stat_function(fun = function(x) power(1.2, 0.1, 0.1, 6*x, 6*(1-x)),   mapping = aes(color = "1.2x, 10%"), size=0.8) +
  #     stat_function(fun = function(x) power(1.2, 0.1, 0.1, 10*x, 10*(1-x)), mapping = aes(color = "1.2x, 10%"), size=0.8) +
  #     stat_function(fun = function(x) power(1.2, 0.1, 0.1, 15*x, 15*(1-x)), mapping = aes(color = "1.2x, 10%"), size=0.8) +
  #     stat_function(fun = function(x) power(1.2, 0.1, 0.1, 25*x, 25*(1-x)), mapping = aes(color = "1.2x, 10%"), size=0.8) +
  #     theme_bw() + xlab("Fraction of replicates cDNA") + ylab("Power") +
  #     scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
  #     scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
  #     scale_color_manual(name = "Effect size", values = c("#FF0000", "#FF8888", "#0000FF", "#8888FF")) +
  #     scale_alpha_manual(name = "Effect size", values = c(0.7)) +
  #     scale_shape_manual(name = "N Replicates", values = c("6"=16, "10"=17, "15"=15, "25"=8), breaks=c("6", "10", "15", "25")) +
  #     theme(panel.grid.major = element_line(size = 1, colour = "grey95"),
  #           plot.title = element_text(size=14, hjust=0.5, face="bold"),
  #           plot.subtitle = element_text(size=12, hjust=0.5)) +
  #     ggtitle(label = powerPlotTitle, subtitle = sprintf("UDP fraction %d%%", as.integer(udpFraction*100)))
  #   p
  # }  
  
  p.replicate_allocation = egg::ggarrange(replicateAllocationPlot(0.01, showTitle=T), replicateAllocationPlot(0.1, showTitle=F),
                                          ncol=1, heights = c(1,1), draw = F)
  
  plot_list = list(cv_plots = p.cv_plots, power = p.power, replicate_allocation = p.replicate_allocation)
  return(plot_list)
}


textPlot = function(text, title = NA) {
  ggThemeBlank = theme_bw() + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.border = element_blank())
  p = ggplot(data.frame(x=1:10, y=1:10), aes(x,y)) +
    geom_blank() + ggThemeBlank +
    annotate("text", x=5, y=5, label = text)
  if (!is.na(title)) {
    p = p + ggtitle(title)
  }
  return(p)
}

getReadUDP = function(seq, profile_chars = character(0)) {
  seq_chars = strsplit(seq, "")[[1]]
  getReadCharsUDP(strsplit(seq, "")[[1]], profile_chars)
}

getReadCharsUDP = function(seq_chars, profile_chars = character(0)) {
  if (length(profile_chars) > 0) {
    if (length(profile_chars) != length(seq_chars)) {
      stop(sprintf("Length of sequence to get UDP (%d chars) should be the same as the HDR profile length (%d chars).", length(seq_chars), length(profile_chars)))
    }
    seq_chars[seq_chars != '*' & profile_chars == "-"] = "-"
  } else {
    seq_chars[seq_chars != '*'] = "-"
  }
  return(str_c(seq_chars, collapse = ""))
}

# When there are "ambiguous" DNA letters in a sequence, we can make a set of
# sequences that represent all possible values. E.g. H indicates A or C or T (but not G).
# If two positions have ambiguous letters, then the combinations multiply.
# To handle this case, I use recursion below, without checking whether the combinations
# may become too large.
getHDRProfiles = function(hdr_profile) {
  hdr_profile_set = c()
  if (grepl("B|D|H|V", hdr_profile)) {
    hits = str_locate_all(hdr_profile, "B|D|H|V")[[1]]
    for (i in 1:nrow(hits)) {
      nt = str_sub(hdr_profile, hits[i,"start"], hits[i,"start"])
      hdr1 <- hdr2 <- hdr3 <- hdr_profile
      if (nt == "B") {
        str_sub(hdr1, hits[i,"start"], hits[i,"start"]) = "C"
        str_sub(hdr2, hits[i,"start"], hits[i,"start"]) = "G"
        str_sub(hdr3, hits[i,"start"], hits[i,"start"]) = "T"
      } else if (nt == "D") {
        str_sub(hdr1, hits[i,"start"], hits[i,"start"]) = "A"
        str_sub(hdr2, hits[i,"start"], hits[i,"start"]) = "G"
        str_sub(hdr3, hits[i,"start"], hits[i,"start"]) = "T"
      } else if (nt == "H") {
        str_sub(hdr1, hits[i,"start"], hits[i,"start"]) = "A"
        str_sub(hdr2, hits[i,"start"], hits[i,"start"]) = "C"
        str_sub(hdr3, hits[i,"start"], hits[i,"start"]) = "T"
      } else if (nt == "V") {
        str_sub(hdr1, hits[i,"start"], hits[i,"start"]) = "A"
        str_sub(hdr2, hits[i,"start"], hits[i,"start"]) = "C"
        str_sub(hdr3, hits[i,"start"], hits[i,"start"]) = "G"
      }
      hdr_profile_set = c(getHDRProfiles(hdr1), getHDRProfiles(hdr2), getHDRProfiles(hdr3))
    }
  } else {
    hdr_profile_set = hdr_profile
  }
  hdr_profile_set
}

# This method was slower unless there are many non-empty "profile chars"
# getReadSiteProfile_OLD = function(seq, profile_chars) {
#   seq_chars = strsplit(seq, "")[[1]]
#   if (length(profile_chars) != length(seq_chars)) {
#     stop(sprintf("Length of sequence to get read profile (%d chars) should be the same as the profile length (%d chars).", length(seq_chars), length(profile_chars)))
#   }
#   out_chars = seq_chars[isDNALetter(profile_chars)]
#   return(paste0(out_chars, collapse = ""))
# }

# profile_positions should have the indices of DNA letters to get in the
# input sequence
getReadSiteProfile = function(seq, profile_positions) {
  out_str = strrep(" ", times = length(profile_positions))
  for (i in 1:length(profile_positions)) {
    substr(out_str, i, i) = substr(seq, profile_positions[i], profile_positions[i])
  }
  return(out_str)
}

# getReadMismatchProfile = function(seq, ref_sequence) {
#   #print(sprintf("Seq:%s\t%s", seq, ref_sequence))
#   seq_chars = strsplit(seq, "")[[1]]
#   ref_chars = strsplit(ref_sequence, "")[[1]]
#   getReadCharsMismatchProfile(seq_chars, ref_chars)
# }
getReadCharsMismatchProfile = function(seq_chars, ref_chars) {
  output = rep("-", length(ref_chars))
  mismatchpos = (seq_chars != ref_chars)
  output[mismatchpos] = seq_chars[mismatchpos]
  return(str_c(output, collapse=""))
}

udpToBinary = function(udp) {
  udp_array = strsplit(udp,"")[[1]]
  binary_array = rep(0, length(udp_array))
  binary_array[udp_array == '*'] = 1
  return(binary_array)
}

strToChars = function(s) {
  strsplit(s,"")[[1]]
}

isDNALetter = function(s_chars) {
  return(s_chars %in% c('A', 'C', 'G', 'T'))
}

getMismatchCount = function(s, ref) {
  s_chars = strToChars(s)
  ref_chars = strToChars(ref)
  letterPos = isDNALetter(s_chars)
  sum(s_chars[letterPos] != ref_chars[letterPos])
}

getMismatchCharsCount = function(s_chars, ref_chars) {
  letterPos = isDNALetter(s_chars)
  sum(s_chars[letterPos] != ref_chars[letterPos])
}

getUDPRatioEstimate = function(udp.counts.df, replicates.df, numerator = "num_reads", denominator = "num_wt_reads", batchCol = NULL, randomEffectsCols = NULL) {
  num_reads = udp.counts.df[, numerator, drop = T]
  num_wt_reads = udp.counts.df[, denominator, drop = T]
  udp.counts.df$num_reads = num_reads
  udp.counts.df$num_wt_reads = num_wt_reads
  
  if (sum(udp.counts.df$num_reads) == 0) {
    # We can't run a model when all counts are zeroes
    return(tibble(effect = 0, effect_sd = NA, effect_confint_lo = NA, effect_confint_hi = NA, df_estimated = NA, pval = NA, method = NA))
  }
  
  numWtZero = sum(udp.counts.df$num_wt_reads == 0)
  if (numWtZero > 0) {
    udp.counts.df = udp.counts.df %>% filter(num_wt_reads != 0)
    warning(sprintf("In getUDPRatioEstimate: removing %d replicates with zero WT reads.", numWtZero))
  }
  if (is.null(randomEffectsCols) & is.null(batchCol)) {
    getUDPRatioEstimateWelch(cDNA_UDP_counts = udp.counts.df$num_reads[udp.counts.df$type == "cDNA"],
                             cDNA_wt_counts = udp.counts.df$num_wt_reads[udp.counts.df$type == "cDNA"],
                             gDNA_UDP_counts = udp.counts.df$num_reads[udp.counts.df$type == "gDNA"],
                             gDNA_wt_counts = udp.counts.df$num_wt_reads[udp.counts.df$type == "gDNA"])
  } else {
    udp.counts.df = udp.counts.df %>%
      dplyr::left_join(replicates.df %>% select(replicate, type, one_of(batchCol, randomEffectsCols)),
                       by = c("type", "replicate"))
    getUDPRatioEstimateLMER(udp.counts.df, batchCol = batchCol, randomEffectsCols = randomEffectsCols)
  }
}

getUDPRatioEstimateLMER = function(udp.counts.df, batchCol = NULL, randomEffectsCols = NULL) {
  formulaStr = "ratio_cmp ~ type"
  if (!is.null(batchCol)) {
    formulaStr = sprintf("%s + (1|%s)", formulaStr, batchCol)
    if (!is.null(randomEffectsCols)) {
      randomEffectsCols = setdiff(randomEffectsCols, batchCol)
    }
  }
  if (!is.null(randomEffectsCols)) {
    for (col in randomEffectsCols) {
      formulaStr = sprintf("%s + (1|%s)", formulaStr, col)
    }
  }
  
  udp.counts.df$ratio_cmp = udp.counts.df$num_reads / udp.counts.df$num_wt_reads
  if (class(udp.counts.df$ratio_cmp) != "numeric") {
    warning("WTF is going on?")
    str1 = sprintf("Class = %s", class(udp.counts.df$num_reads))
    str2 = sprintf("Class = %s", class(udp.counts.df$num_wt_reads))
    str3 = sprintf("Class = %s", class(udp.counts.df$ratio_cmp))
  }
  if (typeof(udp.counts.df$ratio_cmp) != "double") {
    warning("WTF is going on?")
    str1 = sprintf("Type = %s", typeof(udp.counts.df$num_reads))
    str2 = sprintf("Type = %s", typeof(udp.counts.df$num_wt_reads))
    str3 = sprintf("Type = %s", typeof(udp.counts.df$ratio_cmp))
  }
  udp.counts.df$type = factor(udp.counts.df$type, levels = c("gDNA", "cDNA"))
  
  # Run the model 
  lmer_res = lmerTest::lmer(as.formula(formulaStr), data = udp.counts.df)
  lmer_sm = summary(lmer_res)
  
  # Extract effect size and SE from the lmer
  meanRatio  = lmer_sm$coefficients["(Intercept)", "Estimate"]
  meanSE     = lmer_sm$coefficients["(Intercept)", "Std. Error"]
  cDNAEffect = lmer_sm$coefficients["typecDNA", "Estimate"]
  effectSE   = lmer_sm$coefficients["typecDNA", "Std. Error"]
  UNS = (meanRatio + cDNAEffect) / meanRatio
  
  # We have estimates and SEs from the model. I don't have a good way to determine
  # the SE of the ratio (cDNAEffect + meanRatio) / meanRatio. In general, the uncertainty
  # of meanRatio may be higher than that of the cDNAEffect -- but in the ratio this is
  # compensated for (e.g. if meanRatio is higher, then the numerator is also higher),
  # so that the propagation of uncertainty formula won't work, since we don't know the
  # covariance between the numerator and the denominator.
  # Therefore, I just approximate the SE of the ratio as SE(cDNAEffect) / meanRatio.
  UNS_SE = effectSE / meanRatio
  
  conf = confint(lmer_res, method = "Wald")
  cDNAEffectLo = conf["typecDNA", "2.5 %"]
  cDNAEffectHi = conf["typecDNA", "97.5 %"]
  UNS_confint_lo = (cDNAEffectLo + meanRatio) / meanRatio
  UNS_confint_hi = (cDNAEffectHi + meanRatio) / meanRatio
  
  methodStr = sprintf("LMM: %s", str_replace(formulaStr, "ratio_cmp ", ""))
  return(tibble(effect = UNS,
                effect_sd = UNS_SE,
                effect_confint_lo = UNS_confint_lo,
                effect_confint_hi = UNS_confint_hi,
                df_estimated = lmer_sm$coefficients["typecDNA", "df"],
                pval = lmer_sm$coefficients["typecDNA", "Pr(>|t|)"],
                method = methodStr))
}


getUDPRatioEstimateWelch = function(cDNA_UDP_counts, cDNA_wt_counts, gDNA_UDP_counts, gDNA_wt_counts) {
  N_cDNA = length(cDNA_wt_counts)
  N_gDNA = length(gDNA_wt_counts)
  
  # There are two ways that we could calculate the UDP:WT ratio (for both cDNA and gDNA).
  # One is to compute it for each replicate and take the average. The other is to sum
  # UDP counts for replicate, and divide by the summed WT counts for replicate. This latter
  # method effectively weights the ratio by the number of reads in each replicate. The main
  # problem with this method is that it gives a value that will differ from the one we
  # are assuming in our statistical test further down (t test). For example, our confidence
  # interval may include 1 (no effect), but our p value gives a significant result. For this
  # reason I use the unweighted UDP:WT ratios. It has the downside of not using all of the
  # reads across the replicates "efficiently", but if outliers are excluded properly, then
  # it should be a more robust measure of the true effect.
  cDNA_ratio = mean(cDNA_UDP_counts / cDNA_wt_counts)
  gDNA_ratio = mean(gDNA_UDP_counts / gDNA_wt_counts)
  sd_cDNA_ratio = sd(cDNA_UDP_counts / cDNA_wt_counts)
  sd_gDNA_ratio = sd(gDNA_UDP_counts / gDNA_wt_counts)
  SE_cDNA_ratio = sd_cDNA_ratio / sqrt(N_cDNA)
  SE_gDNA_ratio = sd_gDNA_ratio / sqrt(N_gDNA)
  
  effect = cDNA_ratio / gDNA_ratio
  effect_se = ratioUncertaintySD(cDNA_ratio, SE_cDNA_ratio, gDNA_ratio, SE_gDNA_ratio, 0)
  
  # Use a Welch two-sample T test to see whether the UDP:WT ratio differs
  # in cDNA vs. gDNA. This test accounts for there being different variance
  # in cDNA and gDNA.
  pvalue = NaN
  tryres = tryCatch({
    res = t.test(cDNA_UDP_counts / cDNA_wt_counts,
                 gDNA_UDP_counts / gDNA_wt_counts,
                 alternative = "two.sided", var.equal = F)
    pvalue = res$p.value
  }, error = function(e) e, warning = function(w) w)
  if (is(tryres, "error")) {
    warning("Error in t.test")
  }
  # Use a t distribution to determine the confidence interval. We estimate
  # the degrees of freedom using the Welch T test approximation. This gives
  # a wider confidence interval than if we assumed a normal distribution or
  # a standard T distribution.
  df_estimated = welch_df_estimate(sd_cDNA_ratio, N_cDNA, sd_gDNA_ratio, N_gDNA)
  tScore_alpha_0.05 = qt(0.975, df = df_estimated, lower.tail=T)
  effect_confint = tScore_alpha_0.05 * effect_se
  
  return(tibble(effect = effect,
                effect_sd = effect_se,
                effect_confint_lo = effect - tScore_alpha_0.05 * effect_se,
                effect_confint_hi = effect + tScore_alpha_0.05 * effect_se,
                df_estimated = df_estimated,
                pval = pvalue,
                method = "Welch t-test"))
}

getHDRRatioEstimate_WeightedByReads = function(stats.df, colToCompare = "num_hdr_reads") {
  cDNA.df = stats.df %>% dplyr::filter(type == "cDNA") %>% as.data.frame()
  gDNA.df = stats.df %>% dplyr::filter(type == "gDNA") %>% as.data.frame()
  N_cDNA = nrow(cDNA.df)
  N_gDNA = nrow(gDNA.df)
  
  cDNA_ratio = mean(cDNA.df[,colToCompare]) / mean(cDNA.df$num_wt_reads)
  gDNA_ratio = mean(gDNA.df[,colToCompare]) / mean(gDNA.df$num_wt_reads)
  # sd_cDNA_ratio = vecRatioUncertaintySD(cDNA.df[,colToCompare], cDNA.df$num_wt_reads)
  # sd_gDNA_ratio = vecRatioUncertaintySD(gDNA.df[,colToCompare], gDNA.df$num_wt_reads)
  sd_cDNA_ratio = sd(cDNA.df[,colToCompare] / cDNA.df$num_wt_reads)
  sd_gDNA_ratio = sd(gDNA.df[,colToCompare] / gDNA.df$num_wt_reads)
  SE_cDNA_ratio = sd_cDNA_ratio / sqrt(N_cDNA)
  SE_gDNA_ratio = sd_gDNA_ratio / sqrt(N_gDNA)
  
  effect = cDNA_ratio / gDNA_ratio
  effect_se = ratioUncertaintySD(cDNA_ratio, SE_cDNA_ratio, gDNA_ratio, SE_gDNA_ratio, 0)
  
  # I was unsure whether our estimate of this score represents a Z score (normal distribution)
  # or a T score. But assuming that it is a T score corresponds with the result you would
  # get from a linear model, and makes more sense given that we have a finite sample.
  #score = (effect - 1) / effect_se
  #pvalue = 2 * pnorm(-abs(score))
  #pvalue = 2*pt(score, df=38, lower.tail=F)
  
  # Use a Welch two-sample T test to see whether the UDP:WT ratio differs
  # in cDNA vs. gDNA. This test accounts for there being different variance
  # in cDNA and gDNA.
  res = t.test(cDNA.df[,colToCompare] / cDNA.df$num_wt_reads,
               gDNA.df[,colToCompare] / gDNA.df$num_wt_reads,
               alternative = "two.sided", var.equal = F)
  
  # Use a t distribution to determine the confidence interval. We estimate
  # the degrees of freedom using the Welch T test approximation. This gives
  # a wider confidence interval than if we assumed a normal distribution or
  # a standard T distribution.
  df_estimated = welch_df_estimate(sd_cDNA_ratio, N_cDNA, sd_gDNA_ratio, N_gDNA)
  tScore_alpha_0.05 = qt(0.975, df = df_estimated, lower.tail=T)
  effect_confint = tScore_alpha_0.05 * effect_se
  
  return(list(effect = effect, effect_sd = effect_se, effect_confint = effect_confint, pval = res$p.value))
}

welch_df_estimate = function(sdA, nA, sdB, nB) {
  round((sdA^2/nA + sdB^2/nB)^2 / ( (sdA^2/nA)^2 / (nA - 1) + (sdB^2/nB)^2 / (nB - 1) ))
}

getHDRRatioEstimate_lmer = function(stats.df, colToCompare = "num_hdr_reads") {
  # Test using a linear model to determine the effect of "cDNA", i.e.
  # whether it differs from gDNA. We first compute the UDP:WT ratio
  # for each sample (both cDNA and gDNA).
  stats.df$type = factor(stats.df$type, levels=c("gDNA", "cDNA"))
  stats.df$ratio_cmp = stats.df[,colToCompare] / stats.df$num_wt_reads
  
  #lm_res = lm(ratio_cmp ~ type, data = stats.df, weights = stats.df$num_reads)
  lm_res = lm(ratio_cmp ~ type, data = stats.df)
  lm_summary = summary(lm_res)
  #ggplot(stats.df, aes(x=type, y=ratio_cmp)) + geom_boxplot(outlier.shape=NA) + geom_jitter(width=0.25)
  
  # Get the values for parameters from the linear model result
  gDNA_ratio_est = lm_summary$coefficients["(Intercept)", "Estimate"]
  cDNA_ratio_est = gDNA_ratio_est + lm_summary$coefficients["typecDNA", "Estimate"]
  cDNA_gDNA_ratio_est = cDNA_ratio_est / gDNA_ratio_est
  se_gDNA_ratio = lm_summary$coefficients["(Intercept)", "Std. Error"]
  se_cDNA_ratio = lm_summary$coefficients["typecDNA", "Std. Error"]
  se_cDNA_gDNA_ratio = sqrt( cDNA_gDNA_ratio_est^2 * ((se_gDNA_ratio/gDNA_ratio_est)^2 + (se_cDNA_ratio/cDNA_ratio_est)^2) )

  # Propagation of uncertainty formula for calculating the SE of the
  # cDNA:gDNA ratio
  se_cDNA_gDNA = sqrt( (sd(stats.df %>% dplyr::filter(type == "gDNA") %>% .$ratio_cmp)^2 / sum(stats.df$type == "gDNA")) + (sd(stats.df %>% dplyr::filter(type == "cDNA") %>% .$ratio_cmp)^2 / sum(stats.df$type == "cDNA")) )

  
  score = (cDNA_ratio_est - gDNA_ratio_est) / se_cDNA_ratio
  se_cDNA_gDNA_ratio = cDNA_gDNA_ratio_est * (se_cDNA_ratio/cDNA_ratio_est)
  score2 = (cDNA_gDNA_ratio_est - 1) / se_cDNA_gDNA_ratio
  
  pvalue = 2 * pnorm(-abs(score)) # assume normal distribution
  pchisq(score^2, df=38, lower.tail=F) # 
  2*pt(score, df=38, lower.tail=F) # assume T distribution

  # A formula I found earlier for determining the "effective" degrees of freedom
  # for a T test with unequal variances. I think this is basically the same as
  # Welch's T test.
  df = round((se_cDNA_ratio + se_gDNA_ratio)^2 / ( ( (sd_cDNA_ratio^2/nrow(cDNA.df))^2 / (nrow(cDNA.df) - 1) ) + ((sd_gDNA_ratio^2/nrow(gDNA.df))^2 / (nrow(gDNA.df) - 1) ) ))
  
  # Use a mixed model to account for batch effects, e.g. of cDNA generation,
  # barcoding, sequencing.
  stats.df2 = stats.df %>% dplyr::left_join(
    replicates.df %>% dplyr::select(replicate, replicate_PCR, replicate_cDNA, replicate_barcode, replicate_sequencing),
    by="replicate")
  lmer_res = lmer(ratio_cmp ~ type + (1|replicate_cDNA), data = stats.df2)
  drop1(update(lmer_res, REML=F), test="Chisq")
  lmer_res = lmer(ratio_cmp ~ type + (1|replicate_cDNA) + (1|replicate_barcode), data = stats.df2)
  drop1(update(lmer_res, REML=F), test="Chisq")
  lmer_res = lmer(ratio_cmp ~ type + (1|replicate_cDNA) + (1|replicate_sequencing), data = stats.df2)
  drop1(update(lmer_res, REML=F), test="Chisq")
  lmer_res = lmer(ratio_cmp ~ type + (1|replicate_cDNA) + (1|replicate_PCR), data = stats.df2)
  drop1(update(lmer_res, REML=F), test="Chisq")
  lmer_res = lmer(ratio_cmp ~ type + (1|replicate_cDNA) + (1|replicate_barcode), data = stats.df2)
  drop1(update(lmer_res, REML=F), test="Chisq")
  lmer_res = lmer(ratio_cmp ~ type + (1|replicate_cDNA) + (1|replicate_barcode) + (1|replicate_sequencing), data = stats.df2)
  drop1(update(lmer_res, REML=F), test="Chisq")
  lmer_res = lmer(ratio_cmp ~ type + (1|replicate_cDNA) + (1|replicate_barcode) + (1|replicate_sequencing) + (1|replicate_PCR), data = stats.df2)
  drop1(update(lmer_res, REML=F), test="Chisq")
  
  lmer_res = lmer(ratio_cmp ~ type + (1|replicate_cDNA/replicate_PCR/replicate_barcode), data = stats.df2)
  drop1(update(lmer_res, REML=F), test="Chisq")

  cDNA.df = stats.df %>% dplyr::filter(type == "cDNA") %>% as.data.frame()
  gDNA.df = stats.df %>% dplyr::filter(type == "gDNA") %>% as.data.frame()
  
  cDNA_ratio = mean(cDNA.df[,colToCompare]) / mean(cDNA.df$num_wt_reads)
  effect = mean(cDNA.df[,colToCompare] / cDNA.df$num_wt_reads) - mean(gDNA.df[,colToCompare] / gDNA.df$num_wt_reads)
  gDNA_ratio = mean(gDNA.df[,colToCompare]) / mean(gDNA.df$num_wt_reads)
  cDNA_ratio-gDNA_ratio
  sd_cDNA_ratio = vecRatioUncertaintySD(cDNA.df[,colToCompare], cDNA.df$num_wt_reads)
  sd_gDNA_ratio = vecRatioUncertaintySD(gDNA.df[,colToCompare], gDNA.df$num_wt_reads)
  SE_cDNA_ratio = sd_cDNA_ratio / sqrt(nrow(cDNA.df))
  SE_gDNA_ratio = sd_gDNA_ratio / sqrt(nrow(gDNA.df))
  
  effect = cDNA_ratio / gDNA_ratio
  effect = abs(cDNA_ratio - gDNA_ratio)
  effect_se = ratioUncertaintySD(cDNA_ratio, SE_cDNA_ratio, gDNA_ratio, SE_gDNA_ratio, 0)
  effect_se = sqrt(SE_cDNA_ratio^2 + SE_gDNA_ratio^2)
  
  SE_cDNA_ratio = sd(cDNA.df[,colToCompare] / cDNA.df$num_wt_reads) / sqrt(nrow(cDNA.df))
  SE_gDNA_ratio = sd(gDNA.df[,colToCompare] / gDNA.df$num_wt_reads) / sqrt(nrow(gDNA.df))
  
  
  df = round((sd_cDNA_ratio^2/nrow(cDNA.df) + sd_gDNA_ratio^2/nrow(gDNA.df))^2 / ( ( (sd_cDNA_ratio^2/nrow(cDNA.df))^2 / (nrow(cDNA.df) - 1) ) + ((sd_gDNA_ratio^2/nrow(gDNA.df))^2 / (nrow(gDNA.df) - 1) ) ))
  Zscore = (effect) / effect_se
  pchisq(Zscore^2, df=36, lower.tail=F)
  2*pt(Zscore, df=38, lower.tail=F)
  pvalue = 2 * pnorm(-abs(Zscore))
  return(list(effect = effect, effect_sd = effect_se, pval = pvalue))
}

# Standard deviation function that returns zero if there is just one value
sdz = function(vec, na.rm = T) {
  if (length(vec) <= 1) 0 else sd(vec, na.rm = na.rm)
}

vecRatioUncertaintySD = function(A, B) {
  ratioUncertaintySD(mean(A), sdz(A), mean(B), sdz(B), cov(A, B))
}

ratioUncertaintySD = function(meanA, sdA, meanB, sdB, covAB) {
  abs(meanA / meanB) * sqrt((sdA / meanA)^2 + (sdB / meanB)^2 - (2*covAB / (meanA*meanB)))
}

printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}

