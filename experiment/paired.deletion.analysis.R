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
suppressMessages(library(optparse))
suppressMessages(library(variancePartition))
suppressMessages(library(doParallel))
cl <- makeCluster(4)
registerDoParallel(cl)
#library(profvis)
options(stringsAsFactors = F)

opt = NULL
#setwd("/Users/jeremys/work/opentargets/experiment/transcribed/batch1_redo")
#setwd("/Users/jeremys/work/opentargets/experiment/transcribed/batch1")


main = function()
{
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  opt_list <- list(
    make_option(c("--regions"), type="character", default=NULL, metavar="FILE.tsv", help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--replicates"), type="character", default=NULL, metavar="FILE.tsv", help="Path to tab-separated input file listing replicates. Required."),
    make_option(c("--out"), type="character", default=NULL, metavar="outpath", help="Base path for output files. Required."),
    make_option(c("--minMapQ"), type="integer", default=0, metavar="INT", help="Discard reads with mapping quality lower than minMapQ [default %default]"),
    #make_option(c("--subsample"), type="numeric", default=NULL, help=""),
    make_option(c("--max_mismatch_frac"), type="numeric", default=0.05, metavar="FLOAT", help="[default %default]"),
    make_option(c("--viewing_window"), type="integer", default=50, metavar="INT", help="[default %default]"),
    make_option(c("--editing_window"), type="integer", default=10, metavar="INT", help="[default %default]"),
    make_option(c("--min_window_overlap"), type="integer", default=30, metavar="INT", help="[default %default]"),
    make_option(c("--exclude_multiple_deletions"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--exclude_nonspanning_reads"), type="logical", default=T, action="store_true", help="[default %default]"),
    make_option(c("--exclude_nonspanning_deletions"), type="logical", default=T, action="store_true", help="[default %default]"),
    make_option(c("--uns_plot_min_gDNA"), type="integer", default=10, metavar="INT", help="[default %default]"),
    make_option(c("--uns_plot_min_cDNA"), type="integer", default=0, metavar="INT", help="[default %default]"),
    make_option(c("--uns_plot_max_udps"), type="integer", default=40, metavar="INT", help="[default %default]"),
    make_option(c("--no_allele_profile"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--no_site_profile"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--no_udp_profile"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--no_uns_profile"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--no_stats"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--no_replicate_plots"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--variance_analysis"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--variance_analysis_min_count"), type="integer", default=100, action="store_true", help="[default %default]"),
    make_option(c("--power_analysis"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--read_data"), type="character", default=NULL)
  )
  
  parser = OptionParser(usage = "paired.deletion.analysis.R --regions file.tsv --replicates file.tsv --out output_path [options]",
                        option_list = opt_list)
  opt <<- parse_args(parser)
  
  opt <<- list(#input = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1/analysis/udp_replicates.tsv",
    regions = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1_redo/regions.tsv",
    replicates = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1_redo/replicates.tsv",
    out = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1_redo/analysis/batch1_redo.sarah.edit_dels.2bp_window",
    minMapQ = 0,
    subsample = NULL,
    max_mismatch_frac = 0.05,
    #viewing_window = 30, # window around the cut site for viewing
    viewing_window = 40,
    editing_window = 1, # window around the cut site for counting deletions as edits
    min_window_overlap = 30, # minimum number of read bases correctly aligned within the region of interest
    exclude_multiple_deletions = F,
    exclude_nonspanning_reads = T,
    exclude_nonspanning_deletions = T,
    uns_plot_min_gDNA = 10,
    uns_plot_min_cDNA = 0,
    uns_plot_max_udps = 40,
    no_allele_profile = F,
    no_site_profile = F,
    no_udp_profile = F,
    no_uns_profile = F,
    no_stats = F,
    no_replicate_plots = T,
    variance_analysis = T,
    power_analysis = T,
    variance_analysis_min_count = 100,
    read_data = "analysis/batch1_redo.sarah.read_data.rds")
  opt <<- list(
    regions = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1/regions.tsv",
    replicates = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1/replicates.tsv",
    out = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1/analysis/batch1.2bp_window.new",
    minMapQ = 0,
    subsample = NULL,
    max_mismatch_frac = 0.05,
    #viewing_window = 30, # window around the cut site for viewing
    viewing_window = 40,
    editing_window = 1, # window around the cut site for counting deletions as edits
    min_window_overlap = 30, # minimum number of read bases correctly aligned within the region of interest
    exclude_multiple_deletions = F,
    exclude_nonspanning_reads = T,
    exclude_nonspanning_deletions = T,
    uns_plot_min_gDNA = 10,
    uns_plot_min_cDNA = 0,
    uns_plot_max_udps = 40,
    no_allele_profile = F,
    no_site_profile = F,
    no_udp_profile = F,
    no_uns_profile = F,
    no_stats = F,
    no_replicate_plots = T,
    variance_analysis = F,
    power_analysis = T,
    variance_analysis_min_count = 100,
    read_data = "analysis/batch1.read_data.new.rds")
  cat("All options:\n")
  printList(opt)
  
  checkRequiredOpt = function(opt, parser, optName) {
    if (is.null(opt[[optName]])) {
      print_help(parser)
      cat(sprintf("ERROR: argument --%s is required\n\n", optName))
      stop()
    }
  }
  checkRequiredOpt(opt, parser, "regions")
  checkRequiredOpt(opt, parser, "replicates")
  checkRequiredOpt(opt, parser, "out")
  
  regions.df = readr::read_tsv(opt$regions, col_types="ccciiiiccc")
  #regions.df = regions.df %>% dplyr::filter((index %in% c("29")))
  region_names = unique(regions.df$name)
  
  replicates.df = readr::read_tsv(opt$replicates)
  replicates.df$replicate = as.character(replicates.df$replicate)
  # Check that each replicate has a distinct name; otherwise it will
  # be impossible to distinguish replicates in the output files
  for (region_name in region_names) {
    df = replicates.df %>% dplyr::filter(name == region_name)
    if (any(duplicated(df$replicate))) {
      stop("Region %s: not all replicates have distinct names. It will be impossible to distinguish results in output files that relate to specific replicates.")
    }
  }
  
  read_data = NULL
  opt$save_read_data = F
  if (!is.null(opt$read_data)) {
    opt$save_read_data = !file.exists(opt$read_data)
    if (file.exists(opt$read_data)) {
      read_data = readRDS(opt$read_data)
    }
  }
  
  all_region_plots = list()
  all_results = list()
  region_index = 1
  #profvis({
  for (region_name in region_names) {
    df = regions.df %>% dplyr::filter(name == region_name)
    if (nrow(df) > 1) {
      stop("Error: multiple regions with the same region name (%s). Each input line for the --regions file should have a unique name.")
    } else if (nrow(df) == 0) {
      stop("Internal error.")
    }
    cur_region = as.list(df[1,])
    cur_replicates.df = replicates.df %>% dplyr::filter(name == region_name)
    result = doFullDeletionAnalysis(region = cur_region, replicates.df = cur_replicates.df, read_data)
    result$region = region_name
    result$replicate.df = cur_replicates.df
    
    all_results[[region_index]] = result
    all_region_plots[[region_index]] = result$plot_list
    region_index = region_index + 1
  }
  #})
  
  # Save read_data object if the option is set
  if (opt$save_read_data) {
    read_data_list = unlist(lapply(all_results, FUN = function(x) x$read_data),
                            recursive = F)
    saveRDS(read_data_list, file=opt$read_data)
  }
  
  settings.df = data.frame(minMapQ = opt$minMapQ,
                           max_mismatch_frac = opt$max_mismatch_frac,
                           viewing_window = opt$viewing_window, # window around the cut site for viewing
                           editing_window = opt$editing_window, # window around the cut site for counting deletions as edits
                           min_window_overlap = opt$min_window_overlap, # minimum number of read bases correctly aligned within the region of interest
                           exclude_multiple_deletions = opt$exclude_multiple_deletions,
                           exclude_nonspanning_reads = opt$exclude_nonspanning_reads,
                           exclude_nonspanning_deletions = opt$exclude_nonspanning_deletions,
                           uns_plot_min_gDNA = opt$uns_plot_min_gDNA,
                           uns_plot_min_cDNA = opt$uns_plot_min_cDNA,
                           uns_plot_max_udps = opt$uns_plot_max_udps)
  settings.df = data.frame(setting=names(opt), value=as.character(opt)) %>%
    filter(setting %in% c("minMapQ", "max_mismatch_frac", "viewing_window", "editing_window", "min_window_overlap",
                          "exclude_multiple_deletions", "exclude_nonspanning_reads", "exclude_nonspanning_deletions",
                          "uns_plot_min_gDNA", "uns_plot_min_cDNA", "uns_plot_max_udps"))
  ggThemeBlank = theme_bw() + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.border = element_blank())
  myTableTheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.8)),
                                            colhead = list(fg_params=list(cex = 0.8)),
                                            rowhead = list(fg_params=list(cex = 0.8)))
  settingsPlot = ggplot(data.frame(x=1:10, y=1:10), aes(x,y)) + geom_blank() + ggThemeBlank +
    annotation_custom(tableGrob(settings.df, theme = myTableTheme), xmin=1, xmax=9, ymin=1, ymax=9) +
    ggtitle("CRISPR editing analysis settings")
  
  fname = sprintf("%s.plots.test.pdf", opt$out)
  pdf(file = fname, width = 8, height = 7)
  print(settingsPlot)
  print(all_region_plots)
  dev.off()
  
  # fname = sprintf("%s.summary_plots.pdf", opt$out)
  # pdf(file = fname, width = 8, height = 7)
  # print(settingsPlot)
  # for (region_plots in all_region_plots) {
  #   print(region_plots$stats)
  #   print(region_plots$merged_udp)
  #   print(region_plots$merged_del_profile)
  #   print(region_plots$merged_UNS)
  # }
  # dev.off()
  
  if (!opt$no_allele_profile) {
    df = bind_rows(lapply(all_results, function(res) res$wt_hdr.df))
    fname = sprintf("%s.mismatch_profile.tsv", opt$out)
    write.table(df, fname, quote=F, row.names=F, col.names=T, na="", sep="\t")
  }
  if (!opt$no_udp_profile) {
    df = bind_rows(lapply(all_results, function(res) res$replicate.udp.df))
    fname = sprintf("%s.replicate_udps.tsv", opt$out)
    write.table(df, fname, quote=F, row.names=F, col.names=T, na="", sep="\t")
    
    df = bind_rows(lapply(all_results, function(res) res$merged.udp.df))
    fname = sprintf("%s.merged_udps.tsv", opt$out)
    write.table(df, fname, quote=F, row.names=F, col.names=T, na="", na="", sep="\t")
  }
  if (!opt$no_uns_profile) {
    df = bind_rows(lapply(all_results, function(res) res$uns.df))
    fname = sprintf("%s.uns.tsv", opt$out)
    write.table(df, fname, quote=F, row.names=F, col.names=T, sep="\t")
  }
  if (!opt$no_site_profile) {
    df = bind_rows(lapply(all_results, function(res) res$site.profiles.df))
    fname = sprintf("%s.site_profiles.tsv", opt$out)
    write.table(df, fname, quote=F, row.names=F, col.names=T, na="", sep="\t")
  }
  if (!opt$no_stats) {
    df = bind_rows(lapply(all_results, function(res) res$stats.df))
    fname = sprintf("%s.replicate_stats.tsv", opt$out)
    write.table(df, fname, quote=F, row.names=F, col.names=T, na="", sep="\t")
    
    # Also write out a summary of stats per region, rather than per replicate
    #hdr_prop_error_pval = lapply(all_results, function(res) res$stats.df)
    getListItemField = function(l, fieldList) {
      for (field in fieldList) {
        l = l[[field]]
        if (is.null(l)) {
          break
        }
      }
      if (is.null(l)) {NA} else {l}
    }
    getBetaBinomialPval = function(res, field) {
      x = res$stats.summary[[field]]
      if (is.null(x)) {NA} else {summary(x)@Coef$`Pr(> |z|)`[2]}
    }
    stats.summary.df = data.frame(
      name = region_names,
      hdr.error_prop.pval      = sapply(all_results, FUN = function(res) getListItemField(res, c("stats.summary", "hdr_wt.error_prop.res", "pval"))),
      hdr.error_prop.effect    = sapply(all_results, FUN = function(res) getListItemField(res, c("stats.summary", "hdr_wt.error_prop.res", "effect"))),
      hdr.error_prop.effect_sd = sapply(all_results, FUN = function(res) getListItemField(res, c("stats.summary", "hdr_wt.error_prop.res", "effect_sd"))),
      del.error_prop.pval      = sapply(all_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.error_prop.res", "pval"))),
      del.error_prop.effect    = sapply(all_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.error_prop.res", "effect"))),
      del.error_prop.effect_sd = sapply(all_results, FUN = function(res) getListItemField(res, c("stats.summary", "del_wt.error_prop.res", "effect_sd"))),
      hdr.beta_bin.pval = sapply(all_results, FUN = getBetaBinomialPval, "hdr_wt.beta_binomial.res"),
      del.beta_bin.pval = sapply(all_results, FUN = getBetaBinomialPval, "del_wt.beta_binomial.res")
    )
    fname = sprintf("%s.region_stats.tsv", opt$out)
    write.table(stats.summary.df, fname, quote=F, row.names=F, col.names=T, sep="\t")
  }
  
  # if (opt$variance_analysis) {
  #   replicate.udp.df = bind_rows(lapply(all_results, function(res) res$replicate.udp.df))
  #   # Get all columns of the replicate dataframe which begin with "replicate"
  #   replicate_cols = colnames(replicates.df)[grepl("^replicate_", colnames(replicates.df))]
  # }
}


doFullDeletionAnalysis = function(region, replicates.df, read_data = NULL)
{
  if (nrow(replicates.df) < 1) {
    stop(sprintf("doFullDeletionAnalysis: No replicates specified for region %s.", region$name))
  }
  replicate_list = list()
  replicate_udps = list()
  replicate_plots = list()
  replicate_site_profiles = list()
  replicate_wt_hdr = list()
  # stats_list = list()
  locus_name = region$name
  sites = list(start = region$start,
               end = region$end,
               highlight_site = region$highlight_site,
               cut_site = region$cut_site)
  rel_sites = getRelativeSites(sites, opt$viewing_window)
  for (i in 1:nrow(replicates.df)) {
    result = doReplicateDeletionAnalysis(locus_name,
                                         replicates.df[i,]$replicate,
                                         replicates.df[i,]$type,
                                         replicates.df[i,]$bam,
                                         region$sequence_name,
                                         sites,
                                         region$hdr_allele_profile,
                                         region$wt_allele_profile,
                                         region$ref_sequence,
                                         read_data)
    result$name = locus_name
    result$replicate = replicates.df[i,]$replicate
    result$type = replicates.df[i,]$type
    replicate_list[[i]] = result
    
    #stats_list[[i]] = result$stats
    
    if (!opt$no_site_profile) {
      replicate_site_profiles[[i]] = result$sites.profile.df %>% dplyr::mutate(name = locus_name,
                                                                               replicate = replicates.df[i,]$replicate,
                                                                               type = replicates.df[i,]$type) %>%
        dplyr::select(name, replicate, type, everything())
    }
    if (!opt$no_allele_profile) {
      replicate_wt_hdr[[i]] = result$wt_hdr.df %>% dplyr::mutate(name = locus_name,
                                                                 replicate = replicates.df[i,]$replicate,
                                                                 type = replicates.df[i,]$type) %>%
        dplyr::select(name, replicate, type, everything()) %>%
        dplyr::arrange(-num_reads)
    }
    
    replicate_udps[[i]] = result$udp.df %>% dplyr::mutate(name = locus_name,
                                                          replicate = replicates.df[i,]$replicate,
                                                          type = replicates.df[i,]$type) %>%
      dplyr::select(name, replicate, type, everything())
  }
  
  if (!opt$no_allele_profile) {
    wt_hdr.df = bind_rows(replicate_wt_hdr)
  }
  
  if (!opt$no_site_profile) {
    # Make a table of the "site profiles" - combinations of alleles at sites of interest
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
  
  stats_res = getFullReplicateStats(replicate_list)
  
  # Merge together results from all replicates
  replicate.udp.df = bind_rows(replicate_udps)
  
  if (!opt$no_replicate_plots) {
    for (i in 1:nrow(replicates.df)) {
      cur_replicate = replicates.df[i,]$replicate
      cur_type = replicates.df[i,]$type
      plot_title = sprintf("%s replicate %s, %s", locus_name, cur_replicate, cur_type)
      udp.df = replicate.udp.df %>% dplyr::filter(name == locus_name, replicate == cur_replicate, type == cur_type)
      replicate_plots[[i]] = getUDPPlot(udp.df, plot_title, rel_sites)
    }
  }
  
  merged.udp.df = summarise(replicate.udp.df %>% group_by(name, type, udp),
                            num_reads = sum(num_reads),
                            is_hdr_allele = first(is_hdr_allele),
                            is_wt_allele = first(is_wt_allele),
                            has_any_deletion = first(has_any_deletion),
                            has_crispr_deletion = first(has_crispr_deletion),
                            deletion_start = first(deletion_start),
                            deletion_end = first(deletion_end),
                            deletion2_start = first(deletion2_start),
                            deletion2_end = first(deletion2_end),
                            avg_seq_length = mean(avg_seq_length),
                            avg_mismatch_count = mean(avg_mismatch_count)) %>%
    arrange(-num_reads) %>% ungroup()
  
  region_title = sprintf("%s merged replicates", region$name)
  p.udp_gDNA = getUDPPlot(merged.udp.df %>% dplyr::filter(type == "gDNA"), plot_title="gDNA", sites=rel_sites)
  p.udp_cDNA = getUDPPlot(merged.udp.df %>% dplyr::filter(type == "cDNA"), plot_title="cDNA", sites=rel_sites)
  # now add the title
  p.title = ggdraw() + draw_label(region_title, fontface='bold')
  p.cDNA_gDNA = plot_grid(p.udp_gDNA, p.udp_cDNA, nrow=1)
  p.merged_UDP = plot_grid(p.title, p.cDNA_gDNA, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
  
  p.merged_del_profile = getDeletionProfilePlot(replicate.udp.df,
                                                rel_sites,
                                                plot_title = "Merged replicates",
                                                show_average = T, show_replicates = F)
  p.merged_del_profile = p.merged_del_profile + theme(axis.title.x=element_blank())
  p.replicate_del_profile = getDeletionProfilePlot(replicate.udp.df,
                                                   rel_sites,
                                                   plot_title = "Individual replicates",
                                                   show_average = T, show_replicates = T)
  p.del_profile = egg::ggarrange(ggdraw() + draw_label(region$name, fontface='bold'),
                                 p.merged_del_profile, p.replicate_del_profile, 
                                 ncol=1, heights=c(0.1, 0.5, 0.5), draw = F)
  
  p.merged_UNS = getUNSPlot(replicate.udp.df, rel_sites, plot_title = region_title,
                            min_gDNA_count = opt$uns_plot_min_gDNA, min_cDNA_count = opt$uns_plot_min_cDNA,
                            max_udps = opt$uns_plot_max_udps)
  
  uns_res = getUNSData(replicate.udp.df, rel_sites, region_title, min_gDNA_count = opt$uns_plot_min_gDNA,
                       min_cDNA_count = opt$uns_plot_min_cDNA)
  uns.df = NULL
  uns.replicates.df = NULL
  if (!is.null(uns_res)) {
    uns.df = uns_res$udp.dels.df %>% dplyr::mutate(name = region$name) %>% dplyr::select(name, everything())
    uns.replicates.df = uns_res$replicate.dels.df %>% dplyr::mutate(name = region$name) %>% dplyr::select(name, everything())
  }
  
  p.variance_components = NULL
  if (opt$variance_analysis) {
    p.variance_components = getVarianceComponentsPlots(replicate.udp.df, replicates.df, method = "log_reads", min_udp_total_count = opt$variance_analysis_min_count)
  }
  
  p.variance_fit = NULL
  p.power = NULL
  if (opt$power_analysis) {
    power_plots = getPowerPlots(replicate.udp.df, replicates.df, titlestr = region$name, min_udp_total_count = opt$variance_analysis_min_count)
    p.variance_fit = power_plots$cv_plots
    p.power = power_plots$power
  }
  
  plot_list = list(stats = stats_res$p.stats1,
                   merged_udp = p.merged_UDP,
                   merged_del_profile = p.del_profile,
                   merged_UNS = p.merged_UNS,
                   variance_components = p.variance_components,
                   variance_fit = p.variance_fit,
                   power = p.power,
                   replicate_plots = replicate_plots)
  
  read_data = lapply(replicate_list, FUN = function(x) x$read_data)
  names(read_data) = paste(replicates.df$name, replicates.df$replicate, replicates.df$type, sep="_")
  
  result_list = list(replicate_list = replicate_list,
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

getFullReplicateStats = function(replicate_list) {
  stats.df = bind_rows(lapply(replicate_list, function(res) as.data.frame(res$stats)))
  stats.df$HDR_WT_ratio = (stats.df$num_hdr_reads / stats.df$num_wt_reads)
  stats.df$HDR_rate = stats.df$num_hdr_reads / stats.df$num_reads
  stats.df$DEL_WT_ratio = (stats.df$num_deletion_reads / stats.df$num_wt_reads)
  stats.df$DEL_rate = stats.df$num_deletion_reads / stats.df$num_reads
  stats.df$editing_rate = stats.df$num_edit_reads / stats.df$num_reads
  
  # Summary of stats from the propagation of errors method
  hdr_ratio_res = NULL
  del_ratio_res = NULL
  if (sum(stats.df$type == "cDNA") < 2 | sum(stats.df$type == "gDNA") < 2) {
    summary.method.prop = sprintf("Propagation of uncertainty\nUnable to calculate with < 2 replicates")
  } else {
    hdr_ratio_res = getHDRRatioEstimate(stats.df, colToCompare = "num_hdr_reads")
    hdr.conf.interval.str = sprintf("95%% CI: (%.3g, %.3g)",
                                    hdr_ratio_res$effect - 1.96*hdr_ratio_res$effect_sd,
                                    hdr_ratio_res$effect + 1.96*hdr_ratio_res$effect_sd)
    hdr.summary.prop = sprintf("Propagation of uncertainty\ncDNA:gDNA ratio (HDR/WT): %.3g\n%s,    p = %.3g",
                               hdr_ratio_res$effect, hdr.conf.interval.str, hdr_ratio_res$pval)
    
    del_ratio_res = getHDRRatioEstimate(stats.df, colToCompare = "num_deletion_reads")
    del.conf.interval.str = sprintf("95%% CI: (%.3g, %.3g)",
                                    del_ratio_res$effect - 1.96*del_ratio_res$effect_sd,
                                    del_ratio_res$effect + 1.96*del_ratio_res$effect_sd)
    del.summary.prop = sprintf("cDNA:gDNA ratio (DEL/WT): %.3g\n%s,    p = %.3g",
                               del_ratio_res$effect, del.conf.interval.str, del_ratio_res$pval)
    summary.method.prop = paste(hdr.summary.prop, del.summary.prop, sep = "\n")
  }
  
  # Summary of stats from the beta binomial model method
  fm1 = NULL
  fm2 = NULL
  if (sum(stats.df$type == "cDNA") < 1 | sum(stats.df$type == "gDNA") < 1) {
    summary.method.betabin = sprintf("Beta binomial\nUnable to calculate\nOnly %d gDNA and %d cDNA replicates",
                                     sum(stats.df$type == "gDNA"), sum(stats.df$type == "cDNA"))
  } else {
    test.df.hdr = stats.df %>%
      dplyr::mutate(name = paste(name, replicate, type, sep="_"),
                    is_cDNA = as.integer(type == "cDNA")) %>%
      dplyr::select(name, is_cDNA, wt = num_wt_reads, hdr = num_hdr_reads)
    #fm1 = betabin(cbind(hdr, wt) ~ is_cDNA, ~ 1, data = test.df[,2:4])
    #test.df.hdr$is_cDNA[3] = T
    betabin_warn_err = ""
    res_hdr = tryCatch(fm1 <- betabin(cbind(hdr, wt) ~ is_cDNA, ~ 1, data = test.df.hdr[2:3,2:4]),
                       error = function(e) e, warning = function(w) w)
    betabin_hdr_star = ""
    if (is(res_hdr, "warning") | is(res_hdr, "error")) {
      betabin_hdr_star = "**"
      betabin_warn_err = paste(strwrap(paste("**", trimws(res_hdr$message)), width = 50), collapse = "\n")
      if (is(res_hdr, "warning")) {
        # When it's a warning it still produces output, but for some reason
        # when we catch the exception the return value fm1 isn't assigned
        fm1 <- betabin(cbind(hdr, wt) ~ is_cDNA, ~ 1, data = test.df.hdr[2:3,2:4])
      }
    }
    betabin_hdr_pval_str = "**"
    if (!is.null(fm1)) {
      betabin_hdr_pval_str = sprintf("%.3g%s", summary(fm1)@Coef$`Pr(> |z|)`[2], betabin_hdr_star)
    }
    
    test.df.del = stats.df %>%
      dplyr::mutate(name = paste(name, replicate, type, sep="_"),
                    is_cDNA = as.integer(type == "cDNA")) %>%
      dplyr::select(name, is_cDNA, wt = num_wt_reads, del = num_deletion_reads)
    #fm2 = betabin(cbind(del, wt) ~ is_cDNA, ~ 1, data = test.df.del[,2:4])
    res_del = tryCatch(fm2 <- betabin(cbind(del, wt) ~ is_cDNA, ~ 1, data = test.df.del[,2:4]),
                       error = function(e) e, warning = function(w) w)
    betabin_del_star = ""
    if (is(res_del, "warning") | is(res_del, "error")) {
      betabin_del_star = "**"
      betabin_warn_err = paste("**", trimws(res_del$message))
      if (is(res_del, "warning")) {
        # When it's a warning it still produces output, but for some reason
        # when we catch the exception the return value fm1 isn't assigned
        fm2 <- betabin(cbind(del, wt) ~ is_cDNA, ~ 1, data = test.df.del[,2:4])
      }
    }
    betabin_del_pval_str = "**"
    if (!is.null(fm2)) {
      betabin_del_pval_str = sprintf("%.3g%s", summary(fm2)@Coef$`Pr(> |z|)`[2], betabin_del_star)
    }
    
    summary.method.betabin = sprintf("Beta binomial\ncDNA:gDNA ratio (HDR/WT): p = %s\ncDNA:gDNA ratio (DEL/WT): p = %s\n%s",
                                     betabin_hdr_pval_str, betabin_del_pval_str, betabin_warn_err)
  }
  
  # Convert to strings for nice printing (with 3 significant digits)
  stats.plot.df = stats.df %>% dplyr::select(replicate, type, num_udps, HDR_WT_ratio, DEL_WT_ratio,
                                             HDR_rate, DEL_rate, editing_rate, num_reads, num_hdr_reads, num_wt_reads,
                                             num_deletion_reads, reads_excluded_for_minoverlap, reads_excluded_for_mismatches,
                                             reads_excluded_nonspanning, reads_excluded_for_multiple_deletions)
  stats.plot.df = stats.plot.df %>% dplyr::rename("HDR reads" = num_hdr_reads,
                                                  "WT reads" = num_wt_reads,
                                                  "Deletion reads" = num_deletion_reads,
                                                  "excluded-minoverlap" = reads_excluded_for_minoverlap,
                                                  "excluded-mismatches" = reads_excluded_for_mismatches,
                                                  "excluded-nonspanning" = reads_excluded_nonspanning,
                                                  "excluded-mult.deletions" = reads_excluded_for_multiple_deletions)
  stats.plot.df$HDR_WT_ratio = sprintf("%.3g", stats.plot.df$HDR_WT_ratio)
  stats.plot.df$DEL_WT_ratio = sprintf("%.3g", stats.plot.df$DEL_WT_ratio)
  stats.plot.df$HDR_rate = sprintf("%.2f%%", 100 * stats.plot.df$HDR_rate)
  stats.plot.df$DEL_rate = sprintf("%.2f%%", 100 * stats.plot.df$DEL_rate)
  stats.plot.df$editing_rate = sprintf("%.2f%%", 100 * stats.plot.df$editing_rate)
  ggThemeBlank = theme_bw() + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), panel.border = element_blank())
  
  mycex = 0.8
  if (ncol(stats.plot.df) > 12) { mycex = 0.6 }
  myTableTheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = mycex)),
                                            colhead = list(fg_params=list(cex = mycex)),
                                            rowhead = list(fg_params=list(cex = mycex)),
                                            padding = unit(c(2, 4), "mm"))
  p.stats1 = ggplot(data.frame(x=1:10, y=1:10), aes(x,y)) + geom_blank() + ggThemeBlank +
    annotation_custom(tableGrob(t(stats.plot.df), theme = myTableTheme), xmin=1, xmax=9, ymin=1, ymax=8) +
    annotate("text", x=5, y=10, label = sprintf("%s summary", stats.df$name[1]), vjust = 1, fontface = 2, size = 5) +
    annotate("text", x=1, y=9.3, label = summary.method.prop, vjust = 1, hjust = 0, size = 3.5) +
    annotate("text", x=9.5, y=9.3, label = summary.method.betabin, vjust = 1, hjust = 1, size = 3.5)
  if (!opt$no_stats) {
    p.stats1 = p.stats1 + annotate("text", x=1, y=1, label = "Data saved in *.replicate_stats.tsv", vjust = 1, hjust = 0, size = 3)
  }
  
  stats.summary = list(hdr_wt.beta_binomial.res = fm1,
                       hdr_wt.error_prop.res = hdr_ratio_res,
                       del_wt.beta_binomial.res = fm2,
                       del_wt.error_prop.res = del_ratio_res,
                       error_prop.summary = summary.method.prop,
                       beta_binomial.summary = summary.method.betabin)
  
  # mainStr = sprintf("Estimate ('is_cDNA') for whether HDR / WT reads\nacross replicates differs for cDNA vs. gDNA",
  #                   summary(fm1)@Coef$`Pr(> |z|)`[2])
  # summaryStr = paste(capture.output(summary(fm1)), collapse = "\n")
  # p.stats2 = ggplot(data.frame(x=c(0,1), y=c(0,1)), aes(x=x, y=y)) + theme_bw() + ggtitle(sprintf("Region %s stats", stats.df$name[1])) +
  #   coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  #   annotate("text", x=0, y=1, label = mainStr, hjust = 0, vjust = 1, fontface = 2) +
  #   annotate("text", x=0, y=0.85, label = summaryStr, hjust = 0, vjust = 1) +
  #   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
  #         axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())
  
  result_list = list(p.stats1 = p.stats1,
                     stats.df = stats.df,
                     stats.summary = stats.summary)
  return(result_list)
}

######################################################################################

doReplicateDeletionAnalysis = function(name, replicate, type, bam_file, sequence_name, sites, hdr_profile, wt_profile, ref_sequence, read_data_all = NULL)
{
  cat(sprintf("\n\nAnalysing region %s, replicate %s, file %s, %s:%d-%d, highlight site %d, cut site %d\n",
              name, replicate, bam_file, sequence_name, sites$start, sites$end, sites$highlight_site, sites$cut_site))
  cat(sprintf("HDR allele: %s\n", hdr_profile))
  cat(sprintf("WT allele:  %s\n", wt_profile))
  cat(sprintf("REF sequence: %s\n", ref_sequence))
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
    read_data = getRegionReadsFromBam(name, replicate, bam_file, sequence_name, sites$start, sites$end, ref_sequence)
    if (is.null(read_data)) {
      cat(sprintf("\nERROR: No reads retrieved at the specified region for replicate %s from file %s\n", replicate_name, bam_file))
      stop()
    }
  } else {
    read_data = read_data_all[[replicate_name]]
    if (is.null(read_data)) {
      cat(sprintf("\nERROR: Failed to get read data for replicate %s from file %s\n", replicate_name, opt$read_data))
      stop()
    }
  }
  reads.df = read_data$alignedReads
  
  cat(sprintf("%d cDNA reads of %d total (%.1f%%) were soft-clipped\n", read_data$num_softclipped, read_data$num_reads, 100.0 * read_data$num_softclipped / read_data$num_reads))
  cat(sprintf("%d cDNA reads of %d total (%.1f%%) were hard-clipped\n", read_data$num_hardclipped, read_data$num_reads, 100.0 * read_data$num_hardclipped / read_data$num_reads))
  cat(sprintf("%d cDNA reads of %d total (%.1f%%) had insertions\n", read_data$num_insertion, read_data$num_reads, 100.0 * read_data$num_insertion / read_data$num_reads))
  stats$num_softclipped = read_data$num_softclipped
  stats$num_hardclipped = read_data$num_hardclipped
  stats$num_insertion = read_data$num_insertion
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
    stats$reads_excluded_nonspanning = sum(!reads.df$spanning_read)
    cat(sprintf("%d of %d reads (%.2f%%) excluded due to not spanning the site of interest (position %d)\n",
                stats$reads_excluded_nonspanning, stats$num_reads, 100.0 * stats$reads_excluded_nonspanning / stats$num_reads,
                span_site))
    reads.df = reads.df[reads.df$spanning_read, ]
  }
  
  # Exclude reads that don't cover enough of the region of interest
  exclude_for_overlap = (reads.df$seq_length < opt$min_window_overlap)
  stats$reads_excluded_for_minoverlap = sum(exclude_for_overlap)
  cat(sprintf("%d of %d reads (%.2f%%) excluded due to aligning to less than %d bp in the region of interest\n",
              stats$reads_excluded_for_minoverlap, stats$num_reads, 100.0 * stats$reads_excluded_for_minoverlap / stats$num_reads,
              opt$min_window_overlap))
  reads.df = reads.df[!exclude_for_overlap, ]
  
  # Identify unique deletion profile (UDP) for each read
  #reads.df$udp = getReadUDPs(reads.df$region_read, wt_profile_chars)
  reads.df$udp = sapply(reads.df$read_chars, FUN=getReadCharsUDP, wt_profile_chars)
  #reads.df$udp = sapply(reads.df$region_read, FUN=getReadUDP, wt_profile_chars) # SLOWER
  
  reads.df$has_any_deletion = sapply(reads.df$udp, FUN=function(s) grepl("[*]", s))
  rel_cut_site = sites$cut_site - sites$start + 1
  if (opt$exclude_nonspanning_deletions) {
    reads.df$has_crispr_deletion = grepl("[*]", substr(reads.df$udp, rel_cut_site - opt$editing_window, rel_cut_site + opt$editing_window))
    
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
    reads.df$udp[reads.df$has_any_deletion] = sapply(reads.df$udp[reads.df$has_any_deletion], updateUDPEdits, rel_cut_site, opt$editing_window)
  } else {
    # Accept deletions anywhere in the read
    reads.df$has_crispr_deletion = reads.df$has_any_deletion
  }
  # Include deletions anywhere in the read
  reads.df$has_multiple_deletions = reads.df$has_any_deletion
  reads.df$has_multiple_deletions[reads.df$has_multiple_deletions] = 
    sapply(reads.df$udp[reads.df$has_multiple_deletions], FUN=function(s) grepl("[*]+[^*]+[*]+", s))
  
  # Identify which reads are WT or HDR. Technically this only depends on the UDP
  # and not the entire read itself, but we compute it now so that we can accurately
  # count mismatches, without counting the introduced SNPs themselves as mismatches.
  reads.df$is_wt_allele = (reads.df$udp == wt_profile)
  if (any(reads.df$is_wt_allele & reads.df$has_crispr_deletion)) {
    stop("ERROR: unexpected - found read classified as WT allele but having deletion")
  }
  n_wt_vars = sum(wt_profile_chars != '-' & wt_profile_chars != ref_seq_chars)
  if (n_wt_vars > 0) {
    reads.df$mismatch_count[reads.df$is_wt_allele] = reads.df$mismatch_count[reads.df$is_wt_allele] - n_wt_vars
  }
  # make sure we only count as WT those reads which actually cover the HDR site
  reads.df$is_wt_allele[!reads.df$spanning_read] = NA
  
  reads.df$is_hdr_allele = F
  if (hdr_profile != "") {
    if (nchar(hdr_profile) != region_length) {
      stop(sprintf("The HDR allele profile given has length %d, but the region size (end - start + 1) is %d", nchar(hdr_profile), region_length))
    }
    # Check that the HDR profile and WT profile have DNA characters at
    # the same positions
    if (any((hdr_profile_chars == "-") != (wt_profile_chars == "-"))) {
      stop("Error: HDR profile and WT profile should both indicate the expected sequence letters at the same positions")
    }
    reads.df$is_hdr_allele = (reads.df$udp == hdr_profile)
    reads.df$is_hdr_allele = reads.df$is_hdr_allele & !reads.df$has_crispr_deletion
    # We don't want to count the HDR site itself as a mismatch
    n_hdr_vars = sum(hdr_profile_chars != '-' & hdr_profile_chars != ref_seq_chars)
    if (n_hdr_vars > 0) {
      reads.df$mismatch_count[reads.df$is_hdr_allele] = reads.df$mismatch_count[reads.df$is_hdr_allele] - n_hdr_vars
    }
  }
  if (any(reads.df$mismatch_count < 0)) {
    stop("ERROR: something went wrong - got a negative mismatch count.")
  }
  
  # Exclude reads with too many mismatches
  exclude_for_mismatches = (reads.df$mismatch_count / reads.df$seq_length) > opt$max_mismatch_frac
  stats$reads_excluded_for_mismatches = sum(exclude_for_mismatches)
  cat(sprintf("%d of %d reads (%.2f%%) excluded due to having more than %.2f%% of the read being mismatches\n",
              stats$reads_excluded_for_mismatches, stats$num_reads, 100.0 * stats$reads_excluded_for_mismatches / stats$num_reads,
              opt$max_mismatch_frac * 100))
  reads.df = reads.df[!exclude_for_mismatches, ]
  
  stats$num_wt_reads = sum(reads.df$is_wt_allele * reads.df$count, na.rm = T)
  stats$num_hdr_reads = sum(reads.df$is_hdr_allele * reads.df$count, na.rm = T)
  stats$num_deletion_reads = sum(reads.df$has_crispr_deletion * reads.df$count, na.rm = T)
  
  # Aggregate reads according to their UDP
  udp.df = summarise(reads.df %>% group_by(udp),
                     num_reads = sum(count),
                     is_hdr_allele = first(is_hdr_allele),
                     is_wt_allele = first(is_wt_allele),
                     spanning_read = first(spanning_read),
                     has_any_deletion = first(has_any_deletion),
                     has_crispr_deletion = first(has_crispr_deletion),
                     has_multiple_deletions = first(has_multiple_deletions),
                     avg_seq_length = mean(seq_length),
                     avg_mismatch_count = mean(mismatch_count)) %>%
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
  wt_hdr.df = NA
  if (!opt$no_allele_profile) {
    wt_hdr.reads.df = reads.df %>% dplyr::filter(is_wt_allele | is_hdr_allele)
    #wt_hdr.reads.df$mismatch_profile = getReadMismatchProfiles(wt_hdr.reads.df$region_read, ref_seq_chars)
    wt_hdr.reads.df$mismatch_profile = sapply(wt_hdr.reads.df$read_chars, FUN=getReadCharsMismatchProfile, ref_seq_chars)
    
    wt_hdr.df = summarise(wt_hdr.reads.df %>% group_by(mismatch_profile),
                          num_reads = sum(count),
                          mismatch_count = first(mismatch_count),
                          spanning_read = first(spanning_read),
                          is_wt_allele = first(is_wt_allele),
                          is_hdr_allele = first(is_hdr_allele)) %>%
      dplyr::arrange(-num_reads)
  }
  
  # Make a table which has all variations of read sequences at the "profile" sites,
  # i.e. those sites used to identify WT and HDR alleles
  sites.profile.df = NA
  if (!opt$no_site_profile) {
    #reads.df$sites_profile = getReadSiteProfiles(reads.df$region_read, wt_profile_chars)
    profile_positions = which(isDNALetter(wt_profile_chars))
    reads.df$sites_profile = sapply(reads.df$region_read, FUN=getReadSiteProfile, profile_positions)
    sites.profile.df = reads.df %>% group_by(sites_profile) %>%
      summarise(count = sum(count))
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
  read_data = NULL
  if (length(sam_reads) > 0) {
    read_data = getAlignedReads(sam_reads, ref_sequence, start, end)
  }
  return(read_data)
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
  read_contains_deletion = FALSE
  read_status = ""
  expanded_cigar = ""
  cigar_bits = ""
  cursor_read = 1
  
  for (i in 1:length(sam_cigar_parsed$letters)) {
    type = sam_cigar_parsed$letters[i]
    
    if (type == 'M') {
      cigar_bits[i] = substr(sam_read, cursor_read, cursor_read + sam_cigar_parsed$numbers[i] - 1)
      #expanded_cigar = paste0(expanded_cigar, substr(sam_read, cursor_read, cursor_read + sam_cigar_parsed$numbers[i] - 1))
      cursor_read = cursor_read + sam_cigar_parsed$numbers[i]
    }
    else if (type == 'D') {
      cigar_bits[i] = strrep('*', sam_cigar_parsed$numbers[i])
      #expanded_cigar = paste0(expanded_cigar, strrep('*', sam_cigar_parsed$numbers[i]))
      read_contains_deletion = T
    } else if (type == 'S') {
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
  #print registered_read
  return(list("read" = registered_read, "read_ok" = read_ok, "read_status" = read_status, "deletion" = read_contains_deletion))
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
  for (i in 1:nrow(keptReads.df)) {
    # computes the read sequence within the coords of interest
    output = registerRead(keptReads.df$relative_pos[i], keptReads.df$sam_cigar[i], keptReads.df$seq[i], ref_seq_length)
    if (output$read_ok) {
      alignedReads[i] = output$read
    }
    if (output$read_status == "hardclipped") {
      num_hardclipped = num_hardclipped + 1
    } else if (output$read_status == "insertion") {
      num_insertion = num_insertion + 1
    } else if (output$read_status == "softclipped") {
      num_softclipped = num_softclipped + 1
    }
  }
  keptReads.df$region_read = alignedReads
  keptReads.df = keptReads.df %>% dplyr::filter(region_read != "") %>%
    dplyr::select(region_read, count, sam_position, sam_cigar, seq)
  # Some reads may have a distinct sequence yet have the same registered
  # read sequence, because the soft-clipped sequence could differ.
  
  return(list("alignedReads" = keptReads.df, "discardedReads" = discardedReads.df, "num_reads" = num_reads,
              "num_hardclipped" = num_hardclipped, "num_insertion" = num_insertion,
              "num_softclipped" = num_softclipped, "num_outsidewindow" = num_outsidewindow))
}


getUDPPlot = function(udp.df, plot_title, sites) {
  plot_list = list()
  # We allow there to be up to 2 deletions in the UDP
  udp.df = udp.df %>% filter(has_crispr_deletion, !is.na(deletion_start))
  if (nrow(udp.df) == 0) {
    return(egg::ggarrange(textPlot("No UDPs to plot."), top=plot_title, draw = F))
  }
  udp.dels.df = udp.df %>% dplyr::select(deletion_start, deletion_end, deletion2_start, deletion2_end)
  udp.dels.df = udp.dels.df %>% arrange(deletion_start, deletion2_start)
  udp.dels.df$y = 1:nrow(udp.dels.df)
  udp.plot.df = bind_rows(udp.dels.df %>% dplyr::select(y, deletion_start, deletion_end),
                          udp.dels.df %>% dplyr::select(y, deletion_start=deletion2_start, deletion_end=deletion2_end))
  udp.plot.df = udp.plot.df %>% dplyr::filter(!is.na(deletion_start))
  
  xmax = nchar(udp.df$udp[1])
  segment_size = 0.5
  if (nrow(udp.dels.df) > 200) {
    segment_size = 0.3
  }
  simple_theme = theme_bw() + theme(panel.grid.minor = element_line(colour="grey90", size=0.1), panel.grid.major = element_line(colour="grey90", size=0.1))
  p.udp_dist = ggplot(udp.plot.df) +
    geom_segment(aes(x = deletion_start, xend = deletion_end, y = y, yend = y), size = segment_size, color = "dodgerblue3") +
    simple_theme + coord_cartesian(xlim=c(sites$start, sites$end)) + scale_y_reverse() +
    xlab("Nucleotide position") + ylab("Cumulative UDP count")
  
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
  getDelReadCount = function(i) {
    sum(isPositionDel(i) * udp.df$num_reads)
  }
  
  count.plot.df$udp_count = sapply(count.plot.df$x, FUN=getDelCount)
  count.plot.df$read_count = sapply(count.plot.df$x, FUN=getDelReadCount)
  p.udp_count = ggplot(count.plot.df, aes(x=x, y=udp_count)) +
    geom_bar(stat="identity", fill="dodgerblue2") +
    simple_theme + theme(axis.title.x=element_blank()) + ylab("UDP count") +
    coord_cartesian(xlim=c(sites$start, sites$end))
  p.read_count = ggplot(count.plot.df, aes(x=x, y=read_count)) +
    geom_bar(stat="identity", fill="dodgerblue2") +
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
  
  p.full = egg::ggarrange(p.read_count, p.udp_count, p.udp_dist, ncol=1, heights=c(1,1,4),
                          top = plot_title, draw = F)
  return(p.full)
}

getDeletionProfilePlot = function(replicate.udp.df, sites, plot_title = NA, show_average = T, show_replicates = T) {
  if (!show_average & !show_replicates) {
    stop("getDeletionProfilePlot: One of show_average and show_replicates should be true.")
  }
  if (nrow(replicate.udp.df) == 0) {
    return(egg::ggarrange(textPlot("No UDPs to plot."), top=plot_title, draw = F))
  }
  xmax = nchar(replicate.udp.df$udp[1])
  # Get the deletion profile for each replicate separately
  isPositionDel = function(i) {
    (udp.char.matrix[,i] == '*')
  }
  getDelReadCount = function(i) {
    sum(isPositionDel(i) * cur.df$num_reads)
  }
  unique_reps.df = unique(replicate.udp.df %>% dplyr::select(type, replicate))
  counts.list = list()
  for (i in 1:nrow(unique_reps.df)) {
    curType = unique_reps.df[i,]$type
    curReplicate = unique_reps.df[i,]$replicate
    count.df = data.frame(x=1:xmax, type = curType, replicate = curReplicate)
    cur.df = replicate.udp.df %>% dplyr::filter(replicate == curReplicate, type == curType)
    udp.char.matrix = str_split_fixed(cur.df$udp, "", n = nchar(cur.df$udp[1]))
    count.df$del_pct = 100 * sapply(count.df$x, FUN=getDelReadCount) / sum(cur.df$num_reads)
    counts.list = c(counts.list, list(count.df))
  }
  # Combine deletion profiles for replicates
  delpct.plot.df = bind_rows(counts.list)
  delpct.plot.df$replicate = as.character(delpct.plot.df$replicate)
  
  if (show_average) {
    # Merge gDNA replicates together (and similarly for cDNA), and
    # determine deletion profiles
    # merged.udp.df = summarise(replicate.udp.df %>% group_by(udp, type),
    #                           num_reads = sum(num_reads)) %>%
    #   arrange(-num_reads)
    merged.udp.df = summarise(replicate.udp.df %>% group_by(type, udp),
                              num_reads = sum(num_reads)) %>%
      arrange(-num_reads)
    unique_reps.df = unique(merged.udp.df %>% ungroup() %>% dplyr::select(type))
    counts.list = list()
    for (i in 1:nrow(unique_reps.df)) {
      curType = unique_reps.df[i,]$type
      count.df = data.frame(x=1:xmax, type = curType)
      cur.df = merged.udp.df %>% dplyr::filter(type == curType)
      udp.char.matrix = str_split_fixed(cur.df$udp, "", n = nchar(cur.df$udp[1]))
      count.df$del_pct = 100 * sapply(count.df$x, FUN=getDelReadCount) / sum(cur.df$num_reads)
      counts.list = c(counts.list, list(count.df))
    }
    merged.delpct.plot.df = bind_rows(counts.list)
    merged.delpct.plot.df$type[merged.delpct.plot.df$type == "gDNA"] = "gDNA average"
    merged.delpct.plot.df$type[merged.delpct.plot.df$type == "cDNA"] = "cDNA average"
    merged.delpct.plot.df$replicate = merged.delpct.plot.df$type
  }
  
  delpct.plot.df$type = paste(delpct.plot.df$type, "replicate")
  alpha_values = c(`cDNA average`=0.3, `gDNA average`=0.3, `cDNA replicate`=0.9, `gDNA replicate`=0.9)
  color_values = c(`cDNA average`="firebrick1", `gDNA average`="dodgerblue3", `cDNA replicate`="red", `gDNA replicate`="blue")
  size_values = c(`cDNA average`=1.8, `gDNA average`=1.8, `cDNA replicate`=0.3, `gDNA replicate`=0.3)
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
  plot.df$type = factor(plot.df$type, levels = c("cDNA average", "gDNA average", "cDNA replicate", "gDNA replicate"))
  simple_theme = theme_bw() + theme(panel.grid.minor = element_line(colour="grey90", size=0.1), panel.grid.major = element_line(colour="grey90", size=0.1))
  p.gDNA_cDNA = ggplot() +
    geom_line(aes(x=x, y=del_pct, color=type, size=type, alpha=type, group=paste(type, replicate)), data = plot.df) +
    simple_theme + xlab("Nucleotide position") + ylab("Deletion frequency (%)") +
    coord_cartesian(xlim=c(sites$start, sites$end)) +
    scale_color_manual(values = color_values) +
    scale_size_manual(values = size_values) +
    scale_alpha_manual(values = alpha_values)
  
  if (!is.na(plot_title)) {
    p.gDNA_cDNA = p.gDNA_cDNA + ggtitle(plot_title)
  }
  if (!is.na(sites$highlight_site)) {
    p.gDNA_cDNA = p.gDNA_cDNA + geom_vline(xintercept = sites$highlight_site, color="darkgreen", alpha=0.5)
  }
  if (!is.na(sites$cut_site)) {
    p.gDNA_cDNA = p.gDNA_cDNA + geom_vline(xintercept = sites$cut_site, color="grey20", linetype = "longdash", alpha=0.5)
  }
  return(p.gDNA_cDNA)
}


getUNSData = function(replicate.udp.df, sites, region_name, min_gDNA_count = 10, min_cDNA_count = 0) {
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
  getType = function(type_rep) { sapply(type_rep, FUN=function(x) strsplit(x, "^" , fixed=T)[[1]][1]) }
  getReplicate = function(type_rep) { sapply(type_rep, FUN=function(x) strsplit(x, "^" , fixed=T)[[1]][2]) }
  replicate.udp.filled.df = replicate.udp.df %>%
    dplyr::mutate(type_replicate = paste(type, replicate, sep="^")) %>%
    dplyr::select(udp, type_replicate, num_reads) %>%
    tidyr::spread(type_replicate, num_reads, fill = 0) %>%
    tidyr::gather(key="type_replicate", value="num_reads", -udp) %>%
    dplyr::mutate(type = getType(type_replicate),
                  replicate = getReplicate(type_replicate)) %>%
    dplyr::select(udp, type, replicate, num_reads)
  
  # We want to get a dataframe with a column each for gDNA and cDNA counts
  # for each UDP and REPLICATE. But we want to retain the other columns, which
  # are largely duplicated. We do this in a couple of steps.
  # First spread the count columns out
  spread.udp.df = replicate.udp.filled.df %>% dplyr::select(udp, type, replicate, num_reads) %>%
    spread(type, num_reads, fill = 0)
  # Aggregate the remaining columns by UDP
  grouped.udp.df = replicate.udp.df %>% 
    group_by(udp) %>%
    dplyr::summarise(is_hdr_allele = any(is_hdr_allele, na.rm = T),
                     is_wt_allele = any(is_wt_allele, na.rm = T),
                     has_crispr_deletion = any(has_crispr_deletion, na.rm = T),
                     deletion_start = first(deletion_start),
                     deletion_end = first(deletion_end),
                     gDNA_total_count = sum(num_reads[type == "gDNA"]),
                     cDNA_total_count = sum(num_reads[type == "cDNA"]))
  
  # Following this, we have a row per replicate per UDP, with counts for cDNA
  # and gDNA, which are zero if the UDP wasn't observed in the replicate.
  spread.udp.df = spread.udp.df %>% dplyr::left_join(grouped.udp.df, by="udp") %>%
    dplyr::rename(gDNA_count = gDNA, cDNA_count = cDNA) %>%
    dplyr::mutate(total_count = gDNA_count + cDNA_count)
  
  wt.df = spread.udp.df %>%
    dplyr::filter(is_wt_allele) %>%
    group_by(replicate) %>%
    dplyr::summarise(udp = first(udp),
                     gDNA_count = sum(gDNA_count),
                     cDNA_count = sum(cDNA_count),
                     is_hdr_allele = F, is_wt_allele = T) %>%
    dplyr::select(udp, replicate, everything())
  
  # Results not likely to be stable if the WT count is too low
  all_wt_gDNA_count = sum(wt.df$gDNA_count)
  all_wt_cDNA_count = sum(wt.df$cDNA_count)
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
  
  wt_mean_cDNA_count =  mean(wt.df$cDNA_count)
  wt_mean_gDNA_count =  mean(wt.df$gDNA_count)
  
  # Compute the UNS for each UDP and replicate
  replicate.dels.df = spread.udp.df %>% 
    dplyr::filter(is_hdr_allele | is_wt_allele | (has_crispr_deletion & gDNA_total_count >= min_gDNA_count & cDNA_total_count >= min_cDNA_count))
  replicate.dels.df = replicate.dels.df %>% arrange(!is_wt_allele, !is_hdr_allele, -total_count)
  
  replicate.dels.df$uns = (replicate.dels.df$cDNA_count / wt_mean_cDNA_count) / (replicate.dels.df$gDNA_count / wt_mean_gDNA_count)
  
  # Here we compute (for every UDP) the UNS estimate from all replicate, and
  # the associated standard deviation using the propagation of error method.
  # First we need to add columns with the wt counts for cDNA and gDNA for
  # each matched replicate
  replicate.dels.df2 = replicate.dels.df %>%
    dplyr::left_join(wt.df %>% dplyr::select(replicate, wt_gDNA_count = gDNA_count, wt_cDNA_count = cDNA_count),
                     by = c("replicate"))
  udp.dels.df = replicate.dels.df2 %>% group_by(udp) %>%
    summarise(gDNA_total_count = sum(gDNA_count, na.rm = T),
              cDNA_total_count = sum(cDNA_count, na.rm = T),
              mean_gDNA_count = mean(gDNA_count, na.rm = T),
              mean_cDNA_count = mean(cDNA_count, na.rm = T),
              cDNA_ratio = (mean_cDNA_count / wt_mean_cDNA_count),
              sd_cDNA_ratio = vecRatioUncertaintySD(cDNA_count, wt_cDNA_count),
              gDNA_ratio = (mean_gDNA_count / wt_mean_gDNA_count),
              sd_gDNA_ratio = vecRatioUncertaintySD(gDNA_count, wt_gDNA_count),
              is_hdr_allele = first(is_hdr_allele),
              is_wt_allele = first(is_wt_allele),
              has_crispr_deletion = first(has_crispr_deletion),
              deletion_start = first(deletion_start),
              deletion_end = first(deletion_end),
              total_count = sum(total_count),
              uns = (cDNA_ratio / gDNA_ratio),
              sd_uns = ratioUncertaintySD(cDNA_ratio, sd_cDNA_ratio, gDNA_ratio, sd_gDNA_ratio, 0)) %>%
    dplyr::arrange(!is_wt_allele, !is_hdr_allele, -total_count)
  
  return(list(udp.dels.df = udp.dels.df, replicate.dels.df = replicate.dels.df))
}


getUNSPlot = function(replicate.udp.df, sites, plot_title, min_gDNA_count = 10, min_cDNA_count = 0, max_udps = 40) {
  udp_types = unique(replicate.udp.df$type)
  if (!("cDNA" %in% udp_types & "gDNA" %in% udp_types)) {
    return(egg::ggarrange(textPlot("Cannot make UNS plot without both cDNA and gDNA replicates"), top=plot_title, draw = F))
  }
  
  uns_res = getUNSData(replicate.udp.df, sites, region_name=plot_title, min_gDNA_count, min_cDNA_count)
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
                       plot.margin = unit(c(0,0,0,0), "cm"))
  
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
  p.udp = ggplot() + 
    geom_point(aes(x = pos, y = id, colour = `HDR allele`), size = dot_size, shape = 19, alpha = 0.7, data = plot.gather.df[plot.gather.df$udpchar == '-',]) +
    geom_point(aes(x = pos, y = id, colour = `HDR allele`), size = dash_size, shape = '-', alpha = 0.8, data = plot.gather.df[plot.gather.df$udpchar == '*',]) +
    geom_point(aes(x = pos, y = id, colour = `HDR allele`), size = dot_size * 1.5, shape = 19, alpha = 0.8, data = plot.gather.df[plot.gather.df$udpchar == '-' & plot.gather.df$`HDR allele`,]) +
    geom_point(aes(x = pos, y = id, colour = `HDR allele`), size = dash_size, shape = '-', alpha = 0.9, data = plot.gather.df[plot.gather.df$udpchar == '*' & plot.gather.df$`HDR allele`,]) +
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
                       plot.margin = unit(c(0,0,0,0), "cm"))
  
  if (!is.na(sites$cut_site)) {
    p.udp = p.udp + geom_vline(xintercept = sites$cut_site, color="grey20", linetype = "longdash", alpha=0.5)
  }
  
  # Left side of the plot, showing UDP-normalised scores (UNS)
  uns.plot.df = data.frame(y = label(ddata)$x,
                           x = 0,
                           log2uns = log2(udp.dels.df$uns + 1),
                           cDNA_lower = udp.dels.df$uns < 1,
                           gDNA = udp.dels.df$gDNA_total_count,
                           log2uns_conf_hi = log2(udp.dels.df$uns + 1.96 * udp.dels.df$sd_uns + 1),
                           log2uns_conf_lo = log2(udp.dels.df$uns + 1.96 * udp.dels.df$sd_uns + 1))
  uns.plot.df = data.frame(y = label(ddata)$x,
                           x = 0,
                           uns = udp.dels.df$uns,
                           cDNA_lower = udp.dels.df$uns < 1,
                           gDNA = udp.dels.df$gDNA_total_count,
                           uns_conf_hi = udp.dels.df$uns + 1.96 * udp.dels.df$sd_uns,
                           uns_conf_lo = udp.dels.df$uns - 1.96 * udp.dels.df$sd_uns)
  
  uns.replicates.plot.df = udp.del_replicates.df %>% 
    dplyr::select(udp, uns, gDNA_count) %>%
    dplyr::left_join(udp.order.map, by="udp")
  
  maxDotSize = 5.5
  if (numUDPs > 60) {
    maxDotSize = 3.5
  } else if (numUDPs > 30) {
    maxDotSize = 4.5
  }
  
  min_uns_threshold = 0.25
  uns.plot.df$uns = sapply(uns.plot.df$uns, FUN = function(x) max(x, min_uns_threshold))
  uns.plot.df$uns_conf_lo = sapply(uns.plot.df$uns_conf_lo, FUN = function(x) max(x, 0.01))
  max_uns_display = max(1.5, ceiling(max(uns.plot.df$uns)))
  min_uns_display = min(uns.plot.df$uns)
  
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
    scale_color_manual(values=c("red", "blue")) +
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
                       plot.margin = unit(c(0,0,0,0), "cm"))
  
  p.udp_profile = egg::ggarrange(p.dendro, p.udp, p.uns, nrow=1, ncol=3, widths=c(1,4,1), top = plot_title, draw = F)
  #p.udp_profile = plot_grid(p.uns, p.udp, p.dendro, nrow=1, ncol=3, rel_widths = c(2,4,1))
  return(p.udp_profile)
}


fitVarianceComponents = function(replicate.udp.spread.df, replicates.type.df) {
  replicate.udp.spread.df = replicate.udp.spread.df %>% as.data.frame()
  rownames(replicate.udp.spread.df) = replicate.udp.spread.df$udp
  
  # Get all columns of the replicate dataframe which begin with "replicate"
  replicate_cols = colnames(replicates.df)[grepl("^replicate_", colnames(replicates.df))]
  for (col in replicate_cols) {
    replicates.type.df[, col] = as.factor(replicates.type.df[, col, drop=T])
  }
  
  #form <- ~ (1|replicate_cDNA)
  #vp <- fitExtractVarPartModel(replicate.udp.spread.df %>% dplyr::select(-udp), form, replicates.type.df)
  formulaStr = sprintf("~ (1|%s)", replicate_cols[1])
  for (col in replicate_cols[-1]) {
    formulaStr = sprintf("%s + (1|%s)", formulaStr, col)
  }
  vp <- fitExtractVarPartModel(replicate.udp.spread.df %>% dplyr::select(-udp), as.formula(formulaStr), replicates.type.df)
  #View(data.frame(vp))
  vp
}

doVariancePartitionPlot = function(vp, residuals=T, color = NA) {
  vp.df = data.frame(vp) %>% 
    mutate(udp = rownames(vp)) %>%
    gather(variable, variance, -udp) %>%
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
  p = ggplot(vp.df, aes(x=variable, y=pctvariance, fill="grey")) +
    geom_violin(scale="width") +
    geom_boxplot(width=0.2, fill="grey80", outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.6, mapping = aes(color = "black")) +
    theme_bw() +
    theme(legend.position="none") + 
    theme(axis.text.x=element_text(angle=0), axis.title.x=element_blank(), panel.grid=element_blank()) +
    ylab("Variance explained (%)") +
    scale_fill_manual(values = c("gray95"))
  #  scale_fill_manual(values = c(gg_color_hue(numVariables-1), "gray85"))
  if (!is.na(color)) {
    p = p + scale_color_manual(values = c(rep(color, numVariables)))
  }
  return(p)
}

getVarianceComponentsPlots = function(replicate.udp.df, replicates.df, method = "log_reads", min_udp_total_count = 100) {
  replicate_cols = colnames(replicates.df)[grepl("^replicate_", colnames(replicates.df))]
  if (length(replicate_cols) < 1) {
    stop("--variance_analysis option given but no columns beginning with 'replicate_' found in the --replicates file")
  }
  
  filterUDPs = function(udp.df, min_udp_total_count) {
    # Spread the replicates out into columns (cDNA_1, cDNA_2, etc.), with
    # fill = 0 for when a replicate has no value for a given UDP, then
    # gather again.
    udp.filled.df = udp.df %>%
      dplyr::select(udp, type, replicate, num_reads) %>%
      tidyr::spread(replicate, num_reads, fill = 0) %>%
      tidyr::gather(key="replicate", value="num_reads", -udp, -type)
    udp.summary = udp.filled.df %>% group_by(udp, type) %>%
      dplyr::summarise(sum = sum(num_reads))
    
    # Remove UDPs with too few total reads
    udp.filled.df %>% dplyr::inner_join(udp.summary %>% dplyr::filter(sum >= min_udp_total_count), by=c("udp", "type"))
  }
  
  # This function returns a value based on the "num_reads" for each replicate,
  # but which is what we use to determine variance components. Raw number of
  # reads is not really suitable. The recommended value is "deviation_from_poisson",
  # where this represents the distance of a given replicate from the mean number
  # of reads across replicates, ajusted for the expected variability for a
  # poisson process.
  getValueForVariance = function(replicate.udp.df, method) {
    if (method == "reads") {
      value = replicate.udp.df$num_reads
    }
    if (method == "log_reads") {
      value = log(replicate.udp.df$num_reads + 1)
    }
    else if (method == "deviation_from_poisson") {
      mean_reads_per_udp.df= replicate.udp.df %>% dplyr::group_by(udp) %>%
        dplyr::summarise(mean_reads = mean(num_reads))
      replicate.udp.df = replicate.udp.df %>% dplyr::left_join(mean_reads_per_udp.df, by="udp")
      value = abs(replicate.udp.df$num_reads - replicate.udp.df$mean_reads) / sqrt(replicate.udp.df$mean_reads)
    } else {
      stop(sprintf("Unrecognized variance plots method: '%s'", method))
    }
    value
  }
  
  #########################################################################
  dnaType = "cDNA"
  replicate.udp.cDNA.df = filterUDPs(replicate.udp.df %>% dplyr::filter(type == dnaType), min_udp_total_count)
  
  replicate.udp.cDNA.df$value = getValueForVariance(replicate.udp.cDNA.df, method)
  
  # Spread the replicates out into columns (cDNA_1, cDNA_2, etc.), with
  # fill = 0 for when a replicate has no value for a given UDP
  replicate.udp.cDNA.spread.df = replicate.udp.cDNA.df %>%
    dplyr::select(udp, replicate, value) %>%
    tidyr::spread(replicate, value, fill = 0)
  vp = fitVarianceComponents(replicate.udp.cDNA.spread.df,
                             replicates.df %>% dplyr::filter(type == dnaType))
  vp.cDNA.df = data.frame(vp)
  p.cDNA = doVariancePartitionPlot(vp, color="red") + ggtitle(dnaType)
  
  #########################################################################
  # Do the same for gDNA
  dnaType = "gDNA"
  replicate.udp.gDNA.df = filterUDPs(replicate.udp.df %>% dplyr::filter(type == dnaType), min_udp_total_count)
  replicate.udp.gDNA.df$value = getValueForVariance(replicate.udp.gDNA.df, method)
  
  replicate.udp.gDNA.spread.df = replicate.udp.gDNA.df %>%
    dplyr::select(udp, replicate, value) %>%
    tidyr::spread(replicate, value, fill = 0)
  vp = fitVarianceComponents(replicate.udp.gDNA.spread.df,
                             replicates.df %>% dplyr::filter(type == dnaType))
  vp.gDNA.df = data.frame(vp)
  p.gDNA = doVariancePartitionPlot(vp, color="blue") + ggtitle(dnaType)
  
  cDNA.summary.df = variancePartitionSummaryTable(vp.cDNA.df)
  gDNA.summary.df = variancePartitionSummaryTable(vp.gDNA.df)
  
  ggThemeBlank = theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                                    axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
                                    panel.border = element_blank())
  myTableTheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.75)),
                                            colhead = list(fg_params=list(cex = 0.75)),
                                            rowhead = list(fg_params=list(cex = 0.75)))
  
  p.stats1 = ggplot(data.frame(x=1:10, y=1:10), aes(x,y)) + geom_blank() + ggThemeBlank +
    annotation_custom(tableGrob(cDNA.summary.df, theme = myTableTheme), xmin=1, xmax=5, ymin=2, ymax=8.5) +
    annotation_custom(tableGrob(gDNA.summary.df, theme = myTableTheme), xmin=6, xmax=10, ymin=2, ymax=8.5) +
    # annotate("text", x=3.5, y=8, label = "cDNA", hjust = 0.5, fontface = 2, size = 5) +
    # annotate("text", x=7.5, y=8, label = "gDNA", hjust = 0.5, fontface = 2, size = 5) +
    annotate("text", x=5, y=10, label = sprintf("%s variance components", replicates.df$name[1]), vjust = 1, fontface = 2, size = 5) +
    annotate("text", x=5, y=9, label = sprintf("Method = %s", method), vjust = 1, size = 3) +
    annotate("text", x=1.5, y=1, label = sprintf("Based on %d UDPs with read count > %d", nrow(vp.cDNA.df), min_udp_total_count),
             hjust = 0, vjust = 0, size = 3) +
    annotate("text", x=6.5, y=1, label = sprintf("Based on %d UDPs with read count > %d", nrow(vp.gDNA.df), min_udp_total_count),
             hjust = 0, vjust = 0, size = 3)
  
  #egg::ggarrange(p.stats1, grid.arrange(p.cDNA, p.gDNA, ncol=2), nrow = 1)
  plot_grid(p.stats1, plot_grid(p.cDNA, p.gDNA, nrow = 1), nrow = 2, rel_heights = c(0.4, 0.6))
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
                                       top=powerPlotTitle, draw = F))
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
    xlab("UDP read %") + ylab("Power") +
    ggtitle(label = powerPlotTitle,
            subtitle = sprintf("Based on %d cDNA and %d gDNA replicates", n_cDNA_rep, n_gDNA_rep)) +
    theme(panel.grid.major = element_line(size = 1, colour = "grey95"),
          plot.title = element_text(size=14, hjust=0.5, face="bold"),
          plot.subtitle = element_text(size=12, hjust=0.5))
  
  plot_list = list(cv_plots = p.cv_plots, power = p.power)
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

getHDRRatioEstimate = function(stats.df, colToCompare = "num_hdr_reads") {
  cDNA.df = stats.df %>% dplyr::filter(type == "cDNA") %>% as.data.frame()
  gDNA.df = stats.df %>% dplyr::filter(type == "gDNA") %>% as.data.frame()
  
  cDNA_ratio = mean(cDNA.df[,colToCompare]) / mean(cDNA.df$num_wt_reads)
  gDNA_ratio = mean(gDNA.df[,colToCompare]) / mean(gDNA.df$num_wt_reads)
  sd_cDNA_ratio = vecRatioUncertaintySD(cDNA.df[,colToCompare], cDNA.df$num_wt_reads)
  sd_gDNA_ratio = vecRatioUncertaintySD(gDNA.df[,colToCompare], gDNA.df$num_wt_reads)
  SE_cDNA_ratio = sd_cDNA_ratio / sqrt(nrow(cDNA.df))
  SE_gDNA_ratio = sd_gDNA_ratio / sqrt(nrow(gDNA.df))
  
  effect = cDNA_ratio / gDNA_ratio
  effect_se = ratioUncertaintySD(cDNA_ratio, SE_cDNA_ratio, gDNA_ratio, SE_gDNA_ratio, 0)
  
  Zscore = (effect - 1) / effect_se
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


###########################################################################
##' commandArgs parsing
##' 
##' return a named list of command line arguments
##'
##' Usage:
##' call the R script thus
##'   ./myfile.R --args myarg=something
##' or
##'   R CMD BATCH --args myarg=something myfile.R
##'
##' Then in R do
##'   myargs <- getArgs()
##' and myargs will be a named list
##' > str(myargs)
##' List of 2
##' $ file : chr "myfile.R"
##' $ myarg: chr "something"
##'
##' @title getArgs
##' @param verbose print verbage to screen 
##' @param defaults a named list of defaults, optional
##' @return a named list
##' @author Chris Wallace
getArgs = function(verbose=FALSE, defaults=NULL) {
  myargs <- gsub("^--","",commandArgs(TRUE))
  setopts <- !grepl("=",myargs)
  if(any(setopts))
    myargs[setopts] <- paste(myargs[setopts],"=notset",sep="")
  myargs.list <- strsplit(myargs,"=")
  myargs <- lapply(myargs.list,"[[",2 )
  names(myargs) <- lapply(myargs.list, "[[", 1)
  
  ## logicals
  if(any(setopts))
    myargs[setopts] <- TRUE
  
  ## defaults
  if(!is.null(defaults)) {
    defs.needed <- setdiff(names(defaults), names(myargs))
    if(length(defs.needed)) {
      myargs[ defs.needed ] <- defaults[ defs.needed ]
    }
  }
  
  ## verbage
  if(verbose) {
    cat("read",length(myargs),"named args:\n")
    print(myargs)
  }
  myargs
}


###########################################################################

system.time( main() )
