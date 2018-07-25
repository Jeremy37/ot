#!/usr/bin/env Rscript

# This script reads in a VCF file and an "ASECounts" file. The VCF file has the genotypes
# for a single individual in the region of 1 gene. The ASECounts file has allele-specific
# read counts for possibly multiple samples that come from the individual and which cover
# the gene. The script performs statistical tests for allele-specific expression, and
# produces some plots.
library(tidyverse)
library(stringr)
library(egg) # for ggarrange
library(cowplot) # could probably use just one of these two plotting libraries
library(aod) # for betabin
library(ggdendro)
library(optparse)
options(stringsAsFactors = F)

#myargs <- NULL
opt = NULL

setwd("/Users/jeremys/work/opentargets/experiment/transcribed/batch1")


main = function()
{
  option_list <- list(
    make_option(c("--regions"), type="character", default=NULL, help=""),
    make_option(c("--replicates"), type="character", default=NULL, help=""),
    make_option(c("--out"), type="character", default=NULL, help=""),
    make_option(c("--minMapQ"), type="integer", default=0, help="[default %default]"),
    make_option(c("--subsample"), type="numeric", default=NULL, help=""),
    make_option(c("--max_mismatch_frac"), type="numeric", default=0.05, help="[default %default]"),
    make_option(c("--viewing_window"), type="integer", default=40, help="[default %default]"),
    make_option(c("--editing_window"), type="integer", default=10, help="[default %default]"),
    make_option(c("--min_window_overlap"), type="integer", default=30, help="[default %default]"),
    make_option(c("--exclude_multiple_deletions"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--exclude_nonspanning_reads"), type="logical", default=T, action="store_true", help="[default %default]"),
    make_option(c("--exclude_nonspanning_deletions"), type="logical", default=T, action="store_true", help="[default %default]"),
    make_option(c("--uns_plot_min_gDNA"), type="integer", default=10, help="[default %default]"),
    make_option(c("--uns_plot_min_cDNA"), type="integer", default=0, help="[default %default]"),
    make_option(c("--uns_plot_max_udps"), type="integer", default=40, help="[default %default]"),
    make_option(c("--no_allele_profile"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--no_site_profile"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--no_udp_profile"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--no_stats"), type="logical", default=F, action="store_true", help="[default %default]")
  )
  # 
  # opt <<- parse_args(OptionParser(option_list=option_list))

  opt <<- list(#input = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1/analysis/udp_replicates.tsv",
                  regions = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1/analysis/regions.tsv",
                  replicates = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1/analysis/replicates.tsv",
                  out = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1/analysis/batch1.all_dels",
                  minMapQ = 0,
                  subsample = NULL,
                  max_mismatch_frac = 0.05,
                  viewing_window = 30, # window around the cut site for viewing
                  # viewing_window = NULL,
                  editing_window = 10, # window around the cut site for counting deletions as edits
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
                  no_stats = F)
  print(opt)
  
  # if (is.null(opt$minMapQ)) { opt$minMapQ <<- 0 }
  # if (is.null(opt$max_mismatch_frac)) { opt$max_mismatch_frac <<- 0.05 }
  # if (is.null(opt$editing_window)) { opt$editing_window <<- 10 }
  # if (is.null(opt$viewing_window)) { opt$viewing_window <<- 1000 }
  # if (is.null(opt$min_window_overlap)) { opt$min_window_overlap <<- 30 }
  # if (is.null(opt$exclude_multiple_deletions)) { opt$exclude_multiple_deletions <<- F }
  # if (is.null(opt$exclude_nonspanning_reads)) { opt$exclude_nonspanning_reads <<- T }
  # if (is.null(opt$exclude_nonspanning_deletions)) { opt$exclude_nonspanning_deletions <<- T }
  # if (is.null(opt$save_allele_profile)) { opt$save_allele_profile <<- F }
  # if (is.null(opt$save_site_profile)) { opt$save_site_profile <<- F }
  
  regions.df = readr::read_tsv(opt$regions, col_types="ccciiiiccc")
  #regions.df = regions.df %>% dplyr::filter((index %in% c("29")))
  #regions.df = regions.df %>% dplyr::filter((index %in% c("4", "7", "8", "28", "31")))
  #regions.df = regions.df %>% dplyr::filter(index %in% c("22"))
  #regions.df = regions.df %>% dplyr::filter(index %in% c("25", "27"))
  #regions.df = regions.df %>% dplyr::filter(index %in% c("1", "3", "12", "15", "17", "18", "24", "26"))
  #input.df = input.df %>% dplyr::filter(name %in% c("1", "3", "12", "15", "17", "18", "24", "26"))
  
  replicates.df = readr::read_tsv(opt$replicates, col_types="ccicc")
  
  region_names = unique(regions.df$name)
  all_region_plots = list()
  all_results = list()
  region_index = 1
  #profvis({
  for (region_name in region_names) {
    df = regions.df %>% dplyr::filter(name == region_name)
    if (nrow(df) != 1) {
      stop("Error: multiple regions with the same region name (%s). Each input line for the --regions file should have a unique name.")
    }
    cur_region = as.list(df[1,])
    cur_replicates.df = replicates.df %>% dplyr::filter(name == region_name)
    result = doFullDeletionAnalysis(cur_region, cur_replicates.df)
    result$region = region_name
    result$replicate.df = cur_replicates.df
    
    all_results[[region_index]] = result
    all_region_plots[[region_index]] = result$plot_list
    region_index = region_index + 1
  }
  #})
  
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
  ggThemeBlank = theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                                    axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())
  myTableTheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.8)),
                                            colhead = list(fg_params=list(cex = 0.8)),
                                            rowhead = list(fg_params=list(cex = 0.8)))
  settingsPlot = ggplot(data.frame(x=1:10, y=1:10), aes(x,y)) + geom_blank() + ggThemeBlank +
    annotation_custom(tableGrob(settings.df, theme = myTableTheme), xmin=1, xmax=9, ymin=1, ymax=9) +
    ggtitle("CRISPR editing analysis settings")
  
  fname = sprintf("%s.all_plots.pdf", opt$out)
  pdf(file = fname, width = 8, height = 7)
  print(settingsPlot)
  print(all_region_plots)
  dev.off()
  
  fname = sprintf("%s.summary_plots.pdf", opt$out)
  pdf(file = fname, width = 8, height = 7)
  for (region_plots in all_region_plots) {
    print(region_plots$stats)
    print(region_plots$merged_udp)
    print(region_plots$merged_del_profile)
    print(region_plots$merged_UNS)
  }
  dev.off()
  
  getReplicatesDF = function(all_results, dfname) {
    bind_rows(sapply(all_results, 
                     function(res) {
                       sapply(res$replicate_list, function(res2) {
                         df = res2[[dfname]]
                         df$name = res2$name
                         df$replicate = res2$replicate
                         df$type = res2$type
                         df %>% dplyr::select(name, replicate, type, everything())
                       })
                     }) )
  }
  
  if (!opt$no_allele_profile) {
    df = bind_rows(lapply(all_results, function(res) res$wt_hdr.df))
    fname = sprintf("%s.mismatch_profile.tsv", opt$out)
    write.table(df, fname, quote=F, row.names=F, col.names=T, sep="\t")
  }
  if (!opt$no_udp_profile) {
    df = bind_rows(lapply(all_results, function(res) res$replicate.udp.df))
    fname = sprintf("%s.replicate_udps.tsv", opt$out)
    write.table(df, fname, quote=F, row.names=F, col.names=T, sep="\t")
    
    df = bind_rows(lapply(all_results, function(res) res$merged.udp.df))
    fname = sprintf("%s.merged_udps.tsv", opt$out)
    write.table(df, fname, quote=F, row.names=F, col.names=T, sep="\t")
  }
  if (!opt$no_site_profile) {
    df = bind_rows(lapply(all_results, function(res) res$site.profiles.df))
    fname = sprintf("%s.site_profiles.tsv", opt$out)
    write.table(df, fname, quote=F, row.names=F, col.names=T, sep="\t")
  }
  if (!opt$no_stats) {
    df = bind_rows(lapply(all_results, function(res) res$stats.df))
    fname = sprintf("%s.stats.tsv", opt$out)
    write.table(df, fname, quote=F, row.names=F, col.names=T, sep="\t")
  }
}


doFullDeletionAnalysis = function(region, replicates.df)
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
  sites = list(seq_start = region$start,
               seq_end = region$end,
               seq_highlight_site = region$highlight_site,
               seq_cut_site = region$cut_site)
  sites = getRelativeSites(sites, opt$viewing_window)
  for (i in 1:nrow(replicates.df)) {
    result = doReplicateDeletionAnalysis(locus_name,
                                         replicates.df[i,]$replicate,
                                         replicates.df[i,]$type,
                                         replicates.df[i,]$bam,
                                         region$sequence_name,
                                         sites,
                                         region$hdr_allele_profile,
                                         region$wt_allele_profile,
                                         region$ref_sequence)
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
        dplyr::select(name, replicate, type, everything())
    }
    
    replicate_udps[[i]] = result$udp.df %>% dplyr::mutate(name = locus_name,
                                                          replicate = replicates.df[i,]$replicate,
                                                          type = replicates.df[i,]$type) %>%
      dplyr::select(name, replicate, type, everything())
    
    plot_title = sprintf("%s replicate %d, %s", locus_name, replicates.df[i,]$replicate, replicates.df[i,]$type)
    replicate_plots[[i]] = getUDPPlot(result$udp.df, plot_title, highlight_site = sites$highlight_site, cut_site = sites$cut_site)
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
    site.profiles.df = bind_rows(site.profiles.df, site.profile.merged)
  }
  
  stats.df = bind_rows(lapply(replicate_list, function(res) as.data.frame(res$stats)))
  stats.df$HDR_WT_ratio = (stats.df$num_hdr_reads / stats.df$num_wt_reads)
  stats.df$HDR_rate = stats.df$num_hdr_reads / stats.df$num_reads
  stats.df$DEL_WT_ratio = (stats.df$num_deletion_reads / stats.df$num_wt_reads)
  stats.df$DEL_rate = stats.df$num_deletion_reads / stats.df$num_reads
  stats.df$editing_rate = stats.df$num_edit_reads / stats.df$num_reads
  
  # Summary of stats from the propagation of errors method
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
  if (sum(stats.df$type == "cDNA") < 1 | sum(stats.df$type == "gDNA") < 1) {
    summary.method.betabin = sprintf("Beta binomial\nUnable to calculate\nOnly %d gDNA and %d cDNA replicates",
                                     sum(stats.df$type == "gDNA"), sum(stats.df$type == "cDNA"))
  } else {
    test.df.hdr = stats.df %>%
      dplyr::mutate(name = paste(name, replicate, type, sep="_"),
                    is_cDNA = as.integer(type == "cDNA")) %>%
      dplyr::select(name, is_cDNA, wt = num_wt_reads, hdr = num_hdr_reads)
    #fm1 = betabin(cbind(hdr, wt) ~ is_cDNA, ~ 1, data = test.df[,2:4])
    test.df.hdr$is_cDNA[3] = T
    fm1 = NULL
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
    fm2 = NULL
    res_del = tryCatch(fm2 <- betabin(cbind(del, wt) ~ is_cDNA, ~ 1, data = test.df.del[,2:4]),
                       error = function(e) e, warning = function(w) w)
    betabin_del_star = ""
    if (is(res_del, "warning") | is(res_del, "error")) {
      betabin_del_star = "**"
      betabin_warn_err = paste("**", trimws(res_del$message))
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
                                             reads_excluded_for_multiple_deletions)
  stats.plot.df = stats.plot.df %>% dplyr::rename("HDR reads" = num_hdr_reads,
                                                  "WT reads" = num_wt_reads,
                                                  "Deletion reads" = num_deletion_reads,
                                                  "excluded-minoverlap" = reads_excluded_for_minoverlap,
                                                  "excluded-mismatches" = reads_excluded_for_mismatches,
                                                  "excluded-mult.deletions" = reads_excluded_for_multiple_deletions)
  stats.plot.df$HDR_WT_ratio = sprintf("%.3g", stats.plot.df$HDR_WT_ratio)
  stats.plot.df$DEL_WT_ratio = sprintf("%.3g", stats.plot.df$DEL_WT_ratio)
  stats.plot.df$HDR_rate = sprintf("%.3f%%", 100 * stats.plot.df$HDR_rate)
  stats.plot.df$DEL_rate = sprintf("%.3f%%", 100 * stats.plot.df$DEL_rate)
  stats.plot.df$editing_rate = sprintf("%.3f%%", 100 * stats.plot.df$editing_rate)
  ggThemeBlank = theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                                    axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())
  myTableTheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.8)),
                                            colhead = list(fg_params=list(cex = 0.8)),
                                            rowhead = list(fg_params=list(cex = 0.8)))
  p.stats1 = ggplot(data.frame(x=1:10, y=1:10), aes(x,y)) + geom_blank() + ggThemeBlank +
    annotation_custom(tableGrob(t(stats.plot.df), theme = myTableTheme), xmin=1, xmax=9, ymin=1, ymax=8) +
    annotate("text", x=5, y=10, label = sprintf("%s summary", stats.df$name[1]), vjust = 1, fontface = 2, size = 5) +
    annotate("text", x=1, y=9.3, label = summary.method.prop, vjust = 1, hjust = 0, size = 3.5) +
    annotate("text", x=9.5, y=9.3, label = summary.method.betabin, vjust = 1, hjust = 1, size = 3.5)
  #p.stats1
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
  
  # Merge together results from all replicates
  replicate.udp.df = bind_rows(replicate_udps)

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
  p.udp_gDNA = getUDPPlot(merged.udp.df %>% dplyr::filter(type == "gDNA"), "gDNA",
                          highlight_site = sites$highlight_site, cut_site = sites$cut_site)
  p.udp_cDNA = getUDPPlot(merged.udp.df %>% dplyr::filter(type == "cDNA"), "cDNA",
                          highlight_site = sites$highlight_site, cut_site = sites$cut_site)
  # now add the title
  p.title = ggdraw() + draw_label(region_title, fontface='bold')
  p.cDNA_gDNA = plot_grid(p.udp_gDNA, p.udp_cDNA, nrow=1)
  p.merged_UDP = plot_grid(p.title, p.cDNA_gDNA, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
  
  p.merged_del_profile = getDeletionProfilePlot(replicate.udp.df,
                                                plot_title = "Merged replicates",
                                                highlight_site = sites$highlight_site, cut_site = sites$cut_site,
                                                show_average = T, show_replicates = F)
  p.merged_del_profile = p.merged_del_profile + theme(axis.title.x=element_blank())
  p.replicate_del_profile = getDeletionProfilePlot(replicate.udp.df,
                                                   plot_title = "Individual replicates",
                                                   highlight_site = sites$highlight_site, cut_site = sites$cut_site,
                                                   show_average = T, show_replicates = T)
  p.del_profile = ggarrange(ggdraw() + draw_label(region$name, fontface='bold'),
                            p.merged_del_profile, p.replicate_del_profile, 
                            ncol=1, heights=c(0.1, 0.5, 0.5), draw = F)
  
  p.merged_UNS = getUNSPlot(merged.udp.df, plot_title = region_title,
                            min_gDNA_count = opt$uns_plot_min_gDNA, min_cDNA_count = opt$uns_plot_min_cDNA,
                            max_udps = opt$uns_plot_max_udps, cut_site = sites$cut_site)
  
  plot_list = list(stats=p.stats1,
                   merged_udp = p.merged_UDP,
                   merged_del_profile = p.del_profile,
                   merged_UNS = p.merged_UNS,
                   replicate_plots = replicate_plots)
  
  return(list(replicate_list = replicate_list,
              plot_list = plot_list,
              replicate.udp.df = replicate.udp.df,
              merged.udp.df = merged.udp.df,
              site.profiles.df = site.profiles.df,
              wt_hdr.df = wt_hdr.df,
              stats.df = stats.df,
              stats.summary = stats.summary))
}

doReplicateDeletionAnalysis = function(name, replicate, type, bam_file, sequence_name, sites, hdr_profile, wt_profile, ref_sequence)
{
  cat(sprintf("\n\nAnalysing region %s, replicate %d, file %s, %s:%d-%d, highlight site %d, cut site %d\n",
              name, replicate, bam_file, sequence_name, sites$seq_start, sites$seq_end, sites$seq_highlight_site, sites$seq_cut_site))
  cat(sprintf("HDR allele: %s\n", hdr_profile))
  cat(sprintf("WT allele:  %s\n", wt_profile))
  cat(sprintf("REF sequence: %s\n", ref_sequence))
  stats = list(name = name, replicate = replicate, type = type)
  
  read_data = getRegionReadsFromBam(name, replicate, bam_file, sequence_name, sites$seq_start, sites$seq_end, ref_sequence)
  
  cat(sprintf("%d cDNA reads of %d total (%.1f%%) were soft-clipped\n", read_data$num_softclipped, read_data$num_reads, 100.0 * read_data$num_softclipped / read_data$num_reads))
  cat(sprintf("%d cDNA reads of %d total (%.1f%%) were hard-clipped\n", read_data$num_hardclipped, read_data$num_reads, 100.0 * read_data$num_hardclipped / read_data$num_reads))
  cat(sprintf("%d cDNA reads of %d total (%.1f%%) had insertions\n", read_data$num_insertion, read_data$num_reads, 100.0 * read_data$num_insertion / read_data$num_reads))
  stats$num_softclipped = read_data$num_softclipped
  stats$num_hardclipped = read_data$num_hardclipped
  stats$num_insertion = read_data$num_insertion
  stats$num_reads = read_data$num_reads
  
  reads.df = data.frame(region_read = read_data$alignedReads)
  reads.df$read_chars = sapply(reads.df$region_read, FUN=function(s) strsplit(s, ""))
  
  region_length = (sites$seq_end - sites$seq_start + 1)
  if (nchar(wt_profile) != region_length) {
    stop(sprintf("The WT allele profile given has length %d, but the region size (end - start + 1) is %d", nchar(wt_profile), region_length))
  }
  
  if (is.na(hdr_profile)) {
    hdr_profile = ""
  }
  
  sites = getRelativeSites(sites, opt$viewing_window)
  if (!is.null(opt$viewing_window)) {
    # Subset everything to a window around the cut site - including HDR profile, WT profile, ref_sequence, region read
    reads.df$region_read = sapply(reads.df$region_read, FUN=function(s) {substr(s, sites$start, sites$end)})
    reads.df$read_chars = sapply(reads.df$read_chars, FUN=function(x) {list(x[sites$start:sites$end])})
    if (hdr_profile != "") {
      hdr_profile = substr(hdr_profile, sites$start, sites$end)
    }
    wt_profile = substr(wt_profile, sites$start, sites$end)
    ref_sequence = substr(ref_sequence, sites$start, sites$end)
  }
  region_length = (sites$end - sites$start + 1)
  
  hdr_profile_chars = strsplit(hdr_profile, "")[[1]]
  wt_profile_chars = strsplit(wt_profile, "")[[1]]
  ref_seq_chars = strsplit(ref_sequence, "")[[1]]
  
  stats$num_reads = nrow(reads.df)
  
  reads.df$seq_length = sapply(reads.df$region_read, FUN=function(s) str_count(s, "[ACGT]"))
  #reads.df$seq_length = sapply(reads.df$region_read, FUN=function(s) nchar(gsub("[-*]", "", s)))  # SLOWER
  
  reads.df$mismatch_count = sapply(reads.df$read_chars, getMismatchCharsCount, ref_seq_chars)
  #reads.df$mismatch_count = sapply(reads.df$region_read, getMismatchCount, ref_sequence)  # SLOWER
  
  reads.df$spanning_read = T
  span_site = sites$cut_site
  if (!is.na(sites$highlight_site)) {
    span_site = sites$highlight_site
  }
  reads.df$spanning_read = sapply(reads.df$region_read, FUN=function(s) (substring(s, span_site, span_site) != '-'))
  if (opt$exclude_nonspanning_reads) {
    stats$nonspanning_reads = sum(!reads.df$spanning_read)
    cat(sprintf("%d of %d reads (%.2f%%) excluded due to not spanning the site of interest (position %d)\n",
                stats$nonspanning_reads, stats$num_reads, 100.0 * stats$nonspanning_reads / stats$num_reads,
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
  reads.df$udp = sapply(reads.df$read_chars, FUN=getReadCharsUDP, wt_profile_chars)
  #reads.df$udp = sapply(reads.df$region_read, FUN=getReadUDP, wt_profile_chars) # SLOWER
  
  reads.df$has_any_deletion = sapply(reads.df$udp, FUN=function(s) grepl("[*]", s))
  if (opt$exclude_nonspanning_deletions) {
    reads.df$has_crispr_deletion = grepl("[*]", substr(reads.df$udp, sites$cut_site - opt$editing_window, sites$cut_site + opt$editing_window))
    
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
    reads.df$udp = sapply(reads.df$udp, updateUDPEdits, sites$cut_site, opt$editing_window)
  } else {
    # Accept deletions anywhere in the read
    reads.df$has_crispr_deletion = reads.df$has_any_deletion
  }
  # Include deletions anywhere in the read
  reads.df$has_multiple_deletions = sapply(reads.df$udp, FUN=function(s) grepl("[*]+[^*]+[*]+", s))
  
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
  exclude_for_mismatches = reads.df$mismatch_count / reads.df$seq_length > opt$max_mismatch_frac
  stats$reads_excluded_for_mismatches = sum(exclude_for_mismatches)
  cat(sprintf("%d of %d reads (%.2f%%) excluded due to having more than %.2f%% of the read being mismatches\n",
              stats$reads_excluded_for_mismatches, stats$num_reads, 100.0 * stats$reads_excluded_for_mismatches / stats$num_reads,
              opt$max_mismatch_frac * 100))
  reads.df = reads.df[!exclude_for_mismatches, ]
  
  stats$num_wt_reads = sum(reads.df$is_wt_allele, na.rm = T)
  stats$num_hdr_reads = sum(reads.df$is_hdr_allele, na.rm = T)
  stats$num_deletion_reads = sum(reads.df$has_crispr_deletion, na.rm = T)
  
  # Aggregate reads according to their UDP
  udp.df = summarise(reads.df %>% group_by(udp),
                     num_reads = n(),
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
    wt_hdr.reads.df$mismatch_profile = sapply(wt_hdr.reads.df$read_chars, FUN=getReadCharsMismatchProfile, ref_seq_chars)
    wt_hdr.df = summarise(wt_hdr.reads.df %>% group_by(mismatch_profile),
                          num_reads = n(),
                          mismatch_count = first(mismatch_count),
                          spanning_read = first(spanning_read),
                          is_wt_allele = first(is_wt_allele),
                          is_hdr_allele = first(is_hdr_allele))
  }
  
  # Make a table which has all variations of read sequences at the "profile" sites,
  # i.e. those sites used to identify WT and HDR alleles
  sites.profile.df = NA
  if (!opt$no_site_profile) {
    wt_profile_positions = which(isDNALetter(wt_profile_chars))
    reads.df$sites_profile = sapply(reads.df$region_read, getReadSiteProfile, wt_profile_positions)
    sites.profile.df = reads.df %>% group_by(sites_profile) %>%
      summarise(count = n())
  }
  
  return(list(udp.df = udp.df, wt_hdr.df = wt_hdr.df, sites.profile.df = sites.profile.df, stats = stats))
}

getRelativeSites = function(sites, viewing_window = NULL) {
  region_length = (sites$seq_end - sites$seq_start + 1)
  sites$cut_site = sites$seq_cut_site - sites$seq_start + 1
  sites$highlight_site = sites$seq_highlight_site - sites$seq_start + 1
  sites$start = 1
  sites$end = region_length
  if (!is.null(viewing_window)) {
    # Subset everything to a window around the cut site - including HDR profile, WT profile, ref_sequence, region read
    sites$start = sites$cut_site - viewing_window
    if (sites$start < 1) {
      sites$start = 1
    }
    sites$end = sites$cut_site + viewing_window - 1
    if (sites$end > region_length) {
      sites$end = region_length
    }
    sites$cut_site = sites$cut_site - sites$start + 1
    sites$highlight_site = sites$highlight_site - sites$start + 1
  }
  return(sites)
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
  
  read_data = getAlignedReads(sam_reads, ref_sequence, start, end)
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
  alignedReads = vector(mode = "character", length = length(sam_reads)) # Preallocate a vector of aligned reads
  discardedReads = vector(mode = "character", length = length(sam_reads)) # Preallocate a vector of discarded reads
  num_reads = 0
  num_aligned = 0
  num_discarded = 0
  num_softclipped = 0
  num_hardclipped = 0
  num_insertion = 0
  num_outsidewindow = 0
  for (i in 1:length(sam_reads)) {
    line_split = strsplit(sam_reads[i], "\t", fixed = T)[[1]]
    num_reads = num_reads + 1
    
    if (length(line_split) <= 10) { next } # avoid the .SAM header
    sam_position = as.integer(line_split[4])
    relative_pos = sam_position - start
    if (abs(relative_pos) > 1000) {
      num_outsidewindow = num_outsidewindow + 1
      stop("Likely ERROR: read position and reference coordinates differ by more than 1000 bp")
    }
    
    sam_cigar = line_split[6]
    sam_read = line_split[10]
    if (sam_cigar == "*") {
      num_discarded = num_discarded + 1
      discardedReads[num_discarded] = sam_reads[i]
      next # Ignore this read, as it doesn't have a valid CIGAR
    }
    
    # computes the read sequence within the coords of interest
    output = registerRead(relative_pos, sam_cigar, sam_read, ref_seq_length)
    if (output$read_ok) {
      #alignedReads = c(alignedReads, output$read)
      num_aligned = num_aligned + 1
      alignedReads[num_aligned] = output$read
    } else {
      num_discarded = num_discarded + 1
      discardedReads[num_discarded] = sam_reads[i]
    }
    if (output$read_status == "hardclipped") {
      num_hardclipped = num_hardclipped + 1
    } else if (output$read_status == "insertion") {
      num_insertion = num_insertion + 1
    } else if (output$read_status == "softclipped") {
      num_softclipped = num_softclipped + 1
    }
    #print "Registered read:"
    #print(alignedReads[i])			
  }
  alignedReads = alignedReads[alignedReads != ""]
  discardedReads = discardedReads[discardedReads != ""]
  
  return(list("alignedReads" = alignedReads, "discardedReads" = discardedReads, "num_reads" = num_reads,
              "num_hardclipped" = num_hardclipped, "num_insertion" = num_insertion,
              "num_softclipped" = num_softclipped, "num_outsidewindow" = num_outsidewindow))
}

getUDPPlot = function(udp.df, plot_title, highlight_site = NA, cut_site = NA) {
  plot_list = list()
  # We allow there to be up to 2 deletions in the UDP
  udp.df = udp.df %>% filter(has_crispr_deletion, !is.na(deletion_start))
  if (nrow(udp.df) == 0) {
    return(ggarrange(textPlot("No UDPs to plot."), top=plot_title, draw = F))
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
    simple_theme + coord_cartesian(xlim=c(1, xmax)) + scale_y_reverse() +
    xlab("Nucleotide position") + ylab("Cumulative UDP count")
  
  count.plot.df = data.frame(x=1:xmax)
  isPositionDel = function(i) {
    sapply(udp.df$udp, FUN=function(x) substr(x, i, i) == '*')
  }
  getDelCount = function(i) {
    sum(isPositionDel(i))
  }
  getDelReadCount = function(i) {
    sum(isPositionDel(i) * udp.df$num_reads)
  }
  
  count.plot.df$udp_count = sapply(count.plot.df$x, FUN=getDelCount)
  count.plot.df$read_count = sapply(count.plot.df$x, FUN=getDelReadCount)
  p.udp_count = ggplot(count.plot.df, aes(x=x, y=udp_count)) +
    geom_bar(stat="identity", fill="dodgerblue2") +
    simple_theme + theme(axis.title.x=element_blank()) + ylab("UDP count") +
    coord_cartesian(xlim=c(1, xmax))
  p.read_count = ggplot(count.plot.df, aes(x=x, y=read_count)) +
    geom_bar(stat="identity", fill="dodgerblue2") +
    simple_theme + theme(axis.title.x=element_blank()) + ylab("Read count") +
    coord_cartesian(xlim=c(1, xmax))
  
  if (!is.na(highlight_site)) {
    p.udp_dist = p.udp_dist + geom_vline(xintercept = highlight_site, color="darkgreen", alpha=0.5)
    p.udp_count = p.udp_count + geom_vline(xintercept = highlight_site, color="darkgreen", alpha=0.5)
    p.read_count = p.read_count + geom_vline(xintercept = highlight_site, color="darkgreen", alpha=0.5)
  }
  if (!is.na(cut_site)) {
    p.udp_dist = p.udp_dist + geom_vline(xintercept = cut_site, color="grey20", linetype = "longdash", alpha=0.5)
    p.udp_count = p.udp_count + geom_vline(xintercept = cut_site, color="grey20", linetype = "longdash", alpha=0.5)
    p.read_count = p.read_count + geom_vline(xintercept = cut_site, color="grey20", linetype = "longdash", alpha=0.5)
  }
  
  p.full = ggarrange(p.read_count, p.udp_count, p.udp_dist, ncol=1, heights=c(1,1,4),
                     top = plot_title, draw = F)
  return(p.full)
}

getDeletionProfilePlot = function(replicate.udp.df, plot_title = NA, highlight_site = NA, cut_site = NA, show_average = T, show_replicates = T) {
  if (!show_average & !show_replicates) {
    stop("getDeletionProfilePlot: One of show_average and show_replicates should be true.")
  }
  if (nrow(replicate.udp.df) == 0) {
    return(ggarrange(textPlot("No UDPs to plot."), top=plot_title, draw = F))
  }
  xmax = nchar(replicate.udp.df$udp[1])
  # Get the deletion profile for each replicate separately
  getDelReadCount = function(i) {
    del_positions = sapply(cur.df$udp, FUN=function(x) substr(x, i, i) == '*')
    sum(del_positions * cur.df$num_reads)
  }
  unique_reps.df = unique(replicate.udp.df %>% dplyr::select(type, replicate))
  counts.list = list()
  for (i in 1:nrow(unique_reps.df)) {
    curType = unique_reps.df[i,]$type
    curReplicate = unique_reps.df[i,]$replicate
    count.df = data.frame(x=1:xmax, type = curType, replicate = curReplicate)
    cur.df = replicate.udp.df %>% dplyr::filter(replicate == curReplicate, type == curType)
    count.df$del_pct = 100 * sapply(count.df$x, FUN=getDelReadCount) / sum(cur.df$num_reads)
    counts.list = c(counts.list, list(count.df))
  }
  # Combine deletion profiles for replicates
  delpct.plot.df = bind_rows(counts.list)
  delpct.plot.df$replicate = as.character(delpct.plot.df$replicate)
  
  if (show_average) {
    # Merge gDNA replicates together (and similarly for cDNA), and
    # determine deletion profiles
    merged.udp.df = summarise(replicate.udp.df %>% group_by(udp, type),
                              num_reads = sum(num_reads)) %>%
      arrange(-num_reads)
    unique_reps.df = unique(merged.udp.df %>% ungroup() %>% dplyr::select(type))
    counts.list = list()
    for (i in 1:nrow(unique_reps.df)) {
      curType = unique_reps.df[i,]$type
      count.df = data.frame(x=1:xmax, type = curType)
      cur.df = merged.udp.df %>% dplyr::filter(type == curType)
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
    scale_color_manual(values = color_values) +
    scale_size_manual(values = size_values) +
    scale_alpha_manual(values = alpha_values)
  
  if (!is.na(plot_title)) {
    p.gDNA_cDNA = p.gDNA_cDNA + ggtitle(plot_title)
  }
  if (!is.na(highlight_site)) {
    p.gDNA_cDNA = p.gDNA_cDNA + geom_vline(xintercept = highlight_site, color="darkgreen", alpha=0.5)
  }
  if (!is.na(cut_site)) {
    p.gDNA_cDNA = p.gDNA_cDNA + geom_vline(xintercept = cut_site, color="grey20", linetype = "longdash", alpha=0.5)
  }
  return(p.gDNA_cDNA)
}

getUNSPlot = function(udp.df, plot_title, min_gDNA_count = 10, min_cDNA_count = 0, max_udps = 80, cut_site = NA) {
  udp_types = unique(udp.df$type)
  if (!("cDNA" %in% udp_types & "gDNA" %in% udp_types)) {
    return(ggarrange(textPlot("Cannot make UNS plot without both cDNA and gDNA replicates"), top=plot_title, draw = F))
  }
  # We want to get a dataframe with a column each for gDNA and cDNA counts
  # for each UDP. But we want to retain the other columns, which are largely
  # duplicated. We do this in a couple of steps.
  # First spread the count columns out
  spread.udp.df = udp.df %>% dplyr::select(udp, type, num_reads) %>%
    spread(type, num_reads, fill = 0)
  # Aggregate the remaining columns by UDP
  grouped.udp.df = udp.df %>% 
    group_by(udp) %>%
    dplyr::summarise(is_hdr_allele = any(is_hdr_allele, na.rm = T),
                     is_wt_allele = any(is_wt_allele, na.rm = T),
                     has_crispr_deletion = any(has_crispr_deletion, na.rm = T),
                     deletion_start = first(deletion_start),
                     deletion_end = first(deletion_end))
  
  spread.udp.df = spread.udp.df %>% dplyr::left_join(grouped.udp.df, by="udp") %>%
    dplyr::rename(gDNA_count = gDNA, cDNA_count = cDNA) %>%
    dplyr::mutate(total_count = gDNA_count + cDNA_count)
  
  # Results not likely to be stable if the WT count is too low
  wtCount = spread.udp.df %>%
    dplyr::filter(is_wt_allele) %>%
    arrange(-gDNA_count) %>%
    group_by(is_wt_allele) %>%
    dplyr::summarise(gDNA_count = sum(gDNA_count),
                     cDNA_count = sum(cDNA_count))
  if (wtCount$gDNA_count < 100) {
    warning(sprintf("%s wild-type gDNA count (%d) is low and may lead to unstable estimates.", plot_title, wtCount$gDNA_count))
    if (wtCount$gDNA_count < 1) {
      stop(sprintf("%s wild-type gDNA count is zero!", plot_title))
    }
  }
  if (wtCount$cDNA_count < 100) {
    warning(sprintf("%s wild-type cDNA count (%d) is low and may lead to unstable estimates.", plot_title, wtCount$cDNA_count))
    if (wtCount$cDNA_count < 1) {
      stop(sprintf("%s wild-type cDNA count is zero!", plot_title))
    }
  }
  
  dels.df = spread.udp.df %>% 
    dplyr::filter(is_hdr_allele | (has_crispr_deletion & gDNA_count >= min_gDNA_count & cDNA_count >= min_cDNA_count))
  dels.df$overlaps_site = F
  if (!is.na(max_udps)) {
    if (nrow(dels.df) > max_udps) {
      dels.df = dels.df %>% arrange(!is_hdr_allele, -total_count) %>% .[1:max_udps,]
    }
  }
  
  dels.df$uns = (dels.df$cDNA_count / wtCount$cDNA_count) / (dels.df$gDNA_count / wtCount$gDNA_count)
  
  # matrix containing the deletion binary code
  M = as.matrix(do.call(rbind, lapply(as.list(dels.df$udp), udpToBinary)))
  M_chars = as.matrix(do.call(rbind, lapply(as.list(dels.df$udp), strToChars)))
  
  model <- hclust(dist(M))
  dhc <- as.dendrogram(model)
  ddata <- dendro_data(dhc, type = "rectangle")
  dendro_span = max(ddata$segments$y, ddata$segments$yend) - min(ddata$segments$y, ddata$segments$yend)
  
  plot.df = as.data.frame(M_chars[model$order,])
  colnames(plot.df) = as.character(c(1:ncol(plot.df)))
  profileSpan = ncol(plot.df)
  plot.df$id = c(1:nrow(plot.df))
  plot.df[,"HDR allele"] = dels.df[model$order,]$is_hdr_allele
  
  plot.gather.df = plot.df %>% tidyr::gather(key = "position", value = "udpchar", 1:profileSpan)
  plot.gather.df$pos = as.numeric(plot.gather.df$position)
  
  p.dendro = ggplot() + 
    geom_segment(aes(x = segment(ddata)$x, y = segment(ddata)$y, xend = segment(ddata)$xend, yend = segment(ddata)$yend)) +
    coord_cartesian(xlim=c(min(plot.gather.df$id), max(plot.gather.df$id)), ylim=c(min(plot.gather.df$pos), max(plot.gather.df$pos))) +
    ylab("") + scale_y_continuous(expand = c(0, 0)) +
    coord_flip() +
    theme_bw() + theme(axis.title.y = element_blank(),
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       axis.line = element_blank(),
                       panel.border = element_blank(),
                       panel.grid = element_blank(),
                       plot.margin = unit(c(0,0,0,0), "cm"))
  
  dot_size = 1.5
  dash_size = 0.5
  if (max(plot.gather.df$pos > 70)) {
    dot_size = 0.8
    dash_size = 0.25
  }
  p.udp = ggplot() + 
    geom_point(aes(y = pos, x = id, colour = `HDR allele`), size = dash_size, shape = 19, data = plot.gather.df[plot.gather.df$udpchar == '-',]) +
    geom_point(aes(y = pos, x = id, colour = `HDR allele`), size = dot_size, shape = '-', alpha = 0.8, data = plot.gather.df[plot.gather.df$udpchar == '*',]) +
    scale_color_manual(values = c("grey20", "blue")) +
    geom_point(aes(y = pos, x = id, shape = udpchar), color = "black", size = 2.5, data = plot.gather.df[!(plot.gather.df$udpchar == '*' | plot.gather.df$udpchar == '-'),]) + scale_shape_identity() +
    coord_cartesian(xlim=c(min(plot.gather.df$id), max(plot.gather.df$id)), ylim=c(min(plot.gather.df$pos), max(plot.gather.df$pos))) +
    ylab("Position") + scale_y_continuous(expand = c(0.01, 0)) + 
    coord_flip() +
    theme_bw() + theme(legend.position = "none",
                       axis.title.y = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       axis.line = element_blank(),
                       panel.border = element_blank(),
                       plot.margin = unit(c(0,0,0,0), "cm"))
  
  if (!is.na(cut_site)) {
    p.udp = p.udp + geom_hline(yintercept = cut_site, color="grey20", linetype = "longdash", alpha=0.5)
  }
  
  # Left side of the plot, showing UDP-normalised scores (UNS)  
  uns.plot.df = data.frame(x = label(ddata)$x,
                           y = log2(dels.df$uns[model$order] + 1),
                           yend = 0,
                           cDNA_lower = dels.df$uns[model$order] < 1,
                           gDNA = dels.df$gDNA_count[model$order],
                           overlaps_site = dels.df$overlaps_site[model$order])
  p.uns = ggplot(uns.plot.df) +
    annotate("segment", x = uns.plot.df$x, xend = uns.plot.df$x, y = uns.plot.df$y, yend = uns.plot.df$yend, colour = 'gray') +
    geom_point(aes(x = x, y = y, size = log10(gDNA), colour = cDNA_lower)) +
    scale_size(range = c(0.5, 5)) + scale_y_reverse() +
    coord_cartesian(xlim=c(min(plot.gather.df$id), max(plot.gather.df$id))) +
    coord_flip() +
    xlab("Unique deletion profile") + ylab("log2(UNS+1)") +
    theme_bw() + theme(legend.position = "left",
                       axis.line = element_blank(),
                       panel.border = element_blank(),
                       plot.margin = unit(c(0,0,0,0), "cm")) +
    scale_color_discrete(guide=F)
  
  p.udp_profile = ggarrange(p.uns, p.udp, p.dendro, nrow=1, ncol=3, widths=c(1,4,1), top = plot_title, draw = F)
  #p.udp_profile = plot_grid(p.uns, p.udp, p.dendro, nrow=1, ncol=3, rel_widths = c(2,4,1))
  return(p.udp_profile)
}

textPlot = function(text, title = NA) {
  ggThemeBlank = theme_bw() + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())
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
  return(paste0(seq_chars, collapse = ""))
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
  out_chars = character(length = length(profile_positions))
  for (i in 1:length(profile_positions)) {
    out_chars[i] = substr(seq, profile_positions[i], profile_positions[i])
  }
  return(paste0(out_chars, collapse = ""))
}

getReadMismatchProfile = function(seq, ref_sequence) {
  #print(sprintf("Seq:%s\t%s", seq, ref_sequence))
  seq_chars = strsplit(seq, "")[[1]]
  ref_chars = strsplit(ref_sequence, "")[[1]]
  getReadCharsMismatchProfile(seq_chars, ref_chars)
}
getReadCharsMismatchProfile = function(seq_chars, ref_chars) {
  output = rep("-", length(ref_chars))
  mismatchpos = (seq_chars != ref_chars)
  output[mismatchpos] = seq_chars[mismatchpos]
  return(paste0(output, collapse=""))
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

isDNALetter = function(s) {
  return(s == 'A' | s == 'C' | s == 'G' | s == 'T')
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
  
  effect = cDNA_ratio / gDNA_ratio
  effect_sd = ratioUncertaintySD(cDNA_ratio, sd_cDNA_ratio, gDNA_ratio, sd_gDNA_ratio, 0)
  
  Zscore = (effect - 1) / effect_sd
  pvalue = 2 * pnorm(-abs(Zscore))
  return(list(effect = effect, effect_sd = effect_sd, pval = pvalue))
}

vecRatioUncertaintySD = function(A, B) {
  ratioUncertaintySD(mean(A), sd(A), mean(B), sd(B), cov(A, B))
}

ratioUncertaintySD = function(meanA, sdA, meanB, sdB, covAB) {
  abs(meanA / meanB) * sqrt((sdA / meanA)^2 + (sdB / meanB)^2 - (2*covAB / (meanA*meanB)))
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

main()
