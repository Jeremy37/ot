#!/usr/bin/env Rscript

# This script reads in a VCF file and an "ASECounts" file. The VCF file has the genotypes
# for a single individual in the region of 1 gene. The ASECounts file has allele-specific
# read counts for possibly multiple samples that come from the individual and which cover
# the gene. The script performs statistical tests for allele-specific expression, and
# produces some plots.
library(tidyverse)
library(stringr)
library(egg)
library(aod)
options(stringsAsFactors = F)
myargs <- NULL

main = function()
{
  myargs <<- getArgs()
  myargs <<- list(input = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1/analysis/udp_replicates.tsv",
                  out = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1/analysis/batch1",
                  max_mismatch_frac = 0.05,
                  editing_site_zoom_window = 30, # window around the HDR site for displaying a zoom of the editing region
                  min_read_region_overlap = 30, # minimum number of read bases correctly aligned within the region of interest
                  plot_per_locus = F,
                  exclude_nonspanning = F)
  print(myargs)
  
  if (is.na(myargs$exclude_nonspanning)) { myargs$exclude_nonspanning <<- F }
  
  input.df = readr::read_tsv(myargs$input, col_types="cicciiiccc")
  regions = unique(input.df$name)
  all_region_plots = list()
  region_index = 1
  for (region in regions) {
    result = doFullDeletionAnalysis(input.df %>% dplyr::filter(name == region))
    all_region_plots[[region_index]] = result$plot_list
    if (myargs$plot_per_locus) {
      fname = sprintf("%s.%s.plots.pdf", myargs$out, input.df[i,]$name)
      pdf(file = fname, width = 8, height = 7)
      print(result$plot_list)
      dev.off()
    }
    region_index = region_index + 1
  }
  
  if (!myargs$plot_per_locus) {
    fname = sprintf("%s.all_plots.pdf", myargs$out)
    pdf(file = fname, width = 8, height = 7)
    print(all_region_plots)
    dev.off()
    
    fname = sprintf("%s.summary_plots.pdf", myargs$out)
    pdf(file = fname, width = 8, height = 7)
    for (region_plots in all_region_plots) {
      print(region_plots$stat_plots$stats1)
      print(region_plots$merged_plots[[1]])
      print(region_plots$merged_plots[[2]])
    }
    dev.off()
  }
  
}


doFullDeletionAnalysis = function(input.df)
{
  replicate_list = list()
  replicate_plots = list()
  stats_list = list()
  for (i in 1:nrow(input.df)) {
    result = doReplicateDeletionAnalysis(input.df[i,]$name,
                                         input.df[i,]$replicate,
                                         input.df[i,]$cigars_file,
                                         input.df[i,]$chr,
                                         input.df[i,]$start,
                                         input.df[i,]$end,
                                         input.df[i,]$hdr_site,
                                         input.df[i,]$hdr_allele,
                                         input.df[i,]$wt_allele,
                                         input.df[i,]$ref_sequence)
    replicate_list[[i]] = result$udp.df
    replicate_plots[[i]] = result$plots
    stats_list[[i]] = result$stats
    
    fname = sprintf("%s.%s_%d.udps.tsv", myargs$out, input.df[i,]$name, input.df[i,]$replicate)
    print(sprintf("Found %d unique deletion profiles. Saving to file %s.", nrow(result$udp.df), fname))
    write.table(result$udp.df, fname, quote=F, row.names=F, col.names=T, sep="\t")
    
    fname = sprintf("%s.%s_%d.mismatch_profile.tsv", myargs$out, input.df[i,]$name, input.df[i,]$replicate)
    write.table(result$wt_hdf.df, fname, quote=F, row.names=F, col.names=T, sep="\t")
  }
  
  fname = sprintf("%s.%s_%d.stats.tsv", myargs$out, input.df[i,]$name, input.df[i,]$replicate)
  stats.df = bind_rows(stats_list)
  write.table(stats.df, fname, quote=F, row.names=F, col.names=T, sep="\t")
  
  # Analyse all input replicates together
  test.df = data.frame(name = c(paste(stats.df$name, stats.df$replicate, "gDNA", sep="_"), paste(stats.df$name, stats.df$replicate, "cDNA", sep="_")),
                       is_cDNA = c(0,0,0, 1,1,1),
                       wt  = c(stats.df$num_wt_reads_gDNA, stats.df$num_wt_reads_cDNA),
                       mut = c(stats.df$num_hdr_reads_gDNA, stats.df$num_hdr_reads_cDNA))
  fm1 = betabin(cbind(mut, wt) ~ is_cDNA, ~ 1, data = test.df[,2:4])
  
  stats.df$HDR_rate_gDNA = stats.df$num_hdr_reads_gDNA / stats.df$num_gDNA_reads
  stats.df$editing_rate_gDNA = (stats.df$num_edit_reads_gDNA + stats.df$num_hdr_reads_gDNA) / stats.df$num_gDNA_reads
  
  stats.df$cDNA_HDR_WT_ratio = (stats.df$num_hdr_reads_cDNA / stats.df$num_wt_reads_cDNA)
  stats.df$gDNA_HDR_WT_ratio = (stats.df$num_hdr_reads_gDNA / stats.df$num_wt_reads_gDNA)
  stats.df$cDNA_gDNA_HDR_ratio = stats.df$cDNA_HDR_WT_ratio / stats.df$gDNA_HDR_WT_ratio
  stats.plot.df = stats.df %>% dplyr::select(name, replicate, num_udps, cDNA_gDNA_HDR_ratio, cDNA_HDR_WT_ratio, gDNA_HDR_WT_ratio, HDR_rate_gDNA,
                                             editing_rate_gDNA, num_gDNA_reads, num_cDNA_reads,
                                             num_hdr_reads_gDNA, num_hdr_reads_cDNA, num_wt_reads_gDNA, num_wt_reads_cDNA,
                                             reads_excluded_for_minoverlap, reads_excluded_for_mismatches, reads_excluded_for_multiple_deletions)
  stats.plot.df$cDNA_HDR_WT_ratio = sprintf("%.3f", stats.plot.df$cDNA_HDR_WT_ratio)
  stats.plot.df$gDNA_HDR_WT_ratio = sprintf("%.3f", stats.plot.df$gDNA_HDR_WT_ratio)
  stats.plot.df$cDNA_gDNA_HDR_ratio = sprintf("%.3f", stats.plot.df$cDNA_gDNA_HDR_ratio)
  stats.plot.df$HDR_rate_gDNA = sprintf("%.3f%%", 100 * stats.plot.df$HDR_rate_gDNA)
  stats.plot.df$editing_rate_gDNA = sprintf("%.3f%%", 100 * stats.plot.df$editing_rate_gDNA)
  ggThemeBlank = theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                                    axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())
  myTableTheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.8)),
                                            colhead = list(fg_params=list(cex = 0.8)),
                                            rowhead = list(fg_params=list(cex = 0.8)))
  p.stats1 = ggplot(data.frame(x=1:10, y=1:10), aes(x,y)) + geom_blank() + ggThemeBlank +
    annotation_custom(tableGrob(t(stats.plot.df), theme = myTableTheme), xmin=1, xmax=9, ymin=1, ymax=9) +
    annotate("text", x=5, y=10, label = sprintf("Region %s summary", stats.df$name[1]), vjust = 1, fontface = 2, size = 5) +
    annotate("text", x=5, y=9, label = sprintf("Difference in HDR ratio, p = %.3g", summary(fm1)@Coef$`Pr(> |z|)`[2]), vjust = 1, fontface = 2, size = 4)
  p.stats1
  
  mainStr = sprintf("Estimate ('is_cDNA') for whether HDR / WT reads\nacross replicates differs for cDNA vs. gDNA",
                    summary(fm1)@Coef$`Pr(> |z|)`[2])
  summaryStr = paste(capture.output(summary(fm1)), collapse = "\n")
  p.stats2 = ggplot(data.frame(x=c(0,1), y=c(0,1)), aes(x=x, y=y)) + theme_bw() + ggtitle(sprintf("Region %s stats", stats.df$name[1])) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    annotate("text", x=0, y=1, label = mainStr, hjust = 0, vjust = 1, fontface = 2) +
    annotate("text", x=0, y=0.85, label = summaryStr, hjust = 0, vjust = 1) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())
  
  # Merge together results from all replicates and generate new plots
  concat.df = bind_rows(replicate_list)
  merged.udp.df = summarise(concat.df %>% group_by(udp),
            gDNA_count = sum(gDNA_count),
            cDNA_count = sum(cDNA_count),
            total_count = sum(gDNA_count, cDNA_count),
            is_hdr_allele = first(is_hdr_allele),
            is_wt_allele = first(is_wt_allele),
            has_deletion = first(has_deletion),
            deletion_start = first(deletion_start),
            deletion_end = first(deletion_end),
            deletion_length = first(deletion_length),
            avg_seq_length = mean(avg_seq_length),
            avg_mismatch_count = mean(avg_mismatch_count)) %>%
    arrange(-total_count)
  merged_plots_list = getEditingPlots(merged.udp.df, sprintf("Region %s merged replicates", input.df[1,]$name), input.df[1,]$hdr_site - input.df[1,]$start)
  
  plot_list = list(stat_plots=list(stats1=p.stats1, stats2=p.stats2), replicate_plots=replicate_plots, merged_plots=merged_plots_list)
  
  return(list(replicate_list = replicate_list, plot_list = plot_list, stats_list = stats_list))
}


doReplicateDeletionAnalysis = function(name, replicate, cigars_file, chr, start, end, hdr_chr_site, hdr_allele, wt_allele, ref_sequence)
{
  print(sprintf("\n\nAnalysing region %s, replicate %d, file %s, %s:%d-%d, site %d, HDR allele %s, WT allele %s",
                name, replicate, cigars_file, chr, start, end, hdr_chr_site, hdr_allele, wt_allele))
  hdr_site = hdr_chr_site - start
  stats = list(name = name, replicate = replicate)
  
  region.df = readr::read_tsv(cigars_file)
  stats$num_cigars = nrow(region.df)
  stats$num_gDNA_reads = sum(region.df$gDNA_count)
  stats$num_cDNA_reads = sum(region.df$cDNA_count)
  stats$num_reads = stats$num_gDNA_reads + stats$num_cDNA_reads
  
  region.df$total_count = sapply(1:nrow(region.df), FUN=function(i) region.df[i,]$gDNA_count + region.df[i,]$cDNA_count)
  
  region.df$seq_length = sapply(region.df$region_read, FUN=function(s) nchar(gsub("[-*]", "", s)))
  region.df$has_deletion = sapply(region.df$region_read, FUN=function(s) grepl("[*]", s))
  # We only consider a deletion as an "editing deletion" if it covers a region near the HDR site
  del_check_start = max(hdr_site - 30, 1)
  del_check_end = min(hdr_site + 30, nchar(region.df$region_read[1]))
  region.df$has_edit = sapply(substr(region.df$region_read, del_check_start, del_check_end),
                                  FUN=function(s) grepl("[*]", s))
  region.df$has_multiple_deletions = sapply(region.df$region_read, FUN=function(s) grepl("[*]+[^*]+[*]+", s))
  region.df$spanning_read = T
  
  region.df$is_hdr_allele = F
  if (!is.na(hdr_site)) {
    if (is.na(hdr_allele)) {
      stop("If --hdr_site is specified, the --hdr_allele should also be specified")
    }
    region.df$is_hdr_allele = sapply(region.df$region_read, FUN=function(s) {substr(s, hdr_site, hdr_site) == hdr_allele})
    region.df$is_hdr_allele = region.df$is_hdr_allele & !region.df$has_edit
    # We don't want to count the HDR site itself as a mismatch
    if (nchar(hdr_allele) == 1) {
      region.df$mismatch_count[region.df$is_hdr_allele] = region.df$mismatch_count[region.df$is_hdr_allele] - 1
    }
    stats$num_hdr_reads_gDNA = sum(region.df$gDNA_count[region.df$is_hdr_allele])
    stats$num_hdr_reads_cDNA = sum(region.df$cDNA_count[region.df$is_hdr_allele])

    region.df$is_wt_allele = sapply(region.df$region_read, FUN=function(s) {substr(s, hdr_site, hdr_site) == wt_allele})
    region.df$is_wt_allele = region.df$is_wt_allele & !region.df$has_edit
    
    region.df$spanning_read = sapply(region.df$region_read, FUN=function(s) (substring(s, hdr_site, hdr_site) != '-'))
    if (myargs$exclude_nonspanning) {
      stats$nonspanning_reads = sum(region.df$total_count[!region.df$spanning_read])
      stats$nonspanning_cigars = sum(!region.df$spanning_read)
      print(sprintf("%d of %d cigars (%d of %d reads, %.2f%%) excluded due to not spanning the site of interest (position %d)",
                    stats$nonspanning_cigars, stats$num_cigars,
                    stats$nonspanning_reads, stats$num_reads, 100.0 * stats$nonspanning_reads / stats$num_reads,
                    hdr_site))
      region.df = region.df[region.df$spanning_read, ]
    }
  } else {
    region.df$is_wt_allele = !region.df$has_edit
  }
  # make sure we only count as WT those reads which actually cover the HDR site
  region.df$is_wt_allele[!region.df$spanning_read] = NA
  stats$num_wt_reads_gDNA = sum(region.df$gDNA_count[region.df$is_wt_allele], na.rm = T)
  stats$num_wt_reads_cDNA = sum(region.df$cDNA_count[region.df$is_wt_allele], na.rm = T)
  
  # Exclude reads that don't cover enough of the region of interest
  exclude_for_overlap = (region.df$seq_length < myargs$min_read_region_overlap)
  stats$reads_excluded_for_minoverlap = sum(region.df$total_count[exclude_for_overlap])
  stats$cigars_excluded_for_minoverlap = sum(exclude_for_overlap)
  print(sprintf("%d of %d cigars (%d of %d reads, %.2f%%) excluded due to aligning to less than %d bp in the region of interest",
                stats$cigars_excluded_for_minoverlap, stats$num_cigars,
                stats$reads_excluded_for_minoverlap, stats$num_reads, 100.0 * stats$reads_excluded_for_minoverlap / stats$num_reads,
                myargs$min_read_region_overlap))
  region.df = region.df[!exclude_for_overlap, ]
  
  exclude_for_mismatches = region.df$mismatch_count / region.df$seq_length > myargs$max_mismatch_frac
  stats$reads_excluded_for_mismatches = sum(region.df$total_count[exclude_for_mismatches])
  stats$cigars_excluded_for_mismatches = sum(exclude_for_mismatches)
  print(sprintf("%d of %d cigars (%d of %d reads, %.2f%%) excluded due to having more than %f%% of the read being mismatches",
                stats$cigars_excluded_for_mismatches, stats$num_cigars,
                stats$reads_excluded_for_mismatches, stats$num_reads, 100.0 * stats$reads_excluded_for_mismatches / stats$num_reads,
                myargs$max_mismatch_frac * 100))
  region.df = region.df[!exclude_for_mismatches, ]
  
  stats$num_edit_reads_gDNA = sum(region.df$gDNA_count[region.df$has_edit], na.rm = T)
  stats$num_edit_reads_cDNA = sum(region.df$cDNA_count[region.df$has_edit], na.rm = T)
  
  stats$reads_excluded_for_multiple_deletions = sum(region.df$total_count[region.df$has_multiple_deletions])
  stats$cigars_excluded_for_multiple_deletions = sum(region.df$has_multiple_deletions)
  print(sprintf("%d of %d cigars (%d of %d reads, %.2f%%) excluded due to having more than one separate deletion",
                stats$cigars_excluded_for_multiple_deletions, stats$num_cigars,
                stats$reads_excluded_for_multiple_deletions, stats$num_reads, 100.0 * stats$reads_excluded_for_multiple_deletions / stats$num_reads))
  region.df = region.df[!region.df$has_multiple_deletions, ]
  
  # Convert reads to unique deletion profiles (UDP), and merge identical UDPs
  region.df$udp = sapply(region.df$region_read, FUN=getReadUDP, hdr_site)
  
  udp.df = summarise(region.df %>% group_by(udp),
                     gDNA_count = sum(gDNA_count),
                     cDNA_count = sum(cDNA_count),
                     total_count = sum(gDNA_count, cDNA_count),
                     is_hdr_allele = first(is_hdr_allele),
                     is_wt_allele = first(is_wt_allele),
                     has_deletion = first(has_deletion),
                     num_cigars = n(),
                     avg_seq_length = mean(seq_length),
                     avg_mismatch_count = mean(mismatch_count),
                     cigars = paste(cigar, collapse=",")) %>%
    arrange(-total_count)
  
  stats$num_udps = nrow(udp.df)
  
  # Order UDPs by the deletion start position and plot
  getDeletionLoc = function(i) {
    if (udp.df$has_deletion[i]) {
      return(sapply(udp.df$udp[i], FUN=function(udp) str_locate(udp, "[\\*]+")))
    }
    return(NA)
  }
  udp_dels = sapply(1:nrow(udp.df), FUN=getDeletionLoc)
  udp.df$deletion_start = sapply(udp_dels, FUN=function(x) x[1])
  udp.df$deletion_end = sapply(udp_dels, FUN=function(x) x[2]+1)
  udp.df$deletion_length = sapply(udp_dels, FUN=function(x) x[2]-x[1]+1)
  
  plot_list = getEditingPlots(udp.df, sprintf("Region %s replicate %d", name, replicate), hdr_site)
  
  # Make a table which has just the different versions of the WT and HDR alleles
  wt_hdr.df = region.df %>% dplyr::filter(is_wt_allele | is_hdr_allele)
  wt_hdr.df$mismatch_profile = sapply(wt_hdr.df$region_read, FUN=getReadMismatchProfile, ref_sequence)
  wt_hdr.df = wt_hdr.df %>% dplyr::select(mismatch_profile, region_read, cigar, everything())

  return(list(udp.df = udp.df, wt_hdf.df = wt_hdr.df, plots = plot_list, stats = stats))
}


getEditingPlots = function(udp.df, plot_title, hdr_site) {
  plot_list = list()
  udp.plot.df = udp.df %>% filter(!is.na(deletion_start)) %>% dplyr::select(deletion_start, deletion_end, deletion_length)
  udp.plot.df = udp.plot.df %>% arrange(deletion_start)
  udp.plot.df$y = 1:nrow(udp.plot.df)
  xmax = nchar(udp.df$udp[1])
  simple_theme = theme_bw() + theme(panel.grid.minor = element_line(colour="grey90", size=0.1), panel.grid.major = element_line(colour="grey90", size=0.1))
  p.udp_dist = ggplot(udp.plot.df) +
    geom_segment(aes(x = deletion_start, xend = deletion_end, y = y, yend = y), color = "dodgerblue3") +
    geom_vline(xintercept = hdr_site, color="firebrick", alpha=0.5) +
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
    sum(isPositionDel(i) * udp.df$total_count)
  }
  getDelReadCountgDNA = function(i) {
    sum(isPositionDel(i) * udp.df$gDNA_count)
  }
  getDelReadCountcDNA = function(i) {
    sum(isPositionDel(i) * udp.df$cDNA_count)
  }
  
  count.plot.df$udp_count = sapply(count.plot.df$x, FUN=getDelCount)
  count.plot.df$read_count = sapply(count.plot.df$x, FUN=getDelReadCount)
  p.udp_count = ggplot(count.plot.df, aes(x=x, y=udp_count)) +
    geom_bar(stat="identity", color="dodgerblue2") +
    geom_vline(xintercept = hdr_site, color="firebrick", alpha=0.5) +
    simple_theme + theme(axis.title.x=element_blank()) + ylab("UDP count") +
    coord_cartesian(xlim=c(1, xmax))
  p.read_count = ggplot(count.plot.df, aes(x=x, y=read_count)) +
    geom_bar(stat="identity", color="dodgerblue2") +
    geom_vline(xintercept = hdr_site, color="firebrick", alpha=0.5) +
    simple_theme + theme(axis.title.x=element_blank()) + ylab("Read count") +
    coord_cartesian(xlim=c(1, xmax))
  
  p.full = ggarrange(p.read_count, p.udp_count, p.udp_dist, ncol=1, heights=c(1,1,4),
                     top = plot_title, draw = F)
  plot_list[["full"]] = p.full
  #print(p.full)
  
  count.plot.df$gDNA = 100 * sapply(count.plot.df$x, FUN=getDelReadCountgDNA) / sum(udp.df$gDNA_count)
  count.plot.df$cDNA = 100 * sapply(count.plot.df$x, FUN=getDelReadCountcDNA) / sum(udp.df$cDNA_count)
  comparison.plot.df = count.plot.df %>% tidyr::gather(key="type", value="count", gDNA:cDNA)
  # Determine which is higher (gDNA or cDNA)
  #gDNA_higher = count.plot.df$gDNA[hdr_site] > count.plot.df$cDNA[hdr_site]
  p.gDNA_cDNA = ggplot(comparison.plot.df, aes(x=x, y=count, color=type, fill=type)) + geom_line() +
    geom_vline(xintercept = hdr_site, color="firebrick", alpha=0.5) +
    simple_theme + xlab("Nucleotide position") + ylab("Deletion frequency (%)")
  
  if (hdr_site) {
    # Zoom in to get a better view of the region near the HDR site
    w = myargs$editing_site_zoom_window
    comparison.plot.zoom.df = comparison.plot.df %>% dplyr::filter(x >= hdr_site - w & x <= hdr_site + w)
    p.gDNA_cDNA.zoom = ggplot(comparison.plot.zoom.df, aes(x=x, y=count, color=type, fill=type)) + geom_line() +
      simple_theme + geom_vline(xintercept = hdr_site, color="firebrick", alpha=0.5) +
      theme(axis.title.x=element_blank()) + ylab("Deletion frequency (%)")
    
    p.gDNA_cDNA = ggarrange(p.gDNA_cDNA.zoom, p.gDNA_cDNA, ncol=1, heights=c(1,1),
                            top = plot_title, draw = F)
  }
  plot_list[["gDNA_cDNA"]] = p.gDNA_cDNA
  return(plot_list)
}


getReadUDP = function(seq, hdr_site) {
  paste0(gsub("[^*]", "-", substr(seq, 1, hdr_site-1)),
         substr(seq, hdr_site, hdr_site),
         gsub("[^*]", "-", substr(seq, hdr_site+1, nchar(seq))))
}


getReadMismatchProfile = function(seq, ref_sequence) {
  #print(sprintf("Seq:%s\t%s", seq, ref_sequence))
  seq_chars = strsplit(seq, "")[[1]]
  ref_chars = strsplit(ref_sequence, "")[[1]]
  output = rep("-", length(ref_chars))
  output[seq_chars != ref_chars] = seq_chars[seq_chars != ref_chars]
  return(paste0(output, collapse=""))
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
