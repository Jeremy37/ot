#!/usr/bin/env Rscript
suppressMessages(library(optparse))

#source("/Users/jeremys/work/opentargets/src/experiment/genie.functions.R")
source("/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/src/experiment/genie.functions.R")


main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
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
    make_option(c("--qc_plot_max_udps"), type="integer", default=20, metavar="INT", help="[default %default]"),
    make_option(c("--qc_plot_min_udp_fraction"), type="integer", default=0.005, metavar="FLOAT", help="[default %default]"),
    make_option(c("--qc_plot_exclude_wt"), type="logical", default=T, action="store_true", help="[default %default]"),
    make_option(c("--del_span_start"), type="integer", default=NULL, help="[default %default]"),
    make_option(c("--del_span_end"), type="integer", default=NULL, help="[default %default]"),
    make_option(c("--uns_plot_min_gDNA"), type="integer", default=10, metavar="INT", help="[default %default]"),
    make_option(c("--uns_plot_min_cDNA"), type="integer", default=0, metavar="INT", help="[default %default]"),
    make_option(c("--uns_plot_max_udps"), type="integer", default=40, metavar="INT", help="[default %default]"),
    make_option(c("--no_allele_profile"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--no_site_profile"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--no_udp_profile"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--no_uns_profile"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--no_stats"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--replicate_qc_plots"), type="logical", default=T, action="store_true", help="[default %default]"),
    make_option(c("--no_summary_plots"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--no_replicate_plots"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--variance_analysis"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--variance_analysis_min_count"), type="integer", default=100, help="[default %default]"),
    make_option(c("--variance_analysis_min_fraction"), type="numeric", default=0.001, help="[default %default]"),
    make_option(c("--power_analysis"), type="logical", default=F, action="store_true", help="[default %default]"),
    make_option(c("--plot_width"), type="integer", default=8, metavar="INT", help="[default %default]"),
    make_option(c("--plot_height"), type="integer", default=7, metavar="INT", help="[default %default]"),
    make_option(c("--read_data"), type="character", default=NULL)
  )
  
  parser = OptionParser(usage = "paired.deletion.analysis.R --regions file.tsv --replicates file.tsv --out output_path [options]",
                        option_list = option_list)
  opt <<- parse_args(parser)

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
  
  runGenIE(opt)
}


###########################################################################

system.time( main() )
