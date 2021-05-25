source("/Users/jeremys/work/opentargets/src/experiment/genie.functions.R")


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
  min_window_overlap = 50, # minimum number of read bases correctly aligned within the region of interest
  allele_match_window = 5,
  exclude_multiple_deletions = F,
  exclude_nonspanning_reads = T,
  exclude_nonspanning_deletions = T,
  ratio_to_total_reads = F,
  use_cdna_dels_only = F,
  qc_plot_max_udps = 10,
  qc_plot_min_udp_fraction = 0.002,
  qc_plot_exclude_wt = T,
  qc_outlier_threshold = 0,
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
  show_udp_sharing = F,
  replicate_qc_plots = T,
  no_replicate_plots = T,
  variance_analysis = F,
  power_analysis = F,
  plot_width = 8,
  plot_height = 5,
  variance_analysis_min_count = 100,
  variance_analysis_min_fraction = 0.001,
  variance_analysis_split_by_fraction = F)

###############################################################################
opts <<- default_opts
setwd("/Users/jeremys/work/opentargets/experiment/transcribed/power_analysis")

opts$regions = "/Users/jeremys/work/opentargets/experiment/transcribed/power_analysis/power_analysis.regions.tsv"
opts$replicates = "/Users/jeremys/work/opentargets/experiment/transcribed/power_analysis/power_analysis.meta.tsv"
opts$out = "/Users/jeremys/work/opentargets/experiment/transcribed/power_analysis/analysis/power_analysis.2"
opts$read_data = "analysis/power_analysis.read_data.rds"
opts$power_analysis = T
system.time(runGenIE(opts))

opts$regions = "/Users/jeremys/work/opentargets/experiment/transcribed/power_analysis/power_analysis.regions.1.tsv"
opts$replicates = "/Users/jeremys/work/opentargets/experiment/transcribed/power_analysis/power_analysis.meta.tsv"
opts$out = "/Users/jeremys/work/opentargets/experiment/transcribed/power_analysis/analysis/power_analysis.1"
opts$read_data = "analysis/power_analysis.read_data.rds"
opts$power_analysis = T
system.time(runGenIE(opts))

opts$regions = "/Users/jeremys/work/opentargets/experiment/transcribed/power_analysis/power_analysis.regions.ABHD4.tsv"
opts$replicates = "/Users/jeremys/work/opentargets/experiment/transcribed/power_analysis/power_analysis.meta.tsv"
opts$out = "/Users/jeremys/work/opentargets/experiment/transcribed/power_analysis/analysis/power_analysis.ABHD4.2"
opts$read_data = "analysis/power_analysis.read_data.rds"
opts$power_analysis = T
opts$grep_analysis = F
opts$plot_width = 7
opts$plot_height = 4.5
system.time(runGenIE(opts))


###############################################################################
opts <<- default_opts
setwd("/Users/jeremys/work/opentargets/experiment/transcribed/batch1_variance")

opts$regions = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1_variance/batch1_variance.regions.no_SORL1.tsv"
opts$replicates = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1_variance/batch1_variance.meta.newsamples.all.tsv"
opts$out = "/Users/jeremys/work/opentargets/experiment/transcribed/batch1_variance/analysis/batch1_variance.newsamples.all"
opts$read_data = "analysis/batch1_variance.newsamples.all.read_data.rds"
opts$variance_analysis = T
opts$plot_width = 8
opts$plot_height = 5.5
system.time(runGenIE(opts))


###############################################################################
# Editing CLU and CCDC6 SNPs in BOB cell line
opts <<- default_opts
setwd("/Users/jeremys/work/opentargets/experiment/transcribed/BOB_CLU_CCDC6")

opts$regions = "/Users/jeremys/work/opentargets/experiment/transcribed/BOB_CLU_CCDC6/BOB_CLU_CCDC6.regions.tsv"
opts$replicates = "/Users/jeremys/work/opentargets/experiment/transcribed/BOB_CLU_CCDC6/BOB_CLU_CCDC6.genie.meta.tsv"
opts$out = "/Users/jeremys/work/opentargets/experiment/transcribed/BOB_CLU_CCDC6/analysis/BOB_CLU_CCDC6.new2"
opts$read_data = "analysis/BOB_CLU_CCDC6.read_data2.rds"
opts$ratio_to_total_reads = F
opts$viewing_window = 40
opts$editing_window = 20
system.time(runGenIE(opts))


###############################################################################
opts <<- default_opts
setwd("/Users/jeremys/work/opentargets/experiment/transcribed/batch3")

#opts$regions = "/Users/jeremys/work/opentargets/experiment/transcribed/batch3/batch3.genie.regions.tsv"
opts$regions = "/Users/jeremys/work/opentargets/experiment/transcribed/batch3/batch3.genie.regions.TAF1C.tsv"
opts$replicates = "/Users/jeremys/work/opentargets/experiment/transcribed/batch3/batch3.genie.meta.tsv"
opts$out = "/Users/jeremys/work/opentargets/experiment/transcribed/batch3/analysis/batch3.20bp_edit_window.TAF1C"
opts$read_data = "analysis/batch3.read_data.rds"
opts$plot_width = 7
opts$plot_height = 4.5
opts$deletion_analysis = T
system.time(runGenIE(opts))

# SDF4
opts$regions = "/Users/jeremys/work/opentargets/experiment/transcribed/batch3/batch3.genie.regions.exon_primers.tsv"
opts$replicates = "/Users/jeremys/work/opentargets/experiment/transcribed/batch3/batch3.genie.meta.exon_primers.tsv"
opts$max_mismatch_frac = 0.62
opts$out = "/Users/jeremys/work/opentargets/experiment/transcribed/batch3/analysis/batch3.SDF4.exonic"
opts$read_data = "analysis/batch3.read_data.SDF4.rds"
system.time(runGenIE(opts))


###############################################################################
opts <<- default_opts
setwd("/Users/jeremys/work/opentargets/experiment/transcribed/batch4_MUL1_CD33")

opts$regions = "/Users/jeremys/work/opentargets/experiment/transcribed/batch4_MUL1_CD33/batch4_MUL1.genie.regions.tsv"
opts$replicates = "/Users/jeremys/work/opentargets/experiment/transcribed/batch4_MUL1_CD33/batch4_MUL1_CD33.genie.meta.tsv"
opts$out = "/Users/jeremys/work/opentargets/experiment/transcribed/batch4_MUL1_CD33/analysis/batch4_MUL1.20bp_edit_window"
opts$read_data = "analysis/batch4_MUL1_CD33.read_data.rds"
opts$deletion_analysis = T
opts$uns_plot_max_udps = 30
opts$plot_width = 6.8
opts$plot_height = 4.3
system.time(runGenIE(opts))

# Slightly shorter PDF size for UNS plot
opts$plot_width = 7
opts$plot_height = 3.4
opts$out = "/Users/jeremys/work/opentargets/experiment/transcribed/batch4_MUL1_CD33/analysis/batch4_MUL1.20bp_edit_window.2"
system.time(runGenIE(opts))

# Slightly larger PDF size for QC portion of figure
opts$out = "/Users/jeremys/work/opentargets/experiment/transcribed/batch4_MUL1_CD33/analysis/batch4_MUL1.20bp_edit_window.3"
opts$plot_width = 7.4
opts$plot_height = 5
system.time(runGenIE(opts))


###############################################################################
opts <<- default_opts
setwd("/Users/jeremys/work/opentargets/experiment/transcribed/CCDC6_satmut_genie")

opts$regions = "/Users/jeremys/work/opentargets/experiment/transcribed/CCDC6_satmut_genie/ccdc6.satmut.regions.rs1171830.tsv"
opts$replicates = "/Users/jeremys/work/opentargets/experiment/transcribed/CCDC6_satmut_genie/ccdc6.satmut.genie.meta.by_batch.rs1171830.tsv"
opts$out = "/Users/jeremys/work/opentargets/experiment/transcribed/CCDC6_satmut_genie/analysis/CCDC6_satmut_genie.HDR.by_batch.rs1171830"
opts$read_data = "analysis/CCDC6_satmut_genie.indiv.read_data.by_batch.rs1171830.rds"
opts$deletion_analysis = T
opts$plot_width = 7
opts$plot_height = 4.5
system.time(runGenIE(opts))


