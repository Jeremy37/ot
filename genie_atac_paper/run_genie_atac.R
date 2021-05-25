suppressMessages(library(tidyverse))
library(broom)

setwd("/Users/jeremys/work/opentargets/R/rgenie/")
#devtools::document()
devtools::build(vignettes = F)
devtools::install()
#devtools::build_vignettes()
library(rgenie)

dir = "/Users/jeremys/work/opentargets/experiment/atac/"

###############################################################################
# genie_atac_hets
setwd(file.path(dir, "genie_atac_hets"))

genie_het_regions = readr::read_tsv("genie_atac_hets.genie_regions.tsv")
genie_het_replicates = readr::read_tsv("genie_atac_hets.genie_replicates.tsv")

grep_results = grep_analysis(genie_het_regions,
                             genie_het_replicates,
                             required_match_left = 6,
                             required_match_right = 6,
                             min_mapq = 0,
                             quiet = F)

#grep_summary_plot(grep_results[[1]])
write_results(grep_results, base_path = paste0(dir, "genie_atac_hets/analysis/genie_atac_hets.grep_results"))

aln_results = alignment_analysis(genie_het_regions,
                                 genie_het_replicates,
                                 required_match_left = 6,
                                 required_match_right = 6,
                                 crispr_del_window = 100,
                                 min_mapq = 0,
                                 max_mismatch_frac = 0.05,
                                 min_aligned_bases = 50,
                                 exclude_multiple_deletions = FALSE,
                                 exclude_nonspanning_reads = TRUE,
                                 allele_profile = FALSE,
                                 del_span_start = -20,
                                 del_span_end = 20,
                                 quiet = F)

#experiment_summary_plot(grep_results, aln_results)
write_results(aln_results, base_path = paste0(dir, "genie_atac_hets/analysis/genie_atac_hets.alignment_results"))

res$replicate_stats$hdr_frac = res$replicate_stats$num_hdr_reads / (res$replicate_stats$num_wt_reads + res$replicate_stats$num_hdr_reads)
res$replicate_stats$ratio = res$replicate_stats$num_hdr_reads / res$replicate_stats$num_wt_reads
#theme_set(theme_)
ggplot(res$replicate_stats, aes(x=name, y=ratio, fill=type)) +
  geom_boxplot() +
  geom_jitter(alpha=0.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

aln_plots = alignment_analysis_plots(aln_results[[1]],
                                     opts = genie_plot_options(),
                                     variance_components_plot = FALSE,
                                     power_plots = FALSE)


pdf("analysis/genie_atac_hets.pdf", width=7, height=8)
experiment_summary_plot(grep_results, aln_results)

for (i in 1:length(grep_results)) {
  grep_plot = grep_summary_plot(grep_results[[i]])

  aln_plots = alignment_analysis_plots(aln_results[[i]],
                                       opts = genie_plot_options(),
                                       variance_components_plot = FALSE,
                                       power_plots = FALSE)
  print(grep_plot)
  print(aln_plots)
}

dev.off()


###############################################################################
# genie_atac_hets - tagmentation comparison
# Here we compare different options for how gDNA is measured - either amplicon
# PCR, tagmentation forward, tagmentation reverse. We don't really care about
# the ATAC, only the gDNA, so we just use the same ATAC samples to get it
# to run.
setwd(file.path(dir, "genie_atac_hets"))

genie_het_regions = readr::read_tsv("genie_atac_tagment_comparison.regions.tsv")
genie_het_replicates = readr::read_tsv("genie_atac_tagment_comparison.replicates.tsv")

grep_results = grep_analysis(genie_het_regions,
                             genie_het_replicates,
                             required_match_left = 6,
                             required_match_right = 6,
                             min_mapq = 0,
                             quiet = F)

#grep_summary_plot(grep_results[[1]])
write_results(grep_results, base_path = paste0(dir, "genie_atac_hets/analysis/genie_atac_tagment_comparison.grep_results"))
saveRDS(grep_results, paste0(dir, "genie_atac_hets/analysis/genie_atac_tagment_comparison.grep_results.rds"))

aln_results = alignment_analysis(genie_het_regions,
                                 genie_het_replicates,
                                 required_match_left = 6,
                                 required_match_right = 6,
                                 crispr_del_window = 100,
                                 min_mapq = 0,
                                 max_mismatch_frac = 0.05,
                                 min_aligned_bases = 50,
                                 exclude_multiple_deletions = FALSE,
                                 exclude_nonspanning_reads = TRUE,
                                 allele_profile = FALSE,
                                 del_span_start = -6,
                                 del_span_end = 6,
                                 quiet = F)

#experiment_summary_plot(grep_results, aln_results)
write_results(aln_results, base_path = paste0(dir, "genie_atac_hets/analysis/genie_atac_tagment_comparison.alignment_results"))
saveRDS(aln_results, paste0(dir, "genie_atac_hets/analysis/genie_atac_tagment_comparison.alignment_results.rds"))

pdf("analysis/genie_atac_tagment_comparison.pdf", width=7, height=8)
experiment_summary_plot(grep_results, aln_results)

for (i in 1:length(grep_results)) {
  grep_plot = grep_summary_plot(grep_results[[i]])

  aln_plots = alignment_analysis_plots(aln_results[[i]],
                                       opts = genie_plot_options(),
                                       variance_components_plot = FALSE,
                                       power_plots = FALSE)
  print(grep_plot)
  print(aln_plots)
}

dev.off()



###############################################################################
# genie_atac_hets - PCR input comparison
# Here we compare 3 different amounts of ATAC DNA input for the linear PCR.
setwd(file.path(dir, "genie_atac_hets"))

genie_het_regions = readr::read_tsv("genie_atac_pcr_input.regions.tsv")
genie_het_replicates = readr::read_tsv("genie_atac_pcr_input.replicates.tsv")

grep_results = grep_analysis(genie_het_regions,
                             genie_het_replicates,
                             required_match_left = 10,
                             required_match_right = 10,
                             min_mapq = 0,
                             quiet = F)

#grep_summary_plot(grep_results[[1]])
write_results(grep_results, base_path = paste0(dir, "genie_atac_hets/analysis/genie_atac_pcrinput_comparison.grep_results"))

aln_results = alignment_analysis(genie_het_regions,
                                 genie_het_replicates,
                                 required_match_left = 10,
                                 required_match_right = 10,
                                 crispr_del_window = 100,
                                 min_mapq = 0,
                                 max_mismatch_frac = 0.05,
                                 min_aligned_bases = 50,
                                 exclude_multiple_deletions = FALSE,
                                 exclude_nonspanning_reads = TRUE,
                                 allele_profile = FALSE,
                                 del_span_start = -6,
                                 del_span_end = 6,
                                 quiet = F)

#experiment_summary_plot(grep_results, aln_results)
write_results(aln_results, base_path = paste0(dir, "genie_atac_hets/analysis/genie_atac_pcrinput_comparison.alignment_results"))

pdf("analysis/genie_atac_pcrinput_comparison.pdf", width=7, height=8)
experiment_summary_plot(grep_results, aln_results)

for (i in 1:length(grep_results)) {
  grep_plot = grep_summary_plot(grep_results[[i]])

  aln_plots = alignment_analysis_plots(aln_results[[i]],
                                       opts = genie_plot_options(),
                                       variance_components_plot = FALSE,
                                       power_plots = FALSE)
  print(grep_plot)
  print(aln_plots)
}

dev.off()



###############################################################################
# genie_atac_v2 - 8 SNPs

setwd(file.path(dir, "genie_atac_v2"))
genie_atac_v2_regions = readr::read_tsv("genie_atac_v2.regions.noMIR8070.tsv")
genie_atac_v2_regions = readr::read_tsv("genie_atac_v2.regions.tsv")
genie_atac_v2_replicates = readr::read_tsv("genie_atac_v2.meta.tsv")

grep_results = grep_analysis(genie_atac_v2_regions,
                             genie_atac_v2_replicates,
                             required_match_left = 6,
                             required_match_right = 6,
                             min_mapq = 0,
                             quiet = F)

grep_results = readRDS(paste0(dir, "genie_atac_v2/analysis/genie_atac_v2.grep_results.rds"))
saveRDS(grep_results, paste0(dir, "genie_atac_v2/analysis/genie_atac_v2.grep_results.rds"))
write_results(grep_results, base_path = paste0(dir, "genie_atac_v2/analysis/genie_atac_v2.grep_results"))

aln_results = alignment_analysis(genie_atac_v2_regions,
                                 genie_atac_v2_replicates,
                                 required_match_left = 6,
                                 required_match_right = 6,
                                 crispr_del_window = 100,
                                 min_mapq = 0,
                                 max_mismatch_frac = 0.05,
                                 min_aligned_bases = 50,
                                 exclude_multiple_deletions = FALSE,
                                 exclude_nonspanning_reads = TRUE,
                                 allele_profile = FALSE,
                                 del_span_start = -6,
                                 del_span_end = 6,
                                 quiet = F)

aln_results = readRDS(paste0(dir, "genie_atac_v2/analysis/genie_atac_v2.aln_results.rds"))
saveRDS(aln_results, paste0(dir, "genie_atac_v2/analysis/genie_atac_v2.aln_results.rds"))
write_results(aln_results, base_path = paste0(dir, "genie_atac_v2/analysis/genie_atac_v2.alignment_results"))

experiment_summary_plot(grep_results, aln_results)
#experiment_summary_plot(grep_results, NULL)
#experiment_summary_plot(NULL, aln_results)

# alignment_result = aln_results[[1]]
# alignment_summary_plot(alignment_result)
# aln_plots = alignment_analysis_plots(aln_results[[1]],
#                                      opts = genie_plot_options(),
#                                      variance_components_plot = FALSE,
#                                      power_plots = FALSE)

pdf("analysis/genie_atac_v2.2.pdf", width=7, height=8)
experiment_summary_plot(grep_results, aln_results)

for (i in 1:length(grep_results)) {
  grep_plot = grep_summary_plot(grep_results[[i]])

  aln_plots = alignment_analysis_plots(aln_results[[i]],
                                       opts = genie_plot_options(),
                                       variance_components_plot = FALSE,
                                       power_plots = FALSE)
  print(grep_plot)
  print(aln_plots)
}
dev.off()

allele_effect_plot(aln_results[[i]])

###############################################################################
# ctcf_satmut

setwd(file.path(dir, "genie_atac_ctcf_satmut_mar2021"))
genie_ctcf_regions = readr::read_tsv("genie_atac_ctcf_satmut.regions_indiv.tsv")
genie_ctcf_replicates = readr::read_tsv("genie_atac_ctcf_satmut.replicates.tsv")

grep_results = grep_analysis(genie_ctcf_regions,
                             genie_ctcf_replicates,
                             required_match_left = 6,
                             required_match_right = 6,
                             min_mapq = 0,
                             quiet = F)
write_results(grep_results, base_path = paste0(dir, "genie_atac_ctcf_satmut_mar2021/analysis/genie_atac_ctcf_satmut_mar2021.grep_results"))
saveRDS(grep_results, paste0(dir, "genie_atac_ctcf_satmut_mar2021/analysis/genie_atac_ctcf_satmut_mar2021.grep_results.rds"))

experiment_summary_plot(grep_results, NULL)
grep_summary_plot(grep_results[[1]])

start_time = Sys.time()
aln_results = alignment_analysis(genie_ctcf_regions,
                                 genie_ctcf_replicates,
                                 required_match_left = 6,
                                 required_match_right = 6,
                                 crispr_del_window = 100,
                                 min_mapq = 0,
                                 max_mismatch_frac = 0.05,
                                 min_aligned_bases = 50,
                                 exclude_multiple_deletions = FALSE,
                                 exclude_nonspanning_reads = TRUE,
                                 allele_profile = FALSE,
                                 del_span_start = -10,
                                 del_span_end = 10,
                                 quiet = F)
print(Sys.time() - start_time)
write_results(aln_results, base_path = paste0(dir, "genie_atac_ctcf_satmut_mar2021/analysis/genie_atac_ctcf_satmut_mar2021.alignment_results"))
saveRDS(aln_results, paste0(dir, "genie_atac_ctcf_satmut_mar2021/analysis/genie_atac_ctcf_satmut_mar2021.alignment_results.rds"))

plts = alignment_analysis_plots(aln_results,
                                opts = genie_plot_options(),
                                variance_components_plot = FALSE,
                                power_plots = FALSE)

experiment_summary_plot(grep_results, aln_results)
#experiment_summary_plot(grep_results, NULL)
#experiment_summary_plot(NULL, aln_results)

alignment_result = aln_results[[1]]
alignment_summary_plot(alignment_result)

aln_plots = alignment_analysis_plots(aln_results[[1]],
                                     opts = genie_plot_options(),
                                     variance_components_plot = FALSE,
                                     power_plots = FALSE)

pdf("analysis/genie_atac_ctcf_satmut_mar2021.pdf", width=7, height=8)

experiment_summary_plot(grep_results, aln_results)

for (i in 1:length(grep_results)) {
  grep_plot = grep_summary_plot(grep_results[[i]])
  print(grep_plot)
}
aln_plots = alignment_analysis_plots(aln_results[[i]],
                                     opts = genie_plot_options(),
                                     variance_components_plot = FALSE,
                                     power_plots = FALSE)
print(aln_plots)

dev.off()



###############################################################################
# Random


