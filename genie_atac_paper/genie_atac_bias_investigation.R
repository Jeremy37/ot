suppressMessages(library(tidyverse))
library(rgenie)
library(broom)

dir = "/Users/jeremys/work/opentargets/experiment/atac/"

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

grep_summary_plot(grep_results[[1]])

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
                                 del_span_start = -20,
                                 del_span_end = 20,
                                 quiet = F)

setd9_res = alignment_analysis(genie_atac_v2_regions %>% filter(name == "SETD9"),
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
                               del_span_start = -20,
                               del_span_end = 20,
                               quiet = F)
plts = alignment_analysis_plots(setd9_res[[1]],
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

pdf("analysis/genie_atac_v2.pdf", width=7, height=8)
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


# Investigations with genie_atac_v2 data
aln_dfs = bind_results(aln_results)
allele_effect_df = aln_dfs$allele_effect %>%
  left_join(genie_atac_v2_regions %>% select(name, highlight_site))
allele_effect_df = allele_effect_df %>%
  mutate(del_start_dist = deletion_start - highlight_site,
         del_end_dist = deletion_end - highlight_site,
         min_dist = del_start_dist)
allele_effect_df$min_dist[!is.na(allele_effect_df$del_end_dist) & abs(allele_effect_df$del_end_dist) < abs(allele_effect_df$del_start_dist)] = allele_effect_df$del_end_dist[!is.na(allele_effect_df$del_end_dist) & abs(allele_effect_df$del_end_dist) < abs(allele_effect_df$del_start_dist)]

ggplot(allele_effect_df, aes(x=min_dist, y=log(uns))) +
  geom_point()

# A scatterplot of the deletion effect size log(UNS) vs. the deletion end position,
# stratified by edited region.
pdf(file = "analysis/genie_atac_v2.effect_vs_del_end_dist.scatterplot.pdf", width=9, height=8)
ggplot(allele_effect_df %>% filter(total_count > 200, abs(uns) < 10, abs(del_end_dist) < 30), aes(x=del_end_dist, y=uns, alpha=log10(total_count))) +
  geom_point() +
  facet_wrap(~name, scales = "free_y") +
  geom_smooth(method = "lm", mapping = aes(weight = total_count)) +
  scale_y_log10() +
  ggtitle("Del effect size vs. del-SNP distance")

dev.off()

ggplot(allele_effect_df %>% filter(total_count > 200, abs(uns) < 10, abs(del_end_dist) < 30), aes(x=del_end_dist, y=log(uns), alpha=log10(total_count))) +
  geom_point() +
  facet_wrap(~name, scales = "free_y") +
  geom_smooth(mapping = aes(weight = total_count))


# regions = unique(allele_effect_df$name)
# #regions = "NUCB2"
# for (region in regions) {
#   df = allele_effect_df %>% filter(name == region, total_count > 200, abs(uns) < 10, abs(del_end_dist) < 30) %>%
#     na.omit()
#   print(region)
#   print(summary(lm(log(uns) ~ del_end_dist, data = df, weights = total_count)))
# }

allele_effect_test = allele_effect_df %>%
  filter(total_count > 200, abs(uns) < 10, abs(del_end_dist) < 30) %>%
  na.omit() %>%
  nest(regions = -name) %>%
  mutate(lmres = map(regions, ~ lm(log(uns) ~ del_end_dist, data = .x)),
         tidied = map(lmres, tidy)) %>%
  unnest(tidied)

allele_effect_test %>%
  select(name, term, estimate, p.value) %>%
  arrange(term) %>%
  write_tsv("analysis/genie_atac_v2.effect_vs_del_end_dist.lm_result.tsv")


# Look at the fraction of reads with deletions which start/end at a given distance from the cut site
del_effect_df = allele_effect_df %>%
  filter(has_crispr_deletion) %>%
  left_join(genie_atac_v2_regions %>% select(name, cut_site)) %>%
  mutate(del_cut_site_dist = if_else(deletion_end <= cut_site, deletion_end - cut_site,
                                     if_else(deletion_start >= cut_site, deletion_start - cut_site, 0)))

get_num_deleted_bases = function(str) {
  chars = strsplit(str, "")[[1]]
  sum(chars == "*")
}

del_pos_df = del_effect_df %>%
  group_by(name, del_cut_site_dist) %>%
  arrange(desc(gDNA_total_count)) %>%
  summarise(udp = first(udp),
            gDNA_read_count = sum(gDNA_total_count)) %>%
  select(name, udp, del_cut_site_dist, gDNA_read_count) %>%
  group_by(name, udp) %>%
  mutate(num_deleted_bases = get_num_deleted_bases(udp),
         single_base_del = num_deleted_bases == 1)

del_cdf_df = del_pos_df %>%
  group_by(name) %>%
  arrange(del_cut_site_dist) %>%
  mutate(cum_pct = cumsum(gDNA_read_count / sum(gDNA_read_count)))
ecdf(del_pos_df$del_cut_site_dist)

# A plot showing the number of gDNA reads with a deletion that ends at a given
# distance from the cut site, stratified by edited region.
pdf(file = "analysis/genie_atac_v2.gDNA_reads_vs_del_end_dist.scatterplot.pdf", width=9, height=8)

ggplot(del_pos_df, aes(x=del_cut_site_dist, y=gDNA_read_count, col=single_base_del)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ name, ncol = 3) +
  scale_y_log10() +
  ggtitle("Del - cut site distance vs. gDNA read count")

ggplot(del_pos_df, aes(x=del_cut_site_dist, y=gDNA_read_count, col=single_base_del)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ name, ncol = 3) +
  scale_y_log10() +
  coord_cartesian(xlim = c(-35, 35)) +
  ggtitle("(ZOOM) Del - cut site distance vs. gDNA read count")

ggplot(del_cdf_df, aes(x=del_cut_site_dist, y=cum_pct)) +
  geom_line() +
  facet_wrap(~ name, ncol = 3) +
  ggtitle("Cumulative read count by del-cut site distance")

dev.off()
