library(tidyverse)

setwd("/Users/jeremys/work/opentargets/genie_atac/plots")
theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))

# barplot_theme = theme_bw(10) +
#   theme(axis.line = element_line(colour = "black", size=0.25),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(size = 10))

atac_dir = "/Users/jeremys/work/opentargets/experiment/atac"


###############################################################################
# Het experiment

genie_hets = read_tsv(file.path(atac_dir, "genie_atac_hets/analysis", "genie_atac_hets.alignment_results.replicate_stats.saved.tsv"))

genie_hets = genie_hets %>%
  mutate(hdr_frac = num_hdr_reads / (num_wt_reads + num_hdr_reads),
         ratio = num_hdr_reads / num_wt_reads,
         group = paste(rsid, type, direction)) %>%
  rename(Type = type)

allelic_ratios = read_tsv(file.path(atac_dir, "genie_atac_hets/analysis", "het_site_allelic_imbalance.txt")) %>%
  mutate(ratio = alt_reads / ref_reads,
         Type = "AS ratio",
         direction = "",
         group = paste(rsid, Type, direction))

plot_df = bind_rows(genie_hets %>% select(rsid, Type, ratio, group, direction),
                    allelic_ratios %>% select(rsid, Type, ratio, group))

positive_snps = c("rs12212281", "rs10826861", "rs12046948")
negative_snps = c("rs11233493", "rs6450345", "rs113541101")
plot_df$control = "negative controls"
plot_df$control[plot_df$rsid %in% positive_snps] = "positive controls"
plot_df$control = factor(plot_df$control, levels = c("positive controls", "negative controls"))
plot_df = plot_df %>%
  mutate(rsid = factor(rsid, levels = c(positive_snps, negative_snps)))

# genie_hets_pos = plot_df %>%
#   filter(rsid %in% positive_snps)
# genie_hets_neg = plot_df %>%
#   filter(!rsid %in% positive_snps)

p = ggplot(plot_df %>% filter(!(Type == "gDNA" & direction == "rev")),
           aes(x=rsid, y=ratio, fill=Type, group=group)) +
  geom_boxplot() +
  geom_jitter(alpha=0.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  facet_wrap(~control, ncol = 1, scales = "free_x") +
  ylab("Alt / Ref ratio") +
  xlab("Targeted SNP")

pdf("genie_hets.v1.pdf", width = 7.5, height = 6)
print(p)
dev.off()

plot_df$Type[plot_df$Type == "gDNA"] = "PCR"
plot_df$Type[plot_df$Type == "AS ratio"] = "Allele-specific\nATAC reads"
plot_df$Type = factor(as.character(plot_df$Type), levels = c("PCR", "ATAC", "Allele-specific\nATAC reads"))

p = ggplot(plot_df %>% filter(!(Type == "PCR" & direction == "rev")),
           aes(x=Type, y=ratio, fill=Type, group=group)) +
  geom_boxplot() +
  geom_jitter(alpha=0.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  facet_wrap(~rsid, ncol = 3) +
  ylab("Alt / Ref ratio") +
  xlab("Targeted SNP")

pdf("genie_hets.v2.pdf", width = 7.5, height = 6)
print(p)
dev.off()


p = ggplot(plot_df %>% filter(is.na(direction) | direction != "rev"),
           aes(x=Type, y=ratio, fill=Type, group=group)) +
  geom_boxplot(aes(color=Type)) +
  geom_jitter(alpha=0.75, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  facet_wrap(~rsid, ncol = 3) +
  ylab("Alt / Ref ratio") +
  xlab("Targeted SNP")

pdf("genie_hets.fwd.pdf", width = 7.5, height = 6)
print(p)
dev.off()

p = ggplot(plot_df %>% filter(is.na(direction) | direction != "rev"),
           aes(x=Type, y=ratio, fill=Type, group=group, color=Type)) +
  geom_jitter(size=3, alpha=0.75, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  facet_wrap(~rsid, ncol = 3) +
  ylab("Alt / Ref ratio") +
  xlab("Targeted SNP")

pdf("genie_hets.fwd.v2.pdf", width = 7.5, height = 6)
print(p)
dev.off()


p = ggplot(plot_df %>% filter(is.na(direction) | direction == "rev", control == "positive controls"),
           aes(x=Type, y=ratio, fill=Type, group=group)) +
  geom_boxplot(aes(color=Type)) +
  geom_jitter(alpha=0.75, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  facet_wrap(~rsid, ncol = 3) +
  ylab("Alt / Ref ratio") +
  xlab("Targeted SNP")

pdf("genie_hets.rev.pdf", width = 7.5, height = 3.5)
print(p)
dev.off()


###############################################################################
# Het experiment - considering gDNA only - tagmented fwd, reverse, or amplicon

genie_hets = read_tsv(file.path(atac_dir, "genie_atac_hets/analysis", "genie_atac_tagment_comparison.alignment_results.replicate_stats.tsv"))
#genie_hets = read_tsv(file.path(atac_dir, "genie_atac_hets/analysis", "genie_atac_tagment_comparison.grep_results.replicate_stats.tsv"))

genie_hets = genie_hets %>%
  mutate(hdr_frac = num_hdr_reads / (num_wt_reads + num_hdr_reads),
         ratio = num_hdr_reads / num_wt_reads,
         group = paste(name, type)) %>%
  rename(Type = type)

p = ggplot(genie_hets,
           aes(x=name, y=ratio, fill=Type, group=group)) +
  geom_boxplot() +
  geom_jitter(alpha=0.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  ylab("Alt / Ref ratio") +
  xlab("Targeted SNP")
p
#facet_wrap(~rsid, ncol = 3) +

genie_hets$AmpType = "Amplicon"
genie_hets$AmpType[grepl("tag_fwd", genie_hets$name)] = "Tagmented Fwd"
genie_hets$AmpType[grepl("tag_rev", genie_hets$name)] = "Tagmented Rev"

getRsid = function(x) {
  str_split(x, "_")[[1]][2]
}
genie_hets$SNP = sapply(genie_hets$name, getRsid)

p = ggplot(genie_hets %>% filter(Type == "gDNA"), aes(x=AmpType, y=ratio, fill="fill", group=AmpType)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(shape=SNP, group=AmpType, color=SNP), alpha=0.8, position = position_jitter(width = 0.2), size=3) +
  scale_fill_manual(values = c(fill = "grey80"), guide=F) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  ylab("Alt / Ref ratio") +
  xlab("gDNA amplification type")
p

pdf("genie_hets.tagment_comparison.pdf", width = 7.5, height = 6)
print(p)
dev.off()



###############################################################################
# Testing different amounts of DNA input for PCR

genie_rep = read_tsv(file.path(atac_dir, "genie_atac_hets/analysis", "genie_atac_pcrinput_comparison.alignment_results.replicate_stats.tsv"))
#genie_rep = read_tsv(file.path(atac_dir, "genie_atac_hets/analysis", "genie_atac_pcrinput_comparison.grep_results.replicate_stats.tsv"))

genie_rep = genie_rep %>%
  filter(type == "gDNA") %>%
  mutate(hdr_frac = num_hdr_reads / (num_wt_reads + num_hdr_reads),
         ratio = num_hdr_reads / num_wt_reads,
         dna_amount = gsub("ul", " ul", gsub("2_rs12269414_", "", name)),
         group = paste(name, type)) %>%
  rename(Type = type)
genie_rep$dna_amount = factor(genie_rep$dna_amount, levels = c("15 ul", "10 ul", "5 ul"))

p = ggplot(genie_rep,
           aes(x=dna_amount, y=ratio, fill=Type, group=group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("grey80"), guide = F) +
  geom_jitter(size = 3, alpha = 0.75, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  ylab("Alt / Ref ratio") +
  xlab("Amount of input DNA")
p

pdf("genie_hets.pcr_input_comparison.pdf", width = 7.5, height = 6)
print(p)
dev.off()




###############################################################################
# Editing SNPs (v2)



