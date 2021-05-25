library(tidyverse)
library(gridExtra)
library(ggseqlogo)

dir = "/Users/jeremys/work/opentargets/experiment/transcribed/CCDC6_satmut_genie/analysis/"
setwd(dir)

file_path = file.path(dir, "CCDC6_satmut_genie.allHDR.lmer.region_stats.tsv")
file_path = file.path(dir, "CCDC6_satmut_genie.HDR.by_batch.region_stats.tsv")

spliceai.df = readr::read_tsv(file.path(dir, "spliceai.ccdc6_region.tsv")) %>%
  rename(pos_hg37 = pos)
genie.df = readr::read_tsv(file.path(dir, "CCDC6_satmut_genie.allHDR.lmer.region_stats.tsv")) %>%
  separate(name, sep="_", into=c("pos", "ref", "alt")) %>%
  select(-starts_with("del"))
colnames(genie.df) = gsub(".error_prop", "", colnames(genie.df))
genie.df$pos[genie.df$pos == "rs1171830"] = "59906128"
genie.df = genie.df %>% mutate(pos = as.integer(pos)) %>% filter(!is.na(pos))

# Update GenIE values for rs1171830 to reflect what you would get if Kolf2 were
# hom ref (C) and were mutated in the opposite direction
rs1171830.df = genie.df[genie.df$pos == 59906128,]
rs1171830.df$alt[rs1171830.df$alt == "C"] = "A"
rs1171830.df$ref = "C"
rs1171830.df$hdr.effect = 1 / rs1171830.df$hdr.effect
genie.df = bind_rows(genie.df, rs1171830.df)

genie.df$pos_hg37 = genie.df$pos + (61665886 - 59906128)
merged.df = genie.df %>%
  left_join(spliceai.df, by=c("pos_hg37", "ref", "alt")) %>%
  rename(pos_hg38 = pos) %>%
  na.omit()

p1 = ggplot(merged.df, aes(x=DS_DL, y=hdr.effect)) +
  geom_point(alpha = 0.7) +
  geom_errorbar(aes(ymin=hdr.effect_confint_lo, ymax=hdr.effect_confint_hi), width=0.02, alpha=0.5) +
  xlab("SpliceAI donor loss score") + ylab("GenIE HDR effect size") +
  geom_hline(yintercept = 1.0, col = "green") +
  coord_cartesian(ylim=c(0, 8)) +
  ggtitle("GenIE effect size vs. SpliceAI score")

# Or with a log scale
p2 = ggplot(merged.df, aes(x=DS_DL, y=hdr.effect)) +
  geom_point(alpha = 0.7) +
  geom_errorbar(aes(ymin=hdr.effect_confint_lo, ymax=hdr.effect_confint_hi), width=0.02, alpha=0.5) +
  xlab("SpliceAI donor loss score") + ylab("GenIE HDR effect size") +
  geom_hline(yintercept = 1.0, col = "green") +
  ggtitle("GenIE effect size vs. SpliceAI score (log scale)") +
  scale_y_log10()

merged.df = merged.df %>%
  select(chr, pos_hg38, pos_hg37, everything())
write.table(merged.df, file="ccdc6.spliceai_vs_genie.tsv", quote=F, col.names=T, row.names=F, sep="\t")

# Make a barplot of position vs. both Genie effect and SpliceAI score
# to do this nicely we need to transform the values to be on the same
# scale.
max_genie_effect = max(merged.df$hdr.effect) - 1
merged.df$DS_DL_norm = merged.df$DS_DL * max_genie_effect
merged.bypos.df = merged.df %>%
  select(pos_hg38, ref, alt, GenIE = hdr.effect, SpliceAI = DS_DL_norm, hdr.confint_lo = hdr.effect_confint_lo, hdr.confint_hi = hdr.effect_confint_hi) %>%
  mutate(variant = paste0(pos_hg38, "_", ref, ">", alt),
         GenIE = GenIE - 1, hdr.confint_lo = hdr.confint_lo - 1, hdr.confint_hi = hdr.confint_hi - 1) %>%
gather(key=Type, value=value, -variant, -pos_hg38, -ref, -alt, -hdr.confint_lo, -hdr.confint_hi)

merged.bypos.df$hdr.confint_lo[merged.bypos.df$Type == "SpliceAI"] = NA
merged.bypos.df$hdr.confint_hi[merged.bypos.df$Type == "SpliceAI"] = NA

p3 = ggplot(merged.bypos.df, aes(x=variant, y=value, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=hdr.confint_lo, ymax=hdr.confint_hi, group=Type), position=position_dodge(.9), alpha=0.3) +
  theme_bw(11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8),
        legend.justification = c(1, 1), legend.position = c(0.9, 0.95)) +
  ylab("Genie Effect - 1") +
  scale_y_continuous(sec.axis = sec_axis(~ . / max_genie_effect, name="SpliceAI donor loss score")) +
  coord_cartesian(ylim=c(-0.5, 6.5)) +
  ggtitle("GenIE vs. SpliceAI by position")

# Plot a sequence logo representing the effect sizes of the various mutations.
# The idea is to show that our results recapitulate known splice site motifs.
# Since we're measuring reads that aren't spliced, we basically want the inverse
# of the number of reads as a measure of the level of splicing.
ccdc6.seqlogo.df = readr::read_tsv("ccdc6_region.seqlogo_summary.tsv") %>%
  mutate(effect_inv = 1 / hdr_effect)
ccdc6.seqlogo.summary = ccdc6.seqlogo.df %>% group_by(pos) %>% summarise(effect_inv_sum = sum(effect_inv))
ccdc6.seqlogo.df = ccdc6.seqlogo.df %>%
  left_join(ccdc6.seqlogo.summary, by="pos") %>%
  mutate(effect_norm = effect_inv / effect_inv_sum)
plot.df = ccdc6.seqlogo.df


getSeqLogoData = function(merged.df, col) {
  minpos = min(merged.df$pos_hg38)
  maxpos = max(merged.df$pos_hg38)
  genie.logo.data = merged.df %>% select(pos_hg38, alt, one_of(col)) %>%
    spread(key = alt, value = col)
  genie.logo.data = rbind(genie.logo.data,
                          data.frame(pos_hg38 = setdiff(minpos:maxpos, genie.logo.data$pos_hg38),
                                     A = NA, C = NA, G = NA, T = NA)) %>%
    arrange(pos_hg38) %>% as.data.frame()
  rownames(genie.logo.data) = genie.logo.data$pos_hg38
  genie.logo.data %>% select(-pos_hg38) %>% t() %>% as.matrix() %>% round(digits=3)
}
p.genie.seqlogo = ggseqlogo(getSeqLogoData(merged.df, "hdr.effect") - 1, method='custom', seq_type='dna') +
  theme(text = element_text(size=10)) +
  ggtitle("Sequence logo of GenIE effect size")

p.spliceai.seqlogo = ggseqlogo(getSeqLogoData(merged.df, "DS_DL"), method='custom', seq_type='dna') +
  theme(text = element_text(size=10)) +
  ggtitle("SpliceAI DS_DL score")


genie.logo.data = merged.df %>% select(pos_hg38, alt, hdr.effect) %>%
  spread(key = alt, value = hdr.effect) %>%
  filter(59906118 <= pos_hg38, pos_hg38 <= 59906122)
genie.logo.data[is.na(genie.logo.data)] = 1
genie.logo.data = 1 / (genie.logo.data %>% select(-pos_hg38) %>% t() %>% as.matrix())
genie.logo.data = genie.logo.data[, ncol(genie.logo.data):1]
rownames(genie.logo.data) = c("T", "G", "C", "A")

p.genie.importance1 = ggseqlogo(genie.logo.data, method='bits', seq_type='dna') +
  theme(text = element_text(size=10)) +
  ggtitle("GenIE nucleotide importance")
p.genie.importance2 = ggseqlogo(genie.logo.data, method='prob', seq_type='dna') +
  theme(text = element_text(size=10)) +
  ggtitle("GenIE nucleotide importance")


pdf("ccdc6.genie_splicing.pdf", width=7, height=6)
p1
p2
p3
grid.arrange(p.genie.seqlogo,
             p.spliceai.seqlogo, 
             p.genie.importance1, 
             p.genie.importance2, ncol=2)
dev.off()

