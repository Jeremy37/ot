#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

# I got gene positions from UCSC, and then manually added +/- 500 kb to these regions
cd $JS/gwas/AD/data
tabix $JS/datasets/GWAS/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.txt.gz 21:26752861-28043138 > Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.near_APP.txt
tabix $JS/datasets/GWAS/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.txt.gz 1:226558273-227583804 > Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.near_PSEN2.txt
tabix $JS/datasets/GWAS/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.txt.gz 14:73103143-74190399 > Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.near_PSEN1.txt
tabix $JS/datasets/GWAS/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.txt.gz 6:40626244-41630924 > Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.near_TREM2.txt
tabix $JS/datasets/GWAS/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.txt.gz 19:35895303-36899211 > Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.near_TYROBP.txt

F=$JS/datasets/GWAS/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta_v2.sorted.txt.gz
(gzhead 1 $F; tabix $F 6:2576998-3615421) > Alzheimers_disease_Liu_2017_UKBB_GWAX_meta_v2.sorted.near_RIPK1.txt


# Run the following in R
pdf("AD_GWAS_near_APP.pdf", width=8, height=6)
df = readr::read_tsv("/Users/jeremys/work/opentargets/gwas/AD/data/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.near_APP.txt")
ggplot(df, aes(x=BP, y=-log10(META_P))) + geom_point(alpha=0.7) + theme_bw() + ggtitle("AD GWAS near APP")
dev.off()

pdf("AD_GWAS_near_PSEN1.pdf", width=8, height=6)
df = readr::read_tsv("/Users/jeremys/work/opentargets/gwas/AD/data/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.near_PSEN1.txt")
ggplot(df, aes(x=BP, y=-log10(META_P))) + geom_point(alpha=0.7) + theme_bw() + ggtitle("AD GWAS near PSEN1")
dev.off()

pdf("AD_GWAS_near_PSEN2.pdf", width=8, height=6)
df = readr::read_tsv("/Users/jeremys/work/opentargets/gwas/AD/data/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.near_PSEN2.txt")
ggplot(df, aes(x=BP, y=-log10(META_P))) + geom_point(alpha=0.7) + theme_bw() + ggtitle("AD GWAS near PSEN2")
dev.off()

pdf("AD_GWAS_near_TREM2.pdf", width=8, height=6)
df = readr::read_tsv("/Users/jeremys/work/opentargets/gwas/AD/data/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.near_TREM2.txt")
ggplot(df, aes(x=BP, y=-log10(META_P))) + geom_point(alpha=0.7) + theme_bw() + ggtitle("AD GWAS near TREM2")
dev.off()

pdf("AD_GWAS_near_TYROBP.pdf", width=8, height=6)
df = readr::read_tsv("/Users/jeremys/work/opentargets/gwas/AD/data/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.near_TYROBP.txt")
ggplot(df, aes(x=BP, y=-log10(META_P))) + geom_point(alpha=0.7) + theme_bw() + ggtitle("AD GWAS near TYROBP")
dev.off()


pdf("AD_GWAS_near_RIPK1.pdf", width=8, height=6)
df = readr::read_tsv("/Users/jeremys/work/opentargets/gwas/AD/data/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta_v2.sorted.near_RIPK1.txt")
ggplot(df, aes(x=BP, y=-log10(META_P))) + geom_point(alpha=0.7) + theme_bw() + ggtitle("AD GWAS near RIPK1")
dev.off()


# Get variants near AGAP1 in VCF format for running DeepSEA
cd $JS/gwas/PD/data
tabix $JS/datasets/GWAS/Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.sorted.txt.gz 2:236720731-237120731 > $JS/gwas/PD/data/Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.sorted.near_AGAP1.txt

hashJoin.pl --hashFile $JS/gwas/PD/data/Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.sorted.near_AGAP1.txt \
  --scanFile /warehouse/compgen_wh05/js29/sensoryneurons/reference/1000genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz \
  --colHashFile 1 --colScanFile 3 \
  > $JS/gwas/PD/data/AGAP1_variants.txt
#  --scanFile $JS/reference/GRCh37/dbsnp/1000genomes.20130502.sites.cut.gz \


gzhead 1 $JS/datasets/GWAS/PD.proxy.bgen.stats.cut.gz > $JS/gwas/PD/data/PD.proxy.bgen.stats.cut.near_SIPA1L2.txt
tabix $JS/datasets/GWAS/PD.proxy.bgen.stats.cut.gz 1:232333712-232851243 >> $JS/gwas/PD/data/PD.proxy.bgen.stats.cut.near_SIPA1L2.txt


################################################################################
# New variants to check
F=$JS/datasets/GWAS/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta_v2.sorted.txt.gz
(gzhead 1 $F; tabix $F 5:74151084-75151084) > Alzheimers_disease_Liu_2017_UKBB_GWAX_meta_v2.sorted.near_HMGCR.txt
(gzhead 1 $F; zcat $F | grep -P "rs3846662\t") > Alzheimers_disease_Liu_2017_UKBB_GWAX_meta_v2.rs3846662.txt

setwd("/Users/jeremys/work/opentargets/gwas/AD/data")
pdf("AD_GWAS_near_HMGCR.pdf", width=8, height=6)
df = readr::read_tsv("/Users/jeremys/work/opentargets/gwas/AD/data/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta_v2.sorted.near_HMGCR.txt")
ggplot(df, aes(x=BP, y=-log10(META_P))) + geom_point(alpha=0.7) + theme_bw() + ggtitle("AD GWAS near HMGCR")
dev.off()

