#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

# Get variants near AGAP1 in VCF format for running DeepSEA
tabix $JS/datasets/GWAS/Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.sorted.txt.gz 2:236720731-237120731 > $JS/gwas/PD/data/Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.sorted.near_AGAP1.txt

hashJoin.pl --hashFile $JS/gwas/PD/data/Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.sorted.near_AGAP1.txt \
  --scanFile /warehouse/compgen_wh05/js29/sensoryneurons/reference/1000genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz \
  --colHashFile 1 --colScanFile 3 \
  > $JS/gwas/PD/data/AGAP1_variants.txt
#  --scanFile $JS/reference/GRCh37/dbsnp/1000genomes.20130502.sites.cut.gz \

