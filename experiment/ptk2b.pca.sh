#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

#mkdir $JS/macrophage/genotype
cd $JS/macrophage/genotype


# Get genotypes for putative PTK2B causal SNP rs28834970
(zcat $JS/macrophage/genotypes/imputed.86_samples.sorted.filtered.named.vcf.gz | head -n 200 | grep "^#CHROM";
 zcat $JS/macrophage/genotypes/imputed.86_samples.sorted.filtered.named.vcf.gz \
   | grep -P "rs28834970\t") > imputed.86_samples.rs28834970.vcf

