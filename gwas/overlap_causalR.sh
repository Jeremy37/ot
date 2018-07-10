#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
AD=$JS/gwas/AD

cd $AD/genes

# Get genes within a window of AD lead SNPs
echo -e "leadSNP\tlocus\tgeneID\tsymbol\tchr\tgene_start\tgene_end" > AD.leadSNPs.1Mb_window.gene_overlaps.txt
bedtools intersect -a AD.leadSNPs.1Mb_window.chr.bed -b $JS/reference/GRCh37/gencode.v28lift37.genes.pc.bed -wa -wb \
    | awk 'BEGIN{OFS="\t"}{print $4,$5,$10,$11,$6,$7,$8}' | less >> AD.leadSNPs.1Mb_window.gene_overlaps.txt


cd $AD/causalR

Rscript $JS/src/gwas/overlap_causalR.R ../genes/AD.leadSNPs.1Mb_window.gene_overlaps.txt MAYO_ADvCntl_DiffExp_ALL.txt \
   FullCausalRresults_MAYO_ADvCntl.pathlength_1.txt \
   FullCausalRresults_MAYO_ADvCntl.pathlength_2.txt \
   FullCausalRresults_MAYO_ADvCntl.pathlength_3.txt \
   FullCausalRresults_MAYO_ADvCntl.pathlength_4.txt \
   > AD.GWAS_locus_genes.CausalR_annotated.txt

Rscript $JS/src/gwas/overlap_causalR.R ../genes/AD.leadSNPs.1Mb_window.gene_overlaps.txt MAYO_ADvCntl_DiffExp_ALL.txt \
   FullCausalRresults_MAYO_ADvCntl.CER.pathlength_1.txt \
   FullCausalRresults_MAYO_ADvCntl.CER.pathlength_2.txt \
   FullCausalRresults_MAYO_ADvCntl.CER.pathlength_3.txt \
   FullCausalRresults_MAYO_ADvCntl.CER.pathlength_4.txt \
   > AD.GWAS_locus_genes.CausalR_annotated.CER.txt

Rscript $JS/src/gwas/overlap_causalR.R ../genes/AD.leadSNPs.1Mb_window.gene_overlaps.txt MAYO_ADvCntl_DiffExp_ALL.txt \
   FullCausalRresults_MAYO_ADvCntl.TCX.pathlength_1.txt \
   FullCausalRresults_MAYO_ADvCntl.TCX.pathlength_2.txt \
   FullCausalRresults_MAYO_ADvCntl.TCX.pathlength_3.txt \
   FullCausalRresults_MAYO_ADvCntl.TCX.pathlength_4.txt \
   > AD.GWAS_locus_genes.CausalR_annotated.TCX.txt
   
