#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
AD=$JS/gwas/AD

cd $AD/genes

# Get genes within a window of AD lead SNPs
echo -e "leadSNP\tlocus\tgeneID\tsymbol\tchr\tgene_start\tgene_end" > AD.leadSNPs.200kb_window.gene_overlaps.txt
bedtools intersect -a AD.leadSNPs.200kb_window.chr.bed -b $JS/reference/GRCh37/gencode.v28lift37.genes.pc.bed -wa -wb \
    | awk 'BEGIN{OFS="\t"}{print $4,$5,$10,$11,$6,$7,$8}' | less >> AD.leadSNPs.200kb_window.gene_overlaps.txt

echo -e "leadSNP\tlocus\tgeneID\tsymbol\tchr\tgene_start\tgene_end" > AD.leadSNPs.1Mb_window.gene_overlaps.txt
bedtools intersect -a AD.leadSNPs.1Mb_window.chr.bed -b $JS/reference/GRCh37/gencode.v28lift37.genes.pc.bed -wa -wb \
    | awk 'BEGIN{OFS="\t"}{print $4,$5,$10,$11,$6,$7,$8}' | less >> AD.leadSNPs.1Mb_window.gene_overlaps.txt

echo -e "leadSNP\tlocus\tgeneID\tsymbol\tchr\tgene_start\tgene_end" > AD.leadSNPs.2Mb_window.gene_overlaps.txt
bedtools intersect -a AD.leadSNPs.2Mb_window.chr.bed -b $JS/reference/GRCh37/gencode.v28lift37.genes.pc.bed -wa -wb \
    | awk 'BEGIN{OFS="\t"}{print $4,$5,$10,$11,$6,$7,$8}' | less >> AD.leadSNPs.2Mb_window.gene_overlaps.txt

Rscript ~/src/R/aggregateColumn.R --file AD.leadSNPs.200kb_window.gene_overlaps.txt --groupby leadSNP --aggregate symbol > AD.leadSNPs.200kb_window.gene_overlaps.summary.txt
Rscript ~/src/R/aggregateColumn.R --file AD.leadSNPs.1Mb_window.gene_overlaps.txt --groupby leadSNP --aggregate symbol > AD.leadSNPs.1Mb_window.gene_overlaps.summary.txt

Rscript $JS/src/gwas/annotateCandidateGenes.R AD.leadSNPs.1Mb_window.gene_overlaps.txt $JS/reference/tissueRPKM/tissues.selected.rpkm_average.txt > AD.leadSNPs.1Mb_window.gene_overlaps.rpkms.txt


# Get genes within a window of PD lead SNPs
cd $JS/gwas/PD/genes
echo -e "leadSNP\tgeneID\tsymbol\tchr\tgene_start\tgene_end" > PD.leadSNPs.200kb_window.gene_overlaps.txt
bedtools intersect -a PD.leadSNPs.200kb_window.chr.bed -b $JS/reference/GRCh37/gencode.v28lift37.genes.pc.bed -wa -wb \
    | awk 'BEGIN{OFS="\t"}{print $4,$9,$10,$5,$6,$7}' | less >> PD.leadSNPs.200kb_window.gene_overlaps.txt

echo -e "leadSNP\tgeneID\tsymbol\tchr\tgene_start\tgene_end" > PD.leadSNPs.1Mb_window.gene_overlaps.txt
bedtools intersect -a PD.leadSNPs.1Mb_window.chr.bed -b $JS/reference/GRCh37/gencode.v28lift37.genes.pc.bed -wa -wb \
    | awk 'BEGIN{OFS="\t"}{print $4,$9,$10,$5,$6,$7}' | less >> PD.leadSNPs.1Mb_window.gene_overlaps.txt

echo -e "leadSNP\tgeneID\tsymbol\tchr\tgene_start\tgene_end" > PD.leadSNPs.2Mb_window.gene_overlaps.txt
bedtools intersect -a PD.leadSNPs.2Mb_window.chr.bed -b $JS/reference/GRCh37/gencode.v28lift37.genes.pc.bed -wa -wb \
    | awk 'BEGIN{OFS="\t"}{print $4,$9,$10,$5,$6,$7}' | less >> PD.leadSNPs.2Mb_window.gene_overlaps.txt

Rscript ~/src/R/aggregateColumn.R --file PD.leadSNPs.200kb_window.gene_overlaps.txt --groupby leadSNP --aggregate symbol > PD.leadSNPs.200kb_window.gene_overlaps.summary.txt
Rscript ~/src/R/aggregateColumn.R --file PD.leadSNPs.1Mb_window.gene_overlaps.txt --groupby leadSNP --aggregate symbol > PD.leadSNPs.1Mb_window.gene_overlaps.summary.txt

Rscript $JS/src/gwas/annotateCandidateGenes.R PD.leadSNPs.1Mb_window.gene_overlaps.txt $JS/reference/tissueRPKM/tissues.selected.rpkm_average.txt > PD.leadSNPs.1Mb_window.gene_overlaps.rpkms.txt
