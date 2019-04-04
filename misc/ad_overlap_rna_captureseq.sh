#!/bin/bash
OT=/lustre/scratch115/realdata/mdt3/projects/otcoregen
JS=$OT/jeremys
AD=$JS/AD_finemap

cd $AD

# First we need to liftover our SNP positions to GRCh38, since the Mattick
# RNA-CaptureSeq annotations are GRCh38.
submitJobs.py --MEM 1000 -j CrossMap.microglia.GRCh38_to_GRCh37 -q yesterday -n 2 \
    -c "CrossMap.py vcf $JS/software/CrossMap/GRCh38_to_GRCh37.chain.gz microglia.GRCh38.vcf.gz $JS/reference/GRCh37/GRCh37.75.dna.toplevel.fa.gz microglia.GRCh37.vcf"

sed '1d' $AD/annotated/AD.credible_sets.annotated.colocs.prob_gt_1e-4.reordered.txt \
  | awk 'BEGIN{FS="\t";OFS="\t"} {print $2,$3,$3,$4}' > AD.credible_sets.prob_gt_1e-4.snps.GRCh37.bed

CrossMap.py bed $JS/software/CrossMap/GRCh37_to_GRCh38.chain.gz \
    AD.credible_sets.prob_gt_1e-4.snps.GRCh37.bed \
    AD.credible_sets.prob_gt_1e-4.snps.GRCh38.bed


# In principle, not all SNPs may have crosmapped properly, in which case the order
# of SNPs in AD.credible_sets.prob_gt_1e-4.snps.GRCh38.txt might not match up and
# everything would be screwed. But I've checked, and all match currently. Need to
# be careful if running this again. (This would probably better be done as a left
# join.)
#(echo -e "pos_b38\trsID"; cut -f 3-4 AD.credible_sets.prob_gt_1e-4.snps.GRCh38.bed) > AD.credible_sets.prob_gt_1e-4.snps.GRCh38.txt
(echo -e "pos_b38\trsID"; cut -f 3 AD.credible_sets.prob_gt_1e-4.snps.GRCh38.bed) > AD.credible_sets.prob_gt_1e-4.snps.GRCh38.txt

paste <(cut -f 1-3 $AD/annotated/AD.credible_sets.annotated.colocs.prob_gt_1e-4.reordered.txt) \
      AD.credible_sets.prob_gt_1e-4.snps.GRCh38.txt \
      <(cut -f 1-3 --complement $AD/annotated/AD.credible_sets.annotated.colocs.prob_gt_1e-4.reordered.txt) \
      > $AD/annotated/AD.credible_sets.annotated.colocs.prob_gt_1e-4.reordered.grch38.tsv


sed '1d' $AD/annotated/AD.credible_sets.annotated.colocs.prob_gt_1e-4.reordered.grch38.tsv \
  | awk 'BEGIN{FS="\t";OFS="\t"}{print "chr"$2,$4,$4,$5,$12,$1,$13,$14}' \
  | sort -k1,1 -k2,2n > $AD/AD.credset_snps.prob_gt_1e-4.grch38.bed

echo -e "locus\tchr\tpos\trsid\tfinemap_prob\tconsequence\tsymbol\tdistance_to_rna_captureseq" > AD.credset_snps.prob_gt_1e-4.overlap.rna_captureseq.tsv
bedtools closest -D a -t first \
  -a $AD/AD.credset_snps.prob_gt_1e-4.grch38.bed \
  -b $JS/datasets/RNA_CaptureSeq/mattick_neuropsych/GSE118158_filtered_hybrid_transcriptome.bed \
  | awk 'BEGIN{FS="\t";OFS="\t"}{print $6,$1,$2,$4,$5,$7,$8,$21}' >> AD.credset_snps.prob_gt_1e-4.overlap.rna_captureseq.tsv

cat AD.credset_snps.prob_gt_1e-4.overlap.rna_captureseq.tsv | awk 'BEGIN{FS="\t";OFS="\t"}{if ($5 != "" && $5 >= 0.05) print}' > AD.credset_snps.prob_gt_1e-4.overlap.rna_captureseq.p_gt_0.05.tsv

