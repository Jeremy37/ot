#!/bin/bash
# This script is to get candidate causal ATAC QTL SNPs in macrophages which
# are also in H3K27ac peaks in Blueprint macrphages or monocytes, for Erica
# to use in her H3K27ac ChIP-seq CRISPR assays.
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

cd $JS/macrophage/GRCh37/ATAC

MAC_ATAC_QTL=/lustre/scratch117/cellgen/team170/ka8/projects/macrophage-gxe-study/results/ATAC/rasqual/output_fixed_gt/naive_100kb/naive_50kb.eigenMT.txt
wc -l $MAC_ATAC_QTL
# Filter to ATAC peaks with p less than a threshold
cat $MAC_ATAC_QTL | awk '$4 < 1e-3' > macrophage.naive_50kb.eigenMT.p_lt_1e-3.txt

# Get the set of K27ac QTLs from Blueprint. Convert to BED so that we can overlap
# with the macrophage QTLs
zcat $JS/datasets/blueprint/mono_K27AC.gene_minp.txt.gz \
  | perl -ane '@peak=split(/:/,$F[0]); print(join("\t", @peak, @F)."\n");' \
  | gzip > mono_K27AC.gene_minp.bed.gz

# Join the macrophage QTL data with dbSNP in GRCh37 coords (Uses too much memory)
#hashJoin.pl --hashFile $JS/reference/GRCh37/dbsnp/1000genomes.20130502.sites.cut.gz --scanFile macrophage.naive_50kb.eigenMT.p_lt_1e-3.txt \
#  --colHashFile 1 --colScanFile 3 | awk 'BEGIN{OFS="\t"}{print $1,$2,$2,$3,$4,$5,$6,$7,$8,$9,$10}' \
#  > macrophage.naive_50kb.eigenMT.p_lt_1e-3.GRCh37.txt
submitJobs.py -j hashJoin_1000g_macrophage_qtl --MEM 20000 -q yesterday -c "bash ./do_hash_join.sh"


# Intersect the macrophage lead QTL SNP positions with the Blueprint QTL peaks
#echo -e "chr\tpos\tmacrophage_snp\tmacrophage_peak\tmacrophage_pval\tmacrophage_beta\tmonocyte_peak\tmonocyte_snp_chr\tmonocyte_snp_pos\tmonocyte_snp\tmonocyte_pval\tmonocyte_beta" \
#  > macrophage.overlap.mono_K27AC.txt
bedtools intersect -a macrophage.naive_50kb.eigenMT.p_lt_1e-3.GRCh37.txt -b  mono_K27AC.gene_minp.bed.gz -wa -wb \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$5,$7,$9,$15,$16,$17,$18,$19,$20}' \
  > macrophage.overlap.mono_K27AC.txt

#echo -e "chr\tpos\tmacrophage_snp\tmacrophage_peak_GRCh37\tmacrophage_pval\tmacrophage_beta\tmonocyte_peak\tmonocyte_snp_chr\tmonocyte_snp_pos\tmonocyte_snp\tmonocyte_pval\tmonocyte_beta" \
#  > macrophage.overlap.mono_K27AC.peaks.txt
hashJoin.pl --hashFile $JS/macrophage/GRCh37/ATAC_macrophage.peaks.GRCh37.bed --scanFile macrophage.overlap.mono_K27AC.txt \
  --colHashFile 4 --colScanFile 4 | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13":"$14"-"$15}' \
  | sort -k1,1n -k2,2n > macrophage.overlap.mono_K27AC.peaks.txt

echo -e "chr\tpos\tmacrophage_snp\tmacrophage_peak_GRCh37\tmacrophage_peak_GRCh38\tmacrophage_pval\tmacrophage_beta\tmonocyte_peak\tmonocyte_snp_chr\tmonocyte_snp_pos\tmonocyte_snp\tmonocyte_pval\tmonocyte_beta" \
  > macrophage.overlap.mono_K27AC.peaks2.txt
hashJoin.pl --hashFile $JS/macrophage/GRCh38/ATAC_macrophage.peaks.GRCh38.bed --scanFile macrophage.overlap.mono_K27AC.peaks.txt \
  --colHashFile 4 --colScanFile 4 | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$13,$14":"$15"-"$16,$5,$6,$7,$8,$9,$10,$11,$12}' \
  | sort -k1,1n -k2,2n >> macrophage.overlap.mono_K27AC.peaks2.txt
