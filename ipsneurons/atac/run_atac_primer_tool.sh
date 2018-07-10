#!/bin/bash

JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
ATAC=$JS/ipsneurons/GRCh38/ATAC
DIR=$JS/ipsneurons/GRCh38/ATAC/analysis/atac_primer_tool
cd $DIR

mkdir input_bams
cd input_bams
mkdir npc neuron ineuron ipsc
END=bam
END=bam.bai
ln -s $ATAC/4859STDY7028449/4859STDY7028449.$END ineuron/BOB_iN_D9_ATAC_1.$END
ln -s $ATAC/4859STDY7028450/4859STDY7028450.$END ineuron/BOB_iN_D9_ATAC_2.$END
ln -s $ATAC/4859STDY7028451/4859STDY7028451.$END ineuron/BOB_iN_D9_ATAC_3.$END
ln -s $ATAC/4859STDY7028452/4859STDY7028452.$END ineuron/BOB_iN_D11_ATAC_1.$END
ln -s $ATAC/4859STDY7028453/4859STDY7028453.$END ineuron/BOB_iN_D11_ATAC_2.$END
ln -s $ATAC/4859STDY7028454/4859STDY7028454.$END ineuron/BOB_iN_D11_ATAC_3.$END

ln -s $ATAC/4859STDY7028455/4859STDY7028455.$END ipsc/Kolf2_iPS_ATAC_1.$END
ln -s $ATAC/4859STDY7028456/4859STDY7028456.$END ipsc/Kolf2_iPS_ATAC_2.$END
ln -s $ATAC/4859STDY7028457/4859STDY7028457.$END ipsc/Kolf2_iPS_ATAC_3.$END

ln -s $ATAC/4859STDY7079821/4859STDY7079821.$END npc/Kolf2_NPC_ATAC_1.$END
ln -s $ATAC/4859STDY7079822/4859STDY7079822.$END npc/Kolf2_NPC_ATAC_2.$END
ln -s $ATAC/4859STDY7079823/4859STDY7079823.$END npc/Kolf2_NPC_ATAC_3.$END

ln -s $ATAC/4859STDY7079824/4859STDY7079824.$END neuron/Kolf2_Neu_D35_ATAC_1.$END
ln -s $ATAC/4859STDY7079825/4859STDY7079825.$END neuron/Kolf2_Neu_D35_ATAC_2.$END
ln -s $ATAC/4859STDY7079826/4859STDY7079826.$END neuron/Kolf2_Neu_D35_ATAC_3.$END
ln -s $ATAC/4859STDY7079827/4859STDY7079827.$END neuron/Kolf2_Neu_D36_ATAC_4.$END

mkdir all
cp -a npc/* all; cp -a neuron/* all; cp -a ineuron/* all; cp -a ipsc/* all

# ATAC Primer Tool expects peaks file to be unzipped and have exactly 4 cols
zcat $ATAC/peaks/atac_ipsc_peaks.narrowPeak.gz | cut -f 1-4 > $ATAC/peaks/atac_ipsc_peaks.narrowPeak.bed
zcat $ATAC/peaks/atac_npc_peaks.narrowPeak.gz | cut -f 1-4 > $ATAC/peaks/atac_npc_peaks.narrowPeak.bed
zcat $ATAC/peaks/atac_ineuron_peaks.narrowPeak.gz | cut -f 1-4 > $ATAC/peaks/atac_ineuron_peaks.narrowPeak.bed
zcat $ATAC/peaks/atac_neuron_peaks.narrowPeak.gz | cut -f 1-4 > $ATAC/peaks/atac_neuron_peaks.narrowPeak.bed

zcat $ATAC/peaks/atac_ipsc_peaks.narrowPeak.gz | cut -f 1-4 | head -n 100 > $ATAC/peaks/atac_ipsc_peaks.narrowPeak.head100.bed


cd $DIR
################################################################################
# Run 1
# These runs didn't work, since APT is way too slow when given a large number
# of peaks, due to badly written R code
python $JS/software/ATACPrimerTool/pipelines/ATACPrimerTool.py -O $DIR -S APT_IPSC -G hg38 -I input_bams/ipsc -B $ATAC/peaks/atac_ipsc_peaks.narrowPeak.head5.bed --genome_fasta $JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa

submitJobs.py --MEM 4000 -j APT_ipsc -n 2 -q yesterday    -c "python $JS/software/ATACPrimerTool/pipelines/ATACPrimerTool.py -O $DIR -S APT_ipsc -G hg38 -I input_bams/ipsc -B $ATAC/peaks/atac_ipsc_peaks.narrowPeak.head100.bed --genome_fasta $JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa"

submitJobs.py --MEM 4000 -j APT_ipsc -n 2 -q normal    -c "python $JS/software/ATACPrimerTool/pipelines/ATACPrimerTool.py -O $DIR -S APT_ipsc -G hg38 -I input_bams/ipsc -B $ATAC/peaks/atac_ipsc_peaks.narrowPeak.bed --genome_fasta $JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa"
submitJobs.py --MEM 4000 -j APT_ineuron -n 2 -q normal -c "python $JS/software/ATACPrimerTool/pipelines/ATACPrimerTool.py -O $DIR -S APT_ineuron -G hg38 -I input_bams/ineuron -B $ATAC/peaks/atac_ineuron_peaks.narrowPeak.bed --genome_fasta $JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa"
submitJobs.py --MEM 4000 -j APT_neuron -n 2 -q normal  -c "python $JS/software/ATACPrimerTool/pipelines/ATACPrimerTool.py -O $DIR -S APT_neuron -G hg38 -I input_bams/neuron -B $ATAC/peaks/atac_neuron_peaks.narrowPeak.bed --genome_fasta $JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa"
submitJobs.py --MEM 4000 -j APT_npc -n 2 -q normal     -c "python $JS/software/ATACPrimerTool/pipelines/ATACPrimerTool.py -O $DIR -S APT_npc -G hg38 -I input_bams/npc -B $ATAC/peaks/atac_npc_peaks.narrowPeak.bed --genome_fasta $JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa"

# Focusing on a small number of peaks where we have QTLs worked fine
# I made the BED file sn_qtl_in_ipsc_atac_peaks.bed by looking at the Google
# sheet I made for Sarah, "Sensory neuron atac QTLs overlapping iPSC atac peaks",
# and selecting the top peaks where a SNP had PPA > 0.5. I filtered out peaks
# with low R2 for fSNP or rSNP (< 0.8).
submitJobs.py --MEM 4000 -j APT_sn_qtl_in_ipsc_peak -n 2 -q normal \
  -c "python $JS/software/ATACPrimerTool/pipelines/ATACPrimerTool.py -O $DIR -S APT_sn_qtl_in_ipsc_peak -G hg38 -I input_bams/ipsc -B sn_qtl_in_ipsc_atac_peaks.bed --genome_fasta $JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa"


################################################################################
# Run 2
# Here I am re-running APT giving all our BAMs as input so that it can better
# estimate the correlation between # spanning reads and peak size.
submitJobs.py --MEM 4000 -j APT_sn_qtl_in_ipsc_peak -n 2 -q normal \
  -c "python $JS/software/ATACPrimerTool/pipelines/ATACPrimerTool.py -O $DIR -S APT_sn_qtl_in_ipsc_peak -G hg38 -I input_bams/all -B sn_qtl_in_ipsc_atac_peaks.bed --genome_fasta $JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa"

# Here we do the same run, but giving windows of 200 bp around each QTL SNP
submitJobs.py --MEM 4000 -j APT_sn_qtl_in_ipsc_peak_win200 -n 2 -q normal \
  -c "python $JS/software/ATACPrimerTool/pipelines/ATACPrimerTool.py -O $DIR -S APT_sn_qtl_in_ipsc_peak_win200 -G hg38 -I input_bams/all -B sn_qtl_in_ipsc_atac_peaks.win200.bed --genome_fasta $JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa"

# Add back in the information from my full file annotating these peaks
# First we need to cut column 5 from the file into two
cat APT_sn_qtl_in_ipsc_peak/qPCR_regions_corr.txt | perl -ane '@X=split(/_window/, $F[4]); print join("\t", @F[1..3], $X[0], @F[4..6])."\n";' > APT_sn_qtl_in_ipsc_peak/qPCR_regions_corr.split.txt
cat APT_sn_qtl_in_ipsc_peak_win200/qPCR_regions_corr.txt | perl -ane '@X=split(/_window/, $F[4]); print join("\t", @F[1..3], $X[0], @F[4..6])."\n";' > APT_sn_qtl_in_ipsc_peak_win200/qPCR_regions_corr.split.txt

hashJoin.pl --hashFile sn_qtl_in_ipsc_atac_peaks.full.txt --scanFile APT_sn_qtl_in_ipsc_peak/qPCR_regions_corr.split.txt --colHashFile 4 --colScanFile 4 --outerjoin --header --insep $'\t' > sn_qtl_in_ipsc_atac_peaks.output.ann.txt
hashJoin.pl --hashFile sn_qtl_in_ipsc_atac_peaks.win200.full.txt --scanFile APT_sn_qtl_in_ipsc_peak_win200/qPCR_regions_corr.split.txt --colHashFile 4 --colScanFile 4 --outerjoin --header --insep $'\t' > sn_qtl_in_ipsc_atac_peaks.win200.output.ann.txt


################################################################################
# Genotypes
# Get genotypes for SNPs of interest in KOLF_2

sed '1d' sn_qtl_in_ipsc_atac_peaks.win200.full.txt | cut -f 5 > sn_qtl_in_ipsc_atac_peaks.snps.txt
zcat $JS/ipsneurons/GRCh38/genotypes/kolf_2.imputed_phased.20150604.GRCh38.vcf.gz | grep -wF -f sn_qtl_in_ipsc_atac_peaks.snps.txt > sn_qtl_in_ipsc_atac_peaks.kolf2_vcf.txt
echo -e "snpid\tkolf2_genotype" > sn_qtl_in_ipsc_atac_peaks.kolf2_genotype_map.txt
cat sn_qtl_in_ipsc_atac_peaks.kolf2_vcf.txt | perl -ane '@GT=split(/:/, $F[9]); print join("\t", $F[2],$GT[0])."\n";' >> sn_qtl_in_ipsc_atac_peaks.kolf2_genotype_map.txt

hashJoin.pl --hashFile sn_qtl_in_ipsc_atac_peaks.kolf2_genotype_map.txt --scanFile win200.tmp.txt --colHashFile 1 --colScanFile 9 --outerjoin --header --insep $'\t' > win200.new.txt
