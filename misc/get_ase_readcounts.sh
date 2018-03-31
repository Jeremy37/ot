#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
PATH=$JS/software/bcftools-1.6/bin:$PATH
#KOLF_VCF=$JS/ipsneurons/GRCh38/genotypes/kolf_2.imputed_phased.20150604.GRCh38.vcf.gz
#KOLF_VCF=$JS/ipsneurons/GRCh38/genotypes/kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.biallelic.snps.hets.chr.vcf.gz
KOLF_VCF=$JS/ipsneurons/GRCh38/genotypes/kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.biallelic.snps.chr.vcf.gz

cd $JS/ipsneurons/GRCh38/RNA

# It seems essential for ASEReadCounter that you use the same FASTA as was used
# for alignment of the BAM file, otherwise the program runs for hours and only
# dies at the end on some random chromosome (chrUn_KN707615v1_decoy), without
# outputting anything.
#GRCh38_FASTA=$JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa
GRCh38_FASTA=/lustre/scratch117/core/sciops_repository/references/Homo_sapiens/GRCh38_15_plus_hs38d1/all/fasta/Homo_sapiens.GRCh38_15_plus_hs38d1.fa

# First get ASE read counts for all samples
cut -f 1 irods.sample_lanes.txt | head -n 1 | submitJobs.py --MEM 8000 -j bamCountASE -q yesterday \
    -c "python ~/src/utils/bam/bamCountASE.py --indir $JS/ipsneurons/GRCh38/RNA --outdir $JS/ipsneurons/GRCh38/RNA --insuffix .bam --reference $GRCh38_FASTA --sites $KOLF_VCF --Xmx 7500m --execute True"

cut -f 1 irods.sample_lanes.txt | tail -n +2 | submitJobs.py --MEM 8000 -j bamCountASE \
    -c "python ~/src/utils/bam/bamCountASE.py --indir $JS/ipsneurons/GRCh38/RNA --outdir $JS/ipsneurons/GRCh38/RNA --insuffix .bam --reference $GRCh38_FASTA --sites $KOLF_VCF --Xmx 7500m --execute True"
grep "Successfully completed" FarmOut/bamCountASE.68*.txt


# Index the ASECounts files so that we can retrieve subsets
for samp in 486*; do
  bgzip -c $samp/$samp.ASEcounts > $samp/$samp.ASEcounts.gz
  tabix -s 1 -b 2 -e 2 -S 1 $samp/$samp.ASEcounts.gz
done

# Also create a file with ASE counts from all KOLF2 samples together
KOLF2_SAMPLES=( "4860STDY7079828" "4860STDY7079829" "4860STDY7079830" "4860STDY7079831" "4860STDY7079832" "4860STDY7079833" "4860STDY7079834" )
paste <(echo "Sample") <(head -n 1 4860STDY7079828/4860STDY7079828.ASEcounts) > analysis/ipsneurons.kolf2.ASEcounts.nochr
for samp in "${KOLF2_SAMPLES[@]}"; do
  sed '1d' $samp/$samp.ASEcounts | sed 's/^chr//' | perl -sne 'print "$samp\t".$_;' -- -samp=$samp >> analysis/ipsneurons.kolf2.ASEcounts.nochr
done
bgzip analysis/ipsneurons.kolf2.ASEcounts.nochr

# And the same for all samples including BOB iNeurons
paste <(echo "Sample") <(head -n 1 4860STDY7028458/4860STDY7028458.ASEcounts) > analysis/ipsneurons.ASEcounts.nochr
for samp in 486*; do
  sed '1d' $samp/$samp.ASEcounts | sed 's/^chr//' | perl -sne 'print "$samp\t".$_;' -- -samp=$samp >> analysis/ipsneurons.ASEcounts.nochr
done
bgzip analysis/ipsneurons.ASEcounts.nochr


# Get ASE counts for a specific gene across all samples
ENSEMBL_EXONS=$JS/reference/GRCh38/Homo_sapiens.GRCh38.91.exon_start_end.bed
GENE_NAME=PTK2B
GENEID=ENSG00000120899

#mkdir analysis/$GENE_NAME

# First get the exonic regions for the gene from the ensembl file
#cat $JS/reference/GRCh38/Homo_sapiens.GRCh38.91.chr.exon_start_end.bed \
#    | awk -v geneid="$GENEID" '$4 == geneid' | uniq > analysis/$GENE_NAME/tmp.exons.$GENEID.bed
#bedtools merge -i analysis/$GENE_NAME/tmp.exons.$GENEID.bed > analysis/$GENE_NAME/tmp.exons.$GENEID.merged.bed

#paste <(echo "Sample") <(gzhead 1 4860STDY7028458/4860STDY7028458.ASEcounts.gz) > analysis/$GENE_NAME/$GENE_NAME.ASEcounts.txt
#for samp in 486*; do
#  tabix $samp/$samp.ASEcounts.gz -B analysis/$GENE_NAME/tmp.exons.$GENEID.merged.bed | perl -sne 'print "$samp\t".$_;' -- -samp=$samp >> analysis/$GENE_NAME/$GENE_NAME.ASEcounts.txt
#done

# Get the relevant variants from the VCF file, since we need phase information
# to determine ASE levels
#bcftools view --regions-file analysis/$GENE_NAME/tmp.exons.$GENEID.merged.bed --genotype het --output-type z $KOLF_VCF > analysis/$GENE_NAME/$GENE_NAME.kolf_2.hets.vcf.gz
#rm analysis/$GENE_NAME/tmp.exons.*

cd $JS/ipsneurons/GRCh38/RNA/analysis

# No het SNPs in Kolf2
GENE_NAME=APOE
GENE_NAME=TREM2
GENE_NAME=ABI3
GENE_NAME=SORL1
GENE_NAME=CLU
GENE_NAME=ADAM10
GENE_NAME=BIN1
GENE_NAME=CCDC6
GENE_NAME=INPP5D
GENE_NAME=EPHA1

# Low depth
GENE_NAME=ACE
GENE_NAME=SLC24A4
GENE_NAME=CR1
GENE_NAME=PLCG2

# Has ASE SNPs with reasonable depth
GENE_NAME=PTK2B
GENE_NAME=SPPL2A
GENE_NAME=ABCA7
GENE_NAME=CD2AP
GENE_NAME=PICALM



# Run it twice, to get the background ASE for all genes, either using all transcript SNPs
# or only coding transcript SNPs
Rscript $JS/src/misc/plotGeneASE.R --args asefile=ipsneurons.kolf2.ASEcounts.nochr.gz \
                                          vcf=../../genotypes/kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.biallelic.snps.hets.vcf.gz \
                                          gene="$GENE_NAME" \
                                          codingOnly="T" \
                                          sampleSet="ipsneurons" \
                                          allgeneAse="ase/ase.ipsneurons.expressed_genes.ase.txt" \
                                          sampleGenotypeMap=ipsneurons.sample_genotype_map.txt \
                                          sampleGroupsFile=ipsneurons.sampleGroups.txt \
                                          out=ase/ase

GENES=( APOE TREM2 ABI3 SORL1 CLU ADAM10 BIN1 CCDC6 INPP5D EPHA1 ACE SLC24A4 CR1 PLCG2 PTK2B SPPL2A ABCA7 CD2AP PICALM )
for gene in "${GENES[@]}"
do
  submitJobs.py --MEM 6000 -j plotGeneASE.${gene} -q normal -n 2 \
  -c "Rscript $JS/src/misc/plotGeneASE.R --args asefile=ipsneurons.kolf2.ASEcounts.nochr.gz \
                                          vcf=$JS/ipsneurons/GRCh38/genotypes/kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.biallelic.snps.hets.vcf.gz \
                                          gene=$gene \
                                          codingOnly="F" \
                                          sampleSet="ipsneurons" \
                                          allgeneAse="ase/ase.ipsneurons.expressed_genes.ase.txt" \
                                          sampleGenotypeMap=ipsneurons.sample_genotype_map.txt \
                                          sampleGroupsFile=ipsneurons.sampleGroups.txt \
                                          out=ase/ase.ipsneurons"
done


