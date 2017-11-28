newgrp otcoregen # Set primary group to otcoregen so new files have this group by default
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
source $JS/software/ipsneurons.software.sh

###############################################################################
## ATAC-seq - get ASE counts
cd $SEQ/ATAC

# Select certain variants
#KOLF_VCF=$JS/ipsneurons/GRCh38/genotypes/kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.vcf.gz
KOLF_VCF=$JS/ipsneurons/GRCh38/genotypes/kolf_2.imputed_phased.20150604.GRCh38.vcf.gz
zcat $KOLF_VCF | head -n 1000 | grep '^#' > kolf.ase.vcf
zcat $KOLF_VCF | grep -wF -f $JS/ipsneurons/src/ase.snplist.txt >> kolf.ase.vcf


#Use ASEReadCounter to count allele-specific expression
cut -f1 sample_lists/sensoryneurons.samples.v5.new.txt | sed '1d' | submitJobs.py --MEM 4000 --jobname bamCountASE --queue normal --command "python ~/src/utils/bam/bamCountASE.py --indir $STARDIR --outdir $STARDIR --insuffix .Aligned.sortedByCoord.RG.bam --reference $GRCh38_FASTA --sites $VCF_ASE --Xmx 2500m --execute True"

#GRCh38_FASTA=$JS/reference/GRCh38/GRCh38.fa
#GRCh38_FASTA=$JS/reference/GRCh38/GRCh38.chr.fa
GRCh38_FASTA=$JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa


for samp in 485*; do
  samtools view -b $samp/$samp.bam "chr4:89724099-89838315" > $samp/$samp.SNCAregion.bam
  samtools index $samp/$samp.SNCAregion.bam
done

cut -f 1 irods.sample_lanes.txt | tail -n +2 | submitJobs.py --MEM 4000 -j bamCountASE \
    -c "python ~/src/utils/bam/bamCountASE.py --indir $JS/ipsneurons/GRCh38/ATAC --outdir $JS/ipsneurons/GRCh38/ATAC --insuffix .bam --reference $GRCh38_FASTA --sites analysis/kolf.ase.vcf --Xmx 2500m --execute True"

head -n 1 4859STDY7079824/4859STDY7079824.ASEcounts > analysis/ATAC.asecounts.txt
for samp in 485*; do
  paste <(echo "$samp") <(tail -n +2 $samp/$samp.ASEcounts) >> analysis/ATAC.asecounts.txt
done



###############################################################################
## RNA-seq - get ASE counts
cd $SEQ/RNA

GRCh38_FASTA=$JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa

cut -f 1 irods.sample_lanes.txt | tail -n +2 | submitJobs.py --MEM 4000 -j bamCountASE \
    -c "python ~/src/utils/bam/bamCountASE.py --indir $JS/ipsneurons/GRCh38/RNA --outdir $JS/ipsneurons/GRCh38/RNA --insuffix .bam --reference $GRCh38_FASTA --sites analysis/kolf.ase.vcf --Xmx 2500m --execute True"

head -n 1 4860STDY7028458/4860STDY7028458.ASEcounts > analysis/RNAseq.asecounts.txt
for samp in 486*; do
  paste <(echo "$samp") <(tail -n +2 $samp/$samp.ASEcounts) >> analysis/RNAseq.asecounts.txt
done

for samp in 486*; do
  #samtools view -b $samp/$samp.bam "chr4:89724099-89838315" > $samp/$samp.SNCAregion.bam
  samtools index $samp/$samp.SNCAregion.bam
done

