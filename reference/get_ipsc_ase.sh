#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
GRCh37_FASTA=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/reference/GRCh37/hs37d5.fa

INPUT_VCF=$JS/ipsc/GRCh37/merged.hipsci.9samples.wec.gtarray.HumanCoreExome-12_v1_0.imputed_phased.20150604.genotypes.vcf.gz
IPSC_VCF=$JS/ipsc/GRCh37/imputed.9_samples.snps.biallelic.hets.INFO_08.vcf.gz

cd $JS/ipsc/GRCh37

# Make a multi-sample VCF by combining single-sample VCFs
SAMPLES=( HPSI0613i-nukw_1 HPSI0813i-guss_1 HPSI0613i-auim_2 HPSI0613i-febc_2 HPSI1213i-hehd_2 HPSI0214i-pelm_3 HPSI0314i-qaqx_1 HPSI1113i-podx_1 HPSI0114i-joxm_1 )
for samp in "${SAMPLES[@]}"
do
  find /lustre/scratch116/vr/projects/hipsci/releases/data/gtarray/imputed_vcf/ -name "*$samp*"
done
# Extract the relevant samples from HIPSCI VCFs

# Manually create a file with the relevant sample paths, and pass to bcftools merge
mkdir tmp
while read F; do
  #echo "1" | submitJobs.py --MEM 6000 -j extract.9samples.from.vcf -q normal -n 2 \
  #  -c "bcftools view --samples HPSI0613i-nukw_1,HPSI0813i-guss_1,HPSI0613i-auim_2,HPSI0613i-febc_2,HPSI1213i-hehd_2,HPSI0214i-pelm_3,HPSI0314i-qaqx_1,HPSI1113i-podx_1,HPSI0114i-joxm_1 -Oz \
  #     /lustre/scratch116/vr/projects/hipsci/releases/data/gtarray/releases/merged_files/REL-2018-01/$F > tmp/$F.9samples.gz"
  tabix -p vcf tmp/$F.9samples.gz
done < hipsci_vcf_files.perchr.txt

cd tmp
bcftools concat -Oz --file-list hipsci_vcf_files.tmp.perchr.txt > ../merged.hipsci.9samples.wec.gtarray.HumanCoreExome-12_v1_0.imputed_phased.20150604.genotypes.vcf.gz

#bcftools merge -Oz --file-list hipsci_vcf_files.9samples.txt > hipsci.9samples.wec.gtarray.HumanCoreExome-12_v1_0.imputed_phased.20150604.genotypes.vcf.gz
#submitJobs.py --MEM 13000 -j subsetVCF -q yesterday -n 2 \
#    -c "bcftools norm --multiallelics + --fasta-ref $GRCh37_FASTA -Ou $INPUT_VCF \
#    | bcftools view -m2 -M2 -v snps --genotype het -Ou \
#    | bcftools filter -i 'IMP2[1] >= 0.8' \
#    | bcftools sort -Oz > new.imputed.9_samples.snps.biallelic.hets.INFO_08.vcf.gz"

submitJobs.py --MEM 13000 -j subsetVCF -q yesterday -n 2 \
    -c "bcftools norm --multiallelics + --fasta-ref $GRCh37_FASTA -Ou $INPUT_VCF \
    | bcftools view -m2 -M2 -v snps --genotype het -Ou \
    | bcftools filter -i 'INFO >= 0.8' \
    | bcftools sort -Oz > new.imputed.9_samples.snps.biallelic.hets.INFO_08.vcf.gz"
tabix -p vcf new.imputed.9_samples.snps.biallelic.hets.INFO_08.vcf.gz

submitJobs.py --MEM 2000 -j subsetVCF.togenes -q yesterday -n 2 \
    -c "bcftools view --regions-file $JS/reference/GRCh37/ensembldb.91.grch37.exons.bed.gz -Ou new.imputed.9_samples.snps.biallelic.hets.INFO_08.vcf.gz \
    | bcftools sort -Oz > new.imputed.9_samples.snps.biallelic.hets.INFO_08.exons.vcf.gz"
tabix -p vcf new.imputed.9_samples.snps.biallelic.hets.INFO_08.exons.vcf.gz



cd $JS/ipsc/GRCh37/RNA

# Make directories with symbolic links to the HIPSCI RNA-seq BAM files, so that
# we can run ASEReadCounter and save the results here
BAMDIR=/lustre/scratch116/vr/projects/hipsci/releases/data/rnaseq/realigned_bams
SAMPLES=( HPSI0613i-nukw_1 HPSI0813i-guss_1 HPSI0613i-auim_2 HPSI0613i-febc_2 )
for samp in "${SAMPLES[@]}"
do
    mkdir $samp
    ln -s $BAMDIR/EGAS00001000593/$samp/$samp.hs37d5.bwa.realigned.star.markdup.rnaseq.20150415.bam $samp/$samp.bam
    ln -s $BAMDIR/EGAS00001000593/$samp/$samp.hs37d5.bwa.realigned.star.markdup.rnaseq.20150415.bam.bai $samp/$samp.bam.bai
done

SAMPLES=( HPSI1213i-hehd_2 HPSI0214i-pelm_3 HPSI0314i-qaqx_1 HPSI1113i-podx_1 HPSI0114i-joxm_1 )
for samp in "${SAMPLES[@]}"
do
    mkdir $samp
    ln -s $BAMDIR/ERP007111/$samp/$samp.hs37d5.bwa.realigned.star.markdup.rnaseq.20150415.bam $samp/$samp.bam
    ln -s $BAMDIR/ERP007111/$samp/$samp.hs37d5.bwa.realigned.star.markdup.rnaseq.20150415.bam.bai $samp/$samp.bam.bai
done


# Need the reference fasta, with a sequence .dict file, for GATK
submitJobs.py --MEM 13000 -j createDict.hs37d5 -q yesterday -n 2 \
  -c "/software/jre1.8.0_74/bin/java -jar $JS/software/picard/picard.jar CreateSequenceDictionary REFERENCE=$JS/reference/GRCh37/hs37d5.fa OUTPUT=$JS/reference/GRCh37/hs37d5.dict"

while read samp; do
    echo $samp | submitJobs.py --MEM 6000 -j bamCountASE \
        -c "python ~/src/utils/bam/bamCountASE.py --indir $JS/ipsc/GRCh37/RNA --outdir $JS/ipsc/GRCh37/RNA \
        --insuffix .bam --reference $GRCh37_FASTA --sites $JS/ipsc/GRCh37/new.imputed.9_samples.snps.biallelic.hets.INFO_08.vcf.gz \
        --Xmx 2500m --execute True"
done < ipsc.sample_names.txt


RNADIR=$JS/ipsc/GRCh37/RNA
paste <(echo "Sample") <(head -n 1 $RNADIR/HPSI1213i-hehd_2/HPSI1213i-hehd_2.ASEcounts) > ase/ipsc.ASEcounts.9samples
while read samp; do
  sed '1d' $RNADIR/$samp/$samp.ASEcounts | perl -sne 'print "$samp\t".$_;' -- -samp=$samp >> ase/ipsc.ASEcounts.9samples
done < ipsc.sample_names.txt
bgzip ase/ipsc.ASEcounts.9samples

cd $RNADIR
gene=PTK2B
submitJobs.py --MEM 20000 -j plotGeneASE.${gene} -q yesterday -n 2 \
  -c "Rscript $JS/src/misc/plotGeneASE.R --args asefile=ase/ipsc.ASEcounts.9samples.gz \
                                          vcf=$JS/ipsc/GRCh37/new.imputed.9_samples.snps.biallelic.hets.INFO_08.exons.vcf.gz \
                                          gene=${gene} \
                                          codingOnly=F \
                                          grch37=T \
                                          sampleSet=ipsc \
                                          allgeneAse=ase/ipsc.expressed_genes.9.ase.txt \
                                          sampleGenotypeMap=ipsc.sample_genotype_map.txt \
                                          sampleGroupsFile=ipsc.sample_groups.txt \
                                          out=ase/ase.ipsc.9samples"
                                          
GENES=( APOE TREM2 ABI3 SORL1 CLU ADAM10 BIN1 CCDC6 INPP5D EPHA1 ACE SLC24A4 CR1 PLCG2 PTK2B SPPL2A ABCA7 CD2AP PICALM )
for gene in "${GENES[@]}"
do
  submitJobs.py --MEM 20000 -j plotGeneASE.${gene} -q yesterday -n 2 \
    -c "Rscript $JS/src/misc/plotGeneASE.R --args asefile=ase/ipsc.ASEcounts.9samples.gz \
                                          vcf=$JS/ipsc/GRCh37/new.imputed.9_samples.snps.biallelic.hets.INFO_08.exons.vcf.gz \
                                          gene=${gene} \
                                          codingOnly=F \
                                          grch37=T \
                                          sampleSet=ipsc \
                                          allgeneAse=ase/ipsc.expressed_genes.9.ase.txt \
                                          sampleGenotypeMap=ipsc.sample_genotype_map.txt \
                                          sampleGroupsFile=ipsc.sample_groups.txt \
                                          out=ase/ase.ipsc.9samples"
done
