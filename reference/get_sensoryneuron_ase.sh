#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
STARDIR=/warehouse/compgen_wh05/js29/sn_hg38/STAR

cd $JS/sensoryneurons/GRCh38/RNA

GRCh38_FASTA=/lustre/scratch117/cellgen/team170/ka8/annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa
SN_VCF=/warehouse/compgen_wh05/js29/sensoryneurons/rasqual/input/imputed.97_samples.snps_indels.INFO_08.MAF_0.05.RASQUAL.vcf.gz

cd $JS/sensoryneurons/GRCh38

#HPSI1213i-hehd_2
#HPSI0613i-nukw_1
#HPSI0214i-pelm_3
#HPSI0813i-guss_1
#HPSI1113i-podx_1
#HPSI0314i-qaqx_1

submitJobs.py --MEM 13000 -j subsetVCF -q yesterday -n 2 \
    -c "bcftools norm --multiallelics + --fasta-ref $GRCh38_FASTA -Ou $SN_VCF \
    | bcftools view --samples HPSI1213i-hehd_2,HPSI0613i-nukw_1,HPSI0214i-pelm_3 -m2 -M2 -v snps --genotype het -Ou \
    | bcftools sort -Oz > imputed.3_samples.snps.biallelic.hets.INFO_08.MAF_0.05.RASQUAL.vcf.gz"
tabix -p vcf imputed.3_samples.snps.biallelic.hets.INFO_08.MAF_0.05.RASQUAL.vcf.gz

submitJobs.py --MEM 13000 -j subsetVCF -q yesterday -n 2 \
    -c "bcftools norm --multiallelics + --fasta-ref $GRCh38_FASTA -Ou $SN_VCF \
    | bcftools view --samples HPSI1213i-hehd_2,HPSI0613i-nukw_1,HPSI0214i-pelm_3,HPSI0813i-guss_1,HPSI1113i-podx_1,HPSI0314i-qaqx_1,HPSI0613i-auim_2,HPSI0613i-febc_2,HPSI0114i-joxm_1 -m2 -M2 -v snps --genotype het -Ou \
    | bcftools sort -Oz > imputed.9_samples.snps.biallelic.hets.INFO_08.MAF_0.05.RASQUAL.vcf.gz"
tabix -p vcf imputed.9_samples.snps.biallelic.hets.INFO_08.MAF_0.05.RASQUAL.vcf.gz

submitJobs.py --MEM 2000 -j subsetVCF.togenes -q yesterday -n 2 \
    -c "bcftools view --regions-file $JS/reference/GRCh38/ensembldb.91.exons.bed.gz -Ou imputed.9_samples.snps.biallelic.hets.INFO_08.MAF_0.05.RASQUAL.vcf.gz \
    | bcftools sort -Oz > imputed.9_samples.snps.biallelic.hets.INFO_08.MAF_0.05.RASQUAL.exons.vcf.gz"
tabix -p vcf imputed.9_samples.snps.biallelic.hets.INFO_08.MAF_0.05.RASQUAL.exons.vcf.gz




cd $JS/sensoryneurons/GRCh38/RNA

#while read samp; do
#    echo $samp | submitJobs.py --MEM 4000 -j bamCountASE \
#        -c "python ~/src/utils/bam/bamCountASE.py --indir $JS/sensoryneurons/GRCh38/RNA --outdir $JS/sensoryneurons/GRCh38/RNA \
#        --insuffix .Aligned.sortedByCoord.out.bam --reference $GRCh38_FASTA --sites $JS/sensoryneurons/GRCh38/imputed.3_samples.snps.biallelic.hets.INFO_08.MAF_0.05.RASQUAL.vcf.gz \
#        --Xmx 2500m --execute True"
#done < sn.sample_names.txt


paste <(echo "Sample") <(head -n 1 $STARDIR/hehd_2a/hehd_2a.ASEcounts) > ase/sensoryneuron.ASEcounts.3samples
while read samp; do
  sed '1d' $STARDIR/$samp/$samp.ASEcounts | perl -sne 'print "$samp\t".$_;' -- -samp=$samp >> ase/sensoryneuron.ASEcounts.3samples
done < sensoryneuron.sample_names.3.txt
bgzip ase/sensoryneuron.ASEcounts.3samples

paste <(echo "Sample") <(head -n 1 $STARDIR/hehd_2a/hehd_2a.ASEcounts) > ase/sensoryneuron.ASEcounts.9samples
while read samp; do
  sed '1d' $STARDIR/$samp/$samp.ASEcounts | perl -sne 'print "$samp\t".$_;' -- -samp=$samp >> ase/sensoryneuron.ASEcounts.9samples
done < sensoryneuron.sample_names.9.txt
bgzip ase/sensoryneuron.ASEcounts.9samples


cd $JS/sensoryneurons/GRCh38/RNA
gene=PTK2B
submitJobs.py --MEM 20000 -j plotGeneASE.9samples.${gene} -q yesterday -n 2 \
    -c "Rscript $JS/src/misc/plotGeneASE.R --args asefile=ase/sensoryneuron.ASEcounts.9samples.gz \
                                          vcf=$JS/sensoryneurons/GRCh38/imputed.9_samples.snps.biallelic.hets.INFO_08.MAF_0.05.RASQUAL.vcf.gz \
                                          gene=${gene} \
                                          codingOnly=F \
                                          sampleSet=sensoryneuron \
                                          allgeneAse=ase/sensoryneuron.expressed_genes.9.ase.txt \
                                          sampleGenotypeMap=sensoryneuron.sample_genotype_map.txt \
                                          sampleGroupsFile=sensoryneuron.sample_groups.txt \
                                          out=ase/ase.sensoryneuron.9samples"
                                          
GENES=( APOE TREM2 ABI3 SORL1 CLU ADAM10 BIN1 CCDC6 INPP5D EPHA1 ACE SLC24A4 CR1 PLCG2 PTK2B SPPL2A ABCA7 CD2AP PICALM )
for gene in "${GENES[@]}"
do
  submitJobs.py --MEM 6000 -j plotGeneASE.${gene} -q normal -n 2 \
    -c "Rscript $JS/src/misc/plotGeneASE.R --args asefile=ase/sensoryneuron.ASEcounts.9samples.gz \
                                          vcf=$JS/sensoryneurons/GRCh38/imputed.9_samples.snps.biallelic.hets.INFO_08.MAF_0.05.RASQUAL.exons.vcf.gz \
                                          gene=${gene} \
                                          codingOnly=F \
                                          sampleSet=sensoryneuron \
                                          allgeneAse=ase/sensoryneuron.expressed_genes.9.ase.txt \
                                          sampleGenotypeMap=sensoryneuron.sample_genotype_map.txt \
                                          sampleGroupsFile=sensoryneuron.sample_groups.txt \
                                          out=ase/ase.sensoryneuron.9samples"
done
