#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
KAUR_STAR=/lustre/scratch117/cellgen/team170/ka8/projects/macrophage-trQTLs/processed/salmonella/STAR

cd $JS/macrophage/GRCh38/STAR

#SAMPLES=( hehd_A hehd_B hehd_C hehd_D nukw_A nukw_B nukw_C nukw_D pelm_A pelm_B pelm_C pelm_D )

while read samp; do
    mkdir $samp
    ln -s $KAUR_STAR/$samp/$samp.Aligned.sortedByCoord.out.bam $samp/$samp.Aligned.sortedByCoord.out.bam
    ln -s $KAUR_STAR/$samp/$samp.Aligned.sortedByCoord.out.bam.bai $samp/$samp.Aligned.sortedByCoord.out.bam.bai
done < macrophage.sample_names.9.txt

GRCh38_FASTA=/lustre/scratch117/cellgen/team170/ka8/annotations/GRCh38/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa
MACROPHAGE_VCF=/lustre/scratch117/cellgen/team170/ka8/projects/macrophage-gxe-study/genotypes/SL1344/imputed_20151005/imputed.86_samples.snps_only.vcf.gz

cd $JS/macrophage/GRCh38

submitJobs.py --MEM 13000 -j subsetVCF -q yesterday -n 2 \
    -c "bcftools norm --multiallelics + --fasta-ref $GRCh38_FASTA -Ou $MACROPHAGE_VCF \
    | bcftools view -m2 -M2 -v snps --genotype het -Ou \
    | bcftools filter -i 'INFO >= 0.8' \
    | bcftools sort | bgzip > imputed.86_samples.biallelic.snps.hets.vcf.gz"

submitJobs.py --MEM 13000 -j subsetVCF -q yesterday -n 2 \
    -c "bcftools view --samples HPSI1213i-hehd_2,HPSI0613i-nukw_1,HPSI0214i-pelm_3 --genotype het -Oz imputed.86_samples.biallelic.snps.hets.vcf.gz \
    > imputed.3samples.biallelic.snps.hets.vcf.gz"
tabix -p vcf imputed.3samples.biallelic.snps.hets.vcf.gz

submitJobs.py --MEM 13000 -j subsetVCF -q yesterday -n 2 \
    -c "bcftools view --samples HPSI1213i-hehd_2,HPSI0613i-nukw_1,HPSI0214i-pelm_3,HPSI0813i-guss_1,HPSI1113i-podx_1,HPSI0314i-qaqx_1,HPSI0613i-auim_2,HPSI0613i-febc_2,HPSI0114i-joxm_1 --genotype het -Oz imputed.86_samples.biallelic.snps.hets.vcf.gz \
    > imputed.9samples.biallelic.snps.hets.vcf.gz"
tabix -p vcf imputed.9samples.biallelic.snps.hets.vcf.gz

submitJobs.py --MEM 2000 -j subsetVCF.togenes -q yesterday -n 2 \
    -c "bcftools view --regions-file $JS/reference/GRCh38/ensembldb.91.exons.bed.gz -Ou imputed.9samples.biallelic.snps.hets.vcf.gz \
    | bcftools sort -Oz > imputed.9samples.biallelic.snps.hets.exons.vcf.gz"
tabix -p vcf imputed.9samples.biallelic.snps.hets.exons.vcf.gz


cd $JS/macrophage/GRCh38/STAR

while read samp; do
    echo $samp | submitJobs.py --MEM 4000 -j bamCountASE \
        -c "python ~/src/utils/bam/bamCountASE.py --indir $JS/macrophage/GRCh38/STAR --outdir $JS/macrophage/GRCh38/STAR \
        --insuffix .Aligned.sortedByCoord.out.bam --reference $GRCh38_FASTA --sites $JS/macrophage/GRCh38/imputed.9samples.biallelic.snps.hets.vcf.gz --Xmx 2500m --execute True"
done < macrophage.sample_names.9.txt

# Combine ASE counts from multiple samples together
cd $JS/macrophage/GRCh38/STAR
paste <(echo "Sample") <(head -n 1 hehd_A/hehd_A.ASEcounts) > analysis/ase/macrophage.ASEcounts
while read samp; do
  sed '1d' $samp/$samp.ASEcounts | hperl -sne 'print "$samp\t".$_;' -- -samp=$samp >> analysis/ase/macrophage.ASEcounts
done < macrophage.sample_names.txt
bgzip analysis/ase/macrophage.ASEcounts

SAMPLE_SETS=( naive ifng SL ifng_SL )
for sampset in "${SAMPLE_SETS[@]}"
do
  paste <(echo "Sample") <(head -n 1 hehd_A/hehd_A.ASEcounts) > analysis/ase/macrophage.ASEcounts.$sampset
  while read samp; do
    echo $samp
    sed '1d' $samp/$samp.ASEcounts | perl -sne 'print "$samp\t".$_;' -- -samp=$samp >> analysis/ase/macrophage.ASEcounts.$sampset
  done < macrophage.sample_names.$sampset.txt
  bgzip analysis/ase/macrophage.ASEcounts.$sampset
done


cd $JS/macrophage/GRCh38/STAR/analysis
gene=PTK2B
submitJobs.py --MEM 20000 -j plotGeneASE.${gene} -q yesterday -n 2 \
  -c "Rscript $JS/src/misc/plotGeneASE.R --args asefile=ase/macrophage.ASEcounts.3samples.gz \
                                          vcf=$JS/macrophage/GRCh38/imputed.3samples.biallelic.snps.hets.exons.vcf.gz \
                                          gene=${gene} \
                                          codingOnly=F \
                                          sampleSet=macrophage \
                                          allgeneAse=ase/ase.macrophage.expressed_genes.3samples.txt \
                                          sampleGenotypeMap=macrophage.sample_genotype_map.txt \
                                          sampleGroupsFile=macrophage.sample_groups.txt \
                                          out=ase/ase.3samples"


GENES=( APOE TREM2 ABI3 SORL1 CLU ADAM10 BIN1 CCDC6 INPP5D EPHA1 ACE SLC24A4 CR1 PLCG2 PTK2B SPPL2A ABCA7 CD2AP PICALM )
GENES=( PTK2B )
SAMPLE_SETS=( naive ifng SL ifng_SL )

for sampset in "${SAMPLE_SETS[@]}"
do

for gene in "${GENES[@]}"
do
  submitJobs.py --MEM 15000 -j plotGeneASE.$sampset.$gene -q normal -n 2 \
    -c "Rscript $JS/src/misc/plotGeneASE.R --args asefile=ase/macrophage.ASEcounts.$sampset.gz \
                                          vcf=$JS/macrophage/GRCh38/imputed.9samples.biallelic.snps.hets.exons.vcf.gz \
                                          gene=${gene} \
                                          codingOnly=F \
                                          sampleSet=macrophage \
                                          allgeneAse=ase/ase.macrophage.expressed_genes.$sampset.txt \
                                          sampleGenotypeMap=macrophage.sample_genotype_map.txt \
                                          sampleGroupsFile=macrophage.sample_groups.txt \
                                          out=ase/ase.macrophage.$sampset"
done # for gene

done # for sampset