JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

# Run CrossMap to lift the HIPSCI VCF from GRCh37 to GRCh38
HIPSCI_ROOT=/lustre/scratch116/vr/projects/hipsci/releases/data/gtarray/imputed_vcf
KOLF2_VCF_GRCH37=$HIPSCI_ROOT/PRJEB11752/HPSI0114i-kolf_2/HPSI0114i-kolf_2.wec.gtarray.HumanCoreExome-12_v1_0.imputed_phased.20150604.genotypes.vcf.gz

# Testing CrossMap
zcat $KOLF2_VCF_GRCH37 | grep -v "^#" | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3}' | head -n 100 > test.bed
zcat $KOLF2_VCF_GRCH37 | head -n 200 > test.vcf
CrossMap.py bed $JS/software/CrossMap/GRCh37_to_GRCh38.chain.gz test.bed

GRCh38_FASTA=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/reference/GRCh38/GRCh38.fa.gz
CrossMap.py vcf $JS/software/CrossMap/GRCh37_to_GRCh38.chain.gz test.vcf $GRCh38_FASTA test.out.vcf


# Run CrossMap
CrossMap.py vcf $JS/software/CrossMap/GRCh37_to_GRCh38.chain.gz <(zcat $KOLF2_VCF_GRCH37 | head -n 200) $GRCh38_FASTA kolf_2.imputed_phased.20150604.GRCh38.test.vcf

submitJobs.py --MEM 1000 -j CrossMap.kolf2.GRCh37_to_GRCh38 -q yesterday -n 2 \
    -c "CrossMap.py vcf $JS/software/CrossMap/GRCh37_to_GRCh38.chain.gz $KOLF2_VCF_GRCH37 $GRCh38_FASTA kolf_2.imputed_phased.20150604.GRCh38.vcf"

# This command fails with too many temporary files - I guess you can't sort a
# VCF file this big with bcftools
#submitJobs.py --MEM 20000 -j sortVCF -q yesterday -n 2 -c "bcftools view kolf_2.imputed_phased.20150604.GRCh38.vcf | bcftools sort -O z -m 18000 --temp-dir ./tmp > kolf_2.imputed_phased.20150604.GRCh38.vcf.gz"

submitJobs.py --MEM 1000 -j sortVCF -q yesterday -n 2 \
  -c '(cat kolf_2.imputed_phased.20150604.GRCh38.vcf | grep "^#"; cat kolf_2.imputed_phased.20150604.GRCh38.vcf | grep -v "^#" | sort -k1,1 -nk2,2) | bgzip > kolf_2.imputed_phased.20150604.GRCh38.vcf.gz'
rm kolf_2.imputed_phased.20150604.GRCh38.vcf
tabix -p vcf kolf_2.imputed_phased.20150604.GRCh38.vcf.gz

# Subset to SNPs with imputation INFO score >= 0.8
submitJobs.py --MEM 2000 -j filterINFO.0.8 -q yesterday -n 2 -c "bcftools view kolf_2.imputed_phased.20150604.GRCh38.vcf.gz | bcftools filter -i 'IMP2[1] >= 0.8' | bcftools sort -O z -m 4000 > kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.vcf.gz"
submitJobs.py --MEM 13000 -j sortVCF.INFO.0.8 -q yesterday -n 2 -c "bcftools view kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.vcf.gz | bcftools sort -O z -m 12000 > kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.sorted.vcf.gz"

mv kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.vcf.gz kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.vcf.gz.old
mv kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.sorted.vcf.gz kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.vcf.gz

tabix -p vcf kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.vcf.gz

# Subset to SNPs in gene exons
bcftools view kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.vcf.gz --regions-file $JS/reference/GRCh38/Homo_sapiens.GRCh38.89.exons.bed -O z > kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.exons.vcf.gz

###############################################################################
# Get VCF line for variant by position, e.g. rs356168 (GRCh38, 4:89753280, GRCh37, 4:90674431)
tabix kolf_2.imputed_phased.20150604.GRCh38.vcf.gz 4:89753280-89753280
# The variant is not present in KOLF2

