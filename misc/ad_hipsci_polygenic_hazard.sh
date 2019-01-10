JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
HIPSCI_ROOT=/lustre/scratch116/vr/projects/hipsci/releases/data/gtarray/imputed_vcf

#cd $JS
#md polygenic_hazard
cd $JS/polygenic_hazard


grep -wFf phs.snps.txt <(zcat $DBSNP) > phs.snps.dbsnp.txt


KOLF2_VCF_GRCH37=$HIPSCI_ROOT/PRJEB11752/HPSI0114i-kolf_2/HPSI0114i-kolf_2.wec.gtarray.HumanCoreExome-12_v1_0.imputed_phased.20150604.genotypes.vcf.gz

# Get a list of all HIPSCI VCF files
find $HIPSCI_ROOT -name "*.vcf.gz" > hipsci.vcf.files.txt
wc -l hipsci.vcf.files.txt
#1953 hipsci.vcf.files.txt

# Extract the polygenic hazard score SNPs from the HIPSCI VCFs
mkdir tmp
i=1
while read hipsci_vcf; do
  echo "$i: $hipsci_vcf"
  i=$(($i + 1))
  vcfname=$(basename "$hipsci_vcf")
  sampleid=`echo $vcfname | perl -ne '@a=split(/\./); print $a[0]'`
  echo $sampleid
  submitJobs.py -j $sampleid.grep -q normal --nostdin -c "zcat $hipsci_vcf | grep -wFf phs.snps.txt > tmp/$vcfname.snps.txt"
done < hipsci.vcf.files.txt
#done < <(head -n 5 hipsci.vcf.files.txt)

grep "Successfully completed" FarmOut/*.grep.*.txt | wc -l
grep -iP "Failed|TERM|error|exit code" FarmOut/*.grep.*.txt | wc -l


# First get the start of the merged file for KOLF2, but leave out the kolf2
# genotypes since these will be added in the loop later
# Efficient way to find specific SNP IDs in a file
grep -wFf phs.snps.txt <(zcat $KOLF2_VCF_GRCH37) | cut -f 1-9 > hipsci.merged.phs.snps.vcf
echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" > hipsci.merged.phs.header.txt

i=1
while read hipsci_vcf; do
  i=$(($i + 1))
  vcfname=$(basename "$hipsci_vcf")
  sampleid=`echo $vcfname | perl -ne '@a=split(/\./); print $a[0]'`
  #grep -wFf phs.snps.txt <(zcat $hipsci_vcf) > tmp/$vcfname.snps.txt
  
  cp hipsci.merged.phs.snps.vcf hipsci.merged.phs.snps.vcf2
  cp hipsci.merged.phs.header.txt hipsci.merged.phs.header.txt2
  paste hipsci.merged.phs.snps.vcf2 <(cut -f 10 tmp/$vcfname.snps.txt) > hipsci.merged.phs.snps.vcf
  paste hipsci.merged.phs.header.txt2 <(echo $sampleid) > hipsci.merged.phs.header.txt
done < hipsci.vcf.files.txt
#done < <(head -n 50 hipsci.vcf.files.txt)


cp hipsci.merged.phs.snps.vcf hipsci.merged.phs.snps.vcf2
cat hipsci.merged.phs.header.txt hipsci.merged.phs.snps.vcf2 > hipsci.merged.phs.snps.vcf

