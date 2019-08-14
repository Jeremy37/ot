###############################################################################
# This script contains steps for comparing the AD GWAS meta-analysis results
# between our v2 and v3 analyses, as well as the Finemap results from each.
OT=/lustre/scratch115/realdata/mdt3/projects/otcoregen
JS=$OT/jeremys
ROOT=$JS/AD_finemap
STATSv2=$JS/AD_finemap_v1/summary_stats/AD.IGAP1_GWAX_v2.meta.bgz
STATSv3=$JS/AD_finemap/summary_stats/AD.IGAP1_GWAX_v3.meta.bgz
cd $ROOT

mkdir compare_versions

# First make a list of loci that are included in both meta-analysis versions
cp AD.gwsig_loci.tsv AD.compare_loci.tsv
# **** update this file manually to add columns lead_chrpos_1, lead_chrpos_2,
# which are used below to get the appropriate FINEMAP file in each version.
# 
cp AD.compare_loci.tsv AD.compare_loci.finemap.tsv
# Add columns lead_chrpos_1 and lead_chrpos_2, which are used below to get the
# appropriate FINEMAP file in each version

# Edit the AD.compare_loci.tsv file to set ncausal to 1 for all loci except
# PTK2B-CLU and ABCA7, since this is the number of SNPs used in the v2 analysis.
# But in AD.compare_loci.finemap.tsv, we compare different number of causals
# within the v3 meta-analysis, so we leave the new estimated number of causals.


################################################################################
# Make a table of all SNPs at our loci of interest to indicate whether they
# appear in each GWAS:
# Lambert, Kunkle, UKBB GWAX, Meta v2, Meta v3

zcat $OT/AD_PD_finemap/summary_stats/ADD.proxy_v2.bgen.stats.gz | bgzip > $JS/datasets/GWAS/AD.GWAX_v2.bgen.stats.gz
tabix -s 2 -b 3 -e 3 -S 1 $JS/datasets/GWAS/AD.GWAX_v2.bgen.stats.gz

body() {
    IFS= read -r header
    printf '%s\n' "$header"
    "$@"
}
zcat $JS/datasets/GWAS/Kunkle_etal_Stage1_results.txt.gz | tr ' ' '\t' | body sort -k1,1 -k2,2n | bgzip > $JS/datasets/GWAS/Kunkle_etal_Stage1_results.bgzip.txt.gz
tabix -s 1 -b 2 -e 2 -S 1 $JS/datasets/GWAS/Kunkle_etal_Stage1_results.bgzip.txt.gz

# We want to extract SNPs from each meta-analysis (v2 and v3) that are significant
# in either. So first we get a list of SNPs from each, merge together, and then
# extract all of these from each summary stats file.
python $JS/src/gwas/extract_locus_stats.py --locus_file AD.compare_loci.tsv --gwas_file $STATSv2 --window 500000 --out compare_versions/v2.meta.locus_
python $JS/src/gwas/extract_locus_stats.py --locus_file AD.compare_loci.tsv --gwas_file $STATSv3 --window 500000 --out compare_versions/v3.meta.locus_
python $JS/src/gwas/extract_locus_stats.py --locus_file AD.compare_loci.tsv --gwas_file $JS/datasets/GWAS/Alzheimers_disease_Lambert_2013.sorted.txt.gz --window 500000 --out compare_versions/Lambert.locus_
python $JS/src/gwas/extract_locus_stats.py --locus_file AD.compare_loci.tsv --gwas_file $JS/datasets/GWAS/Kunkle_etal_Stage1_results.bgzip.txt.gz --window 500000 --out compare_versions/Kunkle.locus_
python $JS/src/gwas/extract_locus_stats.py --locus_file AD.compare_loci.tsv --gwas_file $JS/datasets/GWAS/AD.GWAX_v2.bgen.stats.gz --window 500000 --out compare_versions/gwax.locus_

(gzhead 1 $STATSv2; zcat compare_versions/v2.meta.locus_*.gz) | bgzip > compare_versions/v2.meta.loci.gz
(gzhead 1 $STATSv3; zcat compare_versions/v3.meta.locus_*.gz) | bgzip > compare_versions/v3.meta.loci.gz
(gzhead 1 $JS/datasets/GWAS/Alzheimers_disease_Lambert_2013.sorted.txt.gz; zcat compare_versions/Lambert.locus_*.gz) | bgzip > compare_versions/Lambert.loci.gz
(gzhead 1 $JS/datasets/GWAS/Kunkle_etal_Stage1_results.bgzip.txt.gz; zcat compare_versions/Kunkle.locus_*.gz) | bgzip > compare_versions/Kunkle.loci.gz
(gzhead 1 $JS/datasets/GWAS/AD.GWAX_v2.bgen.stats.gz; zcat compare_versions/gwax.locus_*.gz) | bgzip > compare_versions/gwax.loci.gz
rm compare_versions/v2.meta.locus_*.gz
rm compare_versions/v3.meta.locus_*.gz
rm compare_versions/Lambert.locus_*.gz
rm compare_versions/Kunkle.locus_*.gz
rm compare_versions/gwax.locus_*.gz



# Merge together finemap output files for each version

# Previously this is where I got v2 finemap stats from, but this isn't what we actually
# used in my spreadsheet before so I've changed to the below.
# head -n 1 $JS/AD_finemap_v1/finemap/output/1_161156033.IGAP1_GWAX.1.snp > compare_versions/v2.meta.finemap.snp
# while read line; do
#   echo $line
#   lead_chrpos=`echo $line | awk '{print $16}'`
#   fname=$JS/AD_finemap_v1/finemap/output/$lead_chrpos.IGAP1_GWAX.1.snp
#   if [ ! -f $fname ]; then
#     fname=$JS/AD_finemap_v1/finemap/output/$lead_chrpos.IGAP1_GWAX.2.snp
#   fi
#   echo $fname
#   sed '1d' $fname >> compare_versions/v2.meta.finemap.snp
# done < <(sed '1d' AD.compare_loci.tsv)

# This also didn't match what was in my previous spreadsheet, argh!
# head -n 1 $JS/AD_finemap_v1/finemap_credset/1_161156033.IGAP1_GWAX_v2.1.snp > compare_versions/v2.meta.finemap.snp
# while read line; do
#   #echo $line
#   lead_chrpos=`echo $line | awk '{print $16}'`
#   fname=$JS/AD_finemap_v1/finemap_credset/$lead_chrpos.IGAP1_GWAX_v2.1.snp
#   echo $fname
#   if [ ! -f $fname ]; then
#     fname=$JS/AD_finemap_v1/finemap_credset/$lead_chrpos.IGAP1_GWAX_v2.2.snp
#   fi
#   if [ ! -f $fname ]; then
#     echo 'FILE NOTE FOUND!'
#   fi
#   # Fix chromosome names which are e.g. "01" to "1"
#   sed '1d' $fname | perl -ane 'if ($F[2]=~/0([\d])/) {$F[2]=$1} print join(" ", @F)."\n"' >> compare_versions/v2.meta.finemap.snp
# done < <(sed '1d' AD.compare_loci.tsv)

cat $JS/AD_finemap_v1/annotated/AD.credible_sets.annotated.colocs.txt | cut -f 3-5,7,8,21 > compare_versions/v2.meta.finemap.snp

head -n 1 $JS/AD_finemap/finemap/out_ncausal_1/1_161155392.IGAP1_GWAX.1.snp > compare_versions/v3.meta.finemap.snp
while read line; do
  echo $line
  lead_chrpos=`echo $line | awk '{print $17}'`
  ncausal=`echo $line | awk '{print $8}'`
  fname=$JS/AD_finemap/finemap/out_ncausal_1/$lead_chrpos.IGAP1_GWAX.1.snp
  if [ "$ncausal" -eq 2 ]; then
    fname=$JS/AD_finemap/finemap/out_ncausal_2/$lead_chrpos.IGAP1_GWAX.2.snp
  fi
  sed '1d' $fname >> compare_versions/v3.meta.finemap.snp
done < <(sed '1d' AD.compare_loci.tsv)

sed -i 's/ /\t/g' compare_versions/v2.meta.finemap.snp
sed -i 's/ /\t/g' compare_versions/v3.meta.finemap.snp

Rscript $JS/src/gwas/compare_gwas_versions.R AD.compare_loci.tsv \
                                compare_versions/v2.meta.loci.gz \
                                compare_versions/v3.meta.loci.gz \
                                compare_versions/Kunkle.loci.gz \
                                compare_versions/Lambert.loci.gz \
                                compare_versions/gwax.loci.gz \
                                compare_versions/v2.meta.finemap.snp \
                                compare_versions/v3.meta.finemap.snp \
                                compare_versions/result



# Compare FINEMAP results on the V3 meta with differing numbers of causal SNPs
# First prepare input files
arr=( 2_106366056 2_127892810 6_40942196 7_143107588 8_27468503 11_60021948 15_59022615 16_81773209 17_61560763 19_1039444 21_27534261 )
head -n 1 finemap/out_ncausal_1/1_161155392.IGAP1_GWAX.1.snp > compare_versions/v3.meta.finemap.ncausal_1.snp
for locus in "${arr[@]}"
do
    echo "$locus" # or do whatever with individual element of the array
    sed '1d' finemap/out_ncausal_1/$locus.IGAP1_GWAX.1.snp >> compare_versions/v3.meta.finemap.ncausal_1.snp
done

head -n 1 finemap/out_ncausal_2/1_161155392.IGAP1_GWAX.2.snp > compare_versions/v3.meta.finemap.ncausal_2.snp
for locus in "${arr[@]}"
do
    echo "$locus" # or do whatever with individual element of the array
    sed '1d' finemap/out_ncausal_2/$locus.IGAP1_GWAX.2.snp >> compare_versions/v3.meta.finemap.ncausal_2.snp
done

arr=( 6_40942196 )
head -n 1 finemap/out_ncausal_3_minfreq_0.001/6_40942196.IGAP1_GWAX.3.snp > compare_versions/v3.meta.finemap.ncausal_3.snp
for locus in "${arr[@]}"
do
    echo "$locus" # or do whatever with individual element of the array
    sed '1d' finemap/out_ncausal_3_minfreq_0.001/$locus.IGAP1_GWAX.3.snp >> compare_versions/v3.meta.finemap.ncausal_3.snp
done

arr=( 6_40942196 )
head -n 1 finemap/out_ncausal_4_minfreq_0.001/6_40942196.IGAP1_GWAX.4.snp > compare_versions/v3.meta.finemap.ncausal_4.snp
for locus in "${arr[@]}"
do
    echo "$locus" # or do whatever with individual element of the array
    sed '1d' finemap/out_ncausal_4_minfreq_0.001/$locus.IGAP1_GWAX.4.snp >> compare_versions/v3.meta.finemap.ncausal_4.snp
done


sed -i 's/ /\t/g' compare_versions/v3.meta.finemap.ncausal_1.snp
sed -i 's/ /\t/g' compare_versions/v3.meta.finemap.ncausal_2.snp
sed -i 's/ /\t/g' compare_versions/v3.meta.finemap.ncausal_3.snp
sed -i 's/ /\t/g' compare_versions/v3.meta.finemap.ncausal_4.snp


zcat $JS/AD_finemap/annotated/AD.IGAP1_GWAX_v3.annotated.selected.tsv.gz | cut -f 11,30 > compare_versions/v3.meta.finemap_nc.tsv

Rscript $JS/src/gwas/compare_finemap_ncausal.R AD.compare_loci.finemap.tsv \
                                compare_versions/v3.meta.finemap.ncausal_1.snp \
                                compare_versions/v3.meta.finemap.ncausal_2.snp \
                                compare_versions/v3.meta.finemap.ncausal_3.snp \
                                compare_versions/v3.meta.finemap.ncausal_4.snp \
                                gcta/output_1e-5/cond/merged_loci.cond.out.flt.tsv \
                                compare_versions/finemap.ncausal_comparison



head -n 1 compare_versions/result.meta_v3_merged.tsv > compare_versions/results.genie_snps.tsv
sed '1d' compare_versions/result.meta_v3_merged.tsv | grep -wF -f <(sed '1d' $JS/experiment/transcribed/all_genie_snps.txt | cut -f 1) \
    >> compare_versions/results.genie_snps.tsv

