JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
ROADMAP=/lustre/scratch117/cellgen/team170/js29/annotations/Roadmap
cd $JS/gwas/AD

# First prepare the input fine mapping file from Toby
# R code
df = read.csv("/lustre/scratch117/cellgen/team170/js29/opentargets/AD/data/TJ.AD.all_causal_v1.csv")
write.table(df, "/lustre/scratch117/cellgen/team170/js29/opentargets/AD/data/TJ.AD.all_causal_v1.orig.txt", quote=F, col.names=T, row.names=F, sep="\t")
##
head -n 1 $JS/gwas/AD/data/TJ.AD.all_causal_v1.orig.txt | cut -f 1 --complement > $JS/gwas/AD/data/TJ.AD.all_causal_v1.txt
sed '1d' $JS/gwas/AD/data/TJ.AD.all_causal_v1.orig.txt | awk 'BEGIN{OFS="\t"}{print $12"_"$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' \
   >> $JS/gwas/AD/data/TJ.AD.all_causal_v1.txt

$JS/gwas/src/doOverlaps.sh $JS/gwas/AD/data/TJ.AD.all_causal_v1.txt
