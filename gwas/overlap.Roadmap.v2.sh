JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
cd $JS/gwas/AD

AD_FINEMAP=$JS/gwas/AD/AD_finemap_output
head -n 1 $AD_FINEMAP/1_207750568.IGAP1_GWAX.1.snp | tr ' ' '\t' > data/IGAP1_GWAX.combined.snp
for F in $AD_FINEMAP/*.IGAP1_GWAX.*.snp; do
  sed '1d' $F | tr ' ' '\t' >> data/IGAP1_GWAX.combined.snp
done

AD_CREDSET=$JS/gwas/AD/AD_credible_sets
head -n 1 $AD_CREDSET/1_207750568.IGAP1_GWAX.1.set | tr ' ' '\t' > data/IGAP1_GWAX.combined.set
for F in $AD_CREDSET/*.IGAP1_GWAX.*.set; do
  sed '1d' $F | tr ' ' '\t' >> data/IGAP1_GWAX.combined.set
done

$JS/src/gwas/merge_Toby_Jimmy_finemap.R data/TJ.AD.all_causal_v1.txt data/IGAP1_GWAX.combined.set data/IGAP1_GWAX.combined.snp

$JS/src/gwas/doOverlaps.sh toby.jimmy.finemap.merged.txt toby.jimmy.finemap.merged

bash $JS/src/gwas/overlap.JEME.sh toby.jimmy.finemap.merged.input.chr.bed
Rscript $JS/src/gwas/merge_JEME_overlaps.R
# produces file toby.jimmy.finemap.annotated.roadmapEnhPromLinks.txt



################################################################################
# PD

OT=/lustre/scratch115/realdata/mdt3/projects/otcoregen
JS=$OT/jeremys
cd $JS/gwas/PD
NAME=PD.finemap

cp -r $OT/AD_PD_finemap/PD_credible_sets .
cp -r $OT/AD_PD_finemap/PD_finemap_output .
mkdir finemap
mkdir ipsneurons
mkdir roadmap
cp archive/ipsneurons/*bedfile.txt ipsneurons
cp archive/roadmap/*bedfilelist*.txt roadmap

PD_CREDSET=$JS/gwas/PD/PD_credible_sets
paste <(echo "locus") <(head -n 1 $PD_CREDSET/pd-gwas+gwax.caucasian.1_226916078.1.set | tr ' ' '\t') > finemap/pd_gwas_gwax.combined.set
for F in $PD_CREDSET/pd-gwas+gwax.caucasian.*.set; do
  locus=`echo $F | perl -ne '@fparts=split(/\./); print join(".", $fparts[2]);'`
  sed '1d' $F | tr ' ' '\t' | perl -sane 'chomp; print join("\t", $locus, @F)."\n";' -- -locus=$locus >> finemap/pd_gwas_gwax.combined.set
done

cat finemap/pd_gwas_gwax.combined.set | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$3,$5,$6,$7,$8,$9,$11,$12}' > finemap/PD.credible_sets.annot_input.set

$JS/src/gwas/doOverlaps.sh finemap/PD.credible_sets.annot_input.set $NAME

bash $JS/src/gwas/overlap.JEME.sh PD.finemap.input.chr.bed

submitJobs.py --MEM 10000 -j PD.annotate.JEME -q yesterday \
    -c "Rscript $JS/src/gwas/merge_JEME_overlaps.R 'gwas/PD' 'PD.finemap.annotated'"


################################################################################
# Next run lines in run_coloc.sh
# Then run annotate.coloc.R
# Finally, filter SNPs on finemap probability (if desired)
(head -n 1 $JS/gwas/PD/PD.finemap.annotated.colocs.txt; \
 sed '1d' $JS/gwas/PD/PD.finemap.annotated.colocs.txt | awk '$8 > 0') \
 > $JS/gwas/PD/PD.finemap.annotated.colocs.snpprob_gt_0.txt

(head -n 1 $JS/gwas/PD/PD.finemap.annotated.colocs.txt; \
 sed '1d' $JS/gwas/PD/PD.finemap.annotated.colocs.txt | awk '$7 < 1e-5 || $8 > 0') \
 > $JS/gwas/PD/PD.finemap.annotated.colocs.p_lt_1e-5.txt

sed '1d' $JS/gwas/PD/PD.finemap.annotated.colocs.txt | awk '$7 >= 1e-5 && $8 > 0' | wc -l



gzip $JS/gwas/PD/PD.finemap.annotated.colocs.txt
