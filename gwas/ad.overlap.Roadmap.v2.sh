JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
cd $JS/gwas/AD

AD_FINEMAP=$JS/../AD_PD_finemap/AD_finemap_output
head -n 1 $AD_FINEMAP/1_207750568.IGAP1_GWAX.1.snp | tr ' ' '\t' > data/IGAP1_GWAX.combined.snp
for F in $AD_FINEMAP/*.IGAP1_GWAX.*.snp; do
  sed '1d' $F | tr ' ' '\t' >> data/IGAP1_GWAX.combined.snp
done

AD_CREDSET=$JS/../AD_PD_finemap/AD_credible_sets
head -n 1 $AD_CREDSET/1_207750568.IGAP1_GWAX.1.set | tr ' ' '\t' > data/IGAP1_GWAX.combined.set
for F in $AD_CREDSET/*.IGAP1_GWAX.*.set; do
  sed '1d' $F | tr ' ' '\t' >> data/IGAP1_GWAX.combined.set
done

../src/merge_Toby_Jimmy_finemap.R data/TJ.AD.all_causal_v1.txt data/IGAP1_GWAX.combined.set data/IGAP1_GWAX.combined.snp

$JS/gwas/src/doOverlaps.sh toby.jimmy.finemap.merged.txt

