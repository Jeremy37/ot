
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
KAUR=/lustre/scratch117/cellgen/team170/ka8/projects/macrophage-gxe-study/results
KAUR_ATAC=$KAUR/ATAC/fastqtl/output
KAUR_EQTL=$KAUR/SL1344/fastqtl/output
KAUR_SQTL=$KAUR/SL1344/leafcutter/fastqtl_output_col_quantile

# PTK2B
SNP=rs2322599
SNP=rs28834970
SNP=rs6987305
zgrep -P "$SNP " --with-filename \
  $KAUR_SQTL/naive_100kb_pvalues.coords.txt.gz $KAUR_SQTL/IFNg_100kb_pvalues.coords.txt.gz \
  $KAUR_SQTL/SL1344_100kb_pvalues.coords.txt.gz $KAUR_SQTL/IFNg_SL1344_100kb_pvalues.coords.txt.gz \
  > $JS/misc/$SNP.$ASSOC.mac_sqtl.fastql.txt

zgrep -P "$SNP\t" --with-filename \
  $KAUR_EQTL/naive_500kb_pvalues.sorted.txt.gz $KAUR_EQTL/IFNg_500kb_pvalues.sorted.txt.gz \
  $KAUR_EQTL/SL1344_500kb_pvalues.sorted.txt.gz $KAUR_EQTL/IFNg_SL1344_500kb_pvalues.sorted.txt.gz \
  > $JS/misc/$SNP.$ASSOC.mac_eqtl.fastql.tsv

zgrep -P "$SNP\t" --with-filename \
  $KAUR_ATAC/naive_500kb_pvalues.sorted.txt.gz $KAUR_ATAC/IFNg_500kb_pvalues.sorted.txt.gz \
  $KAUR_ATAC/SL1344_500kb_pvalues.sorted.txt.gz $KAUR_ATAC/IFNg_SL1344_500kb_pvalues.sorted.txt.gz \
  > $JS/misc/$SNP.$ASSOC.mac_caqtl.fastql.tsv


# SORL1
SNP=rs11218343
# Strangely this doesn't find any tests in Kaur's macrophages. Why wasn't this variant tested?

# FTO lead obesity SNP
SNP=rs9939609

# SPPL2A lead AD SNP
SNP=rs17646025
ASSOC=SPPL2A

# Get splicing associations for the alternatively spliced SPPL2A exon in macrophages
zgrep "clu_29460" --with-filename \
  $KAUR_SQTL/naive_100kb_pvalues.coords.txt.gz $KAUR_SQTL/IFNg_100kb_pvalues.coords.txt.gz \
  $KAUR_SQTL/SL1344_100kb_pvalues.coords.txt.gz $KAUR_SQTL/IFNg_SL1344_100kb_pvalues.coords.txt.gz \
  > $JS/misc/$ASSOC.clu_29460.mac_sqtl.fastql.txt

# Get AD associations for all UKBB variants near the exon, in case relevant variants were
# missed in the meta-analysis
zcat $JS/../AD_PD_finemap/summary_stats/ADD.proxy_v2.bgen.stats.gz | awk '$2 == "15" && $3 > 50919812 && $3 < 51119812' > $JS/misc/rs17646025.within_100kb.ukbb_proxy_v2.txt

# Get AD associations for SPPL2A region in Kunkle et al.
zcat $JS/datasets/GWAS/Kunkle_etal_Stage1_results.txt.gz | awk '$1 == "15" && $2 > 50919812 && $2 < 51119812' > $JS/misc/rs17646025.within_100kb.kunkle.txt

