###############################################################################
# This script sets up some links needed by the AD finemapping annotation script
# located in jeremys/AD_finemap.
OTCOREGEN=/lustre/scratch115/realdata/mdt3/projects/otcoregen
JS=$OTCOREGEN/jeremys
AD_FINEMAP=$JS/AD_finemap 
cd $AD_FINEMAP


# Backup all data provided by Biogen
rsync -a -u $OTCOREGEN/AD_PD_finemap /warehouse/compgen_wh05/js29/opentargets

mkdir $AD_FINEMAP/finemap_credset
cp $OTCOREGEN/AD_PD_finemap/AD_credible_sets/v2/* $AD_FINEMAP/finemap_credset
cp $OTCOREGEN/AD_PD_finemap/AD_finemap_output/v2/* $AD_FINEMAP/finemap_credset

# Create symlinks in AD_finemap dir to UKBB summary stats
mkdir $AD_FINEMAP/summary_stats
ln -s $OTCOREGEN/AD_PD_finemap/summary_stats/ADD.proxy_v2.bgen.stats.gz $AD_FINEMAP/summary_stats/ADD.proxy_v2.bgen.stats.gz
zcat $OTCOREGEN/AD_PD_finemap/summary_stats/AD.IGAP1_GWAX_v3.meta.gz | tr ' ' '\t' | bgzip > $AD_FINEMAP/summary_stats/AD.IGAP1_GWAX_v3.meta.bgz
tabix -s 1 -b 2 -e 2 -S 1 $AD_FINEMAP/summary_stats/AD.IGAP1_GWAX_v3.meta.bgz

ln -s $OTCOREGEN/AD_PD_finemap/ukbb_frequencies/ $AD_FINEMAP/summary_stats/ukbb_frequencies
ln -s $OTCOREGEN/AD_PD_finemap/ukbb_sample_regions/ $AD_FINEMAP/gcta/input/ukbb_sample_regions

ln -s $OTCOREGEN/jeremys/datasets/SpliceAI/spliceai.merge.tsv.bgz $AD_FINEMAP/reference/spliceai.merge.tsv.bgz
ln -s $OTCOREGEN/jeremys/datasets/SpliceAI/spliceai.merge.tsv.bgz.tbi $AD_FINEMAP/reference/spliceai.merge.tsv.bgz.tbi


mkdir $AD_FINEMAP/reference
ln -s $OTCOREGEN/jeremys/reference/dbSNP $AD_FINEMAP/reference/dbSNP
ln -s $OTCOREGEN/jeremys/reference/vep_impact_table.tsv $AD_FINEMAP/reference/vep_impact_table.tsv

mkdir $AD_FINEMAP/reference/annotations
cd $AD_FINEMAP/reference/annotations
ln -s /warehouse/compgen_wh05/js29/annotations/Roadmap Roadmap
cat $OTCOREGEN/jeremys/ipsneurons/GRCh37/ATAC/peaks/atac_ipsc_peaks.narrowPeak | tr -s '\n' > $AD_FINEMAP/reference/annotations/atac_ipsc_peaks.GRCh37.narrowPeak
cat $OTCOREGEN/jeremys/ipsneurons/GRCh37/ATAC/peaks/atac_npc_peaks.narrowPeak | tr -s '\n' > $AD_FINEMAP/reference/annotations/atac_npc_peaks.GRCh37.narrowPeak
cat $OTCOREGEN/jeremys/ipsneurons/GRCh37/ATAC/peaks/atac_neuron_peaks.narrowPeak | tr -s '\n' > $AD_FINEMAP/reference/annotations/atac_neuron_peaks.GRCh37.narrowPeak
cat $OTCOREGEN/jeremys/ipsneurons/GRCh37/ATAC/peaks/atac_ineuron_peaks.narrowPeak | tr -s '\n' > $AD_FINEMAP/reference/annotations/atac_ineuron_peaks.GRCh37.narrowPeak
cp $OTCOREGEN/jeremys/macrophage/GRCh37/ATAC_macrophage.peaks.GRCh37.bed $AD_FINEMAP/reference/annotations/atac_macrophage.peaks.GRCh37.bed
cp $OTCOREGEN/jeremys/datasets/glass_microglia/ATAC/peaks/atac_microglia_peaks.narrowPeak $AD_FINEMAP/reference/annotations/atac_microglia_peaks.GRCh37.narrowPeak
cp $OTCOREGEN/jeremys/sensoryneurons/GRCh37/ATAC_consensus_peaks.bed $AD_FINEMAP/reference/annotations/atac_sensoryneuron_peaks.GRCh37.bed
gzip $AD_FINEMAP/reference/annotations/*.narrowPeak
gzip $AD_FINEMAP/reference/annotations/*.bed

wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/human_permissive_enhancers_phase_1_and_2.bed.gz
mv human_permissive_enhancers_phase_1_and_2.bed.gz fantom.human_permissive_enhancers_phase_1_and_2.bed.gz

md $AD_FINEMAP/gcta/input
ln -s $OTCOREGEN/AD_PD_finemap/ukbb_sample_regions $AD_FINEMAP/gcta/input/ukbb_sample_regions

ln -s $OTCOREGEN/jeremys/annotation/JEME $AD_FINEMAP/reference/JEME

####### GWAS data for coloc
mkdir $AD_FINEMAP/coloc
ln -s $JS/datasets/GWAS $AD_FINEMAP/coloc/GWAS


####### QTL datasets
SNDIR=/warehouse/compgen_wh05/js29/sensoryneurons
md $AD_FINEMAP/coloc/qtl_data/sensoryneurons
cd $AD_FINEMAP/coloc/qtl_data/sensoryneurons
# cp $SNDIR/rasqual/imputed.97_samples.snps_indels.INFO_08.selected_samples.vcf.gz $AD_FINEMAP/coloc/qtl_data/sensoryneurons/sensoryneurons.97_samples.INFO_08.GRCh38.vcf.gz
# cp $SNDIR/fastqtl/output/permutations.10k.allgenes.cis500k.PCs.20.fdr0.1.coords.txt $AD_FINEMAP/coloc/qtl_data/sensoryneurons/
# cp $SNDIR/fastqtl/output/nominals.allgenes.cis500k.ePCs.20.pvalues.notskipped.coords.txt.gz $AD_FINEMAP/coloc/qtl_data/sensoryneurons/

ln -s $OTCOREGEN/jeremys/reference/GRCh37 $AD_FINEMAP/reference/GRCh37
ln -s $SNDIR/datasubmission/biostudies/updated/eqtl.fastqtl.500k.nominals.txt.gz eqtl.fastqtl.500k.nominals.txt.gz
ln -s $SNDIR/datasubmission/biostudies/updated/eqtl.fastqtl.500k.permutations.txt.gz eqtl.fastqtl.500k.permutations.txt.gz
ln -s $SNDIR/datasubmission/biostudies/updated/sqtl.fastqtl.500k.nominals.txt.gz sqtl.fastqtl.500k.nominals.txt.gz
ln -s $SNDIR/datasubmission/biostudies/updated/sqtl.fastqtl.500k.permutations.txt.gz sqtl.fastqtl.500k.permutations.txt.gz
mkdir $AD_FINEMAP/reference/sensoryneurons
ln -s $JS/sensoryneurons/GRCh37/imputed.97_samples.snps_indels.INFO_08.GRCh37.vcf.gz $AD_FINEMAP/reference/sensoryneurons/sensoryneurons.GRCh37.vcf.gz


## Blueprint
BPDIR=$AD_FINEMAP/coloc/qtl_data/blueprint
ln -s $JS/datasets/blueprint/README.QTL $BPDIR/README.QTL
ln -s $JS/datasets/blueprint/mono_gene_nor_combat_peer_10_all_summary.txt.gz $BPDIR/mono_gene_nor_combat_peer_10_all_summary.txt.gz
ln -s $JS/datasets/blueprint/mono_psi_peer_10_all_summary.txt.gz $BPDIR/mono_psi_peer_10_all_summary.txt.gz
ln -s $JS/datasets/blueprint/mono_K27AC_log2rpm_peer_10_all_summary.txt.gz $BPDIR/mono_K27AC_log2rpm_peer_10_all_summary.txt.gz
ln -s $JS/datasets/blueprint/mono_K4ME1_log2rpm_peer_10_all_summary.txt.gz $BPDIR/mono_K4ME1_log2rpm_peer_10_all_summary.txt.gz
ln -s $JS/datasets/blueprint/mono_gene_nor_combat_peer_10_all_summary.txt.gz $BPDIR/mono_gene_nor_combat_peer_10_all_summary.txt.gz
ln -s $JS/datasets/blueprint/mono_psi_peer_10_all_summary.txt.gz $BPDIR/mono_psi_peer_10_all_summary.txt.gz
ln -s $JS/datasets/blueprint/neut_K27AC_log2rpm_peer_10_all_summary.txt.gz $BPDIR/neut_K27AC_log2rpm_peer_10_all_summary.txt.gz
ln -s $JS/datasets/blueprint/neut_K4ME1_log2rpm_peer_10_all_summary.txt.gz $BPDIR/neut_K4ME1_log2rpm_peer_10_all_summary.txt.gz
ln -s $JS/datasets/blueprint/neut_gene_nor_combat_peer_10_all_summary.txt.gz $BPDIR/neut_gene_nor_combat_peer_10_all_summary.txt.gz
ln -s $JS/datasets/blueprint/neut_psi_peer_10_all_summary.txt.gz $BPDIR/neut_psi_peer_10_all_summary.txt.gz
ln -s $JS/datasets/blueprint/tcel_K27AC_log2rpm_peer_10_all_summary.txt.gz $BPDIR/tcel_K27AC_log2rpm_peer_10_all_summary.txt.gz
ln -s $JS/datasets/blueprint/tcel_K4ME1_log2rpm_peer_10_all_summary.txt.gz $BPDIR/tcel_K4ME1_log2rpm_peer_10_all_summary.txt.gz
ln -s $JS/datasets/blueprint/tcel_gene_nor_combat_peer_10_all_summary.txt.gz $BPDIR/tcel_gene_nor_combat_peer_10_all_summary.txt.gz
ln -s $JS/datasets/blueprint/tcel_psi_peer_10_all_summary.txt.gz $BPDIR/tcel_psi_peer_10_all_summary.txt.gz

mv /lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/coloc/qtl_data/blueprint/blueprint* $BPDIR/

## xQTL
mv /lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/coloc/qtl_data/xQTL/xQTL* $AD_FINEMAP/coloc/qtl_data/xQTL/
ln -s $JS/reference/xQTL/eQTLs_all.txt.gz $AD_FINEMAP/coloc/qtl_data/xQTL/eQTLs_all.txt.gz
ln -s $JS/reference/xQTL/mQTLs_all.txt.gz $AD_FINEMAP/coloc/qtl_data/xQTL/mQTLs_all.txt.gz
ln -s $JS/reference/xQTL/haQTLs_all.txt.gz $AD_FINEMAP/coloc/qtl_data/xQTL/haQTLs_all.txt.gz
ln -s $JS/reference/xQTL/maf.gz $AD_FINEMAP/coloc/qtl_data/xQTL/maf.gz


ln -s $JS/reference/GRCh38/Homo_sapiens.GRCh38.91.exon_start_end.simple.bed $JS/AD_finemap/reference/Homo_sapiens.GRCh38.91.exon_start_end.simple.bed
ln -s $JS/reference/GRCh38/gencode.v27.geneid_to_hgnc.txt $JS/AD_finemap/reference/gencode.v27.geneid_to_hgnc.txt
ln -s $JS/reference/GRCh37/GRCh37.75.dna.toplevel.fa.gz $JS/AD_finemap/reference/GRCh37.75.dna.toplevel.fa.gz


## GTEx
ln -s /lustre/scratch118/humgen/resources/GTEx/AnalysisV7/GTEx_Analysis_v7_eQTL_all_associations $AD_FINEMAP/coloc/qtl_data/GTEx/summary

