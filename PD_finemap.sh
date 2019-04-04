###############################################################################
# This script sets up some links needed by the AD finemapping annotation script
# located in jeremys/AD_finemap.
OTCOREGEN=/lustre/scratch115/realdata/mdt3/projects/otcoregen
JS=$OTCOREGEN/jeremys
PD_FINEMAP=$JS/PD_finemap 
cd $PD_FINEMAP


# Backup all data provided by Biogen
rsync -a -u $OTCOREGEN/AD_PD_finemap /warehouse/compgen_wh05/js29/opentargets

mkdir $PD_FINEMAP/finemap_credset
cp $OTCOREGEN/AD_PD_finemap/PD_credible_sets/* $PD_FINEMAP/finemap_credset
cp $OTCOREGEN/AD_PD_finemap/PD_finemap_output/* $PD_FINEMAP/finemap_credset

# Create symlinks in AD_finemap dir to UKBB summary stats
mkdir $PD_FINEMAP/summary_stats
ln -s $OTCOREGEN/AD_PD_finemap/summary_stats/PD.proxy.bgen.stats.gz $PD_FINEMAP/summary_stats/PD.proxy.bgen.stats.gz
ln -s $OTCOREGEN/AD_PD_finemap/summary_stats/pd-gwas+gwax.summarystats_2e6-signalWindows_20180330.txt.gz $PD_FINEMAP/summary_stats/pd-gwas+gwax.summarystats_2e6-signalWindows_20180330.txt.gz
#tabix -s 1 -b 2 -e 2 -S 1 $AD_FINEMAP/summary_stats/AD.IGAP1_GWAX_v3.meta.bgz

mkdir $PD_FINEMAP/reference
ln -s $OTCOREGEN/jeremys/reference/dbSNP $PD_FINEMAP/reference/dbSNP
ln -s $OTCOREGEN/jeremys/reference/vep_impact_table.tsv $PD_FINEMAP/reference/vep_impact_table.tsv

ln -s $OTCOREGEN/jeremys/datasets/SpliceAI/spliceai.merge.tsv.bgz $PD_FINEMAP/reference/spliceai.merge.tsv.bgz
ln -s $OTCOREGEN/jeremys/datasets/SpliceAI/spliceai.merge.tsv.bgz.tbi $PD_FINEMAP/reference/spliceai.merge.tsv.bgz.tbi


mkdir $PD_FINEMAP/reference/annotations
cd $PD_FINEMAP/reference/annotations
ln -s /warehouse/compgen_wh05/js29/annotations/Roadmap Roadmap
zcat $OTCOREGEN/jeremys/ipsneurons/GRCh37/ATAC/peaks/atac_ipsc_peaks.narrowPeak.gz | tr -s '\n' | gzip > $PD_FINEMAP/reference/annotations/atac_ipsc_peaks.GRCh37.narrowPeak.gz
zcat $OTCOREGEN/jeremys/ipsneurons/GRCh37/ATAC/peaks/atac_npc_peaks.narrowPeak.gz | tr -s '\n' | gzip > $PD_FINEMAP/reference/annotations/atac_npc_peaks.GRCh37.narrowPeak.gz
zcat $OTCOREGEN/jeremys/ipsneurons/GRCh37/ATAC/peaks/atac_neuron_peaks.narrowPeak.gz | tr -s '\n' | gzip > $PD_FINEMAP/reference/annotations/atac_neuron_peaks.GRCh37.narrowPeak.gz
zcat $OTCOREGEN/jeremys/ipsneurons/GRCh37/ATAC/peaks/atac_ineuron_peaks.narrowPeak.gz | tr -s '\n' | gzip > $PD_FINEMAP/reference/annotations/atac_ineuron_peaks.GRCh37.narrowPeak.gz
cat $OTCOREGEN/jeremys/macrophage/GRCh37/ATAC_macrophage.peaks.GRCh37.bed | gzip > $PD_FINEMAP/reference/annotations/atac_macrophage.peaks.GRCh37.bed.gz
cat $OTCOREGEN/jeremys/datasets/glass_microglia/ATAC/peaks/atac_microglia_peaks.narrowPeak | gzip > $PD_FINEMAP/reference/annotations/atac_microglia_peaks.GRCh37.narrowPeak.gz
cat $OTCOREGEN/jeremys/sensoryneurons/GRCh37/ATAC_consensus_peaks.bed | gzip > $PD_FINEMAP/reference/annotations/atac_sensoryneuron_peaks.GRCh37.bed.gz

cp $JS/AD_finemap_v1/reference/hgnc.ensembl.map.txt $JS/AD_finemap/reference/hgnc.ensembl.map.txt
cp $JS/AD_finemap/reference/hgnc.ensembl.map.txt $JS/PD_finemap/reference/hgnc.ensembl.map.txt

#wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/human_permissive_enhancers_phase_1_and_2.bed.gz
#mv human_permissive_enhancers_phase_1_and_2.bed.gz fantom.human_permissive_enhancers_phase_1_and_2.bed.gz
ln -s $JS/AD_finemap/reference/annotations/fantom.human_permissive_enhancers_phase_1_and_2.bed.gz fantom.human_permissive_enhancers_phase_1_and_2.bed.gz

ln -s $OTCOREGEN/jeremys/annotation/JEME $PD_FINEMAP/reference/JEME

####### GWAS data for coloc
mkdir $PD_FINEMAP/coloc
ln -s $JS/datasets/GWAS $PD_FINEMAP/coloc/GWAS


####### QTL datasets
SNDIR=/warehouse/compgen_wh05/js29/sensoryneurons
md $PD_FINEMAP/coloc
ln -s $JS/AD_finemap/coloc/qtl_data $PD_FINEMAP/coloc/qtl_data
#md $PD_FINEMAP/coloc/qtl_data/sensoryneurons
#cd $PD_FINEMAP/coloc/qtl_data/sensoryneurons
# cp $SNDIR/rasqual/imputed.97_samples.snps_indels.INFO_08.selected_samples.vcf.gz $PD_FINEMAP/coloc/qtl_data/sensoryneurons/sensoryneurons.97_samples.INFO_08.GRCh38.vcf.gz
# cp $SNDIR/fastqtl/output/permutations.10k.allgenes.cis500k.PCs.20.fdr0.1.coords.txt $PD_FINEMAP/coloc/qtl_data/sensoryneurons/
# cp $SNDIR/fastqtl/output/nominals.allgenes.cis500k.ePCs.20.pvalues.notskipped.coords.txt.gz $PD_FINEMAP/coloc/qtl_data/sensoryneurons/

ln -s $JS/reference/GRCh38/Homo_sapiens.GRCh38.91.exon_start_end.simple.bed $PD_FINEMAP/reference/Homo_sapiens.GRCh38.91.exon_start_end.simple.bed
ln -s $JS/reference/GRCh38/gencode.v27.geneid_to_hgnc.txt $PD_FINEMAP/reference/gencode.v27.geneid_to_hgnc.txt
ln -s $JS/reference/GRCh37/GRCh37.75.dna.toplevel.fa.gz $PD_FINEMAP/reference/GRCh37.75.dna.toplevel.fa.gz


