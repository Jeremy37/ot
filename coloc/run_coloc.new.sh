JS=/Users/jeremys/work/opentargets
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

cd $JS/coloc

GWAS_NAME=AD.meta
GWAS_FILE=$JS/datasets/GWAS/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.txt.gz
GWAS_SIGNALS=$JS/coloc/AD.signals.txt

GWAS_NAME=PD.meta
GWAS_FILE=$JS/datasets/GWAS/Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.sorted.txt.gz
GWAS_SIGNALS=$JS/coloc/PD.signals.txt

PLOT_THRESHOLD=0.5

# Make GWAS signals files
echo -e "chr\tpos\trsid\tlocus\tlocus_name\tp\tlead_snp_prob" > $JS/coloc/AD.signals.txt
sed '1d' $JS/gwas/AD/AD.finemap.leadSNPs.txt | awk 'BEGIN{OFS="\t"}{print $5,$6,$3,$1,$2,$4,$7}' >> $JS/coloc/AD.signals.txt

echo -e "chr\tpos\trsid\tlocus\tp\tlead_snp_prob" > $JS/coloc/PD.signals.txt
sed '1d' $JS/gwas/PD/PD.finemap.leadSNPs.txt | awk 'BEGIN{OFS="\t"}{print $4,$5,$2,$1,$3,$6}' >> $JS/coloc/PD.signals.txt


###############################################################################
# iPSC-derived sensory neurons
QTL_DIR=$JS/coloc/qtl_data/sensory_neuron

submitJobs.py --MEM 3000 -j coloc.$GWAS_NAME.sens_neur.500k -q normal \
    -c "Rscript $JS/src/coloc/qtlColoc.R \
        --qtl_name sens_neur.eqtl.500k \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/sensory_neuron.fastqtl.cis500k.nominals.GRCh38.txt.gz \
        --qtl_signals $QTL_DIR/sensory_neuron.fastqtl.cis500k.qtl_signals.txt  \
        --qtl_variant_info $QTL_DIR/sensory_neuron.variant_info.GRCh37.gz  \
        --qtl_samplesize 97 \
        --outdir $JS/coloc/new \
        --plot_style line --plot_threshold $PLOT_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5"
#--qtl_list


# NEED TO FIX the sQTLs - coloc doesn't work when such a small window is used
# for the sQTLs, because generally the sQTL peak doesn't overlap the GWAS peak.
submitJobs.py --MEM 3000 -j coloc.$GWAS_NAME.sens_neur.sqtl.500k -q normal \
    -c "Rscript $JS/src/coloc/qtlColoc.R \
        --qtl_name sens_neur.sqtl \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/sensory_neuron.sqtl.fastqtl.500k.nominals.GRCh38.txt.gz \
        --qtl_signals $QTL_DIR/sensory_neuron.sqtl.fastqtl.qtl_signals.txt  \
        --qtl_variant_info $QTL_DIR/sensory_neuron.variant_info.GRCh37.gz  \
        --qtl_samplesize 97 \
        --outdir $JS/coloc/new \
        --plot_style line --plot_threshold $PLOT_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5"



###############################################################################
# Blueprint
QTL_DIR=$JS/coloc/qtl_data/blueprint

submitJobs.py --MEM 3000 -j coloc.$GWAS_NAME.blueprint.mono_gene_eQTL -q normal \
    -c "Rscript $JS/src/coloc/qtlColoc.R \
        --qtl_name blueprint.mono_gene_eQTL \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/blueprint.mono_gene_eQTL.nominals.txt.gz \
        --qtl_signals $QTL_DIR/blueprint.mono_gene_eQTL.qtl_signals.txt.gz  \
        --qtl_variant_info $QTL_DIR/blueprint.mono_gene_eQTL.variant_info.txt.gz  \
        --qtl_samplesize 194 \
        --outdir $JS/coloc/new \
        --plot_style line --plot_threshold $PLOT_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5"

submitJobs.py --MEM 3000 -j coloc.$GWAS_NAME.blueprint.mono_K27AC -q normal \
    -c "Rscript $JS/src/coloc/qtlColoc.R \
        --qtl_name blueprint.mono_K27AC \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/blueprint.mono_K27AC.nominals.txt.gz \
        --qtl_signals $QTL_DIR/blueprint.mono_K27AC.qtl_signals.txt.gz  \
        --qtl_variant_info $QTL_DIR/blueprint.mono_K27AC.variant_info.txt.gz  \
        --qtl_samplesize 172 \
        --outdir $JS/coloc/new \
        --plot_style line --plot_threshold $PLOT_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5"

submitJobs.py --MEM 3000 -j coloc.$GWAS_NAME.blueprint.mono_K4ME1 -q normal \
    -c "Rscript $JS/src/coloc/qtlColoc.R \
        --qtl_name blueprint.mono_K4ME1 \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/blueprint.mono_K4ME1.nominals.txt.gz \
        --qtl_signals $QTL_DIR/blueprint.mono_K4ME1.qtl_signals.txt.gz  \
        --qtl_variant_info $QTL_DIR/blueprint.mono_K4ME1.variant_info.txt.gz  \
        --qtl_samplesize 162 \
        --outdir $JS/coloc/new \
        --plot_style line --plot_threshold $PLOT_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5"


###############################################################################
# xQTL
QTL_DIR=$JS/coloc/qtl_data/xQTL

QTL=xQTL_eQTL
submitJobs.py --MEM 8000 -j coloc.$GWAS_NAME.$QTL -q normal \
    -c "Rscript $JS/src/coloc/qtlColoc.R \
        --qtl_name $QTL \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/$QTL.nominals.txt.gz \
        --qtl_signals $QTL_DIR/$QTL.qtl_signals.p_lt_1e-4.txt  \
        --qtl_variant_info $QTL_DIR/$QTL.variant_info.txt.gz  \
        --qtl_samplesize 194 \
        --outdir $JS/coloc/new \
        --plot_style line --plot_threshold $PLOT_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5"

QTL=xQTL_mQTL
submitJobs.py --MEM 8000 -j coloc.$GWAS_NAME.$QTL -q normal \
    -c "Rscript $JS/src/coloc/qtlColoc.R \
        --qtl_name $QTL \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/$QTL.nominals.txt.gz \
        --qtl_signals $QTL_DIR/$QTL.qtl_signals.p_lt_1e-4.txt  \
        --qtl_variant_info $QTL_DIR/$QTL.variant_info.txt.gz  \
        --qtl_samplesize 194 \
        --outdir $JS/coloc/new \
        --plot_style line --plot_threshold $PLOT_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5"

QTL=xQTL_haQTL
submitJobs.py --MEM 8000 -j coloc.$GWAS_NAME.$QTL -q normal \
    -c "Rscript $JS/src/coloc/qtlColoc.R \
        --qtl_name $QTL \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/$QTL.nominals.txt.gz \
        --qtl_signals $QTL_DIR/$QTL.qtl_signals.p_lt_1e-4.txt  \
        --qtl_variant_info $QTL_DIR/$QTL.variant_info.txt.gz  \
        --qtl_samplesize 194 \
        --outdir $JS/coloc/new \
        --plot_style line --plot_threshold $PLOT_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5"


###############################################################################
# iPSC-derived macrophages
QTL_DIR=$JS/coloc/qtl_data/macrophage

STATE=naive
STATE=IFNg
STATE=SL1344
STATE=IFNg_SL1344
QTL_SET="beta_lt_1e-3"
QTL_SET="beta_lt_0.05"

QTL=ips_macrophage.eQTL_$STATE
submitJobs.py --MEM 8000 -j coloc.$GWAS_NAME.$QTL -q normal \
    -c "Rscript $JS/src/coloc/qtlColoc.R \
        --qtl_name $QTL \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/ips_macrophage.eQTL_$STATE.nominals.txt.gz \
        --qtl_signals $QTL_DIR/ips_macrophage.eQTL_$STATE.signals.$QTL_SET.txt \
        --qtl_variant_info $QTL_DIR/ips_macrophage.variant_info.GRCh37.txt.gz \
        --qtl_samplesize 84 \
        --outdir $JS/coloc/new \
        --plot_style line --plot_threshold $PLOT_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5"

QTL=ips_macrophage.leafcutter_QTL_$STATE
submitJobs.py --MEM 8000 -j coloc.$GWAS_NAME.$QTL -q normal \
    -c "Rscript $JS/src/coloc/qtlColoc.R \
        --qtl_name $QTL \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/$QTL.nominals.txt.gz \
        --qtl_signals $QTL_DIR/$QTL.signals.$QTL_SET.txt \
        --qtl_variant_info $QTL_DIR/ips_macrophage.variant_info.GRCh37.txt.gz \
        --qtl_samplesize 84 \
        --outdir $JS/coloc/new \
        --plot_style line --plot_threshold $PLOT_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5"

QTL=ips_macrophage.txrevise_ends_$STATE
submitJobs.py --MEM 8000 -j coloc.$GWAS_NAME.$QTL -q normal \
    -c "Rscript $JS/src/coloc/qtlColoc.R \
        --qtl_name $QTL \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/$QTL.nominals.txt.gz \
        --qtl_signals $QTL_DIR/$QTL.signals.$QTL_SET.txt \
        --qtl_variant_info $QTL_DIR/ips_macrophage.variant_info.GRCh37.txt.gz \
        --qtl_samplesize 84 \
        --outdir $JS/coloc/new \
        --plot_style line --plot_threshold $PLOT_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5"

QTL=ips_macrophage.txrevise_promoters_$STATE
submitJobs.py --MEM 8000 -j coloc.$GWAS_NAME.$QTL -q normal \
    -c "Rscript $JS/src/coloc/qtlColoc.R \
        --qtl_name $QTL \
        --gwas_name $GWAS_NAME --gwas_file $GWAS_FILE --gwas_signals $GWAS_SIGNALS \
        --overlap_dist 5e5 --window 2e5 \
        --qtl_nominal $QTL_DIR/$QTL.nominals.txt.gz \
        --qtl_signals $QTL_DIR/$QTL.signals.$QTL_SET.txt \
        --qtl_variant_info $QTL_DIR/ips_macrophage.variant_info.GRCh37.txt.gz \
        --qtl_samplesize 84 \
        --outdir $JS/coloc/new \
        --plot_style line --plot_threshold $PLOT_THRESHOLD \
        --p1 1e-4 --p2 1e-4 --p12 1e-5"


grep "Successfully completed" FarmOut/coloc.*.180*.txt | wc -l



###############################################################################
# Annotate genes for sQTLs

Rscript $JS/src/coloc/annotate_sqtl_genes.R $JS/coloc/new/coloc.$GWAS_NAME.ips_macrophage.leafcutter_QTL_naive.5e+05.txt 1 \
  $JS/reference/GRCh38/Homo_sapiens.GRCh38.91.exon_start_end.simple.bed \
  $JS/reference/GRCh38/gencode.v27.geneid_to_hgnc.txt

Rscript $JS/src/coloc/annotate_sqtl_genes.R $JS/coloc/new/coloc.$GWAS_NAME.ips_macrophage.leafcutter_QTL_IFNg.5e+05.txt 1 \
  $JS/reference/GRCh38/Homo_sapiens.GRCh38.91.exon_start_end.simple.bed \
  $JS/reference/GRCh38/gencode.v27.geneid_to_hgnc.txt

Rscript $JS/src/coloc/annotate_sqtl_genes.R $JS/coloc/new/coloc.$GWAS_NAME.ips_macrophage.leafcutter_QTL_SL1344.5e+05.txt 1 \
  $JS/reference/GRCh38/Homo_sapiens.GRCh38.91.exon_start_end.simple.bed \
  $JS/reference/GRCh38/gencode.v27.geneid_to_hgnc.txt

Rscript $JS/src/coloc/annotate_sqtl_genes.R $JS/coloc/new/coloc.$GWAS_NAME.ips_macrophage.leafcutter_QTL_IFNg_SL1344.5e+05.txt 1 \
  $JS/reference/GRCh38/Homo_sapiens.GRCh38.91.exon_start_end.simple.bed \
  $JS/reference/GRCh38/gencode.v27.geneid_to_hgnc.txt


Rscript $JS/src/coloc/annotate_sqtl_genes.R $JS/coloc/new/coloc.$GWAS_NAME.sens_neur.sqtl.5e+05.txt 1 \
  $JS/reference/GRCh38/Homo_sapiens.GRCh38.91.exon_start_end.simple.bed \
  $JS/reference/GRCh38/gencode.v27.geneid_to_hgnc.txt


Rscript $JS/src/coloc/annotate.coloc.R
