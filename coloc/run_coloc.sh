JS=/Users/jeremys/work/opentargets
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

cd $JS/coloc

# First prepare variant information for AD and PD. Since we don't have MAF, put
# in a default value. MAF isn't used for coloc anyway, as we have beta + se
zcat $JS/datasets/GWAS/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.txt.gz | sed '1d' \
    | awk 'BEGIN{OFS="\t"}{print $2,$3,$1,$4,"-",0.1}' | gzip > $JS/datasets/GWAS/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.variant_info.txt.gz

# OLD
#zcat $JS/datasets/GWAS/Parkinsons_disease_Nguyen_2017_UKBB_GWAX.sorted.txt.gz | sed '1d' \
#    | awk 'BEGIN{OFS="\t"}{print $2,$3,$1,$4,"-",0.1}' | gzip > $JS/datasets/GWAS/Parkinsons_disease_Nguyen_2017_UKBB_GWAX.variant_info.txt.gz

zcat $JS/datasets/GWAS/Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.sorted.txt.gz | sed '1d' \
    | awk 'BEGIN{OFS="\t"}{print $2,$3,$1,$4,"-",0.1}' | gzip > $JS/datasets/GWAS/Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.variant_info.txt.gz


###############################################################################
# AD IGAP1
$JS/src/coloc/AD_macrophage_run_coloc.R \
        --phenotype featureCounts --gwas AD.IGAP1 \
         --overlapdist 1e5 --window 2e5 \
       --dir $JS/datasets/GWAS \
        --qtl $JS/macrophage/qtltools/output \
        --outdir $JS/coloc/ipsmacrophage \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesizes $JS/coloc/salmonella_coloc_sample_sizes.txt

submitJobs.py --MEM 6000 -j coloc_AD.IGAP1.macrophage.featureCounts -q yesterday \
    -c "$JS/src/coloc/AD_macrophage_run_coloc.R \
        --phenotype featureCounts --gwas AD.IGAP1 \
        --overlapdist 1e5 --window 2e5 \
        --dir $JS/datasets/GWAS \
        --qtl $JS/macrophage/qtltools/output \
        --outdir $JS/coloc/ipsmacrophage \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesizes $JS/coloc/salmonella_coloc_sample_sizes.txt"

submitJobs.py --MEM 6000 -j coloc_AD.IGAP1.macrophage.leafcutter -q yesterday \
    -c "$JS/src/coloc/AD_macrophage_run_coloc.R \
        --phenotype leafcutter --gwas AD.IGAP1 \
        --overlapdist 1e5 --window 2e5 \
        --dir $JS/datasets/GWAS \
        --qtl $JS/macrophage/qtltools/output \
        --outdir $JS/coloc/ipsmacrophage \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesizes $JS/coloc/salmonella_coloc_sample_sizes.txt"

submitJobs.py --MEM 12000 -j coloc_AD.IGAP1.macrophage.reviseAnnotations -q yesterday \
    -c "$JS/src/coloc/AD_macrophage_run_coloc.R \
        --phenotype reviseAnnotations --gwas AD.IGAP1 \
        --overlapdist 1e5 --window 2e5 \
        --dir $JS/datasets/GWAS \
        --qtl $JS/macrophage/qtltools/output \
        --outdir $JS/coloc/ipsmacrophage \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesizes $JS/coloc/salmonella_coloc_sample_sizes.txt"


GWASNAME=AD.meta
GWASNAME=PD.meta

###############################################################################
# AD GWAX meta-analysis
submitJobs.py --MEM 6000 -j coloc.$GWASNAME.macrophage.featureCounts -q normal \
    -c "$JS/src/coloc/AD_macrophage_run_coloc.R \
        --phenotype featureCounts --gwas $GWASNAME \
        --overlapdist 1e5 --window 2e5 \
        --dir $JS/datasets/GWAS \
        --qtl $JS/macrophage/qtltools/output \
        --outdir $JS/coloc/ipsmacrophage \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesizes $JS/coloc/salmonella_coloc_sample_sizes.txt"

submitJobs.py --MEM 6000 -j coloc.$GWASNAME.macrophage.leafcutter -q normal \
    -c "$JS/src/coloc/AD_macrophage_run_coloc.R \
        --phenotype leafcutter --gwas $GWASNAME \
        --overlapdist 1e5 --window 2e5 \
        --dir $JS/datasets/GWAS \
        --qtl $JS/macrophage/qtltools/output \
        --outdir $JS/coloc/ipsmacrophage \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesizes $JS/coloc/salmonella_coloc_sample_sizes.txt"

submitJobs.py --MEM 12000 -j coloc.$GWASNAME.macrophage.reviseAnnotations -q normal \
    -c "$JS/src/coloc/AD_macrophage_run_coloc.R \
        --phenotype reviseAnnotations --gwas $GWASNAME \
        --overlapdist 1e5 --window 2e5 \
        --dir $JS/datasets/GWAS \
        --qtl $JS/macrophage/qtltools/output \
        --outdir $JS/coloc/ipsmacrophage \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesizes $JS/coloc/salmonella_coloc_sample_sizes.txt"


###############################################################################
# Annotate sQTL genes

Rscript $JS/src/coloc/annotate_sqtl_genes.R $JS/coloc/ipsmacrophage/coloc.$GWASNAME.leafcutter.1e+05.txt 2 \
  $JS/reference/GRCh38/Homo_sapiens.GRCh38.91.exon_start_end.simple.bed \
  $JS/reference/GRCh38/gencode.v27.geneid_to_hgnc.txt


###############################################################################
# BLUEPRINT
submitJobs.py --MEM 8000 -j coloc.$GWASNAME.blueprint.mono.eqtl -q yesterday \
    -c "$JS/src/coloc/qtl_gwas_coloc.R \
        --qtlname mono_gene_eQTL --qtltype blueprint --gwasname $GWASNAME \
        --overlapdist 1e5 --window 2e5 \
        --gwasdir $JS/datasets/GWAS \
        --gwasthresh 1e-5 \
        --qtlnominal $JS/datasets/blueprint/mono_gene_eQTL.sorted.txt.gz \
        --qtlgeneminp $JS/datasets/blueprint/mono_gene_eQTL.gene_minp.txt.gz \
        --qtlvariantinfo $JS/datasets/blueprint/mono_gene_eQTL.variant_info.txt.gz \
        --gwasvariantinfo $JS/datasets/blueprint/mono_gene_eQTL.variant_info.txt.gz \
        --outdir $JS/coloc/blueprint \
        --qtlfdrthresh 0.1 \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesize 194"

submitJobs.py --MEM 8000 -j coloc.$GWASNAME.blueprint.mono.h3k27ac -q yesterday \
    -c "$JS/src/coloc/qtl_gwas_coloc.R \
        --qtlname mono_K27AC --qtltype blueprint --gwasname $GWASNAME \
        --overlapdist 1e5 --window 2e5 \
        --gwasdir $JS/datasets/GWAS \
        --gwasthresh 1e-5 \
        --qtlnominal $JS/datasets/blueprint/mono_K27AC.sorted.txt.gz \
        --qtlgeneminp $JS/datasets/blueprint/mono_K27AC.gene_minp.txt.gz \
        --qtlvariantinfo $JS/datasets/blueprint/mono_K27AC.variant_info.txt.gz \
        --outdir $JS/coloc/blueprint \
        --qtlfdrthresh 0.1 \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesize 172"

submitJobs.py --MEM 8000 -j coloc.$GWASNAME.blueprint.mono.h3k4me1 -q yesterday \
    -c "$JS/src/coloc/qtl_gwas_coloc.R \
        --qtlname mono_K4ME1 --qtltype blueprint --gwasname $GWASNAME \
        --overlapdist 1e5 --window 2e5 \
        --gwasdir $JS/datasets/GWAS \
        --gwasthresh 1e-5 \
        --qtlnominal $JS/datasets/blueprint/mono_K4ME1.sorted.txt.gz \
        --qtlgeneminp $JS/datasets/blueprint/mono_K4ME1.gene_minp.txt.gz \
        --qtlvariantinfo $JS/datasets/blueprint/mono_K4ME1.variant_info.txt.gz \
        --outdir $JS/coloc/blueprint \
        --qtlfdrthresh 0.1 \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesize 162"


###############################################################################
# xQTL

# There are about 5000 SNPs tested per gene in the xQTL dataset. They haven't
# done any FDR assessment, so I tried to do it myself using the R qvalue package
# after applying a Bonferroni correction (5000 SNPs) to each gene. This gave me
# an estimate of around 11% FDR at a p value of 1e-5, and FDR of around 1.3%
# when p = 1e-6.

submitJobs.py --MEM 3000 -j coloc.$GWASNAME.xQTL_eQTL -q normal \
    -c "$JS/src/coloc/qtl_gwas_coloc.R \
        --qtlname xQTL_eQTL --qtltype xQTL --gwasname $GWASNAME \
        --overlapdist 1e5 --window 2e5 \
        --gwasdir $JS/datasets/GWAS \
        --gwasthresh 1e-5 \
        --qtlnominal $JS/reference/xQTL/xQTL_eQTL.sorted.txt.gz \
        --qtlgeneminp $JS/reference/xQTL/xQTL_eQTL.gene_minp.txt.gz \
        --qtlvariantinfo $JS/reference/xQTL/xQTL_eQTL.variant_info.txt.gz \
        --outdir $JS/coloc/xQTL \
        --qtlpthresh 1e-5 \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesize 494"

submitJobs.py --MEM 3000 -j coloc.$GWASNAME.xQTL_haQTL -q normal \
    -c "$JS/src/coloc/qtl_gwas_coloc.R \
        --qtlname xQTL_haQTL --qtltype xQTL --gwasname $GWASNAME \
        --overlapdist 1e5 --window 2e5 \
        --gwasdir $JS/datasets/GWAS \
        --gwasthresh 1e-5 \
        --qtlnominal $JS/reference/xQTL/xQTL_haQTL.sorted.txt.gz \
        --qtlgeneminp $JS/reference/xQTL/xQTL_haQTL.gene_minp.txt.gz \
        --qtlvariantinfo $JS/reference/xQTL/xQTL_haQTL.variant_info.txt.gz \
        --outdir $JS/coloc/xQTL \
        --qtlpthresh 1e-5 \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesize 433"

submitJobs.py --MEM 4000 -j coloc.$GWASNAME.xQTL_mQTL -q normal \
    -c "$JS/src/coloc/qtl_gwas_coloc.R \
        --qtlname xQTL_mQTL --qtltype xQTL --gwasname $GWASNAME \
        --overlapdist 1e5 --window 2e5 \
        --gwasdir $JS/datasets/GWAS \
        --gwasthresh 1e-5 \
        --qtlnominal $JS/reference/xQTL/xQTL_mQTL.sorted.txt.gz \
        --qtlgeneminp $JS/reference/xQTL/xQTL_mQTL.gene_minp.txt.gz \
        --qtlvariantinfo $JS/reference/xQTL/xQTL_mQTL.variant_info.txt.gz \
        --outdir $JS/coloc/xQTL \
        --qtlpthresh 1e-5 \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesize 468"


###############################################################################
# iPSC-derived sensory neurons

submitJobs.py --MEM 3000 -j coloc.$GWASNAME.sens_neur.100k -q normal \
    -c "$JS/src/coloc/qtl_gwas_coloc.R \
        --qtlname sens_neur.100k --qtltype sens_neur --gwasname $GWASNAME \
        --overlapdist 1e5 --window 2e5 \
        --gwasdir $JS/datasets/GWAS \
        --gwasthresh 1e-5 \
        --qtlnominal $JS/sensoryneurons/GRCh38/fastqtl.nominals.cis500k.ePCs.20.coords.txt.gz \
        --qtlgeneminp $JS/sensoryneurons/GRCh38/fastqtl.permutations.10k.allgenes.cis100k.PCs.20.fdr0.1.txt  \
        --qtlvariantinfo $JS/sensoryneurons/GRCh38/imputed.97_samples.snps_indels.INFO_08.variant_info.GRCh38.gz \
        --gwasvariantinfo $JS/sensoryneurons/GRCh37/imputed.97_samples.snps_indels.INFO_08.variant_info.GRCh37.gz \
        --outdir $JS/coloc/sensoryneuron \
        --qtlfdrthresh 0.1 \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesize 97"

submitJobs.py --MEM 3000 -j coloc.$GWASNAME.sens_neur.500k -q normal \
    -c "$JS/src/coloc/qtl_gwas_coloc.R \
        --qtlname sens_neur.500k --qtltype sens_neur --gwasname $GWASNAME \
        --overlapdist 1e5 --window 2e5 \
        --gwasdir $JS/datasets/GWAS \
        --gwasthresh 1e-5 \
        --qtlnominal $JS/sensoryneurons/GRCh38/fastqtl.nominals.cis500k.ePCs.20.coords.txt.gz \
        --qtlgeneminp $JS/sensoryneurons/GRCh38/fastqtl.permutations.10k.allgenes.cis500k.PCs.20.fdr0.1.txt  \
        --qtlvariantinfo $JS/sensoryneurons/GRCh38/imputed.97_samples.snps_indels.INFO_08.variant_info.GRCh38.gz \
        --gwasvariantinfo $JS/sensoryneurons/GRCh37/imputed.97_samples.snps_indels.INFO_08.variant_info.GRCh37.gz \
        --outdir $JS/coloc/sensoryneuron \
        --qtlfdrthresh 0.1 \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesize 97"

submitJobs.py --MEM 3000 -j coloc.$GWASNAME.sens_neur.sqtl -q normal \
    -c "$JS/src/coloc/qtl_gwas_coloc.R \
        --qtlname sens_neur.sqtl --qtltype sens_neur --gwasname $GWASNAME \
        --overlapdist 1e5 --window 2e5 \
        --gwasdir $JS/datasets/GWAS \
        --gwasthresh 1e-5 \
        --qtlnominal $JS/sensoryneurons/GRCh38/fastqtl.sqtl.nominals.PCs.5.pvalues.coords.txt.gz \
        --qtlgeneminp $JS/sensoryneurons/GRCh38/fastqtl.sqtl.permutations.10k.PCs.5.fdr0.1.txt  \
        --qtlvariantinfo $JS/sensoryneurons/GRCh38/imputed.97_samples.snps_indels.INFO_08.variant_info.GRCh38.gz \
        --gwasvariantinfo $JS/sensoryneurons/GRCh37/imputed.97_samples.snps_indels.INFO_08.variant_info.GRCh37.gz \
        --outdir $JS/coloc/sensoryneuron \
        --qtlfdrthresh 0.1 \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesize 97"



