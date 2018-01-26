JS=/Users/jeremys/work/opentargets
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys


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



###############################################################################
# AD GWAX meta-analysis
submitJobs.py --MEM 6000 -j coloc_AD.meta.macrophage.featureCounts -q yesterday \
    -c "$JS/src/coloc/AD_macrophage_run_coloc.R \
        --phenotype featureCounts --gwas AD.meta \
        --overlapdist 1e5 --window 2e5 \
        --dir $JS/datasets/GWAS \
        --qtl $JS/macrophage/qtltools/output \
        --outdir $JS/coloc/ipsmacrophage \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesizes $JS/coloc/salmonella_coloc_sample_sizes.txt"

submitJobs.py --MEM 6000 -j coloc_AD.meta.macrophage.leafcutter -q normal \
    -c "$JS/src/coloc/AD_macrophage_run_coloc.R \
        --phenotype leafcutter --gwas AD.meta \
        --overlapdist 1e5 --window 2e5 \
        --dir $JS/datasets/GWAS \
        --qtl $JS/macrophage/qtltools/output \
        --outdir $JS/coloc/ipsmacrophage \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesizes $JS/coloc/salmonella_coloc_sample_sizes.txt"

submitJobs.py --MEM 12000 -j coloc_AD.meta.macrophage.reviseAnnotations -q yesterday \
    -c "$JS/src/coloc/AD_macrophage_run_coloc.R \
        --phenotype reviseAnnotations --gwas AD.meta \
        --overlapdist 1e5 --window 2e5 \
        --dir $JS/datasets/GWAS \
        --qtl $JS/macrophage/qtltools/output \
        --outdir $JS/coloc/ipsmacrophage \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesizes $JS/coloc/salmonella_coloc_sample_sizes.txt"


###############################################################################
# Annotate sQTL genes

Rscript annotate_sqtl_genes.R $JS/coloc/ipsmacrophage/coloc.AD.meta.leafcutter.1e+05.txt 2 \
  $JS/reference/GRCh38/Homo_sapiens.GRCh38.91.gene_start_end.bed \
  $JS/reference/GRCh38/Homo_sapiens.GRCh38.91.exon_start_end.bed \
  > $JS/coloc/ipsmacrophage/coloc.AD.meta.leafcutter.1e+05.ann.txt


###############################################################################

submitJobs.py --MEM 8000 -j coloc_AD.meta.blueprint.mono -q yesterday \
    -c "$JS/src/coloc/qtl_gwas_coloc.blueprint.R \
        --phenotype mono_gene_eQTL --gwas AD.meta \
        --overlapdist 1e5 --window 2e5 \
        --dir $JS/datasets/GWAS \
        --qtl $JS/datasets/blueprint \
        --outdir $JS/coloc/blueprint \
        --p1 1e-3 --p2 1e-4 --p12 5e-5 \
        --samplesize 194"


