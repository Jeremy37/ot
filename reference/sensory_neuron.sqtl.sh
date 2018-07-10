#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

cd $JS/sensoryneurons/GRCh38


(zcat /warehouse/compgen_wh05/js29/sensoryneurons/leafcutter/snqtl.v5.clusters_perind_numers.ratios.txt.gz | head -n 1 | sed 's/chrom/cluster/' | sed 's/.qtl.bam//g';
 zcat /warehouse/compgen_wh05/js29/sensoryneurons/leafcutter/snqtl.v5.clusters_perind_numers.ratios.txt.gz | sed '1d') \
 | gzip > sqtl/sqtl.clusters_perind_numers.ratios.txt.gz

# Extract genotype counts from the VCF file of sensory neuron lead SNPs
python vcfToGenotypeCounts.py --vcf all.leadSNPs.vcf.gz | gzip > all.leadSNPs.genotype_counts.txt.gz

cd $JS/sensoryneurons/GRCh38/sqtl
Rscript $JS/src/misc/plotSensoryneuronSqtl.R --args cluster=clu_55409 snpid=rs4147914 genotypeCounts=../all.leadSNPs.genotype_counts.txt.gz genotypeMap=../sensoryneuron.merged_sample_genotype_map.txt clusterPerind=sqtl.clusters_perind_numers.ratios.txt.gz nominals=sqtl.fastqtl.nominals.txt.gz

Rscript $JS/src/misc/plotSensoryneuronSqtl.R --args cluster=clu_48150 snpid=rs118189385 genotypeCounts=../all.leadSNPs.genotype_counts.txt.gz genotypeMap=../sensoryneuron.merged_sample_genotype_map.txt clusterPerind=sqtl.clusters_perind_numers.ratios.txt.gz nominals=sqtl.fastqtl.nominals.txt.gz

Rscript $JS/src/misc/plotSensoryneuronSqtl.R --args cluster=clu_28666 snpid=rs700828 genotypeCounts=../all.leadSNPs.genotype_counts.txt.gz genotypeMap=../sensoryneuron.merged_sample_genotype_map.txt clusterPerind=sqtl.clusters_perind_numers.ratios.txt.gz nominals=sqtl.fastqtl.nominals.txt.gz

Rscript $JS/src/misc/plotSensoryneuronSqtl.R --args cluster=clu_17169 snpid=rs13432074 genotypeCounts=../all.leadSNPs.genotype_counts.txt.gz genotypeMap=../sensoryneuron.merged_sample_genotype_map.txt clusterPerind=sqtl.clusters_perind_numers.ratios.txt.gz nominals=sqtl.fastqtl.nominals.txt.gz

Rscript $JS/src/misc/plotSensoryneuronSqtl.R --args cluster=clu_9125 snpid=rs35186494 genotypeCounts=../all.leadSNPs.genotype_counts.txt.gz genotypeMap=../sensoryneuron.merged_sample_genotype_map.txt clusterPerind=sqtl.clusters_perind_numers.ratios.txt.gz nominals=sqtl.fastqtl.nominals.txt.gz
