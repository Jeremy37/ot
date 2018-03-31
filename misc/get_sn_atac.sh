#!/bin/bash
# This script is to get candidate causal ATAC QTL SNPs, which can then be
# used for Sarah's CRISPR assays for ATAC-altering variants
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

SN=/warehouse/compgen_wh05/js29/sensoryneurons
SNQTL=$SN/snqtl/results

cd $JS/sensoryneurons/GRCh38/ATAC

gzhead 100 $SNQTL/atac/rasqual.1k.all.with_header.ppa.gz > rasqual.1k.all.with_header.head.txt
(gzhead 1 $SNQTL/atac/rasqual.1k.all.with_header.ppa.gz; \
 zcat $SNQTL/atac/rasqual.1k.all.with_header.ppa.gz | awk '$26 < 0.01 && $27 > 0.0001') \
| gzip > rasqual.1k.pthresh0.01.ppathresh.0.0001.txt

cp $SNQTL/atac/rasqual.1k.leadSNPs.fdr0.1.ann.txt .
cp /warehouse/compgen_wh05/js29/sensoryneurons/atac/rasqual/output/rasqual.1k.fdr0.1.peakcoords.txt .


IPSC_PEAKS=$JS/ipsneurons/GRCh38/ATAC/peaks/atac_ipsc_peaks.narrowPeak
cp /warehouse/compgen_wh05/js29/sensoryneurons/atac/ATAC_consensus_peaks.metadata.txt sensoryneuron_consensus_atac_peaks.GRCh38.txt

Rscript get_sn_atac.R \
  rasqual.1k.leadSNPs.fdr0.1.ann.txt \
  rasqual.1k.pthresh0.01.ppathresh.0.0001.txt \
  sensoryneuron_consensus_atac_peaks.GRCh38.bed \
  $JS/ipsneurons/GRCh38/ATAC/peaks/atac_ipsc_peaks.narrowPeak


# Check the genotypes of these SNPs in kolf2
zcat $JS/ipsneurons/GRCh38/genotypes/kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.biallelic.snps.chr.vcf.gz | grep rs117481827
zcat $JS/ipsneurons/GRCh38/genotypes/kolf_2.imputed_phased.20150604.GRCh38.vcf.gz | grep rs117481827


