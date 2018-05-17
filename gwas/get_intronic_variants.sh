#!/bin/bash

ROOT=/lustre/scratch115/realdata/mdt3/projects/otcoregen
OUTDIR=$ROOT/jeremys/gwas/AD/finemap
cd $OUTDIR

(head -n 1 AD.finemap.max_causal_2.snp; sed '1d' AD.finemap.max_causal_2.snp | awk '$1 ~ /2_127892810/') > AD.finemap.BIN1.snp

BASE=AD.finemap.BIN1
# Get SNP IDs to pass to VEP for annotation
sed '1d' $BASE.snp | cut -f 3 > $BASE.rsids.txt
# Pass this file to VEP - run manually in web browser
# NOTE: HAVE TO CHANGE FROM WINDOWS CR-LF line endings to Unix LF

# Annotate the finemap causal probability for each SNP
hashJoin.pl --hashFile $BASE.snp --scanFile $BASE.vep.txt --colHashFile 3 --colScanFile 1 --header | gzip > $BASE.vep_annotated.txt.gz

# Select intronic variants from the VEP output
(gzhead 1 $BASE.vep_annotated.txt.gz; zcat $BASE.vep_annotated.txt.gz | awk '$4 ~ /intron_variant/') | gzip > $BASE.vep.intronic.txt.gz
(gzhead 1 $BASE.vep_annotated.txt.gz; zcat $BASE.vep_annotated.txt.gz | awk '$4 ~ /splice/') | gzip > $BASE.vep.splice.txt.gz
(gzhead 1 $BASE.vep_annotated.txt.gz; zcat $BASE.vep_annotated.txt.gz | awk '$4 ~ /intron_variant|UTR|synonymous|missense|start_lost|stop_gained|splice/') | gzip > $BASE.vep.transcribed.txt.gz

TYPE=intronic
TYPE=splice
TYPE=transcribed
zcat $BASE.vep.$TYPE.txt.gz | cut -f 1-4,6-12,55-60,62 | awk '$18 > 0' > $BASE.vep.$TYPE.pcausal_gt_0.txt
hashJoin.pl --hashFile ../AD.locusnames.txt --scanFile $BASE.vep.$TYPE.pcausal_gt_0.txt --colHashFile 2 --colScanFile 17 --header > $BASE.vep.$TYPE.pcausal_gt_0.locusnames.txt
#hashJoin.pl --hashFile $ROOT/jeremys/reference/tissueRPKM/tissues.selected.rpkm_average.txt --scanFile AD.finemap.vep.$TYPE.pcausal_gt_0.locusnames.txt --colHashFile 1 --colScanFile 6 --header > AD.finemap.vep.$TYPE.pcausal_gt_0.locusnames.rpkms.txt



################################################################################
# Do the same for PD
ROOT=/lustre/scratch115/realdata/mdt3/projects/otcoregen
OUTDIR=$ROOT/jeremys/gwas/PD/finemap
cd $OUTDIR

sed '1d' pd_gwas_gwax.combined.set | cut -f 3 > PD.finemap.rsids.txt
(head -n 1 pd_gwas_gwax.combined.set; sed '1d' pd_gwas_gwax.combined.set | awk '$8 > 0') > pd_gwas_gwax.combined.pcausal_gt_0.set
sed '1d' pd_gwas_gwax.combined.pcausal_gt_0.set | awk '$8 > 0' | cut -f 3 > PD.finemap.pcausal_gt_0.rsids.txt
# Pass this file to VEP - run manually in web browser
# Select default options, except:
# - include 1000 genomes continental allele frequencies
# - include exon/intron numbers
# - identify canonical transcripts
# - include the 2 possible splicing predictions (dbscSNV, MaxEntScan)
# NOTE: HAVE TO CHANGE FROM WINDOWS CR-LF line endings to Unix LF

# Annotate the finemap causal probability for each SNP
hashJoin.pl --hashFile pd_gwas_gwax.combined.pcausal_gt_0.set --scanFile PD.finemap.pcausal_gt_0.vep.txt --colHashFile 3 --colScanFile 1 --header | gzip > PD.finemap.pcausal_gt_0.vep_annotated.txt.gz
hashJoin.pl --hashFile pd_gwas_gwax.combined.set --scanFile PD.finemap.vep.txt --colHashFile 3 --colScanFile 1 --header | gzip > PD.finemap.vep_annotated.txt.gz

# Select intronic variants from the VEP output
(gzhead 1 PD.finemap.vep_annotated.txt.gz; zcat PD.finemap.vep_annotated.txt.gz | awk '$4 ~ /intron_variant/') | gzip > PD.finemap.vep.intronic.txt.gz
(gzhead 1 PD.finemap.vep_annotated.txt.gz; zcat PD.finemap.vep_annotated.txt.gz | awk '$4 ~ /splice/') > PD.finemap.vep.splice.txt
(gzhead 1 PD.finemap.vep_annotated.txt.gz; zcat PD.finemap.vep_annotated.txt.gz | awk '$4 ~ /intron_variant|UTR|synonymous|missense|start_lost|stop_gained/') | gzip > PD.finemap.vep.transcribed.txt.gz

TYPE=intronic
TYPE=splice
TYPE=transcribed
gzhead 1 PD.finemap.vep.$TYPE.txt.gz > PD.finemap.vep.$TYPE.txt.header
zcat PD.finemap.vep.$TYPE.txt.gz | cut -f 1-4,6-12,56-67 | awk '$23 > 0' > PD.finemap.vep.$TYPE.pcausal_gt_0.txt
#hashJoin.pl --hashFile $ROOT/jeremys/reference/tissueRPKM/tissues.selected.rpkm_average.txt --scanFile PD.finemap.vep.$TYPE.pcausal_gt_0.txt --colHashFile 1 --colScanFile 6 --header > PD.finemap.vep.$TYPE.pcausal_gt_0.rpkms.txt



Rscript $ROOT/jeremys/src/gwas/get_intronic_variants.rpkms.R

