#!/bin/bash

ROOT=/lustre/scratch115/realdata/mdt3/projects/otcoregen
OUTDIR=$ROOT/jeremys/gwas/AD/finemap
cd $OUTDIR

# Get SNP IDs to pass to VEP for annotation
sed '1d' AD.finemap.snp | cut -f 3 > AD.finemap.rsids.txt
# Pass this file to VEP - run manually in web browser
# NOTE: HAVE TO CHANGE FROM WINDOWS CR-LF line endings to Unix LF

# Annotate the finemap causal probability for each SNP
hashJoin.pl --hashFile AD.finemap.snp --scanFile AD.finemap.vep.txt --colHashFile 3 --colScanFile 1 --header | gzip > AD.finemap.vep_annotated.txt.gz

# Select intronic variants from the VEP output
(gzhead 1 AD.finemap.vep_annotated.txt.gz; zcat AD.finemap.vep_annotated.txt.gz | awk '$4 ~ /intron_variant/') | gzip > AD.finemap.vep.intronic.txt.gz
(gzhead 1 AD.finemap.vep_annotated.txt.gz; zcat AD.finemap.vep_annotated.txt.gz | awk '$4 ~ /splice/') | gzip > AD.finemap.vep.splice.txt.gz
(gzhead 1 AD.finemap.vep_annotated.txt.gz; zcat AD.finemap.vep_annotated.txt.gz | awk '$4 ~ /intron_variant|UTR|synonymous|missense|start_lost|stop_gained|splice/') | gzip > AD.finemap.vep.transcribed.txt.gz

TYPE=intronic
TYPE=splice
TYPE=transcribed
zcat AD.finemap.vep.$TYPE.txt.gz | cut -f 1-4,6-12,55-60 | awk '$17 > 0' > AD.finemap.vep.$TYPE.pcausal_gt_0.txt
hashJoin.pl --hashFile ../AD.locusnames.txt --scanFile AD.finemap.vep.$TYPE.pcausal_gt_0.txt --colHashFile 2 --colScanFile 15 --header > AD.finemap.vep.$TYPE.pcausal_gt_0.locusnames.txt
#hashJoin.pl --hashFile $ROOT/jeremys/reference/tissueRPKM/tissues.selected.rpkm_average.txt --scanFile AD.finemap.vep.$TYPE.pcausal_gt_0.locusnames.txt --colHashFile 1 --colScanFile 6 --header > AD.finemap.vep.$TYPE.pcausal_gt_0.locusnames.rpkms.txt


################################################################################
# Do the same for PD
ROOT=/lustre/scratch115/realdata/mdt3/projects/otcoregen
OUTDIR=$ROOT/jeremys/gwas/PD/finemap
cd $OUTDIR

sed '1d' PD.finemap.snp | cut -f 3 > PD.finemap.rsids.txt
# Pass this file to VEP - run manually in web browser

# Annotate the finemap causal probability for each SNP
hashJoin.pl --hashFile PD.finemap.snp --scanFile PD.finemap.vep.txt --colHashFile 3 --colScanFile 1 --header | gzip > PD.finemap.vep_annotated.txt.gz

# Select intronic variants from the VEP output
(gzhead 1 PD.finemap.vep_annotated.txt.gz; zcat PD.finemap.vep_annotated.txt.gz | awk '$4 ~ /intron_variant/') | gzip > PD.finemap.vep.intronic.txt.gz
(gzhead 1 PD.finemap.vep_annotated.txt.gz; zcat PD.finemap.vep_annotated.txt.gz | awk '$4 ~ /splice/') > PD.finemap.vep.splice.txt
(gzhead 1 PD.finemap.vep_annotated.txt.gz; zcat PD.finemap.vep_annotated.txt.gz | awk '$4 ~ /intron_variant|UTR|synonymous|missense|start_lost|stop_gained/') | gzip > PD.finemap.vep.transcribed.txt.gz

TYPE=intronic
TYPE=splice
TYPE=transcribed
zcat PD.finemap.vep.$TYPE.txt.gz | cut -f 1-4,6-12,56-61 | awk '$17 > 0' > PD.finemap.vep.$TYPE.pcausal_gt_0.txt
#hashJoin.pl --hashFile $ROOT/jeremys/reference/tissueRPKM/tissues.selected.rpkm_average.txt --scanFile PD.finemap.vep.$TYPE.pcausal_gt_0.txt --colHashFile 1 --colScanFile 6 --header > PD.finemap.vep.$TYPE.pcausal_gt_0.rpkms.txt



Rscript $ROOT/jeremys/src/gwas/get_intronic_variants.rpkms.R

