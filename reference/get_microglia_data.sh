#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
PATH=~js29/.aspera/connect/bin:$PATH

cd $JS/datasets/microglia/ncbi/dbGaP-15880

# Instructions on downloading SRA data are here:
# https://www.ncbi.nlm.nih.gov/books/NBK36439/

# Following the instructions to use vdb-config graphically did not work,
# but doing this worked:
$JS/software/sratoolkit.2.9.0-ubuntu64/bin/vdb-config --import $JS/datasets/microglia/prj_15880.ngc

# Use SRA toolkit's prefetch to get files specified in the .krt file
KRT_FILE=$JS/datasets/microglia/cart_DAR60784_201804131119.krt
$JS/software/sratoolkit.2.9.0-ubuntu64/bin/prefetch $KRT_FILE

#$JS/software/sratoolkit.2.9.0-ubuntu64/bin/prefetch SRR5955079

#$JS/software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump --gzip -O $JS/datasets/microglia/fastq $JS/datasets/microglia/ncbi/dbGaP-15880/sra/SRR5955082.sra 
#$JS/software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump --gzip SRR5955082

mkdir $JS/datasets/microglia/fastq
submitJobs.py --MEM 1000 -j fastq-dump -q yesterday -c "bash $JS/src/reference/microglia_fastq_dump.sh"


$JS/software/sratoolkit.2.9.0-ubuntu64/bin/vdb-decrypt /lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/datasets/microglia/PhenoGenotypeFiles/RootStudyConsentSet_phs001373.HumanMicroglia.v1.p1.c1.HMB
