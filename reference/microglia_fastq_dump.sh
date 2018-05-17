#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

for f in $JS/datasets/microglia/ncbi/dbGaP-15880/sra/*.sra; do
    filename=$(basename "$f")
    filebase="${filename%.*}"
    if [ ! -f $JS/datasets/microglia/fastq/$filebase.fastq.gz ]; then
        $JS/software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump --gzip -O $JS/datasets/microglia/fastq $f
    fi
done
