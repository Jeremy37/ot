#!/bin/bash

JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

STAR --runThreadN 6 --runMode genomeGenerate \
--genomeDir $JS/reference/GRCh38/STAR_gencode_v27_index \
--genomeFastaFiles $JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa \
--sjdbGTFfile $JS/reference/GRCh38/gencode.v27.annotation.gtf
--sjdbOverhang 149
