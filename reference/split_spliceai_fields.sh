#!/bin/bash

(echo -e 'chr\tpos\tref\talt\tsymbol\tstrand\ttype\tdist\tDS_AG\tDS_AL\tDS_DG\tDS_DL\tDP_AG\tDP_AL\tDP_DG\tDP_DL';
 zcat spliceai.merge.vcf.bgz | grep -v '^#' \
 | perl -ane 'if ($F[7] =~ /SYMBOL=([^;]*);STRAND=([+|-]);TYPE=([E|I]);DIST=([^;]*);DS_AG=([^;]*);DS_AL=([^;]*);DS_DG=([^;]*);DS_DL=([^;]*);DP_AG=([^;]*);DP_AL=([^;]*);DP_DG=([^;]*);DP_DL=(-*.*)/) { print join("\t", $F[0], $F[1], $F[3], $F[4], $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12)."\n"; }') \
 | bgzip > spliceai.merge.tsv.bgz

tabix -s 1 -b 2 -e 2 -S 1 spliceai.merge.tsv.bgz
