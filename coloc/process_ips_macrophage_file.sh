#!/bin/bash
FNAME=$1
OUT=$2

zcat $FNAME \
  | perl -ane 'print join("\t", $F[0],$F[8],$F[9],$F[7],$F[11])."\n";' \
  | sort -k2,2 -k3,3n | uniq \
  | bgzip > $OUT.nominals.txt.gz

tabix -s 2 -b 3 -e 3 $OUT.nominals.txt.gz
