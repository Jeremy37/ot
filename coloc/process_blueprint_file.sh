#!/bin/bash
FNAME=$1
OUT=$2

# For colocalisation analyses we need a file of minp values per gene, a file
# with information on the variants tested, and a file with the nominal p values
# for all tests should be in a particular format:
# geneid chr pos snp_id p.value beta Bonferroni.p.value FDR MAF std.error_of_beta

zcat $FNAME \
  | perl -ane '@fields=split(/:|_/,$F[0]); print join("\t", @fields[0..1], $F[1], @fields[2..3], $F[7])."\n";' \
  | sort -k1,1 -k2,2n | uniq \
  | bgzip > $OUT.variant_info.txt.gz

# Get in P value per gene
getLeadSnps.pl -f $FNAME --genecol 3 --pcol 4 --sep " " | cut -f 10 --complement \
  | perl -ane '@fields=split(/:|_/,$F[0]); print join("\t", $F[2], @fields[0..1], $F[1], @F[3..8])."\n";' \
  | gzip > $OUT.gene_minp.txt.gz

zcat $FNAME \
  | perl -ane '@fields=split(/:|_/,$F[0]); print join("\t", $F[2], @fields[0..1], $F[1], @F[3..8])."\n";' \
  | sort -k2,2 -k3,3n | uniq \
  | bgzip > $OUT.sorted.txt.gz

tabix -s 2 -b 3 -e 3 $OUT.sorted.txt.gz

