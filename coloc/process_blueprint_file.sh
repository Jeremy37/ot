#!/bin/bash
FNAME=$1
OUT=$2

# For colocalisation analyses we need a file of minp values per gene, a file
# with information on the variants tested, and a file with the nominal p values
# for all tests should be in a particular format:
# geneid chr pos snp_id p.value beta Bonferroni.p.value FDR MAF std.error_of_beta

(echo -e "chr\tpos\trsid\tMAF"; \
 zcat $FNAME \
  | perl -ane '@fields=split(/:|_/,$F[0]); print join("\t", @fields[0..1], $F[1], $F[7])."\n";' \
  | sort -k1,1 -k2,2n | uniq) \
  | bgzip > $OUT.variant_info.txt.gz

# Get in P value per gene
(echo -e "feature\tchr\tpos\trsid\tp.value\tbeta\tBonferroni.p.value\tFDR\talt_allele_frequency\tstd.error_of_beta"; \
 getLeadSnps.pl -f $FNAME --genecol 3 --pcol 4 --sep " " | cut -f 10 --complement \
  | perl -ane '@fields=split(/:|_/,$F[0]); print join("\t", $F[2], @fields[0..1], $F[1], @F[3..8])."\n";') \
  | gzip > $OUT.qtl_signals.txt.gz

(gzhead 1 $OUT.qtl_signals.txt.gz; zcat $OUT.qtl_signals.txt.gz | awk '$8 < 0.05') > $OUT.qtl_signals.FDR_0.05.txt

zcat $FNAME \
  | perl -ane '@fields=split(/:|_/,$F[0]); print join("\t", $F[2], @fields[0..1], $F[1], $F[3])."\n";' \
  | sort -k2,2 -k3,3n | uniq \
  | bgzip > $OUT.nominals.txt.gz

tabix -s 2 -b 3 -e 3 $OUT.nominals.txt.gz
