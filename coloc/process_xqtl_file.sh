#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

FNAME_QTL=$1
FNAME_MAF=$2
OUT=$3

# For colocalisation analyses we need a file of minp values per gene, a file
# with information on the variants tested, and a file with the nominal p values
# for all tests should be in a particular format:
# geneid chr pos snp_id p.value beta Bonferroni.p.value FDR MAF std.error_of_beta

# The ref/alt alleles are expected in the variant info file, but they aren't
# actually used. Since we don't have them, we put "-" in those fields.
zcat $FNAME_QTL | sed '1d' \
  | awk 'BEGIN{OFS="\t"}{ print $2,$3,$1 }' \
  | sort -k1,1 -k2,2n | uniq \
  | bgzip > $OUT.variant_info.partial.txt.gz

Rscript $JS/src/coloc/xQTL_add_MAF.R $OUT.variant_info.partial.txt.gz $FNAME_MAF $OUT.variant_info.txt.gz

# Get in P value per gene
zcat $FNAME_QTL | sed '1d' | getLeadSnps.pl -f stdin --genecol 4 --pcol 8 | cut -f 9 --complement \
  | awk 'BEGIN{OFS="\t"}{ print $4,$2,$3,$1,$5,$6,$7,$8 }' \
  | gzip > $OUT.gene_minp.txt.gz

zcat $FNAME_QTL | sed '1d' \
  | awk 'BEGIN{OFS="\t"}{ print $4,$2,$3,$1,$5,$6,$7,$8 }' \
  | sort -k2,2 -k3,3n | uniq \
  | bgzip > $OUT.sorted.txt.gz

tabix -s 2 -b 3 -e 3 $OUT.sorted.txt.gz
#
#FNAME_QTL=eQTLs_all.txt.gz
#FNAME_MAF=maf.gz
#OUT=xQTL_eQTL
#zcat $FNAME_QTL | head -n 1000 | sed '1d' | getLeadSnps.pl -f stdin --genecol 4 --pcol 8 | cut -f 9 --complement \
