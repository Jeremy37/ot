#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

FNAME_QTL=$1
FNAME_MAF=$2
OUT=$3

# For colocalisation analyses we need a file of minp values per gene, a file
# with information on the variants tested, and a file with the nominal p values
# for all tests should be in a particular format:
# geneid chr pos snp_id p.value beta Bonferroni.p.value FDR MAF std.error_of_beta

zcat $FNAME_QTL | sed '1d' \
  | awk 'BEGIN{OFS="\t"}{ print $2,$3,$1 }' \
  | sort -k1,1 -k2,2n | uniq \
  | bgzip > $OUT.variant_info.partial.txt.gz

Rscript $JS/src/coloc/xQTL_add_MAF.R $OUT.variant_info.partial.txt.gz $FNAME_MAF $OUT.variant_info.txt.gz

# Get in P value per gene
(echo -e "feature\tchr\tpos\trsid\tp.value\tfeatureChromosome\tfeaturePositionStart\tSpearmanRho"; \
 zcat $FNAME_QTL | sed '1d' | getLeadSnps.pl -f stdin --genecol 4 --pcol 8 | cut -f 9 --complement \
  | awk 'BEGIN{OFS="\t"}{ print $4,$2,$3,$1,$8,$5,$6,$7 }') \
  | gzip > $OUT.qtl_signals.txt.gz

(gzhead 1 $OUT.qtl_signals.txt.gz; zcat $OUT.qtl_signals.txt.gz | awk '$5 < 1e-4') > $OUT.qtl_signals.p_lt_1e-4.txt


zcat $FNAME_QTL | sed '1d' \
  | awk 'BEGIN{OFS="\t"}{ print $4,$2,$3,$1,$8 }' \
  | sort -k2,2 -k3,3n | uniq \
  | bgzip > $OUT.nominals.txt.gz

tabix -s 2 -b 3 -e 3 $OUT.nominals.txt.gz
#
#FNAME_QTL=eQTLs_all.txt.gz
#FNAME_MAF=maf.gz
#OUT=xQTL_eQTL
#zcat $FNAME_QTL | head -n 1000 | sed '1d' | getLeadSnps.pl -f stdin --genecol 4 --pcol 8 | cut -f 9 --complement \

