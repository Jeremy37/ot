#!/bin/bash

JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
SPLICEAI=$JS/datasets/SpliceAI

cd $SPLICEAI

zcat $SPLICEAI/whole_genome_filtered_spliceai_scores.vcf.gz | bgzip > $SPLICEAI/whole_genome_filtered_spliceai_scores.vcf.bgz
zcat $SPLICEAI/exome_spliceai_scores.vcf.gz | bgzip > $SPLICEAI/exome_spliceai_scores.vcf.bgz

# Merged together whole genome (filtered) and exome (all SNP) predictions
submitJobs.py --MEM 1000 -j merge_genome_exome -q yesterday \
  -c "(zcat $SPLICEAI/whole_genome_filtered_spliceai_scores.vcf.gz | head -n 15; \
  (zcat $SPLICEAI/whole_genome_filtered_spliceai_scores.vcf.gz; zcat $SPLICEAI/exome_spliceai_scores.vcf.gz) \
    | grep -v '^#' | sort -k1,1 -k2,2n | uniq) \
  | bgzip > spliceai.merge.vcf.bgz"

tabix -p vcf spliceai.merge.vcf.bgz

(zcat $SPLICEAI/whole_genome_filtered_spliceai_scores.vcf.gz | head -n 15; \
  (zcat $SPLICEAI/whole_genome_filtered_spliceai_scores.vcf.gz | head -n 1000; zcat $SPLICEAI/exome_spliceai_scores.vcf.gz | head -n 10000) \
    | grep -v '^#' | sort -k1,1 -k2,2n | uniq) \
  | bgzip > spliceai.merge.test.vcf.bgz
  
(echo -e "chr\tpos\tref\talt\tsymbol\tstrand\ttype\tdist\tDS_AG\tDS_AL\tDS_DG\tDS_DL\tDP_AG\tDP_AL\tDP_DG\tDP_DL";
 zcat spliceai.merge.test.vcf.bgz | grep -v "^#" \
 | perl -ane "if ($F[7] =~ /SYMBOL=([^;]*);STRAND=([+|-]);TYPE=([E|I]);DIST=([^;]*);DS_AG=([^;]*);DS_AL=([^;]*);DS_DG=([^;]*);DS_DL=([^;]*);DP_AG=([^;]*);DP_AL=([^;]*);DP_DG=([^;]*);DP_DL=(-*.*)/) { print join(\"\t\", $F[0], $F[1], $F[3], $F[4], $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12).\"\n\"; }") \
 > spliceai.merge.test.tsv


submitJobs.py --MEM 1000 -j split_spliceai_fields -q yesterday \
  -c "$JS/src/reference/split_spliceai_fields.sh"

tabix -s 1 -b 2 -e 2 -S 1 spliceai.merge.tsv.bgz

tabix spliceai.merge.tsv.bgz 2:127792810-127992810 | bgzip > spliceai.merge.test.tsv.bgz
submitJobs.py --MEM 2000 -j get_spliceai_ucsc_tracks -q yesterday \
  -c "python $JS/src/reference/get_spliceai_ucsc_tracks.py"

fetchChromSizes hg19 > chromsizes.hg19.txt
wigToBigWig spliceai.ucsc.ag.wig chromsizes.hg19.txt spliceai.ucsc.ag.bw
wigToBigWig spliceai.ucsc.al.wig chromsizes.hg19.txt spliceai.ucsc.al.bw
wigToBigWig spliceai.ucsc.dg.wig chromsizes.hg19.txt spliceai.ucsc.dg.bw
wigToBigWig spliceai.ucsc.dl.wig chromsizes.hg19.txt spliceai.ucsc.dl.bw
wigToBigWig spliceai.ucsc.merged.wig chromsizes.hg19.txt spliceai.ucsc.merged.bw


tabix $SPLICEAI/spliceai.merge.vcf.bgz 10:1000000-1001000

(gzhead 1 $SPLICEAI/spliceai.merge.tsv.bgz; \
 tabix $SPLICEAI/spliceai.merge.tsv.bgz 15:51016500-51019500) \
 > spliceai/SPPL2A.spliceai.51016500-51019500.tsv


cd $JS/AD_finemap
mkdir spliceai
#echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > ad.spliceai.regions.tsv
echo -e "chr\tpos\tref\talt\tsymbol\tstrand\ttype\tdist\tDS_AG\tDS_AL\tDS_DG\tDS_DL\tDP_AG\tDP_AL\tDP_DG\tDP_DL" > ad.spliceai.regions.tsv
python $JS/src/misc/ad_get_spliceai_annotations.py
gzip spliceai/ad.spliceai.regions.tsv

cd spliceai
# In R:
readr::read_tsv("ad.spliceai.regions.tsv.gz") %>%
  group_by(chr, pos) %>%
  summarise(max_ag = max(DS_AG), max_al = max(DS_AL), max_dg = max(DS_DG), max_dl = max(DS_DL)) %>%
  write.table(file="ad.spliceai.by_pos.tsv2", sep="\t", quote=F, col.names=T, row.names=F)
gzip ad.spliceai.by_pos.tsv

# Convert the TSV file into separate track hubs for UCSC
python $JS/src/misc/ad_get_spliceai.ucsc_tracks.py

wigToBigWig spliceai.ucsc.ag.wig chromsizes.hg19.txt spliceai.ucsc.ag.bw
wigToBigWig spliceai.ucsc.al.wig chromsizes.hg19.txt spliceai.ucsc.al.bw
wigToBigWig spliceai.ucsc.dg.wig chromsizes.hg19.txt spliceai.ucsc.dg.bw
wigToBigWig spliceai.ucsc.dl.wig chromsizes.hg19.txt spliceai.ucsc.dl.bw
wigToBigWig spliceai.ucsc.merged.wig chromsizes.hg19.txt spliceai.ucsc.merged.bw

#echo "track type=wiggle_0 name=SpliceAI-SPPL2A-wig\nvariableStep chrom=chr15 color=255,0,0" > spliceai/spliceai.ad_regions.ag.wig
#zcat ad.spliceai.regions.tsv | sed '1d' | cut -f 2,9 > spliceai/spliceai.ad_regions.ag.wig


