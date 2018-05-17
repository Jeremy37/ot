#!/bin/bash
#newgrp otcoregen # Set primary group to otcoregen so new files have this group by default
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
SEQ=$JS/ipsneurons/GRCh38
SRC=$JS/src/ipsneurons
source $JS/software/ipsneurons.software.sh
kinit # login to iRODS

###############################################################################
## RNA-seq
cd $SEQ/RNA

# Example command to get metadata - this is used within the scripts below
# imeta qu -z seq -d study = 'OpenTargets AD PD finemap ATAC' and target = 1 and manual_qc = 1
python ~/src/utils/irods/irodsGetSamplesInStudy.py --studyName "OpenTargets AD PD finemap RNA" |  cut -f1 -d "." | uniq > irods.lanelets.txt
# Manually extract those corresponding to new samples into irods.lanelets.psen.txt

# Map lanelet ids to sample ids
# imeta ls -d /seq/23846/23846_8#6.cram
python ~/src/utils/irods/irodsFetchMeta.py --irodsList irods.lanelets.psen.txt | sort -k1 > irods.sample_lanes.psen.txt 

# Fetch lanelets in cram format from irods
# iget /seq/23846/23846_8#6.cram
cut -f1 irods.lanelets.psen.txt | python ~/src/utils/irods/fetch-irods.py --dir cram/ --suffix .cram

for f in cram/*.cram; do
  if [ ! -f $f.crai ]; then
    samtools index $f
  fi
done

# Extract the "genome" file for bedtools, i.e. chromosome lengths, from the CRAM header
###samtools view -H "cram/23882_5#1.cram" | ~/src/utils/bam/genomeFileFromBamHeader.py > GRCh38.genome.txt
###samtools view -H "cram/25563_1#1.cram" | ~/src/utils/bam/genomeFileFromBamHeader.py > GRCh38.genome.psen.txt
# diff GRCh38.genome.psen.txt GRCh38.genome.txt # These two are the same

# Convert from CRAM to BigWig format for easier viewing in IGV
cat irods.sample_lanes.psen.txt | submitJobs.py --MEM 2000 --jobname mergeBams \
    --command "~/src/utils/bam/mergeBams.py --indir cram --outdir . --insuffix .cram --outsuffix .bam"

cut -f 1 irods.sample_lanes.psen.txt | submitJobs.py --MEM 1000 --jobname indexBams \
    --command "~/src/utils/bam/indexBams.py --bamdir . --insuffix .bam"

cut -f 1 irods.sample_lanes.psen.txt | submitJobs.py --MEM 3000 --jobname bamCoverageToBigWig \
    --command "~/src/utils/coverage/bam2bigwig.py --genome GRCh38.genome.txt --indir . --split"

# Ideally we would use verifyBamID to ensure that the files match the expected
# sample, but I don't know the genetic background for the PSEN line, and don't
# yet have genotypes for the BOB line


# Count reads mapping to each chromosome
cut -f1 irods.sample_lanes.psen.txt | submitJobs.py --MEM 3000 -n 2 -j countReadsPerChr -c "python ~js29/src/utils/bam/countReadsPerChr.py --indir . --outdir . --insuffix .bam"
echo "4860STDY7352094" | submitJobs.py --MEM 3000 -n 2 -j countReadsPerChr -c "python ~js29/src/utils/bam/countReadsPerChr.py --indir . --outdir . --insuffix .bam"

# Make a file with chromosomes 1...X,Y,MT named chrcounts.all.txt
f=./4860STDY7352093/4860STDY7352093.chr_counts
cat $f | sed "s/^[ ]*//" | perl -ne '@l=split(" ", $_); if ($l[1] =~ "^[chr]*([0-9]+|M|MT|X|Y|\\*)") { print join("\t", @l)."\n"; }' | sortByChrPos.pl --file - --chrcol 2 --poscol 1 > $f.sorted.txt
(echo "chr"; cut -f 2 $f.sorted.txt) > ./chrcounts.all.txt
for f in ./4860STDY7352*/*.chr_counts; do
    cat $f | sed "s/^[ ]*//" | perl -ne '@l=split(" ", $_); if ($l[1] =~ "^[chr]*([0-9]+|M|MT|X|Y|\\*)") { print join("\t", @l)."\n"; }' | sortByChrPos.pl --file - --chrcol 2 --poscol 1 > $f.sorted.txt
    fbname=`basename "$f" | cut -d. -f1`
    cp ./chrcounts.all.txt ./chrcounts.all.tmp
    paste ./chrcounts.all.tmp <(echo $fbname; cut -f 1 $f.sorted.txt) > ./chrcounts.all.txt
done


# Count reads mapping to genes
GENCODE_GTF=$JS/reference/GRCh38/gencode.v27.annotation.gtf
cut -f1 irods.sample_lanes.psen.txt | submitJobs.py --MEM 1000 -j featureCounts -c "~/src/utils/counts/bam2counts.py --sampleDir . --gtf $GENCODE_GTF --strand 2 --countsSuffix .counts.txt --bamSuffix .bam --execute True"

# I noticed a high rate of reads not mapping to any feature (gene). The same was not
# true for my sensory neurons. I thought I would check whether I get a lower rate
# when counting using the same gene annotations as I did for sensory neurons.
GENCODE_BASIC_GTF=/lustre/scratch117/cellgen/team170/js29/sensoryneurons/annotations/Homo_sapiens.GRCh38.79.gencode_basic.gtf
cut -f1 irods.sample_lanes.txt | head -n 1 | submitJobs.py --MEM 1000 -j featureCounts -c "~/src/utils/counts/bam2counts.py --sampleDir . --gtf $GENCODE_BASIC_GTF --strand 2 --countsSuffix .counts.basic.txt --bamSuffix .bam --execute True"
# With the gencode_basic annotations, I actually get slightly worse results, so this
# does not explain the low rate of reads mapping to features.

# Look at the summary of gene read counts - I made this manually into a spreadsheet
cat ./*/*.counts.txt.summary

cd $SEQ/RNA/analysis
# Combine counts for samples into a single file
Rscript $SRC/rna/ipsneurons.combine_counts.R ../irods.sample_lanes.psen.txt .. ipsneurons.psen.counts.txt.gz

cat ../irods.sample_lanes.txt ../irods.sample_lanes.psen.txt > ../irods.sample_lanes.all.txt
Rscript $SRC/rna/ipsneurons.combine_counts.R ../irods.sample_lanes.all.txt .. ipsneurons.all.counts.txt.gz

