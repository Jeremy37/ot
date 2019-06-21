#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
EXP_DIR=$JS/experiment/RNA/CCDC6_clones_new
#md $EXP_DIR
cd $EXP_DIR

STAR_REFERENCE=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/reference/GRCh38/STAR_gencode_v27_index/
BAMDIR=$EXP_DIR/bam

#mkdir cram
python ~/src/utils/irods/irodsGetSamplesInStudy.py --studyName "OpenTargets AD PD finemap RNA" |  cut -f1 -d "." | uniq > irods.lanelets.tsv

# Map lanelet ids to sample ids
# imeta ls -d /seq/23846/23846_8#6.cram
python ~/src/utils/irods/irodsFetchMeta.py --irodsList irods.lanelets.tsv | sort -k1 > irods.sample_lanes.tsv
# Delete lines not corresponding to the latest RNA-seq experiments
cp irods.sample_lanes.tsv irods.samples.tsv
# Edit irods.samples.tsv to put in our own sample names

# Fetch lanelets in cram format from irods
# iget /seq/23846/23846_8#6.cram
cut -f1 irods.lanelets.tsv | python ~/src/utils/irods/fetch-irods.py --dir cram/ --suffix .cram

for f in cram/*.cram; do
  if [ ! -f $f.crai ]; then
    samtools index $f
  fi
done

# Merge the separate CRAM files from different sequencing lanes together for each sample 
cat irods.samples.tsv | submitJobs.py --MEM 2000 --jobname mergeBams \
    --command "~/src/utils/bam/mergeBams.py --indir cram --outdir bam --insuffix .cram --outsuffix .bam"
grep "Successfully completed" FarmOut/mergeBams.*.txt | wc -l
grep -iP "Abort|TERM|error" FarmOut/mergeBams.*.txt | wc -l


cut -f 1 irods.samples.tsv | submitJobs.py --MEM 1000 --jobname indexBams \
    --command "~/src/utils/bam/indexBams.py --bamdir bam --insuffix .bam"
grep "Successfully completed" FarmOut/indexBams.*.txt | wc -l
grep -iP "Abort|Failed|TERM|error" FarmOut/indexBams.*.txt | wc -l


# Make bigwig coverage tracks
samtools view -H "cram/29466_6#1.cram" | ~/src/utils/bam/genomeFileFromBamHeader.py > GRCh38.genome.txt

cut -f 1 irods.samples.tsv | submitJobs.py --MEM 3000 --jobname bamCoverageToBigWig \
    --command "~/src/utils/coverage/bam2bigwig.py --genome GRCh38.genome.txt --indir bam --split"
grep "Successfully completed" FarmOut/bamCoverageToBigWig.*.txt | wc -l
grep "memory" FarmOut/bamCoverageToBigWig.*.txt | wc -l
grep -iP "Failed|TERM|error" FarmOut/bamCoverageToBigWig.*.txt | wc -l

# Count reads mapping to genes
GENCODE_GTF=$JS/reference/GRCh38/gencode.v27.annotation.gtf
cut -f1 irods.samples.tsv | submitJobs.py --MEM 1000 -j featureCounts -c "~/src/utils/counts/bam2counts.py --sampleDir bam --gtf $GENCODE_GTF --strand 2 --countsSuffix .counts.tsv --bamSuffix .bam --execute True"
grep "Successfully completed" FarmOut/featureCounts.*.txt | wc -l
grep -iP "Failed|TERM|error" FarmOut/featureCounts.*.txt | wc -l

# Look at the summary of gene read counts - I made this manually into a spreadsheet
cat bam/*/*.counts.tsv.summary

# Combine counts for samples into a single file
Rscript $JS/src/misc/featureCounts.combine_counts.R irods.samples.tsv bam CCDC6_clones_new.counts.tsv.gz .counts.tsv

paste <(zcat CCDC6_clones_new.counts.tsv.gz) \
      <(zcat ../CCDC6_clones_rs1171830/CCDC6_clones_rs1171830.counts.tsv.gz | cut -f 1-2 --complement) \
      | gzip > CCDC6_clones_all.counts.tsv.gz


# Get a summary of counts across the samples
for f in bam/*/*.counts.tsv.summary; do
  filename=$(basename "$f")
  filename_no_ext="${filename%.*}"
  sample_name=`echo $filename_no_ext | sed 's/.counts.tsv//g'`
  #echo $sample_name >> analysis/CCDC6_clones.batch2.counts.summary.samples.txt
  sed '1d' $f | perl -sne 'print $sample."\t".$_' -- -sample=$sample_name >> analysis/CCDC6_clones.batch2.counts.summary.tsv
done

for f in ../CCDC6_clones_rs1171830/bam/*/*.counts.txt.summary; do
  filename=$(basename "$f")
  filename_no_ext="${filename%.*}"
  sample_name=`echo $filename_no_ext | sed 's/.counts.txt//g'`
  echo $sample_name >> analysis/CCDC6_clones.batch1.counts.summary.samples.txt
  sed '1d' $f | perl -sne 'print $sample."\t".$_' -- -sample=$sample_name >> analysis/CCDC6_clones.batch1.counts.summary.tsv
done

cat analysis/CCDC6_clones.batch1.counts.summary.samples.txt analysis/CCDC6_clones.batch2.counts.summary.samples.txt > analysis/CCDC6_clones.all.counts.summary.samples.txt

(echo -e "sample\tcategory\tcount";
 cat analysis/CCDC6_clones.batch1.counts.summary.tsv;
 cat analysis/CCDC6_clones.batch2.counts.summary.tsv) \
  > analysis/CCDC6_clones.all.counts.summary.tsv

