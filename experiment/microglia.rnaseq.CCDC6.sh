#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
EXP_DIR=$JS/experiment/RNA/microglia_ccdc6
cd $EXP_DIR

STAR_REFERENCE=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/reference/GRCh38/STAR_gencode_v27_index/

#mkdir cram
python ~/src/utils/irods/irodsGetSamplesInStudy.py --studyName "OpenTargets AD PD finemap RNA" |  cut -f1 -d "." | uniq > irods.lanelets.tsv

# Map lanelet ids to sample ids
# imeta ls -d /seq/23846/23846_8#6.cram
python ~/src/utils/irods/irodsFetchMeta.py --irodsList irods.lanelets.new.tsv | sort -k1 > irods.sample_lanes.new.tsv
# Delete lines not corresponding to the latest RNA-seq experiments
cp irods.sample_lanes.tsv irods.sample_lanes.new.tsv
cp irods.sample_lanes.new.tsv irods.samples.tsv
# Edit irods.samples.tsv to put in our own sample names

# Fetch lanelets in cram format from irods
# iget /seq/23846/23846_8#6.cram
cut -f1 irods.lanelets.new.tsv | python ~/src/utils/irods/fetch-irods.py --dir cram/ --suffix .cram

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
samtools view -H "cram/31739_5#31.cram" | ~/src/utils/bam/genomeFileFromBamHeader.py > GRCh38.genome.txt

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

# Link to previously aligned BAMs for microglia
ln -s $JS/experiment/RNA/CCDC6_clones_rs1171830/bam/WT1_KOLF2_micro bam/WT1_KOLF2_micro
ln -s $JS/experiment/RNA/CCDC6_clones_rs1171830/bam/WT2_KOLF2_micro bam/WT2_KOLF2_micro
ln -s $JS/experiment/RNA/CCDC6_clones_rs1171830/bam/WT3_KOLF2_micro bam/WT3_KOLF2_micro
ln -s $JS/experiment/RNA/CCDC6_clones_rs1171830/bam/HOM_HDR_D5_micro bam/HOM_HDR_D5_micro
ln -s $JS/experiment/RNA/CCDC6_clones_rs1171830/bam/HOM_HDR_D8_micro bam/HOM_HDR_D8_micro
ln -s $JS/experiment/RNA/CCDC6_clones_rs1171830/bam/HOM_HDR_H11_micro bam/HOM_HDR_H11_micro

sed '1d' microglia_ccdc6.meta.all.tsv | cut -f 1 > samples.counts_to_combine.tsv
# Combine counts for samples into a single file
Rscript $JS/src/misc/featureCounts.combine_counts.R samples.counts_to_combine.tsv bam microglia.all.counts.tsv.gz .counts.tsv


# Get a summary of counts across the samples
echo -e "sample\tcategory\tcount" > analysis/microglia_ccdc6.all.counts.summary.tsv
for f in bam/*/*.counts.tsv.summary; do
  filename=$(basename "$f")
  filename_no_ext="${filename%.*}"
  sample_name=`echo $filename_no_ext | sed 's/.counts.tsv//g'`
  #echo $sample_name >> analysis/CCDC6_clones.batch2.counts.summary.samples.txt
  sed '1d' $f | perl -sne 'print $sample."\t".$_' -- -sample=$sample_name >> analysis/microglia_ccdc6.all.counts.summary.tsv
done

# In R
summary.df = read_tsv("experiment/RNA/microglia_ccdc6/analysis/microglia_ccdc6.all.counts.summary.tsv") %>%
  dplyr::filter(grepl("Assigned|Ambiguity|NoFeatures|Unmapped", type)) %>%
  mutate(type = factor(as.character(type), levels=c("Unassigned_Ambiguity", "Unassigned_Unmapped", "Unassigned_NoFeatures", "Assigned"))) %>%
  group_by(sampleid) %>% mutate(assigned_count = count[type == "Assigned"])
ggplot(summary.df, aes(x=fct_reorder(sampleid, assigned_count), y=count, fill=type)) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 35, hjust = 1))
