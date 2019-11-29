#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
PATH=$JS/software4/bwa-0.7.17:$PATH
PATH=$JS/software4/biobambam2/2.0.79/bin/:$PATH
PATH=$JS/software/FLASH-1.2.11-Linux-x86_64/:$PATH

EXP_DIR=$JS/experiment/transcribed/BOB_CLU_wt_spike2
#md $EXP_DIR
cd $EXP_DIR

# Make a fasta file with amplicon reference sequences
REFERENCE=$JS/experiment/transcribed/reference/genie_all_amplicons.fa
#submitJobs.py --MEM 1000 -j bwa_index -q yesterday -c "bwa index -a bwtsw $REFERENCE"

# Manually create a metadata file that describes the mapping from fastq filenames
# to samples.
META=BOB_CLU_wt_spike2.genie.meta.tsv

# Create links to fastq files
mkdir fastq
cd fastq
find $JS/experiment/transcribed/fastq/MiSeq_Walkup_217-140498364 -type f -name '*.fastq.gz' -exec ln -s {} . \;
cd $EXP_DIR
mkdir fastq_flash

# Metadata document here:
# https://docs.google.com/spreadsheets/d/1a5NpMXGwGzdxuYvhYNB1_MwkgvFm70Wk5dBaJUJoEeg/edit#gid=54704208

################################################################################
# Use FLASH to stich together read paires prior to aligning. This should greatly
# help the alignment when deletions occur within 20 - 30 bp of the ends of reads.

# We want to pass to Flash parameters for the expected read length and insert size
# to enable calculating the expected read overlap. For this we can merge information
# from the regions file into the metadata file.
Rscript $JS/src/experiment/getFlashInput.R BOB_CLU_wt_spike2.genie.regions.tsv $META 150 > BOB_CLU_wt_spike2.flash_input.tsv

sed '1d' BOB_CLU_wt_spike2.flash_input.tsv | submitJobs.py --MEM 300 --ncores 1 -q normal -j merge_fastq_flash \
    -c "~/src/utils/fastq/flashStitchFastq.py --outdir fastq_flash --indir . --allow_outies"
 
grep "Successfully" FarmOut/merge_fastq_flash*.txt | wc -l
grep -iP "Fail|ERROR|Abort|TERM_" FarmOut/merge_fastq_flash*.txt | wc -l

echo -e "File name\tOutput name\tTotal pairs\tCombined pairs\tUncombined pairs\tPercent combined\tMin overlap\tMax overlap\tMax mismatch dens\tWarning" > flash_output.summary.tsv
for f in FarmOut/merge_fastq_flash*.txt; do
  echo $f
  python $JS/src/experiment/getFlashOutputDetails.py --file $f | sed '1d' >> flash_output.summary.tsv
done

################################################################################
# Align to amplicon sequence
BAMDIR=bam_amplicon
sed '1d' BOB_CLU_wt_spike2.flash_input.tsv | cut -f 1,9 | submitJobs.py --MEM 200 --ncores 2 -q normal -j align_fastq_flash_amplicon \
    -c "~/src/utils/align/bwaMemAlign.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq_flash --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 20 -U 40' --genomeDir $REFERENCE"
sed '1d' BOB_CLU_wt_spike2.flash_input.tsv | grep "CLU_rs867230" | cut -f 1,9 | submitJobs.py --MEM 200 --ncores 2 -q normal -j align_fastq_flash_amplicon \
    -c "~/src/utils/align/bwaMemAlign.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq_flash --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 20 -U 40' --genomeDir $REFERENCE"
grep "Successfully" FarmOut/align_fastq_flash_amplicon*.txt | wc -l
grep -iP "Fail|Error|Abort|TERM_" FarmOut/align_fastq_flash_amplicon*.txt | wc -l


# Sort bams by coordinate
sed '1d' $META | grep "CLU_rs867230" | cut -f3 | submitJobs.py --MEM 1000 --ncores 1 -q normal -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir $BAMDIR --outdir $BAMDIR --noSampleDir"
grep "Successfully" FarmOut/bamsort*.txt | wc -l
grep -iP "Fail|Error|No such|TERM_" FarmOut/bamsort*.txt | wc -l


#Count number of reads for each chromosome
sed '1d' $META | grep "CLU_rs867230" | cut -f3 | submitJobs.py --MEM 200 -q normal -j count_chrs -c "python ~/src/utils/bam/countReadsPerChr.py --indir $BAMDIR --outdir $BAMDIR --insuffix .sortedByCoord.bam --noSampleDir"
grep "Successfully" FarmOut/count_chrs*.txt | wc -l
grep -iP "Fail|Error|Abort|TERM_" FarmOut/count_chrs*.txt | wc -l

ll $BAMDIR/*.chr_counts | sed -e 's/  / /g' | sed -e 's/  / /g'  | sed -e 's/ /\t/g' | cut -f 9 > chr_counts.files.amplicon.txt
Rscript $JS/src/experiment/mergeChrCounts.R chr_counts.files.amplicon.txt > chr_counts.amplicon.all.txt

