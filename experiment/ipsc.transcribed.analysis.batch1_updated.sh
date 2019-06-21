#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
PATH=$PATH:/software/solexa/pkg/bwa/0.7.17

EXP_DIR=$JS/experiment/transcribed/batch1_updated
#md $EXP_DIR
cd $EXP_DIR

# Make a fasta file with amplicon reference sequences
REFERENCE=$JS/experiment/transcribed/batch1_updated/reference/batch1_updated_amplicons.fa
#submitJobs.py --MEM 1000 -j bwa_index -q yesterday -c "bwa index -a bwtsw $REFERENCE"

# Manually create a metadata file that describes the mapping from fastq filenames
# to samples.
META=batch1_updated.meta.tsv


################################################################################
# Use FLASH to stich together read paires prior to aligning. This should greatly
# help the alignment when deletions occur within 20 - 30 bp of the ends of reads.
PATH=$JS/software/FLASH-1.2.11-Linux-x86_64/:$PATH
#mkdir fastq_flash

# We want to pass to Flash parameters for the expected read length and insert size
# to enable calculating the expected read overlap. For this we can merge information
# from the regions file into the metadata file.
Rscript $JS/src/experiment/getFlashInput.R batch1_updated.regions.tsv batch1_updated.meta.tsv > batch1_updated.flash_input.tsv

sed '1d' batch1_updated.flash_input.tsv | submitJobs.py --MEM 1000 --ncores 1 -q normal -j merge_fastq_flash \
    -c "~/src/utils/fastq/flashStitchFastq.py --outdir fastq_flash --indir ."
grep "Successfully" FarmOut/merge_fastq_flash.271*.txt | wc -l
grep -iP "ERROR" FarmOut/merge_fastq_flash.271*.txt | wc -l

echo -e "File name\tOutput name\tTotal pairs\tCombined pairs\tUncombined pairs\tPercent combined\tMin overlap\tMax overlap\tMax mismatch dens\tWarning" > flash_output.summary.tsv
for f in FarmOut/merge_fastq_flash*.txt; do
  echo $f
  python $JS/src/experiment/getFlashOutputDetails.py --file $f | sed '1d' >> flash_output.summary.tsv
done

# Align to amplicon sequence
BAMDIR=bam_amplicon
sed '1d' batch1_updated.flash_input.tsv | cut -f 1,9 | submitJobs.py --MEM 500 --ncores 2 -q normal -j align_fastq_flash_amplicon \
    -c "~/src/utils/align/bwaMemAlign.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq_flash --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 200 -U 40' --genomeDir $REFERENCE"
grep "Successfully" FarmOut/align_fastq_flash_amplicon*.txt | wc -l
grep -iP "Fail|Error" FarmOut/align_fastq_flash_amplicon*.txt | wc -l
    

# Sort bams by coordinate
PATH=/software/solexa/pkg/biobambam/2.0.79/bin/:$PATH
sed '1d' $META | cut -f1 | submitJobs.py --MEM 1000 --ncores 1 -q normal -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir $BAMDIR --outdir $BAMDIR --noSampleDir"
grep "Successfully" FarmOut/bamsort*.txt | wc -l
grep -iP "Fail|Error" FarmOut/bamsort*.txt | wc -l


#Count number of reads for each chromosome
sed '1d' $META | cut -f1 | submitJobs.py --MEM 200 -q normal -j count_chrs -c "python ~/src/utils/bam/countReadsPerChr.py --indir $BAMDIR --outdir $BAMDIR --insuffix .sortedByCoord.bam --noSampleDir"
grep "Successfully" FarmOut/count_chrs*.txt | wc -l
grep -iP "Fail|Error" FarmOut/count_chrs*.txt | wc -l
# Note that I first included separate amplicon sequences in the reference for regions 14 and 15
# but the region 14 amplicon was only 181 bp long, which was incorrect. I updated this and
# re-aligned region 14 and 15 fastq files, but I didn't re-align everything.

ll $BAMDIR/*.chr_counts | sed -e 's/  / /g' | sed -e 's/  / /g'  | sed -e 's/ /\t/g' | cut -f 9 > chr_counts.files.amplicon.txt
Rscript $JS/src/experiment/mergeChrCounts.R chr_counts.files.amplicon.txt > chr_counts.amplicon.all.txt


# Use more lenient soft-clipping parameters for aligning CLU_18. Because it is a short
# amplicon, when there are deletions then the ends of reads go into no-mans-land, and
# the resulting alignments have random insertions, unless they are allowed to be soft-clipped.
sed '1d' batch1_updated.flash_input.tsv | cut -f 1,9 | grep 18_ | submitJobs.py --MEM 500 --ncores 2 -q normal -j align_fastq_flash_amplicon \
    -c "~/src/utils/align/bwaMemAlign.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq_flash --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 20 -U 40' --genomeDir $REFERENCE"
grep "Successfully" FarmOut/align_fastq_flash_amplicon.283*.txt | wc -l
grep -iP "Fail|Error" FarmOut/align_fastq_flash_amplicon.283*.txt | wc -l

PATH=/software/solexa/pkg/biobambam/2.0.79/bin/:$PATH
sed '1d' $META | cut -f1 | grep 18_ | submitJobs.py --MEM 1000 --ncores 1 -q normal -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir $BAMDIR --outdir $BAMDIR --noSampleDir"
grep "Successfully" FarmOut/bamsort.284*.txt | wc -l
grep -iP "Fail|Error" FarmOut/bamsort.284*.txt | wc -l


################################################################################
# Run my analysis pipeline

JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
PATH=$PATH:/software/solexa/pkg/bwa/0.7.17
EXP_DIR=$JS/experiment/transcribed/batch2
cd $EXP_DIR

# For some reason, on the farm my code was crashing deep within dlyr. Installing
# the same (slightly older) version as I had on my Mac seemed to fix this.
# I tried 0.7.4, but on the farm it gives slightly different results! I.e. plots
# aren't properly sorted by UDP start position. I'm not sure what side-effect
# behaviour I'm depending on that would cause this...
# Trying 0.7.5 gave the crash again, but it's intermittent.
require(devtools)
install_version("dplyr", version = "0.7.5", repos = "http://cran.us.r-project.org")

Rscript $JS/src/experiment/paired.deletion.analysis.R --regions regions.test.tsv --replicates replicates.test.tsv --out analysis/batch1_redo.test \
  --viewing_window 40 --exclude_multiple_deletions F --exclude_nonspanning_reads T --exclude_nonspanning_deletions F --editing_window 1

