#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
#PATH=$JS/software4/bwa-0.7.17:$PATH
#PATH=$JS/software4/biobambam2/2.0.79/bin/:$PATH
PATH=$JS/software/FLASH-1.2.11-Linux-x86_64/:$PATH
# farm3
#PATH=/software/solexa/pkg/biobambam/2.0.79/bin/:$PATH
#PATH=$PATH:/software/solexa/pkg/bwa/0.7.17

EXP_DIR=$JS/experiment/transcribed/batch4_MUL1_CD33
#md $EXP_DIR
cd $EXP_DIR

# Make a fasta file with amplicon reference sequences
REFERENCE=$JS/experiment/transcribed/reference/genie_all_amplicons.fa
#submitJobs.py --MEM 1000 -j bwa_index -q yesterday -c "bwa index -a bwtsw $REFERENCE"

# Manually create a metadata file that describes the mapping from fastq filenames
# to samples.
META=batch4_MUL1_CD33.genie.meta.tsv
BAMDIR=bam_amplicon

# Make links to fastq files that are in subfolders
mkdir fastq
cd fastq
find $JS/experiment/transcribed/fastq/MiSeq_Walkup_237-143988857 -type f -name '*.fastq.gz' -exec ln -s {} . \;

# Metadata document here:
# https://docs.google.com/spreadsheets/d/1D2A7OkiOMhMeKWWONFMB2yvyhh7J4fittZRw6wmhjxw/edit#gid=54704208

################################################################################
# Use FLASH to stich together read paires prior to aligning. This should greatly
# help the alignment when deletions occur within 20 - 30 bp of the ends of reads.

# We want to pass to Flash parameters for the expected read length and insert size
# to enable calculating the expected read overlap. For this we can merge information
# from the regions file into the metadata file.
Rscript $JS/src/experiment/getFlashInput.R batch4_MUL1_CD33.genie.regions.tsv $META 150 > batch4_MUL1_CD33.flash_input.tsv

sed '1d' batch4_MUL1_CD33.flash_input.tsv | submitJobs.py --MEM 300 --ncores 1 -q normal -j merge_fastq_flash \
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
sed '1d' batch4_MUL1_CD33.flash_input.tsv | cut -f 1,9 | submitJobs.py --MEM 200 --ncores 2 -q normal -j align_fastq_flash_amplicon \
    -c "~/src/utils/align/bwaMemAlign.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq_flash --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 20 -U 40' --genomeDir $REFERENCE"
grep "Successfully" FarmOut/align_fastq_flash_amplicon*.txt | wc -l
grep -iP "Fail|Error|Abort|TERM_" FarmOut/align_fastq_flash_amplicon*.txt | wc -l


# Sort bams by coordinate (works on farm3 only)
sed '1d' $META | cut -f3 | submitJobs.py --MEM 1000 --ncores 1 -q normal -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir $BAMDIR --outdir $BAMDIR --noSampleDir"
grep "Successfully" FarmOut/bamsort*.txt | wc -l
grep -iP "Fail|Error|No such|TERM_" FarmOut/bamsort*.txt | wc -l


#Count number of reads for each chromosome
sed '1d' $META | cut -f3 | submitJobs.py --MEM 200 -q normal -j count_chrs -c "python ~/src/utils/bam/countReadsPerChr.py --indir $BAMDIR --outdir $BAMDIR --insuffix .sortedByCoord.bam --noSampleDir"
grep "Successfully" FarmOut/count_chrs*.txt | wc -l
grep -iP "Fail|Error|TERM_" FarmOut/count_chrs*.txt | wc -l

ll $BAMDIR/*.chr_counts | sed -e 's/  / /g' | sed -e 's/  / /g'  | sed -e 's/ /\t/g' | cut -f 9 > chr_counts.files.amplicon.txt
Rscript $JS/src/experiment/mergeChrCounts.R chr_counts.files.amplicon.txt > chr_counts.amplicon.all.txt


################################################################################
# Do the same but aligning to an amplicon with just the CD33 exons

# This is really awkward, because in my GenIE functions I don't have a way to
# check cDNA and gDNA against different reference sequences. But because gDNA
# isn't spliced, and I don't count reads with dels / insertions as WT or HDR,
# there's no good way to quantify HDR in both cDNA and gDNA at the same time.
# So I have hacked in a solution to GenIE specific for CD33, which will later
# be removed. I'll provide BAM files for cDNA that are aligned to a reference
# with just the exons, and for gDNA aligned to the full amplicon.
REFERENCE=$JS/experiment/transcribed/reference/genie_cd33_exons.fa
#submitJobs.py --MEM 1000 -j bwa_index -q yesterday -c "bwa index -a bwtsw $REFERENCE"
BAMDIR=bam_cd33_exons
META=batch4_MUL1_CD33.genie.meta.cDNA_exons.tsv

sed '1d' batch4_MUL1_CD33.flash_input.tsv | grep "cDNA" | cut -f 1,9 | submitJobs.py --MEM 200 --ncores 2 -q normal -j align_fastq_flash_cd33exons \
    -c "~/src/utils/align/bwaMemAlign.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq_flash --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 20 -U 40' --genomeDir $REFERENCE"
grep "Successfully" FarmOut/align_fastq_flash_cd33exons*.txt | wc -l
grep -iP "Fail|Error|Abort|TERM_" FarmOut/align_fastq_flash_cd33exons*.txt | wc -l


# Sort bams by coordinate (works on farm3 only)
sed '1d' batch4_MUL1_CD33.genie.meta.cDNA_exons.cDNA_only.tsv | cut -f3 | submitJobs.py --MEM 1000 --ncores 1 -q normal -j bamsort.cd33exons -c "python ~/src/utils/bam/bamSortCoord.py --indir $BAMDIR --outdir $BAMDIR --noSampleDir"
grep "Successfully" FarmOut/bamsort.cd33exons*.txt | wc -l
grep -iP "Fail|Error|No such|TERM_" FarmOut/bamsort.cd33exons*.txt | wc -l


