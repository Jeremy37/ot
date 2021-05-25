#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
PATH=$PATH:$JS/software/bwa-0.7.17
PATH=$JS/software5/biobambam2/bin/:$PATH
PATH=$JS/software/FLASH-1.2.11-Linux-x86_64/:$PATH

EXP_DIR=$JS/experiment/atac/genie_atac_hets
#md $EXP_DIR
cd $EXP_DIR

# Make a fasta file with amplicon reference sequences
REFERENCE=$JS/reference/GRCh38/bwa_index/Homo_sapiens.GRCh38_15.fa
REFERENCE2=$JS/experiment/atac/reference/genie_atac_amplicons.fa
#submitJobs.py --MEM 1000 -j bwa_index -q yesterday -c "bwa index -a bwtsw $REFERENCE2"

# Manually create a metadata file that describes the mapping from fastq filenames
# to samples.
META=genie_atac_hets.genie_replicates.tsv
EXP=genie_atac_hets

################################################################################
# Set up links to fastq files from different runs
cd $EXP_DIR/fastq
ln -s $JS/experiment/transcribed/fastq/Miseq_Walkup_154_ATAC_GenIE Miseq_Walkup_154
ln -s $JS/experiment/transcribed/fastq/Miseq_Walkup_203_ATAC_GenIE Miseq_Walkup_203
ln -s $JS/experiment/transcribed/fastq/Miseq_walkup_298_ATAC_controls Miseq_walkup_298
ln -s $JS/experiment/transcribed/fastq/Miseq_Walkup_307 Miseq_Walkup_307

cd $EXP_DIR
mkdir fastq_flash


################################################################################
# Use FLASH to stich together read paires prior to aligning. This should greatly
# help the alignment when deletions occur within 20 - 30 bp of the ends of reads.

# We want to pass to Flash parameters for the expected read length and insert size
# to enable calculating the expected read overlap. For this we can merge information
# from the regions file into the metadata file.
Rscript $JS/src/experiment/getFlashInput.R $EXP.genie_regions.tsv $META 150 40 > $EXP.flash_input.tsv

sed '1d' $EXP.flash_input.tsv | submitJobs.py --MEM 600 --ncores 1 -q normal -j merge_fastq_flash \
    -c "~/src/utils/fastq/flashStitchFastq.py --outdir fastq_flash --indir . --allow_outies"
 
grep "Successfully" FarmOut/merge_fastq_flash*.txt | wc -l
grep -iP "Fail|ERROR|Abort|TERM_" FarmOut/merge_fastq_flash*.txt | wc -l

echo -e "File name\tOutput name\tTotal pairs\tCombined pairs\tUncombined pairs\tPercent combined\tMin overlap\tMax overlap\tMax mismatch dens\tWarning" > flash_output.summary.tsv
for f in FarmOut/merge_fastq_flash*.txt; do
  echo $f
  python $JS/src/experiment/getFlashOutputDetails.py --file $f | sed '1d' >> flash_output.summary.tsv
done

################################################################################
# Align to genome sequence
BAMDIR=bam_flash
sed '1d' $EXP.flash_input.tsv | cut -f 1,9 | submitJobs.py --MEM 8000 --ncores 2 -q normal -j align_fastq_flash_amplicon \
    -c "~/src/utils/align/bwaMemAlign.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq_flash --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 20 -U 40' --genomeDir $REFERENCE"
grep "Successfully" FarmOut/align_fastq_flash_amplicon*.txt | wc -l
grep -iP "Fail|Error|Abort|TERM_" FarmOut/align_fastq_flash_amplicon*.txt | wc -l

# Sort bams by coordinate
sed '1d' $META | cut -f3 | submitJobs.py --MEM 1000 --ncores 1 -q normal -j bamsort.flash -c "python ~/src/utils/bam/bamSortCoord.py --indir $BAMDIR --outdir $BAMDIR --noSampleDir"
grep "Successfully" FarmOut/bamsort.flash*.txt | wc -l
grep -iP "Fail|Error|TERM_|No such" FarmOut/bamsort.flash*.txt | wc -l

#Count number of reads for each chromosome
sed '1d' $META | cut -f3 | submitJobs.py --MEM 200 -q normal -j count_chrs -c "python ~/src/utils/bam/countReadsPerChr.py --indir $BAMDIR --outdir $BAMDIR --insuffix .sortedByCoord.bam --noSampleDir"
grep "Successfully" FarmOut/count_chrs*.txt | wc -l
grep -iP "Fail|Error|TERM_" FarmOut/count_chrs*.txt | wc -l

ll $BAMDIR/*.chr_counts | sed -e 's/  / /g' | sed -e 's/  / /g'  | sed -e 's/ /\t/g' | cut -f 9 > chr_counts.files.amplicon.flash.txt
Rscript $JS/src/experiment/mergeChrCounts.R chr_counts.files.amplicon.flash.txt > chr_counts.amplicon.flash.all.txt




################################################################################
# Align to amplicon sequence
BAMDIR=bam_amplicon
sed '1d' $EXP.flash_input.tsv | grep "MIR8070" | cut -f 1,9 | submitJobs.py --MEM 8000 --ncores 2 -q normal -j align_fastq_flash_amplicon \
    -c "~/src/utils/align/bwaMemAlign.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq_flash --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 20 -U 40' --genomeDir $REFERENCE2"
grep "Successfully" FarmOut/align_fastq_flash_amplicon*.txt | wc -l
grep -iP "Fail|Error|Abort|TERM_" FarmOut/align_fastq_flash_amplicon*.txt | wc -l

# Sort bams by coordinate
sed '1d' $META | grep "MIR8070" | cut -f3 | submitJobs.py --MEM 1000 --ncores 1 -q normal -j bamsort.flash -c "python ~/src/utils/bam/bamSortCoord.py --indir $BAMDIR --outdir $BAMDIR --noSampleDir"
grep "Successfully" FarmOut/bamsort.flash*.txt | wc -l
grep -iP "Fail|Error|TERM_|No such" FarmOut/bamsort.flash*.txt | wc -l

#Count number of reads for each chromosome
sed '1d' $META | grep "MIR8070" | cut -f3 | submitJobs.py --MEM 200 -q normal -j count_chrs -c "python ~/src/utils/bam/countReadsPerChr.py --indir $BAMDIR --outdir $BAMDIR --insuffix .sortedByCoord.bam --noSampleDir"
grep "Successfully" FarmOut/count_chrs*.txt | wc -l
grep -iP "Fail|Error|TERM_" FarmOut/count_chrs*.txt | wc -l
