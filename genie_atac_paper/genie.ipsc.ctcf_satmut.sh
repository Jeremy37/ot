#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
PATH=$PATH:$JS/software/bwa-0.7.17
PATH=$JS/software/biobambam2/bin/:$PATH
PATH=$JS/software/FLASH-1.2.11-Linux-x86_64/:$PATH

EXP_DIR=$JS/experiment/atac/genie_atac_ctcf_satmut_mar2021
#md $EXP_DIR
cd $EXP_DIR

# Make a fasta file with amplicon reference sequences
REFERENCE=$JS/reference/GRCh38/bwa_index/Homo_sapiens.GRCh38_15.fa
REFERENCE2=$JS/experiment/atac/reference/genie_atac_amplicons.fa
#submitJobs.py --MEM 1000 -j bwa_index -q yesterday -c "bwa index -a bwtsw $REFERENCE2"

# Manually create a metadata file that describes the mapping from fastq filenames
# to samples.
META=genie_atac_ctcf_satmut.replicates.tsv
EXP=genie_atac_ctcf_satmut

# Fastq files for this experiment come from 2 different sequencing runs.
# We want to merge together files with the same barcode.
mkdir fastq
RUN1=../../transcribed/fastq/Miseq_walkup_298_ATAC_controls
RUN2=../../transcribed/fastq/Miseq_Walkup_307

BARCODES=(4 12 20 28 36 44 52 60 68 100 108 116 124 132 140 148 156 164)
#BARCODES=(4 )
for BC in "${BARCODES[@]}"
do
   echo "$BC"
   zcat $RUN1/${BC}_*R1*.fastq.gz $RUN2/${BC}_*R1*.fastq.gz | gzip > fastq/${BC}_R1.fastq.gz
   zcat $RUN1/${BC}_*R2*.fastq.gz $RUN2/${BC}_*R2*.fastq.gz | gzip > fastq/${BC}_R2.fastq.gz
done


################################################################################
# First align paired-end reads to genome
mkdir bam
# mkdir bam_fwd
# mkdir bam_rev

#BAMDIR=bam
#sed '1d' $META | cut -f 3,8,9 | submitJobs.py --MEM 8000 --ncores 2 -q normal -j align_genome_pe \
#    -c "~/src/utils/align/bwaMemAlign.py --outputDir $BAMDIR --noSampleDir --fastqDir . --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 20 -U 40' --genomeDir $REFERENCE"
#grep "Successfully" FarmOut/align_genome_pe.*.txt | wc -l
#grep -iP "Fail|Error|Abort|TERM_" FarmOut/align_genome_pe*.txt | wc -l
#
#
## Sort bams by coordinate
#sed '1d' $META | cut -f3 | submitJobs.py --MEM 1000 --ncores 1 -q normal -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir $BAMDIR --outdir $BAMDIR --noSampleDir"
#grep "Successfully" FarmOut/bamsort.*.txt | wc -l
#grep -iP "Fail|Error|TERM_|No such" FarmOut/bamsort.*.txt | wc -l



################################################################################
# Use FLASH to stich together read paires prior to aligning. This should greatly
# help the alignment when deletions occur within 20 - 30 bp of the ends of reads.
mkdir fastq_flash

# We want to pass to Flash parameters for the expected read length and insert size
# to enable calculating the expected read overlap. For this we can merge information
# from the regions file into the metadata file.
Rscript $JS/src/experiment/getFlashInput.R $EXP.regions_indiv.tsv $EXP.replicates_align.tsv 150 40 > $EXP.flash_input.tsv

sed '1d' $EXP.flash_input.tsv | submitJobs.py --MEM 600 --ncores 1 -q normal -j merge_fastq_flash \
    -c "~/src/utils/fastq/flashStitchFastq.py --outdir fastq_flash --indir . --allow_outies"
 
grep "Successfully" FarmOut/merge_fastq_flash*.txt | wc -l
grep -iP "Fail|ERROR|Abort|TERM_" FarmOut/merge_fastq_flash*.txt | wc -l

echo -e "File name\tOutput name\tTotal pairs\tCombined pairs\tUncombined pairs\tPercent combined\tMin overlap\tMax overlap\tMax mismatch dens\tWarning" > flash_output.summary.tsv
for f in FarmOut/merge_fastq_flash*.txt; do
  echo $f
  python $JS/src/experiment/getFlashOutputDetails.py --file $f | sed '1d' >> flash_output.summary.tsv
done

# Align to genome sequence
BAMDIR=bam_flash
sed '1d' $EXP.flash_input.tsv | cut -f 1,9 | grep 100_ | submitJobs.py --MEM 8000 --ncores 2 -q normal -j align_fastq_flash_amplicon \
    -c "~/src/utils/align/bwaMemAlign.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq_flash --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 20 -U 40' --genomeDir $REFERENCE"
grep "Successfully" FarmOut/align_fastq_flash_amplicon*.txt | wc -l
grep -iP "Fail|Error|Abort|TERM_" FarmOut/align_fastq_flash_amplicon*.txt | wc -l

# Sort bams by coordinate
sed '1d' $EXP.flash_input.tsv | cut -f 1 | grep 100_ | submitJobs.py --MEM 1000 --ncores 1 -q normal -j bamsort.flash -c "python ~/src/utils/bam/bamSortCoord.py --indir $BAMDIR --outdir $BAMDIR --noSampleDir"
grep "Successfully" FarmOut/bamsort.flash*.txt | wc -l
grep -iP "Fail|Error|TERM_|No such" FarmOut/bamsort.flash*.txt | wc -l

#Count number of reads for each chromosome
sed '1d' $EXP.flash_input.tsv | cut -f 1 | submitJobs.py --MEM 200 -q normal -j count_chrs -c "python ~/src/utils/bam/countReadsPerChr.py --indir $BAMDIR --outdir $BAMDIR --insuffix .sortedByCoord.bam --noSampleDir"
grep "Successfully" FarmOut/count_chrs*.txt | wc -l
grep -iP "Fail|Error|TERM_" FarmOut/count_chrs*.txt | wc -l

ll $BAMDIR/*.chr_counts | sed -e 's/  / /g' | sed -e 's/  / /g'  | sed -e 's/ /\t/g' | cut -f 9 > chr_counts.files.amplicon.flash.txt
Rscript $JS/src/experiment/mergeChrCounts.R chr_counts.files.amplicon.flash.txt > chr_counts.amplicon.flash.all.txt



