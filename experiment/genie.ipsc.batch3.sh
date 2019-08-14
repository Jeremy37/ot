#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
PATH=$JS/software4/bwa-0.7.17:$PATH
PATH=$JS/software4/biobambam2/2.0.79/bin/:$PATH
PATH=$JS/software/FLASH-1.2.11-Linux-x86_64/:$PATH

EXP_DIR=$JS/experiment/transcribed/batch3
#md $EXP_DIR
cd $EXP_DIR

# Make a fasta file with amplicon reference sequences
REFERENCE=$JS/experiment/transcribed/reference/batch3_amplicons.fa
#submitJobs.py --MEM 1000 -j bwa_index -q yesterday -c "bwa index -a bwtsw $REFERENCE"

# Manually create a metadata file that describes the mapping from fastq filenames
# to samples.
META=batch3.genie.meta.tsv

# Copy fastq files from Miseq folders all into a single folder
md fastq
for f in ../../fastq/MiSeq_Walkup_205-138408271/FASTQ_Generation_2019-07-12_08_57_52Z-189578145/*/*.fastq.gz; do
  filename=$(basename "$f")
  echo $filename
  ln -s $f $filename
done

cd $EXP_DIR
mkdir fastq_flash


################################################################################
# Use FLASH to stich together read paires prior to aligning. This should greatly
# help the alignment when deletions occur within 20 - 30 bp of the ends of reads.

# We want to pass to Flash parameters for the expected read length and insert size
# to enable calculating the expected read overlap. For this we can merge information
# from the regions file into the metadata file.
Rscript $JS/src/experiment/getFlashInput.R batch3.genie.regions.tsv $META 150 > batch3.flash_input.tsv

sed '1d' batch3.flash_input.tsv | submitJobs.py --MEM 300 --ncores 1 -q normal -j merge_fastq_flash \
    -c "~/src/utils/fastq/flashStitchFastq.py --outdir fastq_flash_no_outies --indir fastq"

sed '1d' batch3.flash_input.tsv | submitJobs.py --MEM 300 --ncores 1 -q normal -j merge_fastq_flash \
    -c "~/src/utils/fastq/flashStitchFastq.py --outdir fastq_flash --indir fastq --allow_outies"

grep "Successfully" FarmOut/merge_fastq_flash.161*.txt | wc -l
grep -iP "Fail|ERROR|Abort|TERM_" FarmOut/merge_fastq_flash.161*.txt | wc -l

echo -e "File name\tOutput name\tTotal pairs\tCombined pairs\tUncombined pairs\tPercent combined\tMin overlap\tMax overlap\tMax mismatch dens\tWarning" > flash_output.summary.tsv
for f in FarmOut/merge_fastq_flash*.txt; do
  echo $f
  python $JS/src/experiment/getFlashOutputDetails.py --file $f | sed '1d' >> flash_output.summary.tsv
done

################################################################################
# Align to amplicon sequence
BAMDIR=bam_amplicon
sed '1d' batch3.flash_input.tsv | cut -f 1,9 | submitJobs.py --MEM 200 --ncores 2 -q normal -j align_fastq_flash_amplicon \
    -c "~/src/utils/align/bwaMemAlign.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq_flash --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 20 -U 40' --genomeDir $REFERENCE"
grep "Successfully" FarmOut/align_fastq_flash_amplicon.161*.txt | wc -l
grep -iP "Fail|Error|Abort|TERM_" FarmOut/align_fastq_flash_amplicon.161*.txt | wc -l


# Sort bams by coordinate
sed '1d' $META | cut -f3 | sort | uniq | submitJobs.py --MEM 1000 --ncores 1 -q normal -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir $BAMDIR --outdir $BAMDIR --noSampleDir"
grep "Successfully" FarmOut/bamsort.161*.txt | wc -l
grep -iP "Fail|Error|No such|TERM_" FarmOut/bamsort.161*.txt | wc -l


#Count number of reads for each chromosome
sed '1d' $META | cut -f3 | sort | uniq | submitJobs.py --MEM 200 -q normal -j count_chrs -c "python ~/src/utils/bam/countReadsPerChr.py --indir $BAMDIR --outdir $BAMDIR --insuffix .sortedByCoord.bam --noSampleDir"
grep "Successfully" FarmOut/count_chrs.161*.txt | wc -l
grep -iP "Fail|Error|TERM_" FarmOut/count_chrs.161*.txt | wc -l

ll $BAMDIR/*.chr_counts | sed -e 's/  / /g' | sed -e 's/  / /g'  | sed -e 's/ /\t/g' | cut -f 9 > chr_counts.files.amplicon.txt
Rscript $JS/src/experiment/mergeChrCounts.R chr_counts.files.amplicon.txt > chr_counts.amplicon.all.txt


################################################################################
# Align the samples with exon-exon primers to the full genome using STAR to be
# able to look at splicing
INPUT=batch3.flash_input.exon_primers.tsv
BAMDIR=bam_genome

Rscript $JS/src/experiment/getFlashInput.R batch3.genie.regions.exon_primers.tsv batch3.genie.meta.exon_primers.tsv > batch3.flash_input.exon_primers.tsv

sed '1d' batch3.flash_input.exon_primers.tsv | submitJobs.py --MEM 300 --ncores 1 -q normal -j merge_fastq_flash \
    -c "~/src/utils/fastq/flashStitchFastq.py --outdir fastq_flash_exon --indir fastq"

grep "Successfully" FarmOut/merge_fastq_flash.109*.txt | wc -l
grep -iP "Fail|ERROR|Abort|TERM_" FarmOut/merge_fastq_flash.109*.txt | wc -l

echo -e "File name\tOutput name\tTotal pairs\tCombined pairs\tUncombined pairs\tPercent combined\tMin overlap\tMax overlap\tMax mismatch dens\tWarning" > flash_output.exonic.summary.tsv
for f in FarmOut/merge_fastq_flash.109*.txt; do
  echo $f
  python $JS/src/experiment/getFlashOutputDetails.py --file $f | sed '1d' >> flash_output.exonic.summary.tsv
done


STAR_REFERENCE=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/reference/GRCh38/STAR_gencode_v27_index/
sed '1d' batch3.flash_input.exon_primers.tsv | cut -f 1,9 | submitJobs.py --MEM 40000 -j align_fastq_flash_genome_STAR --ncores 6 -q yesterday \
    -c "~/src/utils/align/STAR-align.single.py --outputDir $BAMDIR --readFilesInput --genomeDir $STAR_REFERENCE --runThreadN 6 --arg1 'alignEndsType EndToEnd' "
# --arg2 'sjdbFileChrStartEnd bam_genome/STAR.exon_samples.all.SJ.out.tab'"
grep "Successfully" FarmOut/align_fastq_flash_genome_STAR*.txt | wc -l
grep -iP "Fail|ERROR|Abort|TERM_" FarmOut/align_fastq_flash_genome_STAR*.txt | wc -l

# Samtools doesn't seem to work with symlinks, so unfortunately we just copy files
for f in bam_genome/*/*.sortedByCoord.out.bam; do
  filename=$(basename "$f")
  echo $filename
  cp $f bam_genome/$filename
done

# Copy output sorted BAM files back into the bam_genome dir, then...
for f in bam_genome/*.Aligned.sortedByCoord.out.bam; do
	samtools index $f
done

