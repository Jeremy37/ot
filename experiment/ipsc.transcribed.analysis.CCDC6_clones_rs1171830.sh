#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
PATH=$PATH:/software/solexa/pkg/bwa/0.7.17

EXP_DIR=$JS/experiment/transcribed/CCDC6_clones_rs1171830
#md $EXP_DIR
cd $EXP_DIR

# Make a fasta file with amplicon reference sequences
REFERENCE=$JS/experiment/transcribed/reference/batch1_batch2_amplicons.fa
REFERENCE=$JS/experiment/transcribed/reference/CCDC6_insA_amplicon.fa
#submitJobs.py --MEM 1000 -j bwa_index -q yesterday -c "bwa index -a bwtsw $REFERENCE"

# Manually create a metadata file that describes the mapping from fastq filenames
# to samples.
META=CCDC6_clones_rs1171830.meta.tsv
META=CCDC6_clones_rs1171830.meta.exon_primers.tsv

################################################################################
# Use FLASH to stich together read paires prior to aligning. This should greatly
# help the alignment when deletions occur within 20 - 30 bp of the ends of reads.
PATH=$JS/software/FLASH-1.2.11-Linux-x86_64/:$PATH
# mkdir fastq_flash

# We want to pass to Flash parameters for the expected read length and insert size
# to enable calculating the expected read overlap. For this we can merge information
# from the regions file into the metadata file.
Rscript $JS/src/experiment/getFlashInput.R CCDC6_clones_rs1171830.regions.tsv CCDC6_clones_rs1171830.meta.tsv > CCDC6_clones_rs1171830.flash_input.tsv

Rscript $JS/src/experiment/getFlashInput.R CCDC6_clones_rs1171830.regions.tsv CCDC6_clones_rs1171830.meta.exon_primers.tsv > CCDC6_clones_rs1171830.flash_input.exon_primers.tsv

CCDC6_clones_rs1171830.meta.exon_primers.tsv

sed '1d' CCDC6_clones_rs1171830.flash_input.tsv | submitJobs.py --MEM 300 --ncores 1 -q normal -j merge_fastq_flash \
    -c "~/src/utils/fastq/flashStitchFastq.py --outdir fastq_flash --indir ."
    
grep "Successfully" FarmOut/merge_fastq_flash*.txt | wc -l
grep -iP "ERROR" FarmOut/merge_fastq_flash*.txt | wc -l

sed '1d' CCDC6_clones_rs1171830.flash_input.exon_primers.tsv | submitJobs.py --MEM 300 --ncores 1 -q normal -j merge_fastq_flash \
    -c "~/src/utils/fastq/flashStitchFastq.py --outdir fastq_flash --indir ."


echo -e "File name\tOutput name\tTotal pairs\tCombined pairs\tUncombined pairs\tPercent combined\tMin overlap\tMax overlap\tMax mismatch dens\tWarning" > flash_output.summary.tsv
for f in FarmOut/merge_fastq_flash*.txt; do
  echo $f
  python $JS/src/experiment/getFlashOutputDetails.py --file $f | sed '1d' >> flash_output.summary.tsv
done

#### Manually edit CCDC6_clones_rs1171830.flash_input.tsv to remove lines for the
#### 1-bp insertion and move these to a new file, .ins.flash_input.tsv

# Align to amplicon sequence
BAMDIR=bam_amplicon
#mkdir $BAMDIR
INPUT=CCDC6_clones_rs1171830.flash_input.tsv
REFERENCE=$JS/experiment/transcribed/reference/batch1_batch2_amplicons.fa
sed '1d' $INPUT | cut -f 1,9 | submitJobs.py --MEM 200 --ncores 2 -q normal -j align_fastq_flash_amplicon \
    -c "~/src/utils/align/bwaMemAlign.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq_flash --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 200 -U 40' --genomeDir $REFERENCE"

INPUT=CCDC6_clones_rs1171830.ins.flash_input.tsv
REFERENCE_INS=$JS/experiment/transcribed/reference/CCDC6_insA_amplicon.fa
sed '1d' $INPUT | cut -f 1,9 | submitJobs.py --MEM 200 --ncores 2 -q normal -j align_fastq_flash_amplicon \
    -c "~/src/utils/align/bwaMemAlign.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq_flash --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 200 -U 40' --genomeDir $REFERENCE_INS"
grep "Successfully" FarmOut/align_fastq_flash_amplicon*.txt | wc -l
grep -iP "Fail|Error" FarmOut/align_fastq_flash_amplicon*.txt | wc -l

# Sort bams by coordinate
PATH=/software/solexa/pkg/biobambam/2.0.79/bin/:$PATH
sed '1d' $META | cut -f3 | submitJobs.py --MEM 1000 --ncores 1 -q normal -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir $BAMDIR --outdir $BAMDIR --noSampleDir"
grep "Successfully" FarmOut/bamsort*.txt | wc -l
grep -iP "Fail|Error|No such" FarmOut/bamsort*.txt | wc -l


# Align the samples with exon-exon primers to the full genome using STAR to be
# able to look at splicing
INPUT=CCDC6_clones_rs1171830.flash_input.exon_primers.tsv
BAMDIR=bam_genome
REFERENCE=$JS/reference/GRCh38/bwa_index/Homo_sapiens.GRCh38_15.fa

# sed '1d' $INPUT | cut -f 1,9 | submitJobs.py --MEM 8000 --ncores 2 -q normal -j align_fastq_flash_genome_STAR \
#     -c "~/src/utils/align/bwaMemAlign.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq_flash --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 200 -U 40' --genomeDir $REFERENCE"
# grep "Successfully" FarmOut/align_fastq_flash_genome*.txt | wc -l
# grep -iP "Fail|Error" FarmOut/align_fastq_flash_genome*.txt | wc -l

STAR_REFERENCE=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/reference/GRCh38/STAR_gencode_v27_index/
sed '1d' $INPUT | cut -f 1,9 | submitJobs.py --MEM 40000 -j align_fastq_flash_genome_STAR --ncores 6 -q yesterday \
    -c "~/src/utils/align/STAR-align.single.py --outputDir $BAMDIR --readFilesInput --genomeDir $STAR_REFERENCE --runThreadN 6 --arg1 'alignEndsType EndToEnd' --arg2 'sjdbFileChrStartEnd bam_genome/STAR.exon_samples.all.SJ.out.tab'"

# Copy output sorted BAM files back into the bam_genome dir, then...
for f in bam_genome/*.Aligned.sortedByCoord.out.bam; do
	samtools index $f
done

cd bam_genome
samtools merge rs1171830_HDR_E12_cDNA_rep3_exon_all.sorted.bam \
               rs1171830_HDR_E12_cDNA_rep3_exon_pcr7.Aligned.sortedByCoord.out.bam \
               rs1171830_HDR_E12_cDNA_rep3_exon_pcr8.Aligned.sortedByCoord.out.bam \
               rs1171830_HDR_E12_cDNA_rep3_exon_pcr9.Aligned.sortedByCoord.out.bam \
               rs1171830_HDR_E12_cDNA_rep3_exon_pcr10.Aligned.sortedByCoord.out.bam

samtools merge rs1171830_HDR_C11_cDNA_rep3_exon_all.sorted.bam \
               rs1171830_HDR_C11_cDNA_rep3_exon_pcr7.Aligned.sortedByCoord.out.bam \
               rs1171830_HDR_C11_cDNA_rep3_exon_pcr8.Aligned.sortedByCoord.out.bam \
               rs1171830_HDR_C11_cDNA_rep3_exon_pcr9.Aligned.sortedByCoord.out.bam \
               rs1171830_HDR_C11_cDNA_rep3_exon_pcr10.Aligned.sortedByCoord.out.bam

samtools merge rs1171830_HDR_B7_cDNA_rep3_exon_all.sorted.bam \
               rs1171830_HDR_B7_cDNA_rep3_exon_pcr7.Aligned.sortedByCoord.out.bam \
               rs1171830_HDR_B7_cDNA_rep3_exon_pcr8.Aligned.sortedByCoord.out.bam \
               rs1171830_HDR_B7_cDNA_rep3_exon_pcr9.Aligned.sortedByCoord.out.bam \
               rs1171830_HDR_B7_cDNA_rep3_exon_pcr10.Aligned.sortedByCoord.out.bam

samtools merge rs1171830_HDR_A10_cDNA_rep3_exon_all.sorted.bam \
               rs1171830_HDR_A10_cDNA_rep3_exon_pcr7.Aligned.sortedByCoord.out.bam \
               rs1171830_HDR_A10_cDNA_rep3_exon_pcr8.Aligned.sortedByCoord.out.bam \
               rs1171830_HDR_A10_cDNA_rep3_exon_pcr9.Aligned.sortedByCoord.out.bam \
               rs1171830_HDR_A10_cDNA_rep3_exon_pcr10.Aligned.sortedByCoord.out.bam

samtools merge rs1171830_DEL_B5_cDNA_rep3_exon_all.sorted.bam \
               rs1171830_DEL_B5_cDNA_rep3_exon_pcr7.Aligned.sortedByCoord.out.bam \
               rs1171830_DEL_B5_cDNA_rep3_exon_pcr8.Aligned.sortedByCoord.out.bam \
               rs1171830_DEL_B5_cDNA_rep3_exon_pcr9.Aligned.sortedByCoord.out.bam \
               rs1171830_DEL_B5_cDNA_rep3_exon_pcr10.Aligned.sortedByCoord.out.bam

samtools merge rs1171830_KO_C7_cDNA_rep3_exon_all.sorted.bam \
               rs1171830_KO_C7_cDNA_rep3_exon_pcr7.Aligned.sortedByCoord.out.bam \
               rs1171830_KO_C7_cDNA_rep3_exon_pcr8.Aligned.sortedByCoord.out.bam \
               rs1171830_KO_C7_cDNA_rep3_exon_pcr9.Aligned.sortedByCoord.out.bam \
               rs1171830_KO_C7_cDNA_rep3_exon_pcr10.Aligned.sortedByCoord.out.bam

samtools index rs1171830_HDR_E12_cDNA_rep3_exon_all.sorted.bam
samtools index rs1171830_HDR_C11_cDNA_rep3_exon_all.sorted.bam
samtools index rs1171830_HDR_B7_cDNA_rep3_exon_all.sorted.bam
samtools index rs1171830_HDR_A10_cDNA_rep3_exon_all.sorted.bam
samtools index rs1171830_DEL_B5_cDNA_rep3_exon_all.sorted.bam
samtools index rs1171830_KO_C7_cDNA_rep3_exon_all.sorted.bam

submitJobs.py --MEM 60000 -j star_make_index --ncores 6 -q yesterday \
    -c "$JS/src/reference/star_make_index.sh"

submitJobs.py --MEM 60000 -j align_fastq_flash_genome_STAR_pass2 --ncores 6 -q yesterday \
    -c "./run_star_align.sh"
    

#Count number of reads for each chromosome
sed '1d' $META | cut -f3 | submitJobs.py --MEM 200 -q normal -j count_chrs -c "python ~/src/utils/bam/countReadsPerChr.py --indir $BAMDIR --outdir $BAMDIR --insuffix .sortedByCoord.bam --noSampleDir"
grep "Successfully" FarmOut/count_chrs*.txt | wc -l
grep -iP "Fail|Error" FarmOut/count_chrs*.txt | wc -l

ll $BAMDIR/*.chr_counts | sed -e 's/  / /g' | sed -e 's/  / /g'  | sed -e 's/ /\t/g' | cut -f 9 > chr_counts.files.amplicon.txt
Rscript $JS/src/experiment/mergeChrCounts.R chr_counts.files.amplicon.txt > chr_counts.amplicon.all.txt

ll $BAMDIR/*.chr_counts | sed -e 's/  / /g' | sed -e 's/  / /g'  | sed -e 's/ /\t/g' | cut -f 9 > chr_counts.files.amplicon.exon_primers.txt
Rscript $JS/src/experiment/mergeChrCounts.R chr_counts.files.amplicon.exon_primers.txt > chr_counts.amplicon.exon_primers.txt


################################################################################
# Run my analysis pipeline

JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
PATH=$PATH:/software/solexa/pkg/bwa/0.7.17
EXP_DIR=$JS/experiment/transcribed/CCDC6_clones_rs1171830
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


