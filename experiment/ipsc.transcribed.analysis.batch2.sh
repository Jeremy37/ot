#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
PATH=$PATH:/software/solexa/pkg/bwa/0.7.17

EXP_DIR=$JS/experiment/transcribed/batch2
#md $EXP_DIR
cd $EXP_DIR

#cd $JS/reference/GRCh38
#submitJobs.py --MEM 8000 -j bwa_index -c "bwa index -a bwtsw $JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa -p $JS/reference/GRCh38/bwa_index/Homo_sapiens.GRCh38_15.fa"

# Make a fasta file with amplicon reference sequences in reference/batch2_amplicons.fa
submitJobs.py --MEM 8000 -j bwa_index -q yesterday -c "bwa index -a bwtsw $EXP_DIR/reference/batch2_amplicons.fa"


# Manually create a metadata file that describes the mapping from fastq filenames
# to samples.
META=transcribed_experiment.batch2.meta.txt
META=transcribed_experiment.batch2.meta.14_15.txt
META=transcribed_experiment.batch2.meta.9.txt
META=transcribed_experiment.batch2.meta.erica.txt
META=transcribed_experiment.batch2.meta.erica.redo.txt
META=transcribed_experiment.batch2.meta.redo.txt

BAMDIR=bam


# Align to the full human genome
sed '1d' $META | submitJobs.py --MEM 8000 --ncores 2 -q yesterday -j align_fastq \
    -c "~/src/utils/align/bwaAlignPE.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq --nCores 2 --genomeDir $JS/reference/GRCh38/bwa_index/Homo_sapiens.GRCh38_15.fa"
grep "Successfully" FarmOut/align_fastq.917*.txt | wc -l
grep -iP "Fail|Error" FarmOut/align_fastq.917*.txt | wc -l


# ... or align to amplicon sequences
BAMDIR=bam_amplicon
sed '1d' $META | submitJobs.py --MEM 500 --ncores 2 -q normal -j align_fastq_amplicon \
    -c "~/src/utils/align/bwaAlignPE.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq --nCores 2 --genomeDir $EXP_DIR/reference/batch2_amplicons.fa"
grep "Successfully" FarmOut/align_fastq_amplicon.956*.txt | wc -l
grep -iP "Fail|Error" FarmOut/align_fastq_amplicon.956*.txt | wc -l

BAMDIR=bam_amplicon_2
# Change alignment parameters to allow larger deletions
sed '1d' $META | submitJobs.py --MEM 500 --ncores 2 -q normal -j align_fastq_amplicon.2 \
    -c "~/src/utils/align/bwaAlignPE.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq --nCores 2 --params '-O 6,12 -E 0 -A 1 -B 4 -L 150 -T 20 -k 19 -w 120 -d 180' --genomeDir $EXP_DIR/reference/batch2_amplicons.fa"
grep "Successfully" FarmOut/align_fastq_amplicon.2.95*.txt | wc -l
grep -iP "Fail|Error" FarmOut/align_fastq_amplicon.2.95*.txt | wc -l

BAMDIR=bam_amplicon_3
# Change alignment parameters to allow larger deletions
sed '1d' $META | submitJobs.py --MEM 500 --ncores 2 -q normal -j align_fastq_amplicon.3 \
    -c "~/src/utils/align/bwaAlignPE.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq --nCores 2 --params '-O 12,24 -E 1 -A 2 -B 8 -T 35 -k 19 -w 200 -d 300 -L 10 -U 20' --genomeDir $EXP_DIR/reference/batch2_amplicons.fa"
grep "Successfully" FarmOut/align_fastq_amplicon.3.958*.txt | wc -l
grep -iP "Fail|Error" FarmOut/align_fastq_amplicon.3.958*.txt | wc -l

BAMDIR=bam_amplicon_4
# Change alignment parameters to allow larger deletions
sed '1d' $META | submitJobs.py --MEM 500 --ncores 2 -q normal -j align_fastq_amplicon.4 \
    -c "~/src/utils/align/bwaAlignPE.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 20 -U 40' --genomeDir $EXP_DIR/reference/batch2_amplicons.fa"
grep "Successfully" FarmOut/align_fastq_amplicon.4.96*.txt | wc -l
grep -iP "Fail|Error" FarmOut/align_fastq_amplicon.2.96*.txt | wc -l



# Sort bams by coordinate
PATH=/software/solexa/pkg/biobambam/2.0.79/bin/:$PATH
sed '1d' $META | cut -f1 | submitJobs.py --MEM 1000 --ncores 1 -q normal -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir $BAMDIR --outdir $BAMDIR --noSampleDir"
grep "Successfully" FarmOut/bamsort.967*.txt | wc -l
grep -iP "Fail|Error" FarmOut/bamsort.967*.txt | wc -l

#Index bams
#sed '1d' $META | cut -f1 | submitJobs.py --MEM 1000 -j index_bams -c "python ~/src/utils/bam/index-bams.py --bamdir $BAMDIR --insuffix .sortedByCoord.bam --execute True"

#Count number of reads for each chromosome
sed '1d' $META | cut -f1 | submitJobs.py --MEM 200 -q normal -j count_chrs -c "python ~/src/utils/bam/countReadsPerChr.py --indir $BAMDIR --outdir $BAMDIR --insuffix .sortedByCoord.bam --noSampleDir"
grep "Successfully" FarmOut/count_chrs.967*.txt | wc -l
grep -iP "Fail|Error" FarmOut/count_chrs.967*.txt | wc -l
# Note that I first included separate amplicon sequences in the reference for regions 14 and 15
# but the region 14 amplicon was only 181 bp long, which was incorrect. I updated this and
# re-aligned region 14 and 15 fastq files, but I didn't re-align everything.

ll $BAMDIR/*.chr_counts | sed -e 's/  / /g' | sed -e 's/  / /g'  | sed -e 's/ /\t/g' | cut -f 9 > chr_counts.files.amplicon.txt
Rscript $JS/src/experiment/mergeChrCounts.R chr_counts.files.amplicon.txt > chr_counts.amplicon.2.all.txt


################################################################################
# Use FLASH to stich together read paires prior to aligning. This should greatly
# help the alignment when deletions occur within 20 - 30 bp of the ends of reads.
PATH=$JS/software/FLASH-1.2.11-Linux-x86_64/:$PATH
#mkdir fastq_flash

# We want to pass to Flash parameters for the expected read length and insert size
# to enable calculating the expected read overlap. For this we can merge information
# from the regions file into the metadata file.
Rscript $JS/src/experiment/getFlashInput.R batch2.regions.amplicon.tsv transcribed_experiment.batch2.meta.txt > transcribed_experiment.batch2.flash_input.tsv

sed '1d' transcribed_experiment.batch2.flash_input.fixed.tsv | submitJobs.py --MEM 1000 --ncores 1 -q normal -j merge_fastq_flash \
    -c "~/src/utils/fastq/flashStitchFastq.py --outdir fastq_flash --indir fastq"
grep "Successfully" FarmOut/merge_fastq_flash.982*.txt | wc -l
grep -iP "ERROR" FarmOut/merge_fastq_flash.982*.txt | wc -l

echo -e "File name\tOutput name\tTotal pairs\tCombined pairs\tUncombined pairs\tPercent combined\tMin overlap\tMax overlap\tMax mismatch dens\tWarning" > flash_output.summary.tsv
for f in FarmOut/merge_fastq_flash.982*.txt; do
  echo $f
  python $JS/src/experiment/getFlashOutputDetails.py --file $f | sed '1d' >> flash_output.summary.tsv
done

find FarmOut -name 'merge_fastq_flash*.txt' -print0 | xargs -r0 grep -H 'WARNING'
find FarmOut -name 'merge_fastq_flash*.txt' -print0 | xargs -r0 grep -H '5_cDNA_repeat1_pcr4'

# Remove one sample which has a corrupted fastq file
head -n 1 transcribed_experiment.batch2.flash_input.tsv > transcribed_experiment.batch2.flash_input.fixed.tsv
cat transcribed_experiment.batch2.flash_input.tsv | grep "22_gDNA_repeat1_pcr2" >> transcribed_experiment.batch2.flash_input.fixed.tsv

BAMDIR=bam_amplicon_flash
BAMDIR=bam_amplicon_flash_noclip
#mkdir $BAMDIR
sed '1d' transcribed_experiment.batch2.flash_input.tsv | cut -f 1,9 | submitJobs.py --MEM 500 --ncores 2 -q normal -j align_fastq_flash_amplicon \
    -c "~/src/utils/align/bwaMemAlign.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq_flash --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 20 -U 40' --genomeDir $EXP_DIR/reference/batch2_amplicons.fa"
grep "Successfully" FarmOut/align_fastq_flash_amplicon_noclip.986*.txt | wc -l
grep -iP "Fail|Error" FarmOut/align_fastq_flash_amplicon_noclip.986*.txt | wc -l

sed '1d' transcribed_experiment.batch2.flash_input.tsv | cut -f 1,9 | submitJobs.py --MEM 500 --ncores 2 -q normal -j align_fastq_flash_amplicon_noclip \
    -c "~/src/utils/align/bwaMemAlign.py --outputDir $BAMDIR --noSampleDir --fastqDir fastq_flash --nCores 2 --params '-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 200 -U 40' --genomeDir $EXP_DIR/reference/batch2_amplicons.fa"


# Sort bams by coordinate
PATH=/software/solexa/pkg/biobambam/2.0.79/bin/:$PATH
sed '1d' $META | cut -f1 | submitJobs.py --MEM 1000 --ncores 1 -q normal -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir $BAMDIR --outdir $BAMDIR --noSampleDir"
grep "Successfully" FarmOut/bamsort.986*.txt | wc -l
grep -iP "Fail|Error" FarmOut/bamsort.986*.txt | wc -l




################################################################################

# Merge all BAMs together (for each replicate), since each is for a different targeted locus
# First get a file listing the BAM files for each replicate
sets=( cDNA_repeat1_pcr1 cDNA_repeat1_pcr2 cDNA_repeat1_pcr3
       cDNA_repeat2_pcr1 cDNA_repeat2_pcr2 cDNA_repeat2_pcr3
       cDNA_repeat3_pcr1 cDNA_repeat3_pcr2 cDNA_repeat3_pcr3
       gDNA_repeat1_pcr1 gDNA_repeat1_pcr2 gDNA_repeat1_pcr3
       gDNA_repeat2_pcr1 gDNA_repeat2_pcr2 gDNA_repeat2_pcr3
       gDNA_repeat3_pcr1 gDNA_repeat3_pcr2 gDNA_repeat3_pcr3 )
       
sets=( gDNA_repeat1_pcr1 gDNA_repeat1_pcr2 gDNA_repeat1_pcr3
       cDNA_repeat1_pcr1 cDNA_repeat1_pcr2 cDNA_repeat1_pcr3
       cDNA_repeat1_pcr4 cDNA_repeat2_pcr5 cDNA_repeat2_pcr6
       cDNA_repeat2_pcr7 cDNA_repeat2_pcr4 )

for exp in "${sets[@]}"; do
    echo $exp
    sed '1d' $META | perl -sane 'if ($F[0] =~ /$exp/) { print $F[0].".sortedByCoord.bam"."\n" }' -- -exp=$exp > $BAMDIR/transcribed_experiment.batch2.bamlist.$exp.txt
done

cd $BAMDIR
for exp in "${sets[@]}"; do
    submitJobs.py --MEM 1000 -j merge.$exp -c "samtools merge -b transcribed_experiment.batch2.bamlist.$exp.txt batch2.$exp.merged.bam"
done

cd $EXP_DIR
for exp in "${sets[@]}"; do
    echo "batch2.$exp.merged" | submitJobs.py --MEM 4000 --ncores 1 -j bamsort.$exp -c "python ~/src/utils/bam/bamSortCoord.py --indir $BAMDIR --outdir $BAMDIR --noSampleDir"
done

# Count mapped and unmapped reads
# Unmapped reads have the flag 0x4 set
cd $BAMDIR

for f in *_gDNA*.sortedByCoord.bam; do
  unmapped=`samtools view -f 0x4 $f | wc -l`
  mappedUnique=`samtools view -F 0x104 $f | wc -l`
  mappedMulti=`samtools view -F 0x4 -f0x100 $f | wc -l`
  echo -e "$f\t$mappedUnique\t$mappedMulti\t$unmapped" >> ../gDNA.read_counts.amplicon.tsv
done
for f in *_cDNA*.sortedByCoord.bam; do
  unmapped=`samtools view -f 0x4 $f | wc -l`
  mappedUnique=`samtools view -F 0x104 $f | wc -l`
  mappedMulti=`samtools view -F 0x4 -f0x100 $f | wc -l`
  echo -e "$f\t$mappedUnique\t$mappedMulti\t$unmapped" >> ../cDNA.read_counts.amplicon.tsv
done
# Overall things look fine - usually 95% or more of the reads are uniquely mapped
# But a few samples look suspicious - particularly in cDNA, where sometimes a lot
# of reads are unmapped.


################################################################################
# Get counts of read start positions in a bam from each locus to help us
# determine the amplicon coords
echo -e "file\tchr\tpos\tcount" > chr_pos_table.gDNA.min100.tsv
for f in $BAMDIR/*_gDNA*repeat1_pcr1.sortedByCoord.bam; do
  samtools view -F 0x104 $f | cut -f 3,4 | python $JS/src/misc/tabulate_chr_pos.py --prefix $f >> chr_pos_table.gDNA.min100.tsv
done
echo -e "file\tchr\tpos\tcount" > chr_pos_table.cDNA.min5.tsv
for f in $BAMDIR/*_cDNA*repeat1_pcr1.sortedByCoord.bam; do
  samtools view -F 0x104 $f | cut -f 3,4 | python $JS/src/misc/tabulate_chr_pos.py --prefix $f --mincount 5 >> chr_pos_table.cDNA.min5.tsv
done

# in R
df = readr::read_tsv("chr_pos_table.cDNA.min5.tsv")
df$locus = as.integer(sapply(df$file, function(x) strsplit(x, "/|_")[[1]][2]))
df = df %>% dplyr::filter(count > 100) %>% dplyr::arrange(locus, -count)
# View the data.frame and see what positions of highest counts are, which are
# likely amplicon coords


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


submitJobs.py --MEM 8000 -j batch1.del.analysis.1 -q yesterday \
 -c "Rscript $JS/src/experiment/paired.deletion.analysis.R --regions regions.test.tsv --replicates replicates.test.tsv --out analysis/batch1.test \
     --viewing_window 40 --exclude_multiple_deletions F --exclude_nonspanning_reads T --exclude_nonspanning_deletions F"

submitJobs.py --MEM 8000 -j batch1.del.analysis.edit_dels -q yesterday \
 -c "Rscript $JS/src/experiment/paired.deletion.analysis.R --regions regions.test.tsv --replicates replicates.test.tsv --out analysis/batch1.edit_dels \
     --exclude_multiple_deletions F --exclude_nonspanning_reads T --exclude_nonspanning_deletions F"

submitJobs.py --MEM 8000 -j batch1.del.analysis.edit_dels.window_30 -q yesterday \
 -c "Rscript $JS/src/experiment/paired.deletion.analysis.R --regions regions.test.tsv --replicates replicates.test.tsv --out analysis/batch1.edit_dels.60bp_window \
     --viewing_window 30 --exclude_multiple_deletions F --exclude_nonspanning_reads T --exclude_nonspanning_deletions T"

submitJobs.py --MEM 8000 -j batch1.del.analysis.all_dels -q yesterday \
 -c "Rscript $JS/src/experiment/paired.deletion.analysis.R --regions regions.test.tsv --replicates replicates.test.tsv --out analysis/batch1.all_dels \
     --exclude_multiple_deletions F --exclude_nonspanning_reads T --exclude_nonspanning_deletions F"

