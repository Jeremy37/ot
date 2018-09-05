#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
PATH=$PATH:/software/solexa/pkg/bwa/0.7.17

EXP_DIR=$JS/experiment/transcribed/batch1_redo
cd $EXP_DIR

#cd $JS/reference/GRCh38
#submitJobs.py --MEM 8000 -j bwa_index -c "bwa index -a bwtsw $JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa -p $JS/reference/GRCh38/bwa_index/Homo_sapiens.GRCh38_15.fa"

# Find all fastq files from this experiment
find fastq -name *.fastq.gz > all.fastq.files.txt
find fastq -name "*.fastq.gz" -exec cp {} fastq2 \;

# Manually create a metadata file that describes the mapping from fastq filenames
# to samples.
META=transcribed_experiment.batch1_redo.meta.txt
META=transcribed_experiment.batch1_erica.meta.txt
META=transcribed_experiment.batch1_erica.meta.fixed_samples.txt
META=transcribed_experiment.batch1_erica.meta.x.txt

# Align to the full human genome
sed '1d' $META | submitJobs.py --MEM 8000 --ncores 2 -j align_fastq \
    -c "~/src/utils/align/bwaAlignPE.py --outputDir bam --noSampleDir --fastqDir fastq2 --nCores 2 --genomeDir $JS/reference/GRCh38/bwa_index/Homo_sapiens.GRCh38_15.fa"

grep "Successfully" FarmOut/align_fastq.651*.txt | wc -l
grep -i "Fail|Error" FarmOut/align_fastq.651*.txt | wc -l

# Sort bams by coordinate
PATH=/software/solexa/pkg/biobambam/2.0.79/bin/:$PATH
sed '1d' $META | cut -f1 | submitJobs.py --MEM 4000 --ncores 1 -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir bam --outdir bam --noSampleDir"
grep "Successfully" FarmOut/bamsort.660*.txt | wc -l
grep -i "Fail|Error" FarmOut/bamsort.660*.txt | wc -l

#Index bams
#sed '1d' $META | cut -f1 | submitJobs.py --MEM 1000 -j index_bams -c "python ~/src/utils/bam/index-bams.py --bamdir bam --insuffix .sortedByCoord.bam --execute True"

#Count number of reads for each chromosome
sed '1d' $META | cut -f1 | submitJobs.py --MEM 3000 -j count_chrs -c "python ~/src/utils/bam/countReadsPerChr.py --indir bam --outdir bam --insuffix .sortedByCoord.bam --noSampleDir"
grep "Successfully" FarmOut/count_chrs.660*.txt | wc -l
grep -i "Fail|Error" FarmOut/count_chrs.660*.txt | wc -l

ll bam/*.chr_counts | sed -e 's/  / /g' | sed -e 's/ /\t/g' | cut -f 9 > chr_counts.files.txt
Rscript $JS/src/experiment/mergeChrCounts.R chr_counts.files.txt > chr_counts.all.txt

ll bam/*.chr_counts | sed -e 's/  / /g' | sed -e 's/ /\t/g' | cut -f 9 > chr_counts.files.erica.txt
Rscript $JS/src/experiment/mergeChrCounts.R chr_counts.files.erica.txt > chr_counts.erica.txt


# Merge all BAMs together (for each replicate), since each is for a different targeted locus
# First get a file listing the BAM files for each replicate
sets=( cDNA_repeat1_pcr1 cDNA_repeat1_pcr2 cDNA_repeat1_pcr3
       cDNA_repeat2_pcr1 cDNA_repeat2_pcr2 cDNA_repeat2_pcr3
       cDNA_repeat3_pcr1 cDNA_repeat3_pcr2 cDNA_repeat3_pcr3
       gDNA_repeat1_pcr1 gDNA_repeat1_pcr2 gDNA_repeat1_pcr3
       gDNA_repeat2_pcr1 gDNA_repeat2_pcr2 gDNA_repeat2_pcr3
       gDNA_repeat3_pcr1 gDNA_repeat3_pcr2 gDNA_repeat3_pcr3 )
for exp in "${sets[@]}"; do
    echo $exp
    sed '1d' $META | perl -sane 'if ($F[0] =~ /$exp/) { print $F[0].".sortedByCoord.bam"."\n" }' -- -exp=$exp > bam/transcribed_experiment.batch1.bamlist.$exp.txt
done

cd $EXP_DIR/bam
for exp in "${sets[@]}"; do
    submitJobs.py --MEM 1000 -j merge.$exp -c "samtools merge -b transcribed_experiment.batch1.bamlist.$exp.txt batch1.$exp.merged.bam"
done

cd $EXP_DIR
for exp in "${sets[@]}"; do
    echo "batch1.$exp.merged" | submitJobs.py --MEM 4000 --ncores 1 -j bamsort.$exp -c "python ~/src/utils/bam/bamSortCoord.py --indir bam --outdir bam --noSampleDir"
done

# Count mapped and unmapped reads
# Unmapped reads have the flag 0x4 set
cd $EXP_DIR/bam

for f in *_gDNA*.sortedByCoord.bam; do
  unmapped=`samtools view -f 0x4 $f | wc -l`
  mappedUnique=`samtools view -F 0x104 $f | wc -l`
  mappedMulti=`samtools view -F 0x4 -f0x100 $f | wc -l`
  echo -e "$f\t$mappedUnique\t$mappedMulti\t$unmapped" >> ../gDNA.read_counts.tsv
done
for f in *_input*.sortedByCoord.bam; do
  unmapped=`samtools view -f 0x4 $f | wc -l`
  mappedUnique=`samtools view -F 0x104 $f | wc -l`
  mappedMulti=`samtools view -F 0x4 -f0x100 $f | wc -l`
  echo -e "$f\t$mappedUnique\t$mappedMulti\t$unmapped" >> ../gDNA.read_counts.tsv
done

for f in *_cDNA*.sortedByCoord.bam; do
  unmapped=`samtools view -f 0x4 $f | wc -l`
  mappedUnique=`samtools view -F 0x104 $f | wc -l`
  mappedMulti=`samtools view -F 0x4 -f0x100 $f | wc -l`
  echo -e "$f\t$mappedUnique\t$mappedMulti\t$unmapped" >> ../cDNA.read_counts.tsv
done
for f in *_CHIP*.sortedByCoord.bam; do
  unmapped=`samtools view -f 0x4 $f | wc -l`
  mappedUnique=`samtools view -F 0x104 $f | wc -l`
  mappedMulti=`samtools view -F 0x4 -f0x100 $f | wc -l`
  echo -e "$f\t$mappedUnique\t$mappedMulti\t$unmapped" >> ../cDNA.read_counts.tsv
done
# Overall things look fine - usually 95% or more of the reads are uniquely mapped
# But a few samples look suspicious - particularly in cDNA, where sometimes a lot
# of reads are unmapped.


################################################################################
# Run my analysis pipeline

JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
PATH=$PATH:/software/solexa/pkg/bwa/0.7.17
EXP_DIR=$JS/experiment/transcribed/batch1_redo
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


# make_option(c("--regions"), type="character", default=NULL, help=""),
# make_option(c("--replicates"), type="character", default=NULL, help=""),
# make_option(c("--out"), type="character", default=NULL, help=""),
# make_option(c("--minMapQ"), type="integer", default=0, help=""),
# make_option(c("--subsample"), type="numeric", default=NULL, help=""),
# make_option(c("--max_mismatch_frac"), type="numeric", default=0.05, help=""),
# make_option(c("--viewing_window"), type="integer", default=1000, help=""),
# make_option(c("--editing_window"), type="integer", default=10, help=""),
# make_option(c("--min_window_overlap"), type="integer", default=30, help=""),
# make_option(c("--exclude_multiple_deletions"), type="logical", default=F, action="store_true", help=""),
# make_option(c("--exclude_nonspanning_reads"), type="logical", default=T, action="store_true", help=""),
# make_option(c("--exclude_nonspanning_deletions"), type="logical", default=T, action="store_true", help=""),
# make_option(c("--uns_plot_min_gDNA"), type="integer", default=10, help=""),
# make_option(c("--uns_plot_min_cDNA"), type="integer", default=0, help=""),
# make_option(c("--uns_plot_max_udps"), type="integer", default=40, help=""),
# make_option(c("--no_allele_profile"), type="logical", default=F, action="store_true", help=""),
# make_option(c("--no_site_profile"), type="logical", default=F, action="store_true", help="")



################################################################################
# Run analysis for Sarah's PTK2B experiment from August 2018

JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
PATH=$PATH:/software/solexa/pkg/bwa/0.7.17

EXP_DIR=$JS/experiment/transcribed/batch1_redo
cd $EXP_DIR

cp -r fastq/PTK2B fastq_PTK2B
gzip fastq_PTK2B/*.fastq

# Manually create a metadata file that describes the mapping from fastq filenames to samples.
META=sarah_ptk2b.ptk2b.meta.txt

# Align to the full human genome
# **********MAKE sure the fastqDir is correct!
sed '1d' $META | submitJobs.py --MEM 8000 --ncores 2 -j align_fastq.sarah_ptk2b \
    -c "~/src/utils/align/bwaAlignPE.py --outputDir bam --noSampleDir --fastqDir fastq_PTK2B --nCores 2 --genomeDir $JS/reference/GRCh38/bwa_index/Homo_sapiens.GRCh38_15.fa"

grep "Successfully" FarmOut/align_fastq.sarah_ptk2b.6611*.txt | wc -l
grep -i "Fail|Error" FarmOut/align_fastq.sarah_ptk2b.6611*.txt | wc -l

# Sort bams by coordinate
PATH=/software/solexa/pkg/biobambam/2.0.79/bin/:$PATH
sed '1d' $META | cut -f1 | submitJobs.py --MEM 4000 --ncores 1 -j bamsort.sarah_ptk2b -c "python ~/src/utils/bam/bamSortCoord.py --indir bam --outdir bam --noSampleDir"
grep "Successfully" FarmOut/bamsort.sarah_ptk2b.6611*.txt | wc -l
grep -i "Fail|Error" FarmOut/bamsort.sarah_ptk2b.6611*.txt | wc -l


# Use grep to count HDR or WT sequences in fastq files
cd fastq_PTK2B
zcat 249_S249_L001_R1_001.fastq.gz | grep GAATTGCACAACA | wc -l


