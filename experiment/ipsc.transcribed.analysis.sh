#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
PATH=$PATH:/software/solexa/pkg/bwa/0.7.17

cd $JS/experiment/transcribed/batch1

#cd $JS/reference/GRCh38
#submitJobs.py --MEM 8000 -j bwa_index -c "bwa index -a bwtsw $JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa -p $JS/reference/GRCh38/bwa_index/Homo_sapiens.GRCh38_15.fa"

sed '1d' transcribed_experiment.batch1.meta.txt | submitJobs.py --MEM 8000 --ncores 2 -j align_fastq \
    -c "~/src/utils/align/bwaAlignPE.py --outputDir bam_relaxed --noSampleDir --fastqDir fastq --nCores 2 --params '-O 6,7 -E 1,2 -B 6 -L 20' --genomeDir $JS/reference/GRCh38/bwa_index/Homo_sapiens.GRCh38_15.fa"

sed '1d' transcribed_experiment.batch1.meta.txt | grep "26_gDNA_1" | submitJobs.py --MEM 8000 --ncores 2 -j align_fastq \
    -c "~/src/utils/align/bwaAlignPE.py --outputDir bam_relaxed --noSampleDir --fastqDir fastq --nCores 2 --params '-O 10,15 -E 1,4 -A 2 -B 10 -L 80 -d 180' --genomeDir $JS/reference/GRCh38/bwa_index/Homo_sapiens.GRCh38_15.fa"
sed '1d' transcribed_experiment.batch1.meta.txt | grep "26_gDNA_1" | submitJobs.py --MEM 8000 --ncores 2 -j align_fastq \
    -c "~/src/utils/align/bwaAlignPE.py --outputDir bam_relaxed --noSampleDir --fastqDir fastq --nCores 2 --params '-O 6,10 -E 0,2 -B 6 -L 50' --genomeDir $JS/reference/GRCh38/bwa_index/Homo_sapiens.GRCh38_15.fa"
sed '1d' transcribed_experiment.batch1.meta.txt | grep "26_gDNA_1" | submitJobs.py --MEM 8000 --ncores 2 -j align_fastq \
    -c "~/src/utils/align/bwaAlignPE.py --outputDir bam_relaxed --noSampleDir --fastqDir fastq --nCores 2 --params '-O 6,20 -E 0 -B 5 -L 150' --genomeDir $JS/reference/GRCh38/bwa_index/Homo_sapiens.GRCh38_15.fa"
sed '1d' transcribed_experiment.batch1.meta.txt | grep "26_gDNA_1" | submitJobs.py --MEM 8000 --ncores 2 -j align_fastq \
    -c "~/src/utils/align/bwaAlignPE.py --outputDir bam_relaxed --noSampleDir --fastqDir fastq --nCores 2 --params '-O 18 -E 1 -A 3 -B 12 -L 15' --genomeDir $JS/reference/GRCh38/bwa_index/Homo_sapiens.GRCh38_15.fa"
sed '1d' transcribed_experiment.batch1.meta.txt | grep "26_cDNA_1" | submitJobs.py --MEM 8000 --ncores 2 -j align_fastq \
    -c "~/src/utils/align/bwaAlignPE.py --outputDir bam_relaxed --noSampleDir --fastqDir fastq --nCores 2 --params '-O 18 -E 1 -A 3 -B 12 -L 15' --genomeDir $JS/reference/GRCh38/bwa_index/Homo_sapiens.GRCh38_15.fa"

grep "Successfully" FarmOut/align_fastq.*.txt | wc -l
grep -i "Fail|Error" FarmOut/align_fastq.*.txt | wc -l

# Using Kaur's genome index
#sed '1d' transcribed_experiment.batch1.meta.txt | submitJobs.py --MEM 8000 --ncores 2 -j align_fastq \
#    -c "~/src/utils/align/bwaAlignPE.py --outputDir bam --noSampleDir --fastqDir fastq --nCores 2 --params '' --genomeDir /lustre/scratch117/cellgen/team170/ka8/annotations/GRCh38/bwa_index/GRCh38"

#Sort bams by coordinate
PATH=/software/solexa/pkg/biobambam/2.0.79/bin/:$PATH
sed '1d' transcribed_experiment.batch1.meta.txt | cut -f1 | submitJobs.py --MEM 4000 --ncores 1 -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir bam --outdir bam --noSampleDir"
echo "26_gDNA_1" | submitJobs.py --MEM 4000 --ncores 1 -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir bam_relaxed --outdir bam_relaxed --noSampleDir"
echo "26_cDNA_1" | submitJobs.py --MEM 4000 --ncores 1 -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir bam_relaxed --outdir bam_relaxed --noSampleDir"
grep "Successfully" FarmOut/bamsort.*.txt | wc -l
grep -i "Fail|Error" FarmOut/bamsort.*.txt | wc -l

#Index bams
#sed '1d' transcribed_experiment.batch1.meta.txt | cut -f1 | submitJobs.py --MEM 1000 -j index_bams -c "python ~/src/utils/bam/index-bams.py --bamdir bam --insuffix .sortedByCoord.bam --execute True"

#Count number of reads for each chromosome
sed '1d' transcribed_experiment.batch1.meta.txt | cut -f1 | submitJobs.py --MEM 3000 -j count_chrs -c "python ~/src/utils/bam/countReadsPerChr.py --indir bam --outdir bam --insuffix .sortedByCoord.bam --noSampleDir"

ll bam/*.chr_counts | sed -e 's/  / /g' | sed -e 's/ /\t/g' | cut -f 9 > chr_counts.files.txt
Rscript $JS/src/experiment/mergeChrCounts.R chr_counts.files.txt > chr_counts.all.txt


# Merge all BAMs together (for each replicate), since each is for a different targeted locus
# First get a file listing the BAM files for each replicate
sed '1d' transcribed_experiment.batch1.meta.txt | perl -ane 'if ($F[5] =~ /cDNA/ && $F[6] == 1) { print $F[0].".sortedByCoord.bam"."\n" }' > bam/transcribed_experiment.batch1.cDNA.rep1.bamlist.txt
sed '1d' transcribed_experiment.batch1.meta.txt | perl -ane 'if ($F[5] =~ /cDNA/ && $F[6] == 2) { print $F[0].".sortedByCoord.bam"."\n" }' > bam/transcribed_experiment.batch1.cDNA.rep2.bamlist.txt
sed '1d' transcribed_experiment.batch1.meta.txt | perl -ane 'if ($F[5] =~ /cDNA/ && $F[6] == 3) { print $F[0].".sortedByCoord.bam"."\n" }' > bam/transcribed_experiment.batch1.cDNA.rep3.bamlist.txt
sed '1d' transcribed_experiment.batch1.meta.txt | perl -ane 'if ($F[5] =~ /gDNA/ && $F[6] == 1) { print $F[0].".sortedByCoord.bam"."\n" }' > bam/transcribed_experiment.batch1.gDNA.rep1.bamlist.txt
sed '1d' transcribed_experiment.batch1.meta.txt | perl -ane 'if ($F[5] =~ /gDNA/ && $F[6] == 2) { print $F[0].".sortedByCoord.bam"."\n" }' > bam/transcribed_experiment.batch1.gDNA.rep2.bamlist.txt
sed '1d' transcribed_experiment.batch1.meta.txt | perl -ane 'if ($F[5] =~ /gDNA/ && $F[6] == 3) { print $F[0].".sortedByCoord.bam"."\n" }' > bam/transcribed_experiment.batch1.gDNA.rep3.bamlist.txt

cd bam
submitJobs.py --MEM 1000 -j merge_cDNA_1 -c "samtools merge -b transcribed_experiment.batch1.cDNA.rep1.bamlist.txt batch1.cDNA.rep1.merged.bam"
submitJobs.py --MEM 3000 -j merge_cDNA_1 -c "samtools merge -b transcribed_experiment.batch1.cDNA.rep2.bamlist.txt batch1.cDNA.rep2.merged.bam"
submitJobs.py --MEM 3000 -j merge_cDNA_1 -c "samtools merge -b transcribed_experiment.batch1.cDNA.rep3.bamlist.txt batch1.cDNA.rep3.merged.bam"
submitJobs.py --MEM 3000 -j merge_cDNA_1 -c "samtools merge -b transcribed_experiment.batch1.gDNA.rep1.bamlist.txt batch1.gDNA.rep1.merged.bam"
submitJobs.py --MEM 3000 -j merge_cDNA_1 -c "samtools merge -b transcribed_experiment.batch1.gDNA.rep2.bamlist.txt batch1.gDNA.rep2.merged.bam"
submitJobs.py --MEM 3000 -j merge_cDNA_1 -c "samtools merge -b transcribed_experiment.batch1.gDNA.rep3.bamlist.txt batch1.gDNA.rep3.merged.bam"

cd $JS/experiment/transcribed/batch1
echo "batch1.cDNA.rep1.merged" | submitJobs.py --MEM 4000 --ncores 1 -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir bam --outdir bam --noSampleDir"
echo "batch1.cDNA.rep2.merged" | submitJobs.py --MEM 4000 --ncores 1 -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir bam --outdir bam --noSampleDir"
echo "batch1.cDNA.rep3.merged" | submitJobs.py --MEM 4000 --ncores 1 -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir bam --outdir bam --noSampleDir"
echo "batch1.gDNA.rep1.merged" | submitJobs.py --MEM 4000 --ncores 1 -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir bam --outdir bam --noSampleDir"
echo "batch1.gDNA.rep2.merged" | submitJobs.py --MEM 4000 --ncores 1 -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir bam --outdir bam --noSampleDir"
echo "batch1.gDNA.rep3.merged" | submitJobs.py --MEM 4000 --ncores 1 -j bamsort -c "python ~/src/utils/bam/bamSortCoord.py --indir bam --outdir bam --noSampleDir"


# Count mapped and unmapped reads
# Unmapped reads have the flag 0x4 set
cd $JS/experiment/transcribed/batch1/bam

for f in *_gDNA*.sortedByCoord.bam; do
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
# Overall things look fine - usually 95% or more of the reads are uniquely mapped
# But a few samples look suspicious - particularly in cDNA, where sometimes a lot
# of reads are unmapped.


################################################################################
# Run GenERA pipeline

# Convert BAM to SAM files, required as input
cd $JS/experiment/transcribed/batch1

for f in bam/*DNA*.sortedByCoord.bam; do
  NAME=`echo $f | perl -ne '@fparts=split(/\.|\//); print $fparts[1];'`
  echo $NAME
  samtools view $f > sam/$NAME.sam
done


cd $JS/experiment/transcribed/batch1/genera
ln -s $JS/experiment/transcribed/batch1/sam ./DATA/SAM
GenERA_1.py 

bedtools getfasta -fi $JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa -bed regions.bed

cd $JS/experiment/transcribed/batch1
python $JS/src/experiment/paired.deletion.analysis.py --input analysis/input.tsv --out analysis/batch1 --genome $JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa --minMapQ 0 
#--subsample 0.1

python $JS/src/experiment/paired.deletion.analysis.py --input analysis/input.test.tsv --out analysis/batch1_test --genome $JS/reference/GRCh38/Homo_sapiens.GRCh38_15.fa --minMapQ 0 



