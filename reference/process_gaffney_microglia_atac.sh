JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
MG=$JS/datasets/gaffney_microglia/ATAC
PATH=/software/hgi/pkglocal/bwa-0.7.15/bin:$PATH
PATH=/software/solexa/pkg/biobambam/2.0.79/bin/:$PATH

cd $MG

#mkdir bam
cp /nfs/users/nfs_n/nk5/s117/BulkMicroglia/ATAC/Bams/26108/* bam/

mkdir 26108_8_1 26108_8_2 26108_8_3 26108_8_4 26108_8_5


mv 'bam/26108_8#1.bam' 26108_8_1/26108_8_1.bam
mv 'bam/26108_8#2.bam' 26108_8_2/26108_8_2.bam
mv 'bam/26108_8#3.bam' 26108_8_3/26108_8_3.bam
mv 'bam/26108_8#4.bam' 26108_8_4/26108_8_4.bam
mv 'bam/26108_8#5.bam' 26108_8_5/26108_8_5.bam
mv 'bam/26108_8#1.bam.bai' 26108_8_1/26108_8_1.bam.bai
mv 'bam/26108_8#2.bam.bai' 26108_8_2/26108_8_2.bam.bai
mv 'bam/26108_8#3.bam.bai' 26108_8_3/26108_8_3.bam.bai
mv 'bam/26108_8#4.bam.bai' 26108_8_4/26108_8_4.bam.bai
mv 'bam/26108_8#5.bam.bai' 26108_8_5/26108_8_5.bam.bai


sed '1d' microglia.metadata.atac.tsv | cut -f 1 | submitJobs.py --MEM 2000 --jobname sortBams \
    -c "python ~/src/utils/bam/bamSortCoord.py --indir . --outdir . --insuffix .bam --outsuffix .sortedByCoord.bam"
grep "Successfully completed" FarmOut/sortBams.*.txt | wc -l
grep -iP "Failed|TERM|error" FarmOut/sortBams.*.txt | wc -l

# sed '1d' microglia.metadata.atac.tsv | cut -f 1 | submitJobs.py --MEM 1000 --jobname indexBams \
#     -c "python ~/src/utils/bam/indexBams.py --bamdir . --insuffix .sortedByCoord.bam"

# Extract the "genome" file for bedtools, i.e. chromosome lengths, from the CRAM header
samtools view -H "26108_8_1/26108_8_1.bam" | ~/src/utils/bam/genomeFileFromBamHeader.py > GRCh38.genome.txt

# Convert from CRAM to BigWig format for easier viewing in IGV
sed '1d' microglia.metadata.atac.tsv | cut -f 1 | submitJobs.py --MEM 3000 --jobname bamCoverageToBigWig \
    -c "~/src/utils/coverage/bam2bigwig.py --genome GRCh38.genome.txt --indir . --insuffix .sortedByCoord.bam"
grep "Successfully completed" FarmOut/bamCoverageToBigWig.*.txt | wc -l
grep "memory" FarmOut/bamCoverageToBigWig.*.txt | wc -l
grep -iP "Failed|TERM|error" FarmOut/bamCoverageToBigWig.*.txt | wc -l



#mkdir peaks
submitJobs.py --MEM 3000 -j atacMacsCallPeak.microglia -q yesterday \
   -c "bash $JS/src/ipsneurons/atac/atac.macs_call_peaks.samples_file.sh atac.bamfiles.all.txt peaks atac_microglia 0.01"

cd peaks
Rscript ~js29/src/utils/counts/narrowPeakToGFF3.R atac_microglia_peaks.narrowPeak
sed -i -e 's/sequence_feature/exon/g' atac_microglia_peaks.gff3
sed -i -e 's/name/gene_id/g' atac_microglia_peaks.gff3

# How many peaks are there?
wc -l peaks/atac_microglia_peaks.narrowPeak # 72983
# What amount of genome is covered by peaks?
cat peaks/atac_microglia_peaks.narrowPeak | awk '{sum += $3-$2} END {print sum}' # 17,343,538
# Mean peak size is 237 bp

# Convert peaks to bed format
cd $MG
cut -f 1-4 peaks/atac_microglia_peaks.narrowPeak | sed 's/atac_microglia_//g' > peaks/atac_microglia_peaks.bed

sed '1d' microglia.metadata.atac.tsv | cut -f 1 | submitJobs.py --MEM 1000 -j featureCounts -c "~/src/utils/counts/bam2counts.py --sampleDir . --gtf peaks/atac_microglia_peaks.gff3 --strand 0 --countsSuffix .microglia_peaks.counts.tsv --bamSuffix .sortedByCoord.bam --execute True --donotsort False --O True"
grep "Successfully completed" FarmOut/featureCounts*.txt | wc -l
grep -iP "Failed|TERM|error" FarmOut/featureCounts*.txt | wc -l
