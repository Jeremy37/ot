JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
SEQ=$JS/datasets/microglia
PATH=/software/hgi/pkglocal/bwa-0.7.15/bin:$PATH

source $JS/software/ipsneurons.software.sh

###############################################################################
## ATAC-seq
mkdir $SEQ/ATAC
cd $SEQ

# Set up GRCh37 fasta index
submitJobs.py --MEM 8000 -n 2 -q yesterday --jobname bwa_index -c "bwa index -a bwtsw /lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/reference/GRCh37/hs37d5.fa"

sed '1d' glass_microglia.atac.fastq.txt | submitJobs.py --MEM 6000 -n 2 -q normal --jobname bwa_align \
    -c "~/src/utils/align/bwaAlignSE.py --nCores 2 --outdir ATAC --fastqDir fastq --fasta $JS/reference/GRCh37/hs37d5.fa --readFilesInput"
grep -iP "Failed|TERM|error" FarmOut/bwa_align.*.txt | wc -l

# Coordinate sort the BAM files
sed '1d' glass_microglia.metadata.atac.txt | cut -f 1 | submitJobs.py --MEM 2000 --jobname sortBams \
    -c "python ~/src/utils/bam/bamSortCoord.py --indir ATAC --outdir ATAC --insuffix .bam --outsuffix .sortedByCoord.bam"
grep "Successfully completed" FarmOut/sortBams.*.txt | wc -l
grep -iP "Failed|TERM|error" FarmOut/sortBams.*.txt | wc -l

sed '1d' glass_microglia.metadata.atac.txt | cut -f 1 | submitJobs.py --MEM 1000 --jobname indexBams \
    -c "python ~/src/utils/bam/indexBams.py --bamdir ATAC --insuffix .sortedByCoord.bam"

# Extract the "genome" file for bedtools, i.e. chromosome lengths, from the CRAM header
samtools view -H "$SEQ/ATAC/SRR5955079/SRR5955079.bam" | ~/src/utils/bam/genomeFileFromBamHeader.py > ATAC/GRCh37.genome.txt

cd $SEQ/ATAC
# Convert from CRAM to BigWig format for easier viewing in IGV
sed '1d' ../glass_microglia.metadata.atac.txt | cut -f 1 | submitJobs.py --MEM 3000 --jobname bamCoverageToBigWig \
    -c "~/src/utils/coverage/bam2bigwig.py --genome GRCh37.genome.txt --indir . --insuffix .sortedByCoord.bam"
grep "Successfully completed" FarmOut/bamCoverageToBigWig.*.txt | wc -l
grep "memory" FarmOut/bamCoverageToBigWig.*.txt | wc -l
grep -iP "Failed|TERM|error" FarmOut/bamCoverageToBigWig.*.txt | wc -l

##### ATAC-seq QC metrics
sed '1d' ../glass_microglia.metadata.atac.txt | cut -f 1 | submitJobs.py --MEM 3000 --jobname countReadsPerChr \
    -c "python ~/src/utils/bam/countReadsPerChr.py --indir . --outdir . --insuffix .sortedByCoord.bam"

f=./SRR5955079/SRR5955079.chr_counts
cat $f | sed "s/^[ ]*//" | perl -ne '@l=split(" ", $_); if ($l[1] =~ "^[chr]*([0-9]+|M|MT|X|Y|\\*)") { print join("\t", @l)."\n"; }' | sortByChrPos.pl --file - --chrcol 2 --poscol 1 > $f.sorted.txt
(echo "chr"; cut -f 2 $f.sorted.txt) > ./chrcounts.all.txt
for f in ./*/*.chr_counts; do
    cat $f | sed "s/^[ ]*//" | perl -ne '@l=split(" ", $_); if ($l[1] =~ "^[chr]*([0-9]+|M|MT|X|Y|\\*)") { print join("\t", @l)."\n"; }' | sortByChrPos.pl --file - --chrcol 2 --poscol 1 > $f.sorted.txt
    fbname=`basename "$f" | cut -d. -f1`
    cp ./chrcounts.all.txt ./chrcounts.all.tmp
    paste ./chrcounts.all.tmp <(echo $fbname; cut -f 1 $f.sorted.txt) > ./chrcounts.all.txt
done


# Call ATAC-seq peaks across all samples
mkdir peaks
submitJobs.py --MEM 3000 -j atacMacsCallPeak.microglia -q yesterday \
   -c "bash $JS/src/ipsneurons/atac/atac.macs_call_peaks.samples_file.sh atac.bamfiles.all.txt peaks atac_microglia 0.01"
submitJobs.py --MEM 3000 -j atacMacsCallPeak.microglia_exvivo -q yesterday \
   -c "bash $JS/src/ipsneurons/atac/atac.macs_call_peaks.samples_file.sh atac.bamfiles.exvivo.txt peaks atac_microglia_exvivo 0.01"
submitJobs.py --MEM 3000 -j atacMacsCallPeak.microglia_invitro -q yesterday \
   -c "bash $JS/src/ipsneurons/atac/atac.macs_call_peaks.samples_file.sh atac.bamfiles.invitro.txt peaks atac_microglia_invitro 0.01"

# How many peaks are there?
wc -l peaks/atac_microglia_peaks.narrowPeak # 214,253
# What amount of genome is covered by peaks?
cat peaks/atac_microglia_peaks.narrowPeak | awk '{sum += $3-$2} END {print sum}' # 60,269,966
# Mean peak size is 281 bp

wc -l peaks/atac_microglia_exvivo_peaks.narrowPeak # 193,947
cat peaks/atac_microglia_exvivo_peaks.narrowPeak | awk '{sum += $3-$2} END {print sum}' # 55,446,807, mean size 286

wc -l peaks/atac_microglia_invitro_peaks.narrowPeak # 140,413
cat peaks/atac_microglia_invitro_peaks.narrowPeak | awk '{sum += $3-$2} END {print sum}' # 30,605,057, mean size 218


# Convert from .narrowpeak file to gff3
cd peaks
Rscript ~js29/src/utils/counts/narrowPeakToGFF3.R atac_microglia_peaks.narrowPeak
sed -i -e 's/sequence_feature/exon/g' atac_microglia_peaks.gff3
sed -i -e 's/name/gene_id/g' atac_microglia_peaks.gff3


# Count reads mapping to peaks
sed '1d' ../glass_microglia.metadata.atac.txt | cut -f 1 | submitJobs.py --MEM 1000 -j featureCounts -c "~/src/utils/counts/bam2counts.py --sampleDir . --gtf peaks/atac_microglia_peaks.gff3 --strand 0 --countsSuffix .microglia_peaks.counts.txt --bamSuffix .sortedByCoord.bam --execute True --donotsort False --O True"
grep "Successfully completed" FarmOut/featureCounts*.txt | wc -l
grep -iP "Failed|TERM|error" FarmOut/featureCounts*.txt | wc -l

bedtools jaccard -a peaks/atac_microglia_exvivo_peaks.narrowPeak -b peaks/atac_microglia_invitro_peaks.narrowPeak
#intersection	union-intersection	jaccard	n_intersections
#24356822	61695042	0.394794	113237
# For comparison, Jaccard for iNeuron to neuron was 0.323 and for neuron to NPC was 0.469
bedtools jaccard -a peaks/atac_microglia_exvivo_peaks.narrowPeak -b $JS/macrophage/GRCh37/ATAC_macrophage.peaks.GRCh37.nochr.sorted.bed
#intersection	union-intersection	jaccard	n_intersections
#41578790	100015575	0.415723	152156
bedtools jaccard -a $JS/macrophage/GRCh37/ATAC_macrophage.peaks.GRCh37.nochr.sorted.bed -b peaks/atac_microglia_invitro_peaks.narrowPeak
#intersection	union-intersection	jaccard	n_intersections
#28094541	88351100	0.317987	127789
# Interestingly, ex vivo microglia have more overlap with iPSCDMs than with in vitro microglia. That could
# simply be because the iPSDMs have much better coverage than the fairly poor quality in vitro microglia.


cd $SEQ/ATAC/analysis
# Combine counts for samples into a single file
Rscript $JS/src/ipsneurons/rna/ipsneurons.combine_counts.R ../../glass_microglia.metadata.atac.noheader.txt .. microglia.atac.counts.txt.gz .microglia_peaks.counts.txt

# Run microgla.atac.qc.Rmd

