#newgrp otcoregen # Set primary group to otcoregen so new files have this group by default
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
SEQ=$JS/ipsneurons/GRCh38
SRC=$JS/src/ipsneurons
source $JS/software/ipsneurons.software.sh
kinit # login to iRODS

###############################################################################
## ATAC-seq
cd $SEQ/ATAC

# Fetch list of lanelets from IRODS
# imeta qu -z seq -d study = 'OpenTargets AD PD finemap ATAC' and target = 1 and manual_qc = 1
python ~/src/utils/irods/irodsGetSamplesInStudy.py --studyName "OpenTargets AD PD finemap ATAC" |  cut -f1 -d "." | uniq > irods.lanelets.txt

# Map lanelet ids to sample ids
# imeta ls -d /seq/23846/23846_8#6.cram
python ~/src/utils/irods/irodsFetchMeta.py --irodsList irods.lanelets.txt | sort -k1 > irods.sample_lanes.txt 

# Fetch lanelets in cram format from irods
# iget /seq/23846/23846_8#6.cram
cut -f1 irods.lanelets.txt | python ~/src/utils/irods/fetch-irods.py --dir cram/ --suffix .cram

for f in cram/*.cram; do
  if [ ! -f $f.crai ]; then
    samtools index $f
  fi
done

#4859STDY7028452	23520_7#1;23520_8#1
#samtools merge cram/4859STDY7028452.atac.cram cram/23520_7#1.cram cram/23520_8#1.cram
cat irods.sample_lanes.txt | submitJobs.py --MEM 2000 --jobname mergeBams \
    --command "~/src/utils/bam/mergeBams.py --indir cram --outdir bam --insuffix .cram --outsuffix .bam"

grep "Successfully completed" FarmOut/mergeBams.*.txt | wc -l
grep -iP "Failed|TERM|error" FarmOut/mergeBams.*.txt | wc -l

cat irods.sample_lanes.txt | tail -n +8 | head -n 1 | submitJobs.py --MEM 2000 --jobname mergeBams \
    --command "~/src/utils/bam/mergeBams.py --indir cram --outdir bam --insuffix .cram --outsuffix .bam"

cut -f 1 irods.sample_lanes.txt | tail -n +2 | submitJobs.py --MEM 1000 --jobname indexBams \
    --command "~/src/utils/bam/indexBams.py --bamdir bam --insuffix .bam"

# Extract the "genome" file for bedtools, i.e. chromosome lengths, from the CRAM header
samtools view -H "cram/23520_7#1.cram" | ~/src/utils/bam/genomeFileFromBamHeader.py > GRCh38.genome.txt

# Convert from CRAM to BigWig format for easier viewing in IGV

# bamCoverage seems like a good tool, and the option to extend reads to fragments
# is nice. In practice it takes up so much memory (>30 Gb) that I can't get it to
# run for a normal sized BAM, and it seems to use 500 bp bins by default, which is
# way too large. So I favour bedtools genomecov followed by bedGraphToBigWig.
#bamCoverage --extendReads --region 10:456700:891000 -b bam/4859STDY7028449/4859STDY7028449.bam -o 4859STDY7028449.test.bw
#bedtools genomecov -ibam bam/4859STDY7028449/4859STDY7028449.bam -bg > bam/4859STDY7028449/4859STDY7028449.bg
#bedGraphToBigWig bam/4859STDY7028449/4859STDY7028449.bg GRCh38.genome.txt bam/4859STDY7028449/4859STDY7028449.bw

cut -f 1 irods.sample_lanes.txt | submitJobs.py --MEM 3000 --jobname bamCoverageToBigWig \
    --command "~/src/utils/coverage/bam2bigwig.py --genome GRCh38.genome.txt --indir bam"
grep "Successfully completed" FarmOut/bamCoverageToBigWig.*.txt | wc -l
grep "memory" FarmOut/bamCoverageToBigWig.*.txt | wc -l
grep -iP "Failed|TERM|error" FarmOut/bamCoverageToBigWig.*.txt | wc -l

##### ATAC-seq QC metrics
cut -f 1 irods.sample_lanes.txt | submitJobs.py --MEM 3000 --jobname countReadsPerChr --command "python ~/src/utils/bam/countReadsPerChr.py --indir . --outdir . --insuffix .bam"

# Use verifyBamID to ensure that the files match the expected HIPSCI sample
KOLF2_VCF=$SEQ/genotypes/kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.vcf.gz
KOLF2_VCF=$SEQ/genotypes/kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.exons.vcf.gz
cut -f 1 irods.sample_lanes.txt | head -n 1 | submitJobs.py --MEM 1000 --jobname verifyBamID \
    -c  "python ~/src/utils/bam/runVerifyBamID.py --bamdir . --insuffix .bam --vcf $KOLF2_VCF" 
cut -f 1 irods.sample_lanes.txt | submitJobs.py --MEM 1000 --jobname verifyBamID \
    -c  "python ~/src/utils/bam/runVerifyBamID.py --bamdir . --insuffix .bam --vcf $KOLF2_VCF" 
grep "Successfully completed" FarmOut/verifyBamID*.txt | wc -l
grep -iP "Failed|TERM|error" FarmOut/verifyBamID*.txt | wc -l

(echo -ne "Sample\t"; head -n 1 ./4859STDY7028449/4859STDY7028449.verifyBamID.bestSM) > ATAC_samples.verifyBamID.bestSM
for f in ./*/*.verifyBamID.bestSM; do
    (echo -ne "$f\t"; sed '1d' $f) >> ATAC_samples.verifyBamID.bestSM
done

# Count reads mapping to each chromosome
cut -f1 irods.sample_lanes.txt | submitJobs.py --MEM 3000 -n 2 -j countReadsPerChr -c "python ~js29/src/utils/bam/countReadsPerChr.py --indir . --outdir . --insuffix .bam"
# Make a file with chromosomes 1...X,Y,MT named chrcounts.all.txt
f=./4859STDY7028449/4859STDY7028449.chr_counts
cat $f | sed "s/^[ ]*//" | perl -ne '@l=split(" ", $_); if ($l[1] =~ "^[chr]*([0-9]+|M|MT|X|Y|\\*)") { print join("\t", @l)."\n"; }' | sortByChrPos.pl --file - --chrcol 2 --poscol 1 > $f.sorted.txt
(echo "chr"; cut -f 2 $f.sorted.txt) > ./chrcounts.all.txt
for f in ./*/*.chr_counts; do
    cat $f | sed "s/^[ ]*//" | perl -ne '@l=split(" ", $_); if ($l[1] =~ "^[chr]*([0-9]+|M|MT|X|Y|\\*)") { print join("\t", @l)."\n"; }' | sortByChrPos.pl --file - --chrcol 2 --poscol 1 > $f.sorted.txt
    fbname=`basename "$f" | cut -d. -f1`
    cp ./chrcounts.all.txt ./chrcounts.all.tmp
    paste ./chrcounts.all.tmp <(echo $fbname; cut -f 1 $f.sorted.txt) > ./chrcounts.all.txt
done


# Call ATAC-seq peaks across all samples
submitJobs.py --MEM 8000 -j atacMacsCallPeak -q yesterday \
   -c "bash $SRC/atac/atac.macs_call_peaks.sh"

# How many peaks are there?
wc -l atac_multisample_peaks.narrowPeak
#406704
# What amount of genome is covered by peaks?
cat atac_multisample_peaks.narrowPeak | awk '{sum += $3-$2} END {print sum}'
# 112,225,861
# Mean peak size is 276 bp

# Convert from .narrowpeak file to gff3
Rscript ~js29/src/utils/counts/narrowPeakToGFF3.R atac_multisample_peaks.narrowPeak
sed -i -e 's/sequence_feature/exon/g' atac_multisample_peaks.gff3
sed -i -e 's/name/gene_id/g' atac_multisample_peaks.gff3


# Count reads mapping to peaks
cut -f1 irods.sample_lanes.txt | submitJobs.py --MEM 1000 -j featureCounts -c "~/src/utils/counts/bam2counts.py --sampleDir . --gtf atac_multisample_peaks.gff3 --strand 0 --countsSuffix .multisample_peaks.counts.txt --bamSuffix .bam --execute True --donotsort False --O True"

grep "Successfully completed" FarmOut/featureCounts*.txt | wc -l
grep -iP "Failed|TERM|error" FarmOut/featureCounts*.txt | wc -l

# Count reads mapping to peaks on main chromosomes only
cat atac_multisample_peaks.gff3 | perl -ane 'if ($F[0] =~ /^chr[\d|X|Y]+$/) { print }' > atac_multisample_peaks.mainchrs.gff3
wc -l atac_multisample_peaks.mainchrs.gff3
cut -f1 atac_multisample_peaks.mainchrs.gff3 | sort | uniq

cut -f1 irods.sample_lanes.txt | submitJobs.py --MEM 1000 -j featureCounts.mainchrs -c "~/src/utils/counts/bam2counts.py --sampleDir . --gtf atac_multisample_peaks.mainchrs.gff3 --strand 0 --countsSuffix .multisample_peaks.mainchrs.counts.txt --bamSuffix .bam --execute True --donotsort False --O True"

echo "4859STDY7028451" | submitJobs.py --MEM 1000 -j featureCounts.mainchrs -c "~/src/utils/counts/bam2counts.py --sampleDir . --gtf atac_multisample_peaks.mainchrs.gff3 --strand 0 --countsSuffix .multisample_peaks.mainchrs.counts.txt --bamSuffix .bam --execute True --donotsort False --O True"
echo "4859STDY7079823" | submitJobs.py --MEM 1000 -j featureCounts.mainchrs -c "~/src/utils/counts/bam2counts.py --sampleDir . --gtf atac_multisample_peaks.mainchrs.gff3 --strand 0 --countsSuffix .multisample_peaks.mainchrs.counts.txt --bamSuffix .bam --execute True --donotsort False --O True"

grep "Successfully completed" FarmOut/featureCounts.mainchrs*.txt | wc -l
grep -iP "Failed|TERM|error" FarmOut/featureCounts.mainchrs*.txt | wc -l


cd $SEQ/ATAC/analysis
# Combine counts for samples into a single file
Rscript $SRC/rna/ipsneurons.combine_counts.R ../irods.sample_lanes.txt .. ipsneurons.atac.counts.txt.gz .multisample_peaks.counts.txt

Rscript $SRC/rna/ipsneurons.combine_counts.R ../irods.sample_lanes.txt .. ipsneurons.atac.counts.mainchrs.txt.gz .multisample_peaks.mainchrs.counts.txt


submitJobs.py --MEM 6000 -j atacMacsCallPeak.ineuron -q yesterday \
   -c "bash $SRC/atac/atac.macs_call_peaks.samples_file.sh meta/ineuron.bamfiles.txt peaks atac_ineuron 0.01"
submitJobs.py --MEM 6000 -j atacMacsCallPeak.ipsc -q yesterday \
   -c "bash $SRC/atac/atac.macs_call_peaks.samples_file.sh meta/ipsc.bamfiles.txt peaks atac_ipsc 0.01"
submitJobs.py --MEM 6000 -j atacMacsCallPeak.npc -q yesterday \
   -c "bash $SRC/atac/atac.macs_call_peaks.samples_file.sh meta/npc.bamfiles.txt peaks atac_npc 0.01"
submitJobs.py --MEM 6000 -j atacMacsCallPeak.neuron -q yesterday \
   -c "bash $SRC/atac/atac.macs_call_peaks.samples_file.sh meta/neuron.bamfiles.txt peaks atac_neuron 0.01"

Rscript ~js29/src/utils/counts/narrowPeakToGFF3.R peaks/atac_ipsc_peaks.narrowPeak
Rscript ~js29/src/utils/counts/narrowPeakToGFF3.R peaks/atac_npc_peaks.narrowPeak
Rscript ~js29/src/utils/counts/narrowPeakToGFF3.R peaks/atac_neuron_peaks.narrowPeak
Rscript ~js29/src/utils/counts/narrowPeakToGFF3.R peaks/atac_ineuron_peaks.narrowPeak

cat atac_ipsc_peaks.gff3 | perl -ane 'if ($F[0] =~ /^chr[\d|X|Y]+$/) { print }' > atac_ipsc_peaks.mainchrs.gff3
cat atac_npc_peaks.gff3 | perl -ane 'if ($F[0] =~ /^chr[\d|X|Y]+$/) { print }' > atac_npc_peaks.mainchrs.gff3
cat atac_neuron_peaks.gff3 | perl -ane 'if ($F[0] =~ /^chr[\d|X|Y]+$/) { print }' > atac_neuron_peaks.mainchrs.gff3
cat atac_ineuron_peaks.gff3 | perl -ane 'if ($F[0] =~ /^chr[\d|X|Y]+$/) { print }' > atac_ineuron_peaks.mainchrs.gff3

bedtools jaccard -a peaks/atac_ipsc_peaks.narrowPeak -b peaks/atac_npc_peaks.narrowPeak
bedtools jaccard -a peaks/atac_ipsc_peaks.narrowPeak -b peaks/atac_neuron_peaks.narrowPeak
bedtools jaccard -a peaks/atac_ipsc_peaks.narrowPeak -b peaks/atac_ineuron_peaks.narrowPeak
bedtools jaccard -a peaks/atac_npc_peaks.narrowPeak -b peaks/atac_neuron_peaks.narrowPeak
bedtools jaccard -a peaks/atac_npc_peaks.narrowPeak -b peaks/atac_ineuron_peaks.narrowPeak
bedtools jaccard -a peaks/atac_neuron_peaks.narrowPeak -b peaks/atac_ineuron_peaks.narrowPeak

cat meta/ipsc.bamfiles.txt meta/ineuron.bamfiles.txt meta/npc.bamfiles.txt meta/neuron.bamfiles.txt > meta/multisample.bamfiles.txt
# Use relaxed settings to call peaks (FDR 10%)
submitJobs.py --MEM 6000 -j atacMacsCallPeak.ineuron.relaxed -q yesterday \
   -c "bash $SRC/atac/atac.macs_call_peaks.samples_file.sh meta/ineuron.bamfiles.txt peaks atac_ineuron_relaxed 0.1"
submitJobs.py --MEM 6000 -j atacMacsCallPeak.ipsc.relaxed -q yesterday \
   -c "bash $SRC/atac/atac.macs_call_peaks.samples_file.sh meta/ipsc.bamfiles.txt peaks atac_ipsc_relaxed 0.1"
submitJobs.py --MEM 6000 -j atacMacsCallPeak.npc.relaxed -q yesterday \
   -c "bash $SRC/atac/atac.macs_call_peaks.samples_file.sh meta/npc.bamfiles.txt peaks atac_npc_relaxed 0.1"
submitJobs.py --MEM 6000 -j atacMacsCallPeak.neuron.relaxed -q yesterday \
   -c "bash $SRC/atac/atac.macs_call_peaks.samples_file.sh meta/neuron.bamfiles.txt peaks atac_neuron_relaxed 0.1"
submitJobs.py --MEM 8000 -j atacMacsCallPeak.all.relaxed -q yesterday \
   -c "bash $SRC/atac/atac.macs_call_peaks.samples_file.sh meta/multisample.bamfiles.txt peaks atac_multisample_relaxed 0.1"

# For the "relaxed" peak calling settings, we expand peaks by 100 bp on each side
cd peaks
bedsum.sh atac_ineuron_relaxed_peaks.narrowPeak
bedsum.sh atac_ineuron_peaks.narrowPeak

for f in *_relaxed_peaks.narrowPeak; do
    cat $f | gzip > ${f}.orig.gz
    zcat ${f}.orig.gz | awk 'BEGIN{OFS="\t"}{print $1,$2-100,$3+100,$4,$5,$6,$7,$8,$9,$10}' > $f
done


###############################################################################
## RNA-seq
cd $SEQ/RNA

# Example command to get metadata - this is used within the scripts below
imeta qu -z seq -d study = 'OpenTargets AD PD finemap RNA' and target = 1 and manual_qc = 1

#Fetch list of lanelets from IRODS
python ~/src/utils/irods/irodsGetSamplesInStudy.py --studyName "OpenTargets AD PD finemap RNA" |  cut -f1 -d "." | uniq > irods.lanelets.txt

#Map lanelet ids to sample ids
python ~/src/utils/irods/irodsFetchMeta.py --irodsList irods.lanelets.txt | sort -k1 > irods.sample_lanes.txt 

#Fetch lanelets in cram format from irods
cut -f1 irods.lanelets.txt | python ~/src/utils/irods/fetch-irods.py --dir cram/ --suffix .cram
cut -f1 irods.lanelets.new.txt | python ~/src/utils/irods/fetch-irods.py --dir cram/ --suffix .cram

for f in cram/*.cram; do
  if [ ! -f $f.crai ]; then
    samtools index $f
  fi
done

# Extract the "genome" file for bedtools, i.e. chromosome lengths, from the CRAM header
samtools view -H "cram/23882_5#1.cram" | ~/src/utils/bam/genomeFileFromBamHeader.py > GRCh38.genome.txt

# Convert from CRAM to BigWig format for easier viewing in IGV
cat irods.sample_lanes.txt | submitJobs.py --MEM 2000 --jobname mergeBams \
    --command "~/src/utils/bam/mergeBams.py --indir cram --outdir . --insuffix .cram --outsuffix .bam"

grep "Successfully completed" FarmOut/mergeBams.*.txt | wc -l
grep -iP "Failed|TERM|error" FarmOut/mergeBams.*.txt | wc -l

cut -f 1 irods.sample_lanes.txt | submitJobs.py --MEM 1000 --jobname indexBams \
    --command "~/src/utils/bam/indexBams.py --bamdir . --insuffix .bam"

cut -f 1 irods.sample_lanes.txt | submitJobs.py --MEM 3000 --jobname bamCoverageToBigWig \
    --command "~/src/utils/coverage/bam2bigwig.py --genome GRCh38.genome.txt --indir . --split"
echo "4860STDY7028461" | submitJobs.py --MEM 3000 --jobname bamCoverageToBigWig \
    --command "~/src/utils/coverage/bam2bigwig.py --genome GRCh38.genome.txt --indir . --split"

# Use verifyBamID to ensure that the files match the expected HIPSCI sample
KOLF2_VCF=$SEQ/genotypes/kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.vcf.gz
cut -f 1 irods.sample_lanes.txt | submitJobs.py --MEM 1000 --jobname verifyBamID \
    -c  "python ~/src/utils/bam/runVerifyBamID.py --bamdir . --insuffix .bam --vcf $KOLF2_VCF --self" 

(echo -ne "Sample\t"; head -n 1 ./4860STDY7028458/4860STDY7028458.verifyBamID.bestSM) > RNA_samples.verifyBamID.bestSM
for f in ./*/*.verifyBamID.bestSM; do
    (echo -ne "$f\t"; sed '1d' $f) >> RNA_samples.verifyBamID.bestSM
done

# Make a file with chromosomes 1...X,Y,MT named chrcounts.all.txt
f=./4860STDY7028458/4860STDY7028458.chr_counts
cat $f | sed "s/^[ ]*//" | perl -ne '@l=split(" ", $_); if ($l[1] =~ "^[chr]*([0-9]+|M|MT|X|Y|\\*)") { print join("\t", @l)."\n"; }' | sortByChrPos.pl --file - --chrcol 2 --poscol 1 > $f.sorted.txt
(echo "chr"; cut -f 2 $f.sorted.txt) > ./chrcounts.all.txt
for f in ./*/*.chr_counts; do
    cat $f | sed "s/^[ ]*//" | perl -ne '@l=split(" ", $_); if ($l[1] =~ "^[chr]*([0-9]+|M|MT|X|Y|\\*)") { print join("\t", @l)."\n"; }' | sortByChrPos.pl --file - --chrcol 2 --poscol 1 > $f.sorted.txt
    fbname=`basename "$f" | cut -d. -f1`
    cp ./chrcounts.all.txt ./chrcounts.all.tmp
    paste ./chrcounts.all.tmp <(echo $fbname; cut -f 1 $f.sorted.txt) > ./chrcounts.all.txt
done


# Count reads mapping to genes
GENCODE_GTF=$JS/reference/GRCh38/gencode.v27.annotation.gtf
cut -f1 irods.sample_lanes.txt | submitJobs.py --MEM 1000 -j featureCounts -c "~/src/utils/counts/bam2counts.py --sampleDir . --gtf $GENCODE_GTF --strand 2 --countsSuffix .counts.txt --bamSuffix .bam --execute True"
grep "Successfully completed" FarmOut/featureCounts.*.txt | wc -l
grep -iP "Failed|TERM|error" FarmOut/featureCounts.*.txt | wc -l

# I noticed a high rate of reads not mapping to any feature (gene). The same was not
# true for my sensory neurons. I thought I would check whether I get a lower rate
# when counting using the same gene annotations as I did for sensory neurons.
GENCODE_BASIC_GTF=/lustre/scratch117/cellgen/team170/js29/sensoryneurons/annotations/Homo_sapiens.GRCh38.79.gencode_basic.gtf
cut -f1 irods.sample_lanes.txt | head -n 1 | submitJobs.py --MEM 1000 -j featureCounts -c "~/src/utils/counts/bam2counts.py --sampleDir . --gtf $GENCODE_BASIC_GTF --strand 2 --countsSuffix .counts.basic.txt --bamSuffix .bam --execute True"
# With the gencode_basic annotations, I actually get slightly worse results, so this
# does not explain the low rate of reads mapping to features.

# Look at the summary of gene read counts - I made this manually into a spreadsheet
cat ./*/*.counts.txt.summary

cd $SEQ/RNA/analysis
# Combine counts for samples into a single file
Rscript $SRC/rna/ipsneurons.combine_counts.R ../irods.sample_lanes.txt .. ipsneurons.counts.txt.gz

	

