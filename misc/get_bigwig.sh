#!/bin/bash

JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

################### HIPSCI
cd $JS/ipsc/GRCh37/RNA

#for f in $JS/ipsc/GRCh37/RNA/HPSI*/*.bam.bai; do
for f in $JS/ipsc/GRCh37/RNA/HPSI*/*.bam; do
  DIR=`echo $f | perl -ne '@dir=split(/\//); print join("/", @dir[0..($#dir-1)]);'`
  NAME=`echo $f | perl -ne '@dir=split(/\//); @fparts=split(/\./, $dir[$#dir]); print $fparts[0];'`
  mv $f $DIR/$NAME.bam
done

samtools view -H "HPSI0114i-kolf_2/HPSI0114i-kolf_2.hs37d5.bwa.realigned.star.markdup.rnaseq.20160125.bam" | ~/src/utils/bam/genomeFileFromBamHeader.py > GRCh37.genome.txt

# Get bigwig files aware for splicing by using STAR
submitJobs.py --MEM 15000 -n 4 -j kolf_2.bam2bigwig -c "STAR --runMode inputAlignmentsFromBAM --inputBAMfile HPSI0114i-kolf_2/HPSI0114i-kolf_2.bam --outWigType wiggle --outWigStrand Stranded"
submitJobs.py --MEM 15000 -n 4 -j hehd_2.bam2bigwig -c "STAR --runMode inputAlignmentsFromBAM --inputBAMfile HPSI1213i-hehd_2/HPSI1213i-hehd_2.bam --outWigType wiggle --outWigStrand Stranded"
submitJobs.py --MEM 15000 -n 4 -j kolf_2.bam2bigwig -c "STAR --runMode inputAlignmentsFromBAM --inputBAMfile HPSI0114i-kolf_2/HPSI0114i-kolf_2.bam --outWigType wiggle --outWigStrand Stranded"
submitJobs.py --MEM 15000 -n 4 -j kolf_2.bam2bigwig -c "STAR --runMode inputAlignmentsFromBAM --inputBAMfile HPSI0114i-kolf_2/HPSI0114i-kolf_2.bam --outWigType wiggle --outWigStrand Stranded"
submitJobs.py --MEM 15000 -n 4 -j kolf_2.bam2bigwig -c "STAR --runMode inputAlignmentsFromBAM --inputBAMfile HPSI0114i-kolf_2/HPSI0114i-kolf_2.bam --outWigType wiggle --outWigStrand Stranded"
submitJobs.py --MEM 15000 -n 4 -j kolf_2.bam2bigwig -c "STAR --runMode inputAlignmentsFromBAM --inputBAMfile HPSI0114i-kolf_2/HPSI0114i-kolf_2.bam --outWigType wiggle --outWigStrand Stranded"

submitJobs.py --MEM 15000 -n 4 -j kolf_2.bam2bigwig -c "STAR --runMode inputAlignmentsFromBAM --inputBAMfile HPSI0114i-kolf_2/HPSI0114i-kolf_2.bam --outWigType wiggle --outWigStrand Unstranded"
submitJobs.py --MEM 15000 -n 4 -j hehd_2.bam2bigwig -c "STAR --runMode inputAlignmentsFromBAM --inputBAMfile HPSI1213i-hehd_2/HPSI1213i-hehd_2.bam --outWigType wiggle --outWigStrand Unstranded"

echo "HPSI0114i-kolf_2" | submitJobs.py --MEM 10000 -n 2 -j kolf_2.bam2bigwig -q yesterday -c "~/src/utils/coverage/rna_bam2bigwig.py --indir . --insuffix .bam --genome GRCh37.genome.txt --stranded"
echo "HPSI1213i-hehd_2" | submitJobs.py --MEM 10000 -n 2 -j hehd_2.bam2bigwig -q yesterday -c "~/src/utils/coverage/rna_bam2bigwig.py --indir . --insuffix .bam --genome GRCh37.genome.txt"


cut -f 1 ipsc.sample_names.3.txt | submitJobs.py --MEM 10000 -n 2 -j bam2bigwig -q normal -c "~/src/utils/coverage/rna_bam2bigwig.py --indir . --insuffix .bam --genome GRCh37.genome.txt --stranded"
find . -name *.bg -delete

# Wait for the above to finish before running these jobs! (Otherwise temp files collide??)
cut -f 1 ipsc.sample_names.3.txt | submitJobs.py --MEM 10000 -n 2 -j bam2bigwig -q normal -c "~/src/utils/coverage/rna_bam2bigwig.py --indir . --insuffix .bam --genome GRCh37.genome.txt"
find . -name *.bg -delete


cd bigwig
for f in $JS/ipsc/GRCh37/RNA/HPSI*/*.bw; do
  NAME=`echo $f | perl -ne '@dir=split(/\//); @fparts=split(/\./, $dir[$#dir]); print $fparts[0];'`
  ln -s $f $NAME.bw
done

#CrossMap.py bigwig $JS/software/CrossMap/GRCh37_to_GRCh38.chain.gz HPSI0114i-kolf_2.bw HPSI0114i-kolf_2.GRCh38
for f in $JS/ipsc/GRCh37/RNA/bigwig/*.bw; do
  NAME=`echo $f | perl -ne '@dir=split(/\//); @fparts=split(/\./, $dir[$#dir]); print $fparts[0];'`
  submitJobs.py --MEM 10000 -j crossMap.$NAME -q normal -c "CrossMap.py bigwig $JS/software/CrossMap/GRCh37_to_GRCh38.chain.gz $NAME.bw $NAME.GRCh38"
done


################### Macrophages
ln -s /lustre/scratch117/cellgen/team170/ka8/projects/Blood_ATAC/macrophages/bigwig $JS/macrophage/GRCh38/ATAC/bigwig


################### Sensory neurons
SN_RNA=/warehouse/compgen_wh05/js29/sn_hg38/STAR

cd $JS/sensoryneurons/GRCh38/RNA/bigwig
for f in $SN_RNA/*/*.bw; do
  DIR=`echo $f | perl -ne '@dir=split(/\//); print join("/", @dir[0..($#dir-1)]);'`
  FILE=`echo $f | perl -ne '@dir=split(/\//); print $dir[$#dir];'`
  echo $FILE
  ln -s $f $FILE
done

SN_ATAC=/warehouse/compgen_wh05/js29/sn_hg38/ATAC/bwa

cd $JS/sensoryneurons/GRCh38/ATAC/bigwig
for f in $SN_ATAC/*/*.bw; do
  DIR=`echo $f | perl -ne '@dir=split(/\//); print join("/", @dir[0..($#dir-1)]);'`
  FILE=`echo $f | perl -ne '@dir=split(/\//); print $dir[$#dir];'`
  echo $FILE
  ln -s $f $FILE
done


