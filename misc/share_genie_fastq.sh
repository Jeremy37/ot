#!/bin/bash
OT=/lustre/scratch115/realdata/mdt3/projects/otcoregen
JS=$OT/jeremys
DEST=/lustre/scratch117/cellgen/team170/js29/GenIE

cd $DEST/fastq

DIR=$JS/experiment/transcribed
mkdir batch1
for f in $DIR/batch1/fastq/*.fastq.gz; do
  filename=$(basename "$f")
  echo $filename
  cp $f batch1/$filename
done

mkdir batch1_redo
for f in $DIR/batch1_redo/fastq2/*.fastq.gz; do
  filename=$(basename "$f")
  echo $filename
  cp $f batch1_redo/$filename
done

mkdir batch2
for f in $DIR/batch2/fastq/*.fastq.gz; do
  filename=$(basename "$f")
  echo $filename
  cp $f batch2/$filename
done

mkdir clu_sipa1l2
for f in $DIR/clu_sipa1l2/fastq/*.fastq.gz; do
  filename=$(basename "$f")
  echo $filename
  cp $f clu_sipa1l2/$filename
done

