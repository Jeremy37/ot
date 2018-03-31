#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

INPUT_BED=$1

#mkdir JEME
#mkdir JEME/roadmap
for F in $JS/annotation/JEME/encoderoadmap_lasso/*.csv; do
  fname=`basename $F`
  ID=$(echo $fname | cut -f 2 -d '.')
  bedtools intersect -wb -a $INPUT_BED -b <(cat $F | perl -F',' -ane '@x=split(/:|-/, $F[0]); print join("\t", @x, $F[1], $F[2]);') \
    | cut -f 1,3,8,9 > JEME/$ID.overlap.txt
done

