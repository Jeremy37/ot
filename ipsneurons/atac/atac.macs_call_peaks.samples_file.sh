#!/bin/bash
SAMPLE_LIST_FILE=$1
OUTDIR=$2
OUTNAME=$3

FILES=`cat $SAMPLE_LIST_FILE | tr '\n' ' '`

#DIR=`dirname $INFILE`
#FNAME=`basename $INFILE`
#FBASE="${FNAME%.*}"

macs2 callpeak --nomodel --shift -37 --extsize 75 -q 0.01 \
       --outdir $OUTDIR \
       -n $OUTNAME \
       -t $FILES
