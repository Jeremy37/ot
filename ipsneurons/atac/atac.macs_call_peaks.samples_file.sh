#!/bin/bash
SAMPLE_LIST_FILE=$1
OUTDIR=$2
OUTNAME=$3
FDR=$4

FILES=`cat $SAMPLE_LIST_FILE | tr '\n' ' '`

#DIR=`dirname $INFILE`
#FNAME=`basename $INFILE`
#FBASE="${FNAME%.*}"

macs2 callpeak --nomodel --shift -37 --extsize 75 -q $FDR \
       --outdir $OUTDIR \
       -n $OUTNAME \
       -t $FILES
