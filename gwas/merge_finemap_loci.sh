#!/bin/bash

ROOT=/lustre/scratch115/realdata/mdt3/projects/otcoregen
OUTDIR=$ROOT/jeremys/gwas/AD/finemap
cd $OUTDIR

CREDSET=$ROOT/AD_PD_finemap/AD_credible_sets/

fullfile=$CREDSET/1_207750568.IGAP1_GWAX.1.set
filename=$(basename "$fullfile")
extension="${filename##*.}"
filebase=`echo $filename | perl -ne '@a=split(/\./); print $a[0]'`
echo $filebase

# Add "locus" as the first field of the header for the *.set file
# then paste the contents of each *.set file, with each line prefixed by its
# locus name (the root of the *.set filename)
head -n 1 $CREDSET/1_207750568.IGAP1_GWAX.1.set | perl -ane 'chomp; print join("\t", "locus", @F)."\n"' > $OUTDIR/AD.credible_sets.set
for f in $CREDSET/*.set; do
    filename=$(basename "$f")
    filebase=`echo $filename | perl -ne '@a=split(/\./); print $a[0]'`
    sed '1d' $f | perl -sane 'chomp; print join("\t", $locus, @F)."\n"' -- -locus=$filebase >> $OUTDIR/AD.credible_sets.set
done


FINEMAP=$ROOT/AD_PD_finemap/AD_finemap_output
# Add "locus" as the first field to the finemap *.snp output, similar to above
head -n 1 $FINEMAP/1_207750568.IGAP1_GWAX.1.snp | perl -ane 'chomp; print join("\t", "locus", @F)."\n"' > $OUTDIR/AD.finemap.snp
for f in $FINEMAP/*.snp; do
    filename=$(basename "$f")
    filebase=`echo $filename | perl -ne '@a=split(/\./); print $a[0]'`
    sed '1d' $f | perl -sane 'chomp; print join("\t", $locus, @F)."\n"' -- -locus=$filebase >> $OUTDIR/AD.finemap.snp
done

# Do the same for the max_causal_2 runs
head -n 1 $FINEMAP/max_causal_2/1_207750568.IGAP1_GWAX.2.snp | perl -ane 'chomp; print join("\t", "locus", @F)."\n"' > $OUTDIR/AD.finemap.max_causal_2.snp
for f in $FINEMAP/max_causal_2/*.snp; do
    filename=$(basename "$f")
    filebase=`echo $filename | perl -ne '@a=split(/\./); print $a[0]'`
    sed '1d' $f | perl -sane 'chomp; print join("\t", $locus, @F)."\n"' -- -locus=$filebase >> $OUTDIR/AD.finemap.max_causal_2.snp
done


################################################################################
# Do the same for Parkinson's

ROOT=/lustre/scratch115/realdata/mdt3/projects/otcoregen
OUTDIR=$ROOT/jeremys/gwas/PD/finemap
cd $OUTDIR

CREDSET=$ROOT/AD_PD_finemap/PD_finemap_output/

fullfile=$CREDSET/pd-gwas.caucasian.1_155135036.1.set
filename=$(basename "$fullfile")
extension="${filename##*.}"
filebase=`echo $filename | perl -ne '@a=split(/\./); print $a[2]'`
echo $filebase

# Add "locus" as the first field of the header for the *.set file
# then paste the contents of each *.set file, with each line prefixed by its
# locus name (the root of the *.set filename)
head -n 1 $CREDSET/pd-gwas.caucasian.1_155135036.1.set | perl -ane 'chomp; print join("\t", "locus", @F)."\n"' > $OUTDIR/PD.credible_sets.set
for f in $CREDSET/*.set; do
    filename=$(basename "$f")
    filebase=`echo $filename | perl -ne '@a=split(/\./); print $a[2]'`
    sed '1d' $f | perl -sane 'chomp; print join("\t", $locus, @F)."\n"' -- -locus=$filebase >> $OUTDIR/PD.credible_sets.set
done


FINEMAP=$ROOT/AD_PD_finemap/PD_finemap_output
# Add "locus" as the first field to the finemap *.snp output, similar to above
head -n 1 $FINEMAP/pd-gwas.caucasian.1_155135036.snp | perl -ane 'chomp; print join("\t", "locus", @F)."\n"' > $OUTDIR/PD.finemap.snp
for f in $FINEMAP/*.snp; do
    filename=$(basename "$f")
    filebase=`echo $filename | perl -ne '@a=split(/\./); print $a[2]'`
    sed '1d' $f | perl -sane 'chomp; print join("\t", $locus, @F)."\n"' -- -locus=$filebase >> $OUTDIR/PD.finemap.snp
done
