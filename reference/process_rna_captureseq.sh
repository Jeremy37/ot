#!/bin/bash

JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
CAPTURESEQ=$JS/datasets/RNA_CaptureSeq/mattick_neuropsych

cd $CAPTURESEQ

# The value they have put in the bed file "score" field is not very clear.
# Martin Smith replied to my question saying that it was the length of the
# transcript, but this isn't a very useful score. Since bigBed format requires
# the score to be in the range 0-1000, I'm changing all score values to 1000.
cat GSE118158_filtered_hybrid_transcriptome.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,1000,$6,$7,$8,$9,$10,$11,$12}' > GSE118158_filtered_hybrid_transcriptome.score1000.bed

# I copied the standard chromosome details from $JS/reference/GRCh38/GRCh38.chr.fa.fai
# to get the GRCh38.genomesizes.txt file.
# Now convert to bigBed for compatibility with being used as a track hub.
bedToBigBed GSE118158_filtered_hybrid_transcriptome.score1000.bed GRCh38.genomesizes.txt GSE118158_filtered_hybrid_transcriptome.bb --type=bed12


# Convert the main bed file of RNA annotation to one with exons, rather than
# whole genes. This will enable us to annotate the distance of a snp to the
# nearest exon.
python ~/src/python/bed12_to_bed4.py --file GSE118158_filtered_hybrid_transcriptome.bed | sort -k1,1 -k2,2n > GSE118158_filtered_hybrid_transcriptome.exons.chr.tmp.bed
perl ~/src/sortByChrPos.pl --file GSE118158_filtered_hybrid_transcriptome.exons.chr.tmp.bed --chrcol 1 --poscol 2 > GSE118158_filtered_hybrid_transcriptome.exons.chr.bed

bedtools merge -i GSE118158_filtered_hybrid_transcriptome.exons.chr.bed > GSE118158_filtered_hybrid_transcriptome.exons.chr.merged.bed

cat GSE118158_filtered_hybrid_transcriptome.exons.chr.bed | sed 's/^chr//g' > GSE118158_filtered_hybrid_transcriptome.exons.bed
cat GSE118158_filtered_hybrid_transcriptome.exons.chr.merged.bed | sed 's/^chr//g' > GSE118158_filtered_hybrid_transcriptome.exons.merged.bed

# Make a bed file which has just the locations of every exon boundary as a 1-bp feature
cat GSE118158_filtered_hybrid_transcriptome.exons.bed | perl -ane 'print join("\t", $F[0],$F[1],$F[1]+1,$F[3])."\n"; print join("\t", $F[0],$F[2]-1,$F[2],$F[3])."\n";' > GSE118158_filtered_hybrid_transcriptome.spliceBoundaries.bed
cat GSE118158_filtered_hybrid_transcriptome.exons.chr.bed | perl -ane 'print join("\t", $F[0],$F[1],$F[1]+1,$F[3])."\n"; print join("\t", $F[0],$F[2]-1,$F[2],$F[3])."\n";' > GSE118158_filtered_hybrid_transcriptome.spliceBoundaries.chr.bed

cp GSE118158_filtered_hybrid_transcriptome.exons.bed $JS/AD_finemap/reference
cp GSE118158_filtered_hybrid_transcriptome.spliceBoundaries.bed $JS/AD_finemap/reference
