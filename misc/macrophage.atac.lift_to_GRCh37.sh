JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

cd $JS/macrophage/GRCh37

zcat $JS/macrophage/ATAC_peak_metadata.txt.gz | sed '1d' | awk 'BEGIN{OFS="\t"}{print "chr"$4,$5,$6,$1}' > ../GRCh38/ATAC_macrophage.peaks.GRCh38.bed

CrossMap.py bed $JS/software/CrossMap/GRCh38_to_GRCh37.chain.gz \
    ../GRCh38/ATAC_macrophage.peaks.GRCh38.bed \
    ATAC_macrophage.peaks.GRCh37.bed.tmp

# Somehow CrossMap removes the "chr" from the bed file, so we add it back in
cat ATAC_macrophage.peaks.GRCh37.bed.tmp | perl -ne 'print "chr".$_' \
   > ATAC_macrophage.peaks.GRCh37.bed

rm ATAC_macrophage.peaks.GRCh37.bed.tmp
