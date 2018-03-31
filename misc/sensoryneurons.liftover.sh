JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
GRCh37_FASTA=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/reference/GRCh37/GRCh37.p13.genome.fa

cd $JS/sensoryneurons/
cat /warehouse/compgen_wh05/js29/sensoryneurons/atac/ATAC_consensus_peaks.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"peak_"NR}' > GRCh38/ATAC_consensus_peaks.bed

CrossMap.py bed $JS/software/CrossMap/GRCh38_to_GRCh37.chain.gz GRCh38/ATAC_consensus_peaks.bed GRCh37/ATAC_consensus_peaks.bed

bedtools merge -i <(cat GRCh37/ATAC_consensus_peaks.bed | sort -k1,1 -k2,2n) > GRCh37/ATAC_consensus_peaks.merged.bed

