JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
PD=$JS/PD

cd $PD
cat Nalls_2019_genes.clumped.tsv | awk 'BEGIN{FS="\t";OFS="\t"}{print "chr"$1,$2-500000,$3+500000,$4,$5,$6,$7}' > Nalls_2019_genes.clumped.1Mb_window.chr.bed

echo -e "chr\twindowStart\twindowEnd\tlead_snp\tNearest_Gene\tQTL_Nominated_Gene\tlocus_number\tgene_start\tgene_end\tstrand\tgene_id\tsymbol" > Nalls_2019_genes.1Mb_window.gene_overlaps.tsv
bedtools intersect -a Nalls_2019_genes.clumped.1Mb_window.chr.bed -b $JS/AD_finemap/reference/gencode.v29lift37.genes.bed -wa -wb \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13}' >> Nalls_2019_genes.1Mb_window.gene_overlaps.tsv

Rscript $JS/src/misc/get_pd_gene_expression.R --genes Nalls_2019_genes.1Mb_window.gene_overlaps.tsv --expr $JS/reference/tissueExpr/tissues.selected.tpm.tsv.gz \
   > Nalls_2019_genes.1Mb_window.geneTPM.tsv

# Do the same for protein coding genes only
echo -e "chr\twindowStart\twindowEnd\tlead_snp\tNearest_Gene\tQTL_Nominated_Gene\tlocus_number\tgene_start\tgene_end\tstrand\tgene_id\tsymbol" > Nalls_2019_genes.1Mb_window.gene_overlaps.pc.tsv
bedtools intersect -a Nalls_2019_genes.clumped.1Mb_window.chr.bed -b $JS/AD_finemap/reference/gencode.v29lift37.genes.pc.bed -wa -wb \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13}' >> Nalls_2019_genes.1Mb_window.gene_overlaps.pc.tsv

Rscript $JS/src/misc/get_pd_gene_expression.R --genes Nalls_2019_genes.1Mb_window.gene_overlaps.pc.tsv --expr $JS/reference/tissueExpr/tissues.selected.tpm.tsv.gz \
   > Nalls_2019_genes.1Mb_window.geneTPM.pc.tsv

