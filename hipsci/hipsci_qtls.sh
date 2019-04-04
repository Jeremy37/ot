###############################################################################
# This script contains commands used to annotate and fine-map QTLs from the
# HIPSCI project, using summary stats provided by Marc-Jan Bonder.

JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
INPUT=$JS/datasets/hipsci/Input
ROOT=$JS/ipsc/hipsci_finemap

cd $ROOT

# First we need to sort the files provided so that we can extract regions
# and genes of interest
mkdir $ROOT/SplicingLevel

CHR=22
TYPE=SplicingLevel
TYPE=ApaLevel
mkdir $ROOT/$TYPE
for CHR in `seq 1 22`; do
  echo "Chr: $CHR"
  F=$INPUT/$TYPE/full_qtl_results_$CHR.Pval-rescaled.txt.gz
  (gzhead 1 $F; zcat $F | sed '1d' | sort -k1,1 -k2,2n) | sed 's/:/~/g' | bgzip > $ROOT/$TYPE/full_qtl_results_$CHR.rescaled.gene_sorted.txt.gz
  tabix -s 1 -b 2 -e 2 -S 1 $ROOT/$TYPE/full_qtl_results_$CHR.rescaled.gene_sorted.txt.gz
done

TYPE=GeneLevel
mkdir $ROOT/$TYPE
for CHR in `seq 1 22`; do
  echo "Chr: $CHR"
  F=$INPUT/$TYPE/full_qtl_results_$CHR.txt.gz
  (gzhead 1 $F; zcat $F | sed '1d' | sort -k1,1 -k2,2n) | sed 's/:/~/g' | bgzip > $ROOT/$TYPE/full_qtl_results_$CHR.gene_sorted.txt.gz
  tabix -s 1 -b 2 -e 2 -S 1 $ROOT/$TYPE/full_qtl_results_$CHR.gene_sorted.txt.gz
done

time tabix $ROOT/SplicingLevel/full_qtl_results_$CHR.rescaled.gene_sorted.txt.gz "chr22~45316370~45362895~clu_13974:45118934-45122866" | wc -l

# Call R script to go through all QTL associations and print out credible set SNPs
TYPE=SplicingLevel
TYPE=ApaLevel
for CHR in `seq 1 22`; do
  echo "Chr: $CHR"
  F=$ROOT/$TYPE/full_qtl_results_$CHR.rescaled.gene_sorted.txt.gz
  submitJobs.py --MEM 8000 -j finemap_qtl.$TYPE -q normal \
    -c "Rscript $JS/src/hipsci/finemap_qtl_wtccc.R --input $F --out $ROOT/$TYPE/$TYPE.fine.$CHR"
done

TYPE=GeneLevel
for CHR in `seq 1 22`; do
  echo "Chr: $CHR"
  F=$ROOT/$TYPE/full_qtl_results_$CHR.gene_sorted.txt.gz
  submitJobs.py --MEM 8000 -j finemap_qtl.$TYPE -q normal \
    -c "Rscript $JS/src/hipsci/finemap_qtl_wtccc.R --input $F --out $ROOT/$TYPE/$TYPE.fine.$CHR"
done

grep "Successfully" FarmOut/finemap_qtl*.txt | wc -l
grep -iP "Fail|ERROR|Abort|exit code" FarmOut/finemap_qtl*.txt | wc -l


head -n 1 $ROOT/$TYPE/$TYPE.fine.1.summary.tsv > $ROOT/$TYPE/$TYPE.fine.allchrs.summary.tsv
gzhead 1 $ROOT/$TYPE/$TYPE.fine.1.all.tsv > $ROOT/$TYPE/$TYPE.fine.allchrs.allsnps.tsv
for CHR in `seq 1 22`; do
  sed '1d' $ROOT/$TYPE/$TYPE.fine.$CHR.summary.tsv >> $ROOT/$TYPE/$TYPE.fine.allchrs.summary.tsv
  zcat $ROOT/$TYPE/$TYPE.fine.$CHR.all.tsv | sed '1d' >> $ROOT/$TYPE/$TYPE.fine.allchrs.allsnps.tsv
done


body() {
    IFS= read -r header
    printf '%s\n' "$header"
    "$@"
}
cat $ROOT/$TYPE/$TYPE.fine.allchrs.summary.tsv | body sort -k 5,5rg > $ROOT/$TYPE/$TYPE.fine.allchrs.summary.sorted.tsv
gzip $ROOT/$TYPE/$TYPE.fine.allchrs.allsnps.tsv
zcat $ROOT/$TYPE/$TYPE.fine.allchrs.allsnps.tsv.gz | body sort -k 5,5rg | gzip > $ROOT/$TYPE/$TYPE.fine.allchrs.allsnps.sorted.tsv.gz



CHR=22
F=$ROOT/$TYPE/full_qtl_results_$CHR.rescaled.gene_sorted.txt.gz
Rscript $JS/src/hipsci/finemap_qtl_wtccc.R --input $F --out $ROOT/$TYPE/$TYPE.fine.$CHR


TYPE=GeneLevel
F=$ROOT/$TYPE/full_qtl_results_$CHR.gene_sorted.txt.gz
Rscript $JS/src/hipsci/finemap_qtl_wtccc.R --input $F --out $ROOT/$TYPE/$TYPE.fine.$CHR
