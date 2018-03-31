#!/bin/bash

JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
SN=/warehouse/compgen_wh05/js29/sensoryneurons
cd $JS/sensoryneurons/GRCh38

QTL_IN=$SN/fastqtl/output/nominals.allgenes.cis500k.ePCs.20.pvalues.notskipped.coords.txt.gz

zcat $QTL_IN | awk 'BEGIN{OFS="\t"}{print $1,$6,$7,$2,$3,$4,$5}' | sort -k2,2 -k3,3g | bgzip > fastqtl.nominals.cis500k.ePCs.20.coords.txt.gz
tabix -s 2 -b 3 -e 3 fastqtl.nominals.cis500k.ePCs.20.coords.txt.gz

echo -e "phenotype_id\tntests\tmle_shape1\tmle_shape2\tdummy\tsnp_id\tgenedist\tp_val\tbeta\tperm_p\tperm_p_beta\tFDR" > fastqtl.permutations.10k.allgenes.cis100k.PCs.20.fdr0.1.txt 
sed '1d' $SN/fastqtl/output/permutations.10k.allgenes.cis100k.PCs.20.fdr0.1.txt | tr ' ' '\t' >> fastqtl.permutations.10k.allgenes.cis100k.PCs.20.fdr0.1.txt 

echo -e "phenotype_id\tntests\tmle_shape1\tmle_shape2\tdummy\tsnp_id\tgenedist\tp_val\tbeta\tperm_p\tperm_p_beta\tFDR" > fastqtl.permutations.10k.allgenes.cis500k.PCs.20.fdr0.1.txt 
sed '1d' $SN/fastqtl/output/permutations.10k.allgenes.cis500k.PCs.20.fdr0.1.txt | tr ' ' '\t' >> fastqtl.permutations.10k.allgenes.cis500k.PCs.20.fdr0.1.txt 

# For colocalisation analyses we need a file of minp values per gene, a file
# with information on the variants tested, and a file with the nominal p values
# for all tests.
# Instead of the below, we use the permutations output file from FastQTL

#zcat $QTL_IN | getLeadSnps.pl -f stdin --genecol 1 --pcol 4 \
#  | awk 'BEGIN{OFS="\t"}{ print $1,$6,$7,$2,$3,$4,$5 }' \
#  | gzip > fastqtl.nominals.cis500k.ePCs.20.gene_minp.txt.gz


# Also get summary stats for splice QTLs
#echo -e "phenotype_id\tntests\tmle_shape1\tmle_shape2\tdummy\tsnp_id\tdist\tp_val\tbeta\tperm_p\tperm_p_beta\tchr\tpos\tbonf_pval\tcluster_size\tFDR" > fastqtl.sqtl.permutations.10k.PCs.5.fdr0.1.txt
#sed '1d' $SN/leafcutter/fastqtl/output/permutations.10k.PCs.5.fdr0.1.new.txt >> fastqtl.sqtl.permutations.10k.PCs.5.fdr0.1.txt
echo -e "cluster_and_gene\tphenotype_id\tntests\tsnp_id\tchr\tpos\tdist\tp_val\tcluster_size\tFDR\tensembl_ids\tsymbols" > fastqtl.sqtl.permutations.10k.PCs.5.fdr0.1.txt
sed '1d' $SN/leafcutter/fastqtl/sqtls.fdr0.1.genes.annotated.cut.txt \
    | perl -ane 'print join("\t", "$F[10];$F[9];$F[0]", @F)."\n";' \
    >> fastqtl.sqtl.permutations.10k.PCs.5.fdr0.1.txt

zcat $SN/leafcutter/fastqtl/output.new/fastqtl.nominals.PCs.5.pvalues.coords.txt.gz | awk 'BEGIN{OFS="\t"}{print $1,$6,$7,$2,$3,$4,$5}' | sort -k2,2 -k3,3g | bgzip > fastqtl.sqtl.nominals.PCs.5.pvalues.coords.txt.gz
tabix -s 2 -b 3 -e 3 fastqtl.sqtl.nominals.PCs.5.pvalues.coords.txt.gz




# Get sensory neuron variant info for coloc
zcat $SN/rasqual/imputed.97_samples.snps_indels.INFO_08.selected_samples.vcf.gz \
  | grep -v "^#" | cut -f 1-5,8 \
  | perl -ane 'if ($F[5] =~ /.*AC=([0-9]+);AN=([0-9]+)/) { $MAF = $1/$2; print join("\t", @F[0..4], $MAF)."\n" }' \
  | gzip > imputed.97_samples.snps_indels.INFO_08.variant_info.GRCh38.gz


# Update GRCh37 reference fasta to not include "chr" in chromosome names
#zcat $JS/reference/GRCh37/GRCh37.p13.genome.fa.gz | sed -e 's/chr//' | fold -w 80 | bgzip > $JS/reference/GRCh37/GRCh37.p13.genome.nochr.fa.gz


# Run CrossMap
# I had a lot of trouble when running on the full file, because some sequences in the
# CrossMap chain file are not present in the GRCh37 fasta file, e.g. HG1292_PATCH.
# I can't figure out how to get the sequence for HG1292_PATCH, and no existing fasta
# download of the human genome seems to have it. I therefore am trying alternative chain
# files, such as the one from UCSC, hg38ToHg19.over.chain.gz. However, this differs in
# that chromosomes are named e.g. "chr1" rather than "1", and so necessitates renaming
# chromosomes in my VCF file.
cd $JS/sensoryneurons

## This works for the test, but fails for the full VCF file, due to regions which cover
# some "patch" regions of chr1 which aren't present in the Fasta file.
# CrossMap.py vcf $JS/software/CrossMap/GRCh38_to_GRCh37.chain.gz <(zcat GRCh38/imputed.97_samples.snps_indels.INFO_08.cut.vcf.gz | head -n 200) \
#     $JS/reference/GRCh37/GRCh37.p13.genome.nochr.fa.gz \
#     GRCh37/imputed.97_samples.snps_indels.INFO_08.GRCh37.test.vcf
# submitJobs.py --MEM 1000 -j CrossMap.sensoryneurons.GRCh38_to_GRCh37 -q yesterday -n 2 \
#     -c "CrossMap.py vcf $JS/software/CrossMap/GRCh38_to_GRCh37.chain.gz GRCh38/imputed.97_samples.snps_indels.INFO_08.cut.vcf.gz $JS/reference/GRCh37/GRCh37.p13.genome.nochr.fa.gz GRCh37/imputed.97_samples.snps_indels.INFO_08.GRCh37.vcf"

zcat $SN/rasqual/imputed.97_samples.snps_indels.INFO_08.selected_samples.vcf.gz | cut -f 1-9 \
    | perl -ne 'if (/^#/) { print; } else { print "chr".$_; }' \
    | gzip > GRCh38/imputed.97_samples.snps_indels.INFO_08.chr.cut.vcf.gz

CrossMap.py vcf $JS/software/CrossMap/hg38ToHg19.over.chain.gz <(zcat GRCh38/imputed.97_samples.snps_indels.INFO_08.chr.cut.vcf.gz | head -n 200) \
    $JS/reference/GRCh37/GRCh37.p13.genome.fa.gz \
    GRCh37/imputed.97_samples.snps_indels.INFO_08.GRCh37.test.vcf

submitJobs.py --MEM 1000 -j CrossMap.sensoryneurons.GRCh38_to_GRCh37 -q yesterday -n 2 \
    -c "CrossMap.py vcf $JS/software/CrossMap/hg38ToHg19.over.chain.gz GRCh38/imputed.97_samples.snps_indels.INFO_08.chr.cut.vcf.gz $JS/reference/GRCh37/GRCh37.p13.genome.fa.gz GRCh37/imputed.97_samples.snps_indels.INFO_08.GRCh37.chr.vcf"
gzip GRCh37/imputed.97_samples.snps_indels.INFO_08.GRCh37.chr.vcf
zcat GRCh37/imputed.97_samples.snps_indels.INFO_08.GRCh37.chr.vcf.gz | sed -e 's/chr//' | bgzip > GRCh37/imputed.97_samples.snps_indels.INFO_08.GRCh37.vcf.gz



# Get sensory neuron variant info in GRCh37 coords for coloc (needed to match
# positions between GWAS in GRCh37 coords and sensory neurons QTLs in GRCh38).
zcat GRCh37/imputed.97_samples.snps_indels.INFO_08.GRCh37.vcf.gz \
  | grep -v "^#" | cut -f 1-5,8 \
  | perl -ane 'if ($F[5] =~ /.*AC=([0-9]+);AN=([0-9]+)/) { $MAF = $1/$2; print join("\t", @F[0..4], $MAF)."\n" }' \
  | gzip > GRCh37/imputed.97_samples.snps_indels.INFO_08.variant_info.GRCh37.gz

