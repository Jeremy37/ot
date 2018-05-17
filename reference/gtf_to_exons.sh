#!/bin/bash
JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

cd $JS/reference/GRCh38

# I ended up not using this code, as I used R biomaRt to get ensembl transcripts instead.
zcat gencode.v27.annotation.gtf.gz | perl -an -F/'\t'/ -e 'if ($F[2] =~ /gene/ and $F[8] =~ /gene_id "(ENSG[^"]+)"/) { print join("\t", $F[0],$F[3],$F[4],$F[6],$1)."\n"; }' \
    > gencode.v27.genes.bed
zcat gencode.v27.annotation.gtf.gz | perl -an -F/'\t'/ -e 'if ($F[2] =~ /gene/ and $F[8] =~ /gene_id "(ENSG[^"]+)"/ and $F[8] =~ /gene_type "protein_coding"/) { print join("\t", $F[0],$F[3],$F[4],$F[6],$1)."\n"; }' \
    > gencode.v27.genes.pc.bed

zcat gencode.v27.annotation.gtf.gz | perl -an -F/'\t'/ -e 'if ($F[2] =~ /exon/ and $F[8] =~ /gene_id "(ENSG[^"]+)"/) { print join("\t", $F[0],$F[3],$F[4],$F[6],$1)."\n"; }' \
    > gencode.v27.exons.bed


cd $JS/reference/GRCh37

# I ended up not using this code, as I used R biomaRt to get ensembl transcripts instead.
zcat gencode.v28lift37.annotation.gtf.gz | perl -an -F/'\t'/ -e 'if ($F[2] =~ /gene/ and $F[8] =~ /gene_id "(ENSG[^"]+).*gene_name "([^"]+)"/) { print join("\t", $F[0],$F[3],$F[4],$F[6],$1,$2)."\n"; }' \
    > gencode.v28lift37.genes.bed

zcat gencode.v28lift37.annotation.gtf.gz | perl -an -F/'\t'/ -e 'if ($F[2] =~ /gene/ and $F[8] =~ /gene_type "protein_coding"/ and $F[8] =~ /gene_id "(ENSG[^"]+).*gene_name "([^"]+)"/) { print join("\t", $F[0],$F[3],$F[4],$F[6],$1,$2)."\n"; }' \
    > gencode.v28lift37.genes.pc.bed

zcat gencode.v28lift37.annotation.gtf.gz | perl -an -F/'\t'/ -e 'if ($F[2] =~ /exon/ and $F[8] =~ /gene_id "(ENSG[^"]+).*gene_name "([^"]+)"/) { print join("\t", $F[0],$F[3],$F[4],$F[6],$1,$2)."\n"; }' \
    > gencode.v28lift37.exons.bed



