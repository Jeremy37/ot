#!/bin/bash

JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
cd $JS/coloc/qtl_data

################################################################################
# Sensory neurons
SN=/warehouse/compgen_wh05/js29/sensoryneurons

cd $JS/sensoryneurons
OUTDIR=$JS/coloc/qtl_data/sensory_neuron

# Get sensory neuron variant info with MAF, needed for coloc
zcat $SN/rasqual/imputed.97_samples.snps_indels.INFO_08.selected_samples.vcf.gz \
  | grep -v "^#" | cut -f 1-5,8 \
  | perl -ane 'if ($F[5] =~ /.*AC=([0-9]+);AN=([0-9]+)/) { $MAF = $1/$2; print join("\t", @F[0..4], $MAF)."\n" }' \
  | gzip > $OUTDIR/sensory_neuron.variant_info.GRCh38.gz

(echo -e "chr\tpos\trsid\tMAF";
 zcat $JS/sensoryneurons/GRCh37/imputed.97_samples.snps_indels.INFO_08.GRCh37.vcf.gz \
  | grep -v "^#" | cut -f 1-3,8 \
  | perl -ane 'if ($F[3] =~ /.*AC=([0-9]+);AN=([0-9]+)/) { $MAF = $1/$2; print join("\t", @F[0..2], $MAF)."\n" }') \
  | gzip > $OUTDIR/sensory_neuron.variant_info.GRCh37.gz


########### eQTLs
echo -e "feature\tchr\tpos\trsid" >  $JS/coloc/qtl_data/sensory_neuron.fastqtl.cis500k.qtl_signals.txt
sed '1d' $SN/fastqtl/output/permutations.10k.allgenes.cis500k.PCs.20.fdr0.1.coords.txt | awk 'BEGIN{OFS="\t"}{print $1,$13,$14,$6}' >> $OUTDIR/sensory_neuron.fastqtl.cis500k.qtl_signals.txt

cd $JS/sensoryneurons
SN_QTL_IN=$SN/fastqtl/output/nominals.allgenes.cis500k.ePCs.20.pvalues.notskipped.coords.txt.gz
zcat $SN_QTL_IN | awk 'BEGIN{OFS="\t"}{print $1,$6,$7,$2,$3,$4,$5}' | sort -k2,2 -k3,3g | bgzip > GRCh38/fastqtl.nominals.cis500k.ePCs.20.coords.txt.gz
tabix -s 2 -b 3 -e 3 GRCh38/fastqtl.nominals.cis500k.ePCs.20.coords.txt.gz

zcat $SN_QTL_IN | awk 'BEGIN{OFS="\t"}{print $1,$6,$7,$2,$4}' | sort -k2,2 -k3,3g | bgzip > GRCh38/sensory_neuron.fastqtl.cis500k.nominals.GRCh38.txt.gz
tabix -s 2 -b 3 -e 3 GRCh38/sensory_neuron.fastqtl.cis500k.nominals.GRCh38.txt.gz


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

#(echo -e "chr\tpos\trsid\tMAF";
# zcat GRCh37/imputed.97_samples.snps_indels.INFO_08.variant_info.GRCh37.gz | cut -f 1,2,3,6) \
# | bgzip > $JS/coloc/qtl_data/sensory_neuron.variant_info.for_coloc.GRCh37.gz

########### sQTLs
echo -e "feature\tchr\tpos\trsid" > $OUTDIR/sensory_neuron.sqtl.fastqtl.qtl_signals.txt
sed '1d' $JS/sensoryneurons/GRCh38/sqtl/fastqtl.sqtl.permutations.10k.PCs.5.fdr0.1.txt | awk 'BEGIN{OFS="\t"}{print $2,$5,$6,$4}' >> $OUTDIR/sensory_neuron.sqtl.fastqtl.qtl_signals.txt

zcat GRCh38/sqtl/fastqtl.sqtl.nominals.PCs.5.pvalues.coords.txt.gz | cut -f 1,2,3,4,6 | bgzip > $OUTDIR/sensory_neuron.sqtl.fastqtl.nominals.GRCh38.txt.gz
tabix -s 2 -b 3 -e 3 $OUTDIR/sensory_neuron.sqtl.fastqtl.nominals.GRCh38.txt.gz

zcat $SN/leafcutter/fastqtl/output.500k/fastqtl.nominals.500k.PCs.5.pvalues.coords.txt.gz | awk 'BEGIN{OFS="\t"}{print $1,$6,$7,$2,$4}' | sort -k2,2 -k3,3n | bgzip > $OUTDIR/sensory_neuron.sqtl.fastqtl.500k.nominals.GRCh38.txt.gz
tabix -s 2 -b 3 -e 3 $OUTDIR/sensory_neuron.sqtl.fastqtl.500k.nominals.GRCh38.txt.gz


################################################################################
# Blueprint
cd $JS/coloc/qtl_data

submitJobs.py --MEM 500 -j process_bluprint.mono_gene_eQTL -q yesterday \
   -c "$JS/src/coloc/process_blueprint_file.sh $JS/datasets/blueprint/mono_gene_nor_combat_peer_10_all_summary.txt.gz $JS/coloc/qtl_data/blueprint/blueprint.mono_gene_eQTL"

submitJobs.py --MEM 1000 -j process_bluprint.mono_K27AC -q yesterday \
   -c "$JS/src/coloc/process_blueprint_file.sh $JS/datasets/blueprint/mono_K27AC_log2rpm_peer_10_all_summary.txt.gz $JS/coloc/qtl_data/blueprint/blueprint.mono_K27AC"

submitJobs.py --MEM 1000 -j process_bluprint.mono_K4ME1 -q yesterday \
   -c "$JS/src/coloc/process_blueprint_file.sh $JS/datasets/blueprint/mono_K4ME1_log2rpm_peer_10_all_summary.txt.gz $JS/coloc/qtl_data/blueprint/blueprint.mono_K4ME1"


################################################################################
# xQTL
cd $JS/coloc/qtl_data

submitJobs.py --MEM 10000 -j process_xqtl.eQTLs -q normal \
   -c "$JS/src/coloc/process_xqtl_file.sh $JS/reference/xQTL/eQTLs_all.txt.gz $JS/reference/xQTL/maf.gz $JS/coloc/qtl_data/xQTL/xQTL_eQTL"

submitJobs.py --MEM 10000 -j process_xqtl.mQTLs -q normal \
   -c "$JS/src/coloc/process_xqtl_file.sh $JS/reference/xQTL/mQTLs_all.txt.gz $JS/reference/xQTL/maf.gz $JS/coloc/qtl_data/xQTL/xQTL_mQTL"

submitJobs.py --MEM 10000 -j process_xqtl.haQTLs -q normal \
   -c "$JS/src/coloc/process_xqtl_file.sh $JS/reference/xQTL/haQTLs_all.txt.gz $JS/reference/xQTL/maf.gz $JS/coloc/qtl_data/xQTL/xQTL_haQTL"


################################################################################
# iPSC-derived macrophages
KAUR_BASE=/lustre/scratch117/cellgen/team170/ka8/projects/macrophage-trQTLs
KAUR_QTL=$KAUR_BASE/processed/salmonella/qtltools/output

KAUR_VARIANT_INFO=$KAUR_BASE/results/genotypes/salmonella/imputed.86_samples.variant_information.txt.gz
(echo -e "chr\tpos\trsid\tMAF"; \
 zcat $KAUR_VARIANT_INFO \
  | perl -ane 'print join("\t", @F[0..2],$F[6]/$F[7])."\n";' \
  | sort -k1,1 -k2,2n | uniq) \
  | bgzip > $JS/coloc/qtl_data/ips_macrophage.variant_info.GRCh38.txt.gz

KAUR_VARIANT_INFO=$KAUR_BASE/results/genotypes/salmonella/GRCh37/imputed.86_samples.variant_information.GRCh37.txt.gz
(echo -e "chr\tpos\trsid\tMAF"; \
 zcat $KAUR_VARIANT_INFO \
  | perl -ane 'print join("\t", @F[0..2],$F[6]/$F[7])."\n";' \
  | sort -k1,1 -k2,2n | uniq) \
  | bgzip > $JS/coloc/qtl_data/ips_macrophage.variant_info.GRCh37.txt.gz

STATE=naive
STATE=IFNg
STATE=SL1344
STATE=IFNg_SL1344

# Run the below commands for each macrophage states above
(echo -e "feature\tchr\tpos\trsid\tp.value\tbeta_p.value\tpheno_start\tpheno_end\tpheno_strand\tnum_tests\tdist"; \
 zcat $KAUR_QTL/featureCounts/$STATE.permuted.txt.gz | awk 'BEGIN{OFS="\t"}{print $6,$11,$12,$10,$18,$21,$3,$4,$5,$8,$9}' | awk '$6 < 0.001') \
  > $JS/coloc/qtl_data/macrophage/ips_macrophage.eQTL_$STATE.signals.beta_lt_1e-3.txt

(echo -e "feature\tchr\tpos\trsid\tp.value\tbeta_p.value\tpheno_start\tpheno_end\tpheno_strand\tnum_tests\tdist"; \
 zcat $KAUR_QTL/leafcutter/$STATE.permuted.txt.gz | awk 'BEGIN{OFS="\t"}{print $6,$11,$12,$10,$18,$21,$3,$4,$5,$8,$9}' | awk '$6 < 0.001') \
  > $JS/coloc/qtl_data/macrophage/ips_macrophage.leafcutter_QTL_$STATE.signals.beta_lt_1e-3.txt

(echo -e "feature\tchr\tpos\trsid\tp.value\tbeta_p.value\tpheno_start\tpheno_end\tpheno_strand\tnum_tests\tdist"; \
 zcat $KAUR_QTL/txrevise_ends/$STATE.permuted.txt.gz | awk 'BEGIN{OFS="\t"}{print $6,$11,$12,$10,$18,$21,$3,$4,$5,$8,$9}' | awk '$6 < 0.001') \
  > $JS/coloc/qtl_data/macrophage/ips_macrophage.txrevise_ends_$STATE.signals.beta_lt_1e-3.txt

(echo -e "feature\tchr\tpos\trsid\tp.value\tbeta_p.value\tpheno_start\tpheno_end\tpheno_strand\tnum_tests\tdist"; \
 zcat $KAUR_QTL/txrevise_promoters/$STATE.permuted.txt.gz | awk 'BEGIN{OFS="\t"}{print $6,$11,$12,$10,$18,$21,$3,$4,$5,$8,$9}' | awk '$6 < 0.001') \
  > $JS/coloc/qtl_data/macrophage/ips_macrophage.txrevise_promoters_$STATE.signals.beta_lt_1e-3.txt


submitJobs.py --MEM 1000 -j process_ips_macrophage.eQTLs.$STATE -q normal \
   -c "$JS/src/coloc/process_ips_macrophage_file.sh $KAUR_QTL/featureCounts/sorted/$STATE.nominal.sorted.txt.gz $JS/coloc/qtl_data/macrophage/ips_macrophage.eQTL_$STATE"

submitJobs.py --MEM 1000 -j process_ips_macrophage.leafcutter_QTLs.$STATE -q normal \
   -c "$JS/src/coloc/process_ips_macrophage_file.sh $KAUR_QTL/leafcutter/sorted/$STATE.nominal.sorted.txt.gz $JS/coloc/qtl_data/macrophage/ips_macrophage.leafcutter_QTL_$STATE"

submitJobs.py --MEM 1000 -j process_ips_macrophage.txrevise_ends.$STATE -q normal \
   -c "$JS/src/coloc/process_ips_macrophage_file.sh $KAUR_QTL/txrevise_ends/sorted/$STATE.nominal.sorted.txt.gz $JS/coloc/qtl_data/macrophage/ips_macrophage.txrevise_ends_$STATE"

submitJobs.py --MEM 1000 -j process_ips_macrophage.txrevise_promoters.$STATE -q normal \
   -c "$JS/src/coloc/process_ips_macrophage_file.sh $KAUR_QTL/txrevise_promoters/sorted/$STATE.nominal.sorted.txt.gz $JS/coloc/qtl_data/macrophage/ips_macrophage.txrevise_promoters_$STATE"


# Run the below commands for each macrophage states above
(echo -e "feature\tchr\tpos\trsid\tp.value\tbeta_p.value\tpheno_start\tpheno_end\tpheno_strand\tnum_tests\tdist"; \
 zcat $KAUR_QTL/featureCounts/$STATE.permuted.txt.gz | awk 'BEGIN{OFS="\t"}{print $6,$11,$12,$10,$18,$21,$3,$4,$5,$8,$9}' | awk '$6 < 0.05') \
  > $JS/coloc/qtl_data/macrophage/ips_macrophage.eQTL_$STATE.signals.beta_lt_0.05.txt

(echo -e "feature\tchr\tpos\trsid\tp.value\tbeta_p.value\tpheno_start\tpheno_end\tpheno_strand\tnum_tests\tdist"; \
 zcat $KAUR_QTL/leafcutter/$STATE.permuted.txt.gz | awk 'BEGIN{OFS="\t"}{print $6,$11,$12,$10,$18,$21,$3,$4,$5,$8,$9}' | awk '$6 < 0.05') \
  > $JS/coloc/qtl_data/macrophage/ips_macrophage.leafcutter_QTL_$STATE.signals.beta_lt_0.05.txt

(echo -e "feature\tchr\tpos\trsid\tp.value\tbeta_p.value\tpheno_start\tpheno_end\tpheno_strand\tnum_tests\tdist"; \
 zcat $KAUR_QTL/txrevise_ends/$STATE.permuted.txt.gz | awk 'BEGIN{OFS="\t"}{print $6,$11,$12,$10,$18,$21,$3,$4,$5,$8,$9}' | awk '$6 < 0.05') \
  > $JS/coloc/qtl_data/macrophage/ips_macrophage.txrevise_ends_$STATE.signals.beta_lt_0.05.txt

(echo -e "feature\tchr\tpos\trsid\tp.value\tbeta_p.value\tpheno_start\tpheno_end\tpheno_strand\tnum_tests\tdist"; \
 zcat $KAUR_QTL/txrevise_promoters/$STATE.permuted.txt.gz | awk 'BEGIN{OFS="\t"}{print $6,$11,$12,$10,$18,$21,$3,$4,$5,$8,$9}' | awk '$6 < 0.05') \
  > $JS/coloc/qtl_data/macrophage/ips_macrophage.txrevise_promoters_$STATE.signals.beta_lt_0.05.txt


