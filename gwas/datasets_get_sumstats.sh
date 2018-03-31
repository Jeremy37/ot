OT=/lustre/scratch115/realdata/mdt3/projects/otcoregen

cd $OT/jeremys/datasets/GWAS

# Alzheimer's disease
F=$OT/AD_PD_finemap/summary_stats/AD.IGAP1_GWAX.meta.gz
#gzhead 1 $F
#CHR BP SNP A1 A2 GWAS_BETA GWAS_SE GWAS_P GWAX_BETA GWAX_SE GWAX_P META_BETA META_SE META_P I2 HET_P DIRECT

time (echo -e "SNP\tCHR\tBP\tA1\tMAF\tMETA_P\tMETA_BETA\tOR\tlog_OR\tMETA_SE\tz_score\ttrait\tPMID\tused_file"; 
 zcat $F | sed '1d' \
 | perl -ane 'print join("\t", $F[2],$F[0],$F[1],$F[3],"",$F[13],$F[11],"","",$F[12],"","AD","","")."\n";' \
 | sort -k2,2 -k3,3n) \
 | bgzip > Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.txt.gz

tabix -f -s 2 -b 3 -e 3 -S 1 Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.txt.gz


# Get a file of all SNPs below a P value threshold
(gzhead 1 Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.txt.gz; \
 zcat Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.txt.gz | awk 'BEGIN{FS="\t"} {if ($6 < 1e-5) print}') \
 | bgzip > Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.top_hits.txt.gz


################################################################################
# Parkinson's disease
F=$OT/AD_PD_finemap/summary_stats/pd-gwas+gwax.summarystats_2e6-signalWindows_20180330.txt.gz
#gzhead 1 $F
#SNP	CHR	BP	GENPOS	ALLELE1	ALLELE0	A1FREQ	INFO	CHISQ_LINREG	P_LINREG	BETA	SE	CHISQ_BOLT_LMM_INF	P_BOLT_LMM_INF	CHISQ_BOLT_LMM	P_BOLT_LMM

time (echo -e "SNP\tCHR\tBP\tALLELE1\tMAF\tMETA_P\tBETA\tOR\tlog_OR\tSE\tCHISQ_LINREG\ttrait\tPMID\tused_file"; 
 zcat $F | sed '1d' \
 | perl -ane 'print join("\t", $F[2],@F[0..1],$F[3],"",$F[7],$F[5],"","",$F[6],"","PD","","")."\n";' \
 | sort -k2,2 -k3,3n) \
 | bgzip > Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.sorted.txt.gz

tabix -f -s 2 -b 3 -e 3 -S 1 Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.sorted.txt.gz

# Get a file of all SNPs below a P value threshold
(gzhead 1 Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.sorted.txt.gz; \
 zcat Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.sorted.txt.gz | awk 'BEGIN{FS="\t"} {if ($6 < 1e-5) print}') \
 | bgzip > Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.top_hits.1e-5.txt.gz

(gzhead 1 Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.sorted.txt.gz; \
 zcat Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.sorted.txt.gz | awk 'BEGIN{FS="\t"} {if ($6 < 5e-8) print}') \
 | bgzip > Parkinsons_disease_Nguyen_2017_UKBB_GWAX_meta.top_hits.txt.gz


################################################################################
# Parkinson's disease OLD VERSION
F=$OT/AD_PD_finemap/summary_stats/PD.proxy.bgen.stats.gz
#gzhead 1 $F
#SNP	CHR	BP	GENPOS	ALLELE1	ALLELE0	A1FREQ	INFO	CHISQ_LINREG	P_LINREG	BETA	SE	CHISQ_BOLT_LMM_INF	P_BOLT_LMM_INF	CHISQ_BOLT_LMM	P_BOLT_LMM

time (echo -e "SNP\tCHR\tBP\tALLELE1\tMAF\tP_BOLT_LMM\tBETA\tOR\tlog_OR\tSE\tCHISQ_LINREG\ttrait\tPMID\tused_file"; 
 zcat $F | sed '1d' \
 | perl -ane 'print join("\t", @F[0..2],$F[4],"",$F[15],$F[10],"","",$F[11],$F[8],"PD","","")."\n";' \
 | sort -k2,2 -k3,3n) \
 | bgzip > Parkinsons_disease_Nguyen_2017_UKBB_GWAX.sorted.txt.gz

tabix -f -s 2 -b 3 -e 3 -S 1 Parkinsons_disease_Nguyen_2017_UKBB_GWAX.sorted.txt.gz

# Get a file of all SNPs below a P value threshold
(gzhead 1 Parkinsons_disease_Nguyen_2017_UKBB_GWAX.sorted.txt.gz; \
 zcat Parkinsons_disease_Nguyen_2017_UKBB_GWAX.sorted.txt.gz | awk 'BEGIN{FS="\t"} {if ($6 < 1e-5) print}') \
 | bgzip > Parkinsons_disease_Nguyen_2017_UKBB_GWAX.top_hits.1e-5.txt.gz

(gzhead 1 Parkinsons_disease_Nguyen_2017_UKBB_GWAX.sorted.txt.gz; \
 zcat Parkinsons_disease_Nguyen_2017_UKBB_GWAX.sorted.txt.gz | awk 'BEGIN{FS="\t"} {if ($6 < 5e-8) print}') \
 | bgzip > Parkinsons_disease_Nguyen_2017_UKBB_GWAX.top_hits.txt.gz


