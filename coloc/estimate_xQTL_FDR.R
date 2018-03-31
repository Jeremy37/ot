library(qvalue)

# in bash estimate the number of SNPs tested per feature:
# zcat haQTLs_all.txt.gz | cut -f 4 | head -n 1000000 | sort | uniq | wc -l
# 153
# SNPs tested per ha peak ~ 1e6 / 153 =~ 5000

df = readr::read_tsv(file.path(root, "reference", "xQTL", "xQTL_haQTL.gene_minp.minp.txt"), col_names = "pval")
df$pcorr = df$pval * 5000
df$pcorr[df$pcorr > 1] = 1
qv = qvalue::qvalue(df$pcorr)
df$q = qv$qvalues

# in bash estimate the number of SNPs tested per feature:
# zcat mQTLs_all.txt.gz | cut -f 4 | head -n 1000000 | sort | uniq | wc -l
# 1856 
# SNPs tested per me site ~ 1e6 / 1856 =~ 500
df = readr::read_tsv(file.path(root, "reference", "xQTL", "xQTL_mQTL.gene_minp.minp.txt"), col_names = "pval")
df$pcorr = df$pval * 500
df$pcorr[df$pcorr > 1] = 1
qv = qvalue::qvalue(df$pcorr)
df$q = qv$qvalues
