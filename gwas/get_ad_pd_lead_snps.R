#!/usr/bin/env Rscript
library(tidyverse)

root = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys"
#root = "/Users/jeremys/work/opentargets"

####################### AD
finemap.df = readr::read_tsv(file.path(root, "gwas/AD/finemap/AD.credible_sets.set"))

leadsnp.df = finemap.df %>% group_by(locus) %>%
  dplyr::summarise(first(SNP[p == min(p)]), min(p)) %>%
  dplyr::rename(SNP = `first(SNP[p == min(p)])`, p = `min(p)`) %>%
  dplyr::left_join(finemap.df %>% dplyr::select(SNP, Chr, pos, prob))

ad.loci = readr::read_tsv(file.path(root, "gwas/AD/AD.locusnames.txt"))
leadsnp.df = leadsnp.df %>%
  dplyr::left_join(ad.loci %>% dplyr::rename(locus = locus_key)) %>%
  dplyr::select(locus, locus_name, everything()) %>%
  dplyr::arrange(Chr, pos)

write.table(leadsnp.df, file=file.path(root, "gwas/AD/AD.finemap.leadSNPs.txt"), sep="\t", quote=F, row.names=F)

# Write a file of "loci" for getting LD
leadsnp.df = leadsnp.df %>%
  dplyr::mutate(region = paste0(Chr, ":", pos - 249999, "..", pos + 249999))
write.table(leadsnp.df$region, file=file.path(root, "gwas/AD/LD/AD.regions.txt"), sep="\t", quote=F, row.names=F, col.names=F)


####################### PD
finemap.df = readr::read_tsv(file.path(root, "gwas/PD/finemap/pd_gwas_gwax.combined.set"))

leadsnp.df = finemap.df %>% group_by(locus) %>%
  dplyr::summarise(first(SNP[p == min(p)]), min(p)) %>%
  dplyr::rename(SNP = `first(SNP[p == min(p)])`, p = `min(p)`) %>%
  dplyr::left_join(finemap.df %>% dplyr::select(SNP, Chr, pos, prob))

write.table(leadsnp.df, file=file.path(root, "gwas/PD/PD.finemap.leadSNPs.txt"), sep="\t", quote=F, row.names=F)

# Write a file of "loci" for getting LD
leadsnp.df = leadsnp.df %>%
  dplyr::mutate(region = paste0(Chr, ":", pos - 249999, "..", pos + 249999))
write.table(leadsnp.df$region, file=file.path(root, "gwas/PD/LD/PD.regions.txt"), sep="\t", quote=F, row.names=F, col.names=F)

