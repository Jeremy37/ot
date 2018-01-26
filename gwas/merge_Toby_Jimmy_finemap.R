#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
tobyFile = args[1]
finemapSetFile = args[2]
finemapSnpFile = args[3]

toby.df = readr::read_tsv(tobyFile) %>%
  dplyr::rename(Chr = chrom,
                toby.pval = pval,
                toby.pval_cond = pval_cond,
                toby.posterior = posterior) %>%
  dplyr::mutate(chrpos = paste0(Chr, "_", pos))

finemapSet.df = readr::read_tsv(finemapSetFile) %>%
  dplyr::rename(snp = SNP,
                finemap.p = p,
                finemap.setprob = prob,
                finemap.cs95 = cs95,
                finemap.cs99 = cs99) %>%
  dplyr::mutate(chrpos = paste0(Chr, "_", pos))

finemapSnp.df = readr::read_tsv(finemapSnpFile) %>%
  dplyr::rename(finemap.snpprob = snp_prob) %>%
  dplyr::select(snp, finemap.snpprob)

finemapSet.df = finemapSet.df %>% dplyr::left_join(finemapSnp.df, by="snp")

joined.df = toby.df %>% dplyr::left_join(finemapSet.df %>%
                                           dplyr::select(chrpos, snp, finemap.p, finemap.snpprob, finemap.setprob, finemap.cs95, finemap.cs99),
                                         by="chrpos")

finemap.flt.df = finemapSet.df %>% filter(finemap.snpprob > 0.001,
                                       !(chrpos %in% joined.df$chrpos)) %>%
  dplyr::select(Chr, pos, snp, finemap.p, finemap.snpprob, finemap.setprob, finemap.cs95, finemap.cs99) %>%
  dplyr::mutate(signal=NA, ref=NA, alt=NA, toby.pval=NA, toby.pval_cond=NA,
                toby.posterior=NA, consequences=NA, posterior_consequences=NA, region=NA)

merged.df = rbind(joined.df %>% dplyr::select(-chrpos), finemap.flt.df) %>%
  dplyr::arrange(Chr, pos, -finemap.snpprob)

# Wrote this table originally to get the boundaries of the loci, used below
#write.table(merged.df, file="toby.jimmy.finemap.tmp.txt", quote=F, sep="\t", col.names=T, row.names=F)
#View(merged.df)

# In this merged set, we have SNPs where a "signal name" isn't assigned.
# Let's manually fix that up, as it makes ordering by signal and then
# by posterior probability possible.
loci.df = readr::read_tsv("AD.loci.txt")
for (i in 1:nrow(loci.df)) {
  chr = loci.df[i,]$chr
  start = loci.df[i,]$start
  end = loci.df[i,]$end
  signal = loci.df[i,]$signal
  merged.df[merged.df$Chr == loci.df[i,]$chr & merged.df$pos >= start & merged.df$pos <= end,]$signal = signal
}
merged.df %<>% arrange(Chr, signal, -finemap.snpprob)

write.table(merged.df, "toby.jimmy.finemap.merged.txt", quote=F, sep="\t", col.names=T, row.names=F)


