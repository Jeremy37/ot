#!/usr/bin/env Rscript
## This script looks at the credible set files from the AD v2 summary stats
## and evaluates the probabilities of SNPs/indels that were not included in
## FINEMAP (since they weren't in the IGAP1 study).
library(tidyverse)

credset_input_file = "/Users/jeremys/work/opentargets/AD_finemap/annotated/AD.credible_sets.merged.tsv"


main = function() {
  df = readr::read_tsv(credset_input_file) %>%
    dplyr::mutate(F = maf, N = 10000, segment = locus_name)
  df$F[is.na(df$F)] <- 0.20
  
  # First determine 95% credible sets using FINEMAP probs
  segIndices <- getUniqueIndices(df$locus_name)
  
  df$P_Val = df$meta_p
  df$logBF = getBFFromPval(df)
  df$wtccc_meta_ppa = getSegmentNaivePPAs(df)
  
  df$P_Val = df$GWAX_p
  df$logBF = getBFFromPval(df)
  df$wtccc_gwax_ppa = getSegmentNaivePPAs(df)
  
  df$is_indel = nchar(df$UKBB_A0) != nchar(df$UKBB_A1)

  sum(df$wtccc_meta_ppa, na.rm = T)
  sum(df$wtccc_gwax_ppa, na.rm = T)
  sum(df %>% dplyr::filter(is.na(wtccc_meta_ppa)) %>% .$wtccc_gwax_ppa)
  sum(df %>% dplyr::filter(is_indel) %>% .$wtccc_gwax_ppa)
  View(df[is.na(df$wtccc_meta_ppa),])
  View(df %>% dplyr::filter(wtccc_meta_ppa > 0.001 | wtccc_gwax_ppa > 0.001) %>% arrange(locus_name, -is_indel, -wtccc_gwax_ppa))
  x = df %>% dplyr::group_by(locus_name) %>% dplyr::summarise(minp = min(meta_p, na.rm = T))
  
  pdf(file="Missing_snp_probability_sum.pdf", width=6, height = 4)
  indel_ppa.df = df %>% dplyr::filter(is_indel) %>%
    dplyr::group_by(locus_name) %>%
    dplyr::summarise(indel_prob_sum = sum(wtccc_gwax_ppa)) %>%
    arrange(-indel_prob_sum)
  indel_ppa.df$locus_name = factor(indel_ppa.df$locus_name, levels = indel_ppa.df$locus_name)
  ggplot(indel_ppa.df, aes(x=locus_name, y=indel_prob_sum)) + geom_bar(stat="identity") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Causal probability of indels (not in meta-analysis)")
  
  ukbbsnp_ppa.df = df %>% dplyr::filter(is.na(meta_p)) %>%
    dplyr::group_by(locus_name) %>%
    dplyr::summarise(missing_snp_prob_sum = sum(wtccc_gwax_ppa)) %>%
    arrange(-missing_snp_prob_sum)
  ukbbsnp_ppa.df$locus_name = factor(ukbbsnp_ppa.df$locus_name, levels = ukbbsnp_ppa.df$locus_name)
  ggplot(ukbbsnp_ppa.df, aes(x=locus_name, y=missing_snp_prob_sum)) + geom_bar(stat="identity") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Causal prob of all variants not in meta-analysis")
  dev.off()
  
  sum(df$wtccc_prob, na.rm = T)
  sum(df$finemap_prob, na.rm = T)
  x = df %>% dplyr::group_by(locus_name) %>% dplyr::summarise(sumprob = sum(finemap_prob, na.rm = T))
  
  
  # Initialize new columns we're going to add
  # newCols <- c("finemap_cumprob", "finemap_cs95")
  # df[,newCols] <- list(NA_real_, NA_real_)
  # # For each region, determine the SNPs in the FINEMAP credible set
  # for (i in 1:length(segIndices)) {
  #   segID <- as.character(segIndices[[i]][[1]])
  #   startIndex <- segIndices[[i]][[2]]
  #   endIndex <- segIndices[[i]][[3]]
  #   write(segID, stderr())
  #   df.seg <- df[startIndex:endIndex,]
  #   
  #   df.seg$finemap_cumprob = cumsum(df.seg$finemap_prob)
  #   index = which(df.seg$finemap_cumprob >= 0.95)[1]
  #   if (is.na(index)) {
  #     index = nrow(df.seg)
  #   }
  #   df.seg$finemap_cs95 = c(rep(1, index), rep(0, nrow(df.seg) - index))
  #   df[startIndex:endIndex,newCols] <- df.seg[,newCols]
  # }
}


#############################################################################

# First determine the indices of the genes in the table
getUniqueIndices = function(vec)
{
  if (class(vec) == "factor") {
    # For some reason factors are INCREDIBLY slow if used in the code below
    vec <- as.character(vec)
  }
  indicesList <- list()
  if (length(vec) == 0) {
    return(indicesList)
  }
  
  lastVal = vec[1]
  lastIndex = 1
  for (i in 1:length(vec)) {
    if (vec[i] != lastVal) {
      indicesList <- c(indicesList, list(list(lastVal, as.integer(lastIndex), as.integer(i-1))))
      lastVal <- vec[i]
      lastIndex <- i
    }
  }
  indicesList <- c(indicesList, list(list(lastVal, as.integer(lastIndex), as.integer(i))))
  return(indicesList)
}

getBFFromPval = function(df) {
  df[,'Z'] <- pToZStat(df$P_Val)
  df$F[is.na(df$F)] <- 0.10
  df$N[is.na(df$N)] <- 10000
  apply(df[,c('Z','F','N')], 1, function(d) calcLogBF(d['Z'],d['F'],d['N']))
}


pToZStat = function(pVals) {
  sqrt(qchisq(p=pVals, df=1, lower.tail=F))
}

calcLogBF = function(Z, f, N) {
  WW <- 0.1
  V <- approx_v(f, N)
  r <- WW / (V + WW)
  toreturn <- log( sqrt(1-r) ) + (Z*Z*r / 2)
  toreturn
}

approx_v = function(f, N) {
  1 / (2*f*(1-f) * N)
}

getSegmentNaivePPAs = function(df)
{
  segIndices <- getUniqueIndices(df$segment)
  ppas <- vector(mode="double")
  # For each gene/segment...
  for (i in 1:length(segIndices)) {
    startIndex <- segIndices[[i]][[2]]
    endIndex <- segIndices[[i]][[3]]
    ppas[startIndex:endIndex] <- getNaivePPA(df$logBF[startIndex:endIndex])
  }
  ppas
}

getNaivePPA = function(vecLogBF)
{
  logsegbfNaive <- -1000
  for (i in 1:length(vecLogBF)) {
    if (!is.na(vecLogBF[i])) {
      logsegbfNaive <- sumlog(logsegbfNaive, vecLogBF[i])
    }
  }
  vecPPA <- vecLogBF - logsegbfNaive
  exp(vecPPA)
}

sumlog = function(logx, logy)
{
  if (logx > logy) return(logx + log(1 + exp(logy-logx)))
  else return(logy + log(1 + exp(logx-logy)))
}


#############################################################################

main()

