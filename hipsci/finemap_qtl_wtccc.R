#!/usr/bin/env Rscript
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
library(annotables)

DEFAULT_AF = 0.10
DEFAULT_N = 500
rpkmsFile = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/reference/tissueRPKM/tissues.selected.rpkm_average.txt.gz"

# Function that uses WTCCC-style method to compute posterior probabilities of
# causality for each SNP associated with each feature. Since I don't have allele
# frequencies, I just assume all SNPs are a fixed frequency.
getPosteriorProbs = function(df) {
  df$logBF = getBFFromPval(df)
  df$prob = getNaivePPA(df$logBF)
  df
}

getBFFromPval = function(df)
{
  if (is.null(df$F)) {
    df$F = NA
  }
  if (is.null(df$N)) {
    df$N = NA
  }
  df$pvalue[df$pvalue <= 0] = 1e-300
  df[,'Z'] <- pToZStat(df$pvalue)
  df$F[is.na(df$F)] <- DEFAULT_AF
  df$N[is.na(df$N)] <- DEFAULT_N
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

getNaivePPA = function(vecLogBF)
{
  logsegbfNaive <- -1000
  for (i in 1:length(vecLogBF)) {
    logsegbfNaive <- sumlog(logsegbfNaive, vecLogBF[i])
  }
  vecPPA <- vecLogBF - logsegbfNaive
  exp(vecPPA)
}

sumlog = function(logx, logy)
{
  if (logx > logy) return(logx + log(1 + exp(logy-logx)))
  else return(logy + log(1 + exp(logx-logy)))
}



main = function() {
  option_list <- list(
    make_option(c("--input"), type="character", default=NULL, metavar="qtls.tsv.bgz", help="Path to tabix-indexed file with QTL data."),
    make_option(c("--out"), type="character", default=NULL, metavar="outpath", help="Base path for output files. Required.")
  )
  parser = OptionParser(usage = "finemap_qtl_wtccc.R --input file.tsv --out filename_base [options]",
                        option_list = option_list)
  opt <<- parse_args(parser)
  # opt = list(input = "/Users/jeremys/work/opentargets/ipsc/hipsci_finemap/SplicingLevel/full_qtl_results_22.rescaled.gene_sorted.txt.gz",
  #            out = "/Users/jeremys/work/opentargets/ipsc/hipsci_finemap/SplicingLevel/full_qtl_results_22.finemapped")
  # opt = list(input = "/Users/jeremys/work/opentargets/ipsc/hipsci_finemap/ApaLevel/full_qtl_results_13.rescaled.gene_sorted.txt.gz",
  #            out = "/Users/jeremys/work/opentargets/ipsc/hipsci_finemap/ApaLevel/full_qtl_results_13.finemapped")
  checkRequiredOpt = function(opt, parser, optName) {
    if (is.null(opt[[optName]])) {
      print_help(parser)
      cat(sprintf("ERROR: argument --%s is required\n\n", optName))
      stop()
    }
  }
  checkRequiredOpt(opt, parser, "input")
  checkRequiredOpt(opt, parser, "out")
  
  qtl_colnames = c("gene", "pos", "chr", "beta", "se", "ref", "alt", "pvalue")
  qtl.df = readr::read_tsv(opt$input, col_names = T, col_types = "cicddccd") %>%
    rename(pos = snp_pos)
  # Remove ending version (e.g. ".4" in ENSG00123.4)
  qtl.df$gene = gsub("\\.[\\d]+", "", qtl.df$gene, perl=T)
  
  # Filter QTL SNP table to remove genes with best p value worse than 1e-5
  qtl_minp.df = qtl.df %>%
    group_by(gene) %>%
    summarise(gene_minp = min(pvalue, na.rm=T), gene_max_log10p = max(-log10(pvalue), na.rm=T))
  
  qtl.df = qtl.df %>%
    left_join(qtl_minp.df, by="gene") %>%
    filter(gene_max_log10p > 5)
  
  # Get posterior probabilities for SNPs, one gene at a time, using "by"
  res = by(qtl.df, qtl.df$gene, getPosteriorProbs)
  qtl.pp.df = bind_rows(lapply(1:length(res), function(i) res[[i]]))
  
  # Extract one chromosome of data from SpliceAI genome-wide scores
  chr = qtl.pp.df$chr[1]
  spliceaiFname = sprintf("spliceai.merge.%s.tsv.gz", chr)
  if (!file.exists(spliceaiFname)) {
    cmd = sprintf("tabix /lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/datasets/SpliceAI/spliceai.merge.tsv.bgz %s | bgzip > %s", chr, spliceaiFname)
    system(cmd, intern=T)
  }
  
  spliceaiColnames = c("chr", "pos", "ref", "alt", "symbol", "strand", "type", "dist", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL")
  spliceai.df = readr::read_tsv(spliceaiFname, col_names = spliceaiColnames, col_types = "cicccccidddddddd") %>%
    mutate(max_ds = pmax(DS_AG, DS_AL, DS_DG, DS_DL)) %>%
    select(chr, pos, ref, alt, symbol, strand, type, max_ds, everything())
  
  # Add a custom score for each SNP to find good causal candidates
  qtl.new.df = qtl.pp.df %>%
    filter(prob > 1e-3) %>%
    left_join(spliceai.df, by=c("chr", "pos", "ref", "alt")) %>%
    mutate(pval_score = pmin(1, gene_max_log10p / 50)^0.333,
           score = prob * pval_score,
           score_splice = prob * pval_score * max_ds)
  
  max_scores.df = qtl.new.df %>%
    group_by(gene) %>%
    summarise(gene_max_prob = max(prob, na.rm=T),
              gene_max_ds = max(max_ds, na.rm=T),
              gene_max_score = max(score, na.rm=T),
              gene_max_score_splice = max(score_splice, na.rm=T)) %>%
    mutate(gene_max_prob = ifelse(gene_max_prob < 0, NA, gene_max_prob),
           gene_max_ds = ifelse(gene_max_ds < 0, NA, gene_max_ds),
           gene_max_score = ifelse(gene_max_score < 0, NA, gene_max_score),
           gene_max_score_splice = ifelse(gene_max_score_splice < 0, NA, gene_max_score_splice))
  
  
  qtl.out.all.df = qtl.new.df %>%
    left_join(max_scores.df, by="gene") %>%
    arrange(desc(gene_max_score), desc(prob)) %>%
    select(-starts_with("gene_"), starts_with("gene_"))
  
  gzf = gzfile(sprintf("%s.all.tsv.gz", opt$out), "w")
  write.table(qtl.out.all.df, file = gzf,
              col.names=T, row.names=F, sep="\t", quote=F, na="")
  close(gzf)
  
  qtl.out.summary.df = qtl.out.all.df %>%
    group_by(gene) %>%
    summarise(gene_max_prob = first(gene_max_prob), gene_max_ds = first(gene_max_ds), gene_max_log10p = first(gene_max_log10p),
              gene_max_score = first(gene_max_score), gene_max_score_splice = first(gene_max_score_splice),
              prob = first(prob), chr = first(chr), pos = first(pos), ref = first(ref), alt = first(alt),
              symbol = first(symbol), strand = first(strand), type = first(type),
              DS_AG = first(DS_AG), DS_AL = first(DS_AL), DS_DG = first(DS_DG), DS_DL = first(DS_DL),
              DP_AG = first(DP_AG), DP_AL = first(DP_AL), DP_DG = first(DP_DG), DP_DL = first(DP_DL)) %>%
    arrange(desc(gene_max_score))
  
  # Add annotation of gene expression to summary table
  # First get Ensgene ID if it isn't present, using the existing symbol annotation
  grch38_dedup = grch38 %>% filter(!duplicated(ensgene))
  if (!grepl("^ENSG", qtl.out.summary.df$gene[1])) {
    qtl.out.summary.df$splice_event = qtl.out.summary.df$gene
    qtl.out.summary.df = qtl.out.summary.df %>%
      select(-gene) %>%
      left_join(grch38_dedup %>% select(symbol, gene=ensgene)) %>%
      select(splice_event, gene, everything())
  } else {
    qtl.out.summary.df = qtl.out.summary.df %>%
      select(-symbol, -strand) %>%
      left_join(grch38_dedup %>% select(symbol, gene=ensgene, strand), by="gene") %>%
      select(gene:alt, symbol, strand, everything())
  }
  
  rpkm.df = readr::read_tsv(rpkmsFile) %>%
    select(gene_id, IPSC, iNeuron_d11, neuron, macrophage_naive, microglia_gaffney, `gtex.Brain - Hippocampus`)
  qtl.out.summary.df = qtl.out.summary.df %>%
    left_join(rpkm.df, by = c("gene" = "gene_id"))
  
  write.table(qtl.out.summary.df, file = sprintf("%s.summary.tsv", opt$out),
              col.names=T, row.names=F, sep="\t", quote=F, na="")
}


main()


