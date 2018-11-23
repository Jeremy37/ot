#!/usr/bin/env Rscript
library(tidyverse)
library(optparse)
library(coloc)
library(purrrlyr)
library(GenomicRanges)
#source("colocTools.R")

gene_id_map = NULL

#Parse command-line options
option_list <- list(
  make_option(c("--qtl_name"), type="character", default=NULL,
              help="Name of QTLs, used for output."),
  make_option(c("--qtl_nominal"), type="character", default=NULL,
              help="Path to the QTL nominal p values file"),
  make_option(c("--qtl_variant_info"), type="character", default=NULL,
              help="Path to a file with QTL SNP information, i.e. rsID, chr, pos, MAF, which will replace similar columns in the QTL summary files. This is useful for switching between GRCh37/38 coords."),
  make_option(c("--qtl_signals"), type="character", default=NULL,
              help="Path to the QTL min p value per gene file"),
  make_option(c("--qtl_samplesize"), type="integer", default=NULL,
              help="Integer sample size of the QTL study"),
  make_option(c("--qtl_list"), type="character", default=NULL,
              help="Path to a file listing QTL dataset to use for coloc, with each line giving qtl_name, qtl_nominal file, qtl_signals file, and qtl_samplesize"),
  make_option(c("--gwas_name"), type="character", default=NULL,
              help="Name of the GWAS trait"),
  make_option(c("--gwas_file"), type="character", default=NULL,
              help="Path to GWAS summary stats directory."),
  make_option(c("--gwas_signals"), type="character", default=NULL,
              help="Threshold of GWAS p value to use locus for colocalisation"),
  make_option(c("--window"), type="numeric", default=NULL,
              help="Size of the cis window around both the QTL and GWAS lead SNP where summary statistics should be used."),
  make_option(c("--overlap_dist"), type="numeric", default="5e5",
              help="Max distance between QTL lead SNP and GWAS lead SNP."),
  make_option(c("--plot_style"), type="character", default=NULL,
              help="Type of plots to generate as a PDF output - either 'overlay' or 'line' or 'all'"),
  make_option(c("--plot_threshold"), type="numeric", default=0, action="store_true",
              help="Minimum H4 result from coloc for a plot to be generated"),
  make_option(c("--outdir"), type="character", default=NULL,
              help="Path to the output directory."),
  make_option(c("--p1"), type="numeric", default=1e-4,
              help="Prior probability of SNP affecting expression"),
  make_option(c("--p2"), type="numeric", default=1e-4,
              help="Prior probability of SNP affecting GWAS trait"),
  make_option(c("--p12"), type="numeric", default=1e-5,
              help="Prior probability of SNP affecting both traits")
)

importQTLSignals = function(qtlSignalsPath) {
  # This reads a table of significant QTLs, all of which will be tested for coloc
  # with any nearby gwas signal. SNP chr and pos will be mapped using the QTL
  # variant info passed in
  table = readr::read_tsv(qtlSignalsPath)
  expectedCols = c("feature", "chr", "pos", "rsid")
  if (!hasCols(table, expectedCols)) {
    stop(sprintf("In table %s, expected columns with header names %s.",
                 qtlSignalsPath, paste(expectedCols, collapse = ", ")))
  }
  return(table %>% dplyr::select(feature, chr, pos, rsid))
}

importVariantInfo = function(variantInfoPath) {
  table = readr::read_tsv(variantInfoPath, col_names = T, col_types = "cicd")
  expectedCols = c("chr", "pos", "rsid", "MAF")
  if (!hasCols(table, expectedCols)) {
    stop(sprintf("In table %s, expected columns with header names %s.",
                 variantInfoPath, paste(expectedCols, collapse = ", ")))
  }
  table$chr = as.character(table$chr)
  return(table %>% dplyr::select(chr, pos, rsid, MAF))
}

importGWASSignals = function(gwasSignalsPath) {
  table = readr::read_tsv(gwasSignalsPath)
  expectedCols = c("chr", "pos", "rsid")
  if (!hasCols(table, expectedCols)) {
    stop(sprintf("In table %s, expected columns with header names %s.",
                 gwasSignalsPath, paste(expectedCols, collapse = ", ")))
  }
  table$chr = as.character(table$chr)
  return(table %>% dplyr::select(chr, pos, rsid))
}


tabixFetchQTLs = function(phenotype_ranges, tabix_file, qtl_colnames, qtl_coltypes) {
  #Assertions about input
  assertthat::assert_that(class(phenotype_ranges) == "GRanges")
  assertthat::assert_that(assertthat::has_name(GenomicRanges::elementMetadata(phenotype_ranges), "feature"))
  qtl_colnames = c("feature", "chr", "pos", "rsid", "p_nominal")
  qtl_coltypes = c("ccicd")
  result = list()
  for (i in seq_along(phenotype_ranges)) {
    selected_feature = phenotype_ranges[i]$feature
    print(i)
    tabix_table = scanTabixDataFrame(tabix_file, phenotype_ranges[i], 
                                     col_names = qtl_colnames, col_types = qtl_coltypes)[[1]] %>%
      dplyr::filter(feature == selected_feature)
    
    #Add additional columns
    result[[selected_feature]] = tabix_table
  }
  return(result[[1]])
}


main = function() {
  opt <- parse_args(OptionParser(option_list=option_list))
  print(opt)
  # root = "/Users/jeremys/work/opentargets"
  # opt = list(qtl_name = "sens_neur.eqtl.500k", qtl_samplesize = 97, gwas_name = "AD",
  #            qtl_nominal = file.path(root, "coloc/qtl_data/sensory_neuron.fastqtl.cis500k.nominals.GRCh38.for_coloc.txt.gz"),
  #            qtl_variant_info = file.path(root, "coloc/qtl_data/sensory_neuron.variant_info.for_coloc.GRCh37.gz"),
  #            qtl_signals = file.path(root, "coloc/qtl_data/sensory_neuron.fastqtl.cis500k.qtl_signals.txt"),
  #            gwas_file = file.path(root, "datasets/GWAS/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.txt.gz"),
  #            gwas_signals = file.path(root, "coloc/AD.signals.txt"),
  #            window = 2e5, overlap_dist = 5e5,
  #            plot_style = "all", plot_threshold = 0, outdir = file.path(root, "coloc/new"),
  #            p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)

  if (!is.null(opt$plot_style)) {
    if (!(identical(opt$plot_style, "overlay") | identical(opt$plot_style, "line") | identical(opt$plot_style, "all"))) {
      stop(sprintf("Error: unrecognized value for plot_style parameter: %s", opt$plot_style))
    }
  }
  
  gwas_signals = importGWASSignals(opt$gwas_signals) %>%
    dplyr::transmute(chr, gwas_pos = pos, gwas_rsid=rsid) %>%
    dplyr::filter(!is.na(gwas_pos))
  
  qtl_signals = importQTLSignals(opt$qtl_signals) %>%
    dplyr::rename(qtl_pos_orig = pos) %>% dplyr::select(-chr)
  
  qtl_variant_info = importVariantInfo(opt$qtl_variant_info)
  new_qtl_signals = qtl_signals %>%
    dplyr::left_join(qtl_variant_info, by="rsid") %>%
    dplyr::select(-MAF) %>%
    dplyr::rename(qtl_pos = pos, qtl_rsid = rsid)
  if (sum(is.na(new_qtl_signals$qtl_pos)) > 0) {
    warning(sprintf("WARNING: %d of %d QTL signals did not have positions mapped in the QTL variant info", sum(is.na(new_qtl_signals$qtl_pos)), nrow(qtl_signals)))
    new_qtl_signals = new_qtl_signals %>% dplyr::filter(!is.na(qtl_pos))
  }
  qtl_signals = new_qtl_signals
  
  colocPairs = getColocPairs(qtl_signals, gwas_signals, opt$overlap_dist)
  
  # HACK: load table mapping ENSG ID to HGNC symbol, so we can use symbol in plots
  #gene_id_map = readr::read_tsv("/Users/jeremys/work/opentargets/reference/hgnc.ensembl.map.txt") %>%
  gene_id_map = readr::read_tsv("/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/reference/hgnc.ensembl.map.txt") %>%
    dplyr::filter(symbol_type == "symbol")
  
  ###################### Run the coloc for each QTL-GWAS pair! ##################
  write(sprintf("Testing coloc for %d QTL-GWAS signal pairs.", nrow(colocPairs)), stderr())
  
  colocs_list = vector("list", nrow(colocPairs))
  for (i in 1:nrow(colocPairs)) {
      colocs_list[[i]] = colocMolecularQTLs(colocPairs[i,],
                             qtl_summary_path = opt$qtl_nominal,
                             gwas_summary_path = opt$gwas_file,
                             qtl_variant_info = qtl_variant_info,
                             N_qtl = opt$qtl_samplesize,
                             cis_dist = opt$window,
                             p1 = opt$p1, p2 = opt$p2, p12 = opt$p12,
                             plot_style = opt$plot_style, plot_H4_threshold = opt$plot_threshold)
  }
  
  coloc_hits = bind_cols(colocPairs, bind_rows(lapply(colocs_list, FUN = function(x) x$summary)))
  
    # coloc_res = purrrlyr::by_row(colocPairs,
    #                            ~colocMolecularQTLs(.,
    #                                                qtl_summary_path = opt$qtl_nominal,
    #                                                gwas_summary_path = opt$gwas_file, 
    #                                                qtl_variant_info = qtl_variant_info,
    #                                                N_qtl = opt$qtl_samplesize,
    #                                                cis_dist = opt$window,
    #                                                p1 = opt$p1, p2 = opt$p2, p12 = opt$p12)$summary,
    #                            .collate = "rows")
  # coloc_hits = coloc_res %>% dplyr::arrange(chr, gwas_lead)

  coloc_hits$feature = gsub("\\.[\\d]+", "", coloc_hits$feature, perl=T)
  
  coloc_hits = coloc_hits %>%
    dplyr::rename(PP.H0=PP.H0.abf, PP.H1=PP.H1.abf, PP.H2=PP.H2.abf, PP.H3=PP.H3.abf, PP.H4=PP.H4.abf) %>%
    dplyr::select(feature, nsnps, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, qtl_pval, gwas_pval,
                  qtl_lead, gwas_lead, gwas_signal_rsid = gwas_rsid, chr, qtl_pos, gwas_pos)
  
  coloc_output = file.path(opt$outdir, paste("coloc", opt$gwas_name, opt$qtl_name, opt$overlap_dist, "txt", sep = "."))
  write.table(coloc_hits, coloc_output, sep = "\t", quote = FALSE, row.names = FALSE)
  
  if (!is.null(opt$plot_style)) {
    coloc_plots = lapply(colocs_list, FUN = function(x) x$plots)
    pdf(file.path(opt$outdir, sprintf("coloc.%s.%s.%s.%s.plots.pdf", opt$gwas_name, opt$qtl_name, opt$overlap_dist, opt$plot_style)),
        width=8, height=5.5)
    for (colocPlot in coloc_plots) {
      if (!is.null(colocPlot)) {
        print(colocPlot)
      }
    }
    dev.off()
  }
}


###############################################################################
# colocTools

hasCols = function(table, cols) {
  all(sapply(cols, function(x) {x %in% colnames(table)}))
}

#' Determine the set of QTL-GWAS signal pairs to test for colocalisation, based
#' on an overlap window.
#'
#' @param qtl_signals List of data frames with QTL lead pvalues. Each data frame must contain
#' gene_id, snp_id and either p_val or FDR columns, and should not contain other columns.
#' @param gwas_signals Prefix of the GWAS summarystats file
#' @param overlap_dist Max distance between GWAS and QTL variants.
#'
#' @return List of data.frames with phenotype_ids and snp_ids to be tested with coloc.
#' @export
getColocPairs <- function(qtl_signals, gwas_signals, overlap_dist = 5e5) {
  coloc_list = qtl_signals %>%
    dplyr::left_join(gwas_signals, by="chr") %>%
    dplyr::mutate(qtl_gwas_distance = abs(gwas_pos - qtl_pos)) %>%
    dplyr::filter(qtl_gwas_distance < overlap_dist)
  coloc_list
}

colocMolecularQTLs <- function(coloc_pair, qtl_summary_path, gwas_summary_path,
                               qtl_variant_info = NULL, N_qtl = 84, cis_dist = 1e5,
                               p1 = 1e-04, p2 = 1e-04, p12 = 1e-05,
                               plot_style = NULL, plot_H4_threshold = 0) {
  
  expectedCols = c("feature", "chr", "qtl_pos_orig", "gwas_pos", "qtl_rsid", "gwas_rsid")
  assertthat::assert_that(hasCols(coloc_pair, expectedCols))
  assertthat::assert_that(nrow(coloc_pair) == 1)
  
  #Print for debugging
  write(sprintf("Testing coloc for GWAS locus %s and QTL: %s, %s:%d, %s",
                coloc_pair$gwas_rsid, coloc_pair$feature, coloc_pair$chr, coloc_pair$qtl_pos, coloc_pair$qtl_rsid),
        stderr())
  
  coloc_pair = as.list(coloc_pair[1,])
  minPos = min(coloc_pair$qtl_pos, coloc_pair$gwas_pos)
  maxPos = max(coloc_pair$qtl_pos, coloc_pair$gwas_pos)
  diffPos = coloc_pair$gwas_pos - coloc_pair$qtl_pos
  
  result = tryCatch({
    # Make GRanges objects with a distance of cis_dist upstream and downstream
    # of the GWAS and QTL variants. The QTL range should be relative to the coords
    # of the original QTL signal, before coords were replaced
    qtlStart = coloc_pair$qtl_pos_orig - cis_dist
    qtlEnd = coloc_pair$qtl_pos_orig + cis_dist
    if (diffPos > 0) {
      qtlEnd = qtlEnd + diffPos
    } else {
      qtlStart = qtlStart + diffPos
    }
    qtlRange = GRanges(feature = coloc_pair$feature, seqnames = coloc_pair$chr, ranges = IRanges(start = qtlStart, end = qtlEnd), strand = "*")
    gwasRange = GRanges(seqnames = coloc_pair$chr, ranges = IRanges(start = minPos - cis_dist, end = maxPos + cis_dist), strand = "*")
    
    write(sprintf("QTL range: %s:%s:%d-%d\tGWAS range: %s:%d-%d",
                  coloc_pair$feature, coloc_pair$chr, qtlStart, qtlEnd,
                  coloc_pair$chr, minPos - cis_dist, maxPos + cis_dist),
          stderr())
    
    # Fetch QTL summary stats. The code calling colocMolecularQTLs needs to define
    # the function tabixFetchQTLs, which will handle fetching QTLs from potentially
    # differently formatted files.
    qtl_summaries = tabixFetchQTLs(qtlRange, qtl_summary_path)
    if (!is.null(qtl_variant_info)) {
      qtl_summaries = summaryReplaceCoordinates(qtl_summaries, qtl_variant_info)
    }
    
    #Fetch GWAS summary stats
    gwas_summaries = tabixFetchGWASSummary(gwasRange, gwas_summary_path)[[1]]
    
    qtl_min_index = which.min(qtl_summaries$p_nominal)
    gwas_min_index = which.min(gwas_summaries$p_nominal)
    
    #Perform coloc analysis
    coloc_res = colocQtlGWAS(qtl_summaries, gwas_summaries, N_qtl = N_qtl, p1, p2, p12)
    coloc_summary = dplyr::tbl_df(t(data.frame(coloc_res$summary))) %>%
      dplyr::mutate(qtl_pval = qtl_summaries[qtl_min_index,]$p_nominal, gwas_pval = gwas_summaries[gwas_min_index,]$p_nominal,
                    qtl_lead = qtl_summaries[qtl_min_index,]$rsid, gwas_lead = gwas_summaries[gwas_min_index,]$rsid)

    colocPlots = NULL
    if (!is.null(plot_style) & coloc_summary$PP.H4.abf >= plot_H4_threshold) {
      if (identical(plot_style, "overlay")) {
        colocPlots = list(colocOverlayPlot(qtl_summaries, gwas_summaries, coloc_pair, coloc_summary))
      } else if (identical(plot_style, "line")) {
        colocPlots = list(colocLinesPlot(qtl_summaries, gwas_summaries, coloc_pair, coloc_summary))
      } else if (identical(plot_style, "all")) {
        plot1 = colocOverlayPlot(qtl_summaries, gwas_summaries, coloc_pair, coloc_summary)
        plot2 = colocLinesPlot(qtl_summaries, gwas_summaries, coloc_pair, coloc_summary)
        colocPlots = list(plot1, plot2)
      }
    }
    #Summary list
    data_list = list(qtl = qtl_summaries, gwas = gwas_summaries)
    
    result = list(summary = coloc_summary, data = data_list, plots = colocPlots)
  }, error = function(err) {
    print(paste("ERROR:", err))
    result = list(summary = data.frame(nsnps=NA, PP.H0.abf=NA, PP.H1.abf=NA, PP.H2.abf=NA, PP.H3.abf=NA, PP.H4.abf=NA, qtl_pval=NA, gwas_pval=NA, qtl_lead=NA, gwas_lead=NA),
                  data = list(qtl = qtl_summaries, gwas = gwas_summaries), plots = NULL)
  }
  )
  return(result)
}


#' Test colocalisation between molecular QTL and GWAS summary stats
#'
#' @param qtl QTL summary stats (p_nominal, MAF, beta, snp_id)
#' @param gwas GWAS summary stats(beta, se, MAF, log_OR)
#' @param N_qtl Sample size of the QTL mapping study
#'
#' @return coloc.abf result object
#' @export
colocQtlGWAS <- function(qtl, gwas, N_qtl, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05, case_control_proportion = 0.1) {
  #Check that QTL df has all correct names
  assertthat::assert_that(assertthat::has_name(qtl, "rsid"))
  #assertthat::assert_that(assertthat::has_name(qtl, "beta"))
  assertthat::assert_that(assertthat::has_name(qtl, "MAF"))
  assertthat::assert_that(assertthat::has_name(qtl, "p_nominal"))
  
  #Check that GWAS df has all correct names
  assertthat::assert_that(assertthat::has_name(gwas, "beta"))
  assertthat::assert_that(assertthat::has_name(gwas, "se"))
  assertthat::assert_that(assertthat::has_name(gwas, "rsid"))
  assertthat::assert_that(assertthat::has_name(gwas, "log_OR"))
  
  #Count NAs for log_OR and beta
  log_OR_NA_count = length(which(is.na(gwas$log_OR)))
  beta_NA_count = length(which(is.na(gwas$beta)))
  
  #Remove GWAS SNPs with NA std error
  gwas = dplyr::filter(gwas, !is.na(se))
  
  # Subset to SNPs in common between QTL and GWAS
  qtl_in = qtl
  gwas_in = gwas
  qtl = qtl_in %>% dplyr::inner_join(gwas_in %>% dplyr::select(rsid), by="rsid")
  gwas = gwas_in %>% dplyr::inner_join(qtl_in %>% dplyr::select(rsid), by="rsid")
  write(sprintf("%d SNPs in common, of %d QTL SNPs and %d GWAS SNPs", nrow(qtl), nrow(qtl_in), nrow(gwas_in)), stderr())
  if (nrow(qtl) < 2) {
    stop("Too few SNPs in common between QTL and GWAS to run colocalisation.")
  }
  qtl_lead_rsid = qtl_in[which.min(qtl_in$p_nominal),]$rsid
  gwas_lead_rsid = gwas_in[which.min(gwas_in$p_nominal),]$rsid
  if (!qtl_lead_rsid %in% qtl$rsid) {
    write(sprintf("Warning: QTL lead SNP %s is not in set of SNPs in common.", qtl_lead_rsid), stderr())
  }
  if (!gwas_lead_rsid %in% gwas$rsid) {
    write(sprintf("Warning: GWAS lead SNP %s is not in set of SNPs in common.", gwas_lead_rsid), stderr())
  }
  
  # Make QTL dataset object
  df1 = list(pvalues = qtl$p_nominal, 
             N = N_qtl, 
             MAF = qtl$MAF, 
             type = "quant", 
             beta = NA,
             snp = qtl$rsid)
  if ("beta" %in% colnames(qtl)) {
    df1$beta = qtl$beta
  }
  
  #If beta is not specified in the GWAS then use log_OR
  if(beta_NA_count <= log_OR_NA_count){
    assertthat::assert_that(assertthat::has_name(gwas, "MAF"))
    coloc_res = coloc::coloc.abf(dataset1 = df1,
                                 dataset2 = list(beta = gwas$beta, 
                                                 varbeta = gwas$se^2, 
                                                 type = "cc",
                                                 s = case_control_proportion,
                                                 snp = gwas$rsid, 
                                                 MAF = gwas$MAF),
                                 p1 = p1, p2 = p2, p12 = p12)
  } else{
    coloc_res = coloc::coloc.abf(dataset1 = df1,
                                 dataset2 = list(beta = gwas$log_OR, 
                                                 varbeta = gwas$se^2, 
                                                 type = "cc", 
                                                 s = case_control_proportion,
                                                 snp = gwas$rsid),
                                 p1 = p1, p2 = p2, p12 = p12)
  }
  
  return(coloc_res)
}


#' Import a specific region from a tabix-indexed GWAS summary stats file
tabixFetchGWASSummary <- function(granges, summary_path) {
  gwas_col_names = c("rsid", "chr", "pos", "effect_allele", "MAF", 
                     "p_nominal", "beta", "OR", "log_OR", "se", "z_score", "trait", "PMID", "used_file")
  gwas_col_types = c("ccicdddddddccc")
  #  gwas_col_names = c("snp_id", "chr", "pos", "MAF", "p_nominal", "beta", "se")
  #  gwas_col_types = c("ccidddd")
  gwas_pvalues = scanTabixDataFrame(summary_path, granges, col_names = gwas_col_names, col_types = gwas_col_types)
  return(gwas_pvalues)
}


#' A general function to quickly import tabix indexed tab-separated files into data_frame
#'
#' @param tabix_file Path to tabix-indexed text file
#' @param param A instance of GRanges, RangedData, or RangesList 
#' provide the sequence names and regions to be parsed. Passed onto Rsamtools::scanTabix()
#' @param ... Additional parameters to be passed on to readr::read_delim()
#'
#' @return List of data_frames, one for each entry in the param GRanges object.
#' @export
scanTabixDataFrame <- function(tabix_file, param, ...) {
  tabix_list = Rsamtools::scanTabix(tabix_file, param = param)
  df_list = lapply(tabix_list, function(x,...){
    if(length(x) > 0){
      if(length(x) == 1){
        #Hack to make sure that it also works for data frames with only one row
        #Adds an empty row and then removes it
        result = paste(paste(x, collapse = "\n"),"\n",sep = "")
        result = readr::read_delim(result, delim = "\t", ...)[1,]
      }else{
        result = paste(x, collapse = "\n")
        result = readr::read_delim(result, delim = "\t", ...)
      }
    } else{
      #Return NULL if the nothing is returned from tabix file
      result = NULL
    }
    return(result)
  }, ...)
  return(df_list)
}


summaryReplaceCoordinates <- function(summary_df, variant_information) {
  #Remove MAF if it is present
  if ("MAF" %in% colnames(summary_df)) {
    summary_df = dplyr::select(summary_df, -chr, -pos, -MAF) %>%
      dplyr::left_join(variant_information, by = "rsid") %>%
      dplyr::filter(!is.na(pos))
  } else {
    summary_df = dplyr::select(summary_df, -chr, -pos) %>%
      dplyr::left_join(variant_information, by = "rsid") %>%
      dplyr::filter(!is.na(pos))
  }
  return(summary_df)
}

getPlotDF <- function(qtl_summaries, gwas_summaries, coloc_pair) {
  plot.df = bind_rows(qtl_summaries %>% dplyr::select(one_of(c("rsid", "pos", "p_nominal"))) %>% dplyr::mutate(study = "QTL"),
                      gwas_summaries %>% dplyr::select(one_of(c("rsid", "pos", "p_nominal"))) %>% dplyr::mutate(study = "GWAS"))
  
  # Determine which variants are present in both studies
  sharing.df = dplyr::full_join(qtl_summaries %>% dplyr::select(one_of(c("rsid", "pos"))) %>% dplyr::rename(qtl_pos = pos),
                                gwas_summaries %>% dplyr::select(one_of(c("rsid", "pos"))) %>% dplyr::rename(gwas_pos = pos),
                                by="rsid")
  sharing.df$variant = "Shared"
  sharing.df$variant[is.na(sharing.df$gwas_pos) & !is.na(sharing.df$qtl_pos)] = "Not shared"
  sharing.df$variant[is.na(sharing.df$qtl_pos) & !is.na(sharing.df$gwas_pos)] = "Not shared"
  plot.df = plot.df %>% dplyr::left_join(sharing.df %>% dplyr::select(rsid, variant), by="rsid")
  
  plot.df$label = NA
  plot.df$label[plot.df$study == "GWAS" & plot.df$rsid == coloc_pair$gwas_rsid] = coloc_pair$gwas_rsid
  plot.df$label[plot.df$study == "QTL" & plot.df$rsid == coloc_pair$qtl_rsid] = coloc_pair$qtl_rsid
  plot.df
}

colocOverlayPlot <- function(qtl_summaries, gwas_summaries, coloc_pair, coloc_summary) {
  plot.df = getPlotDF(qtl_summaries, gwas_summaries, coloc_pair)
  
  xmax = max(plot.df$pos)
  xrange = xmax - min(plot.df$pos)
  H4ratio = coloc_summary$PP.H4.abf / (coloc_summary$PP.H3.abf + coloc_summary$PP.H4.abf)
  colocLabel = sprintf("H4 = %.3f\nH3 = %.3f\nH4/(H3+H4) = %.3f", coloc_summary$PP.H4.abf, coloc_summary$PP.H3.abf, H4ratio)
  
  if (!is.null(gene_id_map) & grepl("^ENSG", coloc_pair$feature)) {
    coloc_pair$feature = gsub("\\.[\\d]+", "", coloc_pair$feature, perl=T)
    coloc_pair$feature = sprintf("%s / %s", gene_id_map[gene_id_map$ensembl_gene_id == coloc_pair$feature,][1,]$symbol, coloc_pair$feature)
  }
  plotTitle = sprintf("gwas %s, qtl %s: %s", coloc_pair$gwas_rsid, coloc_pair$qtl_rsid, coloc_pair$feature)
  
  colocPlot = ggplot(plot.df, aes(x=pos, y=-log10(p_nominal), color=study, shape=variant)) +
    geom_point(alpha=0.7) + ggtitle(plotTitle) +
    theme_bw() + scale_shape_manual(values=c(4, 19)) +
    geom_text(aes(label=label), color="black", size=3, hjust="left", nudge_x=(xrange/100)) +
    annotate("text", max(plot.df$pos), max(-log10(plot.df$p_nominal)), label=colocLabel, hjust=1, vjust=1, size=3)
  return(colocPlot)
}

colocLinesPlot <- function(qtl_summaries, gwas_summaries, coloc_pair, coloc_summary) {
  plot.df = getPlotDF(qtl_summaries, gwas_summaries, coloc_pair)
  
  plotTitle = sprintf("gwas %s, qtl %s: %s", coloc_pair$gwas_rsid, coloc_pair$qtl_rsid, coloc_pair$feature)
  H4ratio = coloc_summary$PP.H4.abf / (coloc_summary$PP.H3.abf + coloc_summary$PP.H4.abf)
  colocLabel = sprintf("H4 = %.3f\nH3 = %.3f\nH4/(H3+H4) = %.3f", coloc_summary$PP.H4.abf, coloc_summary$PP.H3.abf, H4ratio)
  
  # colocPlot = ggplot(plot.df, aes(x=study, y=-log10(p_nominal), color=rsid, shape=variant, group=rsid)) +
  #   geom_jitter(data=plot.df, mapping=aes(x=study, y=-log10(p_nominal), group=study, color=rsid, shape=variant), alpha=0.7, width=0.1) +
  #   geom_line(alpha=0.3) +
  #   theme_bw() + ggtitle(plotTitle) + scale_shape_manual(values=c(4, 19)) +
  #   geom_text(aes(label=label), color="black", size=3, hjust="left") +
  #   annotate("text", "QTL", max(-log10(plot.df$p_nominal)), label=colocLabel, hjust=1, vjust=1, size=3) +
  #   theme(legend.position = "None")
  
  qtlYMax  = min(300, max(-log10(plot.df %>% dplyr::filter(study == "QTL") %>% .$p_nominal)))
  gwasYMax = min(300, max(-log10(plot.df %>% dplyr::filter(study == "GWAS") %>% .$p_nominal)))
  
  plot.df$log10p = -log10(plot.df$p_nominal)
  plot.df[plot.df$study == "QTL",]$log10p = plot.df[plot.df$study == "QTL",]$log10p * gwasYMax / qtlYMax
  
  colocPlot = ggplot(plot.df, aes(x=study, y=log10p, color=rsid, shape=variant, group=rsid)) +
    geom_jitter(data=plot.df, mapping=aes(x=study, y=log10p, group=study, color=rsid, shape=variant), alpha=0.7, width=0.1) +
    geom_line(alpha=0.3) +
    theme_bw() + ggtitle(plotTitle) + scale_shape_manual(values=c(4, 19)) +
    theme(legend.position = "None") +
    scale_y_continuous(name = "GWAS:  -log10(p)",
                       sec.axis = sec_axis(~ . * qtlYMax / gwasYMax, name = "QTL:  -log10(p)"), limits = c(0, gwasYMax))

  return(colocPlot)
}


###############################################################################

main()

