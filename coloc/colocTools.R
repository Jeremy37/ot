
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
                               p1 = 1e-04, p2 = 1e-04, p12 = 1e-05) {
  
  expectedCols = c("feature", "chr", "qtl_pos", "gwas_pos", "qtl_rsid", "gwas_rsid")
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

    # Fetch QTL summary stats. The code calling colocMolecularQTLs needs to define
    # the function tabixFetchQTLs, which will handle fetching QTLs from potentially
    # differently formatted files.
    qtl_summaries = tabixFetchQTLs(qtlRange, qtl_summary_path)
    #if (!is.null(qtl_variant_info)) {
    #  qtl_summaries = summaryReplaceCoordinates(qtl_summaries, qtl_variant_info)
    #}
    
    #Fetch GWAS summary stats
    gwas_summaries = tabixFetchGWASSummary(gwasRange, gwas_summary_path)[[1]]
    
    qtl_min_index = which.min(qtl_summaries$p_nominal)
    gwas_min_index = which.min(gwas_summaries$p_nominal)
    
    #Perform coloc analysis
    coloc_res = colocQtlGWAS(qtl_summaries, gwas_summaries, N_qtl = N_qtl, p1, p2, p12)
    coloc_summary = dplyr::tbl_df(t(data.frame(coloc_res$summary))) %>%
      dplyr::mutate(qtl_pval = qtl_summaries[qtl_min_index,]$p_nominal, gwas_pval = gwas_summaries[gwas_min_index,]$p_nominal,
                    qtl_lead = qtl_summaries[qtl_min_index,]$rsid, gwas_lead = gwas_summaries[gwas_min_index,]$rsid)
    
    #Summary list
    data_list = list(qtl = qtl, gwas = gwas)
    
    result = list(summary = coloc_summary, data = data_list)
  }, error = function(err) {
    print(paste("ERROR:",err))
    result = list(summary = NULL, data = NULL)
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
  nQTL = nrow(qtl)
  nGWAS = nrow(gwas)
  qtl = qtl %>% dplyr::inner_join(gwas %>% dplyr::select(rsid))
  gwas = gwas %>% dplyr::inner_join(qtl %>% dplyr::select(rsid))
  write(sprintf("%d SNPs in common, of %d QTL SNPs and %d GWAS SNPs", nrow(qtl), nQTL, nGWAS), stderr())
  
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
                                 p1, p2, p12)
  } else{
    coloc_res = coloc::coloc.abf(dataset1 = df1,
                                 dataset2 = list(beta = gwas$log_OR, 
                                                 varbeta = gwas$se^2, 
                                                 type = "cc", 
                                                 s = case_control_proportion,
                                                 snp = gwas$rsid),
                                 p1, p2, p12)
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
    summary_df = dplyr::select(summary_df, -chr, -qtl_pos, -MAF) %>%
      dplyr::left_join(variant_information, by = "rsid") %>%
      dplyr::filter(!is.na(pos)) %>%
      dplyr::rename(qtl_pos = pos)
  } else {
    summary_df = dplyr::select(summary_df, -chr, -qtl_pos) %>%
      dplyr::left_join(variant_information, by = "rsid") %>%
      dplyr::filter(!is.na(pos)) %>%
      dplyr::rename(qtl_pos = pos)
  }
  return(summary_df)
}
