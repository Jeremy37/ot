#!/usr/bin/env Rscript
library("dplyr")
library("tidyr")
library("purrr")
library("coloc")
library("readr")
library("devtools")
library("optparse")

root = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys"
#root = "/Users/jeremys/work/opentargets"
load_all(file.path(root, "software/seqUtils/"))

#Parse command-line options
option_list <- list(
  make_option(c("--phenotype"), type="character", default=NULL,
              help="Type of QTLs used for coloc.", metavar = "type"),
  make_option(c("-w", "--window"), type="character", default=NULL,
              help="Size of the cis window.", metavar = "type"),
  make_option(c("--overlapdist"), type="character", default="1e5",
              help="Max distance between QTL lead SNP and GWAS lead SNP.", metavar = "type"),
  make_option(c("-g", "--gwas"), type="character", default=NULL,
              help="Name of the GWAS trait", metavar = "type"),
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="Path to GWAS summary stats directory.", metavar = "type"),
  make_option(c("-q", "--qtl"), type="character", default=NULL,
              help="Path to the QTL base filename", metavar = "type"),
  make_option(c("--outdir"), type="character", default=NULL,
              help="Path to the output directory.", metavar = "type"),
  make_option(c("--pthresh"), type="numeric", default=NULL,
              help="Considering the minimum nominal p value across SNPs for a molecular feature, pthresh is the maximum minp value at which a feature will be considered a QTL.", metavar = "type"),
  make_option(c("--p1"), type="numeric", default=NULL,
              help="Prior probability of SNP affecting expression", metavar = "type"),
  make_option(c("--p2"), type="numeric", default=NULL,
              help="Prior probability of SNP affecting GWAS trait", metavar = "type"),
  make_option(c("--p12"), type="numeric", default=NULL,
              help="Prior probability of SNP affecting both traits", metavar = "type"),
  make_option(c("-s", "--samplesize"), type="integer", default=NULL,
              help="Integer sample size of the QTL study", metavar = "type")
)
opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

#opt = list(g = "AD.IGAP1", w = "2e5", overlapdist = "1e5", phenotype = "xQTL_eQTL", d = "/Users/jeremys/work/opentargets/datasets/GWAS", o = "/Users/jeremys/work/opentargets/coloc/xQTL/",
#           q = "/Users/jeremys/work/opentargets/reference/xQTL", pthresh = 1e-5, s = 494, p1 = 1e-3, p2 = 1e-4, p12 = 5e-5)

#Extract parameters for CMD options
gwas_id = opt$g
cis_window = as.numeric(opt$w)
overlap_dist = as.numeric(opt$overlapdist)
phenotype = opt$phenotype
gwas_dir = opt$d
qtl_base = paste0(opt$q, "/", phenotype)
outdir = opt$outdir
sample_size = opt$s
p1 = opt$p1
p2 = opt$p2
p12 = opt$p12

gwas_thresh = 1e-5

gwas_dir = file.path(root, "datasets/GWAS")
#qtl_dir = file.path(root, "macrophage/qtltools/output")
#outdir = file.path(root, "coloc/blueprint")

importVariantInfo = function(path) {
  info_col_names = c("chr","pos","snp_id","ref","alt","MAF")
  into_col_types = "cicccd"
  readr::read_delim(path, delim = "\t", col_types = into_col_types, col_names = info_col_names)
}

#Import variant information
GRCh37_variants = importVariantInfo(file.path(root, "reference/xQTL", paste0(phenotype, ".variant_info.txt.gz")))

#Import list of GWAS studies
# gwas_stats_labeled = readr::read_tsv("/lustre/scratch117/cellgen/team170/ka8/projects/macrophage-trQTLs/analysis/data/gwas/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name","type"))
gwas_stats_labeled = readr::read_tsv(file.path(root, "coloc/GWAS_summary_stat_list.labeled.txt"), col_names = c("trait","file_name","type"))

#Construct a new QTL list
sample_sizes_list = list(naive=sample_size)

importQTLs = function(min_pvalue_path) {
  col_names = c("phenotype_id","chr","pos","snp_id","featureChr","featureStart","SpearmanRho","p_val")
  col_types = "cciccidd"
  table = readr::read_delim(min_pvalue_path, col_names = col_names, delim = "\t", col_types = col_types) %>%
    dplyr::filter(!is.na(p_val)) %>%
    dplyr::arrange(p_val)
  return(table)
}

constructQtlList <- function(phenotype, qtl_base, sample_size_list){
  conditions = names(sample_size_list)
  min_pvalues = list()
  qtl_summary_list = list()
  
  #Iterate over conditions and fill lists
  for(condition in conditions){
    min_pvalue_path = paste0(qtl_base, ".gene_minp.txt.gz")
    summary_path = paste0(qtl_base, ".sorted.txt.gz")
    
    min_pvalues[[condition]] = importQTLs(min_pvalue_path) %>% dplyr::select(phenotype_id, snp_id, p_val)
    qtl_summary_list[[condition]] = summary_path
  }
  return(list(min_pvalues = min_pvalues, qtl_summary_list = qtl_summary_list, sample_sizes = sample_size_list))
}

phenotype_values = constructQtlList(phenotype, qtl_base, sample_sizes_list)

#Spcecify the location of the GWAS summary stats file
gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == gwas_id)$file_name
gwas_prefix = file.path(gwas_dir, gwas_file_name)

#Prefilter coloc candidates
qtl_df_list = prefilterColocCandidates(phenotype_values$min_pvalues, gwas_prefix, 
                                       GRCh37_variants, fdr_thresh = NA, p_thresh = opt$pthresh, 
                                       overlap_dist = overlap_dist, gwas_thresh = gwas_thresh)
qtl_pairs = purrr::map_df(qtl_df_list, identity) %>% unique()

#Test for coloc
# coloc_res_list = purrr::map2(phenotype_values$qtl_summary_list, phenotype_values$sample_sizes, 
#                              ~colocMolecularQTLsByRow(qtl_pairs, qtl_summary_path = .x, 
#                                                       gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), 
#                                                       gwas_variant_info = GRCh37_variants,
#                                                       qtl_variant_info = GRCh38_variants, 
#                                                       N_qtl = .y, cis_dist = cis_window))

xQTLTabixFetchQTLs = function(phenotype_ranges, tabix_file) {
  #Assertions about input
  assertthat::assert_that(class(phenotype_ranges) == "GRanges")
  assertthat::assert_that(assertthat::has_name(GenomicRanges::elementMetadata(phenotype_ranges), "phenotype_id"))
  
  qtl_columns = c("phenotype_id", "snp_chr", "snp_pos", "snp_id", "featureChr", "featureStart", "SpearmanRho", "p.value")
  qtl_coltypes = "cciccidd"
  
  result = list()
  for (i in seq_along(phenotype_ranges)){
    selected_phenotype_id = phenotype_ranges[i]$phenotype_id
    print(i)
    tabix_table = scanTabixDataFrame(tabix_file, phenotype_ranges[i], 
                                     col_names = qtl_columns, col_types = qtl_coltypes)[[1]] %>%
      dplyr::filter(phenotype_id == selected_phenotype_id)
    
    #Add additional columns
    result[[selected_phenotype_id]] = tabix_table
  }
  return(result)
}



coloc_res_list = list()
for (qtlType in names(phenotype_values$sample_sizes)) {
  coloc_res = purrrlyr::by_row(qtl_pairs, ~colocMolecularQTLs(., qtl_summary_path = phenotype_values$qtl_summary_list[[qtlType]],
                                                              qtl_type = "xQTL_eQTL",
                                                              gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), 
                                                              gwas_variant_info = GRCh37_variants,
                                                              qtl_variant_info = GRCh37_variants,
                                                              N_qtl = phenotype_values$sample_sizes[[qtlType]], cis_dist = cis_window,
                                                              p1 = p1, p2 = p2, p12 = p12)$summary,
                               .collate = "rows")
  coloc_res_list[[qtlType]] = coloc_res
}

#Export results
coloc_hits = purrr::map_df(coloc_res_list, identity, .id = "condition_name") %>% dplyr::arrange(gwas_lead)
coloc_hits$phenotype_id = gsub("\\.[\\d]+", "", coloc_hits$phenotype_id, perl=T)
  
# Annotate gene names
gene.df = readr::read_tsv(file.path(root, "reference/GRCh38/gencode.v27.geneid_to_hgnc.txt"))
gene.df$ensemblID = gsub("\\.[\\d]+", "", gene.df$ensemblID, perl=T)
coloc_hits2 = coloc_hits %>% dplyr::left_join(gene.df, by=c("phenotype_id"="ensemblID")) %>%
  dplyr::rename(condition=condition_name, ensemblID=phenotype_id, PP.H0=PP.H0.abf, PP.H1=PP.H1.abf, PP.H2=PP.H2.abf, PP.H3=PP.H3.abf, PP.H4=PP.H4.abf) %>%
  dplyr::select(condition, ensemblID, geneSymbol, nsnps, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, qtl_pval, gwas_pval, qtl_lead, gwas_lead) %>%
  dplyr::filter(gwas_pval <= gwas_thresh)

coloc_output = file.path(outdir, paste("coloc", gwas_id, phenotype, overlap_dist, "txt", sep = "."))
write.table(coloc_hits2, coloc_output, sep = "\t", quote = FALSE, row.names = FALSE)

#Debugging example
#colocMolecularQTLs(qtl_pairs[1,], qtl_summary_list$Ctrl, gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), GRCh37_variants, GRCh38_variants, N_qtl = 2e5)
