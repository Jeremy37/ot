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
  make_option(c("--qtlname"), type="character", default=NULL,
              help="Name of QTLs, used for output."),
  make_option(c("--qtltype"), type="character", default=NULL,
              help="Type of QTLs, supported are: sens_neur, xQTL, blueprint, QTLTools"),
  make_option(c("--window"), type="numeric", default=NULL,
              help="Size of the cis window."),
  make_option(c("--overlapdist"), type="numeric", default="1e5",
              help="Max distance between QTL lead SNP and GWAS lead SNP."),
  make_option(c("--gwasname"), type="character", default=NULL,
              help="Name of the GWAS trait"),
  make_option(c("--gwasdir"), type="character", default=NULL,
              help="Path to GWAS summary stats directory."),
  make_option(c("--gwasthresh"), type="numeric", default=NULL,
              help="Threshold of GWAS p value to use locus for colocalisation"),
  make_option(c("--qtlnominal"), type="character", default=NULL,
              help="Path to the QTL nominal p values file"),
  make_option(c("--qtlgeneminp"), type="character", default=NULL,
              help="Path to the QTL min p value per gene file"),
  make_option(c("--qtlvariantinfo"), type="character", default=NULL,
              help="Path to the QTL variant info file"),
  make_option(c("--gwasvariantinfo"), type="character", default=NULL,
              help="Path to the GWAS variant info file. If GWAS uses same coords/genome build as QTLs, this does not need to be provided."),
  make_option(c("--outdir"), type="character", default=NULL,
              help="Path to the output directory."),
  make_option(c("--qtlpthresh"), type="numeric", default=NA,
              help="Considering the minimum nominal p value across SNPs for a molecular feature, qtlpthresh is the maximum minp value at which a feature will be considered a QTL."),
  make_option(c("--qtlfdrthresh"), type="numeric", default=NA,
              help="qtlfdrthresh is the maximum FDR at which a feature will be considered a QTL. Assumes a column named 'fdr' is present in the qtlgeneminp file."),
  make_option(c("--p1"), type="numeric", default=NULL,
              help="Prior probability of SNP affecting expression"),
  make_option(c("--p2"), type="numeric", default=NULL,
              help="Prior probability of SNP affecting GWAS trait"),
  make_option(c("--p12"), type="numeric", default=NULL,
              help="Prior probability of SNP affecting both traits"),
  make_option(c("--samplesize"), type="integer", default=NULL,
              help="Integer sample size of the QTL study")
)
opt <- parse_args(OptionParser(option_list=option_list))

# opt = list(gwasname = "AD.meta", gwasdir = "/Users/jeremys/work/opentargets/datasets/GWAS",
#            gwasthresh = 1e-5, window = 2e5, overlapdist = 1e5,
#            qtlname = "xQTL_eQTL", qtlnominal = "/Users/jeremys/work/opentargets/reference/xQTL/xQTL_eQTL.sorted.txt.gz",
#            qtlgeneminp = "/Users/jeremys/work/opentargets/reference/xQTL/xQTL_eQTL.gene_minp.txt.gz", 
#            qtlvariantinfo = "/Users/jeremys/work/opentargets/reference/xQTL/xQTL_eQTL.variant_info.txt.gz", 
#            outdir = "/Users/jeremys/work/opentargets/coloc/xQTL",
#            qtlpthresh = 1e-5, qtlfdrthresh = NA,
#            p1 = 1e-3, p2 = 1e-4, p12 = 5e-5, samplesize = 494)
#  opt = list(gwasname = "AD.meta", gwasdir = "/Users/jeremys/work/opentargets/datasets/GWAS",
#             gwasthresh = 1e-5, window = 2e5, overlapdist = 1e5,
#             qtlname = "sens_neur.100k", qtltype = "sens_neur",
#             qtlnominal = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/fastqtl.nominals.cis500k.ePCs.20.coords.txt.gz",
#             qtlgeneminp = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/fastqtl.permutations.10k.allgenes.cis100k.PCs.20.fdr0.1.txt",
#             qtlvariantinfo = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/imputed.97_samples.snps_indels.INFO_08.variant_info.GRCh38.gz",
#             gwasvariantinfo = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh37/imputed.97_samples.snps_indels.INFO_08.variant_info.GRCh37.gz",
#             outdir = "/Users/jeremys/work/opentargets/coloc/sensoryneuron",
#             qtlpthresh = NA, qtlfdrthresh = 0.1,
#             p1 = 1e-3, p2 = 1e-4, p12 = 5e-5, samplesize = 97)
# opt = list(gwasname = "AD.meta", gwasdir = "/Users/jeremys/work/opentargets/datasets/GWAS",
#            gwasthresh = 1e-5, window = 2e5, overlapdist = 1e5,
#            qtlname = "sens_neur.sqtl", qtltype = "sens_neur",
#            qtlnominal = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/fastqtl.sqtl.nominals.PCs.5.pvalues.coords.txt.gz",
#            qtlgeneminp = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/fastqtl.sqtl.permutations.10k.PCs.5.fdr0.1.txt",
#            qtlvariantinfo = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/imputed.97_samples.snps_indels.INFO_08.variant_info.GRCh38.gz",
#            gwasvariantinfo = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh37/imputed.97_samples.snps_indels.INFO_08.variant_info.GRCh37.gz",
#            outdir = "/Users/jeremys/work/opentargets/coloc/sensoryneuron",
#            qtlpthresh = NA, qtlfdrthresh = 0.1,
#            p1 = 1e-3, p2 = 1e-4, p12 = 5e-5, samplesize = 97)
# opt = list(gwasname = "AD.meta", gwasdir = "/Users/jeremys/work/opentargets/datasets/GWAS",
#            gwasthresh = 1e-5, window = 2e5, overlapdist = 1e5,
#            qtlname = "mono_gene_eQTL", qtltype = "blueprint",
#            qtlnominal = "/Users/jeremys/work/opentargets/datasets/blueprint/mono_gene_eQTL.sorted.txt.gz",
#            qtlgeneminp = "/Users/jeremys/work/opentargets/datasets/blueprint/mono_gene_eQTL.gene_minp.txt.gz",
#            qtlvariantinfo = "/Users/jeremys/work/opentargets/datasets/blueprint/mono_gene_eQTL.variant_info.txt.gz",
#            gwasvariantinfo = NULL,
#            outdir = "/Users/jeremys/work/opentargets/coloc/blueprint",
#            qtlpthresh = NA, qtlfdrthresh = 0.1,
#            p1 = 1e-3, p2 = 1e-4, p12 = 5e-5, samplesize = 97)
# 
# opt = list(gwasname = "PD.GWAX", gwasdir = "/Users/jeremys/work/opentargets/datasets/GWAS",
#            gwasthresh = 1e-5, window = 2e5, overlapdist = 1e5,
#            qtlname = "mono_gene_eQTL", qtltype = "blueprint",
#            qtlnominal = "/Users/jeremys/work/opentargets/datasets/blueprint/mono_gene_eQTL.sorted.txt.gz",
#            qtlgeneminp = "/Users/jeremys/work/opentargets/datasets/blueprint/mono_gene_eQTL.gene_minp.txt.gz",
#            qtlvariantinfo = "/Users/jeremys/work/opentargets/datasets/blueprint/mono_gene_eQTL.variant_info.txt.gz",
#            gwasvariantinfo = NULL,
#            outdir = "/Users/jeremys/work/opentargets/coloc/blueprint",
#            qtlpthresh = NA, qtlfdrthresh = 0.1,
#            p1 = 1e-3, p2 = 1e-4, p12 = 5e-5, samplesize = 97)

if (is.null(opt$qtltype)) {
   opt$qtltype = opt$qtlname
 }
print(opt)

importVariantInfo = function(path) {
  info_col_names = c("chr","pos","snp_id","ref","alt","MAF")
  into_col_types = "cicccd"
  readr::read_delim(path, delim = "\t", col_types = into_col_types, col_names = info_col_names)
}

# Import variant information. We use SNP rsIDs to match variants between GWAS and QTL studies,
# which enables coloc when the data are based on different genome builds (GRCh37 vs. 38)
qtl_variantinfo = importVariantInfo(opt$qtlvariantinfo)
if (is.null(opt$gwasvariantinfo)) {
  gwas_variantinfo = qtl_variantinfo
} else {
  gwas_variantinfo = importVariantInfo(opt$gwasvariantinfo)
}

#Import list of GWAS studies
# gwas_stats_labeled = readr::read_tsv("/lustre/scratch117/cellgen/team170/ka8/projects/macrophage-trQTLs/analysis/data/gwas/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name","type"))
gwas_stats_labeled = readr::read_tsv(file.path(root, "coloc/GWAS_summary_stat_list.labeled.txt"), col_names = c("trait","file_name","type"))

genericImportQTLs = function(min_pvalue_path) {
  # This reads a table for QTLs which has a header
  table = readr::read_delim(min_pvalue_path, delim = "\t") %>%
    dplyr::filter(!is.na(p_val)) %>%
    dplyr::arrange(FDR)
  return(table)
}

xQTLImportQTLs = function(min_pvalue_path) {
  col_names = c("phenotype_id","chr","pos","snp_id","featureChr","featureStart","SpearmanRho","p_val")
  col_types = "cciccidd"
  table = readr::read_delim(min_pvalue_path, col_names = col_names, delim = "\t", col_types = col_types) %>%
    dplyr::filter(!is.na(p_val)) %>%
    dplyr::arrange(p_val)
  return(table)
}

blueprintImportQTLs = function(min_pvalue_path) {
  col_names = c("phenotype_id","chr","pos","snp_id","p_val","beta","Bonferroni.p.value","FDR","MAF","std.error_of_beta")
  col_types = "ccicdddddd"
  table = readr::read_delim(min_pvalue_path, col_names = col_names, delim = "\t", col_types = col_types) %>%
    dplyr::filter(!is.na(p_val)) %>%
    dplyr::arrange(FDR)
  return(table)
}


tabixFetchQTLs = function(phenotype_ranges, tabix_file, qtl_type) {
  if (qtl_type == "QTLTools") {
    qtl_summaries = qtltoolsTabixFetchPhenotypes(phenotype_ranges, tabix_file)[[1]] %>%
      dplyr::transmute(snp_id, chr = snp_chr, pos = snp_start, p_nominal, beta)
  }
  else if (qtl_type == "FastQTL") {
    qtl_summaries = fastqtlTabixFetchGenes(phenotype_ranges, tabix_file)[[1]]
  }
  else if (qtl_type == "blueprint") {
    qtl_colnames = c("phenotype_id", "chr", "pos", "snp_id", "p_nominal",
                     "beta", "Bonferroni.p.value", "FDR", "MAF", "std.error_of_beta")
    qtl_summaries = tabixFetchQTLsHelper(phenotype_ranges, tabix_file, qtl_colnames, qtl_coltypes = "ccicdddddd")[[1]] %>%
      dplyr::transmute(snp_id, chr, pos, p_nominal, beta)
  }
  else if (qtl_type == "xQTL") {
    qtl_colnames = c("phenotype_id", "chr", "pos", "snp_id", "featureChr", "featureStart", "SpearmanRho", "p_nominal")
    qtl_summaries = tabixFetchQTLsHelper(phenotype_ranges, tabix_file, qtl_colnames, qtl_coltypes = "cciccidd")[[1]] %>%
      dplyr::transmute(snp_id, chr, pos, p_nominal)
  }
  else if (qtl_type == "sens_neur") {
    qtl_colnames = c("phenotype_id", "chr", "pos", "snp_id", "dist", "p_nominal", "beta")
    qtl_summaries = tabixFetchQTLsHelper(phenotype_ranges, tabix_file, qtl_colnames, qtl_coltypes = "ccicidd")[[1]] %>%
      dplyr::transmute(snp_id, chr, pos, p_nominal, beta)
  }
  else if (qtl_type == "hipsci") {
    qtl_colnames = c("chr", "pos", "phenotype_id", "beta", "p_nominal")
    qtl_summaries = tabixFetchQTLsHelper(phenotype_ranges, tabix_file, qtl_colnames, qtl_coltypes = "ccicidd")[[1]] %>%
      dplyr::transmute(snp_id = paste0(chr, "_", pos), chr, pos, p_nominal, beta)
  } else {
    assertthat::assert_that(FALSE, msg="qtl_type not recognized!")
  }
}

tabixFetchQTLsHelper = function(phenotype_ranges, tabix_file, qtl_colnames, qtl_coltypes) {
  #Assertions about input
  assertthat::assert_that(class(phenotype_ranges) == "GRanges")
  assertthat::assert_that(assertthat::has_name(GenomicRanges::elementMetadata(phenotype_ranges), "phenotype_id"))
  
  result = list()
  for (i in seq_along(phenotype_ranges)){
    selected_phenotype_id = phenotype_ranges[i]$phenotype_id
    print(i)
    tabix_table = scanTabixDataFrame(tabix_file, phenotype_ranges[i], 
                                     col_names = qtl_colnames, col_types = qtl_coltypes)[[1]] %>%
      dplyr::filter(phenotype_id == selected_phenotype_id)
    
    #Add additional columns
    result[[selected_phenotype_id]] = tabix_table
  }
  return(result)
}

# This is an anachronism from the coloc code for macrophages, which had multiple conditions
# Here we create a set of lists with one item each, for the "naive" condition
constructQtlList <- function(qtl_type, qtlgeneminp, qtlnominal, samplesize) {
  condition = "naive"
  min_pvalues = list()
  qtl_summary_list = list()
  sample_sizes_list = list()
  
  qtlTable = NULL
  if (qtl_type == "sens_neur" | qtl_type == "hipsci") {
    qtlTable = genericImportQTLs(qtlgeneminp)
  } else if (qtl_type == "blueprint") {
    qtlTable = blueprintImportQTLs(qtlgeneminp)
  } else if (qtl_type == "xQTL") {
    qtlTable = xQTLImportQTLs(qtlgeneminp)
  } else {
    stop(paste0("QTL type ", qtl_type, " not recognized!"))
  }
  if ("FDR" %in% colnames(qtlTable) & !is.null(opt$qtlfdrthresh)) {
    min_pvalues[[condition]] = qtlTable %>% dplyr::select(phenotype_id, snp_id, FDR)
  } else {
    min_pvalues[[condition]] = qtlTable %>% dplyr::select(phenotype_id, snp_id, p_val)
  }
  qtl_summary_list[[condition]] = qtlnominal
  sample_sizes_list[[condition]] = samplesize
  return(list(min_pvalues = min_pvalues, qtl_summary_list = qtl_summary_list, sample_sizes = sample_sizes_list))
}

qtl_values = constructQtlList(opt$qtltype, opt$qtlgeneminp, opt$qtlnominal, opt$samplesize)

#Spcecify the location of the GWAS summary stats file
gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == opt$gwasname)$file_name
gwas_prefix = file.path(opt$gwasdir, gwas_file_name)

#Prefilter coloc candidates
qtl_df_list = prefilterColocCandidates(qtl_values$min_pvalues, gwas_prefix, 
                                       gwas_variantinfo, fdr_thresh = opt$qtlfdrthresh, p_thresh = opt$qtlpthresh, 
                                       overlap_dist = opt$overlapdist, gwas_thresh = opt$gwasthresh)
qtl_pairs = purrr::map_df(qtl_df_list, identity) %>% unique()

#qtl_pairs = qtl_pairs[1:5,]
coloc_res = purrrlyr::by_row(qtl_pairs, ~colocMolecularQTLs(., qtl_summary_path = qtl_values$qtl_summary_list[["naive"]],
                                                              qtl_type = opt$qtltype,
                                                              gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), 
                                                              gwas_variant_info = gwas_variantinfo,
                                                              qtl_variant_info = qtl_variantinfo,
                                                              N_qtl = qtl_values$sample_sizes[["naive"]], cis_dist = opt$window,
                                                              p1 = opt$p1, p2 = opt$p2, p12 = opt$p12)$summary,
                             .collate = "rows")

#Export results
coloc_hits = coloc_res %>% dplyr::arrange(gwas_lead)
coloc_hits$phenotype_id = gsub("\\.[\\d]+", "", coloc_hits$phenotype_id, perl=T)

print(sprintf("%d coloc hits, and %d hits where the gwas p value for overlapping SNPs was less than the threshold %g",
              nrow(coloc_hits), sum(coloc_hits$gwas_pval <= opt$gwasthresh), opt$gwasthresh))

coloc_hits = coloc_hits %>%
  dplyr::rename(PP.H0=PP.H0.abf, PP.H1=PP.H1.abf, PP.H2=PP.H2.abf, PP.H3=PP.H3.abf, PP.H4=PP.H4.abf) %>%
  dplyr::select(phenotype_id, nsnps, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, qtl_pval, gwas_pval, qtl_lead, gwas_lead, gwas_chr, gwas_pos) %>%
  dplyr::filter(gwas_pval <= opt$gwasthresh)

coloc_output = file.path(opt$outdir, paste("coloc", opt$gwasname, opt$qtlname, opt$overlapdist, "txt", sep = "."))
write.table(coloc_hits, coloc_output, sep = "\t", quote = FALSE, row.names = FALSE)

