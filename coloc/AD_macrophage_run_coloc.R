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
              help="Path to the QTL directory.", metavar = "type"),
  make_option(c("--outdir"), type="character", default=NULL,
              help="Path to the output directory.", metavar = "type"),
  make_option(c("--p1"), type="numeric", default=NULL,
              help="Prior probability of SNP affecting expression", metavar = "type"),
  make_option(c("--p2"), type="numeric", default=NULL,
              help="Prior probability of SNP affecting GWAS trait", metavar = "type"),
  make_option(c("--p12"), type="numeric", default=NULL,
              help="Prior probability of SNP affecting both traits", metavar = "type"),
  make_option(c("-s", "--samplesizes"), type="character", default=NULL,
              help="Path to the tab-separated text file with condition names and sample sizes.", metavar = "type")
)
opt <- parse_args(OptionParser(option_list=option_list))
print(opt)
# 
# p1 = 1e-4  # prior prob for SNP associated with QTL
# p2 = 1e-4  # prior prob for SNP associated with GWAS
# p12 = 1e-5 # prior probl for SNP associated with both traits

#Debugging
#opt = list(g = "AD.meta", w = "2e5", overlapdist = "1e5", phenotype = "featureCounts",
# opt = list(g = "PD.GWAX", w = "2e5", overlapdist = "1e5", phenotype = "featureCounts",
#            d = "/Users/jeremys/work/opentargets/datasets/GWAS",
#            outdir = "/Users/jeremys/work/opentargets/coloc/ipsmacrophage/",
#            q = "/Users/jeremys/work/opentargets/macrophage/qtltools/output/",
#            s = "/Users/jeremys/work/opentargets/coloc/salmonella_coloc_sample_sizes.txt",
#            p1 = 1e-3, p2 = 1e-4, p12 = 5e-5)

#Extract parameters for CMD options
gwas_id = opt$g
cis_window = as.numeric(opt$w)
overlap_dist = as.numeric(opt$overlapdist)
phenotype = opt$phenotype
gwas_dir = opt$d
qtl_dir = opt$q
outdir = opt$outdir
sample_size_path = opt$s
p1 = opt$p1  # prior prob for SNP associated with QTL
p2 = opt$p2  # prior prob for SNP associated with GWAS
p12 = opt$p12 # prior probl for SNP associated with both traits

gwas_thresh = 1e-5

#Import variant information
GRCh38_variants = importVariantInformation(file.path(root, "macrophage/genotypes/imputed.86_samples.variant_information.txt.gz"))
GRCh37_variants = importVariantInformation(file.path(root, "macrophage/genotypes/GRCh37/imputed.86_samples.variant_information.GRCh37.vcf.gz"))

#Import list of GWAS studies
# gwas_stats_labeled = readr::read_tsv("/lustre/scratch117/cellgen/team170/ka8/projects/macrophage-trQTLs/analysis/data/gwas/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name","type"))
gwas_stats_labeled = readr::read_tsv(file.path(root, "coloc/GWAS_summary_stat_list.labeled.txt"), col_names = c("trait","file_name","type"))

#Import sample sizes
sample_sizes = readr::read_tsv(sample_size_path, col_names = c("condition_name", "sample_size"), col_types = "ci")
sample_sizes_list = as.list(sample_sizes$sample_size)
names(sample_sizes_list) = sample_sizes$condition_name

tabixFetchQTLs = function(phenotype_ranges, tabix_file, qtl_type) {
  if (qtl_type == "QTLTools") {
    qtl_summaries = qtltoolsTabixFetchPhenotypes(phenotype_ranges, tabix_file)[[1]] %>%
      dplyr::transmute(snp_id, chr = snp_chr, pos = snp_start, p_nominal, beta)
  }
  else if (qtl_type == "FastQTL") {
    qtl_summaries = fastqtlTabixFetchGenes(phenotype_ranges, tabix_file)[[1]]
  } else {
    assertthat::assert_that(FALSE, msg="qtl_type not recognized!")
  }
}

#Construct a new QTL list 
phenotype_values = constructQtlListForColoc(phenotype, qtl_dir, sample_sizes_list)

#Spcecify the location of the GWAS summary stats file
gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == gwas_id)$file_name
gwas_prefix = file.path(gwas_dir, gwas_file_name)

#Prefilter coloc candidates
qtl_df_list = prefilterColocCandidates(phenotype_values$min_pvalues, gwas_prefix, 
                                       GRCh37_variants, fdr_thresh = 0.1, 
                                       overlap_dist = overlap_dist, gwas_thresh = gwas_thresh)
qtl_pairs = purrr::map_df(qtl_df_list, identity) %>% unique()

qtl_pairs = qtl_pairs[1:15, ]

#Test for coloc
# coloc_res_list = purrr::map2(phenotype_values$qtl_summary_list, phenotype_values$sample_sizes, 
#                              ~colocMolecularQTLsByRow(qtl_pairs, qtl_summary_path = .x, 
#                                                       gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), 
#                                                       gwas_variant_info = GRCh37_variants,
#                                                       qtl_variant_info = GRCh38_variants, 
#                                                       N_qtl = .y, cis_dist = cis_window))

print(paste0("Number of QTL pairs: ", nrow(qtl_pairs)))
print(paste0("Number of conditions: ", length(sample_sizes_list)))

coloc_res_list = list()
for (qtlType in names(phenotype_values$sample_sizes)) {
  coloc_res = purrrlyr::by_row(qtl_pairs, ~colocMolecularQTLs(., qtl_type = "QTLTools",
                                                              qtl_summary_path = phenotype_values$qtl_summary_list[[qtlType]], 
                                                              gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), 
                                                              gwas_variant_info = GRCh37_variants,
                                                              qtl_variant_info = GRCh38_variants, 
                                                              N_qtl = phenotype_values$sample_sizes[[qtlType]], cis_dist = cis_window,
                                                              p1 = p1, p2 = p2, p12 = p12)$summary,
                               .collate = "rows")
  coloc_res_list[[qtlType]] = coloc_res
}

#Export results
coloc_hits = purrr::map_df(coloc_res_list, identity, .id = "condition_name") %>% dplyr::arrange(gwas_lead)

if (phenotype == "reviseAnnotations") {
  coloc_hits$reviseAnn = coloc_hits$phenotype_id
  coloc_hits$phenotype_id = sapply(coloc_hits$phenotype_id, function(x) strsplit(x, ".", fixed=T)[[1]][1])

  #coloc_hits = readr::read_tsv(file.path(root, "coloc", "ipsmacrophage", "coloc.AD.meta.reviseAnnotations.1e+05.txt")) %>%
  #  dplyr::rename(phenotype_id = ensemblID) %>% dplyr::select(-geneSymbol)
  #coloc_hits2 = coloc_hits %>% dplyr::left_join(gene.df, by=c("phenotype_id"="ensemblID")) %>%
  #  dplyr::rename(ensemblID=phenotype_id)
  #coloc_hits2 %<>% dplyr::select(condition, reviseAnn, ensemblID, geneSymbol, nsnps, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, qtl_pval, gwas_pval, qtl_lead, gwas_lead)
}

# Annotate gene names
gene.df = readr::read_tsv(file.path(root, "reference/GRCh38/gencode.v27.geneid_to_hgnc.txt"))
gene.df$ensemblID = gsub("\\.[\\d]+", "", gene.df$ensemblID, perl=T)
coloc_hits2 = coloc_hits %>% dplyr::left_join(gene.df, by=c("phenotype_id"="ensemblID")) %>%
  dplyr::rename(condition=condition_name, ensemblID=phenotype_id, PP.H0=PP.H0.abf, PP.H1=PP.H1.abf, PP.H2=PP.H2.abf, PP.H3=PP.H3.abf, PP.H4=PP.H4.abf)

print(sprintf("%d coloc hits, and %d hits where the gwas p value for overlapping SNPs was less than the threshold %g",
              nrow(coloc_hits2), sum(coloc_hits2$gwas_pval <= gwas_thresh), gwas_thresh))

coloc_hits2 = coloc_hits2 %>% dplyr::filter(gwas_pval <= gwas_thresh)
if (phenotype == "reviseAnnotations") {
  # Also select reviseAnn column
  coloc_hits2 = coloc_hits2 %>% dplyr::select(condition, reviseAnn, ensemblID, geneSymbol, nsnps,
                                              PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, qtl_pval, gwas_pval,
                                              qtl_lead, gwas_lead, gwas_chr, gwas_pos)
} else {
  coloc_hits2 = coloc_hits2 %>% dplyr::select(condition, ensemblID, geneSymbol, nsnps,
                                              PP.H0, PP.H1, PP.H2, PP.H3, PP.H4, qtl_pval, gwas_pval,
                                              qtl_lead, gwas_lead, gwas_chr, gwas_pos)
}

coloc_output = file.path(outdir, paste("coloc", gwas_id, phenotype, overlap_dist, "txt", sep = "."))
write.table(coloc_hits2, coloc_output, sep = "\t", quote = FALSE, row.names = FALSE)

#Debugging example
#colocMolecularQTLs(qtl_pairs[1,], qtl_summary_list$Ctrl, gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), GRCh37_variants, GRCh38_variants, N_qtl = 2e5)
