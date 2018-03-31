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
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Path to the output directory.", metavar = "type"),
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
# p1 = 1e-3  # prior prob for SNP associated with QTL
# p2 = 1e-4  # prior prob for SNP associated with GWAS
# p12 = 5e-5 # prior probl for SNP associated with both traits

#Debugging
#opt = list(g = "IBD", w = "2e5", p = "featureCounts", d = "/Volumes/JetDrive/datasets/Inflammatory_GWAS/", o = "results/acLDL/coloc/coloc_lists/",
#           q = "processed/salmonella/qtltools/output/", s = "analysis/data/sample_lists/salmonella_coloc_sample_sizes.txt")

#Extract parameters for CMD options
gwas_id = opt$g
cis_window = as.numeric(opt$w)
overlap_dist = as.numeric(opt$overlapdist)
phenotype = opt$phenotype
gwas_dir = opt$d
qtl_dir = opt$q
outdir = opt$o
sample_size = opt$s
p1 = opt$p1
p2 = opt$p2
p12 = opt$p12

gwas_thresh = 1e-5

# gwas_dir = "/lustre/scratch117/cellgen/team170/ka8/datasets/Inflammatory_GWAS"
# qtl_dir = "/lustre/scratch117/cellgen/team170/ka8/projects/macrophage-trQTLs/processed/salmonella/qtltools/output"
# outdir = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/coloc/AD/ipsmacrophage"
gwas_dir = file.path(root, "datasets/GWAS")
qtl_dir = file.path(root, "macrophage/qtltools/output")
outdir = file.path(root, "coloc/ipsmacrophage")


#Import variant information
GRCh38_variants = importVariantInformation(file.path(root, "macrophage/genotypes/imputed.86_samples.variant_information.txt.gz"))
GRCh37_variants = importVariantInformation(file.path(root, "macrophage/genotypes/GRCh37/imputed.86_samples.variant_information.GRCh37.vcf.gz"))

#Import list of GWAS studies
# gwas_stats_labeled = readr::read_tsv("/lustre/scratch117/cellgen/team170/ka8/projects/macrophage-trQTLs/analysis/data/gwas/GWAS_summary_stat_list.labeled.txt", col_names = c("trait","file_name","type"))
gwas_stats_labeled = readr::read_tsv(file.path(root, "coloc/GWAS_summary_stat_list.labeled.txt"), col_names = c("trait","file_name","type"))

#Construct a new QTL list
sample_sizes_list = list(sample_size)
names(sample_sizes_list) = "naive"
phenotype_values = constructQtlListForColoc(phenotype, qtl_dir, sample_sizes_list)

#Spcecify the location of the GWAS summary stats file
gwas_file_name = dplyr::filter(gwas_stats_labeled, trait == gwas_id)$file_name
gwas_prefix = file.path(gwas_dir, gwas_file_name)

#Prefilter coloc candidates
qtl_df_list = prefilterColocCandidates(phenotype_values$min_pvalues, gwas_prefix, 
                                       GRCh37_variants, fdr_thresh = 0.1, 
                                       overlap_dist = overlap_dist, gwas_thresh = gwas_thresh)
qtl_pairs = purrr::map_df(qtl_df_list, identity) %>% unique()

#Test for coloc
# coloc_res_list = purrr::map2(phenotype_values$qtl_summary_list, phenotype_values$sample_sizes, 
#                              ~colocMolecularQTLsByRow(qtl_pairs, qtl_summary_path = .x, 
#                                                       gwas_summary_path = paste0(gwas_prefix, ".sorted.txt.gz"), 
#                                                       gwas_variant_info = GRCh37_variants,
#                                                       qtl_variant_info = GRCh38_variants, 
#                                                       N_qtl = .y, cis_dist = cis_window))

coloc_res_list = list()
for (qtlType in names(phenotype_values$sample_sizes)) {
  coloc_res = purrrlyr::by_row(qtl_pairs, ~colocMolecularQTLs(., qtl_summary_path = phenotype_values$qtl_summary_list[[qtlType]], 
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
