#!/usr/bin/env Rscript
library("tidyverse")
library("Rsamtools")
library("devtools")
library("cowplot")

#root = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys"
root = "/Users/jeremys/work/opentargets"
load_all(file.path(root, "software/seqUtils/"))

#Parse command-line options
option_list <- list(
  make_option(c("--region1"), type="character", default=NULL,
              help="chr:start-end for region to plot"),
  make_option(c("--region2"), type="character", default=NULL,
              help="chr:start-end for region to plot"),
  make_option(c("--file1"), type="character", default=NULL,
              help="First file of summary stats (tabix-indexed)"),
  make_option(c("--file2"), type="character", default=NULL,
              help="Second file of summary stats (tabix-indexed)"),
  make_option(c("--snpcol1"), type="integer", default=NULL,
              help="1-based column with the file1 position to plot"),
  make_option(c("--pcol1"), type="integer", default=NULL,
              help="1-based column with the file1 statistic to plot"),
  make_option(c("--snpcol2"), type="integer", default=NULL,
              help="1-based column with the file2 position to plot"),
  make_option(c("--pcol2"), type="integer", default=NULL,
              help="1-based column with the file2 statistic to plot"),
  make_option(c("--geneid1"), type="character", default=NULL,
              help="Ensembl ID of the gene in file1 to plot summary stats for (can leave empty)"),
  make_option(c("--geneidcol1"), type="character", default=NULL,
              help=""),
  make_option(c("--matchvariants"), type="character", default=NULL,
              help=""),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Path to the output directory.")
)
opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

if (is.null(opt$file1) | is.null(opt$file2) | is.null(opt$region) | is.null(opt$poscol1) | is.null(opt$statcol1) | is.null(opt$poscol2) | is.null(opt$statcol2) | is.null(opt$outdir)) {
  stop("The following arguments are required: --file1 --file2 --region --col1 --col2 --outdir")
}

# BIN1 locus, hg37: 2:127700000-128000000
# BIN1 locus, hg38: 2:126942424-127242424
opt=list(region1="2:126942424-127242424",
         region2="2:127700000-128000000",
         file1="/Users/jeremys/work/opentargets/macrophage/qtltools/output/featureCounts/sorted/naive.nominal.sorted.txt.gz",
         file2="/Users/jeremys/work/opentargets/datasets/GWAS/Alzheimers_disease_Liu_2017_UKBB_GWAX_meta.sorted.txt.gz",
         file1colnames = "geneid,gene_chr,gene_start,gene_end,strand,n_snps,distance,snp_id,snp_chr,snp_pos1,snp_end,qtl_p,beta,is_lead",
         file2colnames = "snp_id,snp_chr,snp_pos2,A1,MAF,gwas_p,META_BETA,OR,log_OR,META_SE,z_score,trait,PMID,used_file",
         snpcol1="8", snpcol2="1", pcol1="12", pcol2="6",
         geneid1="ENSG00000136717", geneidcol1="1",
         outdir=".")

scanRegion = function(region, file, colnames, selectGeneid = NA) {
  s1 = strsplit(region, ":")[[1]]
  chr = s1[[1]]
  s2 = strsplit(s1[[2]], "-")[[1]]
  start = as.integer(s2[1])
  end = as.integer(s2[2])
  
  gr = GRanges(seqnames = chr, ranges = IRanges(start, end=end))
  df = scanTabixDataFrame(file, param = gr, col_names=colnames)[[1]]
  if (!is.na(selectGeneid)) {
    assertthat::assert_that(assertthat::has_name(df, "geneid"))
    df = df %>% dplyr::filter(geneid == selectGeneid)
  }
  df
}

csvToList = function(csvVals) {
  strsplit(csvVals, ",")[[1]]
}

opt$snpcol1 = as.integer(opt$snpcol1)
opt$snpcol2 = as.integer(opt$snpcol2)
opt$pcol1 = as.integer(opt$pcol1)
opt$pcol2 = as.integer(opt$pcol2)
opt$geneidcol1 = as.integer(opt$geneidcol1)

qtl_df = scanRegion(opt$region1, opt$file1, csvToList(opt$file1colnames), opt$geneid1)
qtlSnpsMissing = sum(!(qtl_df$snp_id %in% gwas_df$snp_id))
if (qtlSnpsMissing > 0) {
  warning(paste0(qtlSnpsMissing, " QTL SNPs missing from GWAS summary stats."))
}

gwas_df = scanRegion(opt$region2, opt$file2, csvToList(opt$file2colnames))
gwasSnpsMissing = sum(!(gwas_df$snp_id %in% qtl_df$snp_id))
if (gwasSnpsMissing > 0) {
  warning(paste0(gwasSnpsMissing, " GWAS SNPs missing from QTL summary stats."))
}

merged_df = qtl_df %>% dplyr::select(snp_id, snp_pos1, qtl_p) %>%
  dplyr::full_join(gwas_df %>% dplyr::select(snp_id, snp_pos2, gwas_p)) %>%
  dplyr::arrange(snp_pos1)
sum(merged_df$snp_pos1 == merged_df$snp_pos2, na.rm=T)

if (opt$matchvariants) {
  #GRCh38_variants = importVariantInformation(file.path(root, "macrophage/genotypes/imputed.86_samples.variant_information.txt.gz"))
  GRCh37_variants = importVariantInformation(file.path(root, "macrophage/genotypes/GRCh37/imputed.86_samples.variant_information.GRCh37.vcf.gz"))
  numStartingSnps = nrow(merged_df)
  merged_df = merged_df %>%
    dplyr::left_join(GRCh37_variants %>% dplyr::select(snp_id, pos, ref, alt, MAF))
  merged_df$pos[is.na(merged_df$pos)] = merged_df$snp_pos2[is.na(merged_df$pos)]
  merged_df$pos[is.na(merged_df$pos)] = merged_df$snp_pos2[is.na(merged_df$pos)]
}

plotAssocs = function(df, title1, title2) {
  qtl_df = df %>% dplyr::filter(!is.na(gwas_p))
  qtl_df_missing = df %>% dplyr::filter(is.na(gwas_p))
  p1 = ggplot(qtl_df, aes(x=pos, y=-log10(qtl_p))) +
    geom_point() + geom_point(data=qtl_df_missing, mapping=aes(x=pos, y=-log10(qtl_p)), color="red") +
    theme_bw() + ggtitle(title1) + scale_colour_discrete(guide=F)
  
  #qtl_df = df %>% dplyr::filter(!is.na(gwas_p))
  #qtl_df_missing = df %>% dplyr::filter(is.na(gwas_p))
  gwas_df = df %>% dplyr::filter(!is.na(qtl_p))
  gwas_df_missing = df %>% dplyr::filter(is.na(qtl_p))
  p2 = ggplot(gwas_df, aes(x=pos, y=-log10(gwas_p))) +
  geom_point() + geom_point(data=gwas_df_missing, mapping=aes(x=pos, y=-log10(gwas_p)), color="red") +
    theme_bw() + ggtitle(title2) + scale_colour_discrete(guide=F)
  
  plot_grid(plotlist=list(p1, p2), ncol=1, align='v')
}

plotAssocs(merged_df, title1 = paste(opt$geneid1, "QTL assoc"), title2 = "AD GWAS")



