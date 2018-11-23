#!/usr/bin/env Rscript
library(tidyverse)

# root = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys"
root = "/Users/jeremys/work/opentargets"
coloc_root = paste0(root, "/AD_finemap/coloc")

colocs.df = readr::read_tsv(file.path(coloc_root, "output", "coloc.AD.meta.eqtl_sqtl_colocs.txt")) %>%
  dplyr::arrange(-H4) %>%
  dplyr::filter(!grepl("mac_txrevise", dataset_short))

mg.colocs.df = colocs.df %>% dplyr::filter(dataset_short == "microglia") %>% dplyr::arrange(locus, -H4) %>%
  dplyr::select(feature, nsnps, chr, gwas_pos, qtl_pval, gwas_pval, qtl_lead, H4, geneSymbol, locus_name)
#View(mg.colocs.df)
mg.sig_colocs.df = mg.colocs.df %>% dplyr::filter(H4 >= 0.9)
mg.top_colocs.df = mg.colocs.df %>% dplyr::filter(H4 >= 0.5, !duplicated(locus_name))

num_colocs.df = colocs.df %>% group_by(dataset_short) %>%
  summarise(num_coloc_tests = n(),
            num_colocs0.5 = sum(H4 >= 0.5),
            num_colocs0.9 = sum(H4 >= 0.9),
            num_PP3gt0.5 = sum(PP.H3 >= 0.5),
            num_PP3gt0.9 = sum(PP.H3 >= 0.9))

num_egenes.df = readr::read_tsv(file.path(coloc_root, "qtl_data", "qtl_num_egenes.tsv"))

num_colocs.df = num_colocs.df %>% dplyr::left_join(num_egenes.df, by="dataset_short")

getPlotDF = function(num_colocs.df) {
  num_colocs.df$is_gtex = grepl("gtex", num_colocs.df$dataset_short)
  
  num_colocs.df$color = 0
  num_colocs.df[num_colocs.df$dataset_short == "microglia",]$color = 1
  num_colocs.df[num_colocs.df$is_gtex,]$color = 2
  num_colocs.df$color = factor(num_colocs.df$color, levels=c(0, 1, 2))
  
  num_colocs.df$text_color = 0
  num_colocs.df[num_colocs.df$dataset_short == "microglia",]$text_color = 1
  num_colocs.df[num_colocs.df$is_gtex,]$text_color = 2
  num_colocs.df$text_color = factor(num_colocs.df$text_color, levels=c(0, 1, 2))
  
  num_colocs.df$size = 0
  num_colocs.df[num_colocs.df$dataset_short == "microglia",]$size = 1
  num_colocs.df$size = factor(num_colocs.df$size, levels=c(0, 1))
  num_colocs.df
}

plot.df = getPlotDF(num_colocs.df)
labelDatasets = c("microglia", "mono_eQTL", "xQTL_eQTL", "mac_eqtl_naive", "mac_eqtl_IFNg_SL1344", "mac_sqtl_naive", "sensneur_eqtl",
                  "gtex.brain_cereb", "gtex.brain_hippo", "gtex.brain_cortex", "gtex.nerve_tibial", "gtex.cells_lcl")
labels.df = plot.df %>% dplyr::filter(dataset_short %in% labelDatasets)

fontSize = 12
pdf(file = file.path(coloc_root, "AD.coloc.comparison.pdf"), width = 7, height = 5)

###############################################################################
# Look at the number of colocs at different thresholds vs. the number of
# eGenes or tests across each dataset
ggplot(plot.df, aes(x=num_egenes, y=num_coloc_tests, label=dataset_short, col=color)) + geom_point() + theme_bw(fontSize) +
  geom_text(mapping = aes(size=size, col=text_color), data = labels.df, hjust=0, vjust=1) +
  scale_size_manual(values=c(3, 5), guide=F) + scale_color_manual(values=c("black", "blue", "grey40"), guide=F) + scale_x_continuous(expand=c(0.1, 0)) +
  ggtitle("Num coloc tests vs. Num eGenes") +
  ylab("Num coloc tests") + xlab("Num eGenes")

ggplot(plot.df, aes(x=num_egenes, y=num_colocs0.5, label=dataset_short, col=color)) + geom_point() + theme_bw(fontSize) +
  geom_text(mapping = aes(size=size, col=text_color), data = labels.df, hjust=0, vjust=1) +
  scale_size_manual(values=c(3, 5), guide=F) + scale_color_manual(values=c("black", "blue", "grey40"), guide=F) + scale_x_continuous(expand=c(0.1, 0)) +
  ggtitle("Num colocs H4 > 0.5 vs. Num eGenes") +
  ylab("Num colocs H4 > 0.5") + xlab("Num eGenes")

ggplot(plot.df, aes(x=num_egenes, y=num_colocs0.9, label=dataset_short, col=color)) + geom_point() + theme_bw(fontSize) +
  geom_text(mapping = aes(size=size, col=text_color), data = labels.df, hjust=0, vjust=1) +
  scale_size_manual(values=c(3, 5), guide=F) + scale_color_manual(values=c("black", "blue", "grey40"), guide=F) + scale_x_continuous(expand=c(0.1, 0)) +
  ggtitle("Num colocs H4 > 0.9 vs. Num eGenes") +
  ylab("Num colocs H4 > 0.9") + xlab("Num eGenes")

ggplot(plot.df, aes(x=num_coloc_tests, y=num_colocs0.5, label=dataset_short, col=color)) + geom_point() + theme_bw(fontSize) +
  geom_text(mapping = aes(size=size, col=text_color), data = labels.df, hjust=0, vjust=1) +
  scale_size_manual(values=c(3, 5), guide=F) + scale_color_manual(values=c("black", "blue", "grey40"), guide=F) + scale_x_continuous(expand=c(0.1, 0)) +
  ggtitle("Num colocs H4 > 0.5 vs. Num coloc tests") +
  ylab("Num colocs H4 > 0.5") + xlab("Num coloc tests")

ggplot(plot.df, aes(x=num_coloc_tests, y=num_colocs0.9, label=dataset_short, col=color)) + geom_point() + theme_bw(fontSize) +
  geom_text(mapping = aes(size=size, col=text_color), data = labels.df, hjust=0, vjust=1) +
  scale_size_manual(values=c(3, 5), guide=F) + scale_color_manual(values=c("black", "blue", "grey40"), guide=F) + scale_x_continuous(expand=c(0.1, 0)) +
  ggtitle("Num colocs H4 > 0.9 vs. Num coloc tests") +
  ylab("Num colocs H4 > 0.9") + xlab("Num coloc tests")


ggplot(plot.df, aes(x=num_PP3gt0.5, y=num_colocs0.5, label=dataset_short, col=color)) + geom_point() + theme_bw(fontSize) +
  geom_text(mapping = aes(size=size, col=text_color), data = labels.df, hjust=0, vjust=1) +
  scale_size_manual(values=c(3, 5), guide=F) + scale_color_manual(values=c("black", "blue", "grey40"), guide=F) + scale_x_continuous(expand=c(0.1, 0)) +
  ggtitle("Num colocs H4 > 0.5 vs. Num indep H3 > 0.5") +
  ylab("Num colocs H4 > 0.5") + xlab("Num indep H3 > 0.5")

ggplot(plot.df, aes(x=num_PP3gt0.9, y=num_colocs0.9, label=dataset_short, col=color)) + geom_point() + theme_bw(fontSize) +
  geom_text(mapping = aes(size=size, col=text_color), data = labels.df, hjust=0, vjust=1) +
  scale_size_manual(values=c(3, 5), guide=F) + scale_color_manual(values=c("black", "blue", "grey40"), guide=F) + scale_x_continuous(expand=c(0.1, 0)) +
  ggtitle("Num colocs H4 > 0.9 vs. Num indep H3 > 0.9") +
  ylab("Num colocs H4 > 0.9") + xlab("Num indep H3 > 0.9")


###############################################################################
# Number of colocs at different thresholds when selecting the single top coloc
# gene per locus
top_colocs.df = colocs.df %>% dplyr::arrange(dataset_short, locus_name, -H4) %>%
  dplyr::filter(!duplicated(paste(dataset_short, locus_name)))
num_colocs.top.df = top_colocs.df %>%
  group_by(dataset_short) %>%
  summarise(num_loci_with_egene = n(),
            num_loci_coloc_0.5 = sum(H4 >= 0.5),
            num_loci_coloc_0.9 = sum(H4 >= 0.9))
num_colocs.top.df = num_colocs.top.df %>%
  dplyr::left_join(num_colocs.df %>% dplyr::select(dataset_short, num_coloc_tests, num_egenes), by="dataset_short")

plot.df = getPlotDF(num_colocs.top.df)
labelDatasets = c("microglia", "mono_eQTL", "xQTL_eQTL", "mac_eqtl_naive", "mac_eqtl_IFNg_SL1344", "mac_sqtl_naive", "sensneur_eqtl",
                  "gtex.brain_cereb", "gtex.brain_hippo", "gtex.brain_cortex", "gtex.nerve_tibial", "gtex.cells_lcl")
labels.df = plot.df %>% dplyr::filter(dataset_short %in% labelDatasets)

ggplot(plot.df, aes(x=num_egenes, y=num_loci_with_egene, label=dataset_short, col=color)) + geom_point() + theme_bw(fontSize) +
  geom_text(mapping = aes(size=size, col=text_color), data = labels.df, hjust=0, vjust=1) +
  scale_size_manual(values=c(3, 5), guide=F) + scale_color_manual(values=c("black", "blue", "grey40"), guide=F) + scale_x_continuous(expand=c(0.1, 0)) +
  ggtitle("Num AD loci with an eGene vs. Num eGenes") +
  ylab("Num AD loci with an eGene") + xlab("Num eGenes")

ggplot(plot.df, aes(x=num_coloc_tests, y=num_loci_with_egene, label=dataset_short, col=color)) + geom_point() + theme_bw(fontSize) +
  geom_text(mapping = aes(size=size, col=text_color), data = labels.df, hjust=0, vjust=1) +
  scale_size_manual(values=c(3, 5), guide=F) + scale_color_manual(values=c("black", "blue", "grey40"), guide=F) + scale_x_continuous(expand=c(0.1, 0)) +
  ggtitle("Num AD loci with an eGene vs. Num colocs tests") +
  ylab("Num AD loci with an eGene") + xlab("Num colocs tests")

ggplot(plot.df, aes(x=num_egenes, y=num_loci_coloc_0.5, label=dataset_short, col=color)) + geom_point() + theme_bw(fontSize) +
  geom_text(mapping = aes(size=size, col=text_color), data = labels.df, hjust=0, vjust=1) +
  scale_size_manual(values=c(3, 5), guide=F) + scale_color_manual(values=c("black", "blue", "grey40"), guide=F) + scale_x_continuous(expand=c(0.1, 0)) +
  ggtitle("Num loci with at least one H4 > 0.5 vs. Num eGenes") +
  ylab("Num AD loci with an eGene H4 > 0.5") + xlab("Num eGenes")

ggplot(plot.df, aes(x=num_egenes, y=num_loci_coloc_0.9, label=dataset_short, col=color)) + geom_point() + theme_bw(fontSize) +
  geom_text(mapping = aes(size=size, col=text_color), data = labels.df, hjust=0, vjust=1) +
  scale_size_manual(values=c(3, 5), guide=F) + scale_color_manual(values=c("black", "blue", "grey40"), guide=F) + scale_x_continuous(expand=c(0.1, 0)) +
  ggtitle("Num loci with at least one H4 > 0.9 vs. Num eGenes") +
  ylab("Num AD loci with an eGene H4 > 0.9") + xlab("Num eGenes")

ggplot(plot.df, aes(x=num_coloc_tests, y=num_loci_coloc_0.5, label=dataset_short, col=color)) + geom_point() + theme_bw(fontSize) +
  geom_text(mapping = aes(size=size, col=text_color), data = labels.df, hjust=0, vjust=1) +
  scale_size_manual(values=c(3, 5), guide=F) + scale_color_manual(values=c("black", "blue", "grey40"), guide=F) + scale_x_continuous(expand=c(0.1, 0)) +
  ggtitle("Num loci with at least one H4 > 0.5 vs. Num coloc tests") +
  ylab("Num AD loci with an eGene H4 > 0.5") + xlab("Num coloc tests")

ggplot(plot.df, aes(x=num_coloc_tests, y=num_loci_coloc_0.9, label=dataset_short, col=color)) + geom_point() + theme_bw(fontSize) +
  geom_text(mapping = aes(size=size, col=text_color), data = labels.df, hjust=0, vjust=1) +
  scale_size_manual(values=c(3, 5), guide=F) + scale_color_manual(values=c("black", "blue", "grey40"), guide=F) + scale_x_continuous(expand=c(0.1, 0)) +
  ggtitle("Num loci with at least one H4 > 0.9 vs. Num coloc tests") +
  ylab("Num AD loci with an eGene H4 > 0.9") + xlab("Num coloc tests")

ggplot(plot.df, aes(x=num_loci_with_egene, y=num_loci_coloc_0.5, label=dataset_short, col=color)) + geom_point() + theme_bw(fontSize) +
  geom_text(mapping = aes(size=size, col=text_color), data = labels.df, hjust=0, vjust=1) +
  scale_size_manual(values=c(3, 5), guide=F) + scale_color_manual(values=c("black", "blue", "grey40"), guide=F) + scale_x_continuous(expand=c(0.1, 0)) +
  ggtitle("Num loci with at least one H4 > 0.5 vs. Num AD loci with an eGene") +
  ylab("Num AD loci with an eGene H4 > 0.5") + xlab("Num AD loci with an eGene")

ggplot(plot.df, aes(x=num_loci_with_egene, y=num_loci_coloc_0.9, label=dataset_short, col=color)) + geom_point() + theme_bw(fontSize) +
  geom_text(mapping = aes(size=size, col=text_color), data = labels.df, hjust=0, vjust=1) +
  scale_size_manual(values=c(3, 5), guide=F) + scale_color_manual(values=c("black", "blue", "grey40"), guide=F) + scale_x_continuous(expand=c(0.1, 0)) +
  ggtitle("Num loci with at least one H4 > 0.9 vs. Num AD loci with an eGene") +
  ylab("Num AD loci with an eGene H4 > 0.9") + xlab("Num AD loci with an eGene")
  

###############################################################################
# Expression of genes with coloc across datasets

# First make sure that we have an ensembl ID for each coloc gene
hgnc.ensgene.df = readr::read_tsv(file.path(root, "reference", "GRCh38", "GRCh38.p12.geneid_hgnc_length.txt.gz")) %>%
  dplyr::filter(grepl("^[0-9]+$|^X|^Y|^MT", chr)) %>%
  dplyr::select(gene_id, hgnc_symbol) %>% dplyr::distinct()

colocs.df = colocs.df %>% dplyr::left_join(hgnc.ensgene.df %>% dplyr::rename(matching_ens_id = gene_id), by=c("geneSymbol" = "hgnc_symbol"))
sum(is.na(colocs.df$ensembl_id))
colocs.df$ensembl_id[is.na(colocs.df$ensembl_id)] = colocs.df$matching_ens_id[is.na(colocs.df$ensembl_id)]
sum(is.na(colocs.df$ensembl_id))

# Get ranks of each gene's expression in each dataset
expr_ranks.df = readr::read_tsv(file.path(root, "reference", "tissueRPKM", "tissues.selected.expr_ranks.txt.gz")) %>%
  dplyr::select(-symbol, -other_symbols, -barres_symbol, -starts_with("Barres"),
                -NPC, -NPC_B, -neuron, -iNeuron_d9, -iNeuron_d11, -NPC_PSEN, -iNeuron_d10_vit)
# Convert ranks to percentiles, and convert NAs into the bottom percentile of
# expression
#expr_ranks.df = expr_ranks.df %>% na.omit()
expr_pct.df = cbind(list(gene_id = expr_ranks.df$gene_id),
                    as.data.frame(lapply(expr_ranks.df[, -1], FUN=function(ranks) ranks / sum(!is.na(ranks)))))
colnames(expr_pct.df) = colnames(expr_ranks.df)
expr_pct.df[is.na(expr_pct.df)] = 1

# Gather the expression data into a pair of columns for joining
expr_pct.df = tidyr::gather(expr_pct.df, key = "dataset", value = "expr_quantile", -gene_id)

expr_dataset_map.df = readr::read_tsv(file.path(coloc_root, "expr_dataset_map.txt")) %>%
  dplyr::select(qtl_dataset_short, expr_dataset)

colocs.df = colocs.df %>% dplyr::left_join(expr_dataset_map.df, by=c("dataset_short" = "qtl_dataset_short"))
colocs.df = colocs.df %>% dplyr::left_join(expr_pct.df, by=c("expr_dataset" = "dataset", "ensembl_id" = "gene_id"))
sum(is.na(colocs.df$expr_quantile))

# A few of the colocs have multiple genes listed in ensembl_id (e.g. when a
# splice QTL overlaps more than one gene), and so these won't have matched
# up their gene expression. For these, we take the max expr of any of the
# listed genes.
colocs.multigene.df = colocs.df[grepl(",", colocs.df$ensembl_id),]
multigene.df = data.frame()
for (i in 1:nrow(colocs.multigene.df)) {
  ids = strsplit(colocs.multigene.df[i,]$ensembl_id, ",", fixed = T)[[1]]
  multigene.df = rbind(multigene.df, data.frame(gene_id = ids,
                                                dataset_short = colocs.multigene.df[i,]$dataset_short,
                                                ensembl_ids = colocs.multigene.df[i,]$ensembl_id,
                                                expr_dataset = colocs.multigene.df[i,]$expr_dataset))
}
multigene.expr.df = multigene.df %>% dplyr::left_join(expr_pct.df, by=c("expr_dataset" = "dataset", "gene_id"))
multigene.expr.df = multigene.expr.df %>% group_by(dataset_short, ensembl_ids) %>%
  dplyr::summarise(expr_rank_multigene = min(expr_quantile, na.rm=T))

colocs.df = colocs.df %>% dplyr::left_join(multigene.expr.df, by = c("dataset_short", "ensembl_id"="ensembl_ids"))
colocs.df[is.na(colocs.df$expr_quantile),]$expr_quantile = colocs.df$expr_rank_multigene[is.na(colocs.df$expr_quantile)]
sum(is.na(colocs.df$expr_quantile))

save.df = colocs.df %>% dplyr::select(-matching_ens_id, -expr_dataset, -expr_rank_multigene)
write.table(save.df, file=file.path(coloc_root, "output", "coloc.AD.meta.eqtl_sqtl_colocs.expr.txt"),
            col.names=T, row.names=F, quote=F, sep="\t", na="")
#mg.save.df = save.df %>% dplyr::filter(dataset_short == "microglia")
#write.table(mg.save.df, file=file.path(coloc_root, "output", "coloc.AD.meta.microglia_colocs.expr.txt"),
#            col.names=T, row.names=F, quote=F, sep="\t", na="")

exprDatasets = c("microglia", "mac_sqtl_naive", "mac_eqtl_naive", "mac_eqtl_IFNg_SL1344", "xQTL_eQTL", "mono_eQTL", "sensneur_eqtl",
                  "gtex.brain_cereb", "gtex.brain_hippo", "gtex.brain_cortex", "gtex.nerve_tibial", "gtex.cells_lcl")
plot.df = colocs.df %>% dplyr::filter(dataset_short %in% exprDatasets) %>%
  dplyr::arrange(locus, dataset, -H4)
plot.df$dataset_short = factor(plot.df$dataset_short, levels = exprDatasets)
plot.df$H4gt05 = plot.df$H4 > 0.5
plot.df$H4gt09 = plot.df$H4 >= 0.9
plot.df$H4gt09 = as.character(plot.df$H4gt09)

# Plots of expression quantiles for all coloc genes
ggplot(plot.df, aes(x=dataset_short, y=expr_quantile, col=H4gt09, alpha=H4gt09, group=dataset_short)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(data = plot.df %>% dplyr::filter(H4 < 0.9), width = 0.25) +
  geom_jitter(data = plot.df %>% dplyr::filter(H4 >= 0.9), width = 0.25) +
  scale_color_manual(values=c("grey40", "blue"), name="PP H4 > 0.9") +
  scale_alpha_manual(values=c(0.5, 1), guide=F) +
  theme_bw(fontSize) + ggtitle("Expression quantiles of all QTL genes at coloc loci") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Expression quantile") + xlab("QTL Dataset")

# ggplot(plot.df %>% dplyr::filter(H4 > 0.5), aes(x=dataset_short, y=expr_quantile)) + geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = 0.25) +
#   scale_color_manual(values=c("grey60", "blue"), name="PP H4 > 0.9") +
#   theme_bw(fontSize) + ggtitle("Expression rank of genes at coloc loci") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggplot(plot.df %>% dplyr::filter(H4 > 0.9), aes(x=dataset_short, y=expr_quantile)) + geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = 0.25) +
#   scale_color_manual(values=c("grey60", "blue"), name="PP H4 > 0.9") +
#   theme_bw(fontSize) + ggtitle("Expression quantiles of genes at coloc loci") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Filter to include only a single top coloc gene per locus / dataset combination
plot.df2 = plot.df %>% dplyr::filter(!duplicated(paste(dataset_short, locus)))

ggplot(plot.df2, aes(x=dataset_short, y=expr_quantile, col=H4gt09, alpha=H4gt09, group=dataset_short)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(data = plot.df2 %>% dplyr::filter(H4 < 0.9), width = 0.25) +
  geom_jitter(data = plot.df2 %>% dplyr::filter(H4 >= 0.9), width = 0.25) +
  scale_color_manual(values=c("grey40", "blue"), name="PP H4 > 0.9") +
  scale_alpha_manual(values=c(0.5, 1), guide=F) +
  theme_bw(fontSize) + ggtitle("Expression quantiles of top coloc gene per locus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Expression quantile") + xlab("QTL Dataset")

# ggplot(plot.df2 %>% dplyr::filter(H4 > 0.9), aes(x=dataset_short, y=expr_quantile)) + geom_boxplot(outlier.shape = NA) +
#   geom_jitter(data = plot.df2 %>% dplyr::filter(H4 >= 0.9), width = 0.25) +
#   scale_color_manual(values=c("grey60", "blue"), name="PP H4 > 0.9") +
#   theme_bw(fontSize) + ggtitle("Expression quantiles of top coloc gene per locus (only H4 > 0.9)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Plot expression according to the rank of each coloc gene within a dataset
plot.df3 = plot.df2 %>% dplyr::arrange(dataset_short, expr_quantile) %>%
  dplyr::group_by(dataset_short) %>%
  dplyr::mutate(index = row_number())

ggplot(plot.df3, aes(x=index, y=expr_quantile, col=H4gt09)) + geom_point() +
  facet_wrap(~dataset_short) +
  scale_color_manual(values=c("grey60", "blue"), name="PP H4 > 0.9") +
  theme_bw(fontSize) + ggtitle("Expression quantiles of top coloc gene per locus") +
  ylab("Expression quantile") + xlab("AD locus index")

# Subset to genes with a coloc having H4 > 0.5
plot.df3 = plot.df2 %>% dplyr::arrange(dataset_short, expr_quantile) %>%
  dplyr::filter(H4 > 0.5) %>%
  dplyr::group_by(dataset_short) %>%
  dplyr::mutate(index = row_number())

ggplot(plot.df3, aes(x=index, y=expr_quantile, col=H4gt09)) + geom_point() +
  facet_wrap(~dataset_short) +
  scale_color_manual(values=c("grey60", "blue"), name="PP H4 > 0.9") +
  theme_bw(fontSize) + ggtitle("Expression quantiles of top coloc gene per locus (only H4 > 0.5)") +
  ylab("Expression quantile") + xlab("AD locus index")


dev.off()

