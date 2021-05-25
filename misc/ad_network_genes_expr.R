library(tidyverse)

ad_dir = "/Users/jeremys/work/opentargets/AD_finemap"
gene_expr_file = paste0(ad_dir, "/reference/tissues.selected.tpm.tsv.gz")

geneExpr.df = readr::read_tsv(gene_expr_file)
network.top = read_tsv(file.path(ad_dir, "network/network.annotated.highly_ranked.tsv"))

network.top = network.top %>%
  left_join(geneExpr.df %>% select(ensgene, primary_microglia_tpm = primary_microglia, ipsc_microglia_tpm = ipsc_microglia, ROSMAP_brain_tpm = ROSMAP_brain),
            by = c("ENSG" = "ensgene"))

roundToString = function(x) {
  if (is.na(x)) { return("") }
  sprintf("%.3g", x)
}

network.numeric = apply(network.top %>% select(minp, minp.quantile, page.rank, page.rank.quantile, Zsco.page.rank.node, Zsco.page.rank.quantile, rankingIte1000.node, primary_microglia_tpm, ipsc_microglia_tpm, ROSMAP_brain_tpm), MARGIN=c(1, 2), FUN=roundToString)

network.top.save = cbind(network.top %>% select(ENSG, gene, biotype, description, degree, chr, start, end),
                        network.numeric)
write_tsv(network.top.save, path = file.path(ad_dir, "network/network.annotated.highly_ranked.with_expression.tsv"), na="")
