#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
genesFile = args[1]
diffexpFile = args[2]
causalRlength1File = args[3]
causalRlength2File = args[4]
causalRlength3File = args[5]
causalRlength4File = args[6]

# root = "/Users/jeremys/work/opentargets/gwas/AD/causalR"
# genesFile = file.path(root, "AD.leadSNPs.1Mb_window.gene_overlaps.txt")
# diffexpFile = file.path(root, "MAYO_ADvCntl_DiffExp_ALL.txt")
# causalRlength1File = file.path(root, "FullCausalRresults_MAYO_ADvCntl.pathlength_1.txt")
# causalRlength2File = file.path(root, "FullCausalRresults_MAYO_ADvCntl.pathlength_2.txt")
# causalRlength3File = file.path(root, "FullCausalRresults_MAYO_ADvCntl.pathlength_3.txt")
# causalRlength4File = file.path(root, "FullCausalRresults_MAYO_ADvCntl.pathlength_4.txt")

genes.df = readr::read_tsv(genesFile)
diffexp.df = readr::read_tsv(diffexpFile) %>% dplyr::rename(symbol=Symbol, geneName=`Entrez Gene Name`, diffexp_log2FC=`Expr Log Ratio`, diffexp_FDR=FDR)
causalRlength1.df = readr::read_tsv(causalRlength1File) %>% dplyr::rename(symbol=NodeName,
                                                                          path1_direction=Regulation,
                                                                          path1_score=Score,
                                                                          path1_correct=Correct,
                                                                          path1_incorrect=Incorrect,
                                                                          path1_ambiguous=Ambiguous,
                                                                          path1_p_val=`p-value`,
                                                                          path1_enrichment_p_val=`Enrichment p-value`)
causalRlength2.df = readr::read_tsv(causalRlength2File) %>% dplyr::rename(symbol=NodeName,
                                                                          path2_direction=Regulation,
                                                                          path2_score=Score,
                                                                          path2_correct=Correct,
                                                                          path2_incorrect=Incorrect,
                                                                          path2_ambiguous=Ambiguous,
                                                                          path2_p_val=`p-value`,
                                                                          path2_enrichment_p_val=`Enrichment p-value`)

causalRlength3.df = readr::read_tsv(causalRlength3File) %>% dplyr::rename(symbol=NodeName,
                                                                          path3_direction=Regulation,
                                                                          path3_score=Score,
                                                                          path3_correct=Correct,
                                                                          path3_incorrect=Incorrect,
                                                                          path3_ambiguous=Ambiguous)
causalRlength4.df = readr::read_tsv(causalRlength4File) %>% dplyr::rename(symbol=NodeName,
                                                                          path4_direction=Regulation,
                                                                          path4_score=Score,
                                                                          path4_correct=Correct,
                                                                          path4_incorrect=Incorrect,
                                                                          path4_ambiguous=Ambiguous)

# Remove duplicate entries from the Differential expression file. It's not clear yet why there
# are multiple entries, Glyn Bradley is looking into it.
diffexp.dedup.df = diffexp.df %>% arrange(diffexp_FDR, diffexp_log2FC) %>% dplyr::filter(!duplicated(symbol))

# Remove duplicate entries from the causalR results. There are two entries for each gene, one
# relating to upregulation and one to downregulation. We take the most significant one.
causalRlength1.df = causalRlength1.df %>%
  arrange(-path1_score, path1_p_val) %>%
  dplyr::filter(!duplicated(symbol)) %>%
  dplyr::mutate(path1_rank = floor(rank(-path1_score)))
causalRlength2.df = causalRlength2.df %>%
  arrange(-path2_score, path2_p_val) %>%
  dplyr::filter(!duplicated(symbol)) %>%
  dplyr::mutate(path2_rank = floor(rank(-path2_score)))

causalRlength3.df = causalRlength3.df %>%
  arrange(-path3_score) %>%
  dplyr::filter(!duplicated(symbol)) %>%
  dplyr::mutate(path3_rank = floor(rank(-path3_score)))
causalRlength4.df = causalRlength4.df %>%
  arrange(-path4_score) %>%
  dplyr::filter(!duplicated(symbol)) %>%
  dplyr::mutate(path4_rank = floor(rank(-path4_score)))

genes.ann.df = genes.df %>%
  dplyr::left_join(diffexp.dedup.df %>% dplyr::select(-ID, -geneName), by="symbol") %>%
  dplyr::left_join(causalRlength1.df, by="symbol") %>%
  dplyr::left_join(causalRlength2.df, by="symbol") %>%
  dplyr::left_join(causalRlength3.df, by="symbol") %>%
  dplyr::left_join(causalRlength4.df, by="symbol")

write.table(genes.ann.df, file="", quote=F, sep="\t", row.names=F, col.names=T, na = "")
