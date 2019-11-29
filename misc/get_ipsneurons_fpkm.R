library(tidyverse)
library(annotables)

dir = "/Users/jeremys/work/opentargets/ipsneurons/GRCh38/RNA"
neurons.df = readr::read_tsv(file.path(dir, "analysis/ipsneurons.fpkm.txt.gz"))

meta = readr::read_tsv(file.path(dir, "sample_id_map.txt"), col_names = c("id", "sample_name")) %>%
  mutate(id = tolower(id))

grch38_ensg_sym = grch38 %>% filter(!duplicated(ensgene)) %>% select(gene_id=ensgene, symbol)
neurons.df = neurons.df %>% left_join(grch38_ensg_sym, by="gene_id")

meta.sel = meta %>% filter(id %in% colnames(neurons.df))

neurons.df = neurons.df %>% select(gene_id, symbol, one_of(meta$id))
rownames(meta.sel) = meta.sel$id
colnames(neurons.df)[3:ncol(neurons.df)] = meta.sel[colnames(neurons.df)[3:ncol(neurons.df)],]$sample_name


counts.fname = file.path("/Users/jeremys/work/opentargets/experiment/RNA/CCDC6_clones_new/CCDC6_clones_all.counts.tsv.gz")
readCounts = function(fname) {
  counts = readtsv(fname)
  # Fix gene names to remove the version number (e.g. ".4" in ENSG00000223972.4)
  gene_ids <- gsub("\\.[\\d]+", "", counts[,1], perl=T)
  counts[,1] = gene_ids
  rownames(counts) = gene_ids
  counts
}
counts.df = readCounts(counts.fname) %>%
  select(gene_id, length, starts_with("WT_"))
gene.meta = counts.df[,c(1,2)]
rpm = apply(counts.df[,-c(1,2)], MARGIN=2, FUN=function(x) x / sum(x)) * 1e6
rpkm.mat = apply(rpm * 1e3, MARGIN=2, FUN=function(x) x / gene.meta$length)
rpkm.df = rpkm.mat %>% as.data.frame() %>% rownames_to_column(var="gene_id")

all.df = neurons.df %>%
  left_join(rpkm.df, by="gene_id") %>%
  group_by(gene_id) %>%
  mutate(BOB_iN_D9_avg = mean(BOB_iN_D9_1, BOB_iN_D9_2, BOB_iN_D9_3),
         BOB_iN_D11_avg = mean(BOB_iN_D11_1, BOB_iN_D11_2, BOB_iN_D11_3),
         Kolf2_NPC_avg = mean(Kolf2_NPC_1, Kolf2_NPC_2, Kolf2_NPC_3),
         Kolf2_Neu_D35_avg = mean(Kolf2_Neu_D35_1, Kolf2_Neu_D35_2, Kolf2_Neu_D35_3),
         Kolf2_iPSC_avg = mean(WT_WT1, WT_WT2, WT_KOLF2, WT_KOLF2_1, WT_KOLF2_2, WT_KOLF2_3)) %>%
  select(gene_id, symbol, ends_with("avg"), everything())

write.table(all.df, file = gzfile(file.path(dir, "ipsneurons.rpkm.all.tsv.gz")), quote=F, col.names=T, row.names=F, sep="\t")

write.table(all.df %>% select(gene_id, symbol, Kolf2_iPSC_avg, ends_with("avg")),
            file = gzfile(file.path(dir, "ipsneurons.rpkm.avgs.tsv.gz")), quote=F, col.names=T, row.names=F, sep="\t")

