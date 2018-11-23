library(tidyverse)

root="/Users/jeremys/work/opentargets/AD_finemap/coloc/"

old_coloc_fname = file.path(root, "output.old", "coloc.AD.meta.eqtl_sqtl_colocs.txt")
new_coloc_fname = file.path(root, "output", "coloc.AD.meta.eqtl_sqtl_colocs.txt")

old.df = readr::read_tsv(old_coloc_fname)
new.df = readr::read_tsv(new_coloc_fname)

cmp.df = new.df %>%
  dplyr::select(feature, dataset_short, locus_name, geneSymbol, H4.corrected = H4) %>%
  dplyr::left_join(old.df %>% dplyr::select(feature, dataset_short, locus_name, H4.old = H4),
                   by=c("feature", "dataset_short", "locus_name"))

ggplot(cmp.df, aes(x=H4.old, y=H4.corrected)) + geom_point(alpha = 0.4) +
  theme_bw()

ggplot(cmp.df, aes(x=H4.old, y=H4.corrected)) + geom_point(alpha = 0.4) +
  theme_bw() +
  coord_cartesian(xlim=c(0.8, 1), ylim=c(0.8, 1))

# Stratify by dataset
ggplot(cmp.df %>% dplyr::filter(!grepl("gtex", dataset_short)), aes(x=H4.old, y=H4.corrected)) +
  geom_point(alpha = 0.4) + theme_bw() +
  facet_wrap(~dataset_short)


cmp.dens.df = cmp.df %>% tidyr::gather(key="version", value="H4", -feature, -dataset_short, -locus_name, -geneSymbol)
ggplot(cmp.dens.df, aes(x=H4, fill=version)) + geom_density(alpha=0.6) + theme_bw()
ggplot(cmp.dens.df, aes(x=H4, fill=version)) + geom_density(alpha=0.6) + theme_bw() +
  coord_cartesian(ylim = c(0,2))

sum(cmp.df$H4.old > 0.9)
sum(cmp.df$H4.corrected > 0.9)
