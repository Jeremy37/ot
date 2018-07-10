
options(stringsAsFactors=F)

df = read.delim("/Users/jeremys/work/opentargets/gwas/AD/finemap.num_variants.txt") %>%
  dplyr::arrange(cs95_variants)
df$locus = factor(as.character(df$locus), levels=as.character(df$locus))
ggplot(df, aes(x=locus, y=cs95_variants)) + geom_bar(stat="identity") +
  theme_bw(16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Locus") + ylab("Variants in\n95% credible set")
