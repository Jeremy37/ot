library("tidyverse")
library("DESeq2")
library("pcaMethods")
options(stringsAsFactors = F)

jsdir = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys"
jsdir = "/Users/jeremys/work/opentargets"
outdir = file.path(jsdir, "experiment/macrophage/PTK2B_DE")
setwd(outdir)

fc_table = read.table(file.path(outdir, "table featurecounts rnaseq macro ptk2b.txt.gz"), header = T, stringsAsFactors = F)
counts_macro_in <- fc_table[,c(7:12)]
counts_macro_in$geneid = rownames(fc_table)

# Read in data for Kaur's samples. They are named like "aipt_A", and we remove the "_A".
kaur_meta = read.table(file.path(jsdir, "macrophage", "RNA_sample_metadata.txt.gz"), header = T, stringsAsFactors = F) %>%
  dplyr::filter(condition_name == "naive")
kaur_counts = read.table(file.path(jsdir, "macrophage", "RNA_count_matrix.txt.gz"), header = T) %>%
  dplyr::select(one_of(kaur_meta$sample_id))
# Remove the "_A" after each sample ID (this indicated naive state)
kaur_sampleids = sapply(colnames(kaur_counts), FUN=function(s) strsplit(s, "_")[[1]][1])
colnames(kaur_counts) = kaur_sampleids
kaur_meta$sample_id = kaur_sampleids
kaur_counts$geneid = rownames(kaur_counts)

# Read in the VCF file for rs28834970 and parse it to get PTK2B genotypes
vcf_file = file.path(jsdir, "macrophage", "genotype", "imputed.86_samples.rs28834970.vcf")
kaur_vcf = read.table(vcf_file, header = T, comment.char = "") %>%
  dplyr::select(starts_with("HPSI"))
gtfield = sapply(t(kaur_vcf[1,]), FUN=function(s) strsplit(s, ":")[[1]][1])

genotypes.df = data.frame(sample_id = sapply(colnames(kaur_vcf), FUN=function(s) strsplit(s, "\\.|_")[[1]][2]),
                          genotype = sapply(gtfield, FUN=function(s) sum(as.integer(strsplit(s, "|", fixed=T)[[1]]))) )
gt_strings.df = data.frame(genotype = c(0, 1, 2), genotype_PTK2B = c("WT", "HET", "HOM"))
genotypes.df = genotypes.df %>% dplyr::left_join(gt_strings.df, by="genotype")
kaur_meta = kaur_meta %>% dplyr::left_join(genotypes.df, by="sample_id")

allcounts.df = counts_macro_in %>% dplyr::inner_join(kaur_counts, by="geneid")
counts_macro = allcounts.df %>% dplyr::select(ends_with("macro"))

macro.coldata = data.frame(row.names = colnames(counts_macro), sample = colnames(counts_macro), genotype_PTK2B = as.factor(rep(c("WT","HOM"), each=3)), pointsize = 3)
kaur.coldata = data.frame(row.names = kaur_meta$sample_id, sample = kaur_meta$sample_id, genotype_PTK2B = kaur_meta$genotype_PTK2B, pointsize = 2)
all.coldata = rbind(macro.coldata, kaur.coldata)

dds <- DESeqDataSetFromMatrix(countData = allcounts.df %>% dplyr::select(-geneid),
                              colData = all.coldata,
                              design = ~ genotype_PTK2B)
dds

#R log transformation to visualise
# Actually use variance stabilizing transform - faster when more samples
# Do it with all samples together
counts.vst <- vst(dds, blind=TRUE)
#View(assay(counts.vst))
colData(counts.vst)

counts.vst.erica = counts.vst[,1:6]

#PCA plot - first just of Erica's 6 samples
# p <- plotPCA(counts.vst.erica, intgroup=c("genotype_PTK2B"))
# p <- p + geom_text(aes_string(label = "name"), color = "black")
# #p <- p + ylim(-8,8) + xlim(-20,20)
# print(p)
# Plot looks the same as if we had done Rlog transform with just those 6 samples.

N_TOP_GENES = 500
# calculate the variance for each gene
rv = rowVars(assay(counts.vst))

# select the ntop genes by variance
select = order(rv, decreasing=TRUE)[seq_len(min(N_TOP_GENES, length(rv)))]

# First do the PCA on Erica's samples and plot it
pcares.macro = pcaMethods::pca(t(assay(counts.vst.erica)[select,]), nPcs = 5)
pca.df = cbind(scores(pcares.macro), macro.coldata, sample = rownames(scores(pcares.macro)))

pdf(file = "erica_ptk2b_pca_projection.pdf", width=9, height=7)

ggplot(pca.df, aes(x=PC1, y=PC2)) + geom_point(aes(col=genotype_PTK2B, size=pointsize)) +
  geom_text(aes(label = sample), size=3.5, alpha=0.8) +
  xlab(paste0("PC1: ", round(pcares.macro@R2[1] * 100),"% variance")) +
  ylab(paste0("PC2: ", round(pcares.macro@R2[2] * 100),"% variance")) +
  scale_size_continuous(range=c(2,3), guide=F) +
  scale_x_continuous(expand = c(0.1, 0.1)) +
  ggtitle("PCA on just Erica's samples") +
  theme_bw()

# Now do the PCA with all samples to see what we get
pcares.all = pcaMethods::pca(t(assay(counts.vst)[select,]), nPcs = 5)
pca.df.all = cbind(scores(pcares.all), all.coldata, sample = rownames(scores(pcares.all)))

ggplot(pca.df.all, aes(x=PC1, y=PC2)) + geom_point(aes(col=genotype_PTK2B, size=pointsize)) +
  geom_text(aes(label = sample), size=3.5, alpha=0.8) +
  xlab(paste0("PC1: ", round(pcares.all@R2[1] * 100),"% variance")) +
  ylab(paste0("PC2: ", round(pcares.all@R2[2] * 100),"% variance")) +
  scale_size_continuous(range=c(2,3), guide=F) +
  scale_x_continuous(expand = c(0.1, 0.1)) +
  ggtitle("PCA on all samples") +
  theme_bw()


# Now project Kaur's samples onto the axes defined by PCA on Erica's samples
erica.mat.norm = assay(counts.vst.erica)[select,] - rowMeans(assay(counts.vst.erica)[select,])
erica.loadings = t(erica.mat.norm) %*% pcares.macro@loadings
erica.loadings.df = cbind(erica.loadings, macro.coldata)

# ggplot(erica.loadings.df, aes(x=PC1, y=PC2)) + geom_point(aes(col=genotype_PTK2B, size=pointsize)) +
#   geom_text(aes(label = sample), size=3.5, alpha=0.8) +
#   xlab(paste0("PC1: ", round(pcares@R2[1] * 100),"% variance")) +
#   ylab(paste0("PC2: ", round(pcares@R2[2] * 100),"% variance")) +
#   scale_size_continuous(range=c(2,3), guide=F) +
#   scale_x_continuous(expand = c(0.1, 0.1)) +
#   ggtitle("PCA on Erica's samples - projected using her PCA axes") +
#   theme_bw()
# Looks the same as the PCA done just on Erica's samples, and so just confirms
# that we are multiplying by the PC loadings properly

all.mat.norm = assay(counts.vst)[select,] - rowMeans(assay(counts.vst.erica)[select,])
all.loadings = t(all.mat.norm) %*% pcares.macro@loadings
all.loadings.df = cbind(all.loadings, all.coldata)
ggplot(all.loadings.df, aes(x=PC1, y=PC2)) + geom_point(aes(col=genotype_PTK2B, size=pointsize)) +
  geom_text(aes(label = sample), size=3.5, alpha=0.8) +
  scale_size_continuous(range=c(2,3), guide=F) +
  scale_x_continuous(expand = c(0.1, 0.1)) +
  ggtitle("All samples projected onto PCA of Erica's samples") +
  theme_bw()

dev.off()

