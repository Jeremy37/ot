library(pheatmap)
library(pcaMethods)

# To use these scripts you must also source the pca.projections.R file:
# source("~/js29/src/R/pca.projection.R")

# Function that produces a heatmap comparing samples with reference samples.
#  ref.mat: a gene expression matrix (rpkms) for reference samples, with gene IDs as rownames
#  sample.mat: gene expression matrix (rpkms) for user samples with gene IDs as rownames
#  ref.meta: a data.frame with columns ID (sample IDs) and group (for samples of the same type)
#  sample.meta: a data.frame with columsn ID and group
#  averageGroups: if True, then samples within a group are averaged before plotting
heatmapWithReference = function(ref.mat, ref.meta, sample.mat, sample.meta, averageGroups = T, title=NULL) {
  if (!averageGroups) {
    ref.cmp.df = as.data.frame(ref.mat)
    ref.cmp.df$gene = toupper(rownames(ref.mat))
  } else {
    ref.cmp.df = data.frame(gene = toupper(rownames(ref.mat)))
    for (group in unique(barres.meta$group)) {
      ref.cmp.df[, group] = rowMeans(ref.mat[, ref.meta$group == group, drop=F])
    }
  }
  
  if (ncol(ref.cmp.df) > 50) {
    warning("More than 50 reference samples to compare - heatmap will be very large")
  } else if (ncol(ref.cmp.df) > 100) {
    stop("More than 100 reference samples to compare. Infeasible to produce heatmap. Aborting.")
  }
  
  sample.df = sample.mat %>% as.data.frame()
  sample.df$gene = rownames(sample.df)
  rpkm.cmp = sample.df %>% dplyr::inner_join(ref.cmp.df, by="gene") %>%
    dplyr::select(-gene)
  
  #cor.pearson = cor(log2(rpkm.cmp + 0.1))
  #pheatmap(cor.pearson, main="Heatmap of gene expression correlations (Pearson)", fontsize_row = 9, fontsize_col = 8, display_numbers = T, number_format = "%.2f", fontsize_number=7)
  
  cor.spearman = cor(rpkm.cmp, method="spearman")
  if (is.null(title)) {
    title = "Heatmap of gene expression correlations (Spearman)"
  }
  fontsize_number = 7
  if (ncol(cor.spearman) > 30) {
    fontsize_number == 5
  }
  pheatmap(cor.spearman, main=title, fontsize_row = 9, fontsize_col = 8, display_numbers = T, number_format = "%.2f", fontsize_number=fontsize_number)
}

# Function to project samples onto a PCA defined by reference samples.
#  ref.mat: a gene expression matrix (rpkms) for reference samples, with gene IDs as rownames
#  sample.mat: gene expression matrix (rpkms) for user samples with gene IDs as rownames
#  ref.meta: a data.frame with columns "ID" (sample IDs) and "group" (for samples of the same type)
#  sample.meta: a data.frame with columns "ID" and "group"
#  averageGroups: if True, then samples within a group are averaged before plotting
#  maxPC: the highest PC to be shown in a plot. E.g. if maxPC=3, then the following plots
#         will be shown: PC1xPC2, PC1xPC3, PC2xPC3
doPcaProjection = function(ref.mat, ref.meta, sample.mat, sample.meta, maxPC=2) {
  if (maxPC < 2) {
    return()
  }
  if (ncol(ref.mat) > 500) {
    pcaProj = bigpcaProjection(log2(as.matrix(ref.mat)+1), gtex.pca.meta, log2(as.matrix(sample.mat)+1), sample.meta)
  } else {
    pcaProj = pcaProjection(log2(as.matrix(ref.mat)+1), ref.meta, log2(as.matrix(sample.mat)+1), sample.meta)
  }
  pca.s = summary(pcaProj$pcaRes)
  
  groupLoc = pcaGroupMeans(pcaProj$sampleCoords, "PC1", "PC2")
  groupLoc$group = rownames(groupLoc)
  xRange = max(groupLoc$x) - min(groupLoc$x)
  xlim1 = min(groupLoc$x) - 0.05*xRange
  xlim2 = max(groupLoc$x) + 0.05*xRange
  p = ggplot(pcaProj$sampleCoords, aes(x=PC1, y=PC2, col=group)) +
    geom_point(alpha=0.7) + theme_bw(12) +
    geom_text(mapping=aes(x=x, y=y, label=group, col=group), data=groupLoc, size=3) +
    scale_colour_discrete(guide = FALSE) +
    coord_cartesian(xlim=c(xlim1, xlim2)) +
    xlab(sprintf("PC1 (%.1f%%)", pca.s$importance["Proportion of Variance","PC1"]*100)) +
    ylab(sprintf("PC2 (%.1f%%)", pca.s$importance["Proportion of Variance","PC2"]*100))
  print(p)
  
  if (maxPC >= 3) {
    groupLoc = pcaGroupMeans(pcaProj$sampleCoords, "PC1", "PC3")
    groupLoc$group = rownames(groupLoc)
    xRange = max(groupLoc$x) - min(groupLoc$x)
    xlim1 = min(groupLoc$x) - 0.05*xRange
    xlim2 = max(groupLoc$x) + 0.05*xRange
    p = ggplot(pcaProj$sampleCoords, aes(x=PC1, y=PC3, col=group)) +
      geom_point(alpha=0.7) + theme_bw(12) +
      geom_text(mapping=aes(x=x, y=y, label=group, col=group), data=groupLoc, size=3) +
      scale_colour_discrete(guide = FALSE) +
      coord_cartesian(xlim=c(xlim1, xlim2)) +
      xlab(sprintf("PC1 (%.1f%%)", pca.s$importance["Proportion of Variance","PC1"]*100)) +
      ylab(sprintf("PC3 (%.1f%%)", pca.s$importance["Proportion of Variance","PC3"]*100))
    print(p)
    
    groupLoc = pcaGroupMeans(pcaProj$sampleCoords, "PC2", "PC3")
    groupLoc$group = rownames(groupLoc)
    xRange = max(groupLoc$x) - min(groupLoc$x)
    xlim1 = min(groupLoc$x) - 0.05*xRange
    xlim2 = max(groupLoc$x) + 0.05*xRange
    p = ggplot(pcaProj$sampleCoords, aes(x=PC2, y=PC3, col=group)) +
      geom_point(alpha=0.7) + theme_bw(12) +
      geom_text(mapping=aes(x=x, y=y, label=group, col=group), data=groupLoc, size=3) +
      scale_colour_discrete(guide = FALSE) +
      coord_cartesian(xlim=c(xlim1, xlim2)) +
      xlab(sprintf("PC2 (%.1f%%)", pca.s$importance["Proportion of Variance","PC2"]*100)) +
      ylab(sprintf("PC3 (%.1f%%)", pca.s$importance["Proportion of Variance","PC3"]*100))
    print(p)
  }
  
  if (maxPC >= 4) {
    groupLoc = pcaGroupMeans(pcaProj$sampleCoords, "PC1", "PC4")
    groupLoc$group = rownames(groupLoc)
    xRange = max(groupLoc$x) - min(groupLoc$x)
    xlim1 = min(groupLoc$x) - 0.05*xRange
    xlim2 = max(groupLoc$x) + 0.05*xRange
    p = ggplot(pcaProj$sampleCoords, aes(x=PC1, y=PC4, col=group)) +
      geom_point(alpha=0.7) + theme_bw(12) +
      geom_text(mapping=aes(x=x, y=y, label=group, col=group), data=groupLoc, size=3) +
      scale_colour_discrete(guide = FALSE) +
      coord_cartesian(xlim=c(xlim1, xlim2)) +
      xlab(sprintf("PC1 (%.1f%%)", pca.s$importance["Proportion of Variance","PC1"]*100)) +
      ylab(sprintf("PC4 (%.1f%%)", pca.s$importance["Proportion of Variance","PC4"]*100))
    print(p)
    
    groupLoc = pcaGroupMeans(pcaProj$sampleCoords, "PC2", "PC4")
    groupLoc$group = rownames(groupLoc)
    xRange = max(groupLoc$x) - min(groupLoc$x)
    xlim1 = min(groupLoc$x) - 0.05*xRange
    xlim2 = max(groupLoc$x) + 0.05*xRange
    p = ggplot(pcaProj$sampleCoords, aes(x=PC2, y=PC4, col=group)) +
      geom_point(alpha=0.7) + theme_bw(12) +
      geom_text(mapping=aes(x=x, y=y, label=group, col=group), data=groupLoc, size=3) +
      scale_colour_discrete(guide = FALSE) +
      coord_cartesian(xlim=c(xlim1, xlim2)) +
      xlab(sprintf("PC2 (%.1f%%)", pca.s$importance["Proportion of Variance","PC2"]*100)) +
      ylab(sprintf("PC4 (%.1f%%)", pca.s$importance["Proportion of Variance","PC4"]*100))
    print(p)
    
    groupLoc = pcaGroupMeans(pcaProj$sampleCoords, "PC3", "PC4")
    groupLoc$group = rownames(groupLoc)
    xRange = max(groupLoc$x) - min(groupLoc$x)
    xlim1 = min(groupLoc$x) - 0.05*xRange
    xlim2 = max(groupLoc$x) + 0.05*xRange
    p = ggplot(pcaProj$sampleCoords, aes(x=PC3, y=PC4, col=group)) +
      geom_point(alpha=0.7) + theme_bw(12) +
      geom_text(mapping=aes(x=x, y=y, label=group, col=group), data=groupLoc, size=3) +
      scale_colour_discrete(guide = FALSE) +
      coord_cartesian(xlim=c(xlim1, xlim2)) +
      xlab(sprintf("PC3 (%.1f%%)", pca.s$importance["Proportion of Variance","PC3"]*100)) +
      ylab(sprintf("PC4 (%.1f%%)", pca.s$importance["Proportion of Variance","PC4"]*100))
    print(p)
  }
}

ensemblToHGNCIDs = function(ensemblIDs, txRefPath) {
  txdb = readRDS(txRefPath)
  genedb = txdb %>% dplyr::select(ensembl_gene_id, external_gene_name) %>% 
    unique() %>% dplyr::rename("gene_id" = "ensembl_gene_id", "gene_name" = "external_gene_name")
  
  gene.df = data.frame(gene_id = ensemblIDs)
  gene.df %<>% dplyr::left_join(genedb %>% dplyr::select(gene_id, gene_name), by="gene_id")
  if (sum(!is.na(gene.df$gene_name)) == 0) {
    warning("ensemblToHGNCIDs: no rownames matched HGNC gene names")
  }
  gene.df
}

HGNCToEnsemblIDs = function(gene_names, txRefPath) {
  txdb = readRDS(txRefPath)
  genedb = txdb %>% dplyr::select(ensembl_gene_id, external_gene_name) %>%
    unique() %>% dplyr::rename("gene_id" = "ensembl_gene_id", "gene_name" = "external_gene_name")
  
  gene.df = data.frame(gene_name = gene_names)
  gene.df %<>% dplyr::left_join(genedb %>% dplyr::select(gene_id, gene_name), by="gene_name")
  if (sum(!is.na(gene.df$gene_id)) == 0) {
    warning("HGNCToEnsemblIDs: no rownames matched Ensembl gene IDs")
  }
  gene.df
}

ensemblToHGNCRownames = function(df, txRefPath) {
  ids.df = ensemblToHGNCIDs(rownames(df), txRefPath)
  # Subset to the genes where we could find an HGNC name, and remove duplicates
  ids.df = ids.df[!is.na(ids.df$gene_name) & !duplicated(ids.df$gene_name),]
  df.new = df[ids.df$gene_id,]
  rownames(df.new) = ids.df$gene_name
  df.new
}

HGNCToEnsemblRownames = function(df, txRefPath) {
  ids.df = HGNCToEnsemblIDs(rownames(df), txRefPath)
  # Subset to the genes where we could find an Ensembl ID, and remove duplicates
  ids.df = ids.df[!is.na(ids.df$gene_id) & !duplicated(ids.df$gene_id),]
  df.new = df[ids.df$gene_name,]
  rownames(df.new) = ids.df$gene_id
  df.new
}

