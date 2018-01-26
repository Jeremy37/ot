library(DESeq2)
library(readr)

readtsv = function(fname, ...) { readr::read_tsv(fname, ...) %>% as.data.frame() }

readCounts = function(fname) {
  counts = readtsv(fname)
  # Fix gene names to remove the version number (e.g. ".4" in ENSG00000223972.4)
  gene_ids <- gsub("\\.[\\d]+", "", counts[,1], perl=T)
  counts[,1] = gene_ids
  rownames(counts) = gene_ids
  counts
}

# Takes the gene read counts output from featureCounts (which has the first two
# cols as gene_id and gene_length), and gets stabilized counts using DESeq2,
# and then returns a matrix of RPKMs, with geneIDs as the rownames.
countsToRPKM = function(counts, gene.meta=NULL, useDeseq=T) {
  if (is.null(gene.meta)) {
    gene.meta = counts[,c(1,2)]
    counts.mat = counts[, -c(1,2)]
    rownames(gene.meta) = gene.meta$gene_id
    rownames(counts.mat) = gene.meta$gene_id
  } else {
    counts.mat = counts
  }
  
  if (!all(rownames(counts.mat) == gene.meta$gene_id)) {
    stop("countsToRPKM: Error - rownames of counts.mat matrix should be the same as the gene IDs in the gene metadata.")
  }

  if (useDeseq) {
    coldata = data.frame(sampleID = colnames(counts.mat))
    dds = DESeqDataSetFromMatrix(counts.mat, coldata, ~0)
    colnames(dds) = colnames(counts.mat)
    vst = varianceStabilizingTransformation(dds)
    dds = estimateSizeFactors(dds)
    #View(assay(vst))
    counts.norm = counts(dds, normalized=T)
    colnames(counts.norm) = colnames(counts.mat)
  } else {
    counts.norm = counts.mat
  }
  
  rpm = apply(counts.norm, MARGIN=2, FUN=function(x) x / sum(x)) * 1e6
  rpkm = rpm * 1e3
  rpkm = apply(rpkm, MARGIN=2, FUN=function(x) x / gene.meta$length)
  rpkm
}
