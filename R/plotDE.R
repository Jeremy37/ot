library(ggplot2)
library(dplyr)

plotDE = function(deseq.res.df, title, sigThreshold=0.01, xlim=NULL, ylim=NULL, maxLabels=100, labelSize=2.2) {
  p = ggplot(deseq.res.df %>% mutate(sig=(padj < sigThreshold)), aes(x=baseMean, y=log2FoldChange, col=sig)) +
    geom_point(size=0.5) +
    scale_color_manual(values=c("black", "red", "dodgerblue")) +
    theme_bw(14) +
    scale_x_log10() + xlab("DESeq2 baseMean expression") +
    theme(legend.position="none") +
    ggtitle(title)
  
  # Determine the points to label by fitting a spline to the points and
  # adjusting it until we get the desired number of points beyond the line
  deseq.res.sig = deseq.res.df %>% filter(padj < sigThreshold) %>%
    arrange(baseMean)
  deseq.res.sig.neg = deseq.res.sig %>% filter(log2FoldChange < 0)
  deseq.res.sig.pos = deseq.res.sig %>% filter(log2FoldChange > 0)
  fit.neg <- smooth.spline(log10(deseq.res.sig.neg$baseMean), deseq.res.sig.neg$log2FoldChange, df=7)
  fit.pos <- smooth.spline(log10(deseq.res.sig.pos$baseMean), deseq.res.sig.pos$log2FoldChange, df=7)
  
  if (nrow(deseq.res.sig) <= maxLabels) {
    deseq.res.sig.plot = deseq.res.sig
  } else {
    # Fit splines to the significant genes, and adjust these to get the
    # desired number of points labelled
    factor = 0.9
    offset = -0.2
    numlabels = sum(deseq.res.sig.pos$log2FoldChange > (predict(fit.pos, log10(deseq.res.sig.pos$baseMean))$y * factor + offset)) +
      sum(deseq.res.sig.neg$log2FoldChange < (predict(fit.neg, log10(deseq.res.sig.neg$baseMean))$y * factor - offset))
    while (numlabels > maxLabels) {
      offset = offset + 0.1
      factor = factor * 1.04
      numlabels = sum(deseq.res.sig.pos$log2FoldChange > (predict(fit.pos, log10(deseq.res.sig.pos$baseMean))$y * factor + offset)) +
        sum(deseq.res.sig.neg$log2FoldChange < (predict(fit.neg, log10(deseq.res.sig.neg$baseMean))$y * factor - offset))
    }
    deseq.res.sig.plot = rbind(deseq.res.sig.pos[deseq.res.sig.pos$log2FoldChange > (predict(fit.pos, log10(deseq.res.sig.pos$baseMean))$y * factor + offset),],
                               deseq.res.sig.neg[deseq.res.sig.neg$log2FoldChange < (predict(fit.neg, log10(deseq.res.sig.neg$baseMean))$y * factor - offset),])
  }
  p = p + annotate(geom="text", x=deseq.res.sig.plot$baseMean, y=deseq.res.sig.plot$log2FoldChange,
                   label=deseq.res.sig.plot$gene_name, col="blue", hjust = -0.2, size=labelSize) 
  #p = p + geom_line(aes(10^x, y*factor - offset), data=as.data.frame(fit.neg[c("x","y")]), col="grey90") +
  #  geom_line(aes(10^x, y*factor + offset), data=as.data.frame(fit.pos[c("x","y")]), col="grey90")
  if (!is.null(xlim) | !is.null(ylim)) {
    p = p + coord_cartesian(xlim=xlim, ylim=ylim)
  }
  return(p)
}
