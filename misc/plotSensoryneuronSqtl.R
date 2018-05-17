#!/usr/bin/env Rscript
library(tidyverse)
library(gridExtra)
library(grid)
myargs <- NULL

# myargs = list()
# myargs$cluster = "clu_55409"
# myargs$snpid = "rs4147914"
# myargs$genotypeCounts = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/all.leadSNPs.genotype_counts.txt.gz"
# myargs$genotypeMap = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/sensoryneuron.merged_sample_genotype_map.txt"
# myargs$clusterPerind = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/sqtl/sqtl.clusters_perind_numers.ratios.txt.gz"
# myargs$nominals = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/sqtl/sqtl.fastqtl.nominals.txt.gz"

main = function()
{
  myargs <<- getArgs()
  print(myargs)
  
  if (is.null(myargs$cluster)) {
    stop("Missing parameter 'cluster'")
  }
  if (is.null(myargs$snpid)) {
    stop("Missing parameter 'snpid'")
  }
  if (is.null(myargs$genotypeCounts)) {
    stop("Missing parameter 'genotypecounts'")
  }
  if (is.null(myargs$genotypeMap)) {
    stop("Missing parameter 'genotypeMap'")
  }
  if (is.null(myargs$clusterPerind)) {
    stop("Missing parameter 'clusterPerind'")
  }
  if (is.null(myargs$output)) {
    myargs$output = paste0("sqtl_plot.", myargs$cluster, ".", myargs$snpid)
  }
  
  gtc.df = read.delim(myargs$genotypeCounts)
  colnames(gtc.df)[1] = "chr"
  gtc.df = gtc.df %>% dplyr::filter(ID == myargs$snpid)
  if (nrow(gtc.df) == 0) {
    stop(sprintf("SNP %s is not in table of genotype counts.", myargs$snpid))
  } else if (nrow(gtc.df) > 1) {
    stop(sprintf("Something is wrong - SNP %s appears multiple times in the table of genotype counts.", myargs$snpid))
  }
  gtc.snp.df = gtc.df %>% dplyr::select(-c(chr, POS, REF, ALT, ID)) %>%
    tidyr::gather(key = sample, value = alt_allele_count)
  gtc.snp.df$sample = gsub(".", "-", gtc.snp.df$sample, fixed=T)
  
  genotypemap.df = read.delim(myargs$genotypeMap)
  rownames(genotypemap.df) = genotypemap.df$sample
  
  cluster.df = readr::read_tsv(myargs$clusterPerind)
  sampleNames = colnames(cluster.df)[2:ncol(cluster.df)]
  cluster.df$clusterName = sapply(cluster.df$cluster, function(s) { strsplit(s, ":", fixed=T)[[1]][[4]]})
  cluster.df$chr = sapply(cluster.df$cluster, function(s) { strsplit(s, ":", fixed=T)[[1]][[1]]})
  cluster.df$start = sapply(cluster.df$cluster, function(s) { strsplit(s, ":", fixed=T)[[1]][[2]]})
  cluster.df$end = sapply(cluster.df$cluster, function(s) { strsplit(s, ":", fixed=T)[[1]][[3]]})
  # Select the relevant cluster from the clusters file
  cluster.df = cluster.df %>% dplyr::filter(clusterName == myargs$cluster) %>%
    dplyr::select(-cluster, everything()) %>%
    dplyr::rename(junction = cluster)
  if (nrow(cluster.df) == 0) {
    stop(sprintf("Cluster %s not found in clusters input file", myargs$cluster))
  }
  
  # Map the sample names in the clusters file to the genotype IDs
  if (!all(sampleNames %in% genotypemap.df$sample)) {
    stop("Not all sample names from the clusters file are present in the genotype map.")
  }
  colnames(cluster.df)[1:length(sampleNames)] = as.character(genotypemap.df[sampleNames,]$genotypeid)
  
  # Reshape the clusters for plotting
  cluster.plot.df = cluster.df %>% tidyr::gather(key = sample, value = ratio, -c(clusterName:junction))
  cluster.plot.df = cluster.plot.df %>% dplyr::left_join(gtc.snp.df, by="sample")
  
  # Make a separate plot for each 
  pdf(file=paste0(myargs$output, ".pdf"))
  
  p = ggplot(cluster.plot.df, aes(x=alt_allele_count, y=ratio, group=alt_allele_count)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.2, alpha=0.7) +
    facet_wrap(~junction) +
    theme_bw() + scale_x_continuous(breaks = c(0,1,2)) +
    ggtitle(paste0(myargs$cluster, "\nFraction of reads by intron junction and genotype"))
  print(p)
  
  if (!is.null(myargs$nominals)) {
    snp.nominals.str <- system(sprintf("gzip -cd %s | grep %s", myargs$nominals, myargs$snpid), intern=TRUE)
    if (length(snp.nominals.str) == 0) {
      write(sprintf("Did not find any p values in nominals file for SNP %s", myargs$snpid), stderr())
    } else {
      snp.nominals.con <- textConnection(snp.nominals.str)
      snp.nominals.df <- read.delim(snp.nominals.con, header=FALSE, stringsAsFactors = F, col.names=c("cluster", "snp", "dist_to_intron_centre", "pval", "slope", "snp_chr", "snp_pos")) %>%
        dplyr::select(cluster, snp, pval)
      snp.nominals.df$clusterName = sapply(snp.nominals.df$cluster, function(s) { strsplit(s, ":", fixed=T)[[1]][[4]]})
      snp.nominals.df = snp.nominals.df %>% dplyr::filter(clusterName == myargs$cluster)
      plot.new()
      grid.table(snp.nominals.df, theme = gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.8))))
    }
  }
  
  dev.off()
}


###########################################################################
##' commandArgs parsing
##' 
##' return a named list of command line arguments
##'
##' Usage:
##' call the R script thus
##'   ./myfile.R --args myarg=something
##' or
##'   R CMD BATCH --args myarg=something myfile.R
##'
##' Then in R do
##'   myargs <- getArgs()
##' and myargs will be a named list
##' > str(myargs)
##' List of 2
##' $ file : chr "myfile.R"
##' $ myarg: chr "something"
##'
##' @title getArgs
##' @param verbose print verbage to screen 
##' @param defaults a named list of defaults, optional
##' @return a named list
##' @author Chris Wallace
getArgs = function(verbose=FALSE, defaults=NULL) {
  myargs <- gsub("^--","",commandArgs(TRUE))
  setopts <- !grepl("=",myargs)
  if(any(setopts))
    myargs[setopts] <- paste(myargs[setopts],"=notset",sep="")
  myargs.list <- strsplit(myargs,"=")
  myargs <- lapply(myargs.list,"[[",2 )
  names(myargs) <- lapply(myargs.list, "[[", 1)
  
  ## logicals
  if(any(setopts))
    myargs[setopts] <- TRUE
  
  ## defaults
  if(!is.null(defaults)) {
    defs.needed <- setdiff(names(defaults), names(myargs))
    if(length(defs.needed)) {
      myargs[ defs.needed ] <- defaults[ defs.needed ]
    }
  }
  
  ## verbage
  if(verbose) {
    cat("read",length(myargs),"named args:\n")
    print(myargs)
  }
  myargs
}


###########################################################################

main()
