#!/usr/bin/env Rscript

# This script reads in a VCF file and an "ASECounts" file. The VCF file has the genotypes
# for a single individual in the region of 1 gene. The ASECounts file has allele-specific
# read counts for possibly multiple samples that come from the individual and which cover
# the gene. The script performs statistical tests for allele-specific expression, and
# produces some plots.
library(tidyverse)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(wiggleplotr)
library(gridExtra)
options(stringsAsFactors = F)
myargs <- NULL

#root = "/Users/jeremys/work/opentargets"
root = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys"
txdb_file = file.path(root, "reference/GRCh38/ensembldb.91.rds")
txmeta_file = file.path(root, "reference/GRCh38/ensembldb.91.transcripts.rds")

# myargs <<- list(asefile = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/ipsneurons/GRCh38/RNA/analysis/ipsneurons.kolf2.ASEcounts.nochr.gz",
#                 vcf = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/ipsneurons/GRCh38/genotypes/kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.biallelic.snps.hets.vcf.gz",
#                 allgeneAse = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/ipsneurons/GRCh38/RNA/analysis/ase/ase.ipsneurons.expressed_genes.ase.txt",
#                 sampleSet = "ipsneurons",
#                 gene = "SNCA",
#                 sampleGenotypeMap = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/ipsneurons/GRCh38/RNA/analysis/ipsneurons.sample_genotype_map.txt",
#                 sampleGroupsFile = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/ipsneurons/GRCh38/RNA/analysis/ipsneurons.sampleGroups.txt",
#                 out = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/ipsneurons/GRCh38/RNA/analysis/ase/ase")

# myargs <<- list(asefile = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/RNA/ase/sensoryneuron.ASEcounts.9samples.gz",
#                 vcf = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/imputed.9_samples.snps.biallelic.hets.INFO_08.MAF_0.05.RASQUAL.vcf.gz",
#                 allgeneAse = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/RNA/ase/sensoryneuron.expressed_genes.9.ase.txt",
#                 sampleSet = "sensoryneurons",
#                 gene = "PTK2B",
#                 sampleGenotypeMap = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/RNA/sensoryneuron.sample_genotype_map.txt",
#                 sampleGroupsFile = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/RNA/sensoryneuron.sample_groups.txt",
#                 codingOnly = "F",
#                 out = "/Users/jeremys/work/opentargets/sensoryneurons/GRCh38/RNA/ase/ase.9samples")

main = function()
{
  myargs <<- getArgs()
  # myargs <<- list(asefile = "/Users/jeremys/work/opentargets/ipsneurons/GRCh38/RNA/analysis/ipsneurons.kolf2.ASEcounts.nochr.gz",
  #                 vcf = "/Users/jeremys/work/opentargets/ipsneurons/GRCh38/genotypes/kolf_2.imputed_phased.20150604.GRCh38.INFO.0.8.biallelic.snps.hets.vcf.gz",
  #                 allgeneAse = "/Users/jeremys/work/opentargets/ipsneurons/GRCh38/RNA/analysis/ase/ase.ipsneurons.expressed_genes.ase.new.txt",
  #                 sampleSet = "ipsneurons",
  #                 gene = "SNCA",
  #                 sampleGenotypeMap = "/Users/jeremys/work/opentargets/ipsneurons/GRCh38/RNA/analysis/ipsneurons.sample_genotype_map.txt",
  #                 sampleGroupsFile = "/Users/jeremys/work/opentargets/ipsneurons/GRCh38/RNA/analysis/ipsneurons.sampleGroups.txt",
  #                 codingOnly = "T",
  #                 out = "/Users/jeremys/work/opentargets/ipsneurons/GRCh38/RNA/analysis/ase/ase")
  print(myargs)
  
  if (is.null(myargs$asefile)) {
    stop("Missing parameter 'asefile'")
  }
  if (is.null(myargs$vcf)) {
    stop("Missing parameter 'vcf'")
  }
  if (is.null(myargs$allgeneAse)) {
    # This parameter indicates the file where summary data for all genes will be stored.
    # This needs to be computed once for each sample set, and in subsequent runs looking
    # at ASE for a specific gene, will be loaded from an RDS file.
    stop("Missing parameter 'allgeneAse'")
  }
  if (is.null(myargs$gene)) {
    stop("Missing parameter 'gene'")
  }
  if (is.null(myargs$sampleGenotypeMap)) {
    stop("Missing parameter 'sampleGenotypeMap'")
  }
  sampleGenotypeMap <- read.delim(myargs$sampleGenotypeMap)
  #, header = F, col.names = c("Sample", "genotypeID", "displayName")
  
  if (is.null(myargs$minhetprob)) {
    # Minimum genotype probability of heterozygote
    myargs$minhetprob <<- 0.9
  }
  myargs$minhetprob <<- as.numeric(myargs$minhetprob)
  
  if (!is.null(myargs$sampleGroupsFile)) {
    sampleGroups <- read.delim(myargs$sampleGroupsFile)
    sampleGroups$displayName = factor(as.character(sampleGroups$displayName), levels=sampleGroups$displayName)
  }
  if (!is.null(myargs$codingOnly)) {
    codingOnly = as.logical(myargs$codingOnly)
    if (!is.na(codingOnly) & codingOnly == F) {
      myargs$codingOnly <<- FALSE
    }
  } else {
    myargs$codingOnly <<- TRUE
  }
  if (is.null(myargs$out)) {
    myargs$out <<- paste0(myargs$gene, ".", "ase")
  }
  if (!is.null(myargs$grch37)) {
    if (myargs$grch37 != "F") {
      txdb_file <<- file.path(root, "reference/GRCh37/ensembldb.91.grch37.rds")
      txmeta_file <<- file.path(root, "reference/GRCh37/ensembldb.91.grch37.transcripts.rds")
    }
  }
  
  vcf.all.df = readr::read_tsv(myargs$vcf, comment = "##") %>%
    dplyr::rename(chr = `#CHROM`, pos=POS, id=ID, ref=REF, alt=ALT) %>%
    dplyr::filter(id != ".") %>%
    dplyr::filter(!duplicated(id))
  
  # Read in the table of allele-specific counts
  ase.all.df = readr::read_tsv(myargs$asefile, col_types="cciccciiiiiiii") %>%
    dplyr::rename(chr=contig, pos=position) %>%
    dplyr::mutate(variantLabel = paste(variantID, refAllele, altAllele, sep="_")) %>%
    dplyr::filter(variantID %in% vcf.all.df$id) %>%
    dplyr::left_join(sampleGenotypeMap, by="Sample")

  if (!is.null(myargs$sampleGroupsFile)) {
    ase.all.df$Sample = factor(as.character(ase.all.df$Sample), levels=sampleGroups$Sample)
  }
  
  # Subset to the gene of interest
  # First get gene exon coordinates
  print(paste0("Loading transcripts file: ", txdb_file))
  txdb = loadDb(txdb_file)
  exons = exonsBy(txdb, by = "tx", use.names = TRUE)
  cdss = cdsBy(txdb, by = "tx", use.names = TRUE)
  #rtracklayer::export.bed(reduce(unlist(exons)), paste0(myargs$out, ".ensembldb.91.exons.bed"))
  
  txmeta = readRDS(txmeta_file)
  selected_transcripts = txmeta %>% dplyr::filter(gene_name == myargs$gene)
  if (myargs$codingOnly) {
    selected_transcripts = selected_transcripts %>% dplyr::filter(transcript_biotype == "protein_coding")
  }
  tx_ids = selected_transcripts$transcript_id
  gene_id = selected_transcripts$gene_id[1]
  gene.exons = exons[tx_ids]

  # Harmonize chromosome names
  if (grepl("chr", ase.all.df$chr[1])) {
    removeChr = function(x) { gsub("^chr", "", x, perl=T)}
    ase.all.df$chr = sapply(ase.all.df$chr, removeChr)
  }
  if (grepl("chr", vcf.all.df$chr[1])) {
    removeChr = function(x) { gsub("^chr", "", x, perl=T)}
    vcf.all.df$chr = sapply(vcf.all.df$chr, removeChr)
  }
  
  vcf.gene.df = subsetByChrPos(vcf.all.df, gene.exons)
  ase.gene.df = subsetByChrPos(ase.all.df, gene.exons)
  print(ase.gene.df)
  if (nrow(ase.gene.df) < 1) {
    stop(paste0("No heterozygous SNPs within gene ", myargs$gene, " to use for ASE."))
  }
  genotype.df = getGenotypeDF(vcf.gene.df)

  # Add in genotypes
  ase.gene.df = ase.gene.df %>%
    dplyr::left_join(genotype.df %>% dplyr::select(genotypeID, id, gt), by=c("genotypeID", "variantID" = "id")) 
  
  variantLabels = vcf.gene.df %>% dplyr::select(pos, id, ref, alt) %>%
    arrange(pos) %>%
    dplyr::mutate(variantLabel = paste(id, ref, alt, sep="_"))

  if (!is.null(myargs$sampleGroupsFile)) {
    # Replace sample names with display names from the sample groups file
    rownames(sampleGroups) = sampleGroups$Sample
    ase.gene.df$displayName = factor(as.character(sampleGroups[ase.gene.df$Sample,]$displayName), levels=sampleGroups$displayName)
  }
  
  # Determine the haplotype-specific read counts
  ase.gene.df$haplotype1 = NA
  ase.gene.df$haplotype1[!is.na(ase.gene.df$gt) & ase.gene.df$gt == "0|1"] = ase.gene.df$refCount[!is.na(ase.gene.df$gt) & ase.gene.df$gt == "0|1"]
  ase.gene.df$haplotype1[!is.na(ase.gene.df$gt) & ase.gene.df$gt == "1|0"] = ase.gene.df$altCount[!is.na(ase.gene.df$gt) & ase.gene.df$gt == "1|0"]
  ase.gene.df$haplotype2 = NA
  ase.gene.df$haplotype2[!is.na(ase.gene.df$gt) & ase.gene.df$gt == "0|1"] = ase.gene.df$altCount[!is.na(ase.gene.df$gt) & ase.gene.df$gt == "0|1"]
  ase.gene.df$haplotype2[!is.na(ase.gene.df$gt) & ase.gene.df$gt == "1|0"] = ase.gene.df$refCount[!is.na(ase.gene.df$gt) & ase.gene.df$gt == "1|0"]
  
  ase.gene.df$variantLabel = factor(as.character(ase.gene.df$variantLabel), levels = as.character(variantLabels$variantLabel))
  ase.gene.df = ase.gene.df %>% dplyr::filter(!is.na(haplotype1) & !is.na(haplotype2))

  # Before we generate plots for ASE, we need the distribution of ASE across
  # all genes and samples, which we precompute once. It takes a while to get
  # genotypes for all SNPs, so we only do this for all SNPs if we haven't
  # already precomputed gene ASE across all genes
  if (!file.exists(myargs$allgeneAse)) {
    # Precompute the distribution of ASE across all genes and samples
    print("Computing ASE for all genes")
    genestats.df = getASEAllGenes(ase.all.df, vcf.all.df)
  } else {
    print(paste("Loading ASE for all genes from file:", myargs$allgeneAse))
    genestats.df = readr::read_tsv(myargs$allgeneAse)
  }
  
  #exons.df = readr::read_tsv("/Users/jeremys/work/opentargets/reference/GRCh38/Homo_sapiens.GRCh38.91.chr.exon_start_end.bed",
  #                           col_names = c("chr", "start", "end", "geneid")) %>%
  #  dplyr::filter(geneid == myargs$gene)
  # exons.gr = reduce(GRanges(seqnames = exons.df$chr,
  #                           ranges = IRanges(exons.df$start, exons.df$end),
  #                           strand = NA,
  #                           exons.df[,4]))
  
  # In determining how to plot the alleles, we want to plot counts for all
  # haplotype 1 alleles (which each may be ref or alt) and counts for
  # haplotype 2 alleles.
  ase.plot.df = ase.gene.df %>%
    dplyr::select(displayName, variantLabel, haplotype1, haplotype2) %>%
    tidyr::gather(allele, count, haplotype1:haplotype2)
  ase.plot.df$allele = factor(ase.plot.df$allele, levels = c("haplotype1", "haplotype2"))
  
  pdf(file = paste0(myargs$out, ".", myargs$gene, ".pdf"), width = 8, height = 7)
  
  p1 = ggplot(ase.plot.df, aes(x=variantLabel, y=count)) +
    geom_bar(aes(fill = allele), position = "dodge", stat="identity") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~displayName) +
    ggtitle(paste0("ASE per sample - ", myargs$gene)) + xlab("Variant") +
    ylim(0, max(5, max(ase.plot.df$count)))
  print(p1)
  
  ase.byvariant = ase.plot.df %>% group_by(variantLabel, allele) %>%
    summarise(count = sum(count))
  p2 = ggplot(ase.byvariant, aes(x=variantLabel, y=count)) +
    geom_bar(aes(fill = allele), position = "dodge", stat="identity") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste0("ASE per variant (across all samples) - ", myargs$gene)) + xlab("Variant") +
    ylim(0, max(5, max(ase.byvariant$count)))
  print(p2)

  # Display the positions of the SNPs used for ASE on top of the
  # transcript definitions
  gene.variants = ase.gene.df %>%
    dplyr::select(chr, pos, variantID, ref=refAllele, alt=altAllele, totalCount) %>%
    dplyr::group_by(variantID) %>%
    dplyr::summarise(chr=first(chr), pos=first(pos), ref=first(ref), alt=first(alt), totalCount=sum(totalCount)) %>%
    dplyr::arrange(pos)
  gene.exon.ranges = ranges(unlist(gene.exons))
  minx = min(gene.exon.ranges@start)
  maxx = max(gene.exon.ranges@start + gene.exon.ranges@width)
  minx = minx - (maxx - minx) * 0.05
  maxx = maxx + (maxx - minx) * 0.05
  cds_ids = tx_ids[tx_ids %in% names(cdss)]
  # ptr = plotTranscripts(exons[tx_ids], cdss[cds_ids], txmeta, rescale_introns = F) +
  #   geom_segment(aes(x=pos, xend=pos, y=0.5, yend=numtranscripts), data=gene.variants, colour="blue", size=0.3, alpha=0.5) +
  #   geom_text(aes(x=pos, y=0.5, label=variantID, angle=45, vjust="top", hjust="right"), data=gene.variants, size=3) +
  #   coord_cartesian(ylim=c(0.5,numtranscripts+0.3), xlim=c(minx, maxx)) +
  #   ggtitle("Transcripts, and variants used for ASE")
  gene.variants$numtranscripts = length(tx_ids)
  ptr = plotTranscripts(exons[tx_ids], cdss[cds_ids], txmeta, rescale_introns = F) +
    geom_segment(aes(x=pos, xend=pos, y=0.5, yend=numtranscripts), data=gene.variants, colour="blue", size=0.3, alpha=0.5) +
    geom_text(aes(x=pos, y=0.5, label=variantID, angle=45, vjust="top", hjust="right"), data=gene.variants, size=3) +
    coord_cartesian(ylim=c(0.5,length(tx_ids)+0.3), xlim=c(minx, maxx)) +
    ggtitle("Transcripts, and variants used for ASE")
  plot(ptr)
  gene.variants = gene.variants %>% dplyr::select(-numtranscripts)
  
  # Display details of the SNPs used for ASE as a table
  plot.new()
  tt = ttheme_default(core = list(fg_params=list(fontsize=9)),
                      colhead = list(fg_params=list(fontsize=9)),
                      rowhead = list(fg_params=list(fontsize=9)))
  grid.table(gene.variants, theme=tt)
  
  ase.bysample = ase.plot.df %>% group_by(displayName, allele) %>%
    summarise(count = sum(count, na.rm=T))
  p3 = ggplot(ase.bysample, aes(x=displayName, y=count)) +
    geom_bar(aes(fill = allele), width=0.8, position = "dodge", stat="identity") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste0("ASE per sample - ", myargs$gene)) + xlab("Sample")
  print(p3)
  
  if (!is.null(myargs$sampleGroupsFile)) {
    ase.bysample = ase.bysample %>% dplyr::left_join(sampleGroups, by="displayName")
    p4 = ggplot(ase.bysample, aes(x=displayName, y=count)) +
      geom_bar(aes(fill = Group, shape=allele), width=0.8, position = "dodge", stat="identity") +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(paste0("ASE per sample - ", myargs$gene)) + xlab("Sample")
    print(p4)
  }
  
  # Statistically test for ASE in individual samples
  ase.bysample.spread = ase.bysample %>% ungroup() %>% tidyr::spread(allele, count)
  ase.bysample.spread$binomial_p = NA
  #ase.bysample.spread$betabinomial_p = NA
  for (i in 1:nrow(ase.bysample.spread)) {
    hap1count = ase.bysample.spread[i,]$haplotype1
    hap2count = ase.bysample.spread[i,]$haplotype2
    if (hap1count + hap2count > 0) {
      ase.bysample.spread[i,]$binomial_p = pbinom(min(hap1count, hap2count), hap1count + hap2count, 0.5) * 2
    }
  }
  # This can happen only when you have exactly equally balanced counts,
  # e.g. 3 of 6 for each allele, in which case p(allele count <= 3) > 0.5
  ase.bysample.spread$binomial_p[ase.bysample.spread$binomial_p > 1] = 1
  
  # Plot for each sample the quantile of gene ASE this gene falls at
  # Note that some genes will not be in the genestats.df, since they
  # had too few ASE SNPs of sufficient depth
  genestats.gene.df = genestats.df %>% dplyr::filter(geneid == gene_id)
  
  ase.bysample.spread = ase.bysample.spread %>%
    dplyr::left_join(genestats.gene.df)
  
  # Display a summary of the stats as a plot
  plot.new()
  tt = ttheme_default(core = list(fg_params=list(fontsize=9)),
                      colhead = list(fg_params=list(fontsize=9)),
                      rowhead = list(fg_params=list(fontsize=9)))
  grid.table(ase.bysample.spread %>%
               dplyr::select(Sample=displayName, Group, haplotype1, haplotype2, binomial_p, ase, depthBin, binAseQuantile),
             theme=tt)
  if (nrow(genestats.gene.df) > 0) {
    genestats.display.df = genestats.df %>% dplyr::left_join(sampleGroups)
    #genestats.display.df = genestats.display.df %>%
    #  dplyr::inner_join(ase.bysample.spread %>% dplyr::select(Sample, depthBin), by=c("Sample", "depthBin"))
    p5 = ggplot(genestats.display.df %>% dplyr::filter(Sample %in% unique(ase.bysample.spread$Sample)),
           aes(x=ase)) +
      geom_histogram(binwidth=0.01) +
      geom_vline(aes(xintercept = ase, color="red"), data=ase.bysample.spread) +
      facet_wrap(~displayName) +
      theme_bw() + ggtitle("Distribution of gene ASE per sample (rel. to all genes)") +
      scale_colour_discrete(guide=F)
    print(p5)
  }
  
  dev.off()
}

getASEAllGenes = function(ase.df, vcf.df) {
  genotype.df = getGenotypeDF(vcf.df)
  # Add in genotypes
  ase.df = ase.df %>%
    dplyr::left_join(genotype.df %>% dplyr::select(genotypeID, id, gt), by=c("genotypeID", "variantID" = "id")) 
  # Determine the haplotype-specific read counts
  ase.df$haplotype1 = NA
  ase.df$haplotype1[!is.na(ase.df$gt) & ase.df$gt == "0|1"] = ase.df$refCount[!is.na(ase.df$gt) & ase.df$gt == "0|1"]
  ase.df$haplotype1[!is.na(ase.df$gt) & ase.df$gt == "1|0"] = ase.df$altCount[!is.na(ase.df$gt) & ase.df$gt == "1|0"]
  ase.df$haplotype2 = NA
  ase.df$haplotype2[!is.na(ase.df$gt) & ase.df$gt == "0|1"] = ase.df$altCount[!is.na(ase.df$gt) & ase.df$gt == "0|1"]
  ase.df$haplotype2[!is.na(ase.df$gt) & ase.df$gt == "1|0"] = ase.df$refCount[!is.na(ase.df$gt) & ase.df$gt == "1|0"]
  
  ase.gr = GRanges(seqnames = ase.df$chr, ranges = IRanges(ase.df$pos, ase.df$pos), strand = NA,
                   ase.df[, "Sample"], ase.df[, 4:ncol(ase.df)])

  txdb = loadDb(txdb_file)
  exons = exonsBy(txdb, by = "tx", use.names = TRUE)
  txmeta = readRDS(txmeta_file)
  if (myargs$codingOnly) {
    txmeta = txmeta %>% dplyr::filter(transcript_biotype == "protein_coding")
  }
  tx_ids = txmeta$transcript_id
  exons.gr = exons[tx_ids]
  #exons.df = readr::read_tsv(file.path(root, "reference/GRCh38/Homo_sapiens.GRCh38.91.exon_start_end.bed"),
  #                           col_types = "ciic", col_names = c("chr", "start", "end", "geneid"))
  #exons.gr = GRanges(seqnames = exons.df$chr, ranges = IRanges(exons.df$start, exons.df$end), strand = NA,
  #                   exons.df[,4])
  
  ase.genehits = as.data.frame(findOverlaps(ase.gr, exons.gr)) %>% dplyr::filter(!duplicated(queryHits))
  ase.overlaps = ase.df[ase.genehits$queryHits,]
  ase.overlaps$geneid = txmeta[ase.genehits$subjectHits,]$gene_id
  ase.overlaps = ase.overlaps %>% dplyr::filter(!is.na(haplotype1) & !is.na(haplotype2))
  
  gene.counts = ase.overlaps %>% group_by(geneid) %>%
    summarise(totalCounts = sum(totalCount),
              hap1counts = sum(haplotype1),
              hap2counts = sum(haplotype2),
              ase = min(hap1counts / totalCounts, hap2counts / totalCounts))
  
  numSamples = length(unique(ase.overlaps$Sample))
  # Only include genes where we have at least 10 total reads over het SNPs
  # per sample on average
  countPerSampleThreshold = 10
  selected.gene.counts = gene.counts %>% dplyr::filter(totalCounts >= numSamples * countPerSampleThreshold)
  selectedGenes = selected.gene.counts$geneid
  ase.overlaps.sel = ase.overlaps %>% dplyr::filter(geneid %in% selectedGenes)
  
  pdf(file = paste0(myargs$out, ".", myargs$sampleSet, ".expressed_genes.pdf"), width = 8, height = 7)
  
  p1 = ggplot(gene.counts, aes(x=log2(totalCounts))) + geom_histogram(bins=50) +
    theme_bw() + ggtitle("Total coverage of het SNPs per gene")
  print(p1)
  
  p2 = ggplot(gene.counts %>% dplyr::filter(geneid %in% selectedGenes), aes(x=log2(totalCounts))) + geom_histogram(bins=50) +
    theme_bw() + ggtitle("Total coverage of het SNPs per gene - selected genes")
  print(p2)
  
  addBinomialStats = function(selected.gene.counts, sampleid) {
    getBinomialP = function(i) {
      hap1 = selected.gene.counts[i,]$hap1counts
      hap2 = selected.gene.counts[i,]$hap2counts
      # Use the smaller haplotype count to get a binomial p value
      # and multiply by 2, since we want to know how likely it is
      # to get an imbalance more than this in either direction
      pbinom(min(hap1, hap2), hap1 + hap2, 0.5) * 2
    }
    selected.gene.counts$binomialP = sapply(1:nrow(selected.gene.counts), getBinomialP)
    selected.gene.counts$aseQuantile = rank(selected.gene.counts$ase) / nrow(selected.gene.counts)
    selected.gene.counts$sigQuantile = rank(selected.gene.counts$binomialP) / nrow(selected.gene.counts)
    
    # Determine quantiles of the number of reads for each gene
    numQuantiles = 10
    countQuantiles = as.integer(quantile(selected.gene.counts$totalCounts, probs=seq(0, 1, 1/numQuantiles)))
    quantileLabels = sapply(1:numQuantiles, function(i) { sprintf("%d-%d", countQuantiles[i], countQuantiles[i+1])})
    selected.gene.counts$depthBin = NA
    selected.gene.counts$binAseQuantile = NA
    for (i in 1:numQuantiles) {
      binSelect = (selected.gene.counts$totalCounts >= countQuantiles[i] & selected.gene.counts$totalCounts <= countQuantiles[i+1])
      selected.gene.counts[binSelect,]$depthBin = quantileLabels[i]
      selected.gene.counts[binSelect,]$binAseQuantile = rank(selected.gene.counts[binSelect,]$ase) / sum(binSelect)
    }
    selected.gene.counts
  }
  
  selected.gene.counts = addBinomialStats(selected.gene.counts)
  p3 = ggplot(selected.gene.counts %>% dplyr::filter(sigQuantile > 0.05), aes(log10(binomialP))) +
    geom_histogram(bins=50) + theme_bw() +
    ggtitle("Distribution of binomial P values for gene ASE\n(where ASE counts aggregated across all samples)")
  print(p3)
  
  dev.off()
  
  genestats.list = list(allsamples = data.frame(Sample = "all samples", selected.gene.counts))
  
  ###################################################################
  # Now determine ASE quantiles per sample, rather than globally
  sampleids = unique(ase.overlaps$Sample)
  for (sampleid in sampleids) {
    selected.gene.counts = ase.overlaps %>%
      dplyr::filter(Sample == sampleid) %>%
      group_by(geneid) %>%
      summarise(totalCounts = sum(totalCount),
                hap1counts = sum(haplotype1),
                hap2counts = sum(haplotype2),
                ase = min(hap1counts / totalCounts, hap2counts / totalCounts)) %>%
      dplyr::filter(totalCounts >= countPerSampleThreshold)
    selected.gene.counts = addBinomialStats(selected.gene.counts)
    genestats.list = c(genestats.list, list(data.frame(Sample = sampleid, selected.gene.counts)))
  }
  
  genestats.df = bind_rows(genestats.list)
  write.table(genestats.df, file=myargs$allgeneAse, sep="\t", quote=F, row.names=F, col.names=T)
  
  #ggplot(genestats.df %>% dplyr::filter(sigQuantile > 0.05), aes(log10(binomialP), fill=Sample)) +
  #  geom_density(alpha=0.3) + theme_bw()
  genestats.df
}

getGenotypeDF = function(vcf.df) {
  # Extract sample genotype. This isn't really correct, because here we assume
  # that the GT field is at the same position for each variant, which it may not be
  #format = "GT:GP:GC:IA:IB:BAF:LRR:PS"
  #sample = "1|0:0,1,0:0.9113:0.446:0.408:0.5689:0.4171:15981712"
  #formatFields = strsplit(format, split=":")[[1]]
  formatFields = strsplit(vcf.df$FORMAT[1], split=":")[[1]]
  gtpos = which(formatFields == "GT")[1]
  gppos = which(formatFields == "GP")[1]
  getGenotype = function(samp) {
    sampFields = strsplit(samp, split=":")[[1]]
    return(sampFields[gtpos])
  }
  getHetProb = function(samp) {
    sampFields = strsplit(samp, split=":")[[1]]
    gtProbs = strsplit(sampFields[gppos], split=",")[[1]]
    return(as.numeric(gtProbs[2]))
  }
  
  # Create a table with columns:
  # genotypeID  rsID  genotype
  # This has samples underneath each other, rather than a column per sample,
  # to facilitate joining on this DF later
  genotype.df = data.frame(genotypeID=character(), id=character(), genotype=character())
  genotype.df = vcf.df %>% dplyr::select(id, 10:ncol(vcf.df))
  #vcf.df.test = vcf.df
  #vcf.df.test[,11] = vcf.df.test[,10]
  #genotype.df = vcf.df.test %>% dplyr::select(id, 10:ncol(vcf.df.test))
  genotype.df = genotype.df %>%
    tidyr::gather(key=genotypeID, value=info, 2:ncol(genotype.df))
  genotype.df$gt = sapply(genotype.df$info, getGenotype)
  genotype.df$hetprob = sapply(genotype.df$info, getHetProb)
  # Subset to heterozygous genotypes
  genotype.df = genotype.df %>% dplyr::filter(hetprob > myargs$minhetprob)
  genotype.df
}

subsetByChrPos = function(df, gene.exons) {
  df.gr = GRanges(seqnames = df$chr, ranges = IRanges(df$pos, df$pos), strand = NA)
  hits = as.data.frame(findOverlaps(df.gr, gene.exons)) %>% dplyr::filter(!duplicated(queryHits))
  df = df[hits$queryHits,]
  df
}

openFile = function(fname) {
  f <- NULL
  if (is.null(fname)) {
    f <- file("stdin")
  } else {
    f <- file(fname, open="r")
  }
  f
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
