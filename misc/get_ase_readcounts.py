#!/usr/bin/env python
import sys
import os
import argparse
import fileinput
import subprocess

parser = argparse.ArgumentParser(description = "Use GATK ASEReadCounter to count allele-specific expression for specific genes.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--indir", type=str, help = "Base directory of the input BAM files.")
parser.add_argument("--outdir", type=str, help = "Base directory of the output files.")
parser.add_argument("--insuffix", type=str, help = "Suffix of the input bam file.", default = ".bam")
parser.add_argument("--outsuffix", type=str, help = "Suffix of the output ASE count file.", default = ".ASEcounts")
parser.add_argument("--reference", type=file, help = "Path to indexed reference FASTA file.")
parser.add_argument("--vcf", type=str, help = "Path to the VCF file containing SNP coordinates and phased variants.")
parser.add_argument("--geneid", type=str, help = "Ensembl ID of gene to get ASE counts for")
parser.add_argument("--geneidfile", type=file, help = "Path to file with Ensembl gene IDs, one per line")
parser.add_argument("--ensemblfile", type=file, help = "Path to Ensembl gene metadata file", default = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/reference/GRCh38/Homo_sapiens.GRCh38.91.exon_start_end.simple.bed")
parser.add_argument("--execute", help = "Execute the script", default = "False")
parser.add_argument("--Xmx", help = "Memory allocated to the Java process.", default = "1500m")
parser.add_argument("--java_path", help = "Path to the Java executable.", default = "/software/java/bin/java")
parser.add_argument("--gatk_path", help = "Path to the GATK executable.", default = "~ka8/software/GenomeAnalysisTK.jar")
args = parser.parse_args()

# If any gene IDs are specific, then we first subset the VCF to include only
# sites overlapping those genes
geneidHash = {}
if args.geneid != None:
	geneidHash[args.geneid] = True
elif args.geneidfile != None:
	for linestr in args.geneidfile:
		geneidHash[linestr.strip()] = True

vcfFile = args.vcf
if len(geneidHash) > 0:
	cmd = "bcftools view --regions-file args.ensemblfile --output-type z " + args.vcf + " > tmp.vcf.gz"

#Construct file names
for line in fileinput.input("-"):
		sample_name = line.rstrip()
		path_in = os.path.join(args.indir, sample_name, sample_name + args.insuffix)
		path_out = os.path.join(args.outdir, sample_name, sample_name + args.outsuffix)
		gatk_path = " ".join([args.java_path, "-jar", "-Xmx"+args.Xmx, args.gatk_path, "-T ASEReadCounter"])
		flags = "-U ALLOW_N_CIGAR_READS -dt NONE --minMappingQuality 10 -rf MateSameStrand"
		command = " ".join([gatk_path, "-R", args.reference, "-I", path_in, "-o", path_out, "-sites", args.sites, flags])
		print(command)
		if args.execute == "True":
			subprocess.call(['bash','-c',command])
