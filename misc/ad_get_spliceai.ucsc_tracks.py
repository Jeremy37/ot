#!/usr/bin/env python
import argparse
import os
import gzip

AG_file = open("spliceai.ucsc.ag.wig",'w')
AL_file = open("spliceai.ucsc.al.wig",'w')
DG_file = open("spliceai.ucsc.dg.wig",'w')
DL_file = open("spliceai.ucsc.dl.wig",'w')
merged_file = open("spliceai.ucsc.merged.wig",'w')

ag_color = "255,0,0"
al_color = "255,0,0"
dg_color = "0,110,165"
dl_color = "50,80,255"
merged_color = "0,0,200"

print >>AG_file, "track type=wiggle_0 name=\"SpliceAI Acceptor Gain\" color={}".format(ag_color)
print >>AL_file, "track type=wiggle_0 name=\"SpliceAI Acceptor Loss\" color={}".format(al_color)
print >>DG_file, "track type=wiggle_0 name=\"SpliceAI Donor Gain\" color={}".format(dg_color)
print >>DL_file, "track type=wiggle_0 name=\"SpliceAI Donor Loss\" color={}".format(dl_color)
print >>merged_file, "track type=wiggle_0 name=\"SpliceAI Max\" color={}".format(merged_color)

lastchr = "X"
for linestr in gzip.open("ad.spliceai.by_pos.tsv.gz",'r'):
	linevals = linestr.split()
	chr = linevals[0]
	pos = linevals[1]
	# Skip header line
	if chr == "chr":
		continue
	
	if chr != lastchr:
		# Write wiggle variableStep line
		print >>AG_file, "variableStep chrom=chr{}".format(chr)
		print >>AL_file, "variableStep chrom=chr{}".format(chr)
		print >>DG_file, "variableStep chrom=chr{}".format(chr)
		print >>DL_file, "variableStep chrom=chr{}".format(chr)
		print >>merged_file, "variableStep chrom=chr{}".format(chr)
		lastchr = chr
	
	print >>AG_file, "{}\t{}".format(pos, linevals[2])
	print >>AL_file, "{}\t{}".format(pos, linevals[3])
	print >>DG_file, "{}\t{}".format(pos, linevals[4])
	print >>DL_file, "{}\t{}".format(pos, linevals[5])
	print >>merged_file, "{}\t{}".format(pos, max(float(linevals[2]), float(linevals[3]), float(linevals[4]), float(linevals[5])) )
