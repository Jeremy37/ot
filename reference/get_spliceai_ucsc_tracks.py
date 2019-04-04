#!/usr/bin/env python
import argparse
import os
import gzip

AG_file = open("spliceai.ucsc.ag.wig",'w')
AL_file = open("spliceai.ucsc.al.wig",'w')
DG_file = open("spliceai.ucsc.dg.wig",'w')
DL_file = open("spliceai.ucsc.dl.wig",'w')
merged_file = open("spliceai.ucsc.merged.wig",'w')

def main():
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

	lastpos = -1
	lastchr = "Z"
	pos_lines = []
	for linestr in gzip.open("spliceai.merge.tsv.bgz",'r'):
		linevals = linestr.split()
		chr = linevals[0]
		pos = linevals[1]
		# Skip header line
		if chr == "chr":
			continue

		if (pos != lastpos and lastpos != -1):
			# Write wiggle variableStep lines
			processPos(chr, lastpos, pos_lines)
			pos_lines = [linevals]
		else:
			pos_lines.append(linevals)
		lastpos = pos
		
		if chr != lastchr:
			# Write wiggle variableStep line
			print >>AG_file, "variableStep chrom=chr{}".format(chr)
			print >>AL_file, "variableStep chrom=chr{}".format(chr)
			print >>DG_file, "variableStep chrom=chr{}".format(chr)
			print >>DL_file, "variableStep chrom=chr{}".format(chr)
			print >>merged_file, "variableStep chrom=chr{}".format(chr)
			lastchr = chr
		


def processPos(chr, pos, lines):
	max_ag = max_al = max_dg = max_dl = 0
	for line in lines:
		ag = float(line[8])
		al = float(line[9])
		dg = float(line[10])
		dl = float(line[11])
		if ag > max_ag:
			max_ag = ag
		if al > max_al:
			max_al = al
		if dg > max_dg:
			max_dg = dg
		if dl > max_dl:
			max_dl = dl
	print >>AG_file, "{}\t{}".format(pos, max_ag)
	print >>AL_file, "{}\t{}".format(pos, max_al)
	print >>DG_file, "{}\t{}".format(pos, max_dg)
	print >>DL_file, "{}\t{}".format(pos, max_dl)
	print >>merged_file, "{}\t{}".format(pos, max(max_ag, max_al, max_dg, max_dl) )
	

main()

