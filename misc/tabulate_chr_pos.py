#!/usr/bin/env python
import sys
import os
import argparse
import fileinput

parser = argparse.ArgumentParser(description = "Use GATK ASEReadCounter to count allele-specific expression for specific genes.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--mincount", type=int, help = "Min count of reads at position for being printed to output", default = 100)
parser.add_argument("--prefix", type=str, help = "Prefix for each line that is output", default = None)
args = parser.parse_args()

chrPosDict = {}

for linestr in fileinput.input("-"):
	chrPosStr = linestr.strip('\n')
	if chrPosStr in chrPosDict:
		chrPosDict[chrPosStr] = chrPosDict[chrPosStr] + 1
	else:
		chrPosDict[chrPosStr] = 1

for chrPosStr, count in chrPosDict.iteritems():
	if count > args.mincount:
		if args.prefix is not None:
			print "\t".join([args.prefix, chrPosStr, str(count)])
		else:
			print "\t".join([chrPosStr, str(count)])
