#!/usr/bin/env python

import argparse
import sys
import os
import os.path
import subprocess
import uuid
from pprint import pprint

parser = argparse.ArgumentParser(description="Analyzes the deletion patterns in gDNA and cDNA from a CRISPR-Cas9 experiment.")
parser.add_argument("--input", type=file, required=True, metavar='FILE', help="Input file name, with each line listing relevant BAM files and regions")
parser.add_argument("--out", type=str, required=True, metavar='FILE', help="Base path for output files")
parser.add_argument("--genome", type=str, required=True, metavar='FILE', help="Path to indexed fasta file for the reference genome")
parser.add_argument("--saveall", action='store_true', help="If true, then extra files are saved with region-aligned and discarded reads for gDNA and cDNA")
parser.add_argument("--minMapQ", type=int, default=0, help="Minimum mapping quality for read to be included")
parser.add_argument("--subsample", type=float, default=1, help="Subsample reads to the specified fraction (e.g. for testing)")
parser.add_argument("--sep", type=str, metavar='SEP[\\t]', default="\t", help="field separator (default: tab)")
parser.add_argument("--debug", action='store_true')

args = parser.parse_args()


def main():
	statsOutputFname = args.out + ".stats.tsv"
# 	with open(statsOutputFname, "w") as outfile:
# 		outfile.write("\t".join(['region_name', 'gDNA_num_reads', 'gDNA_num_softclipped', 'gDNA_num_hardclipped', 'gDNA_num_insertion', 'gDNA_aligned', 'gDNA_discarded', 'cDNA_num_reads', 'cDNA_num_softclipped', 'cDNA_num_hardclipped', 'cDNA_num_insertion', 'cDNA_aligned', 'cDNA_discarded', 'num_cigars']) + "\n")
# 
	refSequences = {}
	linenum = 1
	for linestr in args.input:
		if linestr.startswith("#"):
			continue
		lineVals = linestr.strip('\n').split('\t')
		region_name = lineVals[0]
		gdna_file = lineVals[1]
		cdna_file = lineVals[2]
		chr = lineVals[3]
		roi_start = int(lineVals[4])
		roi_end = int(lineVals[5])
		roi_site = int(lineVals[6])
		hdr_allele = lineVals[7]
		wt_allele = lineVals[8]
		ref_sequence = lineVals[9]
		print "Input line:"
		print linestr
		
		if not ref_sequence:
			# Get the reference genome sequence for the ROI
			region_bedfilename = region_name + "." + str(uuid.uuid4()) + ".bed"
			region_bedfile = open(region_bedfilename, "w")
			region_bedfile.write("%s\t%d\t%d" % (chr, roi_start, roi_end))
			region_bedfile.close()
		
			cmd = "bedtools getfasta -fi %s -bed %s" % (args.genome, region_bedfilename)
			ref_sequence = subprocess.check_output(cmd, shell=True, executable='/bin/bash').splitlines()[1]
			refSequences[region_name] = ref_sequence
			os.remove(region_bedfilename)
		print "ref_sequence:"
		print ref_sequence
		
		if len(ref_sequence) != (roi_end - roi_start):
			print "ERROR: length of reference sequence should be the same as the span from ROI start to end coordinates"
			exit()
		
		#stats = analyzeRegion(region_name, gdna_file, cdna_file, chr, roi_start, roi_end, roi_site, ref_sequence)
		
# 		with open(statsOutputFname, "a") as outfile:
# 			outfile.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" %
# 			              (region_name, stats['gDNA_num_reads'], stats['gDNA_num_softclipped'], stats['gDNA_num_hardclipped'], stats['gDNA_num_insertion'],
# 			               stats['gDNA_aligned'], stats['gDNA_discarded'], stats['cDNA_num_reads'], stats['cDNA_num_softclipped'], stats['cDNA_num_hardclipped'],
# 			               stats['cDNA_num_insertion'], stats['cDNA_aligned'], stats['cDNA_discarded'], stats['num_cigars']))
		
		linenum = linenum + 1
	
	refSequenceOutputFname = args.out + ".ref_sequence.tsv"
	with open(refSequenceOutputFname, "w") as refseq_outfile:
		refseq_outfile.write("\t".join(['region_name', 'ref_sequence']) + "\n")
		for key, value in refSequences.items():
			refseq_outfile.write("\t".join([key, value]) + "\n")


def analyzeRegion(region_name, gdna_file, cdna_file, chr, roi_start, roi_end, roi_site, ref_sequence):
	stats = {}
	# Extract reads from gDNA BAM file
	alignmentFlag = "-F 0x904 " # only include primary alignment for each read, and must be mapped
	mapqFlag = "-q %d " % args.minMapQ
	subsampleStr = ""
	if args.subsample and args.subsample < 1:
		subsampleStr = "-s %f " % args.subsample
	regionStr = chr + ":%d-%d" % (roi_start, roi_end)
	cmd = "samtools view " + alignmentFlag + mapqFlag + subsampleStr + gdna_file + " " + regionStr
	#cmd = "samtools view " + alignmentFlag + mapqFlag + subsampleStr + gdna_file + " " + regionStr + " | head -n 50"
	print cmd
	sam_output = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
	#print "output:"
	#print sam_output
	
	ref_coords = (roi_start, roi_end)
	ZONE = (0, len(ref_sequence))
	gDNA_data = getAlignedReads(sam_output, ref_sequence, ref_coords)
	gDNA_aligned = gDNA_data['alignedReads']
	print "%d gDNA reads of %d total (%.1f%%) were soft-clipped" % (gDNA_data['num_softclipped'], gDNA_data['num_reads'], 100.0 * gDNA_data['num_softclipped'] / gDNA_data['num_reads'])
	print "%d gDNA reads of %d total (%.1f%%) were hard-clipped" % (gDNA_data['num_hardclipped'], gDNA_data['num_reads'], 100.0 * gDNA_data['num_hardclipped'] / gDNA_data['num_reads'])
	print "%d gDNA reads of %d total (%.1f%%) had insertions" % (gDNA_data['num_insertion'], gDNA_data['num_reads'], 100.0 * gDNA_data['num_insertion'] / gDNA_data['num_reads'])
	stats['gDNA_num_softclipped'] = gDNA_data['num_softclipped']
	stats['gDNA_num_hardclipped'] = gDNA_data['num_hardclipped']
	stats['gDNA_num_insertion'] = gDNA_data['num_insertion']
	stats['gDNA_num_reads'] = gDNA_data['num_reads']
	
	# Now that we have the reads all starting at the same position, we can determine
	# their cigar strings relative to this region
	gDNA_cigars = getZoneCigars(gDNA_aligned, ref_sequence)
	
	stats['gDNA_aligned'] = len(gDNA_aligned)
	stats['gDNA_discarded'] = len(gDNA_data['discardedReads'])
	if args.saveall:
		# Save aligned output to file
		alignedFile = args.out + "." + region_name + ".gDNA.aligned.txt"
		print "Saving %d gDNA aligned reads to: %s" % (len(gDNA_aligned), alignedFile)
		with open(alignedFile, "w") as outFile:
			for i in xrange(len(gDNA_aligned)):
				print>>outFile, ("%s\t%s" % (gDNA_aligned[i], gDNA_cigars[i]))
		discardFile = args.out + "." + region_name + ".gDNA.discarded.txt"
		print "Saving %d gDNA discarded reads to: %s" % (len(gDNA_data['discardedReads']), discardFile)
		with open(discardFile, "w") as outFile:
			for line in gDNA_data['discardedReads']:
				print>>outFile, line
	
	# Extract reads from cDNA BAM file
	cmd = "samtools view " + alignmentFlag + mapqFlag + subsampleStr + cdna_file + " " + regionStr
	#cmd = "samtools view " + alignmentFlag + mapqFlag + subsampleStr + cdna_file + " " + regionStr + " | head -n 50"
	print cmd
	sam_output = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
	#print "output:"
	#print sam_output
	cDNA_data = getAlignedReads(sam_output, ref_sequence, ref_coords)
	cDNA_aligned = cDNA_data['alignedReads']
	print "%d cDNA reads of %d total (%.1f%%) were soft-clipped" % (cDNA_data['num_softclipped'], cDNA_data['num_reads'], 100.0 * cDNA_data['num_softclipped'] / cDNA_data['num_reads'])
	print "%d cDNA reads of %d total (%.1f%%) were hard-clipped" % (cDNA_data['num_hardclipped'], cDNA_data['num_reads'], 100.0 * cDNA_data['num_hardclipped'] / cDNA_data['num_reads'])
	print "%d cDNA reads of %d total (%.1f%%) had insertions" % (cDNA_data['num_insertion'], cDNA_data['num_reads'], 100.0 * cDNA_data['num_insertion'] / cDNA_data['num_reads'])
	stats['cDNA_num_softclipped'] = cDNA_data['num_softclipped']
	stats['cDNA_num_hardclipped'] = cDNA_data['num_hardclipped']
	stats['cDNA_num_insertion'] = cDNA_data['num_insertion']
	stats['cDNA_num_reads'] = cDNA_data['num_reads']
	
	cDNA_cigars = getZoneCigars(cDNA_aligned, ref_sequence)
	
	stats['cDNA_aligned'] = len(cDNA_aligned)
	stats['cDNA_discarded'] = len(cDNA_data['discardedReads'])
	if args.saveall:
		alignedFile = args.out + "." + region_name + ".cDNA.aligned.txt"
		print "Saving %d cDNA aligned reads to: %s" % (len(cDNA_aligned), alignedFile)
		with open(alignedFile, "w") as outFile:
			for i in xrange(len(cDNA_aligned)):
				print>>outFile, ("%s\t%s" % (cDNA_aligned[i], cDNA_cigars[i]))
		discardFile = args.out + "." + region_name + ".cDNA.discarded.txt"
		print "Saving %d cDNA discarded reads to: %s" % (len(cDNA_data['discardedReads']), discardFile)
		with open(discardFile, "w") as outFile:
			for line in cDNA_data['discardedReads']:
				print>>outFile, line
	
	gDNA_cigars_table = getCigarsTable(gDNA_cigars, gDNA_aligned)
	cDNA_cigars_table = getCigarsTable(cDNA_cigars, cDNA_aligned)
	cigar_data = matchCigars(gDNA_cigars_table, cDNA_cigars_table)
	merged_cigars = cigar_data['CIGARs']
	merged_cigar_reads = cigar_data['CIGAR_reads']
	cigar_count_gDNA = cigar_data['CIGARs_COUNT_gDNA']
	cigar_count_cDNA = cigar_data['CIGARs_COUNT_cDNA']
	cigar_mutation_count = getCigarMutationCounts(cigar_data['CIGARs'])
	
	outputFname = args.out + "." + region_name + ".cigars.tsv"
	# Get indices for order sorted by decreasing gDNA cigar count
	sortedIndices = sorted(range(len(cigar_count_gDNA)), key=lambda k: -cigar_count_gDNA[k])
	with open(outputFname, "w") as outfile:
		outfile.write("cigar\tregion_read\tgDNA_count\tcDNA_count\tmismatch_count\n")
		for i in sortedIndices:
			outfile.write("\t".join([merged_cigars[i], merged_cigar_reads[i], str(cigar_count_gDNA[i]), str(cigar_count_cDNA[i]), str(cigar_mutation_count[i])]) + "\n")
	
	stats['num_cigars'] = len(merged_cigars)
	return(stats)



# Goes through each SAM file read and gets its sequence with respect to the reference
# zone of interest, where "-" indicates that the read does not cover the position, and
# * indicates a deletion.
def getAlignedReads(reads_input, ref_sequence, ref_coords):
	# reads_input: array of lines from a SAM file
	# ref_sequence: sequence of the amplicon
	# ref_coords: coordinates of the region of interest in the reference
	
	ref_seq_length = len(ref_sequence) # size of the ref genome
	alignedReads = list()
	discardedReads = list()
	num_reads = 0
	num_softclipped = 0
	num_hardclipped = 0
	num_insertion = 0
	num_outsidewindow = 0
	for i,line in enumerate(reads_input.splitlines()): # reads line by line
		line_split = line.split("\t") # converts the string "line" to an array of string
		num_reads += 1
		
		if len(line_split) > 10: # avoid the .SAM header
			sam_position = int(line_split[3]) - 1  # python indexes start at 0
			relative_pos = sam_position - ref_coords[0]
			if abs(relative_pos) > 1000:
				num_outsidewindow += 1
				print "Likely ERROR: read position and reference coordinates differ by more than 1000 bp"
				exit()
			
			sam_cigar = line_split[5]
			sam_read = line_split[9]
			if sam_cigar == "*":
				discardedReads.append(line)
				continue # Ignore this read, as it doesn't have a valid CIGAR
			
			# computes the read sequence within the coords of interest
			output = registerRead(relative_pos, sam_cigar, sam_read, ref_seq_length)
			if output['read_ok']:
				alignedReads.append(output['read'])
			else:
				discardedReads.append(line)
			if output['read_status'] == "hardclipped":
				num_hardclipped += 1
			elif output['read_status'] == "insertion":
				num_insertion += 1
			elif output['read_status'] == "softclipped":
				num_softclipped += 1
			#print "Registered read:"
			#print(alignedReads[i])			
	
	return {'alignedReads':alignedReads, 'discardedReads':discardedReads, 'num_reads':num_reads,
			'num_hardclipped':num_hardclipped, 'num_insertion':num_insertion,
			'num_softclipped':num_softclipped, 'num_outsidewindow':num_outsidewindow}


# this function takes a .sam read and registers to the reference genome, relative
# to the region of interest in the reference
# '-' not sequenced nucleotides
# '*' deleted nucleotides
# if the cigars contains H or I the variable 'read_ok' will be False
# the returned read should be the same length as the reference genome
def registerRead(relative_position, sam_cigar, sam_read, region_length):
	sam_cigar_parsed = parseCigar(sam_cigar)
	read_ok = True
	read_contains_deletion = False
	read_status = ""
	expanded_cigar = ""
	cursor_read = 0
	
	for i in range(0, len(sam_cigar_parsed['list_letter'])):
		
		type = sam_cigar_parsed['list_letter'][i]
		
		if type == 'M':
			expanded_cigar += sam_read[cursor_read: cursor_read + sam_cigar_parsed['list_num'][i]]
			cursor_read += sam_cigar_parsed['list_num'][i]
		elif type == 'D':
			expanded_cigar += '*' * sam_cigar_parsed['list_num'][i]
			read_contains_deletion = True
		elif type == 'S':
			cursor_read += sam_cigar_parsed['list_num'][i]
			read_status = "softclipped"
		elif type == 'H':
			read_status = "hardclipped"
			read_ok = False
			break
		elif type == 'I':
			read_status = "insertion"
			read_ok = False
			break
		else:
			read_ok = False
			break
	
	registered_read = ""
	if read_ok:
		if relative_position < 0:
			if (len(expanded_cigar) - relative_position) > 0:
				registered_read = expanded_cigar[-relative_position:]
		else:
			registered_read = '-' * relative_position
			registered_read += expanded_cigar
	
		if len(registered_read) > region_length:
			registered_read = registered_read[0:region_length]
		else:
			registered_read += '-' * (region_length - len(registered_read))
	
	#print registered_read
	return{'read':registered_read, 'read_ok':read_ok, 'read_status':read_status, 'deletion':read_contains_deletion}


# this function takes a string cigar "10M20D7M" and converts it to arrays of letters ['M','D','M'] and number [10,20,7]
def parseCigar(input_cigar):
	#print(input_cigar)
	X = list(input_cigar) # splits the string into array of char
	
	list_num = []
	list_letter = []
	
	idx_last_letter = -1;
	idx_current = 0;
	
	for x in X:
		if not is_number(x):
			list_letter.append(x)
			idx_1 = idx_last_letter + 1
			idx_2 = idx_current
			list_num.append(int(''.join(X[idx_1:idx_2])))
			idx_last_letter = idx_current
		idx_current += 1
	
	return {'list_letter': list_letter, 'list_num': list_num, 'canonical': set(list_letter) == set(['M', 'D', 'M'])}


# Goes through each "registered" read, which has been altered to begin at the start
# of the reference coords, and builds a custom cigar string relative to this reference
def getZoneCigars(registered_reads, ref_sequence):
	return [mreCRISPR_cigar_in_zone(ref_sequence, regread) for regread in registered_reads]


def getCigarsTable(zone_cigars, registered_reads):
	CIGARs = list() # list of unique CIGARs
	CIGARs_count = list() # corresponding count for each CIGAR
	CIGAR_reads = list() # read corresponding to cigar
	
	for i,cigar in enumerate(zone_cigars):
		# computes custom cigar for the registered read:
		if cigar in CIGARs:
			CIGARs_count[CIGARs.index(cigar)] += 1
		else:
			CIGARs.append(cigar)
			CIGARs_count.append(1)
			CIGAR_reads.append(registered_reads[i])
	
	return{'CIGARs':CIGARs,
		   'CIGARs_count':CIGARs_count,
		   'CIGAR_reads':CIGAR_reads}


# This function computes the cigar in the zone of interest
# example output:
# cigar ['+','D','-']
# cigar_count [10,7,23]
def mreCRISPR_cigar_in_zone(ref_zone_seq, registered_read_zone):
	cigar = list() # cigar letters
	cigar_count = list() # cigar counts
	N = len(ref_zone_seq)
	type = ''
	count = 0
	for i in range(0,N):
		if registered_read_zone[i] == ref_zone_seq[i]:
			if type == '':
				type = '+'
				count = 1
			elif type == '+':
				count += 1
			else:
				cigar.append(type)
				cigar_count.append(count)
				type = '+'
				count = 1
		elif registered_read_zone[i] == '*':
			if type == '':
				type = 'D'
				count = 1
			elif type == 'D':
				count += 1
			else:
				cigar.append(type)
				cigar_count.append(count)
				type = 'D'
				count = 1
		elif registered_read_zone[i] == '-':
			if type == '':
				type = '#'
				count = 1
			elif type == '#':
				count += 1
			else:
				cigar.append(type)
				cigar_count.append(count)
				type = '#'
				count = 1
		else:
			if type == '':
				type = '-'
				count = 1
			elif type == '-':
				count += 1
			else:
				cigar.append(type)
				cigar_count.append(count)
				type = '-'
				count = 1
	cigar.append(type)
	cigar_count.append(count)
	
	cigar_txt = ''
	full_read = True
	
	for i in range(0,len(cigar)):
		cigar_txt += str(cigar_count[i]) + cigar[i]

		if cigar[i] == '#':
			full_read = False
	
	# print(cigar,cigar_count,cigar_txt)
# 	return{'cigar':cigar,
# 		   'cigar_count':cigar_count,
# 		   'cigar_txt':cigar_txt,
# 		   'full_read':full_read}
	return cigar_txt


# test if a character is a number
def is_number(s):
	try:
		int(s)
		return True
	except ValueError:
		return False


# this function pairs cDNA and gDNA and computes mirScore
def matchCigars(gDNA_cigar_table, cDNA_cigar_table):
	# computes overlap and sort by either mirScore or SeqDepth
	CIGARs = list() # contains the list of all unique cigars over the zone
	CIGAR_reads = list()
	CIGARs_COUNT_gDNA = list() 
	CIGARs_COUNT_cDNA = list()
	
	# finds the complete set:
	TOTAL_CIGARs = []
	for i in range(0,len(gDNA_cigar_table['CIGARs'])):
		if not gDNA_cigar_table['CIGARs'][i] in CIGARs:
			CIGARs.append(gDNA_cigar_table['CIGARs'][i])
			CIGAR_reads.append(gDNA_cigar_table['CIGAR_reads'][i])
	
	for i in range(0,len(cDNA_cigar_table['CIGARs'])):
		if not cDNA_cigar_table['CIGARs'][i] in CIGARs:
			CIGARs.append(cDNA_cigar_table['CIGARs'][i])
			CIGAR_reads.append(cDNA_cigar_table['CIGAR_reads'][i])
	
	print('total cigars:%d\tunique cigars:%d' % (len(gDNA_cigar_table['CIGARs']) + len(cDNA_cigar_table['CIGARs']), len(CIGARs)))
	
	# fills out counts
	for i in range(0,len(CIGARs)):
		if CIGARs[i] in gDNA_cigar_table['CIGARs']:
			idx_gDNA = gDNA_cigar_table['CIGARs'].index(CIGARs[i])
			CIGARs_COUNT_gDNA.append(gDNA_cigar_table['CIGARs_count'][idx_gDNA])
		else:
			CIGARs_COUNT_gDNA.append(0)
		
		if CIGARs[i] in cDNA_cigar_table['CIGARs']:
			idx_cDNA = cDNA_cigar_table['CIGARs'].index(CIGARs[i])
			CIGARs_COUNT_cDNA.append(cDNA_cigar_table['CIGARs_count'][idx_cDNA])
		else:
			CIGARs_COUNT_cDNA.append(0)
	
	return{'CIGARs':CIGARs,
		   'CIGAR_reads':CIGAR_reads,
		   'CIGARs_COUNT_gDNA':CIGARs_COUNT_gDNA,
		   'CIGARs_COUNT_cDNA':CIGARs_COUNT_cDNA}


def cigarMutationCount(cigar):
	mutation_count = 0
	cigarVals = parseCigar(cigar)
	cigarLetters = cigarVals['list_letter']
	cigarNumbers = cigarVals['list_num']
	for i in xrange(len(cigarLetters)):
		if cigarLetters[i] == '-':
			mutation_count += cigarNumbers[i]
	return mutation_count

def getCigarMutationCounts(cigars):
	mutation_counts = list()
	for cigar in cigars:
		mutation_counts.append(cigarMutationCount(cigar))
	return mutation_counts



main()

