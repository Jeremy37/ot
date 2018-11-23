#!/usr/bin/env python
import os
import argparse
import re

parser = argparse.ArgumentParser(description = "Extract parameters and output results from running Flash, and output them as a single tab-separated line of values", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--file", type=file, required=True, metavar='FILE', help = "File with output from running Flash")
args = parser.parse_args()

# Lines should have trhee tab-separated fields; the first is the output sample name,
# and the second and third fields are the names of fastq file 1 and fastq file 2,
# which must have paired reads in the same order (i.e. matching line by line).
# Optional fields 4,5,6 give the expected fragment length, read length, and
# fragment length standard deviation.
filedata = args.file.read()

res = re.search(r'/([^/]+)\.extendedFrags.fastq.gz', filedata)
output_name = res.group(1)

res = re.search(r'Min overlap:[\s]*([\d]+)', filedata)
min_overlap = res.group(1)

res = re.search(r'Max overlap:[\s]*([\d]+)', filedata)
max_overlap = res.group(1)

res = re.search(r'Max mismatch density:[\s]*([\d\.]+)', filedata)
max_mismatch_dens = res.group(1)

res = re.search(r'Total pairs:[\s]*([\d]+)', filedata)
total_pairs = res.group(1)

res = re.search(r'Combined pairs:[\s]*([\d]+)', filedata)
combined_pairs = res.group(1)

res = re.search(r'Uncombined pairs:[\s]*([\d]+)', filedata)
uncombined_pairs = res.group(1)

res = re.search(r'Percent combined:[\s]*([\d\.%]+)', filedata)
pct_combined = res.group(1)

res = re.search(r'WARNING:[\s]*([^[]+)', filedata)
warning = ""
if (res):
	warning = res.group(1)
	warning = re.sub("\n", " ", warning)

print("\t".join(["File name", "Output name", "Total pairs", "Combined pairs", "Uncombined pairs", "Percent combined", "Min overlap", "Max overlap", "Max mismatch dens", "Warning"]))
print("\t".join([args.file.name, output_name, total_pairs, combined_pairs, uncombined_pairs, pct_combined, min_overlap, max_overlap, max_mismatch_dens, warning]))
