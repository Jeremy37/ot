#!/usr/bin/env python
import argparse
import os

filepath = "/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/datasets/SpliceAI/spliceai.merge.tsv.bgz"

for linestr in open("AD.IGAP1_GWAX.gwsig_indep.hits",'r'):
    linevals = linestr.split()
    if linevals[0] == "Chr":
        continue
    chrom = linevals[0]
    pos = linevals[2]
    os.system("tabix {} {}:{}-{} >> spliceai/ad.spliceai.regions.tsv".format(filepath, chrom, (int(pos) - 100000), (int(pos) + 100000)) )
