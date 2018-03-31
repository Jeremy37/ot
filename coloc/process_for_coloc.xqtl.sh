#!/bin/bash

JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
cd $JS/reference/xQTL

submitJobs.py --MEM 10000 -j process_xqtl.eQTLs -q yesterday \
   -c "$JS/src/coloc/process_xqtl_file.sh eQTLs_all.txt.gz maf.gz xQTL_eQTL"

submitJobs.py --MEM 10000 -j process_xqtl.mQTLs -q yesterday \
   -c "$JS/src/coloc/process_xqtl_file.sh mQTLs_all.txt.gz maf.gz xQTL_mQTL"

submitJobs.py --MEM 10000 -j process_xqtl.haQTLs -q yesterday \
   -c "$JS/src/coloc/process_xqtl_file.sh haQTLs_all.txt.gz maf.gz xQTL_haQTL"
