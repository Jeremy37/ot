JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

cd $JS/datasets/blueprint

submitJobs.py --MEM 500 -j process_bluprint.mono_gene -q yesterday \
   -c "$JS/src/coloc/process_blueprint_file.sh mono_gene_nor_combat_peer_10_all_summary.txt.gz mono_gene_eQTL"

submitJobs.py --MEM 500 -j process_bluprint.neut_gene -q yesterday \
   -c "$JS/src/coloc/process_blueprint_file.sh neut_gene_nor_combat_peer_10_all_summary.txt.gz neut_gene_eQTL"

submitJobs.py --MEM 500 -j process_bluprint.tcel_gene -q yesterday \
   -c "$JS/src/coloc/process_blueprint_file.sh tcel_gene_nor_combat_peer_10_all_summary.txt.gz tcel_gene_eQTL"


submitJobs.py --MEM 1000 -j process_bluprint.mono_K27AC -q yesterday \
   -c "$JS/src/coloc/process_blueprint_file.sh mono_K27AC_log2rpm_peer_10_all_summary.txt.gz mono_K27AC"

submitJobs.py --MEM 1000 -j process_bluprint.mono_K4ME1 -q yesterday \
   -c "$JS/src/coloc/process_blueprint_file.sh mono_K4ME1_log2rpm_peer_10_all_summary.txt.gz mono_K4ME1"

