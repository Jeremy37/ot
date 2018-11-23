OTCOREGEN=/lustre/scratch115/realdata/mdt3/projects/otcoregen
JS=$OTCOREGEN/jeremys

JS_SRC=/nfs/users/nfs_j/js29/src

md $JS/experiment/macrophage/sequence_data

python $JS_SRC/utils/irods/irodsGetSamplesInStudy.py --studyName "OpenTargets AD PD finemap RNA" > cram.file_list.tsv
cat cram.file_list.tsv | sed 's/.cram//g' > irods.lanelets.txt

# Map lanelet ids to sample ids
python ~/src/utils/irods/irodsFetchMeta.py --irodsList irods.lanelets.txt | sort -k1 > irods.sample_lanes.txt 

# Fetch lanelets in cram format from irods
cut -f1 irods.lanelets.txt | python ~/src/utils/irods/fetch-irods.py --dir cram/ --suffix .cram
