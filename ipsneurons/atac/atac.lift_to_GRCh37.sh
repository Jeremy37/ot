JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys

cd $JS/ipsneurons/GRCh37/ATAC

CrossMap.py bed $JS/software/CrossMap/GRCh38_to_GRCh37.chain.gz \
    $JS/ipsneurons/GRCh38/ATAC/peaks/atac_multisample_peaks.narrowPeak \
    $JS/ipsneurons/GRCh37/ATAC/peaks/atac_multisample_peaks.narrowPeak.tmp

# Somehow CrossMap removes the "chr" from the bed file, so we add it back in
cat $JS/ipsneurons/GRCh37/ATAC/peaks/atac_multisample_peaks.narrowPeak.tmp | perl -ne 'print "chr".$_."\n"' \
   > $JS/ipsneurons/GRCh37/ATAC/peaks/atac_multisample_peaks.narrowPeak


####### Do the same for subsets of the ATAC samples
CrossMap.py bed $JS/software/CrossMap/GRCh38_to_GRCh37.chain.gz $JS/ipsneurons/GRCh38/ATAC/peaks/atac_ipsc_peaks.narrowPeak $JS/ipsneurons/GRCh37/ATAC/peaks/atac_ipsc_peaks.narrowPeak.tmp
cat $JS/ipsneurons/GRCh37/ATAC/peaks/atac_ipsc_peaks.narrowPeak.tmp | perl -ne 'print "chr".$_."\n"' \
   > $JS/ipsneurons/GRCh37/ATAC/peaks/atac_ipsc_peaks.narrowPeak

CrossMap.py bed $JS/software/CrossMap/GRCh38_to_GRCh37.chain.gz $JS/ipsneurons/GRCh38/ATAC/peaks/atac_npc_peaks.narrowPeak $JS/ipsneurons/GRCh37/ATAC/peaks/atac_npc_peaks.narrowPeak.tmp
cat $JS/ipsneurons/GRCh37/ATAC/peaks/atac_npc_peaks.narrowPeak.tmp | perl -ne 'print "chr".$_."\n"' \
   > $JS/ipsneurons/GRCh37/ATAC/peaks/atac_npc_peaks.narrowPeak

CrossMap.py bed $JS/software/CrossMap/GRCh38_to_GRCh37.chain.gz $JS/ipsneurons/GRCh38/ATAC/peaks/atac_neuron_peaks.narrowPeak $JS/ipsneurons/GRCh37/ATAC/peaks/atac_neuron_peaks.narrowPeak.tmp
cat $JS/ipsneurons/GRCh37/ATAC/peaks/atac_neuron_peaks.narrowPeak.tmp | perl -ne 'print "chr".$_."\n"' \
   > $JS/ipsneurons/GRCh37/ATAC/peaks/atac_neuron_peaks.narrowPeak

CrossMap.py bed $JS/software/CrossMap/GRCh38_to_GRCh37.chain.gz $JS/ipsneurons/GRCh38/ATAC/peaks/atac_ineuron_peaks.narrowPeak $JS/ipsneurons/GRCh37/ATAC/peaks/atac_ineuron_peaks.narrowPeak.tmp
cat $JS/ipsneurons/GRCh37/ATAC/peaks/atac_ineuron_peaks.narrowPeak.tmp | perl -ne 'print "chr".$_."\n"' \
   > $JS/ipsneurons/GRCh37/ATAC/peaks/atac_ineuron_peaks.narrowPeak

rm $JS/ipsneurons/GRCh37/ATAC/peaks/*.tmp
