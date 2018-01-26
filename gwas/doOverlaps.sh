#!/bin/bash

F=$1
fullname=`basename $F`
NAME="${fullname%.*}"

cat $F | sed '1d' | awk 'BEGIN{OFS="\t"}{print $2,$3-1,$3,$2"_"$3}' > $NAME.input.bed
cat $F | sed '1d' | awk 'BEGIN{OFS="\t"}{print "chr"$2,$3-1,$3,$2"_"$3}' > $NAME.input.chr.bed

multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.dnase.txt --output roadmap/dnase -vv
multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.dnase.brain.txt --output roadmap/dnase.brain -vv
multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.dnase.blood_tcell.txt --output roadmap/dnase.blood_tcell -vv
multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.dnase.hsc_bcell.txt --output roadmap/dnase.hsc_bcell -vv

multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.H3K27ac.txt --output roadmap/H3K27ac -vv
multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.H3K4me3.txt --output roadmap/H3K4me3 -vv
multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.FantomEnh.txt --output roadmap/FantomEnh -vv

multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.RoadmapEnh.txt --output roadmap/RoadmapEnh -vv
multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.RoadmapEnh.brain.txt --output roadmap/RoadmapEnh.brain -vv
multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.RoadmapEnh.blood_tcell.txt --output roadmap/RoadmapEnh.blood_tcell -vv
multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.RoadmapEnh.hsc_bcell.txt --output roadmap/RoadmapEnh.hsc_bcell -vv

multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist ./ipsneurons/ipsc.atac.bedfile.txt --output ipsneurons/ipsc.atac -vv
multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist ./ipsneurons/npc.atac.bedfile.txt --output ipsneurons/npc.atac -vv
multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist ./ipsneurons/neuron.atac.bedfile.txt --output ipsneurons/neuron.atac -vv
multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist ./ipsneurons/ineuron.atac.bedfile.txt --output ipsneurons/ineuron.atac -vv

multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist ./ipsneurons/ipsMacrophage.atac.bedfile.txt --output ipsneurons/ipsMacrophage.atac -vv


# Merge together all the results
paste <(head -n 1 $F) \
      <(echo -e "iPSC\tipsMacrophage\tNPC\tNeuron\tiNeuron\tDNase\tBrain DNase\tTcell DNase\tBcell DNase\tFantom Enh\tRoadmapEnh\tBrain Enh\tTcell Enh\tBcell Enh") \
      <(echo -e "DNase overlaps\tBrain DNase overlaps\tTcell DNase overlaps\tBcell DNase overlaps\tFantom Enh overlaps\tRoadmapEnh overlaps\tBrain Enh overlaps\tTcell Enh overlaps\tBcell Enh overlaps") \
      > $NAME.annotated.txt
paste <(sed '1d' $F) \
      <(sed '1d' ipsneurons/ipsc.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' ipsneurons/ipsMacrophage.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' ipsneurons/npc.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' ipsneurons/neuron.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' ipsneurons/ineuron.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' roadmap/dnase.overlap.summary.txt | cut -f 4) \
      <(sed '1d' roadmap/dnase.brain.overlap.summary.txt | cut -f 4) \
      <(sed '1d' roadmap/dnase.blood_tcell.overlap.summary.txt | cut -f 4) \
      <(sed '1d' roadmap/dnase.hsc_bcell.overlap.summary.txt | cut -f 4) \
      <(sed '1d' roadmap/FantomEnh.overlap.summary.txt | cut -f 4) \
      <(sed '1d' roadmap/RoadmapEnh.overlap.summary.txt | cut -f 4) \
      <(sed '1d' roadmap/RoadmapEnh.brain.overlap.summary.txt | cut -f 4) \
      <(sed '1d' roadmap/RoadmapEnh.blood_tcell.overlap.summary.txt | cut -f 4) \
      <(sed '1d' roadmap/RoadmapEnh.hsc_bcell.overlap.summary.txt | cut -f 4) \
      <(sed '1d' roadmap/dnase.overlap.summary.txt | cut -f 5) \
      <(sed '1d' roadmap/dnase.brain.overlap.summary.txt | cut -f 5) \
      <(sed '1d' roadmap/dnase.blood_tcell.overlap.summary.txt | cut -f 5) \
      <(sed '1d' roadmap/dnase.hsc_bcell.overlap.summary.txt | cut -f 5) \
      <(sed '1d' roadmap/FantomEnh.overlap.summary.txt | cut -f 5) \
      <(sed '1d' roadmap/RoadmapEnh.overlap.summary.txt | cut -f 5) \
      <(sed '1d' roadmap/RoadmapEnh.brain.overlap.summary.txt | cut -f 5) \
      <(sed '1d' roadmap/RoadmapEnh.blood_tcell.overlap.summary.txt | cut -f 5) \
      <(sed '1d' roadmap/RoadmapEnh.hsc_bcell.overlap.summary.txt | cut -f 5) \
      >> $NAME.annotated.txt
