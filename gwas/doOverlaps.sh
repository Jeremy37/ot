#!/bin/bash

F=$1
NAME=$2

fullname=`basename $F`
#NAME="${fullname%.*}"
cat $F | sed '1d' | awk 'BEGIN{OFS="\t"}{print $2,$3-1,$3,$2"_"$3}' > $NAME.input.bed
cat $F | sed '1d' | awk 'BEGIN{OFS="\t"}{print "chr"$2,$3-1,$3,$2"_"$3}' > $NAME.input.chr.bed

#multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.dnase.txt --output overlaps/roadmap.dnase -vv
#multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.dnase.brain.txt --output overlaps/roadmap.dnase.brain -vv
#multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.dnase.blood_tcell.txt --output overlaps/roadmap.dnase.blood_tcell -vv
#multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.dnase.hsc_bcell.txt --output overlaps/roadmap.dnase.hsc_bcell -vv
#
#multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.H3K27ac.txt --output overlaps/roadmap.H3K27ac -vv
#multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.H3K4me3.txt --output overlaps/roadmap.H3K4me3 -vv
#multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.FantomEnh.txt --output overlaps/FantomEnh -vv
#
#multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.RoadmapEnh.txt --output overlaps/RoadmapEnh -vv
#multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.RoadmapEnh.brain.txt --output overlaps/RoadmapEnh.brain -vv
#multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.RoadmapEnh.blood_tcell.txt --output overlaps/RoadmapEnh.blood_tcell -vv
#multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.RoadmapEnh.hsc_bcell.txt --output overlaps/RoadmapEnh.hsc_bcell -vv
#
#multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist overlaps/ipsc.atac.bedfile.txt --output overlaps/ipsc.atac -vv
#multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist overlaps/npc.atac.bedfile.txt --output overlaps/npc.atac -vv
#multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist overlaps/neuron.atac.bedfile.txt --output overlaps/neuron.atac -vv
#multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist overlaps/ineuron.atac.bedfile.txt --output overlaps/ineuron.atac -vv
#
#multi_bedintersect.py --snpbed $NAME.input.chr.bed --bedfilelist overlaps/ipsMacrophage.atac.bedfile.txt --output overlaps/ipsMacrophage.atac -vv
multi_bedintersect.py --snpbed $NAME.input.bed --bedfilelist overlaps/ipsSensoryNeuron.atac.bedfile.txt --output overlaps/ipsSensoryNeuron.atac -vv

multi_bedintersect.py --snpbed $NAME.input.bed --bedfilelist overlaps/microglia.atac.bedfile.txt --output overlaps/microglia.atac -vv


# Merge together all the results
paste <(head -n 1 $F) \
      <(echo -e "iPSC\tmicroglia\tipsMacrophage\tNPC\tNeuron\tiNeuron\tipsSensNeuron\tDNase\tBrain DNase\tTcell DNase\tBcell DNase\tFantom Enh\tRoadmapEnh\tBrain Enh\tTcell Enh\tBcell Enh") \
      <(echo -e "DNase overlaps\tBrain DNase overlaps\tTcell DNase overlaps\tBcell DNase overlaps\tFantom Enh overlaps\tRoadmapEnh overlaps\tBrain Enh overlaps\tTcell Enh overlaps\tBcell Enh overlaps") \
      > $NAME.annotated.txt
paste <(sed '1d' $F) \
      <(sed '1d' overlaps/ipsc.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/microglia.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/ipsMacrophage.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/npc.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/neuron.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/ineuron.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/ipsSensoryNeuron.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/roadmap.dnase.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/roadmap.dnase.brain.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/roadmap.dnase.blood_tcell.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/roadmap.dnase.hsc_bcell.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/FantomEnh.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/RoadmapEnh.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/RoadmapEnh.brain.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/RoadmapEnh.blood_tcell.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/RoadmapEnh.hsc_bcell.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/roadmap.dnase.overlap.summary.txt | cut -f 5) \
      <(sed '1d' overlaps/roadmap.dnase.brain.overlap.summary.txt | cut -f 5) \
      <(sed '1d' overlaps/roadmap.dnase.blood_tcell.overlap.summary.txt | cut -f 5) \
      <(sed '1d' overlaps/roadmap.dnase.hsc_bcell.overlap.summary.txt | cut -f 5) \
      <(sed '1d' overlaps/FantomEnh.overlap.summary.txt | cut -f 5) \
      <(sed '1d' overlaps/RoadmapEnh.overlap.summary.txt | cut -f 5) \
      <(sed '1d' overlaps/RoadmapEnh.brain.overlap.summary.txt | cut -f 5) \
      <(sed '1d' overlaps/RoadmapEnh.blood_tcell.overlap.summary.txt | cut -f 5) \
      <(sed '1d' overlaps/RoadmapEnh.hsc_bcell.overlap.summary.txt | cut -f 5) \
      >> $NAME.annotated.txt
