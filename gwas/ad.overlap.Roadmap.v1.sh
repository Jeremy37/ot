JS=/lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys
ROADMAP=/lustre/scratch117/cellgen/team170/js29/annotations/Roadmap
cd $JS/gwas/AD

# First prepare the input fine mapping file from Toby
# R code
df = read.csv("/lustre/scratch117/cellgen/team170/js29/opentargets/AD/data/TJ.AD.all_causal_v1.csv")
write.table(df, "/lustre/scratch117/cellgen/team170/js29/opentargets/AD/data/TJ.AD.all_causal_v1.orig.txt", quote=F, col.names=T, row.names=F, sep="\t")
##
head -n 1 $JS/opentargets/AD/data/TJ.AD.all_causal_v1.orig.txt > $JS/opentargets/AD/data/TJ.AD.all_causal_v1.txt
sed '1d' $JS/opentargets/AD/data/TJ.AD.all_causal_v1.orig.txt | awk 'BEGIN{OFS="\t"}{print $1,$12"_"$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' \
   >> $JS/opentargets/AD/data/TJ.AD.all_causal_v1.txt
   
F=$JS/opentargets/AD/data/TJ.AD.all_causal_v1.txt
#cat $JS/opentargets/AD/data/TJ.AD.all_causal_v1.tmp.txt | grep -v "HLA" | awk '$1 != "2566"' > $F
#rm $JS/opentargets/AD/data/TJ.AD.all_causal_v1.tmp.txt

cat $F | sed '1d' | awk 'BEGIN{OFS="\t"}{print $3,$4-1,$4,$3"_"$4}' > TJ.AD.input.bed
cat $F | sed '1d' | awk 'BEGIN{OFS="\t"}{print "chr"$3,$4-1,$4,$3"_"$4}' > TJ.AD.input.chr.bed


multi_bedintersect.py --snpbed TJ.AD.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.dnase.txt --output roadmap/dnase -vv
multi_bedintersect.py --snpbed TJ.AD.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.dnase.brain.txt --output roadmap/dnase.brain -vv
multi_bedintersect.py --snpbed TJ.AD.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.dnase.blood_tcell.txt --output roadmap/dnase.blood_tcell -vv
multi_bedintersect.py --snpbed TJ.AD.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.dnase.hsc_bcell.txt --output roadmap/dnase.hsc_bcell -vv

multi_bedintersect.py --snpbed TJ.AD.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.H3K27ac.txt --output roadmap/H3K27ac -vv
multi_bedintersect.py --snpbed TJ.AD.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.H3K4me3.txt --output roadmap/H3K4me3 -vv
multi_bedintersect.py --snpbed TJ.AD.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.FantomEnh.txt --output roadmap/FantomEnh -vv

multi_bedintersect.py --snpbed TJ.AD.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.RoadmapEnh.txt --output roadmap/RoadmapEnh -vv
multi_bedintersect.py --snpbed TJ.AD.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.RoadmapEnh.brain.txt --output roadmap/RoadmapEnh.brain -vv
multi_bedintersect.py --snpbed TJ.AD.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.RoadmapEnh.blood_tcell.txt --output roadmap/RoadmapEnh.blood_tcell -vv
multi_bedintersect.py --snpbed TJ.AD.input.chr.bed --bedfilelist roadmap/roadmap.bedfilelist.RoadmapEnh.hsc_bcell.txt --output roadmap/RoadmapEnh.hsc_bcell -vv

multi_bedintersect.py --snpbed TJ.AD.input.chr.bed --bedfilelist ./ipsneurons/ipsc.atac.bedfile.txt --output ipsneurons/ipsc.atac -vv
multi_bedintersect.py --snpbed TJ.AD.input.chr.bed --bedfilelist ./ipsneurons/npc.atac.bedfile.txt --output ipsneurons/npc.atac -vv
multi_bedintersect.py --snpbed TJ.AD.input.chr.bed --bedfilelist ./ipsneurons/neuron.atac.bedfile.txt --output ipsneurons/neuron.atac -vv
multi_bedintersect.py --snpbed TJ.AD.input.chr.bed --bedfilelist ./ipsneurons/ineuron.atac.bedfile.txt --output ipsneurons/ineuron.atac -vv

# Merge together all the results
paste <(head -n 1 data/TJ.AD.all_causal_v1.txt) \
      <(echo -e "iPSC\tNPC\tNeuron\tiNeuron\tDNase\tBrain DNase\tTcell DNase\tBcell DNase\tFantom Enh\tRoadmapEnh\tBrain Enh\tTcell Enh\tBcell Enh") \
      <(echo -e "DNase overlaps\tBrain DNase overlaps\tTcell DNase overlaps\tBcell DNase overlaps\tFantom Enh overlaps\tRoadmapEnh overlaps\tBrain Enh overlaps\tTcell Enh overlaps\tBcell Enh overlaps") \
      > TJ.AD.all_causal_v1.annotated.txt
paste <(sed '1d' data/TJ.AD.all_causal_v1.txt) \
      <(sed '1d' ipsneurons/ipsc.atac.overlap.summary.txt | cut -f 4) \
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
      >> TJ.AD.all_causal_v1.annotated.txt
