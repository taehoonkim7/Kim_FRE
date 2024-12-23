#!/bin/bash 

# Load modules
module load gcc/5.2.0
module load intervene/0.6.2
module load bedtools/2.30.0

# Generate master regions
cat FG1_other_CHH_hypo.bed FG2_other_CHH_hypo.bed FG5_other_CHH_hyper.bed FG6_other_CHH_hyper.bed > master.bed
bedtools sort -i master.bed > master.sort.bed
bedtools merge -d 150 -i master.sort.bed > master.DMR.bed #because individual DMRs are merged by 150 bp 

# Filter DMRs by master DMRs
bedtools intersect -a master.DMR.bed -b FG1_other_CHH_hypo.bed -wa > FG1_other_CHH_hypo.master.bed
bedtools intersect -a master.DMR.bed -b FG2_other_CHH_hypo.bed -wa > FG2_other_CHH_hypo.master.bed
bedtools intersect -a master.DMR.bed -b FG5_other_CHH_hyper.bed -wa > FG5_other_CHH_hyper.master.bed
bedtools intersect -a master.DMR.bed -b FG6_other_CHH_hyper.bed -wa > FG6_other_CHH_hyper.master.bed

sort -k1,1n -k2,2n FG1_other_CHH_hypo.master.bed | uniq > FG1_other_CHH_hypo.filter.bed
sort -k1,1n -k2,2n FG2_other_CHH_hypo.master.bed | uniq > FG2_other_CHH_hypo.filter.bed
sort -k1,1n -k2,2n FG5_other_CHH_hyper.master.bed | uniq > FG5_other_CHH_hyper.filter.bed
sort -k1,1n -k2,2n FG6_other_CHH_hyper.master.bed | uniq > FG6_other_CHH_hyper.filter.bed

# Generate Venn diagram 
intervene venn -i FG1_other_CHH_hypo.filter.bed FG2_other_CHH_hypo.filter.bed \
FG5_other_CHH_hyper.filter.bed FG6_other_CHH_hyper.filter.bed \
-o intervene_DMR --type genomic --names FG1_hypo,FG2_hypo,FG5_hyper,FG6_hyper \
--colors "#b28600","#ffd966","#a9d18e","#4e8133" \
--fontsize 7 --figsize 2 2 --save-overlaps

# Get FG-dDMR.bed
cat ./sets/0111_FG2_hypo_FG5_hyper_FG6_hyper.bed ./sets/1011_FG1_hypo_FG5_hyper_FG6_hyper.bed \
./sets/1101_FG1_hypo_FG2_hypo_FG6_hyper.bed ./sets/1110_FG1_hypo_FG2_hypo_FG5_hyper.bed \
./sets/1111_FG1_hypo_FG2_hypo_FG5_hyper_FG6_hyper.bed > FG-dDMR.bed
sort -k1,1n -k2,2n FG-dDMR.bed > FG-dDMR.sort.bed
mv FG-dDMR.sort.bed FG-dDMR.bed
