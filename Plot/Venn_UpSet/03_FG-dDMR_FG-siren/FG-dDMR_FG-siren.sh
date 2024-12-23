#!/bin/bash 

# Load modules
module load gcc/5.2.0
module load intervene/0.6.2
module load bedtools/2.30.0

# Generate master regions
cat FG-dDMR.bed FG-siren.bed > master.bed
bedtools sort -i master.bed > master.sort.bed
bedtools merge -d 0 -i master.sort.bed > master.DMR_siren.bed

bedtools intersect -a master.DMR_siren.bed -b FG-dDMR.bed -wa | uniq > FG-dDMR.filter.bed
bedtools intersect -a master.DMR_siren.bed -b FG-siren.bed -wa | uniq > FG-siren.filter.bed

# Generate Venn diagram 
intervene venn -i FG-dDMR.filter.bed FG-siren.filter.bed \
-o intervene_DMR_siren --type genomic --names FG-dDMR,FG-siren \
--colors "#f17272","#ffd966" \
--fontsize 5 --figsize 2 2 --save-overlaps 

# Get bed files for each region
mv 11_FG-dDMR_FG-siren.bed FG-dDMR_FG-siren.bed
mv 10_FG-dDMR.bed FG-dDMR-only.bed
mv 01_FG-siren.bed FG-siren-only.bed

cat FG-dDMR_FG-siren.bed FG-dDMR-only.bed FG-siren-only.bed > all.bed
bedtools sort -i all.bed > all.sort.bed

bedtools shuffle -i all.sort.bed -g OsChr.txt -chrom > Shuffle.bed 
