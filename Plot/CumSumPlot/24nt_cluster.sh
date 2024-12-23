#!/bin/bash 

# Load modules
module load samtools/1.18
module load bedtools/2.30.0

# Merge replicates of 24-nt siRNA bam files 
samtools merge FG1.24nt.bam FG1_1.uniq.24nt_siRNA.bam FG1_2.uniq.24nt_siRNA.bam #FG1
mv FG2_1.uniq.24nt_siRNA.bam FG2.24nt.bam #FG2
mv FG3_1.uniq.24nt_siRNA.bam FG3.24nt.bam #FG3
samtools merge FG4.24nt.bam FG4_1.uniq.24nt_siRNA.bam FG4_2.uniq.24nt_siRNA.bam #FG4
samtools merge FG5.24nt.bam FG5_1.uniq.24nt_siRNA.bam FG5_2.uniq.24nt_siRNA.bam #FG5
samtools merge FG6.24nt.bam FG6_1.uniq.24nt_siRNA.bam FG6_2.uniq.24nt_siRNA.bam #FG6
samtools merge FL.24nt.bam FL_1.uniq.24nt_siRNA.bam FL_2.uniq.24nt_siRNA.bam #FL

# Make 24-nt siRNA bed files
for sample in FG{1..6} FL
do
bedtools bamtobed -i ${sample}.24nt.bam > ${sample}.24nt.bed #bam-to-bed
done

# Make master cluster
cat FG{1..6}.24nt.bed FL.24nt.bed > master.bed
bedtools sort -i master.bed > master.sort.bed
bedtools merge -d 75 -i master.sort.bed > master.cluster.bed #merge clusters within 75-bp

grep -v -e "ChrUn" -e "ChrSy" master.cluster.bed > temp.bed
mv temp.bed > master.cluster.bed

# Count 24-nt reads by master cluster
for sample in FG{1..6} FL
do
bedtools coverage -a master.cluster.bed -b ${sample}.24nt.bed -counts > ${sample}.cluster.count.bed
done

# Filter 24-nt siRNA clutsers (RPM >= 10)
awk '$4 >= 29' FG1.cluster.count.bed > FG1.cluster.filtered.bed
awk '$4 >= 24' FG2.cluster.count.bed > FG2.cluster.filtered.bed
awk '$4 >= 39' FG3.cluster.count.bed > FG3.cluster.filtered.bed
awk '$4 >= 68' FG4.cluster.count.bed > FG4.cluster.filtered.bed
awk '$4 >= 68' FG5.cluster.count.bed > FG5.cluster.filtered.bed
awk '$4 >= 59' FG6.cluster.count.bed > FG6.cluster.filtered.bed
awk '$4 >= 11' FL.cluster.count.bed > FL.cluster.filtered.bed

# Make master cluster with expression in at least one sample
cat FG{1..6}.cluster.filtered.bed FL.cluster.filtered.bed > master.cluster.filtered.bed
bedtools sort -i master.cluster.filtered.bed > master.cluster.filtered.sort.bed
bedtools merge -d 0 -i master.cluster.filtered.sort.bed > master.cluster.filtered.bed

# Make final cluster 
for sample in FG{1..6} FL
do
bedtools map -a master.cluster.filtered.bed -b ${sample}.cluster.filtered.bed \
-c 4 -o collapse -null 0 > ${sample}.cluster.final.bed
done
