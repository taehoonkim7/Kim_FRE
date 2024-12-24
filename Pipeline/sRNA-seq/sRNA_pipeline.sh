#!/bin/bash

# Pre-cleaning QC
module load fastqc/0.11.7
fastqc FG2_rep1.fq.gz

# Read cleaning
module load trim_galore/0.6.10
trim_galore --quality 20 --length 10 --small_rna -o . ./FG2_rep1.fq.gz

# Alignment (uniquely mapped reads)
module load bowtie/1.2.3
bowtie -v 0 -q -m 1 -S Osativa FG2_rep1.fq.gz FG2_rep1.uniq.sam

# sam-to-bam conversion
module load samtools/1.18
samtools sort -o FG2_rep1.uniq.bam FG2_rep1.uniq.sam
samtools view -b -F 0x4 FG2_rep1.uniq.bam -o FG2_rep1.map.bam #retain only mapped reads

# Filtering other ncRNAs
module load bedtools/2.30.0
bedtools intersect -v -a FG2_rep1.map.bam \
-b rRNA.gff3 tRNA.gff3 snRNA.gff3 snoRNA.gff3 miR_osa.gff3 > FG2_rep1.siRNA.bam

# Getting 24-nt siRNA reads
samtools view -h FG2_rep1.siRNA.bam | awk 'length($10) == 24 || $1 ~ /^@/' | \
samtools view -bS - > FG2_rep1.24nt.bam #retain only 24-nt siRNA reads

# bigWig conversion
module load deeptools/3.5.2
bamCoverage -b FG2_rep1.24nt.bam -o FG2_rep1.24nt.bigWig --binSize 1 --normalizeUsing CPM
