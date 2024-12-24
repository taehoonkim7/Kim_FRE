#!/bin/bash

# FG5-1_S13_L002_R1_001.fastq.gz - read1 from paired-end sequencing
# FG5-1_S13_L002_R2_001.fastq.gz - read2 from paired-end sequencing

# Pre-cleaning QC
module load fastqc/0.11.7

fastqc FG5-1_S13_L002_R1_001.fastq.gz
fastqc FG5-1_S13_L002_R2_001.fastq.gz


# Read cleaning
module load trimmomatic/0.39

trimmomatic PE -phred33 FG5-1_S13_L002_R1_001.fastq.gz FG5-1_S13_L002_R2_001.fastq.gz \
FG5-1_S13_R1_1P.fastq.gz FG5-1_S13_R1_1U.fastq.gz FG5-1_S13_R2_2P.fastq.gz FG5-1_S13_R2_2U.fastq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 

# Alignment
module load hisat2/2.2.1

hisat2 -p 8 --rg-id=${base} --rg PL:ILLUMINA -x Rice_index \
--dta --rna-strandness RF -1 FG5-1_S13_R1_1P.fastq.gz -2 FG5-1_S13_R2_2P.fastq.gz \
-S FG5-1_S13.sam

# Sam-to-bam conversion
module load samtools/1.18
samtools sort -@ 8 -o FG5-1_S13.bam FG5-1_S13.sam

# Read counting
module load htseq/0.11.2
htseq-count --format bam --order pos --mode union --stranded reverse --minaqual 1 \
--type exon --idattr gene_id FG5-1_S13.bam Osativa_323_v7.0.gene_exons.gtf \
> FG5-1.gene.tsv

# TPM conversion
python3 TPM_calculator.py FG5-1.gene.tsv FG5_rep1.TPM.tsv \
Osativa_323_v7.0.gene_exons.gff3 Osativa_323_v7.0.GeneID_PrimaryTranscript.txt
