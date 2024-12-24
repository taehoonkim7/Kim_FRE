#!/bin/bash

# FG2-1_S3_L004_R1_001.fastq.gz - read1 from paired-end sequencing
# FG2-1_S3_L004_R2_001.fastq.gz - read2 from paired-end sequencing

# Pre-cleaning QC
module load fastqc/0.11.7

fastqc FG2-1_S3_L004_R1_001.fastq.gz
fastqc FG2-1_S3_L004_R2_001.fastq.gz


# Read cleaning
module load trimmomatic/0.39

trimmomatic PE -phred33 FG2-1_S3_L004_R1_001.fastq.gz FG2-1_S3_L004_R2_001.fastq.gz \
$FG2-1_S3_R1_1P.fastq.gz FG2-1_S3_R1_1U.fastq.gz FG2-1_S3_R2_2P.fastq.gz FG2-1_S3_R2_2U.fastq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 


# Alignment
module load bowtie2/2.4.5
module load bismark/0.22.3

bismark ./genome -I 50 -N 1 --multicore 8 -q \
-1 FG2-1_S3_R1_1P.fastq.gz -2 FG2-1_S3_R2_2P.fastq.gz


# Deduplication 
deduplicate_bismark --bam -p \
-o FG2-1_S3_R1_1P_bismark_bt2_pe.bam FG2-1_S3_R1_1P_bismark_bt2_pe.bam


# Methylation extraction
bismark_methylation_extractor -p --multicore 8 FG2-1_S3_R1_1P_bismark_bt2_pe.deduplicated.bam

bismark2bedGraph --CX -o FG2-1_S3 \
CHG_OB_FG2-1_S3_1P_bismark_bt2_pe.deduplicated.txt CHG_OT_FG2-1_S3_1P_bismark_bt2_pe.deduplicated.txt \
CHH_OB_FG2-1_S3_1P_bismark_bt2_pe.deduplicated.txt CHH_OT_FG2-1_S3_1P_bismark_bt2_pe.deduplicated.txt \
CpG_OB_FG2-1_S3_1P_bismark_bt2_pe.deduplicated.txt CpG_OT_FG2-1_S3_1P_bismark_bt2_pe.deduplicated.txt

coverage2cytosine --CX --genome_folder ./genome -o FG2-1_S3_deduplicate_each_cytocine FG2-1_S3.gz.bismark.cov.gz
## output: FG2-1_S3_deduplicate_each_cytocine.CX_report.txt