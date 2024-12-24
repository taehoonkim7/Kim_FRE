#!/bin/bash

# Genome preparation
module load gffread/0.12.7
gffread Osativa_323_v7.0.gene_exons.gff3 -T -o Osativa_323_v7.0.gene_exons.gtf

module load hisat2/2.2.1
hisat2_extract_splice_sites.py Osativa_323_v7.0.gene_exons.gtf > splicesites.tsv
hisat2_extract_exons.py Osativa_323_v7.0.gene_exons.gtf > exons.tsv 
hisat2-build -p 8 --ss splicesites.tsv --exon exons.tsv Osativa_323_v7.0.fa Rice_index
