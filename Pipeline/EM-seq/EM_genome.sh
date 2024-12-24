#!/bin/bash

# Genome preparation
module load bowtie2/2.4.5
module load bismark/0.22.3

bismark_genome_preparation --parallel 8 ./genome
## ./genome/Osativa_323_v7.0.fa