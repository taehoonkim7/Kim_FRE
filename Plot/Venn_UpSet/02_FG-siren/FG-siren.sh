#!/bin/bash 

# Load modules
module load gcc/5.2.0
module load intervene/0.6.2

#6-way upset
intervene upset -i FG{1..6}.cumSum.bed \
-o intervene_siren --type genomic --names FG1,FG2,FG3,FG4,FG5,FG6 \
--figsize 6 3 --ninter 50 --save-overlaps

# Get FG-siren.bed
cp ./sets/111111_FG1_FG2_FG3_FG4_FG5_FG6.bed FG-siren.bed
sort -k1,1n -k2,2n FG-siren.bed > FG-siren.sort.bed
mv FG-siren.sort.bed FG-siren.bed
