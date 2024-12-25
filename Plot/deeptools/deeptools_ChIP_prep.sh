#!/bin/bash 

# Load modules
module load deeptools/3.5.2

# Data prepartion
## for PCA plot
multiBigwigSummary bins \
-b GSM8121726_CLSY3_GFP_Panicle_Rep1_5bp_RPKM.bigwig GSM8121727_CLSY3_GFP_Panicle_Rep2_5bp_RPKM.bigwig \
GSM8121728_CLSY3_MYC_Panicle_Rep1_5bp_RPKM.bigwig GSM8121729_CLSY3_MYC_Panicle_Rep2_5bp_RPKM.bigwig \
GSM8121730_INPUT_Rep1_5bp_RPKM.bigwig GSM8121731_INPUT_Rep2_5bp_RPKM.bigwig \
-l GFP1 GFP2 MYC1 MYC2 Input1 Input2 -o bw_summary_ChIP.npz 

plotPCA -in bw_summary_ChIP.npz -o null.png --plotFileFormat png \
--rowCenter --log2 --outFileNameData PCA_ChIP_osclsy3.tab --ntop 200 

## for Profile plot & Heatmap
bigwigCompare -b1 GSM8121726_CLSY3_GFP_Panicle_Rep1_5bp_RPKM.bigwig \
-b2 GSM8121730_INPUT_Rep1_5bp_RPKM.bigwig -bs 5 -o GFP_1_log2.bw

bigwigCompare -b1 GSM8121727_CLSY3_GFP_Panicle_Rep2_5bp_RPKM.bigwig \
-b2 GSM8121731_INPUT_Rep2_5bp_RPKM.bigwig -bs 5 -o GFP_2_log2.bw

bigwigCompare -b1 GSM8121728_CLSY3_MYC_Panicle_Rep1_5bp_RPKM.bigwig \
-b2 GSM8121730_INPUT_Rep1_5bp_RPKM.bigwig -bs 5 -o MYC_1_log2.bw

bigwigCompare -b1 GSM8121729_CLSY3_MYC_Panicle_Rep2_5bp_RPKM.bigwig \
-b2 GSM8121731_INPUT_Rep2_5bp_RPKM.bigwig -bs 5 -o MYC_2_log2.bw
