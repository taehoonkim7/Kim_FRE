#!/bin/bash 

# Load modules
module load deeptools/3.5.2

# plotHeatmap
## mCG
computeMatrix reference-point \
-R ./bed/FG-dDMR_FG-siren.bed ./bed/FG-dDMR-only.bed ./bed/FG-siren-only.bed \
-S ./bigWig_CpG/FG{1..6}_CpG.bw ./bigWig_CpG/FL_CpG.bw \
-b 2000 -a 2000 --binSize 200 \
--referencePoint center --sortRegions keep \
-o matrix_CpG_FG-dDMR_FG-siren.gz 

plotHeatmap -m matrix_CpG_FG-dDMR_FG-siren.gz \
-o plotHeatmap_CpG_FG-dDMR_FG-siren.pdf --plotFileFormat pdf \
--samplesLabel FG{1..6} FL \
--regionsLabel "FG-dDMR_FG-siren" "FG-dDMR-only" "FG-siren-only" \
--refPointLabel "" --xAxisLabel "" \
--sortRegions descend --whatToShow 'heatmap and colorbar' \
--legendLocation none \
--heatmapHeight 7 --heatmapWidth 1 --zMin 0 --yMax 100 \
--missingDataColor '#FFFFE0' \
--colorList '#FFFFE0,#8B0000' #lightyellow,red4

## mCHG
computeMatrix reference-point \
-R ./bed/FG-dDMR_FG-siren.bed ./bed/FG-dDMR-only.bed ./bed/FG-siren-only.bed \
-S ./bigWig_CHG/FG{1..6}_CHG.bw ./bigWig_CHG/FL_CHG.bw \
-b 2000 -a 2000 --binSize 200 \
--referencePoint center --sortRegions keep \
-o matrix_CHG_FG-dDMR_FG-siren.gz 

plotHeatmap -m matrix_CHG_FG-dDMR_FG-siren.gz \
-o plotHeatmap_CHG_FG-dDMR_FG-siren.pdf --plotFileFormat pdf \
--samplesLabel FG{1..6} FL \
--regionsLabel "FG-dDMR_FG-siren" "FG-dDMR-only" "FG-siren-only" \
--refPointLabel "" --xAxisLabel "" \
--sortRegions descend --whatToShow 'heatmap and colorbar' \
--legendLocation none \
--heatmapHeight 7 --heatmapWidth 1 --zMin 0 --yMax 100 \
--missingDataColor '#FFFFE0' \
--colorList '#FFFFE0,#BC4d00' #lightyellow,reddish-brown

## mCHH
computeMatrix reference-point \
-R ./bed/FG-dDMR_FG-siren.bed ./bed/FG-dDMR-only.bed ./bed/FG-siren-only.bed \
-S ./bigWig_CHH/FG{1..6}_CHH.bw ./bigWig_CHH/FL_CHH.bw \
-b 2000 -a 2000 --binSize 200 \
--referencePoint center --sortRegions keep \
-o matrix_CHH_FG-dDMR_FG-siren.gz 

plotHeatmap -m matrix_CHH_FG-dDMR_FG-siren.gz \
-o plotHeatmap_CHH_FG-dDMR_FG-siren.pdf --plotFileFormat pdf \
--samplesLabel FG{1..6} FL \
--regionsLabel "FG-dDMR_FG-siren" "FG-dDMR-only" "FG-siren-only" \
--refPointLabel "" --xAxisLabel "" \
--sortRegions descend --whatToShow 'heatmap and colorbar' \
--legendLocation none \
--heatmapHeight 7 --heatmapWidth 1 --zMin 0 --yMax 100 \
--missingDataColor '#FFFFE0' \
--colorList '#FFFFE0,#EE9A00' #lightyellow,orange2

## 24-nt siRNA
computeMatrix reference-point \
-R ./bed/FG-dDMR_FG-siren.bed ./bed/FG-dDMR-only.bed ./bed/FG-siren-only.bed \
-S ./bigWig_24nt/FG{1..6}.24nt.bw ./bigWig_24nt/FL.24nt.bw \
-b 2000 -a 2000 --binSize 200 \
--referencePoint center --sortRegions keep \
-o matrix_24nt_FG-dDMR_FG-siren.gz 

plotHeatmap -m matrix_24nt_FG-dDMR_FG-siren.gz \
-o plotHeatmap_24nt_FG-dDMR_FG-siren.pdf --plotFileFormat pdf \
--samplesLabel FG{1..6} FL \
--regionsLabel "FG-dDMR_FG-siren" "FG-dDMR-only" "FG-siren-only" \
--refPointLabel "" --xAxisLabel "" \
--sortRegions descend --whatToShow 'heatmap and colorbar' \
--legendLocation none \
--heatmapHeight 7 --heatmapWidth 1 \
--missingDataColor '#FFFFE0' \
--colorList '#FFFFE0,#0000ff' #lightyellow,blue
