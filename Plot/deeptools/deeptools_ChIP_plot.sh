#!/bin/bash 

# Load modules
module load deeptools/3.5.2

# Plot
## FG-dDMR & FG-siren
computeMatrix scale-regions \
-R ./bed/FG-dDMR_FG-siren.bed ./bed/FG-dDMR-only.bed ./bed/FG-siren-only.bed ./bed/Shuffle.bed \
-S ./bigWig_ChIP/GFP_[12]_log2.bw ./bigWig_ChIP/MYC_[12]_log2.bw \
-b 2000 -a 2000 --binSize 10 -m 1000 \
--sortRegions descend \
--averageTypeBins mean \
-o matrix_ChIP_FG-dDMR_FG-siren.gz 

### Profile plot
plotProfile -m matrix_ChIP_FG-dDMR_FG-siren.gz \
-o profile_ChIP_FG-dDMR_FG-siren.pdf --plotFileFormat pdf \
--regionsLabel "FG-dDMR_FG-siren" "FG-dDMR-only" "FG-siren-only" "Shuffle" \
--averageType "mean" \
--samplesLabel GFP1 GFP2 MYC1 MYC2 \
--colors "#f8a66c" "#f17272" "#ffd966" "#cccccc" \
--plotHeight 5 --plotWidth 4 \
--numPlotsPerRow 4 \
--refPointLabel "" --yAxisLabel "" \
--legendLocation best 

### Heatmap
plotHeatmap -m matrix_ChIP_FG-dDMR_FG-siren.gz \
-o heatmap_ChIP_FG-dDMR_FG-siren.pdf --plotFileFormat pdf \
--samplesLabel GFP1 GFP2 MYC1 MYC2 \
--regionsLabel "FG-dDMR_FG-siren" "FG-dDMR-only" "FG-siren-only" "Shuffle" \
--startLabel "" --endLabel "" --xAxisLabel "" \
--sortRegions descend --whatToShow 'heatmap and colorbar' \
--legendLocation none \
--heatmapHeight 25 --heatmapWidth 1 --zMin 0 --zMax 0.8 \
--colorMap viridis


## FRE-adjacent regions
computeMatrix reference-point \
-R ./bed/FRE_{0..3}.bed \
-S ./bigWig_ChIP/GFP_[12]_log2.bw ./bigWig_ChIP/MYC_[12]_log2.bw \
-b 2000 -a 2000 --binSize 10 \
--referencePoint center --sortRegions descend \
--averageTypeBins mean \
-o matrix_ChIP_FRE.gz 

### Profile plot
plotProfile -m matrix_ChIP_FRE.gz \
-o profile_ChIP_FRE.pdf --plotFileFormat pdf \
--regionsLabel "FRE_0" "FRE_1" "FRE_2" "FRE_3" \
--averageType "mean" \
--samplesLabel GFP1 GFP2 MYC1 MYC2 \
--plotHeight 5 --plotWidth 4 \
--numPlotsPerRow 4 \
--refPointLabel "" --yAxisLabel "" \
--legendLocation best 

### Heatmap
plotHeatmap -m matrix_ChIP_FRE.gz \
-o heatmap_ChIP_FRE.pdf --plotFileFormat pdf \
--samplesLabel GFP1 GFP2 MYC1 MYC2 \
--regionsLabel "FRE_0" "FRE_1" "FRE_2" "FRE_3" \
--refPointLabel "" --xAxisLabel "" \
--sortRegions descend --whatToShow 'heatmap and colorbar' \
--legendLocation none \
--heatmapHeight 100 --heatmapWidth 1 --zMin 0 --zMax 1 \
--colorMap viridis

#### zooming FRE_0
computeMatrixOperations subset -m matrix_ChIP_FRE.gz \
-o matrix_ChIP_FRE_0.gz --groups "FRE_0.bed"

plotHeatmap -m matrix_ChIP_FRE_0.gz \
-o heatmap_ChIP_FRE_0.pdf --plotFileFormat pdf \
--samplesLabel GFP1 GFP2 MYC1 MYC2 \
--regionsLabel "FRE_0" \
--refPointLabel "" --xAxisLabel "" \
--sortRegions descend --whatToShow 'heatmap and colorbar' \
--legendLocation none \
--heatmapHeight 5 --heatmapWidth 1 --zMin 0 --zMax 1 \
--colorMap viridis
