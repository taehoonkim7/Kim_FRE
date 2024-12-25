#!/bin/bash 

# Load modules
module load deeptools/3.5.2

# Plot
## 24nt-siRNA
computeMatrix reference-point \
-R ./bed/FRE_{0..3}.bed \
-S ./bigWig_24nt/Seedling.24nt.bw ./bigWig_24nt/FL.24nt.bw \
./bigWig_24nt/Panicle.24nt.bw ./bigWig_24nt/Stamen.24nt.bw \
./bigWig_24nt/FG[16].24nt.bw ./bigWig_24nt/Pistil_Nip.24nt.bw \
./bigWig_24nt/Embryo_Nip.24nt.bw ./bigWig_24nt/Endosperm_Nip.24nt.bw \
-b 2000 -a 2000 --binSize 10 \
--referencePoint center --sortRegions descend \
--averageTypeBins mean --missingDataAsZero \
-o matrix_24nt_FRE.gz 

### for FRE_0
plotProfile -m matrix_24nt_FRE.gz --plotType heatmap \
-o profileHeatmap_24nt_FRE_0.pdf --plotFileFormat pdf \
--regionsLabel "FRE_0" "FRE_1" "FRE_2" "FRE_3" \
--averageType "mean" \
--samplesLabel Seedling FL Panicle Stamen FG1 FG6 Pistil Embryo Endosperm \
--plotHeight 5 --plotWidth 8 \
--numPlotsPerRow 1 \
--yMin 0 --yMax 100 \
--refPointLabel "" --yAxisLabel "" \
--legendLocation best --perGroup \
--colors Blues

### for FRE_1
plotProfile -m matrix_24nt_FRE.gz --plotType heatmap \
-o profileHeatmap_24nt_FRE_1.pdf --plotFileFormat pdf \
--regionsLabel "FRE_0" "FRE_1" "FRE_2" "FRE_3" \
--averageType "mean" \
--samplesLabel Seedling FL Panicle Stamen FG1 FG6 Pistil Embryo Endosperm \
--plotHeight 5 --plotWidth 8 \
--numPlotsPerRow 1 \
--yMin 0 --yMax 4 \
--refPointLabel "" --yAxisLabel "" \
--legendLocation best --perGroup \
--colors Blues

### for FRE_2/3
plotProfile -m matrix_24nt_FRE.gz --plotType heatmap \
-o profileHeatmap_24nt_FRE_2.pdf --plotFileFormat pdf \
--regionsLabel "FRE_0" "FRE_1" "FRE_2" "FRE_3" \
--averageType "mean" \
--samplesLabel Seedling FL Panicle Stamen FG1 FG6 Pistil Embryo Endosperm \
--plotHeight 5 --plotWidth 8 \
--numPlotsPerRow 1 \
--yMin 0 --yMax 0.4 \
--refPointLabel "" --yAxisLabel "" \
--legendLocation best --perGroup \
--colors Blues

## mCG
computeMatrix reference-point \
-R ./bed/FRE_{0..3}.bed \
-S ./bigWig_CpG/Seedling_CG.bw ./bigWig_CpG/FL_CpG.bw \
./bigWig_CpG/Panicle_CG.bw ./bigWig_CpG/Stamen_CG.bw \
./bigWig_CpG/FG[16]_CpG.bw ./bigWig_CpG/Pistil_CG.bw \
./bigWig_CpG/Embryo_CG.bw ./bigWig_CpG/Endosperm_CG.bw \
-b 2000 -a 2000 --binSize 10 \
--referencePoint center --sortRegions keep \
--averageTypeBins mean --skipZeros \
-o matrix_CpG_FRE.gz 

plotProfile -m matrix_CpG_FRE.gz \
--plotType heatmap \
-o profileHeatmap_CpG_FRE.pdf --plotFileFormat pdf \
--regionsLabel "FRE_0" "FRE_1" "FRE_2" "FRE_3" \
--averageType "mean" \
--samplesLabel Seedling FL Panicle Stamen FG1 FG6 Pistil Embryo Endosperm \
--plotHeight 5 --plotWidth 8 \
--numPlotsPerRow 1 \
--yMin 0 \
--refPointLabel "" --yAxisLabel "" \
--legendLocation best --perGroup \
--colors Reds

## mCHG
computeMatrix reference-point \
-R ./bed/FRE_{0..3}.bed \
-S ./bigWig_CHG/Seedling_CHG.bw ./bigWig_CHG/FL_CHG.bw \
./bigWig_CHG/Panicle_CHG.bw ./bigWig_CHG/Stamen_CHG.bw \
./bigWig_CHG/FG[16]_CHG.bw ./bigWig_CHG/Pistil_CHG.bw \
./bigWig_CHG/Embryo_CHG.bw ./bigWig_CHG/Endosperm_CHG.bw \
-b 2000 -a 2000 --binSize 10 \
--referencePoint center --sortRegions keep \
--averageTypeBins mean --skipZeros \
-o matrix_CHG_FRE.gz 

plotProfile -m matrix_CHG_FRE.gz \
--plotType heatmap \
-o profileHeatmap_CHG_FRE.pdf --plotFileFormat pdf \
--regionsLabel "FRE_0" "FRE_1" "FRE_2" "FRE_3" \
--averageType "mean" \
--samplesLabel Seedling FL Panicle Stamen FG1 FG6 Pistil Embryo Endosperm \
--plotHeight 5 --plotWidth 8 \
--numPlotsPerRow 1 \
--yMin 0 \
--refPointLabel "" --yAxisLabel "" \
--legendLocation best --perGroup \
--colors YlOrBr

## mCHH
computeMatrix reference-point \
-R ./bed/FRE_{0..3}.bed \
-S ./bigWig_CHH/Seedling_CHH.bw ./bigWig_CHH/FL_CHH.bw \
./bigWig_CHH/Panicle_CHH.bw ./bigWig_CHH/Stamen_CHH.bw \
./bigWig_CHH/FG[16]_CHH.bw ./bigWig_CHH/Pistil_CHH.bw \
./bigWig_CHH/Embryo_CHH.bw ./bigWig_CHH/Endosperm_CHH.bw \
-b 2000 -a 2000 --binSize 10 \
--referencePoint center --sortRegions keep \
--averageTypeBins mean --skipZeros \
-o matrix_CHH_FRE.gz 

plotProfile -m matrix_CHH_FRE.gz \
--plotType heatmap \
-o profileHeatmap_CHH_FRE.pdf --plotFileFormat pdf \
--regionsLabel "FRE_0" "FRE_1" "FRE_2" "FRE_3" \
--averageType "mean" \
--samplesLabel Seedling FL Panicle Stamen FG1 FG6 Pistil Embryo Endosperm \
--plotHeight 5 --plotWidth 8 \
--numPlotsPerRow 1 \
--yMin 0 --yMax 60 \
--refPointLabel "" --yAxisLabel "" \
--legendLocation best --perGroup \
--colors Oranges
