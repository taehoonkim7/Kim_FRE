#!/bin/bash 

# Load modules
module load deeptools/3.5.2

# plotHeatmap
## mCG
computeMatrix reference-point \
-R ./bed/FG-dDMR_FG-siren.bed ./bed/FG-dDMR-only.bed \
-S ./bigWig_CpG/Seedling_CG.bw ./bigWig_CpG/FL_CpG.bw ./bigWig_CpG/Panicle_CG.bw \
./bigWig_CpG/Stamen_CG.bw ./bigWig_CpG/FG[16]_CpG.bw ./bigWig_CpG/Pistil_CG.bw \
./bigWig_CpG/Pistil_2-3_CG.bw ./bigWig_CpG/Pistil_2-6_CG.bw ./bigWig_CpG/EC_CG.bw \
./bigWig_CpG/CC_CG.bw ./bigWig_CpG/Embryo_CG.bw ./bigWig_CpG/Endosperm_CG.bw \
-b 2000 -a 2000 --binSize 200 \
--referencePoint center --sortRegions keep \
--averageTypeBins mean --skipZeros \
-o matrix_CpG_tissue.gz 

plotHeatmap -m matrix_CpG_tissue.gz \
-o plotHeatmap_CpG_tissue.pdf --plotFileFormat pdf \
--samplesLabel Seedling FL Panicle Stamen FG1 FG6 \
Pistil rdr2-3 rdr2-6 EC CC Embryo Endosperm \
--regionsLabel "FG-dDMR_FG-siren" "FG-dDMR-only" \
--refPointLabel "" --xAxisLabel "" \
--sortRegions descend --whatToShow 'heatmap and colorbar' \
--legendLocation none \
--heatmapHeight 6 --heatmapWidth 1 --zMin 0 --yMax 100 \
--missingDataColor '#FFFFE0' \
--colorList '#FFFFE0,#8B0000' #lightyellow,red4

## mCHG
computeMatrix reference-point \
-R ./bed/FG-dDMR_FG-siren.bed ./bed/FG-dDMR-only.bed \
-S ./bigWig_CHG/Seedling_CHG.bw ./bigWig_CHG/FL_CHG.bw ./bigWig_CHG/Panicle_CHG.bw \
./bigWig_CHG/Stamen_CHG.bw ./bigWig_CHG/FG[16]_CHG.bw ./bigWig_CHG/Pistil_CHG.bw \
./bigWig_CHG/Pistil_2-3_CHG.bw ./bigWig_CHG/Pistil_2-6_CHG.bw ./bigWig_CHG/EC_CHG.bw \
./bigWig_CHG/CC_CHG.bw ./bigWig_CHG/Embryo_CHG.bw ./bigWig_CHG/Endosperm_CHG.bw \
-b 2000 -a 2000 --binSize 200 \
--referencePoint center --sortRegions keep \
--averageTypeBins mean --skipZeros \
-o matrix_CHG_tissue.gz 

plotHeatmap -m matrix_CHG_tissue.gz \
-o plotHeatmap_CHG_tissue.pdf --plotFileFormat pdf \
--samplesLabel Seedling FL Panicle Stamen FG1 FG6 \
Pistil rdr2-3 rdr2-6 EC CC Embryo Endosperm \
--regionsLabel "FG-dDMR_FG-siren" "FG-dDMR-only" \
--refPointLabel "" --xAxisLabel "" \
--sortRegions descend --whatToShow 'heatmap and colorbar' \
--legendLocation none \
--heatmapHeight 6 --heatmapWidth 1 --zMin 0 --yMax 100 \
--missingDataColor '#FFFFE0' \
--colorList '#FFFFE0,#BC4d00' #lightyellow,reddish-brown

## mCHH
computeMatrix reference-point \
-R ./bed/FG-dDMR_FG-siren.bed ./bed/FG-dDMR-only.bed \
-S ./bigWig_CHH/Seedling_CHH.bw ./bigWig_CHH/FL_CHH.bw ./bigWig_CHH/Panicle_CHH.bw \
./bigWig_CHH/Stamen_CHH.bw ./bigWig_CHH/FG[16]_CHH.bw ./bigWig_CHH/Pistil_CHH.bw \
./bigWig_CHH/Pistil_2-3_CHH.bw ./bigWig_CHH/Pistil_2-6_CHH.bw ./bigWig_CHH/EC_CHH.bw \
./bigWig_CHH/CC_CHH.bw ./bigWig_CHH/Embryo_CHH.bw ./bigWig_CHH/Endosperm_CHH.bw \
-b 2000 -a 2000 --binSize 200 \
--referencePoint center --sortRegions keep \
--averageTypeBins mean --skipZeros \
-o matrix_CHH_tissue.gz 

plotHeatmap -m matrix_CHH_tissue.gz \
-o plotHeatmap_CHH_tissue.pdf --plotFileFormat pdf \
--samplesLabel Seedling FL Panicle Stamen FG1 FG6 \
Pistil rdr2-3 rdr2-6 EC CC Embryo Endosperm \
--regionsLabel "FG-dDMR_FG-siren" "FG-dDMR-only" \
--refPointLabel "" --xAxisLabel "" \
--sortRegions descend --whatToShow 'heatmap and colorbar' \
--legendLocation none \
--heatmapHeight 6 --heatmapWidth 1 --zMin 0 --yMax 100 \
--missingDataColor '#FFFFE0' \
--colorList '#FFFFE0,#EE9A00' #lightyellow,orange2
