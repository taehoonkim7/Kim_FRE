# A FRE guides de novo DNA methylation during rice gametophyte development

## Introduction

- Codes and datasets for [Preprint](https://doi.org/10.1101/2024.12.03.626448) (to be updated once published)
- FRE: female-specific RdDM-associated element

## Table of Contents

- [Pipeline](#pipeline)
  - [EM-seq pipeline](#em-seq-pipeline)
  - [sRNA-seq pipeline](#srna-seq-pipeline)
  - [mRNA-seq pipeline](#mrna-seq-pipeline)
- [Plot](#plot)
  - [DensityPlot](#densityplot)
  - [Venn Diagram & UpSet Plot](#venn-diagram--upset-plot)
  - [Cumulative Sum Plot](#cumulative-sum-plot)
  - [Metaplot](#metaplot)
  - [Boxplot](#boxplot)
  - [Plots using *deeptools* module](#plots-using-deeptools-module)
  - [Heatmap](#heatmap)
  - [Expression & DNA methylation](#expression--dna-methylation)
  - [Other plots](#other-plots)

## Pipeline

### EM-seq pipeline

- `EM_genome.sh`: prepare genome index for alignment using bismark
- `EM_pipeline.sh`: run EM-seq pipeline
  - Read QC
  - Read cleaning
  - Alingment
  - Deduplication
  - Methylation extraction
  - output: CX_report file for downstream analyses

### sRNA-seq pipeline

- `sRNA_genome.sh`: prepare genome indexing for alignment using bowtie
- `sRNA_pipeline.sh`: run sRNA-seq pipeline
  - Read QC
  - Read cleaning
  - Alingment
  - Filtering other non-coding RNAs
    - rRNA, tRNA, snRNA, snoRNA (based on [RNAcentral](https://rnacentral.org))
    - miRNA (based on [miRBase](https://www.mirbase.org))
  - Getting 24-nt siRNA reads
  - Conversion into bigWig

### mRNA-seq pipeline

- `mRNA_genome.sh`: prepare genome indexing for alignment using hisat2
- `mRNA_pipeline.sh`: run mRNA-seq pipeline
  - Read QC
  - Read cleaning
  - Alingment
  - Read counting
  - TPM conversion using `TPM_calculator.py`
- `TPM_calculator.py`  
  - Usage: `python TPM_calculator.py FG5-1.gene.tsv FG5_rep1.TPM.tsv Osativa_323_v7.0.gene_exons.gff3 Osativa_323_v7.0.GeneID_PrimaryTranscript.txt`
    - `FG5-1.gene.tsv`: read count matrix from `htseq-count`
    - `FG5_rep1.TPM.tsv`: output file
    - `Osativa_323_v7.0.gene_exons.gff3`: GFF3 gene annotation file
    - `Osativa_323_v7.0.GeneID_PrimaryTranscript.txt`: list of primary transcripts of each gene

## Plot

### DensityPlot

> Density plots: **Fig. 1c**, **Extended Data Fig. 2a,b**

- `process_CX_density.py`:
  - Process CX_report into bedGraph files (1-bp or 50-bp windows)
  - Usage: `python process_CX_density.py OsChr.txt FG2 FG2-1_S3_deduplicate_each_cytocine.CX_report.txt FG2-2_S4_deduplicate_each_cytocine.CX_report.txt`
    - `OsChr.txt`: information of rice chromosomes and lengths
    - `FG2`: outfile header
    - `FG2-1_S3_deduplicate_each_cytocine.CX_report.txt FG2-2_S4_deduplicate_each_cytocine.CX_report.txt`: multiple CX_report files to merge
- `DensityPlot.R`:
  - Calculate fractional differences based on filtered 50-bp windows.
  - Generate csv files with fractional differences.
  - Generate density plots.

### Venn Diagram & UpSet Plot

> Venn diagram for FG-dDMRs: **Fig. 1e**

- `FG-dDMR.sh`:
  - Generate Venn diagram of intersections among followings:
    - CHH hypo-DMRs in FG1 compared to other FGs.
    - CHH hypo-DMRs in FG2 compared to other FGs.
    - CHH hyper-DMRs in FG5 compared to other FGs.
    - CHH hyper-DMRs in FG6 compared to other FGs.

> UpSet plot for FG-*siren* loci: **Fig. 2b**

- `FG-siren.sh`:
  - Generate UpSet plot of intersections among FG1- to FG6-*siren* loci.

> Venn diagram for FG-dDMRs and FG-*siren* loci : **Fig. 2c**

- `FG-dDMR_FG-siren.sh`:
  - Generate Venn diagram of intersections between FG-dDMRs and FG-*siren* loci.
  - Get Shuffle regions using `bedtools shuffle`.

### Cumulative Sum Plot

> Cumulative Sum Plot: **Fig. 2a**

- `24nt_cluster.sh`:
  - Make 24-nt siRNA clusters from each sample (≥ 10 reads per million).
  - Generate input files for `CumSumPlot.R`.
    - `FG1.cluster.final.bed`
    - `FG2.cluster.final.bed`
    - `FG3.cluster.final.bed`
    - `FG4.cluster.final.bed`
    - `FG5.cluster.final.bed`
    - `FG6.cluster.final.bed`
    - `FL.cluster.final.bed`
- `CumSumPlot.R`:
  - Generate cumulative sum plot of 24-nt siRNA expression.
    - *x*-axis: 24-nt siRNA clusters arranged by the descending order of 24-nt siRNA expression
    - *y*-axis: Cumulative sum of 24-nt siRNA expression
  - Identify *siren* loci based on Knee point.
  - Required input files:
    - output `.bed` files from `24nt_cluster.sh`
    - `Sample_Description.csv`: description on all input files

### Metaplot

> Metaplot for genome-wide methylation: **Extended Data Fig. 1a**

- `process_CX_meta_genome.py`:
  - Process CX_report into bedGraph files (1-bp or 50-bp windows)
  - Usage: `python process_CX_meta_genome.py FG2-1_S3_deduplicate_each_cytocine.CX_report.txt OsChr.txt FG2-1 4 500000 100000`
    - `FG2-1_S3_deduplicate_each_cytocine.CX_report.txt`: CX_report file to process
    - `OsChr.txt`: information of rice chromosomes and lengths
    - `FG2-1`: outfile header
    - `4`: minimal mapped reads to keep cytosine
    - `500000`: width of the sliding windows (500-kb)
    - `100000`: shift of the sliding windows (100-kb)
- `Metaplot_genome.R`:
  - Generate metaplot of CG, CHG, CHH methylation across all rice chromosomes.
  - Required input files:
    - output `.txt` files from `process_CX_meta_genome.py` (in `./Input/`)
    - `Sample_Description.csv`: description on all input files
    - `Centromere.csv`: the location of centromeric regions of each chromosome

> Metaplot for genes & TEs: **Extended Data Fig. 1b,g**

- `process_CX_meta_gene.py`:
  - Process CX_report for metaplot of 2-kb flanking regions of genes
  - Usage: `python process_CX_meta_gene.py all.locus_brief_info.txt FG2-1_S3_deduplicate_each_cytocine.CX_report.txt FG2-1 2000 100 20`
    - `all.locus_brief_info.txt` — All locus information in the rice genome
    - `FG2-1_S3_deduplicate_each_cytocine.CX_report.txt`: CX_report file to process
    - `FG2-1`: outfile header
    - `2000`: 2-kb upstream/downstream regions
    - `100`: 100-bp bins for upstream and downstream regions
    - `20`: 20 proportional bins for gene bodies
- `process_CX_meta_TE.py`:
  - Process CX_report for metaplot of 2-kb flanking regions of TEs by families
  - Usage: `python process_CX_meta_TE.py all.con_Chr1-12.txt.mod.EDTA.TEanno.split.gff3 FG2-1_S3_deduplicate_each_cytocine.CX_report.txt FG2-1 2000 100 20`
    - `all.con_Chr1-12.txt.mod.EDTA.TEanno.split.gff3` — GFF3 file for TEs in the rice genome (generated using [EDTA](https://github.com/oushujun/EDTA))
    - `FG2-1_S3_deduplicate_each_cytocine.CX_report.txt`: CX_report file to process
    - `FG2-1`: outfile header
    - `2000`: 2-kb upstream/downstream regions
    - `100`: 100-bp bins for upstream and downstream regions
    - `20`: 20 proportional bins for TE bodies
- `Metaplot_GeneTE.R`:
  - Generate metaplot of CG, CHG, CHH methylation of genes and TEs.
  - Required input files:
    - output `.txt` files from `process_CX_meta_gene.py` (in `./Gene/`)
    - output `.txt` files from `process_CX_meta_TE.py` (in `./TE_split/`)
    - `Sample_Description.csv`: description on all input files

### Boxplot

> Boxplot for gene & TE bodies: **Extended Data Fig. 1c,h**

- `process_CX_box_gene.py`:
  - Process CX_report for boxplot of gene bodies
  - Usage: `python process_CX_box_gene.py all.locus_brief_info.txt FG2-1_S3_deduplicate_each_cytocine.CX_report.txt FG2-1`
    - `all.locus_brief_info.txt` — All locus information in the rice genome
    - `FG2-1_S3_deduplicate_each_cytocine.CX_report.txt`: CX_report file to process
    - `FG2-1`: outfile header
- `process_CX_box_TE.py`:
  - Process CX_report for boxplot of TE bodies
  - Usage: `python process_CX_box_TE.py all.con_Chr1-12.txt.mod.EDTA.TEanno.split.gff3 FG2-1_S3_deduplicate_each_cytocine.CX_report.txt FG2-1`
    - `all.con_Chr1-12.txt.mod.EDTA.TEanno.split.gff3` — GFF3 file for TEs in the rice genome (generated using [EDTA](https://github.com/oushujun/EDTA))
    - `FG2-1_S3_deduplicate_each_cytocine.CX_report.txt`: CX_report file to process
    - `FG2-1`: outfile header
- `Boxplot_GeneTE.R`:
  - Generate boxplots of CG, CHG, CHH methylation of gene bodies and TE bodies.
  - Required input files:
    - output `.txt` files from `process_CX_box_gene.py` (in `./Gene/`)
    - output `.txt` files from `process_CX_box_TE.py` (in `./TE_split/`)
    - `Sample_Description.csv`: description on all input files

> Boxplot for regions of interests (DNA methylation)
>
> - FG-dDMRs: **Fig. 1f**
> - FG-dDMRs & FG-*siren*: **Fig. 2g**, **Extended Data Fig. 3c,d**
> - FG-dDMRs & FG-*siren* across various tissues: **Fig. 3a**, **Extended Data Fig. 4a,b**

- `process_CX_box_region.py`:
  - Process CX_report for boxplot of regions of interest
  - Usage: `python process_CX_box_region.py FG-dDMR.bed FG2-1_S3_deduplicate_each_cytocine.CX_report.txt FG2-1_S3.FG-dDMR.txt`
    - `FG-dDMR.bed` — bed file for genomic regions of interest
    - `FG2-1_S3_deduplicate_each_cytocine.CX_report.txt`: CX_report file to process
    - `FG2-1_S3.FG-dDMR.txt`: outfile
- `Boxplot_FG-dDMR.R`:
  - Generate boxplots of CG, CHG, CHH methylation of FG-dDMRs.
  - Required input files:
    - output `.txt` files from `process_CX_box_region.py` (in `./FG-dDMR/`)
    - `Sample_Description_FG-dDMR.csv`: description on all input files
- `Boxplot_FG-dDMR_FG-siren.R`:
  - Generate boxplots of CG, CHG, CHH methylation of four regions of interests (FG-dDMR $\cap$ FG-*siren*, FG-dDMR-only, FG-*siren*-only, and Shuffle).
  - Required input files:
    - output `.txt` files from `process_CX_box_region.py` (in `./FG-dDMR_FG-siren/`)
    - `Sample_Description_FG-dDMR_FG-siren.csv`: description on all input files
- `Boxplot_tissue.R`:
  - Generate boxplots of CG, CHG, CHH methylation of four regions of interests (FG-dDMR $\cap$ FG-*siren*, FG-dDMR-only, FG-*siren*-only, and Shuffle) across different tissues.
  - Perform Kruskal-Wallis test and Dunn's test for multiple comparison
  - Required input files:
    - output `.txt` files from `process_CX_box_region.py` (in `./tissue/`)
    - `Sample_Description_tissue.csv`: description on all input files

> Boxplot for FG-dDMRs & FG-*siren* (24-nt siRNA): **Fig. 2f**

- `Boxplot_24nt.R`:
  - Generate boxplot of 24-nt siRNA expression of four regions of interests (FG-dDMR $\cap$ FG-*siren*, FG-dDMR-only, FG-*siren*-only, and Shuffle) across different tissues.
  - Required input files:
    - `RPM.24nt.csv`: 24-nt siRNA expression levels (normalized by reads per million) for each region.

### Plots using *deeptools* module

### Heatmap

### Expression & DNA methylation

### Other plots
