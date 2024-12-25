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

> Boxplot for regions of interest (DNA methylation)
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
  - Generate boxplots of CG, CHG, CHH methylation of four regions of interest (FG-dDMR $\cap$ FG-*siren*, FG-dDMR-only, FG-*siren*-only, and Shuffle).
  - Required input files:
    - output `.txt` files from `process_CX_box_region.py` (in `./FG-dDMR_FG-siren/`)
    - `Sample_Description_FG-dDMR_FG-siren.csv`: description on all input files
- `Boxplot_tissue.R`:
  - Generate boxplots of CG, CHG, CHH methylation of four regions of interest (FG-dDMR $\cap$ FG-*siren*, FG-dDMR-only, FG-*siren*-only, and Shuffle) across different tissues.
  - Perform Kruskal-Wallis test and Dunn's test for multiple comparison
  - Required input files:
    - output `.txt` files from `process_CX_box_region.py` (in `./tissue/`)
    - `Sample_Description_tissue.csv`: description on all input files

> Boxplot for FG-dDMRs & FG-*siren* (24-nt siRNA): **Fig. 2f**

- `Boxplot_24nt.R`:
  - Generate boxplot of 24-nt siRNA expression of four regions of interest (FG-dDMR $\cap$ FG-*siren*, FG-dDMR-only, FG-*siren*-only, and Shuffle) across different tissues.
  - Required input files:
    - `RPM.24nt.csv`: 24-nt siRNA expression levels (normalized by reads per million) for each region.

### Plots using *deeptools* module

> Heatmap for regions of interests (24-nt siRNA and DNA methylation)
>
> - FG-dDMRs & FG-*siren*: **Fig. 2d,e**, **Extended Data Fig. 3a,b**
> - FG-dDMRs & FG-*siren* across various tissues: **Fig. 3b**, **Extended Data Fig. 4c,d**

- `deeptools_FG-dDMR_FG-siren.sh`:
  - Generate heatmap for 24-nt siRNA expression and DNA methylation levels of regions of interest (FG-dDMR $\cap$ FG-*siren*, FG-dDMR-only, and FG-*siren*-only).
  - Required input files:
    - bigWig files for 24-nt siRNA and DNA methylation (CG, CHG, and CHH)
    - bed files for regions of interest (FG-dDMR $\cap$ FG-*siren*, FG-dDMR-only, and FG-*siren*-only)
- `deeptools_tissue.sh`:
  - Generate heatmap for DNA methylation levels of regions of interest (FG-dDMR $\cap$ FG-*siren* and FG-dDMR-only) across different tissues.
  - Required input files:
    - bigWig files for DNA methylation (CG, CHG, and CHH)
    - bed files for regions of interest (FG-dDMR $\cap$ FG-*siren* and FG-dDMR-only)

> Heatmap for FRE-adjacent regions (24-nt siRNA and DNA methylation)
>
> - 24-nt siRNA expression: **Fig. 4d**, **Extended Data Fig. 6a**
> - CG, CHG, and CHH methylation: **Extended Data Fig. 6b–d**

- `deeptools_FRE.sh`:
  - Generate heatmap for 24-nt siRNA expression and DNA methylation levels of FRE-adjacent regions across different tissues.
  - Required input files:
    - bigWig files for 24-nt siRNA and DNA methylation (CG, CHG, and CHH)
    - bed files for FRE in the rice genome categorized by the number of mismatches

> Plots using OsCLSY3 ChIP-seq data ([Pal et al., 2024](https://doi.org/10.1038/s41467-024-52239-z))
>
> - PCA plot for different ChIP-seq samples: **Extended Data Fig. 7a**
> - Profile plots for FG-dDMRs & FG-*siren*: **Fig. 5b**, **Extended Data Fig. 7b**
> - Profile plots for FRE-adjacent regions: **Fig. 5c**, **Extended Data Fig. 7c**
> - Heatmaps for FG-dDMRs & FG-*siren*: **Extended Data Fig. 7d**
> - Heatmaps for FRE-adjacent regions: **Extended Data Fig. 7e**

- `deeptools_ChIP_prep.sh`:
  - Prepare input files for `PCA_ChIP_osclsy3.R`
  - Prepare input files for `deeptools_ChIP_plot.sh`
  - Required input files:
    - bigWig files for ChIP-seq samples (OsCLSY3-GFP, OsCLSY3-2xMYC, and Input)
- `PCA_ChIP_osclsy3.R`:
  - Generate PCA plot for different ChIP-seq samples.
  - Required input files:
    - `.tab` file generated by `deeptools_ChIP_prep.sh`
- `deeptools_ChIP_plot.sh`:
  - Generate Profile plots and heatmaps for ChIP-seq signals compared to input.
  - Required input files:
    - processed bigWig files using `deeptools_ChIP_prep.sh`
    - bed files for regions of interest (FG-dDMR $\cap$ FG-*siren* and FG-dDMR-only)
    - bed files for FRE in the rice genome categorized by the number of mismatches

### Expression & DNA methylation

> Boxplots for TPM values for genes categorized by expression levels: **Extended Data Fig. 1d**

- `Boxplot_TPM.R`:
  - Generate boxplot of TPM values of genes categorized based on expression levels.
  - Generate `.txt` files containing genes and their expression categories for each sample.
  - Generate `gene_group_exp.csv` containing genes and their expression categories for all samples.
  - Required input files:
    - `FG.TPM.txt`: TPM values of all genes in all samples analyzed

> Mataplots for DNA methylation levels for genes categorized by expression levels: **Extended Data Fig. 1e**

- `process_CX_meta_exp.py`:
  - Process CX_report for metaplot of 2-kb flanking regions of genes by expression categories.
  - Usage: `python process_CX_meta_exp.py all.locus_brief_info.txt FG2_exp.txt FG2-1_S3_deduplicate_each_cytocine.CX_report.txt FG-1 2000 100 20`
    - `all.locus_brief_info.txt`: All locus information in the rice genome
    - `FG2_exp.txt`: genes and their categories based on expression levels (generated using `Boxplot_TPM.R`)
    - `FG2-1_S3_deduplicate_each_cytocine.CX_report.txt`: CX_report file to process
    - `FG2-1`: outfile header
    - `2000`: 2-kb upstream/downstream regions
    - `100`: 100-bp bins for upstream and downstream regions
    - `20`: 20 proportional bins for gene bodies
- `Metaplot_gene_exp.R`:
  - Generate metaplot of CG, CHG, CHH methylation of genes by expression categories.
  - Required input files:
    - output `.txt` files from `process_CX_meta_exp.py` (in `./Input/`)
    - `Sample_Description_meta.csv`: description on all input files

> Boxplots for DNA methylation levels for genes categorized by expression levels: **Extended Data Fig. 1f**

- `Boxplot_gene_exp.R`:
  - Generate boxplots of CG, CHG, CHH methylation of gene bodies by expression categories.
  - Required input files:
    - output `.txt` files from `process_CX_box_gene.py` (in `./Gene/`)
      - Same as input files used in `Boxplot_GeneTE.R`
    - `gene_group_exp.csv`: output file from `Boxplot_TPM.R`
    - `Sample_Description.csv`: description on all input files

### Other plots

> Barplots for DMR analysis in CG, CHG, and CHH contexts: **Fig. 1d**, **Extended Data Fig. 2d,e**

- `DMR.R`:
  - Generate barplots of the number of DMRs by sequence contexts (CG, CHG, and CHH).
  - Required input files:
    - `DMR.csv`: the number of DMRs

> Genomic elements overlapping with FG-dDMRs
>
> - Barplot for genic elements and TEs overlapping with FG-dDMRs: **Fig. 1g**
> - Barplot for TE families overlapping with FG-dDMRs: **Fig. 1h**

- `Element_FG-dDMR.R`:
  - Generate barplots of the number of genic elements and TEs overlapping with FG-dDMRs.
  - Required input files:
    - `Element_FG-dDMR.csv`: the number of genic elements and TEs
- `TE_FG-dDMR.R`:
  - Generate barplots of the ratio of observed and expected counts of TEs overlapping with FG-dDMRs by families.
  - Generate table of Fisher's exact test results.
  - Required input files:
    - `TE_FG-dDMR.csv`: the number of observed TEs on FG-dDMRs and total TEs in the genome by families

> DNA motif enrichment
>
> - Heatmap for motif enrichment in FG-dDMRs & FG-*siren*: **Extended Data Fig. 5b,c**
> - Heatmap for motif enrichment in TE families: **Extended Data Fig. 5d**
> - Barplot for FRE enrichment in FG-dDMRs & FG-*siren*: **Fig. 4b**
> - Heatmap for 24-nt siRNA expression and FRE occurrence in FG-dDMRs & FG-*siren*: **Fig. 4c**

- `heatmap_motif_region.R`:
  - Generate heatmap for motif enrichment in regions of interest (FG-dDMR $\cap$ FG-*siren*, FG-dDMR-only, FG-*siren*-only, and Shuffle).
  - Generate heatmap for *p* values from Fisher's exact test to compare pairwise motif enrichment in regions of interest (FG-dDMR $\cap$ FG-*siren*, FG-dDMR-only, FG-*siren*-only, and Shuffle).
  - Required input files:
    - `motif_enrich_region.csv`: the number of motif occurrence in regions of interest
- `heatmap_motif_TE.R`:
  - Generate heatmap for motif enrichment in TE families.
  - Required input files:
    - `motif_enrich_TE.csv`: the number of motif occurrence in TE families
- `barplot_FRE_region.R`:
  - Generate barplot for FRE enrichment in regions of interest (FG-dDMR $\cap$ FG-*siren*, FG-dDMR-only, FG-*siren*-only, and Shuffle).
  - Required input files:
    - `motif_enrich_region.csv`: the number of motif occurrence in regions of interest (the same input file as used in `heatmap_motif_region.R`)
- `heatmap_24nt_FRE_region.R`:
  - Generate heatmap for 24-nt siRNA expression levels in a descending order in regions of interest (FG-dDMR $\cap$ FG-*siren*, FG-dDMR-only, and FG-*siren*-only).
  - Generate heatmap for FRE occurrence arranged based on 24-nt siRNA expression in regions of interest (FG-dDMR $\cap$ FG-*siren*, FG-dDMR-only, and FG-*siren*-only).
  - Required input files:
    - `24nt_total.csv`: total 24-nt siRNA reads per sample
    - `count.FG-dDMR_FG-siren.txt`, `count.FG-dDMR_only.txt`, `count.FG-siren_only.txt`: 24-nt siRNA read counts in regions of interest
    - `FRE.FG-dDMR_FG-siren.bed`, `FRE.FG-dDMR_only.bed`, `FRE.FG-siren_only.bed`: FRE occurrence in regions of interest

> Rice *CLASSY* genes
>
> - Heatmap for transcriptional levels of *OsCLSY* genes: **Fig. 5a**
> - Lineplots for GC rates (DNA) and structural prediction scores (protein) of *OsCLSY3*: **Extended Data Fig. 8b**

- `expression_OsCLSY.R`:
  - Generate heatmap for expression of rice CLASSY genes.
  - Required input files:
    - `FG_OsCLSY.csv`: TPM values of OsCLSY genes in developing ovaries (FG1–FG6) and flag leaf
    - `EvoRepro_OsCLSY.csv`: TPM values of OsCLSY genes retrieved from [CoNekT database](https://evorepro.sbs.ntu.edu.sg)
- `get_GC_rate.py`:
  - Calculate GC rate within 25-bp windows.
  - Usage: `python get_GC_rate.py OsCLSY3_dna.fa 25 OsCLSY3_GC.csv`
    - `OsCLSY3_dna.fa`: fasta file of OsCLSY3 DNA sequence
    - `25`: 25-bp window to calculate GC rates
    - `OsCLSY3_GC.csv`: output file
- `plot_GC_rate.R`:
  - Generate lineplot for GC rates of OsCLSY3 DNA sequence.
  - Required input files:
    - `OsCLSY3_GC.csv`: output file generated using `get_GC_rate.py`
- `plot_pLDDT.R`:
  - Generate lineplot for pLDDT (confidence score) of OsCLSY3 protein structure predicted using AlphaFold3.
  - Required input files:
    - `confidence.csv`: confidence score of OsCLSY3 protein structure (exported from PyMOL)
