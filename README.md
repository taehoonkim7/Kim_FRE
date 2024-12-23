# A FRE guides de novo DNA methylation during rice gametophyte development

- Codes and datasets for <https://doi.org/10.1101/2024.12.03.626448> (to be updated once published)
- FRE: female-specific RdDM-associated element

## Pipeline

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

### Venn_UpSet

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

- `process_CX_meta_TE.py`:
