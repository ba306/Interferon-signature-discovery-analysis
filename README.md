# Expression signatures with specificity for type I and II IFN response and relevance for autoimmune diseases and cancer 

### Authors: 
Bogac Aybey, Benedikt Brors, and Eike Staub
> **Status**: Manuscript accepted

## Overview

This repository contains all scripts and data necessary to reproduce the analyses presented in our study on discovering interferon (IFN) response signatures. The study focuses on identifying gene expression signatures specific to Type I and Type II IFN responses. We evaluate the relevance our IFN-I and IFN-II and published IFN signatures in autoimmune diseases and cancer along multiple contexts e.g., platforms, cell types, and experimental setups.

## Repository Structure

- `data/`: Raw and processed datasets.
- `discovery/`: Scripts for identifying IFN-specific gene signatures.
- `functions/`: Custom R functions used throughout the analyses.
- `greenland_scRNAseq_analysis/`: Analysis scripts for the Greenland scRNA-seq dataset.
- `kartha_scRNAseq_analysis/`: Analysis scripts for the Kartha scRNA-seq dataset.
- `published_signatures/`: Published IFN gene signatures used for comparison.
- `SLE_analysis/`: Scripts for analyzing SLE datasets.
- `Tumor_bulk_analysis/`: Scripts for bulk RNA-seq analysis of tumor samples (TCGA and ICI trials).


## Data Availability

Raw and processed data are included in the `data/` folder or can be downloaded using links provided in the scripts.

## Packages and Versions

This project uses the following R packages with the specified versions:

| Package              | Version    |   | Package              | Version    |
|----------------------|------------|---|----------------------|------------|
| dplyr                | 1.1.4      |   | ggplot2              | 3.4.4      |
| foreach              | 1.5.2      |   | gplots               | 3.1.1      |
| qs                   | 0.25.5     |   | limma                | 3.50.3     |
| Seurat               | 4.0.5      |   | netmeta              | 3.2-0      |
| SeuratDisk           | 0.0.0.9019 |   | parallel             | 4.1.2      |
| data.table           | 1.14.2     |   | stringr              | 1.5.1      |
| DESeq2               | 1.34.0     |   | tidyverse            | 2.0.0      |
| ComplexHeatmap       | 2.10.0     |   | randomForest         | 4.7-1.2    |
| cowplot              | 1.1.1      |   | emmeans              | 1.8.6      |
| grid                 | 4.1.2      |   | ggrepel              | 0.9.3      |
| gridtext             | 0.1.5      |   | patchwork            | 1.1.1      |
| RColorBrewer         | 1.1-2      |   | uwot                 | 0.1.11     |
| ggpubr               | 0.6.0      |   | ggbeeswarm           | 0.7.2      |
| ggsci                | 3.2.0      |   | readxl               | 1.4.3      |
| gridExtra            | 2.3        |   | GEOquery             | 2.62.2     |
| reshape2             | 1.4.4      |   | hgu133plus2.db       | 3.13.0     |
| rstatix              | 0.7.2      |   | illuminaHumanv4.db   | 1.26.0     |
| varhandle            | 2.0.6      |   | pROC                 | 1.18.5     |
| doParallel           | 1.0.17     |   | circlize             | 0.4.16     |
| base                 | 4.1.2      |   | tibble               | 3.2.1      |
| tidyr                | 1.3.1      |   | correlation          | 0.8.7      |

The following packages are proprietary and not publicly available:

- `xopggexpdatasets`
- `xopdata`
- `ProfilerAPI2`

These packages are maintained internally and include code that connects to internal SQL databases and internal data systems. As such, they **cannot be shared or installed externally**. To fully reproduce the analysis, access to these internal packages and the corresponding data infrastructure is required.


## Contact

Bogac Aybey
[bogac.aybey@stud.uni-heidelberg.de](mailto:bogac.aybey@stud.uni-heidelberg.de)

