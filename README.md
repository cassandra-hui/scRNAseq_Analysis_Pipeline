# scRNAseq Analysis Pipeling

## Overview

This repository contains scripts and templates for analyzing single-cell RNA sequencing (scRNA-seq) data. The pipeline processes 10X Genomics or Parse data through quality control, normalization, dimensionality reduction, clustering, cell type annotation, and differential expression analysis.


## Analysis Workflow 

Start in the `03_Run` folder and follow the sctipts in `00_scripts`. This will walk through:

0. Data Loading: Import 10X Genomics or Parse data 
1. Quality Control: Filter cells based on quality metrics
2. Normalization: Using SCTransform or LogNormalize methods
3. Dimensionality Reduction: PCA and UMAP visualization
4. Integration (Optional): Batch effect correction
5. Cell Type Annotation: Manual or automated (SingleR, scMRMA)
6. Merging or Breaking Clusters (Optional): Based on reasercher's goals
7. Differential Expression Analysis: Compare expression between groups
8. R Shiny App Setup: The scRNA_shiny app for exporation of your data 

## Report Generation 

The `scRNAseq_Project.Rmd` file serves as a template for creating comprehensive analysis reports. It pulls R objects created from the scripts and displays the information and figures in a rendered HTML file for the reasercher. 

