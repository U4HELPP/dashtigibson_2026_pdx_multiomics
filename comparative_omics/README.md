# Comparative Multi-Omics Analysis

An R workflow for comparative multi-omics analysis of bulk RNA-seq and proteomics data. This repository includes scripts for dimensionality reduction (PCA and PLS-DA) and comparative Gene Set Variation Analysis (GSVA) across MSigDB Hallmark pathways, annotated by cancer type.

## Repository Structure

```text
.
├── 3d_plsda_pca.R               # 3D PLS-DA and PCA visualization for RNA-seq and proteomics data
├── gsva_correlation_heatmap.R   # GSVA, pathway correlation, & heatmap export
├── sample_metadata.csv          # Metadata mapping sample IDs to cancer types
└── results/                     # Directory for exported CSV matrices and visual plots
