# High-Throughput Drug Screening Analysis

An R workflow for high-throughput drug screening analysis across PDX models. This directory includes scripts for differential breast cancer drug response analysis by genetic ancestry and cross-cancer lineage statistical modeling across 1 µM and 10 µM screening concentrations.

## Repository Structure

```text
.
├── drugscreen_ancestry_analysis.R  # Differential breast cancer drug response analysis (EUR vs. AFR) & volcano plots
├── drugscreen_lineage_analysis.R   # Cross-lineage drug response ANOVA, Tukey HSD, Z-score heatmaps, & boxplots
└── sample_metadata.csv             # Metadata mapping sample IDs to cancer types and genetic ancestry groups
