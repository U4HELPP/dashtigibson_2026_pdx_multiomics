# Bulk RNA-Seq Processing — README

**Repo/paths referenced:**
- Alignment pipeline: [jds-nf-workflows (Jackson Laboratory)](https://github.com/TheJacksonLaboratory/jds-nf-workflows)
- Protocol: [RNA-seq Step-by-Step Guide (VCU Bioinformatics Core)](https://github.com/VCU-Bioinformatics-Core/bioinformatics-epicenter-pipelines/wiki/RNA-seq-Step-by-Step-Guide)
- Merge script: [`bulk-rnaseq/merge_counts`](https://github.com/U4HELPP/dashtigibson_2026_pdx_multiomics/tree/main/bulk-rnaseq/merge_counts) *(scripts taken from SOP v1.1 release)*
- Downstream analysis code: [`bulk-rnaseq/analysis`](https://github.com/U4HELPP/dashtigibson_2026_pdx_multiomics/tree/main/bulk-rnaseq/analysis)

## Purpose

Documents the processing steps used for bulk RNA-seq samples in the PDX multi-omics manuscript (*Pan-cancer comparisons of multi-omics data for therapeutic selection*, Dashti-Gibson et al.).

## Pipeline Overview

| Step | Tool/Script | Input | Output |
|---|---|---|---|
| 1 | [jds-nf-workflows](https://github.com/TheJacksonLaboratory/jds-nf-workflows) *(JAX PDX Nextflow pipeline; see [protocol](https://github.com/VCU-Bioinformatics-Core/bioinformatics-epicenter-pipelines/wiki/RNA-seq-Step-by-Step-Guide) for step-by-step usage)* | Raw bulk RNA-seq FastQ files | QC reports, alignments, human/mouse gene and isoform count files |
| 2 | [`merge_counts`](https://github.com/U4HELPP/dashtigibson_2026_pdx_multiomics/tree/main/bulk-rnaseq/merge_counts) *(SOP v1.1 release)* | JAX per-sample count files (see directory structure below) | Single merged count matrix |
| 3 | [`analysis`](https://github.com/U4HELPP/dashtigibson_2026_pdx_multiomics/tree/main/bulk-rnaseq/analysis) | Merged count matrix | PAM50 subtyping, PCA plots, PLS-DA analysis, downstream results |

## Merge Step — Required Directory Structure

The merge script requires JAX pipeline outputs to be consolidated into the following structure via symlinks:

```
jaxresults/
    ├── human_bam/
    ├── human_gene_exp/
    ├── human_isoform_exp/
    ├── mouse_bam/
    ├── mouse_gene_exp/
    ├── mouse_isoform_exp/
```

> **Note:** The cluster-specific SOP used to generate this directory structure and submit the merge job is institution-specific and cannot be shared publicly. A template submission script (SOP v1.1) is included in the repo, stripped of all cluster-specific information, to illustrate the general submission pattern.

## Merge Step — Script Descriptions

| Script | Role |
|---|---|
| `link_jax_results.sh` | Run first. Creates the symlinked directory structure above from the raw JAX results, then generates and executes `merge_rsem_data_250108.sh`. |
| `merge_rsem_data_250108.sh` | Template script (cluster-specific info stripped) containing the merge commands for each file type (gene counts, isoform counts, TPM, etc.), called automatically at the end of `link_jax_results.sh`. |
| [`merge_jax_counts.R`](https://github.com/AmyOlex/WCCTR_RNASeq_Pipeline/blob/master/WorkstationScripts/merge_jax_counts.R) *(WCCTR_RNASeq_Pipeline repo)* | Performs the actual merging of per-sample count files into the single merged count matrix; called by `merge_rsem_data_250108.sh`. |

## Analysis
 
| Step | Tool/Script | Description |
|---|---|---|
| 1 | [`06_PAM50_subtyping.R`](https://github.com/AmyOlex/WCCTR_RNASeq_Pipeline/blob/master/WorkstationScripts/06_PAM50_subtyping.R) *(WCCTR_RNASeq_Pipeline repo)* | PAM50 classification |
 

## Referenced URLs

| Reference | URL |
|---|---|
| jds-nf-workflows (Jackson Laboratory) | https://github.com/TheJacksonLaboratory/jds-nf-workflows |
| RNA-seq Step-by-Step Guide (VCU Bioinformatics Core) | https://github.com/VCU-Bioinformatics-Core/bioinformatics-epicenter-pipelines/wiki/RNA-seq-Step-by-Step-Guide |
| `bulk-rnaseq/merge_counts` (SOP v1.1 release) | https://github.com/U4HELPP/dashtigibson_2026_pdx_multiomics/tree/main/bulk-rnaseq/merge_counts |
| `bulk-rnaseq/analysis` | https://github.com/U4HELPP/dashtigibson_2026_pdx_multiomics/tree/main/bulk-rnaseq/analysis |
| `merge_jax_counts.R` (WCCTR_RNASeq_Pipeline repo) | https://github.com/AmyOlex/WCCTR_RNASeq_Pipeline/blob/master/WorkstationScripts/merge_jax_counts.R |
| `06_PAM50_subtyping.R` (WCCTR_RNASeq_Pipeline repo) | https://github.com/AmyOlex/WCCTR_RNASeq_Pipeline/blob/master/WorkstationScripts/06_PAM50_subtyping.R |
