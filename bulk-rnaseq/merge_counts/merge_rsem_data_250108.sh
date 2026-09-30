#!/bin/bash

## Developer: Amy Olex
## Date: 1/8/25

## This is NOT a SLURM script

module load R/4.4.1

## Merge Human Gene Expression Data
Rscript /path/to/WCCTR_RNASeq_Pipeline/WorkstationScripts/merge_jax_counts.R -r 25.01.08_HUMAN_U54_AllSamples_nf-JaxPDXNet_RSEM_Gene \
-d /path/to/bulkRNASeq/jax_pipeline_results/human_gene_exp/ \
-o /path/to/bulkRNASeq/jax_pipeline_results/merged_data/ \
-s PDX \
-t GENE \
-i gene_id > 25.01.08_human_gene_merge.log


## Merge Mouse Gene Expression Data
Rscript /path/to/WCCTR_RNASeq_Pipeline/WorkstationScripts/merge_jax_counts.R -r 25.01.08_MOUSE_U54_AllSamples_nf-JaxPDXNet_RSEM_Gene \
-d /path/to/bulkRNASeq/jax_pipeline_results/mouse_gene_exp/ \
-o /path/to/bulkRNASeq/jax_pipeline_results/merged_data/ \
-s PDX \
-t GENE \
-i gene_id > 25.01.08_mouse_gene_merge.log

## Merge Human Isoform Expression Data
Rscript /path/to/WCCTR_RNASeq_Pipeline/WorkstationScripts/merge_jax_counts.R -r 25.01.08_HUMAN_U54_AllSamples_nf-JaxPDXNet_RSEM_Isoform \
-d /path/to/bulkRNASeq/jax_pipeline_results/human_isoform_exp/ \
-o /path/to/bulkRNASeq/jax_pipeline_results/merged_data/ \
-s PDX \
-t ISOFORM \
-i transcript_id > 25.01.08_human_isoform_merge.log


## Merge Mouse Isoform Expression Data
Rscript /path/to/WCCTR_RNASeq_Pipeline/WorkstationScripts/merge_jax_counts.R -r 25.01.08_MOUSE_U54_AllSamples_nf-JaxPDXNet_RSEM_Isoform \
-d /path/to/bulkRNASeq/jax_pipeline_results/mouse_isoform_exp/ \
-o /path/to/bulkRNASeq/jax_pipeline_results/merged_data/ \
-s PDX \
-t ISOFORM \
-i transcript_id > 25.01.08_mouse_isoform_merge.log
