# ==============================================================================
# Script: Comparative Hallmark GSVA Analysis (RNA-seq vs. Proteomics)
# Description: Preprocesses log2-transformed expression matrices, performs GSVA 
#              scoring across Hallmark pathways, calculates gene set coverage & 
#              RNA-protein correlations, and exports annotated heatmap plots.
# Author: Nina Dashti-Gibson
# ==============================================================================

# --- 1. Dependencies & Package Setup ---
suppressPackageStartupMessages({
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!require("GSVA", quietly = TRUE)) BiocManager::install("GSVA")
  if (!require("msigdbr", quietly = TRUE)) BiocManager::install("msigdbr")
  if (!require("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
  if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!require("tibble", quietly = TRUE)) install.packages("tibble")
  if (!require("readr", quietly = TRUE)) install.packages("readr")
  if (!require("circlize", quietly = TRUE)) install.packages("circlize")
})

library(GSVA)
library(msigdbr)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)
library(readr)

# --- 2. Configuration & Paths ---
RNA_FILE   <- "25.11.20_HUMAN_U54_AllSamples_nf-JaxPDXNet_RSEM_Gene_TPM_U54manuscript.csv"
PROT_FILE  <- "120325_DetectedProteins_MaxSum_RepAvgImpute_RepAveraged_renamed.csv"
META_FILE  <- "sample_metadata.csv"
OUTPUT_DIR <- "results"

RNA_GENE_COL  <- "geneSym"
PROT_GENE_COL <- "Gene"

CANCER_COLORS <- c(
  "Breast"      = "#cb027e", "Colorectal"  = "#0081cb",
  "Gallbladder" = "#ff924d", "Lung"        = "#74cc00",
  "Melanoma"    = "#784803", "Ovarian"     = "#ca6ce7",
  "Pancreatic"  = "#face20"
)

TYPE_COLORS <- c(RNA = "#E69F00", Protein = "#56B4E9")

# Ensure output directory exists
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# --- 3. Helper & Preprocessing Functions ---

#' Preprocess Expression Data (Aggregation & Log2 Transformation)
preprocess_expression_data <- function(df, type_name, gene_col) {
  if (!gene_col %in% colnames(df)) {
    stop(sprintf("The column '%s' is missing in %s dataset.", gene_col, type_name))
  }
  
  df <- df %>% rename(Gene = all_of(gene_col))
  df <- df[!is.na(df$Gene) & df$Gene != "", ]
  
  # Deduplicate genes via mean aggregation across numeric columns
  df <- df %>% 
    group_by(Gene) %>% 
    summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
    ungroup()
  
  mat <- df %>% column_to_rownames("Gene") %>% as.matrix()
  
  # Log2 transform data (log2(x + 1) handles zero expression)
  mat <- log2(mat + 1)
  return(mat)
}

#' Match & Filter Matrix Sample Columns Against Metadata
align_to_metadata <- function(mat, meta_samples, dataset_name) {
  cols <- colnames(mat)
  mapped_ids <- sapply(cols, function(cn) {
    if (cn %in% meta_samples) return(cn)
    matches <- meta_samples[sapply(meta_samples, function(m) startsWith(cn, m))]
    if (length(matches) > 0) return(matches[1]) else return(NA)
  })
  
  valid_cols <- !is.na(mapped_ids)
  if (sum(valid_cols) == 0) {
    stop(sprintf("No matching samples found in metadata for %s", dataset_name))
  }
  
  mat_sub <- mat[, valid_cols, drop = FALSE]
  colnames(mat_sub) <- mapped_ids[valid_cols]
  
  # Deduplicate sample columns by averaging replicates matching the same Sample ID
  unique_ids <- unique(colnames(mat_sub))
  if (length(unique_ids) < ncol(mat_sub)) {
    mat_sub <- do.call(cbind, lapply(unique_ids, function(id) {
      cols_for_id <- mat_sub[, colnames(mat_sub) == id, drop = FALSE]
      rowMeans(cols_for_id, na.rm = TRUE)
    }))
    colnames(mat_sub) <- unique_ids
  }
  
  return(mat_sub)
}

#' Calculate Gene Set Coverage Percentage
calculate_coverage <- function(pathway_name, detected_genes, gene_sets) {
  set_genes  <- gene_sets[[pathway_name]]
  n_in_set   <- length(set_genes)
  n_detected <- sum(set_genes %in% detected_genes)
  return((n_detected / n_in_set) * 100)
}

# --- 4. Data Loading & Preprocessing ---
message("Loading dataset and metadata files...")
if (!file.exists(RNA_FILE))  stop(paste("RNA file not found:", RNA_FILE))
if (!file.exists(PROT_FILE)) stop(paste("Proteomics file not found:", PROT_FILE))
if (!file.exists(META_FILE)) stop(paste("Metadata file not found:", META_FILE))

rna_data  <- read.csv(RNA_FILE, check.names = FALSE, stringsAsFactors = FALSE)
prot_data <- read.csv(PROT_FILE, check.names = FALSE, stringsAsFactors = FALSE)
meta_data <- read.csv(META_FILE, stringsAsFactors = FALSE)

# Format and log2-transform expression matrices
rna_mat  <- preprocess_expression_data(rna_data, "RNA", gene_col = RNA_GENE_COL)
prot_mat <- preprocess_expression_data(prot_data, "Protein", gene_col = PROT_GENE_COL)

# Filter matrices to include only samples listed in metadata
valid_metadata_samples <- unique(meta_data$Sample)
rna_mat  <- align_to_metadata(rna_mat, valid_metadata_samples, "RNA")
prot_mat <- align_to_metadata(prot_mat, valid_metadata_samples, "Protein")

# Intersect common samples across both modalities and metadata
common_samples <- intersect(colnames(rna_mat), colnames(prot_mat))
message(sprintf("Identified %d matched samples between RNA-seq, Proteomics, and metadata.", length(common_samples)))

if (length(common_samples) == 0) {
  stop("No common samples found across RNA, Protein, and metadata files.")
}

rna_mat  <- rna_mat[, common_samples, drop = FALSE]
prot_mat <- prot_mat[, common_samples, drop = FALSE]

# --- 5. Hallmark Gene Sets & GSVA Calculation ---
message("Retrieving Hallmark gene sets from msigdbr...")
all_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
gs_list       <- split(x = all_gene_sets$gene_symbol, f = all_gene_sets$gs_name)

message("Calculating RNA GSVA scores...")
rna_params <- gsvaParam(exprData = rna_mat, geneSets = gs_list, kcdf = "Gaussian")
rna_gsva   <- gsva(rna_params)

message("Calculating Protein GSVA scores...")
prot_params <- gsvaParam(exprData = prot_mat, geneSets = gs_list, kcdf = "Gaussian")
prot_gsva   <- gsva(prot_params)

if (!is.matrix(rna_gsva))  rna_gsva  <- assay(rna_gsva)
if (!is.matrix(prot_gsva)) prot_gsva <- assay(prot_gsva)

# --- 6. Export Full GSVA Matrices & Correlation Summary ---
date_str <- format(Sys.Date(), "%Y-%m-%d")

common_pathways <- intersect(rownames(rna_gsva), rownames(prot_gsva))

rna_gsva_out  <- rna_gsva[common_pathways, common_samples, drop = FALSE]
prot_gsva_out <- prot_gsva[common_pathways, common_samples, drop = FALSE]

# Export raw GSVA score tables
rna_gsva_df <- tibble::rownames_to_column(as.data.frame(rna_gsva_out, check.names = FALSE), var = "Hallmark_Signature")
prot_gsva_df <- tibble::rownames_to_column(as.data.frame(prot_gsva_out, check.names = FALSE), var = "Hallmark_Signature")

readr::write_csv(rna_gsva_df,  file.path(OUTPUT_DIR, sprintf("%s_Hallmark_GSVA_RNA_AllSignatures.csv", date_str)))
readr::write_csv(prot_gsva_df, file.path(OUTPUT_DIR, sprintf("%s_Hallmark_GSVA_Protein_AllSignatures.csv", date_str)))

# Compute statistical correlation & coverage for all Hallmark signatures
rna_genes_detected  <- rownames(rna_mat)
prot_genes_detected <- rownames(prot_mat)

all_pathway_metrics <- lapply(common_pathways, function(pw) {
  r_vec <- rna_gsva[pw, common_samples]
  p_vec <- prot_gsva[pw, common_samples]
  
  test_res <- tryCatch(
    cor.test(r_vec, p_vec, method = "pearson", use = "pairwise.complete.obs"),
    error = function(e) list(estimate = NA, p.value = NA)
  )
  
  pct_rna  <- calculate_coverage(pw, rna_genes_detected, gs_list)
  pct_prot <- calculate_coverage(pw, prot_genes_detected, gs_list)
  
  data.frame(
    Hallmark_Signature     = pw,
    Pearson_R              = round(as.numeric(test_res$estimate), 4),
    P_Value                = as.numeric(test_res$p.value),
    RNA_Gene_Coverage_Pct  = round(pct_rna, 2),
    Prot_Gene_Coverage_Pct = round(pct_prot, 2),
    stringsAsFactors       = FALSE
  )
})

correlation_df <- do.call(rbind, all_pathway_metrics) %>% arrange(desc(Pearson_R))
readr::write_csv(correlation_df, file.path(OUTPUT_DIR, sprintf("%s_Hallmark_Pathways_Correlation_and_Coverage.csv", date_str)))

message("Exported full GSVA score matrices and correlation summary to /", OUTPUT_DIR)

# --- 7. Select Heatmap Signatures ---
manual_sets <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_NOTCH_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_APICAL_JUNCTION",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HYPOXIA",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_MYOGENESIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_PROTEIN_SECRETION",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
)

if (!is.null(manual_sets)) {
  selected_pathways <- intersect(manual_sets, common_pathways)
  if (length(selected_pathways) == 0) stop("None of the manually selected pathways were found in the data.")
} else {
  top_n <- 20
  selected_pathways <- correlation_df$Hallmark_Signature[1:top_n]
}

# --- 8. Prepare Heatmap Matrices & Annotations ---
processed_data_list <- list()

for (i in seq_along(selected_pathways)) {
  pw <- selected_pathways[i]
  
  rna_scores  <- rna_gsva[pw, ]
  prot_scores <- prot_gsva[pw, ]
  
  pct_rna  <- calculate_coverage(pw, rna_genes_detected, gs_list)
  pct_prot <- calculate_coverage(pw, prot_genes_detected, gs_list)
  r_val    <- cor(rna_scores, prot_scores, method = "pearson", use = "pairwise.complete.obs")
  if (is.na(r_val)) r_val <- 0 
  
  processed_data_list[[length(processed_data_list) + 1]] <- list(
    ScoreVector = rna_scores, RowName = paste0(pw, " (RNA)"),
    Type = "RNA", Pathway = pw, Coverage = pct_rna, Correlation = r_val
  )
  processed_data_list[[length(processed_data_list) + 1]] <- list(
    ScoreVector = prot_scores, RowName = paste0(pw, " (Prot)"),
    Type = "Protein", Pathway = pw, Coverage = pct_prot, Correlation = r_val
  )
}

combined_scores <- do.call(rbind, lapply(processed_data_list, function(x) x$ScoreVector))
rownames(combined_scores) <- sapply(processed_data_list, function(x) x$RowName)

row_annotation_type <- sapply(processed_data_list, function(x) x$Type)
row_split_vector     <- sapply(processed_data_list, function(x) x$Pathway)
coverage_pct_vector  <- sapply(processed_data_list, function(x) x$Coverage)
cor_vector           <- sapply(processed_data_list, function(x) x$Correlation)

# Align Metadata Order
meta_filtered <- meta_data[meta_data$Sample %in% colnames(combined_scores), ]
meta_filtered$Cancer_type <- tools::toTitleCase(meta_filtered$Cancer_type)
meta_filtered <- meta_filtered %>% arrange(Cancer_type, Sample)
combined_scores <- combined_scores[, meta_filtered$Sample]

cor_col_fun <- colorRamp2(c(0, 0.5, 1), c("gold", "turquoise", "darkmagenta"))

# Top Annotation (Cancer Type)
ha_top <- HeatmapAnnotation(
  `Cancer type` = meta_filtered$Cancer_type,
  col = list(`Cancer type` = CANCER_COLORS),
  show_annotation_name = FALSE,
  annotation_name_side = "left",
  simple_anno_size = unit(0.4, "cm"),
  annotation_name_gp = gpar(fontsize = 12),
  annotation_legend_param = list(
    `Cancer type` = list(
      nrow = 2, 
      title_position = "topleft",
      title_gp = gpar(fontsize = 11),
      labels_gp = gpar(fontsize = 10)
    )
  )
)

# Left Annotation (Correlation & Data Type)
ha_left <- rowAnnotation(
  Correlation = cor_vector,
  `Data type` = row_annotation_type,
  col = list(`Data type` = TYPE_COLORS, Correlation = cor_col_fun),
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  annotation_name_rot = 90, 
  annotation_name_gp = gpar(fontsize = 10),
  simple_anno_size = unit(0.4, "cm"),
  annotation_legend_param = list(
    Correlation = list(
      direction = "horizontal", 
      title_position = "topcenter", 
      legend_width = unit(3, "cm"),
      title_gp = gpar(fontsize = 11),
      labels_gp = gpar(fontsize = 10)
    ),
    `Data type` = list(
      title_gp = gpar(fontsize = 11),
      labels_gp = gpar(fontsize = 10)
    )
  )
)

# Right Annotation (Gene Set Coverage Barplot)
bar_fill_colors <- ifelse(row_annotation_type == "RNA", TYPE_COLORS["RNA"], TYPE_COLORS["Protein"])

ha_right <- rowAnnotation(
  `Features present (%)` = anno_barplot(
    coverage_pct_vector,
    gp = gpar(fill = bar_fill_colors, col = NA),
    width = unit(3, "cm"),
    border = TRUE,
    axis_param = list(
      side = "top", 
      at = c(0, 50, 100), 
      labels_rot = 0,
      gp = gpar(fontsize = 8)
    ) 
  ),
  annotation_name_side = "top", 
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 8)
)

# --- 9. Export Heatmap Plot ---
png_path <- file.path(OUTPUT_DIR, sprintf("%s_Hallmark_RNA_vs_Protein_Annotated_Barplot.png", date_str))

png(png_path, width = 14, height = 12, units = "in", res = 300)

ht <- Heatmap(
  combined_scores,
  name = "GSVA score",
  row_split = factor(row_split_vector, levels = selected_pathways),
  cluster_row_slices = FALSE,
  
  row_title_side = "right",
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 11),
  
  row_gap = unit(2.5, "mm"),
  top_annotation = ha_top,
  left_annotation = ha_left,
  right_annotation = ha_right,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,   
  
  show_row_names = FALSE,
  
  show_column_names = TRUE,
  column_names_side = "top", 
  column_names_gp = gpar(fontsize = 10),
  
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  
  heatmap_legend_param = list(
    direction = "horizontal", 
    title_position = "topcenter", 
    legend_width = unit(3, "cm"),
    title_gp = gpar(fontsize = 11),
    labels_gp = gpar(fontsize = 10)
  )
)

draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

dev.off()

message("Analysis complete! All outputs successfully generated in /", OUTPUT_DIR)
