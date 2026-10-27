# ==============================================================================
# Script: Multi-Omics 3D PCA and PLS-DA Analysis (RNA-seq & Proteomics)
# Description: Performs dimensionality reduction (PCA & PLS-DA) on bulk RNA-seq 
#              and Proteomics datasets and exports 3D interactive and static plots.
# ==============================================================================

# --- 1. Dependencies & Package Setup ---
suppressPackageStartupMessages({
  library(mixOmics)
  library(plotly)
  library(dplyr)
  library(tools)
  library(htmlwidgets)
  library(webshot2)
})

# --- 2. Configuration & Paths ---
RNA_EXPR_FILE  <- "25.11.20_HUMAN_U54_AllSamples_nf-JaxPDXNet_RSEM_Gene_TPM_U54manuscript.csv"
PROT_EXPR_FILE <- "120325_PCA_Normalization_Averaged_renamed.csv"
META_FILE      <- "sample_metadata.csv"
OUTPUT_DIR     <- "results"
N_COMP         <- 3

PROT_OUTLIERS  <- c("VCU-PC-064", "VCU-CO-023", "VCU-MC-130")

CANCER_COLORS  <- c(
  "Breast"      = "#cb027e", "Colorectal"  = "#0081cb",
  "Gallbladder" = "#ff924d", "Lung"        = "#74cc00",
  "Melanoma"    = "#784803", "Ovarian"     = "#ca6ce7",
  "Pancreatic"  = "#face20"
)

# Ensure output directory exists
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# --- 3. Shared Helper Functions ---

#' Extract Explained Variance Vector
get_prop_var <- function(model_obj) {
  var_raw <- model_obj$prop_expl_var
  if (is.list(var_raw)) {
    return(if (!is.null(var_raw$X)) var_raw$X else unlist(var_raw))
  }
  return(var_raw)
}

#' Generate Standardized 3D Plotly Scatter Plot
plot_3d_scatter <- function(df, x_col, y_col, z_col, group_col, sample_col, 
                            title_text, x_label, y_label, z_label, color_palette) {
  plot_ly(
    data = df, 
    x = as.formula(paste0("~", x_col)), 
    y = as.formula(paste0("~", y_col)), 
    z = as.formula(paste0("~", z_col)), 
    color = as.formula(paste0("~", group_col)), 
    colors = color_palette, 
    text = as.formula(paste0("~paste('Sample:', ", sample_col, ", '<br>Type:', ", group_col, ")")),
    type = "scatter3d", 
    mode = "markers",
    marker = list(size = 10, opacity = 0.9, line = list(width = 1, color = "DarkSlateGrey"))
  ) %>%
    layout(
      title = title_text,
      scene = list(
        xaxis = list(title = x_label),
        yaxis = list(title = y_label),
        zaxis = list(title = z_label)
      )
    )
}

#' Helper to export interactive HTML widgets and static PNG snapshots
export_visualization <- function(fig, prefix, date_str) {
  html_path <- file.path(OUTPUT_DIR, sprintf("%s_%s_Interactive.html", date_str, prefix))
  png_path  <- file.path(OUTPUT_DIR, sprintf("%s_%s_Static.png", date_str, prefix))
  
  saveWidget(fig, html_path, selfcontained = TRUE)
  webshot2::webshot(html_path, file = png_path, vwidth = 1200, vheight = 900, delay = 2)
}

#' Master Pipeline Runner for Dimensionality Reduction
run_dim_reduction_pipeline <- function(X, metadata_df, datatype_label, export_prefix, date_str) {
  message(sprintf("\n--- Processing %s ---", datatype_label))
  
  # Format response vector Y
  metadata_df$Cancer_type <- tools::toTitleCase(metadata_df$Cancer_type)
  Y                       <- as.factor(metadata_df$Cancer_type)
  
  # Drop zero-variance features
  X <- X[, apply(X, 2, var) > 0]
  
  # Run PLS-DA & PCA
  message("Running PLS-DA...")
  plsda_result <- plsda(X, Y, ncomp = N_COMP)
  
  message("Running PCA...")
  pca_result <- mixOmics::pca(X, ncomp = N_COMP, center = TRUE, scale = TRUE)
  
  # Variance Check
  plsda_var <- get_prop_var(plsda_result)
  pca_var   <- get_prop_var(pca_result)
  
  message(sprintf("PLS-DA PC1: %.2f%%", plsda_var[1] * 100))
  message(sprintf("PCA    PC1: %.2f%%", pca_var[1] * 100))
  
  # Build Plotly Objects
  plsda_df <- data.frame(
    Comp1 = plsda_result$variates$X[, 1],
    Comp2 = plsda_result$variates$X[, 2],
    Comp3 = plsda_result$variates$X[, 3],
    Group = Y, 
    SampleID = metadata_df$Sample
  )
  
  fig_plsda <- plot_3d_scatter(
    df = plsda_df, x_col = "Comp1", y_col = "Comp2", z_col = "Comp3",
    group_col = "Group", sample_col = "SampleID", 
    title_text = sprintf("3D PLS-DA: %s", datatype_label),
    x_label = sprintf("Comp 1 (%.1f%%)", plsda_var[1] * 100),
    y_label = sprintf("Comp 2 (%.1f%%)", plsda_var[2] * 100),
    z_label = sprintf("Comp 3 (%.1f%%)", plsda_var[3] * 100),
    color_palette = CANCER_COLORS
  )
  
  pca_df <- data.frame(
    PC1 = pca_result$variates$X[, 1],
    PC2 = pca_result$variates$X[, 2],
    PC3 = pca_result$variates$X[, 3],
    Group = Y, 
    SampleID = metadata_df$Sample
  )
  
  fig_pca <- plot_3d_scatter(
    df = pca_df, x_col = "PC1", y_col = "PC2", z_col = "PC3",
    group_col = "Group", sample_col = "SampleID", 
    title_text = sprintf("3D PCA: %s", datatype_label),
    x_label = sprintf("PC 1 (%.1f%%)", pca_var[1] * 100),
    y_label = sprintf("PC 2 (%.1f%%)", pca_var[2] * 100),
    z_label = sprintf("PC 3 (%.1f%%)", pca_var[3] * 100),
    color_palette = CANCER_COLORS
  )
  
  # Export Visualizations
  message(sprintf("Exporting %s plots to '%s'...", datatype_label, OUTPUT_DIR))
  export_visualization(fig_plsda, sprintf("PLSDA_%s", export_prefix), date_str)
  export_visualization(fig_pca, sprintf("PCA_%s", export_prefix), date_str)
}

# --- 4. Main Execution ---
date_str <- format(Sys.Date(), "%Y-%m-%d")

if (!file.exists(META_FILE)) stop(paste("Metadata file not found:", META_FILE))
metadata_raw <- read.csv(META_FILE, stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# Part A: Bulk RNA-seq Pipeline
# ------------------------------------------------------------------------------
if (file.exists(RNA_EXPR_FILE)) {
  message("\n==========================================")
  message("Processing Bulk RNA-seq Data")
  message("==========================================")
  
  rna_df <- read.csv(RNA_EXPR_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(rna_df) <- rna_df$gene_id
  X_rna <- t(rna_df[, 3:ncol(rna_df)])
  
  # Prefix sample matching against metadata
  meta_idx_rna <- sapply(rownames(X_rna), function(s_name) {
    idx <- which(startsWith(s_name, metadata_raw$Sample))
    if (length(idx) > 0) return(idx[1]) else return(NA)
  })
  
  keep_rna <- !is.na(meta_idx_rna)
  X_rna        <- X_rna[keep_rna, ]
  meta_rna_df  <- metadata_raw[meta_idx_rna[keep_rna], ]
  
  # Log2 transform expression matrix
  X_rna <- log2(X_rna + 1)
  
  run_dim_reduction_pipeline(X_rna, meta_rna_df, "RNA-seq", "RNA", date_str)
} else {
  warning(paste("RNA expression file not found, skipping:", RNA_EXPR_FILE))
}

# ------------------------------------------------------------------------------
# Part B: Proteomics Pipeline
# ------------------------------------------------------------------------------
if (file.exists(PROT_EXPR_FILE)) {
  message("\n==========================================")
  message("Processing Proteomics Data")
  message("==========================================")
  
  prot_df <- read.csv(PROT_EXPR_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(prot_df) <- prot_df$Protein_ID
  X_prot <- t(prot_df[, 2:ncol(prot_df)])
  
  # Outlier removal
  if (length(PROT_OUTLIERS) > 0) {
    X_prot <- X_prot[!rownames(X_prot) %in% PROT_OUTLIERS, ]
    message(sprintf("Removed %d proteomics outliers: %s", 
                    length(PROT_OUTLIERS), paste(PROT_OUTLIERS, collapse = ", ")))
  }
  
  # Direct sample matching against metadata
  X_prot       <- X_prot[rownames(X_prot) %in% metadata_raw$Sample, ]
  meta_prot_df <- metadata_raw[match(rownames(X_prot), metadata_raw$Sample), ]
  
  run_dim_reduction_pipeline(X_prot, meta_prot_df, "Proteomics", "protein", date_str)
} else {
  warning(paste("Proteomics expression file not found, skipping:", PROT_EXPR_FILE))
}

message("\nMulti-omics 3D PCA and PLS-DA pipeline completed successfully.")
