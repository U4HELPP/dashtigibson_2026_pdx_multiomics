# ==============================================================================
# Script: Proteomics 3D PCA and PLS-DA Analysis
# Description: Performs dimensionality reduction (PCA & PLS-DA) on proteomics 
#              data and exports interactive and static 3D plots.
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
EXPR_FILE         <- "120325_PCA_Normalization_Averaged_renamed.csv"
META_FILE         <- "sample_metadata.csv"
OUTPUT_DIR        <- "results"
N_COMP            <- 3
SAMPLES_TO_REMOVE <- c("VCU-PC-064", "VCU-CO-023", "VCU-MC-130")

CANCER_COLORS <- c(
  "Breast"      = "#cb027e", "Colorectal"  = "#0081cb",
  "Gallbladder" = "#ff924d", "Lung"        = "#74cc00",
  "Melanoma"    = "#784803", "Ovarian"     = "#ca6ce7",
  "Pancreatic"  = "#face20"
)

# Ensure output directory exists
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# --- 3. Helper Functions ---

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

# --- 4. Data Loading & Preprocessing ---
message("Loading expression and metadata files...")
if (!file.exists(EXPR_FILE)) stop(paste("Expression file not found:", EXPR_FILE))
if (!file.exists(META_FILE)) stop(paste("Metadata file not found:", META_FILE))

expression_df <- read.csv(EXPR_FILE, check.names = FALSE, stringsAsFactors = FALSE)
metadata_df   <- read.csv(META_FILE, stringsAsFactors = FALSE)

# Format expression matrix (Proteomics uses Protein_ID in Col 1)
rownames(expression_df) <- expression_df$Protein_ID
expression_matrix       <- expression_df[, 2:ncol(expression_df)]
X                       <- t(expression_matrix)

# Explicitly remove specific outlier samples if defined
if (length(SAMPLES_TO_REMOVE) > 0) {
  X <- X[!rownames(X) %in% SAMPLES_TO_REMOVE, ]
  message(sprintf("Removed %d specified samples: %s", 
                  length(SAMPLES_TO_REMOVE), paste(SAMPLES_TO_REMOVE, collapse = ", ")))
}

# Filter matrix to match samples in metadata
X           <- X[rownames(X) %in% metadata_df$Sample, ]
metadata_df <- metadata_df[match(rownames(X), metadata_df$Sample), ]

if (!all(rownames(X) == metadata_df$Sample)) {
  stop("Error: Sample names between expression matrix and metadata do not match!")
}

# Process response labels
metadata_df$Cancer_type <- tools::toTitleCase(metadata_df$Cancer_type)
Y                       <- as.factor(metadata_df$Cancer_type)

# Drop zero-variance proteins
X <- X[, apply(X, 2, var) > 0]

# --- 5. Run Analyses ---
message("Running PLS-DA on Proteomics data...")
plsda_result <- plsda(X, Y, ncomp = N_COMP)

message("Running PCA on Proteomics data...")
pca_result <- mixOmics::pca(X, ncomp = N_COMP, center = TRUE, scale = TRUE)

# --- 6. Variance Verification ---
plsda_var <- get_prop_var(plsda_result)
pca_var   <- get_prop_var(pca_result)

message("--- VARIANCE CHECK ---")
message(sprintf("PLS-DA PC1: %.2f%%", plsda_var[1] * 100))
message(sprintf("PCA    PC1: %.2f%%", pca_var[1] * 100))

if (round(plsda_var[1], 4) == round(pca_var[1], 4)) {
  warning("PLS-DA and PCA variances are identical. Verify model inputs.")
} else {
  message("Confirmed: PLS-DA and PCA variances are distinct.")
}

# --- 7. Generate 3D Visualizations ---
# PLS-DA Plot
plsda_df <- data.frame(
  Comp1 = plsda_result$variates$X[, 1],
  Comp2 = plsda_result$variates$X[, 2],
  Comp3 = plsda_result$variates$X[, 3],
  Group = Y, 
  SampleID = rownames(X)
)

fig_plsda <- plot_3d_scatter(
  df = plsda_df, x_col = "Comp1", y_col = "Comp2", z_col = "Comp3",
  group_col = "Group", sample_col = "SampleID", title_text = "3D PLS-DA: Proteomics",
  x_label = sprintf("Comp 1 (%.1f%%)", plsda_var[1] * 100),
  y_label = sprintf("Comp 2 (%.1f%%)", plsda_var[2] * 100),
  z_label = sprintf("Comp 3 (%.1f%%)", plsda_var[3] * 100),
  color_palette = CANCER_COLORS
)

# PCA Plot
pca_df <- data.frame(
  PC1 = pca_result$variates$X[, 1],
  PC2 = pca_result$variates$X[, 2],
  PC3 = pca_result$variates$X[, 3],
  Group = Y, 
  SampleID = rownames(X)
)

fig_pca <- plot_3d_scatter(
  df = pca_df, x_col = "PC1", y_col = "PC2", z_col = "PC3",
  group_col = "Group", sample_col = "SampleID", title_text = "3D PCA: Proteomics",
  x_label = sprintf("PC 1 (%.1f%%)", pca_var[1] * 100),
  y_label = sprintf("PC 2 (%.1f%%)", pca_var[2] * 100),
  z_label = sprintf("PC 3 (%.1f%%)", pca_var[3] * 100),
  color_palette = CANCER_COLORS
)

# --- 8. File Exports ---
date_str <- format(Sys.Date(), "%Y-%m-%d")

export_visualization <- function(fig, prefix) {
  html_path <- file.path(OUTPUT_DIR, sprintf("%s_%s_Interactive.html", date_str, prefix))
  png_path  <- file.path(OUTPUT_DIR, sprintf("%s_%s_Static.png", date_str, prefix))
  
  saveWidget(fig, html_path, selfcontained = TRUE)
  webshot2::webshot(html_path, file = png_path, vwidth = 1200, vheight = 900, delay = 2)
}

message("Exporting plots to '", OUTPUT_DIR, "' directory...")
export_visualization(fig_plsda, "PLSDA_protein")
export_visualization(fig_pca, "PCA_protein")

message("Analysis complete. Output files successfully written to /", OUTPUT_DIR)
