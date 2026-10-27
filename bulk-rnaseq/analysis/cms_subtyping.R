
# ==============================================================================
# Script: Consensus Molecular Subtype (CMS) Classification for Colorectal Cancer
# Description: Isolates 'VCU-CO' samples from bulk RNA-seq TPM data, maps gene 
#              symbols to Entrez IDs, log2-transforms expression values, and 
#              classifies samples using the CMS Random Forest (RF) predictor.
# ==============================================================================

# --- 1. Dependencies & Setup ---
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  if (!requireNamespace("CMSclassifier", quietly = TRUE)) {
    devtools::install_github("Sage-Bionetworks/CMSclassifier")
  }
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install("org.Hs.eg.db")
  }
  
  library(CMSclassifier)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(tibble)
})

OUTPUT_DIR <- "results"
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# --- 2. Load & Preprocess Expression Data ---
file_path <- "25.11.20_HUMAN_U54_AllSamples_nf-JaxPDXNet_RSEM_Gene_TPM_U54manuscript.csv"

if (!file.exists(file_path)) {
  stop(paste("Expression file not found:", file_path))
}

cat("Loading TPM expression data from:", file_path, "\n")
df <- read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)

# Filter out missing/blank Gene Symbols
df <- df[!is.na(df$geneSym) & df$geneSym != "NA" & df$geneSym != "", ]

# Isolate sample columns matching the 'VCU-CO' prefix
vcu_co_cols <- grep("^VCU-CO", colnames(df), value = TRUE)

if (length(vcu_co_cols) == 0) {
  stop("No sample columns starting with 'VCU-CO' were found in the expression matrix.")
}

cat(sprintf("Found %d samples matching prefix 'VCU-CO'.\n", length(vcu_co_cols)))

# Aggregate duplicate Gene Symbols by mean across VCU-CO sample columns
df_agg <- df %>%
  group_by(geneSym) %>%
  summarise(across(all_of(vcu_co_cols), ~ mean(as.numeric(.x), na.rm = TRUE))) %>%
  ungroup()

# --- 3. Map Gene Symbols to Entrez IDs ---
cat("Mapping gene symbols to Entrez IDs...\n")
entrez_mapping <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys = df_agg$geneSym, 
  columns = "ENTREZID", 
  keytype = "SYMBOL"
)

# Clean mapping: remove unmapped genes and duplicate symbols
entrez_mapping <- entrez_mapping[!is.na(entrez_mapping$ENTREZID), ]
entrez_mapping <- entrez_mapping[!duplicated(entrez_mapping$SYMBOL), ]

# Merge Entrez IDs back to expression data
df_merged <- merge(df_agg, entrez_mapping, by.x = "geneSym", by.y = "SYMBOL")

# Deduplicate Entrez IDs if multiple symbols mapped to the same Entrez ID
df_unique <- df_merged[!duplicated(df_merged$ENTREZID), ]

# --- 4. Prepare Matrix & Log2 Transform ---
# Extract numeric matrix and force explicit numeric storage mode
expr_matrix <- as.matrix(df_unique[, vcu_co_cols])
storage.mode(expr_matrix) <- "numeric"
rownames(expr_matrix) <- as.character(df_unique$ENTREZID)

# Apply log2 transformation: log2(TPM + 1)
cat("Applying log2(x + 1) transformation to expression matrix...\n")
expr_matrix <- log2(expr_matrix + 1)

# --- 5. Run CMS Random Forest Predictor ---
cat("Running CMSclassifier RF method...\n")
cms_rf_results <- classifyCMS.RF(expr_matrix)

# --- 6. Export Results ---
date_str <- format(Sys.Date(), "%Y-%m-%d")
rf_file  <- file.path(OUTPUT_DIR, sprintf("%s_VCU_CO_CMS_RF_Results.csv", date_str))

write.csv(cms_rf_results, rf_file, row.names = TRUE)
cat(sprintf("\nRF results successfully saved to '%s'\n", rf_file))

# --- 7. Console Summary ---
cat("\nSummary of RF Subtypes (VCU-CO Cohort):\n")
print(table(cms_rf_results$RF.predicted, useNA = "ifany"))