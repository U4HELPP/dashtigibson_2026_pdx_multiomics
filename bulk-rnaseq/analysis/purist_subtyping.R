
# ==============================================================================
# Script: PurIST Subtyping Analysis for PDAC
# Description: Isolates 'VCU-PC' samples from bulk RNA-seq TPM data, aligns
#              gene symbols, and runs the PurIST single-sample classifier.
# ==============================================================================

# --- 1. Dependencies & Setup ---
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})

# --- 2. Configuration & Relative Paths ---
EXPR_FILE      <- "25.11.20_HUMAN_U54_AllSamples_nf-JaxPDXNet_RSEM_Gene_TPM_U54manuscript.csv"
MODEL_FILE     <- "fitteds_public_2019-02-12.Rdata"
FUNCTIONS_FILE <- "functions.R"
OUTPUT_DIR     <- "results"
SAMPLE_PREFIX  <- "^VCU-PC"

# Ensure output directory exists
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# --- 3. Data & Model Validation ---
message("Validating required model and data files...")
if (!file.exists(MODEL_FILE))     stop(paste("Model file missing:", MODEL_FILE))
if (!file.exists(FUNCTIONS_FILE)) stop(paste("Functions file missing:", FUNCTIONS_FILE))
if (!file.exists(EXPR_FILE))      stop(paste("Expression file missing:", EXPR_FILE))

# Load model classifier and helper functions
load(MODEL_FILE)
source(FUNCTIONS_FILE)
classifier <- classifs[[1]]

# --- 4. Load & Preprocess Expression Data ---
message("Loading RNA-seq expression matrix...")
rna_raw <- read.csv(EXPR_FILE, check.names = FALSE, stringsAsFactors = FALSE)

if (!"geneSym" %in% colnames(rna_raw)) {
  stop("The required column 'geneSym' was not found in the expression file.")
}

# Format gene symbols and average duplicate gene rows
rna_processed <- rna_raw %>%
  filter(!is.na(geneSym) & geneSym != "") %>%
  group_by(geneSym) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  ungroup() %>%
  column_to_rownames("geneSym")

# Isolate target sample columns matching the 'VCU-PC' prefix
vcu_pc_cols <- grep(SAMPLE_PREFIX, colnames(rna_processed), value = TRUE)

if (length(vcu_pc_cols) == 0) {
  stop(sprintf("No sample columns matching prefix '%s' were found.", SAMPLE_PREFIX))
}

data_vcu_pc <- rna_processed[, vcu_pc_cols, drop = FALSE]
message(sprintf("Successfully extracted %d target samples matching '%s'.", ncol(data_vcu_pc), SAMPLE_PREFIX))

# --- 5. Run PurIST Classifier & Export ---
message("Running PurIST classifier...")
predictions <- apply_classifier(data = data_vcu_pc, classifier = classifier)
predictions_clean <- na.omit(predictions)

# Save output with date prefix into results directory
date_str    <- format(Sys.Date(), "%Y-%m-%d")
output_file <- file.path(OUTPUT_DIR, sprintf("%s_VCU_PC_PurIST_predictions.csv", date_str))

write.csv(predictions_clean, output_file, row.names = TRUE)
message(sprintf("Analysis complete. PurIST predictions saved to '%s'", output_file))
