
# ==============================================================================
# R Script: Breast Cancer Drug Response Comparison (EUR vs AFR Ancestry)
# Exports analysis tables and volcano plots
# Author: Nina Dashti-Gibson
# ==============================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, dplyr, tidyr, ggplot2, ggrepel, svglite, patchwork)

# ------------------------------------------------------------------------------
# 0. Setup Output Directory and Date Prefix
# ------------------------------------------------------------------------------
results_dir <- "results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

date_stamp <- format(Sys.Date(), "%Y%m%d")

# Helper function to generate output paths: results/YYYYMMDD_filename.ext
get_output_path <- function(filename) {
  file.path(results_dir, paste0(date_stamp, "_", filename))
}

# ------------------------------------------------------------------------------
# 1. Load and Filter Metadata
# ------------------------------------------------------------------------------
metadata_file <- "sample_metadata.csv"
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

bc_metadata <- metadata %>%
  filter(startsWith(PDX, "VCU-BC") & Ancestry %in% c("EUR", "AFR"))

eur_samples <- bc_metadata %>% filter(Ancestry == "EUR") %>% pull(PDX)
afr_samples <- bc_metadata %>% filter(Ancestry == "AFR") %>% pull(PDX)

cat("--- Sample Summary ---\n")
cat("Metadata File:", metadata_file, "\n")
print(table(bc_metadata$Ancestry))
cat("\n")

# ------------------------------------------------------------------------------
# 2. Differential Drug Response Analysis
# ------------------------------------------------------------------------------
analyze_drug_response <- function(excel_path, conc_label, bc_meta) {
  if (!file.exists(excel_path)) {
    warning(paste("File not found:", excel_path))
    return(NULL)
  }
  
  cat(paste0("Processing ", conc_label, " dataset: ", excel_path, "\n"))
  df <- read_excel(excel_path)
  
  eur_cols <- intersect(eur_samples, colnames(df))
  afr_cols <- intersect(afr_samples, colnames(df))
  
  long_df <- df %>%
    pivot_longer(
      cols = dplyr::all_of(c(eur_cols, afr_cols)),
      names_to = "PDX",
      values_to = "Inhibition"
    ) %>%
    left_join(bc_meta %>% select(PDX, Ancestry), by = "PDX") %>%
    filter(!is.na(Inhibition))
  
  results <- long_df %>%
    group_by(Drug, Class) %>%
    group_modify(~ {
      sub_data <- .x
      eur_vals <- sub_data$Inhibition[sub_data$Ancestry == "EUR"]
      afr_vals <- sub_data$Inhibition[sub_data$Ancestry == "AFR"]
      
      mean_eur <- mean(eur_vals, na.rm = TRUE)
      mean_afr <- mean(afr_vals, na.rm = TRUE)
      diff_afr_eur <- mean_afr - mean_eur
      
      p_ttest  <- NA_real_
      p_wilcox <- NA_real_
      
      if (length(eur_vals) >= 2 && length(afr_vals) >= 2) {
        p_ttest  <- tryCatch(t.test(Inhibition ~ Ancestry, data = sub_data, var.equal = FALSE)$p.value, error = function(e) NA)
        p_wilcox <- tryCatch(wilcox.test(Inhibition ~ Ancestry, data = sub_data)$p.value, error = function(e) NA)
      }
      
      tibble(
        N_EUR = length(eur_vals),
        N_AFR = length(afr_vals),
        Mean_Inhibition_EUR = mean_eur,
        Mean_Inhibition_AFR = mean_afr,
        Diff_AFR_minus_EUR = diff_afr_eur,
        p_value_ttest = p_ttest,
        p_value_wilcox = p_wilcox
      )
    }) %>%
    ungroup() %>%
    mutate(
      Concentration = conc_label,
      adj_p_value_FDR = p.adjust(p_value_ttest, method = "BH")
    ) %>%
    arrange(p_value_ttest)
  
  return(results)
}

# ------------------------------------------------------------------------------
# 3. Run Analysis & Export CSV Results
# ------------------------------------------------------------------------------
file_1uM  <- "20251113_All_U54_1uM_inhibition.xlsx"
file_10uM <- "20251113_All_U54_10uM_inhibition.xlsx"

res_1uM  <- analyze_drug_response(file_1uM, "1uM", bc_metadata)
res_10uM <- analyze_drug_response(file_10uM, "10uM", bc_metadata)

combined_results <- bind_rows(res_1uM, res_10uM)

if (!is.null(res_1uM)) {
  write.csv(res_1uM, get_output_path("breast_cancer_ancestry_comparison_1uM.csv"), row.names = FALSE)
}
if (!is.null(res_10uM)) {
  write.csv(res_10uM, get_output_path("breast_cancer_ancestry_comparison_10uM.csv"), row.names = FALSE)
}
if (!is.null(combined_results)) {
  write.csv(combined_results, get_output_path("breast_cancer_ancestry_comparison_combined.csv"), row.names = FALSE)
}

# ------------------------------------------------------------------------------
# 4. Volcano Plot Construction
# ------------------------------------------------------------------------------
make_volcano_plot <- function(res_df, min_effect_size = 20, p_cutoff = 0.05, top_n_labels = 10) {
  if (is.null(res_df)) return(NULL)
  
  df_plot <- res_df %>%
    filter(!is.na(p_value_ttest)) %>%
    mutate(
      neg_log10_p = -log10(p_value_ttest),
      Category = case_when(
        p_value_ttest < p_cutoff & abs(Diff_AFR_minus_EUR) >= min_effect_size ~ paste0("p < ", p_cutoff, " & |Diff| >= ", min_effect_size, "%"),
        p_value_ttest < p_cutoff ~ paste0("p < ", p_cutoff, " only"),
        TRUE ~ "Not Significant"
      )
    )
  
  hits_to_label <- df_plot %>%
    filter(p_value_ttest < p_cutoff & abs(Diff_AFR_minus_EUR) >= min_effect_size) %>%
    slice_max(order_by = abs(Diff_AFR_minus_EUR), n = top_n_labels)
  
  if (nrow(hits_to_label) == 0) {
    hits_to_label <- df_plot %>%
      filter(p_value_ttest < p_cutoff) %>%
      slice_max(order_by = abs(Diff_AFR_minus_EUR), n = top_n_labels)
  }
  
  p <- ggplot(df_plot, aes(x = Diff_AFR_minus_EUR, y = neg_log10_p)) +
    geom_point(aes(color = Category), alpha = 0.8, size = 3.0) +
    scale_color_manual(values = setNames(
      c("#DC143C", "orange", "gray60"),
      c(paste0("p < ", p_cutoff, " & |Diff| >= ", min_effect_size, "%"),
        paste0("p < ", p_cutoff, " only"),
        "Not Significant")
    )) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_vline(xintercept = c(-min_effect_size, min_effect_size), linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray30", linewidth = 0.6) +
    geom_text_repel(
      data = hits_to_label,
      aes(label = Drug),
      fontface = "bold",
      size = 4.8,
      max.overlaps = 25,
      box.padding = 0.5
    ) +
    labs(
      x = "Mean Inhibition Difference: AFR - EUR (%)",
      y = "-log10(p-value)",
      color = "Threshold"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      axis.title = element_text(size = 16, face = "bold"),
      axis.text  = element_text(size = 14, color = "black"),
      legend.title = element_text(size = 15, face = "bold"),
      legend.text  = element_text(size = 14),
      legend.position = "bottom"
    )
  
  return(p)
}

p_1uM  <- make_volcano_plot(res_1uM, min_effect_size = 20, p_cutoff = 0.05)
p_10uM <- make_volcano_plot(res_10uM, min_effect_size = 20, p_cutoff = 0.05)

# ------------------------------------------------------------------------------
# 5. Figure Composition and Export
# ------------------------------------------------------------------------------
if (!is.null(p_1uM) && !is.null(p_10uM)) {
  side_by_side_plot <- (p_1uM | p_10uM) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = NULL,
      subtitle = NULL,
      tag_levels = 'A'
    ) &
    theme(
      legend.position = "bottom",
      plot.tag = element_text(size = 20, face = "bold")
    )
  
  png_out <- get_output_path("volcano_plots_side_by_side.png")
  svg_out <- get_output_path("volcano_plots_side_by_side.svg")
  
  ggsave(png_out, plot = side_by_side_plot, width = 10, height = 4.5, dpi = 300)
  ggsave(svg_out, plot = side_by_side_plot, width = 10, height = 4.5, device = "svg")
  
  cat("Successfully generated and saved results to '", results_dir, "/':\n", sep = "")
  cat(" - ", basename(png_out), "\n", sep = "")
  cat(" - ", basename(svg_out), "\n", sep = "")
  
  print(side_by_side_plot)
}
