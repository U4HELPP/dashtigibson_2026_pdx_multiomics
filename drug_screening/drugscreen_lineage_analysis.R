# ==============================================================================
# Pipeline: ANOVA, Tukey's Post-Hoc, Box Plots, and Z-Score Heatmaps
# Description: Evaluates drug response and drug-class efficacy across cancer
#              lineages/types using screen inhibition metrics.
# Author: Nina Dashti-Gibson
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Package Dependencies
# ------------------------------------------------------------------------------
required_packages <- c(
  "readxl", "readr", "dplyr", "tidyr", "purrr", "broom", 
  "ggplot2", "ggsignif", "svglite", "stringr", "tibble", 
  "grid", "ComplexHeatmap", "circlize"
)

missing_pkgs <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(paste("Missing required packages:", paste(missing_pkgs, collapse = ", "), 
             "\nPlease install them prior to running this script."))
}

library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)
library(ggsignif)
library(stringr)
library(tibble)
library(grid)
library(ComplexHeatmap)
library(circlize)

# ------------------------------------------------------------------------------
# 2. Pipeline Configuration
# ------------------------------------------------------------------------------
CONFIG <- list(
  output_dir       = "results",
  sample_meta_path = "sample_metadata.csv",
  fdr_threshold    = 0.05,
  inputs           = list(
    "1uM"  = "20251113_All_U54_1uM_inhibition.xlsx",
    "10uM" = "20251113_All_U54_10uM_inhibition.xlsx"
  ),
  cancer_colors    = c(
    "breast"     = "#cb027e",
    "colorectal" = "#0081cb",
    "pancreatic" = "#face20"
  )
)

today_str <- format(Sys.Date(), "%Y-%m-%d")
if (!dir.exists(CONFIG$output_dir)) dir.create(CONFIG$output_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 3. Load & Validate Metadata
# ------------------------------------------------------------------------------
if (!file.exists(CONFIG$sample_meta_path)) {
  stop(paste("Metadata file not found at:", CONFIG$sample_meta_path))
}

metadata <- read_csv(CONFIG$sample_meta_path, show_col_types = FALSE) %>%
  select(PDX, CancerType) %>%
  filter(!is.na(PDX) & !is.na(CancerType)) %>%
  mutate(CancerType = tolower(CancerType))

valid_samples <- metadata$PDX

# ------------------------------------------------------------------------------
# 4. Pipeline Execution Function
# ------------------------------------------------------------------------------
#' Run Drug Screen Analysis Pipeline
#'
#' @param data_path String. Path to input Excel data file.
#' @param concentration_label String. Label for the dataset run (e.g., "1uM").
#' @return List containing ANOVA and Tukey post-hoc results for drugs and classes.
run_drug_anova_pipeline <- function(data_path, concentration_label) {
  
  message("\n--------------------------------------------------")
  message(" Processing: ", concentration_label)
  message("--------------------------------------------------")
  
  if (!file.exists(data_path)) {
    warning(paste("File not found:", data_path, "- skipping analysis."))
    return(NULL)
  }
  
  # A. Read and Format Raw Data
  raw_df <- read_excel(data_path)
  colnames(raw_df)[1:2] <- c("Drug", "Class")
  
  long_df <- raw_df %>%
    pivot_longer(
      cols = -c(Drug, Class),
      names_to = "PDX",
      values_to = "Inhibition"
    ) %>%
    filter(PDX %in% valid_samples) %>%
    inner_join(metadata, by = "PDX") %>%
    mutate(
      CancerType = as.factor(CancerType),
      Inhibition = as.numeric(Inhibition)
    )
  
  # B. Individual Drug Level Statistics (ANOVA + Tukey)
  drug_fits <- long_df %>%
    group_by(Drug, Class) %>%
    nest() %>%
    mutate(
      valid = map_lgl(data, ~ length(unique(.x$CancerType[!is.na(.x$Inhibition)])) > 1),
      fit = map2(data, valid, function(d, v) {
        if (!v) return(NULL)
        tryCatch(aov(Inhibition ~ CancerType, data = d), error = function(e) NULL)
      }),
      anova_tidy = map(fit, ~ if (!is.null(.x)) tidy(.x) else NULL)
    )
  
  drug_anova_summary <- drug_fits %>%
    filter(!map_lgl(anova_tidy, is.null)) %>%
    unnest(anova_tidy) %>%
    filter(term == "CancerType") %>%
    ungroup() %>%
    select(Drug, Class, df, sumsq, meansq, statistic, p.value) %>%
    rename(F_statistic = statistic, ANOVA_p_value = p.value) %>%
    mutate(FDR_q_value = p.adjust(ANOVA_p_value, method = "BH")) %>%
    arrange(ANOVA_p_value)
  
  significant_drugs <- drug_anova_summary %>%
    filter(FDR_q_value < CONFIG$fdr_threshold) %>%
    select(Drug, Class)
  
  drug_tukey_summary <- drug_fits %>%
    inner_join(significant_drugs, by = c("Drug", "Class")) %>%
    mutate(
      tukey_tidy = map(fit, function(f) {
        if (is.null(f)) return(NULL)
        tryCatch({
          tukey_res <- TukeyHSD(f)
          as.data.frame(tukey_res$CancerType) %>%
            tibble::rownames_to_column(var = "Comparison")
        }, error = function(e) NULL)
      })
    ) %>%
    filter(!map_lgl(tukey_tidy, is.null)) %>%
    unnest(tukey_tidy) %>%
    ungroup()
  
  if (nrow(drug_tukey_summary) > 0) {
    drug_tukey_summary <- drug_tukey_summary %>%
      select(Drug, Class, Comparison, diff, lwr, upr, `p adj`) %>%
      rename(Tukey_p_adj = `p adj`, Difference = diff, CI_lower = lwr, CI_upper = upr) %>%
      arrange(Tukey_p_adj)
  } else {
    drug_tukey_summary <- tibble(
      Drug = character(), Class = character(), Comparison = character(),
      Difference = numeric(), CI_lower = numeric(), CI_upper = numeric(), Tukey_p_adj = numeric()
    )
  }
  
  # C. Heatmap Generation
  if (nrow(significant_drugs) > 0) {
    message(paste0("Generating Z-score heatmap for ", nrow(significant_drugs), " significant drug(s)..."))
    
    sig_drug_data <- long_df %>%
      inner_join(significant_drugs, by = c("Drug", "Class"))
    
    heat_mat <- sig_drug_data %>%
      select(Drug, PDX, Inhibition) %>%
      pivot_wider(names_from = PDX, values_from = Inhibition) %>%
      column_to_rownames("Drug") %>%
      as.matrix()
    
    sample_annotation <- metadata %>%
      arrange(CancerType, PDX) %>%
      column_to_rownames("PDX")
    
    heat_mat <- heat_mat[, rownames(sample_annotation), drop = FALSE]
    
    # Calculate Z-scores with zero-variance protection
    heat_mat_z <- t(scale(t(heat_mat)))
    heat_mat_z[is.na(heat_mat_z)] <- 0
    
    col_fun <- colorRamp2(
      c(-4, -2, 0, 2, 4),
      c("#021480", "#8267b6", "#ffffff", "#e66f5f", "#c60000")
    )
    
    ha_top <- HeatmapAnnotation(
      `Cancer type` = sample_annotation$CancerType,
      col = list(`Cancer type` = CONFIG$cancer_colors),
      show_annotation_name = FALSE,
      show_legend = TRUE,
      annotation_legend_param = list(
        `Cancer type` = list(
          title = "Cancer type",
          title_gp = gpar(fontsize = 11, fontface = "bold"),
          labels_gp = gpar(fontsize = 10),
          grid_height = unit(4, "mm"),
          grid_width = unit(4, "mm")
        )
      )
    )
    
    ht <- Heatmap(
      heat_mat_z,
      name = "z-score",
      col = col_fun,
      top_annotation = ha_top,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_row_dend = TRUE,
      show_column_dend = FALSE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = gpar(fontsize = 10),
      column_names_gp = gpar(fontsize = 9),
      column_names_rot = 45,
      rect_gp = gpar(col = NA),
      heatmap_legend_param = list(
        title = "z-score",
        direction = "horizontal",
        title_gp = gpar(fontsize = 11, fontface = "bold"),
        labels_gp = gpar(fontsize = 9),
        at = c(-4, -2, 0, 2, 4),
        legend_width = unit(4, "cm")
      )
    )
    
    out_png <- file.path(CONFIG$output_dir, paste0(today_str, "_Heatmap_Significant_Drugs_", concentration_label, ".png"))
    out_svg <- file.path(CONFIG$output_dir, paste0(today_str, "_Heatmap_Significant_Drugs_", concentration_label, ".svg"))
    
    png(out_png, width = 7, height = 6, units = "in", res = 300)
    draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = FALSE)
    dev.off()
    
    svg(out_svg, width = 7, height = 6)
    draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = FALSE)
    dev.off()
  } else {
    message(paste0("No individual drugs met FDR < ", CONFIG$fdr_threshold, " at ", concentration_label, "."))
  }
  
  # D. Drug Class Level Statistics (ANOVA + Tukey)
  class_aggregated_df <- long_df %>%
    group_by(Class, PDX, CancerType) %>%
    summarise(Mean_Inhibition = mean(Inhibition, na.rm = TRUE), .groups = "drop")
  
  class_fits <- class_aggregated_df %>%
    group_by(Class) %>%
    nest() %>%
    mutate(
      valid = map_lgl(data, ~ length(unique(.x$CancerType[!is.na(.x$Mean_Inhibition)])) > 1),
      fit = map2(data, valid, function(d, v) {
        if (!v) return(NULL)
        tryCatch(aov(Mean_Inhibition ~ CancerType, data = d), error = function(e) NULL)
      }),
      anova_tidy = map(fit, ~ if (!is.null(.x)) tidy(.x) else NULL)
    )
  
  class_anova_summary <- class_fits %>%
    filter(!map_lgl(anova_tidy, is.null)) %>%
    unnest(anova_tidy) %>%
    filter(term == "CancerType") %>%
    ungroup() %>%
    select(Class, df, sumsq, meansq, statistic, p.value) %>%
    rename(F_statistic = statistic, ANOVA_p_value = p.value) %>%
    mutate(FDR_q_value = p.adjust(ANOVA_p_value, method = "BH")) %>%
    arrange(ANOVA_p_value)
  
  significant_classes <- class_anova_summary %>%
    filter(FDR_q_value < CONFIG$fdr_threshold) %>%
    select(Class, ANOVA_p_value, FDR_q_value)
  
  class_tukey_summary <- class_fits %>%
    inner_join(significant_classes, by = "Class") %>%
    mutate(
      tukey_tidy = map(fit, function(f) {
        if (is.null(f)) return(NULL)
        tryCatch({
          tukey_res <- TukeyHSD(f)
          as.data.frame(tukey_res$CancerType) %>%
            tibble::rownames_to_column(var = "Comparison")
        }, error = function(e) NULL)
      })
    ) %>%
    filter(!map_lgl(tukey_tidy, is.null)) %>%
    unnest(tukey_tidy) %>%
    ungroup()
  
  if (nrow(class_tukey_summary) > 0) {
    class_tukey_summary <- class_tukey_summary %>%
      select(Class, Comparison, diff, lwr, upr, `p adj`) %>%
      rename(Tukey_p_adj = `p adj`, Difference = diff, CI_lower = lwr, CI_upper = upr) %>%
      arrange(Tukey_p_adj)
  } else {
    class_tukey_summary <- tibble(
      Class = character(), Comparison = character(),
      Difference = numeric(), CI_lower = numeric(), CI_upper = numeric(), Tukey_p_adj = numeric()
    )
  }
  
  # E. Class Box Plot Visualizations
  if (nrow(significant_classes) > 0) {
    message(paste0("Generating box plots for ", nrow(significant_classes), " significant drug class(es)..."))
    
    sig_class_data <- class_aggregated_df %>%
      inner_join(significant_classes, by = "Class") %>%
      mutate(CancerType = str_to_title(CancerType))
    
    tukey_facet_df <- class_tukey_summary %>%
      filter(Tukey_p_adj < 0.05) %>%
      separate(Comparison, into = c("group1", "group2"), sep = "-") %>%
      mutate(
        group1 = str_to_title(group1),
        group2 = str_to_title(group2),
        stars = case_when(
          Tukey_p_adj <= 0.001 ~ "***",
          Tukey_p_adj <= 0.01  ~ "**",
          Tukey_p_adj <= 0.05  ~ "*",
          TRUE                 ~ "ns"
        )
      ) %>%
      inner_join(
        sig_class_data %>%
          group_by(Class) %>%
          summarise(
            y_max = max(Mean_Inhibition, na.rm = TRUE),
            y_min = min(Mean_Inhibition, na.rm = TRUE),
            .groups = "drop"
          ),
        by = "Class"
      ) %>%
      group_by(Class) %>%
      mutate(
        y_range = y_max - y_min,
        step = if_else(y_range == 0, 5, y_range * 0.16),
        y_position = y_max + step * row_number()
      ) %>%
      ungroup()
    
    n_classes   <- nrow(significant_classes)
    n_rows      <- min(n_classes, 3)
    n_cols      <- ceiling(n_classes / n_rows)
    plot_width  <- max(7, n_cols * 3.8)
    plot_height <- max(8, n_rows * 3.2)
    
    cancer_colors_title <- setNames(
      as.character(CONFIG$cancer_colors), 
      str_to_title(names(CONFIG$cancer_colors))
    )
    
    p_facet <- ggplot(sig_class_data, aes(x = CancerType, y = Mean_Inhibition, fill = CancerType)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black", width = 0.5) +
      geom_jitter(width = 0.15, size = 1.0, alpha = 0.85, color = "black") +
      facet_wrap(~ Class, scales = "free_y", nrow = 3) +
      scale_fill_manual(values = cancer_colors_title) +
      theme_bw(base_size = 15) +
      theme(
        aspect.ratio = 0.85,
        strip.background = element_rect(fill = "#f0f0f0"),
        strip.text = element_text(size = 14, face = "bold"),
        title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 13, face = "plain"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 13, color = "black"),
        legend.position = "none"
      ) +
      labs(
        title = paste0("Drug Class Response across Cancer Types (", concentration_label, ")"),
        subtitle = "Significant classes following ANOVA FDR correction (q < 0.05); stars indicate Tukey HSD p < 0.05",
        x = "Cancer Type",
        y = "Mean Growth Inhibition (%)"
      )
    
    if (nrow(tukey_facet_df) > 0) {
      p_facet <- p_facet +
        geom_signif(
          data = tukey_facet_df,
          aes(
            xmin = group1,
            xmax = group2,
            annotations = stars,
            y_position = y_position
          ),
          manual = TRUE,
          inherit.aes = FALSE,
          tip_length = 0.02,
          textsize = 7,
          vjust = 0.3
        )
    }
    
    ggsave(
      filename = file.path(CONFIG$output_dir, paste0(today_str, "_Significant_Classes_Boxplots_", concentration_label, ".png")),
      plot = p_facet, width = plot_width, height = plot_height, dpi = 300
    )
    ggsave(
      filename = file.path(CONFIG$output_dir, paste0(today_str, "_Significant_Classes_Boxplots_", concentration_label, ".svg")),
      plot = p_facet, width = plot_width, height = plot_height, device = "svg"
    )
    
    for (cls in significant_classes$Class) {
      cls_data  <- sig_class_data %>% filter(Class == cls)
      cls_tukey <- class_tukey_summary %>% filter(Class == cls)
      
      annotations <- cls_tukey %>%
        separate(Comparison, into = c("group1", "group2"), sep = "-") %>%
        mutate(
          group1 = str_to_title(group1),
          group2 = str_to_title(group2),
          stars = case_when(
            Tukey_p_adj <= 0.001 ~ "***",
            Tukey_p_adj <= 0.01  ~ "**",
            Tukey_p_adj <= 0.05  ~ "*",
            TRUE                 ~ "ns"
          )
        )
      
      comp_list   <- lapply(1:nrow(annotations), function(i) c(annotations$group1[i], annotations$group2[i]))
      labels_list <- annotations$stars
      
      p_ind <- ggplot(cls_data, aes(x = CancerType, y = Mean_Inhibition, fill = CancerType)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
        geom_jitter(width = 0.15, size = 1.2, alpha = 0.85, color = "black") +
        scale_fill_manual(values = cancer_colors_title) +
        theme_classic(base_size = 14) +
        theme(
          legend.position = "none",
          plot.title = element_text(face = "bold", size = 14)
        ) +
        labs(
          title = paste0(cls, " (", concentration_label, ")"),
          subtitle = paste0("ANOVA FDR q = ", sprintf("%.3e", significant_classes$FDR_q_value[significant_classes$Class == cls])),
          x = "Cancer Type",
          y = "Mean Growth Inhibition (%)"
        )
      
      if (length(comp_list) > 0) {
        y_max <- max(cls_data$Mean_Inhibition, na.rm = TRUE)
        y_range <- y_max - min(cls_data$Mean_Inhibition, na.rm = TRUE)
        step <- if_else(y_range == 0, 5, y_range * 0.16)
        y_positions <- seq(from = y_max + step, length.out = length(comp_list), by = step)
        
        p_ind <- p_ind +
          geom_signif(
            comparisons = comp_list,
            annotations = labels_list,
            y_position = y_positions,
            tip_length = 0.02,
            textsize = 6,
            vjust = 0.3
          )
      }
      
      safe_cls_name <- gsub("[^A-Za-z0-9]", "_", cls)
      ggsave(
        filename = file.path(CONFIG$output_dir, paste0(today_str, "_Boxplot_Class_", safe_cls_name, "_", concentration_label, ".png")),
        plot = p_ind, width = 6, height = 5.5, dpi = 300
      )
    }
  } else {
    message(paste0("No drug classes met FDR < ", CONFIG$fdr_threshold, " at ", concentration_label, "."))
  }
  
  # F. Save CSV Output Summaries
  write_csv(drug_anova_summary, file.path(CONFIG$output_dir, paste0(today_str, "_ANOVA_Individual_Drugs_", concentration_label, ".csv")))
  write_csv(drug_tukey_summary, file.path(CONFIG$output_dir, paste0(today_str, "_Tukey_Individual_Drugs_", concentration_label, ".csv")))
  write_csv(class_anova_summary, file.path(CONFIG$output_dir, paste0(today_str, "_ANOVA_Drug_Classes_", concentration_label, ".csv")))
  write_csv(class_tukey_summary, file.path(CONFIG$output_dir, paste0(today_str, "_Tukey_Drug_Classes_", concentration_label, ".csv")))
  
  return(list(
    drug_anova = drug_anova_summary,
    drug_tukey = drug_tukey_summary,
    class_anova = class_anova_summary,
    class_tukey = class_tukey_summary
  ))
}

# ------------------------------------------------------------------------------
# 5. Pipeline Orchestration
# ------------------------------------------------------------------------------
results <- map2(CONFIG$inputs, names(CONFIG$inputs), run_drug_anova_pipeline)

message("\nAnalysis Complete! Outputs exported to: '", CONFIG$output_dir, "/'")