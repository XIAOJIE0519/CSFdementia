#!/usr/bin/env Rscript
# Meta-Analysis Step 1: Individual Study Differential Expression Analysis
# Performs differential expression analysis for each study separately

cat("================================================================================\n")
cat("META-ANALYSIS STEP 1: INDIVIDUAL STUDY DIFFERENTIAL EXPRESSION\n")
cat("================================================================================\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# Create output directory
output_dir <- "F:/1a-EOD-CSF-protein/meta"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read data
cat("Reading combined expression data...\n")
data <- read.csv("F:/1a-EOD-CSF-protein/combine/combined_expression_matrices.csv", 
                 stringsAsFactors = FALSE, check.names = FALSE)

cat("  Total samples:", nrow(data), "\n")
cat("  Total variables:", ncol(data) - 17, "\n\n")

# Define comparison groups
# EOD includes: EOAD, EODLB, EODSD, EOFTD, EOOD
# LOD includes: LOAD, LODLB, LOFTD, LOOD
comparisons <- list(
  EOAD_vs_CN = list(group1 = "EOAD", group2 = "CN"),
  EOAD_vs_LOAD = list(group1 = "EOAD", group2 = "LOAD"),
  LOAD_vs_CN = list(group1 = "LOAD", group2 = "CN"),
  EOFTD_vs_CN = list(group1 = "EOFTD", group2 = "CN"),
  EOFTD_vs_LOFTD = list(group1 = "EOFTD", group2 = "LOFTD"),
  LOFTD_vs_CN = list(group1 = "LOFTD", group2 = "CN"),
  EODSD_vs_CN = list(group1 = "EODSD", group2 = "CN"),
  EOOD_vs_LOOD = list(group1 = "EOOD", group2 = "LOOD"),
  EOOD_vs_CN = list(group1 = "EOOD", group2 = "CN"),
  LOOD_vs_CN = list(group1 = "LOOD", group2 = "CN"),
  EODLB_vs_CN = list(group1 = "EODLB", group2 = "CN"),
  EODLB_vs_LODLB = list(group1 = "EODLB", group2 = "LODLB"),
  LODLB_vs_CN = list(group1 = "LODLB", group2 = "CN"),
  EOD_vs_CN = list(group1 = c("EOAD", "EODLB", "EODSD", "EOFTD", "EOOD"), group2 = "CN"),
  LOD_vs_CN = list(group1 = c("LOAD", "LODLB", "LOFTD", "LOOD"), group2 = "CN"),
  EOD_vs_LOD = list(group1 = c("EOAD", "EODLB", "EODSD", "EOFTD", "EOOD"), 
                    group2 = c("LOAD", "LODLB", "LOFTD", "LOOD"))
)

# Function to perform differential expression analysis
perform_de_analysis <- function(study_data, comparison_name, group1, group2, study_name) {
  
  # Filter samples for the two groups
  group1_samples <- study_data[study_data$Diagnosis_Derived %in% group1, ]
  group2_samples <- study_data[study_data$Diagnosis_Derived %in% group2, ]
  
  # Check if both groups have samples
  if (nrow(group1_samples) == 0 || nrow(group2_samples) == 0) {
    return(NULL)
  }
  
  # Extract protein data (columns 18 onwards, after Diagnosis_Derived)
  protein_cols <- 18:ncol(study_data)
  protein_names <- colnames(study_data)[protein_cols]
  
  results <- data.frame()
  
  for (i in seq_along(protein_cols)) {
    protein_name <- protein_names[i]
    col_idx <- protein_cols[i]
    
    # Get values for both groups
    group1_values <- as.numeric(group1_samples[, col_idx])
    group2_values <- as.numeric(group2_samples[, col_idx])
    
    # Remove NA values
    group1_values <- group1_values[!is.na(group1_values)]
    group2_values <- group2_values[!is.na(group2_values)]
    
    # Check if we have enough samples
    n_total <- length(group1_values) + length(group2_values)
    if (length(group1_values) < 2 || length(group2_values) < 2 || n_total < 5) {
      next
    }
    
    # Calculate log2 fold change
    mean1 <- mean(group1_values)
    mean2 <- mean(group2_values)
    log2fc <- mean1 - mean2
    
    # Perform t-test
    tryCatch({
      t_test <- t.test(group1_values, group2_values, var.equal = FALSE)
      
      # Calculate standard error
      se <- abs(log2fc / t_test$statistic)
      
      # Calculate F-value (for ANOVA-like interpretation)
      f_value <- t_test$statistic^2
      
      # Store results
      results <- rbind(results, data.frame(
        Protein = protein_name,
        Log2FC = log2fc,
        SE = se,
        F_value = f_value,
        P_value = t_test$p.value,
        Comparison = comparison_name,
        N = n_total,
        Study = study_name,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      # Skip proteins that cause errors
    })
  }
  
  # Calculate FDR if we have results
  if (nrow(results) > 0) {
    results$FDR <- p.adjust(results$P_value, method = "BH")
    
    # Reorder columns
    results <- results[, c("Protein", "Log2FC", "SE", "F_value", "P_value", "FDR", "Comparison", "N", "Study")]
  }
  
  return(results)
}

# Get unique studies
studies <- unique(data$Study)
cat("Found", length(studies), "studies:\n")
print(studies)
cat("\n")

# Perform analysis for each study and comparison
total_analyses <- 0
successful_analyses <- 0

for (study in studies) {
  cat("Processing", study, "...\n")
  
  # Filter data for this study
  study_data <- data[data$Study == study, ]
  
  # Check available diagnoses
  available_diagnoses <- unique(study_data$Diagnosis_Derived)
  cat("  Available diagnoses:", paste(available_diagnoses, collapse = ", "), "\n")
  
  for (comp_name in names(comparisons)) {
    total_analyses <- total_analyses + 1
    
    comp <- comparisons[[comp_name]]
    group1 <- comp$group1
    group2 <- comp$group2
    
    # Check if both groups are available
    has_group1 <- any(group1 %in% available_diagnoses)
    has_group2 <- any(group2 %in% available_diagnoses)
    
    if (!has_group1 || !has_group2) {
      cat("    ", comp_name, ": Skipped (missing groups)\n")
      next
    }
    
    # Perform analysis
    results <- perform_de_analysis(study_data, comp_name, group1, group2, study)
    
    if (is.null(results) || nrow(results) == 0) {
      cat("    ", comp_name, ": No results\n")
      next
    }
    
    # Save results
    output_file <- file.path(output_dir, paste0(study, "_", comp_name, ".csv"))
    write.csv(results, output_file, row.names = FALSE)
    
    successful_analyses <- successful_analyses + 1
    cat("    ", comp_name, ": ", nrow(results), " proteins analyzed -> ", 
        basename(output_file), "\n")
  }
  
  cat("\n")
}

cat("================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n")
cat("Total analyses attempted:", total_analyses, "\n")
cat("Successful analyses:", successful_analyses, "\n")
cat("Output directory:", output_dir, "\n")
cat("\n")

# List all output files
output_files <- list.files(output_dir, pattern = "^study_.*\\.csv$", full.names = FALSE)
cat("Generated", length(output_files), "result files\n")

cat("\n================================================================================\n")
cat("STEP 1 COMPLETE\n")
cat("================================================================================\n")
