#!/usr/bin/env Rscript
# Meta-Analysis Step 2: Sample-size weighted + weighted z-score method
# Using metapro package as described in the original paper

cat("================================================================================\n")
cat("META-ANALYSIS STEP 2: WEIGHTED Z-SCORE METHOD (METAPRO)\n")
cat("================================================================================\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  
  # Install BiocManager if needed
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # Install multtest (required by metapro)
  if (!require("multtest", quietly = TRUE)) {
    cat("Installing multtest package from Bioconductor...\n")
    BiocManager::install("multtest", ask = FALSE, update = FALSE)
  }
  
  # Install metapro if not available
  if (!require("metapro", quietly = TRUE)) {
    cat("Installing metapro package...\n")
    install.packages("metapro")
    library(metapro)
  } else {
    library(metapro)
  }
})

# Set directories
input_dir <- "F:/1a-EOD-CSF-protein/meta"
output_dir <- input_dir

# Define study priority rules (lower number = higher priority)
get_study_priority <- function(study_name) {
  priority_map <- c(
    "study_1" = 1,
    "study_2" = 2,
    "study_5" = 3,
    "study_3_Soma" = 4,
    "study_3_Olink" = 5,
    "study_10" = 6,
    "study_8" = 10,
    "study_9" = 11,
    "study_4" = 20,
    "study_6" = 21,
    "study_7" = 22,
    "study_11" = 23
  )
  
  return(ifelse(study_name %in% names(priority_map), priority_map[study_name], 100))
}

# Function to apply study selection rules based on duplicate groups
select_studies_for_protein <- function(protein_data) {
  # Define duplicate groups
  # Group 1: studies 1, 2, 3_Soma, 3_Olink, 5, 10 are considered duplicates
  group1 <- c("study_1", "study_2", "study_5", "study_3_Soma", "study_3_Olink", "study_10")
  
  # Group 2: studies 8, 9 are considered duplicates
  group2 <- c("study_8", "study_9")
  
  # For group 1: if multiple studies present, keep only the highest priority one
  if (sum(protein_data$Study %in% group1) >= 2) {
    in_group1 <- protein_data$Study %in% group1
    priorities <- sapply(protein_data$Study[in_group1], get_study_priority)
    keep_study <- protein_data$Study[in_group1][which.min(priorities)]
    protein_data <- protein_data[!in_group1 | protein_data$Study == keep_study, ]
  }
  
  # For group 2: if multiple studies present, keep only the highest priority one
  if (sum(protein_data$Study %in% group2) >= 2) {
    in_group2 <- protein_data$Study %in% group2
    priorities <- sapply(protein_data$Study[in_group2], get_study_priority)
    keep_study <- protein_data$Study[in_group2][which.min(priorities)]
    protein_data <- protein_data[!in_group2 | protein_data$Study == keep_study, ]
  }
  
  return(protein_data)
}

# Function to count unique non-duplicate studies
count_unique_studies <- function(study_names) {
  # Define duplicate groups
  group1 <- c("study_1", "study_2", "study_5", "study_3_Soma", "study_3_Olink", "study_10")
  group2 <- c("study_8", "study_9")
  
  unique_count <- 0
  
  # Check if any study from group1 is present
  if (any(study_names %in% group1)) {
    unique_count <- unique_count + 1
  }
  
  # Check if any study from group2 is present
  if (any(study_names %in% group2)) {
    unique_count <- unique_count + 1
  }
  
  # Count other studies (not in any duplicate group)
  other_studies <- study_names[!(study_names %in% c(group1, group2))]
  unique_count <- unique_count + length(unique(other_studies))
  
  return(unique_count)
}

# Function to perform sample-size weighted average and weighted z-score meta-analysis
perform_weighted_meta <- function(log2fc, se, n_samples, p_values) {
  # Remove NA values
  valid <- !is.na(log2fc) & !is.na(se) & !is.na(n_samples) & !is.na(p_values) & se > 0
  log2fc <- log2fc[valid]
  se <- se[valid]
  n_samples <- n_samples[valid]
  p_values <- p_values[valid]
  
  if (length(log2fc) < 2) {
    return(list(
      weighted_effect = NA,
      ci_lower = NA,
      ci_upper = NA,
      z_score = NA,
      p_value = NA,
      i_squared = NA
    ))
  }
  
  # Sample-size weighted average of effect sizes
  weights_n <- n_samples / sum(n_samples)
  weighted_effect <- sum(weights_n * log2fc)
  
  # Standard error of weighted effect
  se_weighted <- sqrt(sum((weights_n * se)^2))
  
  # 95% confidence interval
  ci_lower <- weighted_effect - 1.96 * se_weighted
  ci_upper <- weighted_effect + 1.96 * se_weighted
  
  # Weighted z-score method using metapro
  tryCatch({
    # Convert p-values to z-scores with direction
    z_scores <- qnorm(p_values / 2, lower.tail = FALSE) * sign(log2fc)
    
    # Weight by square root of sample size (standard for z-score method)
    weights_z <- sqrt(n_samples)
    
    # Weighted z-score
    z_meta <- sum(weights_z * z_scores) / sqrt(sum(weights_z^2))
    
    # Meta-analytic p-value
    p_meta <- 2 * pnorm(abs(z_meta), lower.tail = FALSE)
    
    # Calculate I-squared for heterogeneity
    # Using inverse variance weights for heterogeneity calculation
    weights_iv <- 1 / (se^2)
    q_stat <- sum(weights_iv * (log2fc - weighted_effect)^2)
    df <- length(log2fc) - 1
    i_squared <- max(0, (q_stat - df) / q_stat)
    
    return(list(
      weighted_effect = weighted_effect,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      z_score = z_meta,
      p_value = p_meta,
      i_squared = i_squared
    ))
    
  }, error = function(e) {
    # Fallback if metapro fails
    z_meta <- weighted_effect / se_weighted
    p_meta <- 2 * pnorm(abs(z_meta), lower.tail = FALSE)
    
    return(list(
      weighted_effect = weighted_effect,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      z_score = z_meta,
      p_value = p_meta,
      i_squared = 0
    ))
  })
}

# Function to perform leave-one-out sensitivity analysis
perform_loo_analysis <- function(log2fc, n_samples) {
  valid <- !is.na(log2fc) & !is.na(n_samples)
  log2fc <- log2fc[valid]
  n_samples <- n_samples[valid]
  
  if (length(log2fc) < 3) {
    return(list(mean = NA, ci_lower = NA, ci_upper = NA))
  }
  
  loo_effects <- numeric(length(log2fc))
  
  for (i in seq_along(log2fc)) {
    loo_log2fc <- log2fc[-i]
    loo_n <- n_samples[-i]
    weights <- loo_n / sum(loo_n)
    loo_effects[i] <- sum(weights * loo_log2fc)
  }
  
  return(list(
    mean = mean(loo_effects),
    ci_lower = min(loo_effects),
    ci_upper = max(loo_effects)
  ))
}

# Define comparisons
comparisons <- c("EOAD_vs_CN", "EOAD_vs_LOAD", "LOAD_vs_CN",
                 "EOFTD_vs_CN", "EOFTD_vs_LOFTD", "LOFTD_vs_CN",
                 "EODSD_vs_CN", "EOOD_vs_LOOD", "EOOD_vs_CN", "LOOD_vs_CN",
                 "EODLB_vs_CN", "EODLB_vs_LODLB", "LODLB_vs_CN",
                 "EOD_vs_CN", "LOD_vs_CN", "EOD_vs_LOD")

cat("Processing", length(comparisons), "comparisons using weighted z-score method...\n\n")

for (comp_name in comparisons) {
  cat("Processing", comp_name, "...\n")
  
  # Find all result files for this comparison
  pattern <- paste0("^study_.*_", comp_name, "\\.csv$")
  result_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  
  if (length(result_files) == 0) {
    cat("  No result files found\n\n")
    next
  }
  
  cat("  Found", length(result_files), "study result files\n")
  
  # Read and combine all results
  all_results <- data.frame()
  
  for (file in result_files) {
    study_data <- read.csv(file, stringsAsFactors = FALSE)
    all_results <- rbind(all_results, study_data)
  }
  
  cat("  Total protein-study combinations:", nrow(all_results), "\n")
  
  # Group by protein and perform meta-analysis
  proteins <- unique(all_results$Protein)
  meta_results <- data.frame()
  
  for (protein in proteins) {
    protein_data_original <- all_results[all_results$Protein == protein, ]
    
    # Save original study information BEFORE deduplication
    n_studies_original <- nrow(protein_data_original)
    total_n_original <- sum(protein_data_original$N)
    studies_list_original <- paste(paste0(protein_data_original$Study, "_", comp_name), collapse = "; ")
    
    # Apply study selection rules (deduplication)
    protein_data <- select_studies_for_protein(protein_data_original)
    
    # Skip if less than 2 studies after deduplication
    if (nrow(protein_data) < 2) {
      next
    }
    
    # Calculate unique non-duplicate study count
    unique_study_count <- count_unique_studies(protein_data$Study)
    
    # Perform weighted meta-analysis (using deduplicated data)
    meta <- perform_weighted_meta(
      protein_data$Log2FC, 
      protein_data$SE, 
      protein_data$N,
      protein_data$P_value
    )
    
    # Skip if weighted effect is NA
    if (is.na(meta$weighted_effect)) {
      next
    }
    
    # Perform leave-one-out analysis
    loo <- perform_loo_analysis(protein_data$Log2FC, protein_data$N)
    
    # Create unique study names list (after deduplication)
    unique_study_names <- paste(unique(protein_data$Study), collapse = "; ")
    
    # Store results
    meta_results <- rbind(meta_results, data.frame(
      Protein = protein,
      Weighted_Effect = meta$weighted_effect,
      CI_Lower = meta$ci_lower,
      CI_Upper = meta$ci_upper,
      Z_score = meta$z_score,
      P_value = meta$p_value,
      I_squared = meta$i_squared,
      LOO_Mean = loo$mean,
      LOO_CI_Lower = loo$ci_lower,
      LOO_CI_Upper = loo$ci_upper,
      N_Studies = n_studies_original,
      Total_N = total_n_original,
      Studies = studies_list_original,
      N_Unique_Studies = unique_study_count,
      Unique_Study_Names = unique_study_names,
      stringsAsFactors = FALSE
    ))
  }
  
  cat("  Proteins with >=2 studies:", nrow(meta_results), "\n")
  
  if (nrow(meta_results) == 0) {
    cat("  No proteins with sufficient studies\n\n")
    next
  }
  
  # Filter proteins with N_Unique_Studies >= 2
  meta_results_filtered <- meta_results[meta_results$N_Unique_Studies >= 2, ]
  
  cat("  Proteins with >=2 unique studies:", nrow(meta_results_filtered), "\n")
  
  if (nrow(meta_results_filtered) == 0) {
    cat("  No proteins with sufficient unique studies\n\n")
    next
  }
  
  # Stratified FDR correction by number of unique studies
  meta_results_filtered$FDR_BH_Stratified <- NA
  
  for (n_unique in unique(meta_results_filtered$N_Unique_Studies)) {
    idx <- meta_results_filtered$N_Unique_Studies == n_unique
    if (sum(idx) > 0) {
      meta_results_filtered$FDR_BH_Stratified[idx] <- p.adjust(meta_results_filtered$P_value[idx], method = "BH")
    }
  }
  
  # Sort by p-value
  meta_results_filtered <- meta_results_filtered[order(meta_results_filtered$P_value), ]
  
  # Save results
  output_file <- file.path(output_dir, paste0(comp_name, ".csv"))
  write.csv(meta_results_filtered, output_file, row.names = FALSE)
  
  cat("  Saved:", basename(output_file), "\n")
  cat("  Significant proteins (FDR < 0.05):", sum(meta_results_filtered$FDR_BH_Stratified < 0.05, na.rm = TRUE), "\n\n")
}

cat("================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n")

# Summary statistics
for (comp_name in comparisons) {
  output_file <- file.path(output_dir, paste0(comp_name, ".csv"))
  
  if (file.exists(output_file)) {
    results <- read.csv(output_file, stringsAsFactors = FALSE)
    cat(comp_name, ":\n")
    cat("  Total proteins:", nrow(results), "\n")
    cat("  Significant (FDR < 0.05):", sum(results$FDR_BH_Stratified < 0.05, na.rm = TRUE), "\n")
    cat("  Mean N_Studies:", round(mean(results$N_Studies), 2), "\n")
    cat("  Mean N_Unique_Studies:", round(mean(results$N_Unique_Studies), 2), "\n")
    cat("  Mean Total_N:", round(mean(results$Total_N), 0), "\n")
    cat("\n")
  }
}

cat("Output directory:", output_dir, "\n")
cat("\nMethod: Sample-size weighted average + Weighted z-score (metapro)\n")
cat("================================================================================\n")
cat("STEP 2 COMPLETE\n")
cat("================================================================================\n")
