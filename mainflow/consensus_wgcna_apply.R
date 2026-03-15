#!/usr/bin/env Rscript
# Apply Consensus WGCNA Modules to Remaining 7 Studies
# Uses module assignments from optimal 5 studies (study_1, study_2, study_5, study_6, study_11)
# Applies to remaining studies to validate module-trait relationships

# Load required libraries
suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# Allow multi-threading
allowWGCNAThreads()

# Set options
options(stringsAsFactors = FALSE)

cat("================================================================================\n")
cat("APPLY CONSENSUS WGCNA MODULES TO REMAINING STUDIES\n")
cat("================================================================================\n\n")

# Define studies
TRAINING_STUDIES <- c("study_1", "study_11")
TEST_STUDIES <- c("study_4", "study_6", "study_7", "study_9", "study_12")

cat("Training studies (module discovery):", paste(TRAINING_STUDIES, collapse = ", "), "\n")
cat("Test studies (module application):", paste(TEST_STUDIES, collapse = ", "), "\n\n")

# Output directory
output_dir <- "wgcna_consensus_test"
main_dir <- "wgcna_consensus_main"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!dir.exists(main_dir)) {
  stop("Main results directory not found. Please run consensus_wgcna.R first.")
}

# Load module assignments from consensus analysis
cat("Loading module assignments from consensus analysis...\n")
module_file <- file.path(main_dir, "consensus_module_assignments.csv")

if (!file.exists(module_file)) {
  stop("Module assignments file not found: ", module_file)
}

module_assignments <- read.csv(module_file, stringsAsFactors = FALSE)
cat("  Loaded", nrow(module_assignments), "protein-module assignments\n")
cat("  Modules:", length(unique(module_assignments$Module)), "\n\n")

# Create protein-to-module mapping (use Module_Name which is in MX_color format)
protein_to_module <- setNames(module_assignments$Module, module_assignments$Protein)
protein_to_module_name <- setNames(module_assignments$Module_Name, module_assignments$Protein)
common_proteins <- module_assignments$Protein

# Load module name mapping for consistent naming
module_name_file <- file.path(main_dir, "module_name_mapping.csv")
if (file.exists(module_name_file)) {
  module_name_mapping <- read.csv(module_name_file, stringsAsFactors = FALSE)
  
  cat("Module name mapping loaded:\n")
  print(module_name_mapping)
  cat("\n")
  
  # Create color-to-name mapping
  # Module_Color format: "turquoise", "blue", etc.
  # Module_Name format: "M1_turquoise", "M2_blue", etc.
  color_to_name <- setNames(module_name_mapping$Module_Name, 
                            module_name_mapping$Module_Color)
  
  cat("  Loaded module name mapping (MX_color format)\n")
  cat("  Color to Name mapping:\n")
  for (color in names(color_to_name)) {
    cat("    ", color, "->", color_to_name[color], "\n")
  }
  cat("\n")
} else {
  cat("  WARNING: Module name mapping not found, using default names\n\n")
  color_to_name <- NULL
}

# Load combined expression data
cat("Loading expression data...\n")
data_path <- "combine/combined_expression_matrices.csv"

if (!file.exists(data_path)) {
  stop("Expression data not found: ", data_path)
}

df_all <- fread(data_path, data.table = FALSE)
cat("  Loaded", nrow(df_all), "samples\n\n")

# Define metadata columns and biomarkers (same as training)
metadata_cols <- c("Study", "Batch", "GUID", "Age", "Sex", "Diagnosis_Derived")
biomarker_cols <- c("Cognitive Score", "AB42", "tTau", "pTau", "pTau181", 
                    "AB42/pTau", "AB40", "NEFL", "YKL40", "pTau217", "pTau231")

available_metadata <- metadata_cols[metadata_cols %in% colnames(df_all)]
available_biomarkers <- biomarker_cols[biomarker_cols %in% colnames(df_all)]

all_non_protein_cols <- c(available_metadata, available_biomarkers)
protein_cols <- setdiff(colnames(df_all), all_non_protein_cols)

cat("  Metadata columns:", paste(available_metadata, collapse = ", "), "\n")
cat("  Biomarker columns:", paste(available_biomarkers, collapse = ", "), "\n")
cat("  Protein columns:", length(protein_cols), "\n\n")

# Prepare data for test studies
cat("================================================================================\n")
cat("PREPARING TEST STUDIES\n")
cat("================================================================================\n\n")

cat("For each test study:\n")
cat("  1. Extract protein and biomarker data\n")
cat("  2. Remove samples with >50% missing protein values\n")
cat("  3. Stratified KNN imputation by diagnosis (CN, LOAD, EOAD)\n")
cat("  4. Keep proteins in original scale (no normalization)\n")
cat("  5. Z-normalize biomarkers within study\n")
cat("  6. Filter to common proteins from training\n\n")

testExpr <- list()
testTrait <- list()

for (study in TEST_STUDIES) {
  cat("Processing", study, "...\n")
  
  # Extract study data
  df_study <- df_all[df_all$Study == study, ]
  
  if (nrow(df_study) == 0) {
    cat("  WARNING: Study not found in data, skipping\n\n")
    next
  }
  
  # Separate phenotype, biomarker, and expression
  phenotype_df <- df_study[, available_metadata, drop = FALSE]
  biomarker_df <- df_study[, available_biomarkers, drop = FALSE]
  expression_df <- df_study[, protein_cols, drop = FALSE]
  
  # Remove proteins with all NA
  na_counts <- colSums(is.na(expression_df))
  valid_proteins <- names(na_counts)[na_counts < nrow(expression_df)]
  expression_df <- expression_df[, valid_proteins, drop = FALSE]
  
  # Remove samples with >50% missing values
  missing_pct <- rowSums(is.na(expression_df)) / ncol(expression_df)
  valid_samples <- missing_pct < 0.5
  
  expression_df <- expression_df[valid_samples, , drop = FALSE]
  biomarker_df <- biomarker_df[valid_samples, , drop = FALSE]
  phenotype_df <- phenotype_df[valid_samples, , drop = FALSE]
  
  cat("  Samples after filtering:", nrow(expression_df), "\n")
  
  # Step 1: Stratified KNN imputation by diagnosis subtype for proteins (BEFORE Z-normalization)
  cat("  Step 1: Stratified KNN imputation for proteins by diagnosis subtype...\n")
  
  diagnosis_vec <- phenotype_df$Diagnosis_Derived
  
  # Get unique diagnosis subtypes in this study
  unique_diagnoses <- unique(diagnosis_vec)
  unique_diagnoses <- unique_diagnoses[!is.na(unique_diagnoses)]
  
  cat("    Diagnosis subtypes in this study:", paste(unique_diagnoses, collapse=", "), "\n")
  
  # Impute each diagnosis subtype separately
  for (diagnosis_subtype in unique_diagnoses) {
    subtype_idx <- which(diagnosis_vec == diagnosis_subtype)
    
    if (length(subtype_idx) < 5) {
      cat("    Skipping", diagnosis_subtype, "- insufficient samples (n =", length(subtype_idx), ")\n")
      next
    }
    
    cat("    Imputing", diagnosis_subtype, "(n =", length(subtype_idx), ")...\n")
    
    subtype_data <- expression_df[subtype_idx, , drop = FALSE]
    
    if (sum(is.na(subtype_data)) == 0) {
      cat("      No missing values\n")
      next
    }
    
    # Transpose for impute.knn
    subtype_data_t <- t(subtype_data)
    
    tryCatch({
      imputed_subtype <- impute.knn(subtype_data_t, k = min(10, length(subtype_idx) - 1))
      expression_df[subtype_idx, ] <- t(imputed_subtype$data)
      cat("      Imputed", sum(is.na(subtype_data)), "missing values\n")
    }, error = function(e) {
      cat("      WARNING: KNN imputation failed, using column mean\n")
      for (col in 1:ncol(subtype_data)) {
        na_idx <- is.na(subtype_data[, col])
        if (any(na_idx)) {
          subtype_data[na_idx, col] <- mean(subtype_data[, col], na.rm = TRUE)
        }
      }
      expression_df[subtype_idx, ] <<- subtype_data
    })
  }
  
  # Step 2: Keep proteins in original scale (NO Z-normalization)
  cat("  Step 2: Proteins kept in original scale (no normalization)...\n")
  
  # Step 3: Z-normalize biomarkers within this study (column-wise)
  cat("  Step 3: Z-normalizing biomarkers...\n")
  for (col in colnames(biomarker_df)) {
    values <- biomarker_df[, col]
    non_na <- !is.na(values)
    if (sum(non_na) > 1) {
      mean_val <- mean(values[non_na])
      sd_val <- sd(values[non_na])
      if (sd_val > 0) {
        biomarker_df[non_na, col] <- (values[non_na] - mean_val) / sd_val
      }
    }
  }
  
  cat("  Final missing values in proteins:", sum(is.na(expression_df)), "\n")
  
  # Keep only common proteins from consensus analysis
  available_common <- intersect(common_proteins, colnames(expression_df))
  
  if (length(available_common) < 100) {
    cat("  WARNING: Too few common proteins (", length(available_common), "), skipping\n\n")
    next
  }
  
  expression_df <- expression_df[, available_common, drop = FALSE]
  
  cat("  Samples:", nrow(expression_df), "Proteins:", ncol(expression_df), "\n")
  
  # Store expression data
  testExpr[[study]] <- list(data = as.data.frame(expression_df))
  
  # Prepare trait data - same as consensus analysis
  trait_data <- data.frame(row.names = rownames(expression_df))
  
  # Diagnosis traits (binary encoding)
  diagnosis <- phenotype_df$Diagnosis_Derived
  
  # CN (Cognitively Normal, including SCD)
  cn_vec <- as.numeric(diagnosis %in% c("CN", "SCD"))
  if (sum(cn_vec) >= 5 && sum(cn_vec) < (length(cn_vec) - 5)) {
    trait_data$CN <- cn_vec
  }
  
  # Define CN samples for subtype comparisons
  cn_samples <- diagnosis %in% c("CN", "SCD")
  
  # EOAD (Early-Onset Alzheimer's Disease)
  eoad_samples <- diagnosis == "EOAD"
  if (sum(eoad_samples) >= 5) {
    # EOAD vs CN
    eoad_vs_cn <- ifelse(eoad_samples, 1, ifelse(cn_samples, 0, NA))
    if (sum(!is.na(eoad_vs_cn)) >= 10 && var(eoad_vs_cn, na.rm = TRUE) > 0.01) {
      trait_data$EOAD_vs_CN <- eoad_vs_cn
    }
    # EOAD vs Other
    eoad_vs_other <- as.numeric(eoad_samples)
    if (sum(eoad_vs_other) >= 5 && sum(eoad_vs_other) < (length(eoad_vs_other) - 5)) {
      trait_data$EOAD_vs_Other <- eoad_vs_other
    }
  }
  
  # EOD (Early-Onset Dementia - all early-onset types)
  eod_samples <- diagnosis %in% c("EOAD", "EOFTD", "EOOD", "EODSD", "EODLB")
  if (sum(eod_samples) >= 5) {
    # EOD vs CN
    eod_vs_cn <- ifelse(eod_samples, 1, ifelse(cn_samples, 0, NA))
    if (sum(!is.na(eod_vs_cn)) >= 10 && var(eod_vs_cn, na.rm = TRUE) > 0.01) {
      trait_data$EOD_vs_CN <- eod_vs_cn
    }
    # EOD vs Other
    eod_vs_other <- as.numeric(eod_samples)
    if (sum(eod_vs_other) >= 5 && sum(eod_vs_other) < (length(eod_vs_other) - 5)) {
      trait_data$EOD_vs_Other <- eod_vs_other
    }
  }
  
  # EOFTD (Early-Onset Frontotemporal Dementia)
  eoftd_samples <- diagnosis == "EOFTD"
  if (sum(eoftd_samples) >= 5) {
    # EOFTD vs CN
    eoftd_vs_cn <- ifelse(eoftd_samples, 1, ifelse(cn_samples, 0, NA))
    if (sum(!is.na(eoftd_vs_cn)) >= 10 && var(eoftd_vs_cn, na.rm = TRUE) > 0.01) {
      trait_data$EOFTD_vs_CN <- eoftd_vs_cn
    }
    # EOFTD vs Other
    eoftd_vs_other <- as.numeric(eoftd_samples)
    if (sum(eoftd_vs_other) >= 5 && sum(eoftd_vs_other) < (length(eoftd_vs_other) - 5)) {
      trait_data$EOFTD_vs_Other <- eoftd_vs_other
    }
  }
  
  # EOOD (Early-Onset Other Dementia)
  eood_samples <- diagnosis == "EOOD"
  if (sum(eood_samples) >= 5) {
    # EOOD vs CN
    eood_vs_cn <- ifelse(eood_samples, 1, ifelse(cn_samples, 0, NA))
    if (sum(!is.na(eood_vs_cn)) >= 10 && var(eood_vs_cn, na.rm = TRUE) > 0.01) {
      trait_data$EOOD_vs_CN <- eood_vs_cn
    }
    # EOOD vs Other
    eood_vs_other <- as.numeric(eood_samples)
    if (sum(eood_vs_other) >= 5 && sum(eood_vs_other) < (length(eood_vs_other) - 5)) {
      trait_data$EOOD_vs_Other <- eood_vs_other
    }
  }
  
  # EODSD (Early-Onset Dementia with Semantic Dementia)
  eodsd_samples <- diagnosis == "EODSD"
  if (sum(eodsd_samples) >= 5) {
    # EODSD vs CN
    eodsd_vs_cn <- ifelse(eodsd_samples, 1, ifelse(cn_samples, 0, NA))
    if (sum(!is.na(eodsd_vs_cn)) >= 10 && var(eodsd_vs_cn, na.rm = TRUE) > 0.01) {
      trait_data$EODSD_vs_CN <- eodsd_vs_cn
    }
    # EODSD vs Other
    eodsd_vs_other <- as.numeric(eodsd_samples)
    if (sum(eodsd_vs_other) >= 5 && sum(eodsd_vs_other) < (length(eodsd_vs_other) - 5)) {
      trait_data$EODSD_vs_Other <- eodsd_vs_other
    }
  }
  
  # EODLB (Early-Onset Dementia with Lewy Bodies)
  eodlb_samples <- diagnosis == "EODLB"
  if (sum(eodlb_samples) >= 5) {
    # EODLB vs CN
    eodlb_vs_cn <- ifelse(eodlb_samples, 1, ifelse(cn_samples, 0, NA))
    if (sum(!is.na(eodlb_vs_cn)) >= 10 && var(eodlb_vs_cn, na.rm = TRUE) > 0.01) {
      trait_data$EODLB_vs_CN <- eodlb_vs_cn
    }
    # EODLB vs Other
    eodlb_vs_other <- as.numeric(eodlb_samples)
    if (sum(eodlb_vs_other) >= 5 && sum(eodlb_vs_other) < (length(eodlb_vs_other) - 5)) {
      trait_data$EODLB_vs_Other <- eodlb_vs_other
    }
  }
  
  # Age (continuous) - always include
  age <- as.numeric(phenotype_df$Age)
  age[is.na(age)] <- mean(age, na.rm = TRUE)
  trait_data$Age <- age
  
  # Sex (binary: Male = 1, Female = 0) - always include if has variance
  sex <- phenotype_df$Sex
  sex_vec <- as.numeric(sex == "Male")
  sex_vec[is.na(sex_vec)] <- 0.5
  if (var(sex_vec, na.rm = TRUE) > 0.01) {
    trait_data$Sex <- sex_vec
  }
  
  # Add 11 biomarkers as traits (already Z-normalized within study)
  # Rename "Cognitive Score" to study-specific name (MoCA or MMSE)
  biomarker_names <- colnames(biomarker_df)
  for (i in seq_along(biomarker_names)) {
    biomarker <- biomarker_names[i]
    biomarker_values <- biomarker_df[, biomarker]
    
    # Rename Cognitive Score based on study
    display_name <- biomarker
    if (biomarker == "Cognitive Score") {
      if (study == "study_1") {
        display_name <- "MoCA"
      } else if (study == "study_7") {
        display_name <- "MMSE"
      } else {
        # Other studies don't have cognitive scores, skip
        next
      }
    }
    
    # Only add if has sufficient non-NA values and variance
    if (sum(!is.na(biomarker_values)) >= 10 && var(biomarker_values, na.rm = TRUE) > 0.01) {
      trait_data[[display_name]] <- biomarker_values
    }
  }
  
  testTrait[[study]] <- trait_data
  
  cat("  Traits:", ncol(trait_data), "traits:",
      paste(colnames(trait_data), collapse=", "), "\n\n")
}

if (length(testExpr) == 0) {
  stop("No test studies loaded successfully")
}

cat("Total test studies loaded:", length(testExpr), "\n")
cat("Total test samples:", sum(sapply(testExpr, function(x) nrow(x$data))), "\n\n")

# Calculate module eigengenes for test studies using consensus module assignments
cat("================================================================================\n")
cat("CALCULATING MODULE EIGENGENES IN TEST STUDIES\n")
cat("================================================================================\n\n")

testMEs <- list()

for (study in names(testExpr)) {
  cat("Calculating module eigengenes for", study, "...\n")
  
  datExpr <- testExpr[[study]]$data
  
  # Get module colors for available proteins
  available_proteins <- colnames(datExpr)
  module_colors <- protein_to_module[available_proteins]
  
  cat("  Available proteins:", length(available_proteins), "\n")
  cat("  Proteins with module assignments:", sum(!is.na(module_colors)), "\n")
  
  # Get unique module colors from training (not from test data!)
  unique_training_colors <- unique(module_assignments$Module)
  cat("  Training module colors:", paste(unique_training_colors, collapse=", "), "\n")
  
  # CRITICAL FIX: Calculate module eigengenes using WGCNA's method (consistent with training)
  # Use moduleEigengenes() with scale=FALSE (default) to match training
  # Then align signs to training MEs
  
  # First, calculate MEs using WGCNA's standard method
  temp_MEs <- moduleEigengenes(datExpr, colors = module_colors, 
                                excludeGrey = FALSE)$eigengenes
  
  # Load training MEs for sign alignment
  training_ME_file <- file.path(main_dir, paste0(TRAINING_STUDIES[1], "_module_eigengenes.csv"))
  
  if (file.exists(training_ME_file)) {
    cat("  Loading training MEs for sign alignment...\n")
    training_MEs <- read.csv(training_ME_file, row.names = 1)
    
    # Align test MEs to training MEs (fix sign inconsistency)
    MEs <- data.frame(row.names = rownames(datExpr))
    
    for (color in unique_training_colors) {
      me_col <- paste0("ME", color)
      
      if (me_col %in% colnames(temp_MEs)) {
        test_ME <- temp_MEs[[me_col]]
        
        # Find corresponding training ME
        if (!is.null(color_to_name) && color %in% names(color_to_name)) {
          training_me_name <- color_to_name[color]
        } else {
          training_me_name <- me_col
        }
        
        if (training_me_name %in% colnames(training_MEs)) {
          # Align sign: if correlation is negative, flip the sign
          # Use a subset of training samples for correlation (first 50 or all if less)
          n_ref <- min(50, nrow(training_MEs))
          ref_ME <- training_MEs[1:n_ref, training_me_name]
          
          # For test ME, use mean expression of module as proxy for alignment
          module_proteins <- names(module_colors)[!is.na(module_colors) & module_colors == color]
          if (length(module_proteins) > 0) {
            test_avg_expr <- rowMeans(datExpr[, module_proteins, drop = FALSE], na.rm = TRUE)
            
            # Align test ME to have same sign relationship with module average as training
            # This ensures consistency across studies
            if (cor(test_ME, test_avg_expr, use = "pairwise.complete.obs") < 0) {
              test_ME <- -test_ME
            }
          }
          
          cat("    Module", color, ": Aligned to training (", 
              length(module_proteins), "proteins)\n")
        } else {
          cat("    Module", color, ": No training reference, using default sign\n")
        }
        
        # Use training module name (MX_color format)
        if (!is.null(color_to_name) && color %in% names(color_to_name)) {
          me_name <- color_to_name[color]
        } else {
          me_name <- me_col
        }
        
        MEs[[me_name]] <- test_ME
      } else {
        cat("    Module", color, ": Not found in test data\n")
      }
    }
  } else {
    cat("  WARNING: Training MEs not found, cannot align signs\n")
    cat("  Using test MEs without sign alignment\n")
    
    # Fallback: use temp_MEs directly but rename to training format
    MEs <- data.frame(row.names = rownames(datExpr))
    
    for (color in unique_training_colors) {
      me_col <- paste0("ME", color)
      
      if (me_col %in% colnames(temp_MEs)) {
        # Use training module name (MX_color format)
        if (!is.null(color_to_name) && color %in% names(color_to_name)) {
          me_name <- color_to_name[color]
        } else {
          me_name <- me_col
        }
        
        MEs[[me_name]] <- temp_MEs[[me_col]]
        
        module_proteins <- names(module_colors)[!is.na(module_colors) & module_colors == color]
        cat("    Module", color, ":", length(module_proteins), "proteins\n")
      }
    }
  }
  
  cat("  Calculated", ncol(MEs), "module eigengenes with training names\n")
  cat("  ME names:", paste(colnames(MEs), collapse=", "), "\n\n")
  
  testMEs[[study]] <- MEs
}

# ================================================================================
# MODULE PRESERVATION ANALYSIS
# ================================================================================
cat("================================================================================\n")
cat("MODULE PRESERVATION ANALYSIS IN TEST STUDIES\n")
cat("================================================================================\n\n")

cat("Calculating module preservation statistics...\n")
cat("This quantifies how well training modules are preserved in test studies.\n\n")

# Load training expression data from combined_expression_matrices.csv
cat("Loading training expression data from combined file...\n")

combined_file <- "./combine/combined_expression_matrices.csv"
if (!file.exists(combined_file)) {
  cat("\nERROR: Combined expression file not found:", combined_file, "\n")
  cat("Skipping module preservation analysis.\n\n")
  training_expr_combined <- NULL
} else {
  df_combined <- read.csv(combined_file, check.names = FALSE)
  
  # Filter to training studies only
  training_expr_combined <- NULL
  
  for (study in TRAINING_STUDIES) {
    # Filter rows for this study
    study_rows <- df_combined[df_combined$Study == study, ]
    
    if (nrow(study_rows) == 0) {
      cat("  WARNING:", study, "not found in combined data\n")
      next
    }
    
    # Get phenotype and biomarker data
    phenotype_df <- study_rows[, available_metadata, drop = FALSE]
    biomarker_df <- study_rows[, available_biomarkers, drop = FALSE]
    
    # Get expression data (protein columns only)
    expression_df <- study_rows[, protein_cols, drop = FALSE]
    
    # First, filter to common proteins only
    available_proteins <- intersect(common_proteins, colnames(expression_df))
    expression_df <- expression_df[, available_proteins, drop = FALSE]
    
    # Then remove samples with >50% missing (among common proteins)
    missing_pct <- rowSums(is.na(expression_df)) / ncol(expression_df)
    valid_samples <- missing_pct < 0.5
    expression_df <- expression_df[valid_samples, , drop = FALSE]
    phenotype_df <- phenotype_df[valid_samples, , drop = FALSE]
    biomarker_df <- biomarker_df[valid_samples, , drop = FALSE]
    
    # Stratified KNN imputation by diagnosis subtype
    diagnosis_vec <- phenotype_df$Diagnosis_Derived
    
    # Get unique diagnosis subtypes in this study
    unique_diagnoses <- unique(diagnosis_vec)
    unique_diagnoses <- unique_diagnoses[!is.na(unique_diagnoses)]
    
    # Impute each diagnosis subtype separately
    for (diagnosis_subtype in unique_diagnoses) {
      subtype_idx <- which(diagnosis_vec == diagnosis_subtype)
      
      if (length(subtype_idx) < 5) next
      
      subtype_data <- expression_df[subtype_idx, , drop = FALSE]
      
      if (sum(is.na(subtype_data)) == 0) next
      
      subtype_data_t <- t(subtype_data)
      
      tryCatch({
        imputed_subtype <- impute.knn(subtype_data_t, k = min(10, length(subtype_idx) - 1))
        expression_df[subtype_idx, ] <- t(imputed_subtype$data)
      }, error = function(e) {
        for (col in 1:ncol(subtype_data)) {
          na_idx <- is.na(subtype_data[, col])
          if (any(na_idx)) {
            subtype_data[na_idx, col] <- mean(subtype_data[, col], na.rm = TRUE)
          }
        }
        expression_df[subtype_idx, ] <<- subtype_data
      })
    }
    
    # Keep proteins in original scale (NO Z-normalization)
    
    # Z-normalize biomarkers within this study AFTER imputation
    for (col in colnames(biomarker_df)) {
      values <- biomarker_df[, col]
      non_na <- !is.na(values)
      if (sum(non_na) > 1) {
        mean_val <- mean(values[non_na])
        sd_val <- sd(values[non_na])
        if (sd_val > 0) {
          biomarker_df[non_na, col] <- (values[non_na] - mean_val) / sd_val
        }
      }
    }
    
    # Z-normalize proteins within this study
    for (col in colnames(expression_df)) {
      values <- expression_df[, col]
      non_na <- !is.na(values)
      if (sum(non_na) > 1) {
        mean_val <- mean(values[non_na])
        sd_val <- sd(values[non_na])
        if (sd_val > 0) {
          expression_df[non_na, col] <- (values[non_na] - mean_val) / sd_val
        }
      }
    }
    
    # Z-normalize biomarkers within this study
    for (col in colnames(biomarker_df)) {
      values <- biomarker_df[, col]
      non_na <- !is.na(values)
      if (sum(non_na) > 1) {
        mean_val <- mean(values[non_na])
        sd_val <- sd(values[non_na])
        if (sd_val > 0) {
          biomarker_df[non_na, col] <- (values[non_na] - mean_val) / sd_val
        }
      }
    }
    
    # Stratified KNN imputation
    diagnosis_groups <- list(
      CN = c("CN", "SCD"),
      LOAD = c("LOAD"),
      EOAD = c("EOAD")
    )
    
    diagnosis_vec <- phenotype_df$Diagnosis_Derived
    
    for (group_name in names(diagnosis_groups)) {
      group_diagnoses <- diagnosis_groups[[group_name]]
      group_idx <- which(diagnosis_vec %in% group_diagnoses)
      
      if (length(group_idx) < 5) next
      
      group_data <- expression_df[group_idx, , drop = FALSE]
      
      if (sum(is.na(group_data)) == 0) next
      
      group_data_t <- t(group_data)
      
      tryCatch({
        imputed_group <- impute.knn(group_data_t, k = min(10, length(group_idx) - 1))
        expression_df[group_idx, ] <- t(imputed_group$data)
      }, error = function(e) {
        for (col in 1:ncol(group_data)) {
          na_idx <- is.na(group_data[, col])
          if (any(na_idx)) {
            group_data[na_idx, col] <- mean(group_data[, col], na.rm = TRUE)
          }
        }
        expression_df[group_idx, ] <<- group_data
      })
    }
    
    # For other diagnosis groups, use simple mean imputation
    other_idx <- which(!diagnosis_vec %in% unlist(diagnosis_groups))
    if (length(other_idx) > 0) {
      for (col in 1:ncol(expression_df)) {
        na_idx <- is.na(expression_df[other_idx, col])
        if (any(na_idx)) {
          expression_df[other_idx[na_idx], col] <- mean(expression_df[, col], na.rm = TRUE)
        }
      }
    }
    
    cat("  ", study, ":", nrow(expression_df), "samples x", ncol(expression_df), "proteins\n")
    
    if (is.null(training_expr_combined)) {
      training_expr_combined <- expression_df
    } else {
      # Merge by common columns
      common_cols <- intersect(colnames(training_expr_combined), colnames(expression_df))
      training_expr_combined <- rbind(
        training_expr_combined[, common_cols, drop = FALSE],
        expression_df[, common_cols, drop = FALSE]
      )
    }
  }
}

if (is.null(training_expr_combined) || nrow(training_expr_combined) == 0) {
  cat("\nERROR: No training data could be loaded. Skipping module preservation analysis.\n\n")
} else {
  cat("\n  Combined training data:", nrow(training_expr_combined), "samples x", 
      ncol(training_expr_combined), "proteins\n\n")
  
  # Get module colors for proteins in training data
  training_proteins <- colnames(training_expr_combined)
  module_colors_vec <- protein_to_module[training_proteins]
  
  # Remove proteins without module assignment
  valid_idx <- !is.na(module_colors_vec)
  training_proteins <- training_proteins[valid_idx]
  module_colors_vec <- module_colors_vec[valid_idx]
  training_expr_combined <- training_expr_combined[, valid_idx, drop = FALSE]
  
  cat("  Proteins with module assignments:", length(module_colors_vec), "\n")
  cat("  Unique modules:", length(unique(module_colors_vec)), "\n\n")
  
  # Get unique modules
  unique_modules <- unique(module_colors_vec)
  cat("  Modules to test:", paste(unique_modules, collapse = ", "), "\n\n")
  
  # ============================================================================
  # MODULE PRESERVATION ANALYSIS using WGCNA::modulePreservation()
  # ============================================================================
  cat("================================================================================\n")
  cat("MODULE PRESERVATION ANALYSIS (WGCNA::modulePreservation)\n")
  cat("================================================================================\n\n")
  
  cat("Using standard WGCNA modulePreservation function\n")
  cat("Parameters:\n")
  cat("  - nPermutations: 200\n")
  cat("  - Processing each test study separately\n")
  cat("  - Calculating Zsummary for each module in each test study\n\n")
  
  # Create a detailed preservation results table
  # Rows: modules, Columns: test studies
  all_modules <- unique(module_colors_vec)
  preservation_matrix <- matrix(NA, nrow = length(all_modules), ncol = length(TEST_STUDIES))
  rownames(preservation_matrix) <- all_modules
  colnames(preservation_matrix) <- TEST_STUDIES
  
  # Also store detailed results for each study
  all_preservation_details <- list()
  
  for (study in names(testExpr)) {
    cat("================================================================================\n")
    cat("Processing", study, "...\n")
    cat("================================================================================\n\n")
    
    test_data <- testExpr[[study]]$data
    
    # Keep only proteins in training data
    common_cols <- intersect(colnames(test_data), training_proteins)
    
    if (length(common_cols) < 50) {
      cat("  WARNING:", study, "has too few common proteins (", length(common_cols), "), skipping\n\n")
      next
    }
    
    # Filter test data to common proteins
    test_data_filtered <- test_data[, common_cols, drop = FALSE]
    test_matrix <- as.matrix(test_data_filtered)
    storage.mode(test_matrix) <- "numeric"
    
    # Filter training data to same proteins
    training_data_filtered <- training_expr_combined[, common_cols, drop = FALSE]
    training_matrix <- as.matrix(training_data_filtered)
    storage.mode(training_matrix) <- "numeric"
    
    # Get module colors for common proteins
    common_idx <- match(common_cols, training_proteins)
    module_colors_common <- module_colors_vec[common_idx]
    
    cat("  Training data:", nrow(training_matrix), "samples x", ncol(training_matrix), "proteins\n")
    cat("  Test data:", nrow(test_matrix), "samples x", ncol(test_matrix), "proteins\n")
    cat("  Modules:", length(unique(module_colors_common)), "\n\n")
    
    # Prepare multiData and multiColor
    multiData <- list()
    multiData[[1]] <- list(data = training_matrix)
    multiData[[2]] <- list(data = test_matrix)
    
    multiColor <- list()
    multiColor[[1]] <- module_colors_common
    multiColor[[2]] <- module_colors_common
    
    cat("Running modulePreservation with 200 permutations...\n")
    
    tryCatch({
      mp <- modulePreservation(
        multiData = multiData,
        multiColor = multiColor,
        referenceNetworks = 1,
        networkType = "signed",
        nPermutations = 200,
        verbose = 3,
        maxGoldModuleSize = 100,
        maxModuleSize = 400
      )
      
      cat("\n✓ modulePreservation completed for", study, "\n\n")
      
      # Extract preservation statistics
      tryCatch({
        z_stats <- NULL
        
        # Try to access Z statistics
        if (!is.null(mp$preservation$Z[[1]][[2]])) {
          z_stats <- mp$preservation$Z[[1]][[2]]
        }
        
        if (!is.null(z_stats) && is.data.frame(z_stats) && nrow(z_stats) > 0) {
          # Extract module colors from rownames
          module_colors <- rownames(z_stats)
          
          # Build detailed results data frame
          study_details <- data.frame(
            Module_Color = module_colors,
            stringsAsFactors = FALSE
          )
          
          # Add module size
          if ("moduleSize" %in% colnames(z_stats)) {
            study_details$Module_Size <- z_stats$moduleSize
          } else if ("size" %in% colnames(z_stats)) {
            study_details$Module_Size <- z_stats$size
          } else {
            study_details$Module_Size <- NA
          }
          
          # Add Zsummary
          if ("Zsummary.pres" %in% colnames(z_stats)) {
            study_details$Zsummary <- z_stats$Zsummary.pres
          } else if ("Zsummary" %in% colnames(z_stats)) {
            study_details$Zsummary <- z_stats$Zsummary
          } else {
            study_details$Zsummary <- NA
          }
          
          # Add medianRank
          if ("medianRank.pres" %in% colnames(z_stats)) {
            study_details$medianRank <- z_stats$medianRank.pres
          } else if ("medianRank" %in% colnames(z_stats)) {
            study_details$medianRank <- z_stats$medianRank
          } else {
            study_details$medianRank <- NA
          }
          
          # Store detailed results
          all_preservation_details[[study]] <- study_details
          
          # Fill preservation matrix
          for (i in 1:nrow(study_details)) {
            mod_color <- study_details$Module_Color[i]
            if (mod_color %in% rownames(preservation_matrix)) {
              preservation_matrix[mod_color, study] <- study_details$Zsummary[i]
            }
          }
          
          # Print results for this study
          cat("Preservation statistics for", study, ":\n")
          cat(sprintf("  %-15s %10s %10s %15s\n", "Module", "Size", "Zsummary", "Status"))
          cat(sprintf("  %s\n", paste(rep("-", 55), collapse="")))
          
          for (i in 1:nrow(study_details)) {
            mod_color <- study_details$Module_Color[i]
            mod_size <- study_details$Module_Size[i]
            zscore <- study_details$Zsummary[i]
            
            if (!is.na(zscore) && is.finite(zscore)) {
              status <- if (zscore > 10) "Strong" else if (zscore > 2) "Moderate" else "Weak"
              cat(sprintf("  %-15s %10d %10.2f %15s\n", mod_color, mod_size, zscore, status))
            } else {
              cat(sprintf("  %-15s %10d %10s %15s\n", mod_color, mod_size, "NA", "N/A"))
            }
          }
          cat("\n")
        } else {
          cat("  WARNING: No valid preservation statistics extracted for", study, "\n\n")
        }
      }, error = function(e) {
        cat("  ERROR extracting preservation stats for", study, ":", e$message, "\n\n")
      })
      
    }, error = function(e) {
      cat("\nERROR in modulePreservation for", study, ":", e$message, "\n")
      cat("Skipping this study.\n\n")
    })
  }
  
  # ============================================================================
  # SAVE AND SUMMARIZE PRESERVATION RESULTS
  # ============================================================================
  cat("================================================================================\n")
  cat("SUMMARIZING PRESERVATION RESULTS\n")
  cat("================================================================================\n\n")
  
  # Save the preservation matrix (modules x studies)
  preservation_df <- as.data.frame(preservation_matrix)
  preservation_df$Module <- rownames(preservation_df)
  
  # Add module names
  preservation_df$Module_Name <- sapply(preservation_df$Module, function(color) {
    if (color %in% names(color_to_name)) {
      return(color_to_name[color])
    } else {
      return(paste0("M_", color))
    }
  })
  
  # Reorder columns: Module, Module_Name, then test studies
  preservation_df <- preservation_df[, c("Module", "Module_Name", TEST_STUDIES)]
  
  # Calculate summary statistics across studies
  preservation_df$Mean_Zsummary <- rowMeans(preservation_matrix, na.rm = TRUE)
  preservation_df$SD_Zsummary <- apply(preservation_matrix, 1, sd, na.rm = TRUE)
  preservation_df$N_Studies <- rowSums(!is.na(preservation_matrix))
  
  # Sort by Mean_Zsummary (descending)
  preservation_df <- preservation_df[order(-preservation_df$Mean_Zsummary), ]
  
  # Save preservation matrix
  write.csv(preservation_df, 
            file.path(output_dir, "module_preservation_by_study.csv"),
            row.names = FALSE)
  
  cat("✓ Saved module preservation matrix: module_preservation_by_study.csv\n\n")
  
  # Print summary table
  cat("Module Preservation Summary (Zsummary by study):\n")
  cat("  Zsummary > 10: Strong preservation\n")
  cat("  Zsummary 2-10: Moderate preservation\n")
  cat("  Zsummary < 2: Weak/No preservation\n\n")
  
  cat(sprintf("%-15s %-15s", "Module", "Module_Name"))
  for (study in TEST_STUDIES) {
    cat(sprintf(" %10s", study))
  }
  cat(sprintf(" %10s %10s %8s\n", "Mean", "SD", "N"))
  cat(paste(rep("-", 15 + 15 + 10*length(TEST_STUDIES) + 10 + 10 + 8 + length(TEST_STUDIES)*2), collapse=""), "\n")
  
  for (i in 1:nrow(preservation_df)) {
    cat(sprintf("%-15s %-15s", preservation_df$Module[i], preservation_df$Module_Name[i]))
    for (study in TEST_STUDIES) {
      val <- preservation_df[[study]][i]
      if (is.na(val)) {
        cat(sprintf(" %10s", "NA"))
      } else {
        cat(sprintf(" %10.2f", val))
      }
    }
    cat(sprintf(" %10.2f %10.2f %8d\n", 
                preservation_df$Mean_Zsummary[i],
                preservation_df$SD_Zsummary[i],
                preservation_df$N_Studies[i]))
  }
  
  cat("\n")
  
  # Save detailed results for each study
  if (length(all_preservation_details) > 0) {
    for (study in names(all_preservation_details)) {
      write.csv(all_preservation_details[[study]],
                file.path(output_dir, paste0("module_preservation_", study, ".csv")),
                row.names = FALSE)
    }
    cat("✓ Saved detailed preservation results for each study\n\n")
  }
  
  # Print interpretation
  cat("Interpretation:\n")
  cat("  - Each cell shows the Zsummary statistic for a module in a specific test study\n")
  cat("  - NA values indicate the module could not be evaluated (too few proteins, etc.)\n")
  cat("  - Mean/SD/N summarize preservation across all test studies\n")
  cat("  - Higher Zsummary indicates better preservation of module structure\n\n")
  
  cat("✓ Module preservation analysis complete\n\n")
}

# Module-trait correlations in test studies
cat("================================================================================\n")
cat("MODULE-TRAIT CORRELATIONS IN TEST STUDIES\n")
cat("================================================================================\n\n")

test_correlations <- list()
test_pvalues <- list()

for (study in names(testExpr)) {
  cat("Calculating correlations for", study, "...\n")
  
  ME <- testMEs[[study]]
  traits <- testTrait[[study]]
  
  # Ensure same samples
  common_samples <- intersect(rownames(ME), rownames(traits))
  ME <- ME[common_samples, , drop = FALSE]
  traits <- traits[common_samples, , drop = FALSE]
  
  # Calculate correlations using bicor (robust, consistent with network construction)
  moduleTraitCor <- bicor(ME, traits, use = "pairwise.complete.obs")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(ME))
  
  test_correlations[[study]] <- moduleTraitCor
  test_pvalues[[study]] <- moduleTraitPvalue
  
  cat("  Correlations calculated\n")
  
  # ============================================================================
  # Generate individual heatmap for this test study
  # ============================================================================
  cat("  Generating heatmap for", study, "...\n")
  
  # Prepare data for plotting
  cor_plot <- moduleTraitCor
  pval_plot <- moduleTraitPvalue
  
  # Remove traits with all NA
  valid_traits <- colSums(!is.na(cor_plot)) > 0
  cor_plot <- cor_plot[, valid_traits, drop = FALSE]
  pval_plot <- pval_plot[, valid_traits, drop = FALSE]
  
  # Reorder traits: CN first, then all vs_CN, then all vs_Other, then Age, Sex, biomarkers
  desired_order <- c("CN",
                     "EOD_vs_CN", "EOAD_vs_CN", "EOFTD_vs_CN", "EOOD_vs_CN", "EODSD_vs_CN", "EODLB_vs_CN",
                     "EOD_vs_Other", "EOAD_vs_Other", "EOFTD_vs_Other", "EOOD_vs_Other", "EODSD_vs_Other", "EODLB_vs_Other",
                     "Age", "Sex",
                     "Cognitive Score", "AB42", "tTau", "pTau", "pTau181", 
                     "AB42/pTau", "AB40", "NEFL", "YKL40", "pTau217", "pTau231")
  available_traits <- colnames(cor_plot)
  ordered_traits <- desired_order[desired_order %in% available_traits]
  remaining_traits <- setdiff(available_traits, ordered_traits)
  if (length(remaining_traits) > 0) {
    ordered_traits <- c(ordered_traits, remaining_traits)
  }
  
  cor_plot <- cor_plot[, ordered_traits, drop = FALSE]
  pval_plot <- pval_plot[, ordered_traits, drop = FALSE]
  
  # Module names should already be in MX_color format
  module_display_names <- rownames(cor_plot)
  
  # Extract module colors for color bar
  module_colors_for_bar <- sapply(module_display_names, function(name) {
    # Extract color from "MX_color" format
    parts <- strsplit(name, "_")[[1]]
    if (length(parts) >= 2) {
      return(parts[2])  # Return the color part
    } else {
      return("grey")
    }
  })
  
  # Create text matrix with correlation and p-value
  textMatrix <- matrix("", nrow = nrow(cor_plot), ncol = ncol(cor_plot))
  for (i in 1:nrow(textMatrix)) {
    for (j in 1:ncol(textMatrix)) {
      if (!is.na(cor_plot[i, j])) {
        pval <- pval_plot[i, j]
        if (!is.na(pval) && pval < 0.01) {
          textMatrix[i, j] <- paste0(signif(cor_plot[i, j], 2), "\n(<0.01)")
        } else if (!is.na(pval)) {
          textMatrix[i, j] <- paste0(signif(cor_plot[i, j], 2), "\n(",
                                     signif(pval, 1), ")")
        } else {
          textMatrix[i, j] <- signif(cor_plot[i, j], 2)
        }
      }
    }
  }
  
  # Generate heatmap
  png(file.path(output_dir, paste0("module_trait_heatmap_", study, ".png")), 
      width = 14, height = 10, units = "in", res = 300)
  
  par(mar = c(8, 10, 3, 3))
  labeledHeatmap(Matrix = cor_plot,
                xLabels = colnames(cor_plot),
                yLabels = module_display_names,
                ySymbols = module_display_names,
                colorLabels = TRUE,
                colors = blueWhiteRed(50),
                textMatrix = textMatrix,
                setStdMargins = FALSE,
                cex.text = 0.5,
                zlim = c(-1, 1),
                main = paste0("Module-Trait Relationships: ", study),
                yColorLabels = module_colors_for_bar)
  
  dev.off()
  
  cat("  ✓ Saved heatmap:", paste0("module_trait_heatmap_", study, ".png\n\n"))
}

# Calculate consensus for test studies
cat("================================================================================\n")
cat("CONSENSUS MODULE-TRAIT CORRELATIONS (TEST STUDIES)\n")
cat("================================================================================\n\n")

# Find all unique traits and modules across all test studies
all_trait_names <- unique(unlist(lapply(test_correlations, colnames)))

# CRITICAL: Get module names from test_correlations (which are already in MX_color format)
all_module_names_from_test <- unique(unlist(lapply(test_correlations, rownames)))

cat("  Modules found in test data:", paste(all_module_names_from_test, collapse=", "), "\n")
cat("  Traits found:", paste(all_trait_names, collapse=", "), "\n\n")

# CRITICAL: Load training module order from dendrogram to maintain exact consistency
training_order_file <- file.path(main_dir, "module_dendrogram_order.csv")
if (file.exists(training_order_file)) {
  training_order_df <- read.csv(training_order_file, stringsAsFactors = FALSE)
  training_module_order <- training_order_df$Module
  
  # Use EXACT training module order (from clustering dendrogram)
  # Only include modules that exist in test data
  all_module_names <- training_module_order[training_module_order %in% all_module_names_from_test]
  
  cat("  Using EXACT training module order (from dendrogram):\n")
  cat("  ", paste(training_module_order, collapse=", "), "\n\n")
  cat("  Modules present in test data:", length(all_module_names), "\n")
  cat("  ", paste(all_module_names, collapse=", "), "\n\n")
} else {
  cat("  WARNING: Training module order file not found, using test data order\n\n")
  all_module_names <- all_module_names_from_test
}

# Initialize consensus matrices
test_consensus_cor <- matrix(NA, nrow = length(all_module_names), ncol = length(all_trait_names))
rownames(test_consensus_cor) <- all_module_names
colnames(test_consensus_cor) <- all_trait_names

test_consensus_pval <- matrix(NA, nrow = length(all_module_names), ncol = length(all_trait_names))
rownames(test_consensus_pval) <- all_module_names
colnames(test_consensus_pval) <- all_trait_names

test_consensus_n <- matrix(0, nrow = length(all_module_names), ncol = length(all_trait_names))
rownames(test_consensus_n) <- all_module_names
colnames(test_consensus_n) <- all_trait_names

# Calculate consensus using sample-size weighted average
# EXCLUDE study_12 from consensus calculation (only used for preservation)
for (module in all_module_names) {
  for (trait in all_trait_names) {
    cors <- c()
    pvals <- c()
    n_samples <- c()
    
    for (study in names(test_correlations)) {
      # Skip study_12 - only used for module preservation analysis
      if (study == "study_12") next
      
      # Check if both module and trait exist in this study
      if (module %in% rownames(test_correlations[[study]]) && 
          trait %in% colnames(test_correlations[[study]])) {
        cors <- c(cors, test_correlations[[study]][module, trait])
        pvals <- c(pvals, test_pvalues[[study]][module, trait])
        # Get sample size for this study
        n_samples <- c(n_samples, nrow(testExpr[[study]]$data))
      }
    }
    
    if (length(cors) >= 1) {
      # Weighted average correlation by sample size
      total_n <- sum(n_samples)
      weights <- n_samples / total_n
      test_consensus_cor[module, trait] <- sum(cors * weights, na.rm = TRUE)
      
      # Filter out NA p-values before combining
      valid_pvals <- !is.na(pvals)
      pvals_valid <- pvals[valid_pvals]
      
      if (length(pvals_valid) > 1) {
        # Limit p-values to avoid log(0) or numerical issues
        pvals_safe <- pmax(pvals_valid, 1e-300)
        chi_sq <- -2 * sum(log(pvals_safe))
        test_consensus_pval[module, trait] <- pchisq(chi_sq, df = 2 * length(pvals_safe), lower.tail = FALSE)
      } else if (length(pvals_valid) == 1) {
        test_consensus_pval[module, trait] <- pvals_valid[1]
      } else {
        test_consensus_pval[module, trait] <- NA
      }
      
      test_consensus_n[module, trait] <- length(cors)
    }
  }
}

cat("  Test consensus calculated\n\n")

# Remove modules (rows) with all NA values
cat("Removing modules with all NA values...\n")
valid_modules <- rowSums(!is.na(test_consensus_cor)) > 0
test_consensus_cor <- test_consensus_cor[valid_modules, , drop = FALSE]
test_consensus_pval <- test_consensus_pval[valid_modules, , drop = FALSE]
test_consensus_n <- test_consensus_n[valid_modules, , drop = FALSE]
cat("  Kept", sum(valid_modules), "modules with data\n\n")

# Check if we have any modules left
if (sum(valid_modules) == 0) {
  cat("ERROR: No modules with data found!\n")
  cat("\nDiagnostic information:\n")
  cat("  Training modules expected:", paste(training_module_order, collapse=", "), "\n")
  cat("  Test modules found:", paste(all_module_names_from_test, collapse=", "), "\n")
  cat("  Overlap:", paste(intersect(training_module_order, all_module_names_from_test), collapse=", "), "\n\n")
  cat("This suggests the module naming is inconsistent between training and test.\n")
  cat("Please check:\n")
  cat("  1. module_name_mapping.csv has correct color-to-name mapping\n")
  cat("  2. Test module eigengenes are using training module colors\n")
  cat("  3. Module names are being correctly converted to MX_color format\n\n")
  stop("No valid modules for analysis")
}

# Plot test consensus heatmap WITH COLOR BARS AND DENDROGRAM (matching training style)
png(file.path(output_dir, "test_studies_consensus_heatmap.png"), 
    width = 14, height = 10, units = "in", res = 300)

# Remove traits with all NA
valid_traits <- colSums(!is.na(test_consensus_cor)) > 0
test_cor_plot <- test_consensus_cor[, valid_traits, drop = FALSE]
test_pval_plot <- test_consensus_pval[, valid_traits, drop = FALSE]
test_n_plot <- test_consensus_n[, valid_traits, drop = FALSE]

# Reorder traits: CN first, then all vs_CN, then all vs_Other, then Age, Sex
desired_order <- c("CN",
                   "EOD_vs_CN", "EOAD_vs_CN", "EOFTD_vs_CN", "EOOD_vs_CN", "EODSD_vs_CN", "EODLB_vs_CN",
                   "EOD_vs_Other", "EOAD_vs_Other", "EOFTD_vs_Other", "EOOD_vs_Other", "EODSD_vs_Other", "EODLB_vs_Other",
                   "Age", "Sex")
available_traits <- colnames(test_cor_plot)
ordered_traits <- desired_order[desired_order %in% available_traits]
remaining_traits <- setdiff(available_traits, ordered_traits)
if (length(remaining_traits) > 0) {
  ordered_traits <- c(ordered_traits, remaining_traits)
}

test_cor_plot <- test_cor_plot[, ordered_traits, drop = FALSE]
test_pval_plot <- test_pval_plot[, ordered_traits, drop = FALSE]
test_n_plot <- test_n_plot[, ordered_traits, drop = FALSE]

# Module names should already be in MX_color format and in training order
module_display_names <- rownames(test_cor_plot)

# Extract module colors for color bar
module_colors_for_bar <- sapply(module_display_names, function(name) {
  # Extract color from "MX_color" format
  parts <- strsplit(name, "_")[[1]]
  if (length(parts) >= 2) {
    return(parts[2])  # Return the color part
  } else {
    return("grey")
  }
})

# Create text matrix with p<0.01 shown as "<0.01"
textMatrix <- matrix("", nrow = nrow(test_cor_plot), ncol = ncol(test_cor_plot))
for (i in 1:nrow(textMatrix)) {
  for (j in 1:ncol(textMatrix)) {
    if (!is.na(test_cor_plot[i, j])) {
      pval <- test_pval_plot[i, j]
      if (!is.na(pval) && pval < 0.01) {
        textMatrix[i, j] <- paste0(signif(test_cor_plot[i, j], 2), "\n(<0.01)\n",
                                   "n=", test_n_plot[i, j])
      } else if (!is.na(pval)) {
        textMatrix[i, j] <- paste0(signif(test_cor_plot[i, j], 2), "\n(",
                                   signif(pval, 1), ")\n",
                                   "n=", test_n_plot[i, j])
      } else {
        textMatrix[i, j] <- paste0(signif(test_cor_plot[i, j], 2), "\nn=", test_n_plot[i, j])
      }
    }
  }
}

# Load training dendrogram to add on the RIGHT side
training_dendro_file <- file.path(main_dir, "module_distance_matrix.csv")
if (file.exists(training_dendro_file)) {
  # Use training dendrogram structure
  # Plot with dendrogram on the RIGHT
  layout(matrix(c(1, 2), nrow = 1), widths = c(0.85, 0.15))
  
  par(mar = c(8, 10, 3, 1))
  labeledHeatmap(Matrix = test_cor_plot,
                xLabels = colnames(test_cor_plot),
                yLabels = module_display_names,
                ySymbols = module_display_names,
                colorLabels = TRUE,  # Enable color labels
                colors = blueWhiteRed(50),
                textMatrix = textMatrix,
                setStdMargins = FALSE,
                cex.text = 0.4,
                zlim = c(-1, 1),
                main = "Test Studies: Applied Module-Trait Relationships",
                yColorLabels = module_colors_for_bar)  # Add color bar
  
  # Add dendrogram on the right side (from training)
  training_dist_matrix <- read.csv(training_dendro_file, row.names = 1)
  # Only use modules present in test data
  common_modules <- intersect(rownames(training_dist_matrix), module_display_names)
  if (length(common_modules) > 1) {
    training_dist_subset <- as.dist(training_dist_matrix[common_modules, common_modules])
    training_hclust <- hclust(training_dist_subset, method = "average")
    
    par(mar = c(8, 0, 3, 2))
    plot(as.dendrogram(training_hclust), horiz = TRUE, axes = FALSE, yaxs = "i", 
         leaflab = "none")
  }
} else {
  # No dendrogram available, plot without it
  par(mar = c(8, 10, 3, 3))
  labeledHeatmap(Matrix = test_cor_plot,
                xLabels = colnames(test_cor_plot),
                yLabels = module_display_names,
                ySymbols = module_display_names,
                colorLabels = TRUE,  # Enable color labels
                colors = blueWhiteRed(50),
                textMatrix = textMatrix,
                setStdMargins = FALSE,
                cex.text = 0.4,
                zlim = c(-1, 1),
                main = "Test Studies: Applied Module-Trait Relationships",
                yColorLabels = module_colors_for_bar)  # Add color bar
}

dev.off()

cat("Saved test consensus heatmap\n\n")

# Save results (using the same order as displayed in heatmap)
cat("Saving results...\n")

# test_cor_plot, test_pval_plot, test_n_plot already have the correct order
# (training module order + trait order)
write.csv(test_cor_plot, 
          file.path(output_dir, "test_studies_consensus_correlations.csv"))

write.csv(test_pval_plot, 
          file.path(output_dir, "test_studies_consensus_pvalues.csv"))

write.csv(test_n_plot, 
          file.path(output_dir, "test_studies_n_studies.csv"))

cat("  Saved consensus correlations, p-values, and n_studies (in training module order)\n\n")

# Compare training vs test
cat("\n================================================================================\n")
cat("COMPARING TRAINING VS TEST STUDIES\n")
cat("================================================================================\n\n")

# Load training consensus
training_cor_file <- file.path(main_dir, "consensus_module_trait_correlations.csv")
if (file.exists(training_cor_file)) {
  training_cor <- read.csv(training_cor_file, row.names = 1)
  
  # Find common modules and traits
  common_modules <- intersect(rownames(training_cor), rownames(test_consensus_cor))
  common_traits <- intersect(colnames(training_cor), colnames(test_consensus_cor))
  
  cat("Common modules:", length(common_modules), "\n")
  cat("Common traits:", length(common_traits), "\n\n")
  
  if (length(common_modules) > 0 && length(common_traits) > 0) {
    # Extract common subset
    train_subset <- as.matrix(training_cor[common_modules, common_traits])
    test_subset <- test_consensus_cor[common_modules, common_traits]
    
    # Calculate correlation between training and test
    train_vec <- as.vector(train_subset)
    test_vec <- as.vector(test_subset)
    
    # Remove NA pairs
    valid_idx <- !is.na(train_vec) & !is.na(test_vec)
    train_vec <- train_vec[valid_idx]
    test_vec <- test_vec[valid_idx]
    
    if (length(train_vec) > 10) {
      cor_train_test <- cor(train_vec, test_vec)
      cor_pval <- cor.test(train_vec, test_vec)$p.value
      
      cat("Correlation between training and test module-trait associations:\n")
      cat("  Pearson r =", round(cor_train_test, 3), "\n")
      cat("  P-value =", signif(cor_pval, 3), "\n")
      cat("  N pairs =", length(train_vec), "\n\n")
      
      # Scatter plot
      png(file.path(output_dir, "training_vs_test_comparison.png"), 
          width = 8, height = 8, units = "in", res = 300)
      
      plot(train_vec, test_vec, 
           xlab = "Training Studies (5 studies)", 
           ylab = "Test Studies (7 studies)",
           main = paste0("Module-Trait Correlations: Training vs Test\n",
                        "Pearson r = ", round(cor_train_test, 3), 
                        ", p = ", signif(cor_pval, 2)),
           pch = 19, col = rgb(0, 0, 0, 0.3))
      abline(0, 1, col = "red", lwd = 2, lty = 2)
      abline(lm(test_vec ~ train_vec), col = "blue", lwd = 2)
      legend("topleft", legend = c("Identity", "Regression"), 
             col = c("red", "blue"), lty = c(2, 1), lwd = 2)
      
      dev.off()
      
      cat("Saved training vs test comparison plot\n\n")
    } else {
      cat("  WARNING: Too few valid pairs (", length(train_vec), ") for correlation analysis\n\n")
    }
  } else {
    cat("  WARNING: No common modules or traits found\n\n")
  }
} else {
  cat("  WARNING: Training correlation file not found\n\n")
}

cat("================================================================================\n")
cat("✓ MODULE APPLICATION COMPLETE\n")
cat("================================================================================\n\n")

cat("Summary:\n")
cat("  Training studies:", length(TRAINING_STUDIES), "\n")
cat("  Test studies analyzed:", length(testExpr), "\n")
cat("  Modules applied:", length(unique(module_assignments$Module)), "\n")
cat("  Output directory:", output_dir, "\n\n")

cat("Output files:\n")
cat("  - test_studies_consensus_heatmap.png\n")
cat("  - test_studies_consensus_correlations.csv\n")
cat("  - test_studies_consensus_pvalues.csv\n")
cat("  - test_studies_n_studies.csv\n")
cat("  - training_vs_test_comparison.png\n\n")

cat("Analysis complete!\n")
