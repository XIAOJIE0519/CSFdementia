#!/usr/bin/env Rscript
# Consensus WGCNA Analysis - Optimal 4 Studies
# Studies: study_1, study_2, study_5, study_11
# Data integration: Simple merging (no batch correction)
# Imputation: KNN stratified by diagnosis (CN, LOAD, EOAD)
# Analyzes module-trait relationships with EOAD, LOAD subtypes

# Load required libraries
suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(impute)
  library(BRETIGEA)
})

# Allow multi-threading
allowWGCNAThreads()

# Set options
options(stringsAsFactors = FALSE)

cat("================================================================================\n")
cat("CONSENSUS WGCNA ANALYSIS - 4 STUDIES WITH SIMPLE DATA MERGING\n")
cat("================================================================================\n\n")

# Define optimal studies
OPTIMAL_STUDIES <- c("study_1", "study_11")

cat("Using optimal study combination:\n")
cat("  Studies:", paste(OPTIMAL_STUDIES, collapse = ", "), "\n")
cat("  Training: study_1, study_11\n")
cat("  Testing: study_4, study_6, study_7, study_9\n\n")

# Create output directory
output_dir <- "wgcna_consensus_main"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load combined expression data
cat("Loading expression data...\n")
data_path <- "combine/combined_expression_matrices.csv"

if (!file.exists(data_path)) {
  stop("Expression data not found: ", data_path)
}

df_all <- fread(data_path, data.table = FALSE)
cat("  Loaded", nrow(df_all), "samples\n\n")

# ================================================================================
# STEP 1: Identify metadata columns and biomarkers
# ================================================================================
cat("================================================================================\n")
cat("STEP 1: IDENTIFYING METADATA COLUMNS AND BIOMARKERS\n")
cat("================================================================================\n\n")

# Define metadata columns (non-protein columns)
metadata_cols <- c("Study", "Batch", "GUID", "Age", "Sex", "Diagnosis_Derived")

# Define 11 biomarkers
biomarker_cols <- c("Cognitive Score", "AB42", "tTau", "pTau", "pTau181", 
                    "AB42/pTau", "AB40", "NEFL", "YKL40", "pTau217", "pTau231")

# Check which columns exist in data
available_metadata <- metadata_cols[metadata_cols %in% colnames(df_all)]
available_biomarkers <- biomarker_cols[biomarker_cols %in% colnames(df_all)]

cat("  Metadata columns:", paste(available_metadata, collapse = ", "), "\n")
cat("  Biomarker columns:", paste(available_biomarkers, collapse = ", "), "\n\n")

# Get protein columns (everything except metadata and biomarkers)
all_non_protein_cols <- c(available_metadata, available_biomarkers)
protein_cols <- setdiff(colnames(df_all), all_non_protein_cols)

cat("  Total protein columns:", length(protein_cols), "\n\n")

# ================================================================================
# STEP 2: Process each study independently - Z-normalization and imputation
# ================================================================================
cat("================================================================================\n")
cat("STEP 2: PROCESSING EACH STUDY INDEPENDENTLY\n")
cat("================================================================================\n\n")

cat("For each study:\n")
cat("  1. Extract protein and biomarker data\n")
cat("  2. Remove samples with >50% missing protein values\n")
cat("  3. Stratified KNN imputation by diagnosis (CN, LOAD, EOAD)\n")
cat("  4. Keep proteins in original scale (no normalization)\n")
cat("  5. Z-normalize biomarkers within study\n\n")

# Store processed data for each study
study_data_list <- list()
study_biomarker_list <- list()
study_phenotype_list <- list()

for (study in OPTIMAL_STUDIES) {
  cat("Processing", study, "...\n")
  
  df_study <- df_all[df_all$Study == study, ]
  if (nrow(df_study) == 0) {
    cat("  WARNING: Study not found, skipping\n\n")
    next
  }
  
  # Extract phenotype data
  phenotype_df <- df_study[, available_metadata, drop = FALSE]
  
  # Extract protein data
  expression_df <- df_study[, protein_cols, drop = FALSE]
  
  # Extract biomarker data
  biomarker_df <- df_study[, available_biomarkers, drop = FALSE]
  
  # First, remove proteins with all NA in this study
  na_counts <- colSums(is.na(expression_df))
  valid_proteins <- names(na_counts)[na_counts < nrow(expression_df)]
  expression_df <- expression_df[, valid_proteins, drop = FALSE]
  
  cat("  Proteins with data in this study:", ncol(expression_df), "\n")
  
  # Then, remove samples with >50% missing protein values (among available proteins)
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
  
  # Create unique sample names using row indices to avoid duplicates
  sample_names <- paste(study, seq_len(nrow(phenotype_df)), sep = ".")
  rownames(expression_df) <- sample_names
  rownames(biomarker_df) <- sample_names
  rownames(phenotype_df) <- sample_names
  
  # Store processed data
  study_data_list[[study]] <- expression_df
  study_biomarker_list[[study]] <- biomarker_df
  study_phenotype_list[[study]] <- phenotype_df
  
  cat("  Final missing values in proteins:", sum(is.na(expression_df)), "\n")
  cat("  Final missing values in biomarkers:", sum(is.na(biomarker_df)), "\n\n")
}

# ================================================================================
# STEP 3: Identify common proteins across training studies
# ================================================================================
cat("================================================================================\n")
cat("STEP 3: IDENTIFYING COMMON PROTEINS ACROSS TRAINING STUDIES\n")
cat("================================================================================\n\n")

# Get proteins present in each study (after processing)
protein_presence <- list()
for (study in names(study_data_list)) {
  # Protein is "present" if it has <100% missing values
  present_proteins <- colnames(study_data_list[[study]])[colSums(!is.na(study_data_list[[study]])) > 0]
  protein_presence[[study]] <- present_proteins
  cat("  ", study, ":", length(present_proteins), "proteins\n")
}

# Get common proteins (intersection)
common_proteins <- Reduce(intersect, protein_presence)
cat("\n  Common proteins (intersection):", length(common_proteins), "\n\n")

# Filter each study to common proteins only
for (study in names(study_data_list)) {
  study_data_list[[study]] <- study_data_list[[study]][, common_proteins, drop = FALSE]
}

# ================================================================================
# STEP 4: Merge processed data from training studies
# ================================================================================
cat("================================================================================\n")
cat("STEP 4: MERGING PROCESSED DATA FROM TRAINING STUDIES\n")
cat("================================================================================\n\n")

# Combine all studies by simple row binding
normalized_data <- do.call(rbind, study_data_list)
normalized_biomarkers <- do.call(rbind, study_biomarker_list)

# Combine phenotype data
combined_traits_list <- list()
for (study in names(study_phenotype_list)) {
  phenotype_df <- study_phenotype_list[[study]]
  combined_traits_list[[study]] <- data.frame(
    Study = study,
    Diagnosis = phenotype_df$Diagnosis_Derived,
    Age = phenotype_df$Age,
    Sex = phenotype_df$Sex,
    row.names = rownames(phenotype_df)
  )
}
normalized_traits <- do.call(rbind, combined_traits_list)

cat("  Combined protein data:", nrow(normalized_data), "samples x", 
    ncol(normalized_data), "proteins\n")
cat("  Combined biomarker data:", nrow(normalized_biomarkers), "samples x",
    ncol(normalized_biomarkers), "biomarkers\n\n")

# Note: Imputation already done in STEP 2 for each study independently
imputed_data <- normalized_data

# ================================================================================
# STEP 5: Prepare data for WGCNA
# ================================================================================
cat("================================================================================\n")
cat("STEP 5: PREPARING DATA FOR WGCNA\n")
cat("================================================================================\n\n")

# Remove low-variance proteins
protein_vars <- apply(imputed_data, 2, var, na.rm = TRUE)
high_var_proteins <- names(protein_vars)[protein_vars > 0.01]
expression_data <- imputed_data[, high_var_proteins, drop = FALSE]

cat("  Proteins after variance filtering:", ncol(expression_data), "\n")

# Split by study for consensus analysis
multiExpr <- list()
multiTrait <- list()

for (study in OPTIMAL_STUDIES) {
  # Extract samples for this study based on Study prefix in rownames
  study_idx <- grep(paste0("^", study, "\\."), rownames(expression_data))
  
  if (length(study_idx) == 0) next
  
  study_expr <- expression_data[study_idx, , drop = FALSE]
  study_traits <- normalized_traits[study_idx, , drop = FALSE]
  study_biomarkers <- normalized_biomarkers[study_idx, , drop = FALSE]
  
  multiExpr[[study]] <- list(data = as.data.frame(study_expr))
  
  cat("  ", study, ":", nrow(study_expr), "samples x", ncol(study_expr), "proteins\n")
  
  # Prepare trait data
  trait_data <- data.frame(row.names = rownames(study_expr))
  
  diagnosis <- study_traits$Diagnosis
  
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
  age <- as.numeric(study_traits$Age)
  age[is.na(age)] <- mean(age, na.rm = TRUE)
  trait_data$Age <- age
  
  # Sex (binary: Male = 1, Female = 0) - always include if has variance
  sex <- study_traits$Sex
  sex_vec <- as.numeric(sex == "Male")
  sex_vec[is.na(sex_vec)] <- 0.5
  if (var(sex_vec, na.rm = TRUE) > 0.01) {
    trait_data$Sex <- sex_vec
  }
  
  # Add 11 biomarkers as traits (already Z-normalized within study)
  # Rename "Cognitive Score" to study-specific name (MoCA or MMSE)
  biomarker_names <- colnames(study_biomarkers)
  for (i in seq_along(biomarker_names)) {
    biomarker <- biomarker_names[i]
    biomarker_values <- study_biomarkers[, biomarker]
    
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
  
  multiTrait[[study]] <- trait_data
  
  cat("    Traits:", ncol(trait_data), "-", paste(colnames(trait_data), collapse=", "), "\n")
}

cat("\n")

# Check if we have enough studies
if (length(multiExpr) < 2) {
  stop("Need at least 2 studies for consensus analysis")
}

cat("Total studies loaded:", length(multiExpr), "\n")
cat("Total samples:", sum(sapply(multiExpr, function(x) nrow(x$data))), "\n\n")

# Check data structure
cat("Checking data structure...\n")
exprSize <- checkSets(multiExpr)
cat("  Expression data check passed\n")
cat("  All studies have", exprSize$nGenes, "genes\n\n")

# Choose soft-thresholding power
cat("================================================================================\n")
cat("CHOOSING SOFT-THRESHOLDING POWER\n")
cat("================================================================================\n\n")

powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))

# Calculate power for each study
powerTables <- list()

for (study in names(multiExpr)) {
  cat("Analyzing soft threshold for", study, "...\n")
  
  datExpr <- multiExpr[[study]]$data
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
  powerTables[[study]] <- sft$fitIndices
  
  # Find optimal power
  sft_df <- sft$fitIndices
  optimal_power <- sft_df$Power[which.max(sft_df$SFT.R.sq)]
  
  cat("  Optimal power:", optimal_power, 
      "(R^2 =", round(max(sft_df$SFT.R.sq), 3), ")\n\n")
}

# Use consensus power (median of optimal powers)
optimal_powers <- sapply(powerTables, function(x) x$Power[which.max(x$SFT.R.sq)])
consensus_power <- round(median(optimal_powers))

cat("Consensus soft-thresholding power:", consensus_power, "\n\n")

# Plot soft threshold analysis
png(file.path(output_dir, "soft_threshold_analysis.png"), 
    width = 12, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))

# Scale independence
plot(powerTables[[1]]$Power, powerTables[[1]]$SFT.R.sq, type = "n",
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (R^2)",
     main = "Scale Independence")
for (i in seq_along(powerTables)) {
  points(powerTables[[i]]$Power, powerTables[[i]]$SFT.R.sq, 
         col = i, pch = 19)
}
abline(h = 0.8, col = "red", lty = 2)
abline(v = consensus_power, col = "blue", lty = 2)
legend("bottomright", legend = names(powerTables), col = seq_along(powerTables), 
       pch = 19, cex = 0.8)

# Mean connectivity
plot(powerTables[[1]]$Power, powerTables[[1]]$mean.k., type = "n",
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean Connectivity")
for (i in seq_along(powerTables)) {
  points(powerTables[[i]]$Power, powerTables[[i]]$mean.k., 
         col = i, pch = 19)
}
abline(v = consensus_power, col = "blue", lty = 2)

dev.off()

cat("Saved soft threshold analysis plot\n\n")

# ================================================================================
# CONSENSUS NETWORK CONSTRUCTION AND MODULE DETECTION
# Using WGCNA official blockwiseConsensusModules() function
# ================================================================================
cat("================================================================================\n")
cat("CONSENSUS NETWORK CONSTRUCTION AND MODULE DETECTION\n")
cat("Using WGCNA official blockwiseConsensusModules() function\n")
cat("================================================================================\n\n")

cat("Running blockwiseConsensusModules...\n")
cat("  Network type: signed\n")
cat("  Correlation: bicor (robust)\n")
cat("  Power:", consensus_power, "\n")
cat("  Min module size: 10\n")
cat("  Merge cut height: 0.1 (correlation > 0.9)\n")
cat("  Deep split: 4\n")
cat("  minKMEtoStay: 0.30\n")
cat("  Reassignment P-value: 0.05\n")
cat("  pamRespectsDendro: TRUE\n\n")

# Run consensus module detection
consensusNet <- blockwiseConsensusModules(
  multiExpr,
  power = consensus_power,
  minModuleSize = 10,
  deepSplit = 4,
  pamRespectsDendro = TRUE,
  mergeCutHeight = 0.1,
  numericLabels = FALSE,
  minKMEtoStay = 0.30,
  reassignThreshold = 0.05,
  saveTOMs = FALSE,
  verbose = 3,
  maxBlockSize = 10000,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  consensusQuantile = 0.3
)

# Extract module colors
dynamicColors <- consensusNet$colors
cat("\n  Identified", length(unique(dynamicColors)), "modules (before kME-based merging)\n\n")

# Extract consensus dendrogram and TOM
conTree <- consensusNet$dendrograms[[1]]

# ================================================================================
# STEP: Module merging using kME (more robust than ME correlation alone)
# ================================================================================
cat("================================================================================\n")
cat("MODULE MERGING USING kME\n")
cat("================================================================================\n\n")

cat("Merging modules based on kME similarity (not just ME correlation)...\n")
cat("  This is more robust than merging based on ME similarity alone\n\n")

# Calculate module eigengenes for the first dataset
datExpr_for_merge <- multiExpr[[1]]$data
MEList_for_merge <- moduleEigengenes(datExpr_for_merge, colors = dynamicColors)
MEs_for_merge <- MEList_for_merge$eigengenes

# Calculate kME for all genes
cat("  Calculating kME for all genes...\n")
kME_all_genes <- bicor(datExpr_for_merge, MEs_for_merge, use = "pairwise.complete.obs")

# For each gene, check if it should be reassigned based on kME
cat("  Checking for genes that should be reassigned based on kME...\n")
for (i in 1:length(dynamicColors)) {
  current_color <- dynamicColors[i]
  if (current_color == "grey") next
  
  current_ME <- paste0("ME", current_color)
  if (current_ME %in% colnames(kME_all_genes)) {
    current_kME <- kME_all_genes[i, current_ME]
    
    # Find the module with highest kME
    max_kME <- max(kME_all_genes[i, ], na.rm = TRUE)
    best_ME <- colnames(kME_all_genes)[which.max(kME_all_genes[i, ])]
    best_color <- gsub("ME", "", best_ME)
    
    # If another module has much higher kME, consider reassignment
    if (!is.na(max_kME) && !is.na(current_kME) && 
        max_kME > current_kME + 0.1 && max_kME > 0.3) {
      dynamicColors[i] <- best_color
    }
  }
}

# Now merge similar modules based on ME correlation
cat("  Merging similar modules based on ME correlation...\n")
MEList_after_reassign <- moduleEigengenes(datExpr_for_merge, colors = dynamicColors)
MEs_after_reassign <- MEList_after_reassign$eigengenes
MEDiss <- 1 - bicor(MEs_after_reassign, use = "pairwise.complete.obs")
METree <- hclust(as.dist(MEDiss), method = "average")

# Merge modules with correlation > 0.9 (height < 0.1)
merge_kME <- mergeCloseModules(datExpr_for_merge, dynamicColors, 
                                cutHeight = 0.1, verbose = 3)

# Update colors after kME-based merging
dynamicColors <- merge_kME$colors
cat("\n  After kME-based merging:", length(unique(dynamicColors)), "modules\n\n")

# Plot dendrogram with module colors
png(file.path(output_dir, "consensus_dendrogram_modules.png"), 
    width = 12, height = 6, units = "in", res = 300)
plotDendroAndColors(conTree, dynamicColors, "Consensus Modules",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus Gene Dendrogram and Module Colors")
dev.off()

cat("Saved dendrogram plot\n\n")

# Calculate module eigengenes for each study
cat("Calculating module eigengenes for each study...\n")

MEs <- list()
for (study in names(multiExpr)) {
  datExpr <- multiExpr[[study]]$data
  MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
  MEs[[study]] <- MEList$eigengenes
}

cat("  Module eigengenes calculated for all studies\n\n")

# ================================================================================
# MODULE CLEANING: Remove low kME genes and reassign grey module genes
# ================================================================================
cat("================================================================================\n")
cat("MODULE CLEANING\n")
cat("================================================================================\n\n")

cat("Step 1: Remove genes with low intramodular connectivity (kME < 0.3)\n")
cat("Step 2: Reassign grey module genes with kME > 0.25 to their best module\n\n")

# Use the first dataset for cleaning
datExpr_clean <- multiExpr[[1]]$data

# Calculate kME for all genes
cat("  Calculating module membership (kME) for all genes...\n")
kME_matrix <- bicor(datExpr_clean, MEs[[1]], use = "pairwise.complete.obs")

# Initialize cleaned colors
cleaned_colors <- dynamicColors
n_removed <- 0
n_reassigned <- 0

# Step 1: Remove genes with low intramodular kME
cat("  Step 1: Removing genes with kME < 0.3 in their assigned module...\n")
for (i in 1:length(cleaned_colors)) {
  if (cleaned_colors[i] == "grey") next  # Skip grey module for now
  
  module_color <- cleaned_colors[i]
  ME_col <- paste0("ME", module_color)
  
  if (ME_col %in% colnames(kME_matrix)) {
    gene_kME <- kME_matrix[i, ME_col]
    
    if (!is.na(gene_kME) && gene_kME < 0.3) {
      cleaned_colors[i] <- "grey"
      n_removed <- n_removed + 1
    }
  }
}
cat("    Removed", n_removed, "genes with low kME (< 0.3)\n\n")

# Step 2: Reassign grey module genes with high kME
cat("  Step 2: Reassigning grey module genes with kME > 0.25...\n")
grey_indices <- which(cleaned_colors == "grey")

for (i in grey_indices) {
  # Find the module with highest kME for this gene
  gene_kMEs <- kME_matrix[i, ]
  
  # Remove grey module from consideration
  gene_kMEs <- gene_kMEs[!grepl("MEgrey", names(gene_kMEs))]
  
  if (length(gene_kMEs) > 0) {
    max_kME <- max(gene_kMEs, na.rm = TRUE)
    
    if (!is.na(max_kME) && max_kME > 0.25) {
      best_ME <- names(gene_kMEs)[which.max(gene_kMEs)]
      best_module <- gsub("ME", "", best_ME)
      cleaned_colors[i] <- best_module
      n_reassigned <- n_reassigned + 1
    }
  }
}
cat("    Reassigned", n_reassigned, "grey genes to modules (kME > 0.25)\n\n")

# Update module colors
dynamicColors <- cleaned_colors

cat("  Final module counts after cleaning:\n")
module_table <- table(dynamicColors)
for (mod in names(module_table)) {
  cat("    ", mod, ":", module_table[mod], "genes\n")
}
cat("\n")

# Rename modules by size (M1 = largest, M2 = second largest, etc.)
cat("Renaming modules by size...\n")
module_sizes <- table(dynamicColors)
module_sizes <- sort(module_sizes[names(module_sizes) != "grey"], decreasing = TRUE)

# Create mapping: color -> M number
color_to_M <- setNames(paste0("M", seq_along(module_sizes), "_", names(module_sizes)), names(module_sizes))
if ("grey" %in% dynamicColors) {
  grey_size <- sum(dynamicColors == "grey")
  color_to_M["grey"] <- paste0("M", length(module_sizes) + 1, "_grey")
}

# Store original colors for later use
original_colors <- dynamicColors

cat("  Module size ranking:\n")
for (i in seq_along(module_sizes)) {
  cat("    ", color_to_M[names(module_sizes)[i]], ":", names(module_sizes)[i], "(", module_sizes[i], "proteins)\n")
}
if ("grey" %in% dynamicColors) {
  cat("    ", color_to_M["grey"], ": grey (", sum(dynamicColors == "grey"), "proteins)\n")
}
cat("\n")

# Recalculate module eigengenes after cleaning
cat("  Recalculating module eigengenes after cleaning...\n")
MEs <- list()
for (study in names(multiExpr)) {
  datExpr <- multiExpr[[study]]$data
  MEList <- moduleEigengenes(datExpr, colors = original_colors)
  MEs_temp <- MEList$eigengenes
  
  # Rename ME columns to M format (without ME prefix)
  new_colnames <- colnames(MEs_temp)
  for (color in names(color_to_M)) {
    old_name <- paste0("ME", color)
    if (old_name %in% new_colnames) {
      new_colnames[new_colnames == old_name] <- color_to_M[color]
    }
  }
  colnames(MEs_temp) <- new_colnames
  
  MEs[[study]] <- MEs_temp
}
cat("  Module eigengenes recalculated\n\n")

# Correlate modules with traits
cat("================================================================================\n")
cat("MODULE-TRAIT CORRELATIONS\n")
cat("================================================================================\n\n")

# For each study, calculate module-trait correlations
all_correlations <- list()
all_pvalues <- list()

for (study in names(multiExpr)) {
  cat("Calculating correlations for", study, "...\n")
  
  ME <- MEs[[study]]
  traits <- multiTrait[[study]]
  
  # Ensure same samples
  common_samples <- intersect(rownames(ME), rownames(traits))
  ME <- ME[common_samples, , drop = FALSE]
  traits <- traits[common_samples, , drop = FALSE]
  
  # Calculate correlations using bicor (robust, consistent with network construction)
  nMEs <- ncol(ME)
  nTraits <- ncol(traits)
  
  moduleTraitCor <- bicor(ME, traits, use = "pairwise.complete.obs")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(ME))
  
  all_correlations[[study]] <- moduleTraitCor
  all_pvalues[[study]] <- moduleTraitPvalue
  
  # Plot heatmap for this study
  png(file.path(output_dir, paste0(study, "_module_trait_heatmap.png")), 
      width = 12, height = 8, units = "in", res = 300)
  
  # Reorder traits: CN first, then all vs_CN, then all vs_Other, then Age, Sex
  desired_order <- c("CN",
                     "EOD_vs_CN", "EOAD_vs_CN", "EOFTD_vs_CN", "EOOD_vs_CN", "EODSD_vs_CN", "EODLB_vs_CN",
                     "EOD_vs_Other", "EOAD_vs_Other", "EOFTD_vs_Other", "EOOD_vs_Other", "EODSD_vs_Other", "EODLB_vs_Other",
                     "Age", "Sex")
  available_traits <- colnames(traits)
  ordered_traits <- desired_order[desired_order %in% available_traits]
  remaining_traits <- setdiff(available_traits, ordered_traits)
  if (length(remaining_traits) > 0) {
    ordered_traits <- c(ordered_traits, remaining_traits)
  }
  
  moduleTraitCor_ordered <- moduleTraitCor[, ordered_traits, drop = FALSE]
  moduleTraitPvalue_ordered <- moduleTraitPvalue[, ordered_traits, drop = FALSE]
  
  # Create text matrix with correlations and p-values (p<0.01 shown as "<0.01")
  textMatrix <- matrix("", nrow = nrow(moduleTraitCor_ordered), ncol = ncol(moduleTraitCor_ordered))
  for (i in 1:nrow(textMatrix)) {
    for (j in 1:ncol(textMatrix)) {
      pval <- moduleTraitPvalue_ordered[i, j]
      if (pval < 0.01) {
        textMatrix[i, j] <- paste0(signif(moduleTraitCor_ordered[i, j], 2), "\n(<0.01)")
      } else {
        textMatrix[i, j] <- paste0(signif(moduleTraitCor_ordered[i, j], 2), "\n(",
                                   signif(pval, 1), ")")
      }
    }
  }
  
  par(mar = c(8, 10, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor_ordered,
                xLabels = colnames(moduleTraitCor_ordered),
                yLabels = gsub("^ME", "", colnames(ME)),
                ySymbols = gsub("^ME", "", colnames(ME)),
                colorLabels = FALSE,
                colors = blueWhiteRed(50),
                textMatrix = textMatrix,
                setStdMargins = FALSE,
                cex.text = 0.5,
                zlim = c(-1, 1),
                main = paste("Module-Trait Relationships -", study))
  
  dev.off()
  
  cat("  Saved heatmap for", study, "\n\n")
}

# Calculate consensus module-trait correlations
cat("Calculating consensus module-trait correlations...\n")

# Find all unique traits across all studies
all_trait_names <- unique(unlist(lapply(all_correlations, colnames)))
all_module_names <- rownames(all_correlations[[1]])

cat("  All traits found:", paste(all_trait_names, collapse=", "), "\n")

# Initialize consensus matrices
consensus_cor <- matrix(NA, nrow = length(all_module_names), ncol = length(all_trait_names))
rownames(consensus_cor) <- all_module_names
colnames(consensus_cor) <- all_trait_names

consensus_pval <- matrix(NA, nrow = length(all_module_names), ncol = length(all_trait_names))
rownames(consensus_pval) <- all_module_names
colnames(consensus_pval) <- all_trait_names

consensus_n_studies <- matrix(0, nrow = length(all_module_names), ncol = length(all_trait_names))
rownames(consensus_n_studies) <- all_module_names
colnames(consensus_n_studies) <- all_trait_names

# For each module-trait pair, calculate consensus from studies that have that trait
# Use sample-size weighted average instead of simple mean
for (module in all_module_names) {
  for (trait in all_trait_names) {
    # Collect correlations, p-values, and sample sizes from studies that have this trait
    cors <- c()
    pvals <- c()
    n_samples <- c()
    
    for (study in names(all_correlations)) {
      if (trait %in% colnames(all_correlations[[study]])) {
        cors <- c(cors, all_correlations[[study]][module, trait])
        pvals <- c(pvals, all_pvalues[[study]][module, trait])
        # Get sample size for this study
        n_samples <- c(n_samples, nrow(multiExpr[[study]]$data))
      }
    }
    
    # Calculate consensus if at least 1 study has this trait
    if (length(cors) >= 1) {
      # Weighted average correlation by sample size
      total_n <- sum(n_samples)
      weights <- n_samples / total_n
      consensus_cor[module, trait] <- sum(cors * weights, na.rm = TRUE)
      
      # Fisher's method for combining p-values (only if multiple studies)
      if (length(pvals) > 1) {
        chi_sq <- -2 * sum(log(pvals))
        consensus_pval[module, trait] <- pchisq(chi_sq, df = 2 * length(pvals), lower.tail = FALSE)
      } else {
        # Single study: use its p-value directly
        consensus_pval[module, trait] <- pvals[1]
      }
      
      # Record number of studies
      consensus_n_studies[module, trait] <- length(cors)
    }
  }
}

cat("  Consensus calculated\n")
cat("  Traits in consensus:", sum(colSums(!is.na(consensus_cor)) > 0), "\n\n")

# Plot consensus heatmap (only traits with data)
png(file.path(output_dir, "consensus_module_trait_heatmap.png"), 
    width = 14, height = 10, units = "in", res = 300)

# Remove traits with all NA
valid_traits <- colSums(!is.na(consensus_cor)) > 0
consensus_cor_plot <- consensus_cor[, valid_traits, drop = FALSE]
consensus_pval_plot <- consensus_pval[, valid_traits, drop = FALSE]
consensus_n_plot <- consensus_n_studies[, valid_traits, drop = FALSE]

# Reorder traits: CN first, then all vs_CN, then all vs_Other, then Age, Sex
desired_order <- c("CN",
                   "EOD_vs_CN", "EOAD_vs_CN", "EOFTD_vs_CN", "EOOD_vs_CN", "EODSD_vs_CN", "EODLB_vs_CN",
                   "EOD_vs_Other", "EOAD_vs_Other", "EOFTD_vs_Other", "EOOD_vs_Other", "EODSD_vs_Other", "EODLB_vs_Other",
                   "Age", "Sex")
available_traits <- colnames(consensus_cor_plot)
ordered_traits <- desired_order[desired_order %in% available_traits]

# Add any remaining traits not in desired order
remaining_traits <- setdiff(available_traits, ordered_traits)
if (length(remaining_traits) > 0) {
  ordered_traits <- c(ordered_traits, remaining_traits)
}

consensus_cor_plot <- consensus_cor_plot[, ordered_traits, drop = FALSE]
consensus_pval_plot <- consensus_pval_plot[, ordered_traits, drop = FALSE]
consensus_n_plot <- consensus_n_plot[, ordered_traits, drop = FALSE]

# Perform hierarchical clustering on modules
module_dist <- dist(consensus_cor_plot)
module_hclust <- hclust(module_dist, method = "average")

# Reorder modules by clustering
module_order <- module_hclust$order
consensus_cor_plot <- consensus_cor_plot[module_order, , drop = FALSE]
consensus_pval_plot <- consensus_pval_plot[module_order, , drop = FALSE]
consensus_n_plot <- consensus_n_plot[module_order, , drop = FALSE]

# Create text matrix with correlations, p-values (p<0.01 shown as "<0.01"), and n_studies
textMatrix <- matrix("", nrow = nrow(consensus_cor_plot), ncol = ncol(consensus_cor_plot))
for (i in 1:nrow(textMatrix)) {
  for (j in 1:ncol(textMatrix)) {
    if (!is.na(consensus_cor_plot[i, j])) {
      pval <- consensus_pval_plot[i, j]
      if (pval < 0.01) {
        textMatrix[i, j] <- paste0(signif(consensus_cor_plot[i, j], 2), "\n(<0.01)\n",
                                   "n=", consensus_n_plot[i, j])
      } else {
        textMatrix[i, j] <- paste0(signif(consensus_cor_plot[i, j], 2), "\n(",
                                   signif(pval, 1), ")\n",
                                   "n=", consensus_n_plot[i, j])
      }
    }
  }
}

# Plot with dendrogram on the RIGHT
layout(matrix(c(1, 2), nrow = 1), widths = c(0.85, 0.15))

par(mar = c(8, 10, 3, 1))
labeledHeatmap(Matrix = consensus_cor_plot,
              xLabels = colnames(consensus_cor_plot),
              yLabels = rownames(consensus_cor_plot),
              ySymbols = rownames(consensus_cor_plot),
              colorLabels = FALSE,
              colors = blueWhiteRed(50),
              textMatrix = textMatrix,
              setStdMargins = FALSE,
              cex.text = 0.4,
              zlim = c(-1, 1),
              main = "Consensus Module-Trait Relationships")

# Add dendrogram on the right side
par(mar = c(8, 0, 3, 2))
plot(as.dendrogram(module_hclust), horiz = TRUE, axes = FALSE, yaxs = "i", 
     leaflab = "none")


dev.off()

cat("Saved consensus heatmap\n\n")

# Save module clustering dendrogram information
cat("Saving module clustering dendrogram data...\n")

# Create a comprehensive dendrogram table
n_modules <- length(module_hclust$labels)
n_merges <- nrow(module_hclust$merge)

# Initialize detailed dendrogram table
dendrogram_table <- data.frame(
  Step = integer(),
  Cluster1_ID = integer(),
  Cluster1_Name = character(),
  Cluster1_Type = character(),
  Cluster2_ID = integer(),
  Cluster2_Name = character(),
  Cluster2_Type = character(),
  Merge_Height = numeric(),
  New_Cluster_ID = integer(),
  New_Cluster_Size = integer(),
  stringsAsFactors = FALSE
)

# Track cluster membership at each step
cluster_members <- list()
for (i in 1:n_modules) {
  # Use the module names from rownames (already in MX_color format from color_to_M)
  cluster_members[[i]] <- rownames(consensus_cor_plot)[i]
}

# Process each merge step
for (step in 1:n_merges) {
  cluster1_id <- module_hclust$merge[step, 1]
  cluster2_id <- module_hclust$merge[step, 2]
  merge_height <- module_hclust$height[step]
  new_cluster_id <- n_modules + step
  
  # Determine cluster 1 info
  if (cluster1_id < 0) {
    # Negative means it's an original module
    cluster1_type <- "Module"
    cluster1_name <- rownames(consensus_cor_plot)[-cluster1_id]
    cluster1_members <- cluster1_name
  } else {
    # Positive means it's a merged cluster
    cluster1_type <- "Merged_Cluster"
    cluster1_name <- paste0("Cluster_", cluster1_id)
    cluster1_members <- cluster_members[[cluster1_id]]
  }
  
  # Determine cluster 2 info
  if (cluster2_id < 0) {
    cluster2_type <- "Module"
    cluster2_name <- rownames(consensus_cor_plot)[-cluster2_id]
    cluster2_members <- cluster2_name
  } else {
    cluster2_type <- "Merged_Cluster"
    cluster2_name <- paste0("Cluster_", cluster2_id)
    cluster2_members <- cluster_members[[cluster2_id]]
  }
  
  # Store new cluster members
  cluster_members[[new_cluster_id]] <- c(cluster1_members, cluster2_members)
  new_cluster_size <- length(cluster_members[[new_cluster_id]])
  
  # Add to table
  dendrogram_table <- rbind(dendrogram_table, data.frame(
    Step = step,
    Cluster1_ID = cluster1_id,
    Cluster1_Name = cluster1_name,
    Cluster1_Type = cluster1_type,
    Cluster2_ID = cluster2_id,
    Cluster2_Name = cluster2_name,
    Cluster2_Type = cluster2_type,
    Merge_Height = merge_height,
    New_Cluster_ID = new_cluster_id,
    New_Cluster_Size = new_cluster_size,
    stringsAsFactors = FALSE
  ))
}

# Add cluster members as a separate column (comma-separated)
dendrogram_table$Cluster_Members <- sapply(dendrogram_table$New_Cluster_ID, function(id) {
  paste(cluster_members[[id]], collapse = "; ")
})

# Save comprehensive dendrogram table
write.csv(dendrogram_table, 
          file.path(output_dir, "module_dendrogram_detailed.csv"),
          row.names = FALSE)

# Also save the simple order information
module_dendrogram_order <- data.frame(
  Module = rownames(consensus_cor_plot),
  Cluster_Order = 1:length(rownames(consensus_cor_plot)),
  stringsAsFactors = FALSE
)

write.csv(module_dendrogram_order, 
          file.path(output_dir, "module_dendrogram_order.csv"),
          row.names = FALSE)

# Save distance matrix used for clustering
module_dist_matrix <- as.matrix(module_dist)
rownames(module_dist_matrix) <- rownames(consensus_cor_plot)
colnames(module_dist_matrix) <- rownames(consensus_cor_plot)
write.csv(module_dist_matrix,
          file.path(output_dir, "module_distance_matrix.csv"))

cat("  Saved detailed module dendrogram data:\n")
cat("    - module_dendrogram_detailed.csv: Complete merge history with cluster members\n")
cat("    - module_dendrogram_order.csv: Final module order in dendrogram\n")
cat("    - module_distance_matrix.csv: Distance matrix between modules\n\n")

# Save detailed results
cat("Saving detailed results...\n")

# Save CSV files with clustered order
write.csv(consensus_cor_plot, 
          file.path(output_dir, "consensus_module_trait_correlations.csv"))

write.csv(consensus_pval_plot, 
          file.path(output_dir, "consensus_module_trait_pvalues.csv"))

write.csv(consensus_n_plot, 
          file.path(output_dir, "consensus_module_trait_n_studies.csv"))

# Save module assignments with MX_color format (using color_to_M, NOT clustered order)
module_assignments <- data.frame(
  Protein = high_var_proteins,
  Module = original_colors,
  Module_Name = color_to_M[original_colors],
  stringsAsFactors = FALSE
)

write.csv(module_assignments, 
          file.path(output_dir, "consensus_module_assignments.csv"),
          row.names = FALSE)

# Save module color mapping (for reference only, not used in heatmap)
module_color_mapping <- data.frame(
  Module_Color = names(color_to_M),
  Module_Name = color_to_M,
  stringsAsFactors = FALSE
)

write.csv(module_color_mapping,
          file.path(output_dir, "module_name_mapping.csv"),
          row.names = FALSE)

# Save module eigengenes for each study
for (study in names(MEs)) {
  write.csv(MEs[[study]], 
            file.path(output_dir, paste0(study, "_module_eigengenes.csv")))
}

# Create summary report
cat("\n================================================================================\n")
cat("CONSENSUS WGCNA SUMMARY\n")
cat("================================================================================\n\n")

cat("Studies analyzed:", length(multiExpr), "\n")
cat("  Studies:", paste(names(multiExpr), collapse = ", "), "\n\n")

cat("Network parameters:\n")
cat("  Common proteins:", length(common_proteins), "\n")
cat("  Soft-thresholding power:", consensus_power, "\n")
cat("  Number of modules:", length(unique(dynamicColors)), "\n\n")

cat("Module sizes:\n")
module_sizes <- table(dynamicColors)
for (mod in names(module_sizes)) {
  cat("  ", mod, ":", module_sizes[mod], "proteins\n")
}

cat("\nSignificant module-trait associations (consensus p < 0.01):\n")
sig_associations <- which(consensus_pval < 0.01, arr.ind = TRUE)
if (nrow(sig_associations) > 0) {
  for (i in 1:nrow(sig_associations)) {
    row_idx <- sig_associations[i, 1]
    col_idx <- sig_associations[i, 2]
    module <- rownames(consensus_cor)[row_idx]
    trait <- colnames(consensus_cor)[col_idx]
    cor_val <- consensus_cor[row_idx, col_idx]
    pval <- consensus_pval[row_idx, col_idx]
    
    cat(sprintf("  %s - %s: r = %.3f, p = %.2e\n", 
                module, trait, cor_val, pval))
  }
} else {
  cat("  No significant associations at p < 0.01\n")
}

cat("\nOutput directory:", output_dir, "\n")

cat("\n================================================================================\n")
cat("â CONSENSUS WGCNA ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")

# Save session info
sink(file.path(output_dir, "session_info.txt"))
cat("Consensus WGCNA Analysis\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Session Info:\n")
print(sessionInfo())
sink()

cat("Analysis complete!\n")

# ================================================================================
# BRETIGEA CELL TYPE ENRICHMENT ANALYSIS (Fisher Exact Test)
# ================================================================================
cat("\n================================================================================\n")
cat("BRETIGEA CELL TYPE ENRICHMENT ANALYSIS\n")
cat("================================================================================\n\n")

cat("Using BRETIGEA package for Fisher exact test of module-cell type enrichment\n")
cat("Testing enrichment of 6 brain cell types (top 200 markers each)\n\n")

# Load BRETIGEA markers (top 200 genes for each cell type)
tryCatch({
  # Get cell type markers from BRETIGEA
  data("markers_df_human_brain", package = "BRETIGEA")
  
  # Extract top 200 markers for each cell type
  # BRETIGEA uses abbreviated cell type names
  cell_types_map <- c(
    "neurons" = "neu",
    "astrocytes" = "ast", 
    "oligodendrocytes" = "oli",
    "microglia" = "mic",
    "endothelial" = "end",
    "OPCs" = "opc"
  )
  
  bretigea_markers <- list()
  
  for (cell_type in names(cell_types_map)) {
    cell_abbr <- cell_types_map[cell_type]
    markers <- markers_df_human_brain[markers_df_human_brain$cell == cell_abbr, ]
    # Take top 200 markers
    bretigea_markers[[cell_type]] <- head(markers$markers, 200)
    cat("  ", cell_type, ": ", length(bretigea_markers[[cell_type]]), " markers\n", sep="")
  }
  cat("\n")
  
  # Extract all gene symbols from common proteins for background
  all_genes <- sapply(strsplit(common_proteins, "\\|"), function(x) x[1])
  
  # Perform Fisher exact test for each module
  fisher_results <- list()
  
  for (module_color in unique(original_colors)) {
    # Include grey module in cell type enrichment analysis
    
    module_proteins <- common_proteins[original_colors == module_color]
    
    # Extract gene symbols from protein IDs (format: SYMBOL|UNIPROT)
    module_genes <- sapply(strsplit(module_proteins, "\\|"), function(x) x[1])
    
    cat("Testing module", module_color, "(", length(module_genes), "proteins)...\n", sep="")
    
    for (cell_type in names(bretigea_markers)) {
      markers <- bretigea_markers[[cell_type]]
      
      # Create contingency table
      # In module & in markers
      in_both <- sum(module_genes %in% markers)
      # In module & not in markers
      in_module_only <- length(module_genes) - in_both
      # Not in module & in markers
      in_markers_only <- sum(markers %in% all_genes) - in_both
      # Not in module & not in markers
      background <- length(all_genes) - length(module_genes) - in_markers_only
      
      # Fisher exact test
      contingency <- matrix(c(in_both, in_module_only, in_markers_only, background), nrow = 2)
      fisher_test <- fisher.test(contingency, alternative = "greater")
      
      # Store results
      result_key <- paste(module_color, cell_type, sep = "_")
      fisher_results[[result_key]] <- data.frame(
        Module = module_color,
        Module_Name = color_to_M[module_color],
        Cell_Type = cell_type,
        Overlap = in_both,
        Module_Size = length(module_genes),
        Marker_Size = sum(markers %in% all_genes),
        OR = fisher_test$estimate,
        P_value = fisher_test$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine all results
  fisher_df <- do.call(rbind, fisher_results)
  rownames(fisher_df) <- NULL
  
  # Add FDR correction
  fisher_df$FDR <- p.adjust(fisher_df$P_value, method = "fdr")
  
  # Reorder columns
  fisher_df <- fisher_df[, c("Module", "Module_Name", "Cell_Type", "Overlap", "Module_Size", 
                              "Marker_Size", "OR", "P_value", "FDR")]
  
  # Sort by FDR
  fisher_df <- fisher_df[order(fisher_df$FDR), ]
  
  # Save results
  write.csv(fisher_df, 
            file.path(output_dir, "module_celltype_fisher_enrichment.csv"),
            row.names = FALSE)
  
  cat("\nâ?Fisher exact test completed\n")
  cat("  Results saved to: module_celltype_fisher_enrichment.csv\n\n")
  
  # Print significant enrichments (FDR < 0.05)
  sig_enrichments <- fisher_df[fisher_df$FDR < 0.05, ]
  if (nrow(sig_enrichments) > 0) {
    cat("Significant cell type enrichments (FDR < 0.05):\n")
    for (i in 1:nrow(sig_enrichments)) {
      cat(sprintf("  %s - %s: OR=%.2f, overlap=%d/%d, FDR=%.2e\n",
                  sig_enrichments$Module_Name[i],
                  sig_enrichments$Cell_Type[i],
                  sig_enrichments$OR[i],
                  sig_enrichments$Overlap[i],
                  sig_enrichments$Module_Size[i],
                  sig_enrichments$FDR[i]))
    }
  } else {
    cat("No significant enrichments at FDR < 0.05\n")
  }
  cat("\n")
  
}, error = function(e) {
  cat("ERROR in BRETIGEA analysis:", e$message, "\n")
  cat("Skipping cell type enrichment analysis\n\n")
})

# ================================================================================
# HUB GENE IDENTIFICATION (kME and GS)
# ================================================================================
cat("================================================================================\n")
cat("HUB GENE IDENTIFICATION\n")
cat("================================================================================\n\n")

cat("Identifying hub genes based on:\n")
cat("  - kME (module membership): correlation with module eigengene\n")
cat("  - GS (gene significance): correlation with traits\n\n")

# Calculate kME for all genes in all modules
cat("Calculating module membership (kME)...\n")

# Use combined expression data
datExpr_combined <- multiExpr[[1]]$data
ME_combined <- MEs[[1]]

# Calculate kME (correlation between gene expression and ME)
kME_all <- cor(datExpr_combined, ME_combined, use = "pairwise.complete.obs")

# Calculate GS for key traits
cat("Calculating gene significance (GS) for key traits...\n")

# Select key traits for GS calculation
key_traits <- c("EOAD_vs_CN", "LOAD_vs_CN", "EOD_vs_CN", "LOD_vs_CN", "Age")
available_key_traits <- intersect(key_traits, colnames(multiTrait[[1]]))

if (length(available_key_traits) > 0) {
  # Calculate GS (correlation between gene expression and trait)
  trait_data_combined <- multiTrait[[1]][, available_key_traits, drop = FALSE]
  
  # Ensure same samples
  common_samples <- intersect(rownames(datExpr_combined), rownames(trait_data_combined))
  datExpr_for_GS <- datExpr_combined[common_samples, ]
  trait_for_GS <- trait_data_combined[common_samples, , drop = FALSE]
  
  GS_all <- cor(datExpr_for_GS, trait_for_GS, use = "pairwise.complete.obs")
  
  cat("  Calculated GS for", ncol(GS_all), "traits\n\n")
} else {
  cat("  WARNING: No key traits available for GS calculation\n\n")
  GS_all <- NULL
}

# Identify hub genes for each module
hub_genes_list <- list()

for (module_color in unique(original_colors)) {
  if (module_color == "grey") next  # Skip grey module
  
  # Use high_var_proteins instead of common_proteins
  module_proteins <- high_var_proteins[original_colors == module_color]
  module_name <- color_to_M[module_color]
  
  # Get kME for this module - use module_name directly (already in M1_turquoise format)
  if (module_name %in% colnames(kME_all)) {
    kME_module <- kME_all[module_proteins, module_name]
    
    # Create hub gene data frame
    hub_df <- data.frame(
      Protein = module_proteins,
      Module = module_color,
      Module_Name = module_name,
      kME = kME_module,
      stringsAsFactors = FALSE
    )
    
    # Add GS if available
    if (!is.null(GS_all)) {
      for (trait in colnames(GS_all)) {
        hub_df[[paste0("GS_", trait)]] <- GS_all[module_proteins, trait]
      }
    }
    
    # Sort by kME (descending)
    hub_df <- hub_df[order(-hub_df$kME), ]
    
    # Mark top 10 as hub genes
    hub_df$Is_Hub <- FALSE
    hub_df$Is_Hub[1:min(10, nrow(hub_df))] <- TRUE
    
    hub_genes_list[[module_color]] <- hub_df
  }
}

# Combine all hub gene results
all_hub_genes <- do.call(rbind, hub_genes_list)
rownames(all_hub_genes) <- NULL

# Save hub genes
write.csv(all_hub_genes,
          file.path(output_dir, "module_hub_genes.csv"),
          row.names = FALSE)

cat("â?Hub gene identification completed\n")
cat("  Results saved to: module_hub_genes.csv\n\n")

# Print top hub genes for each module
cat("Top 5 hub genes per module (by kME):\n")
for (module_color in unique(original_colors)) {
  if (module_color == "grey") next
  
  module_name <- color_to_M[module_color]
  module_hubs <- all_hub_genes[all_hub_genes$Module == module_color, ]
  
  if (!is.null(module_hubs) && is.data.frame(module_hubs) && nrow(module_hubs) > 0) {
    top_hubs <- head(module_hubs, 5)
    cat("\n  ", module_name, ":\n", sep="")
    for (i in 1:nrow(top_hubs)) {
      cat(sprintf("    %d. %s (kME=%.3f)\n", i, top_hubs$Protein[i], top_hubs$kME[i]))
    }
  }
}

cat("\n\n================================================================================\n")
cat("â?ALL ANALYSES COMPLETE\n")
cat("================================================================================\n\n")
