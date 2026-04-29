#!/usr/bin/env Rscript
# Consensus WGCNA Analysis - New Strategy
# Four separate WGCNA analyses with preservation testing

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(impute)
  library(BRETIGEA)
})

allowWGCNAThreads()
options(stringsAsFactors = FALSE)

cat("================================================================================\n")
cat("NEW CONSENSUS WGCNA STRATEGY - 4 SEPARATE ANALYSES\n")
cat("================================================================================\n\n")

# Define 6 analysis scenarios based on diagnosis groups
SCENARIOS <- list(
  list(name = "LOAD",
       studies = c("study_1", "study_11"),
       diagnoses = list(study_1 = c("CN", "LOAD"), study_11 = c("CN", "LOAD")),
       is_consensus = TRUE),
  list(name = "EODSD",
       studies = c("study_4"),
       diagnoses = list(study_4 = c("CN", "EODSD")),
       is_consensus = FALSE),
  list(name = "EOFTD",
       studies = c("study_6"),
       diagnoses = list(study_6 = c("CN", "EOFTD")),
       is_consensus = FALSE),
  list(name = "LOFTD",
       studies = c("study_6"),
       diagnoses = list(study_6 = c("CN", "LOFTD")),
       is_consensus = FALSE),
  list(name = "EODLB",
       studies = c("study_9"),
       diagnoses = list(study_9 = c("CN", "EODLB")),
       is_consensus = FALSE),
  list(name = "LODLB",
       studies = c("study_9"),
       diagnoses = list(study_9 = c("CN", "LODLB")),
       is_consensus = FALSE)
)

output_base <- "wgcna_consensus_new"
if (!dir.exists(output_base)) {
  dir.create(output_base, recursive = TRUE)
}

# Load data
cat("Loading expression data...\n")
data_path <- "combine/combined_expression_matrices.csv"
df_all <- fread(data_path, data.table = FALSE)
cat("  Loaded", nrow(df_all), "samples\n\n")

metadata_cols <- c("Study", "Batch", "GUID", "Age", "Sex", "Diagnosis_Derived")
biomarker_cols <- c("Cognitive Score", "AB42", "tTau", "pTau", "pTau181", 
                    "AB42/pTau", "AB40", "NEFL", "YKL40", "pTau217", "pTau231")
available_metadata <- metadata_cols[metadata_cols %in% colnames(df_all)]
available_biomarkers <- biomarker_cols[biomarker_cols %in% colnames(df_all)]
all_non_protein_cols <- c(available_metadata, available_biomarkers)
protein_cols <- setdiff(colnames(df_all), all_non_protein_cols)

# Function: Process study with stratified imputation and diagnosis filtering
process_study <- function(study_name, df_all, protein_cols, available_metadata, available_biomarkers, diagnoses_to_keep = NULL) {
  cat("Processing", study_name, "...\n")
  df_study <- df_all[df_all$Study == study_name, ]
  if (nrow(df_study) == 0) return(NULL)
  
  # Filter by diagnosis if specified
  if (!is.null(diagnoses_to_keep)) {
    df_study <- df_study[df_study$Diagnosis_Derived %in% diagnoses_to_keep, ]
    cat("  Filtered to diagnoses:", paste(diagnoses_to_keep, collapse = ", "), "\n")
    if (nrow(df_study) == 0) {
      cat("  WARNING: No samples after diagnosis filtering\n")
      return(NULL)
    }
  }
  
  phenotype_df <- df_study[, available_metadata, drop = FALSE]
  expression_df <- df_study[, protein_cols, drop = FALSE]
  biomarker_df <- df_study[, available_biomarkers, drop = FALSE]
  
  na_counts <- colSums(is.na(expression_df))
  valid_proteins <- names(na_counts)[na_counts < nrow(expression_df)]
  expression_df <- expression_df[, valid_proteins, drop = FALSE]
  
  missing_pct <- rowSums(is.na(expression_df)) / ncol(expression_df)
  valid_samples <- missing_pct < 0.5
  expression_df <- expression_df[valid_samples, , drop = FALSE]
  biomarker_df <- biomarker_df[valid_samples, , drop = FALSE]
  phenotype_df <- phenotype_df[valid_samples, , drop = FALSE]
  
  cat("  Samples:", nrow(expression_df), "\n")
  
  # Stratified KNN imputation by diagnosis subtype
  diagnosis_vec <- phenotype_df$Diagnosis_Derived
  unique_diagnoses <- unique(diagnosis_vec[!is.na(diagnosis_vec)])
  
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
  
  # Z-normalize biomarkers
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
  
  sample_names <- paste(study_name, seq_len(nrow(phenotype_df)), sep = ".")
  rownames(expression_df) <- sample_names
  rownames(biomarker_df) <- sample_names
  rownames(phenotype_df) <- sample_names
  
  return(list(expression = expression_df, biomarker = biomarker_df, phenotype = phenotype_df))
}

# Function: Prepare traits
prepare_traits <- function(phenotype_df, biomarker_df, study_name) {
  trait_data <- data.frame(row.names = rownames(phenotype_df))
  diagnosis <- phenotype_df$Diagnosis_Derived
  
  cn_vec <- as.numeric(diagnosis %in% c("CN", "SCD"))
  if (sum(cn_vec) >= 5 && sum(cn_vec) < (length(cn_vec) - 5)) {
    trait_data$CN <- cn_vec
  }
  
  cn_samples <- diagnosis %in% c("CN", "SCD")
  
  disease_types <- c("EOAD", "EOFTD", "EOOD", "EODSD", "EODLB",
                     "LOAD", "LOFTD", "LOOD", "LODLB")
  for (dx_type in disease_types) {
    dx_samples <- diagnosis == dx_type
    if (sum(dx_samples) >= 5) {
      dx_vs_cn <- ifelse(dx_samples, 1, ifelse(cn_samples, 0, NA))
      if (sum(!is.na(dx_vs_cn)) >= 10 && var(dx_vs_cn, na.rm = TRUE) > 0.01) {
        trait_data[[paste0(dx_type, "_vs_CN")]] <- dx_vs_cn
      }
    }
  }
  
  age <- as.numeric(phenotype_df$Age)
  age[is.na(age)] <- mean(age, na.rm = TRUE)
  trait_data$Age <- age
  
  sex <- phenotype_df$Sex
  sex_vec <- as.numeric(sex == "Male")
  sex_vec[is.na(sex_vec)] <- 0.5
  if (var(sex_vec, na.rm = TRUE) > 0.01) {
    trait_data$Sex <- sex_vec
  }
  
  # Add biomarkers with study-specific cognitive score names
  biomarker_names <- colnames(biomarker_df)
  for (i in seq_along(biomarker_names)) {
    biomarker <- biomarker_names[i]
    biomarker_values <- biomarker_df[, biomarker]
    
    # Rename Cognitive Score based on study
    display_name <- biomarker
    if (biomarker == "Cognitive Score") {
      if (study_name == "study_1") {
        display_name <- "MoCA"
      } else if (study_name == "study_7") {
        display_name <- "MMSE"
      } else {
        # Other studies don't have cognitive scores, skip
        next
      }
    }
    
    if (sum(!is.na(biomarker_values)) >= 10 && var(biomarker_values, na.rm = TRUE) > 0.01) {
      trait_data[[display_name]] <- biomarker_values
    }
  }
  
  return(trait_data)
}

# Function: Calculate ORA with consensus modules
calculate_ora_with_consensus <- function(new_module_colors, new_proteins, consensus_file) {
  # Load consensus module assignments
  if (!file.exists(consensus_file)) {
    cat("  WARNING: Consensus module file not found\n")
    return(NULL)
  }
  
  consensus_modules <- read.csv(consensus_file, stringsAsFactors = FALSE)
  
  # Get all proteins in universe (union of new and consensus)
  all_proteins_universe <- unique(c(new_proteins, consensus_modules$Protein))
  universe_size <- length(all_proteins_universe)
  
  # Perform ORA for each new module against each consensus module
  ora_results <- list()
  
  for (new_mod in unique(new_module_colors)) {
    # Include grey module in ORA analysis
    
    new_mod_proteins <- new_proteins[new_module_colors == new_mod]
    new_mod_size <- length(new_mod_proteins)
    
    for (cons_mod in unique(consensus_modules$Module)) {
      # Include grey module in consensus comparison
      
      cons_mod_proteins <- consensus_modules$Protein[consensus_modules$Module == cons_mod]
      cons_mod_size <- length(cons_mod_proteins)
      
      # Calculate overlap
      overlap <- sum(new_mod_proteins %in% cons_mod_proteins)
      
      # Hypergeometric test
      # q = overlap - 1, m = consensus module size, n = universe - consensus size, k = new module size
      p_value <- phyper(overlap - 1, cons_mod_size, universe_size - cons_mod_size, new_mod_size, lower.tail = FALSE)
      
      # Calculate odds ratio
      a <- overlap
      b <- new_mod_size - overlap
      c <- cons_mod_size - overlap
      d <- universe_size - new_mod_size - cons_mod_size + overlap
      
      if (b > 0 && c > 0) {
        odds_ratio <- (a * d) / (b * c)
      } else {
        odds_ratio <- Inf
      }
      
      result_key <- paste(new_mod, cons_mod, sep = "_vs_")
      ora_results[[result_key]] <- data.frame(
        New_Module = new_mod,
        Consensus_Module = cons_mod,
        Overlap = overlap,
        New_Module_Size = new_mod_size,
        Consensus_Module_Size = cons_mod_size,
        Universe_Size = universe_size,
        P_value = p_value,
        Odds_Ratio = odds_ratio,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(ora_results) == 0) return(NULL)
  
  ora_df <- do.call(rbind, ora_results)
  rownames(ora_df) <- NULL
  ora_df$FDR <- p.adjust(ora_df$P_value, method = "fdr")
  
  return(ora_df)
}

cat("\nStarting analysis loop...\n\n")

# Load consensus module assignments for ORA
consensus_module_file <- "wgcna_consensus_main/consensus_module_assignments.csv"

# MAIN ANALYSIS LOOP
for (scenario_idx in 1:length(SCENARIOS)) {
  scenario <- SCENARIOS[[scenario_idx]]
  cat("\n================================================================================\n")
  cat("SCENARIO", scenario_idx, ":", scenario$name, "\n")
  cat("Studies:", paste(scenario$studies, collapse = ", "), "\n")
  cat("================================================================================\n\n")
  
  output_dir <- file.path(output_base, scenario$name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Process training studies with diagnosis filtering
  cat("Processing training studies...\n")
  study_data_list <- list()
  for (study in scenario$studies) {
    diagnoses_for_study <- scenario$diagnoses[[study]]
    result <- process_study(study, df_all, protein_cols, available_metadata, available_biomarkers, diagnoses_for_study)
    if (!is.null(result)) {
      study_data_list[[study]] <- result
    }
  }
  
  if (length(study_data_list) == 0) {
    cat("ERROR: No training studies loaded\n")
    next
  }
  
  # Find common proteins
  protein_presence <- lapply(study_data_list, function(x) colnames(x$expression))
  common_proteins <- Reduce(intersect, protein_presence)
  cat("Common proteins:", length(common_proteins), "\n\n")
  
  for (study in names(study_data_list)) {
    study_data_list[[study]]$expression <- study_data_list[[study]]$expression[, common_proteins, drop = FALSE]
  }
  
  # Prepare multiExpr and multiTrait
  multiExpr <- list()
  multiTrait <- list()
    
    for (study in names(study_data_list)) {
    multiExpr[[study]] <- list(data = as.data.frame(
      study_data_list[[study]]$expression[, common_proteins, drop = FALSE]
    ))
    multiTrait[[study]] <- prepare_traits(
      study_data_list[[study]]$phenotype,
      study_data_list[[study]]$biomarker,
      study
    )
    cat("  ", study, ":", nrow(multiExpr[[study]]$data), "samples x",
        ncol(multiExpr[[study]]$data), "proteins |  Traits:",
        ncol(multiTrait[[study]]), "\n")
  }
  cat("\n")
  
  # Choose soft-thresholding power (per-study, then take consensus median)
  cat("Choosing soft-thresholding power...\n")
  powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
  powerTables <- list()
  for (study in names(multiExpr)) {
    sft <- pickSoftThreshold(multiExpr[[study]]$data, powerVector = powers, verbose = 0)
    powerTables[[study]] <- sft$fitIndices
    opt <- sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
    cat("  ", study, "optimal power:", opt,
        "(R^2 =", round(max(sft$fitIndices$SFT.R.sq), 3), ")\n")
  }
  optimal_powers <- sapply(powerTables, function(x) x$Power[which.max(x$SFT.R.sq)])
  consensus_power <- round(median(optimal_powers))
  cat("  Consensus power:", consensus_power, "\n\n")
  
  if (scenario$is_consensus && length(multiExpr) > 1) {
    # ---- True consensus WGCNA (blockwiseConsensusModules) ----
    cat("Running blockwiseConsensusModules...\n")
    # Verify data structure
    checkSets(multiExpr)
    
    consensusNet <- blockwiseConsensusModules(
      multiExpr,
      power             = consensus_power,
      minModuleSize     = 10,
      deepSplit         = 4,
      pamRespectsDendro = TRUE,
      mergeCutHeight    = 0.1,
      numericLabels     = FALSE,
      minKMEtoStay      = 0.30,
      reassignThreshold = 0.05,
      saveTOMs          = FALSE,
      verbose           = 3,
      maxBlockSize      = 10000,
      networkType       = "signed",
      TOMType           = "signed",
      corType           = "bicor",
      consensusQuantile = 0.3
    )
    dynamicColors <- consensusNet$colors
    names(dynamicColors) <- colnames(multiExpr[[1]]$data)
    cat("  Identified", length(unique(dynamicColors)), "consensus modules\n\n")
    
    # kME-based reassignment (consistent with consensus_wgcna.R)
    datExpr_ref <- multiExpr[[1]]$data
    MEList_ref  <- moduleEigengenes(datExpr_ref, colors = dynamicColors)
    MEs_ref     <- MEList_ref$eigengenes
    kME_ref     <- bicor(datExpr_ref, MEs_ref, use = "pairwise.complete.obs")
    
    for (i in seq_along(dynamicColors)) {
      cur_col <- dynamicColors[i]
      if (cur_col == "grey") next
      cur_ME  <- paste0("ME", cur_col)
      if (cur_ME %in% colnames(kME_ref)) {
        cur_kME <- kME_ref[i, cur_ME]
        max_kME <- max(kME_ref[i, ], na.rm = TRUE)
        best_ME <- colnames(kME_ref)[which.max(kME_ref[i, ])]
        best_col <- gsub("ME", "", best_ME)
        if (!is.na(max_kME) && !is.na(cur_kME) &&
            max_kME > cur_kME + 0.1 && max_kME > 0.3) {
          dynamicColors[i] <- best_col
        }
      }
    }
    
    # Merge similar modules (ME correlation > 0.9)
    merge_res     <- mergeCloseModules(datExpr_ref, dynamicColors, cutHeight = 0.1, verbose = 3)
    dynamicColors <- merge_res$colors
    cat("  After kME-reassign + merge:", length(unique(dynamicColors)), "modules\n\n")
    
  } else {
    # ---- Single-study WGCNA (blockwiseModules) ----
  cat("Running blockwiseModules...\n")
  net <- blockwiseModules(
      multiExpr[[1]]$data,
      power             = consensus_power,
      minModuleSize     = 10,
      deepSplit         = 4,
    pamRespectsDendro = TRUE,
      mergeCutHeight    = 0.1,
      numericLabels     = FALSE,
      minKMEtoStay      = 0.30,
    reassignThreshold = 0.05,
      saveTOMs          = FALSE,
      verbose           = 3,
      maxBlockSize      = 10000,
      networkType       = "signed",
      TOMType           = "signed",
      corType           = "bicor"
    )
  dynamicColors <- net$colors
    names(dynamicColors) <- colnames(multiExpr[[1]]$data)
    cat("  Identified", length(unique(dynamicColors)), "modules\n\n")
  }
  
  # Module cleaning (kME < 0.3 -> grey; grey with kME > 0.25 -> reassign)
  cat("Module cleaning...\n")
  datExpr_clean <- multiExpr[[1]]$data
  MEList <- moduleEigengenes(datExpr_clean, colors = dynamicColors)
  MEs_clean <- MEList$eigengenes
  kME_matrix <- bicor(datExpr_clean, MEs_clean, use = "pairwise.complete.obs")
  cleaned_colors <- dynamicColors
  
  for (i in seq_along(cleaned_colors)) {
    if (cleaned_colors[i] == "grey") next
    ME_col <- paste0("ME", cleaned_colors[i])
    if (ME_col %in% colnames(kME_matrix)) {
      gene_kME <- kME_matrix[i, ME_col]
      if (!is.na(gene_kME) && gene_kME < 0.3) cleaned_colors[i] <- "grey"
    }
  }
  grey_indices <- which(cleaned_colors == "grey")
  for (i in grey_indices) {
    gene_kMEs <- kME_matrix[i, !grepl("MEgrey", colnames(kME_matrix))]
    if (length(gene_kMEs) > 0) {
      max_kME <- max(gene_kMEs, na.rm = TRUE)
      if (!is.na(max_kME) && max_kME > 0.25) {
        best_ME  <- names(gene_kMEs)[which.max(gene_kMEs)]
        cleaned_colors[i] <- gsub("ME", "", best_ME)
      }
    }
  }
  dynamicColors <- cleaned_colors
  cat("  Final modules:", length(unique(dynamicColors)), "\n\n")
  
  # Rename modules by size (M1 = largest, M2 = second largest, etc.)
  cat("Renaming modules by size...\n")
  module_sizes <- table(dynamicColors)
  module_sizes <- sort(module_sizes[names(module_sizes) != "grey"], decreasing = TRUE)
  
  # Create mapping: color -> M number with color name
  color_to_M <- setNames(paste0("M", seq_along(module_sizes), "_", names(module_sizes)), 
                         names(module_sizes))
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
  
  # Recalculate MEs per study (consistent with consensus_wgcna.R)
  MEs_list <- list()
  for (study in names(multiExpr)) {
    MEList_s <- moduleEigengenes(multiExpr[[study]]$data, colors = original_colors)
    MEs_s    <- MEList_s$eigengenes
    new_colnames <- colnames(MEs_s)
    for (color in names(color_to_M)) {
      old_name <- paste0("ME", color)
      if (old_name %in% new_colnames)
        new_colnames[new_colnames == old_name] <- color_to_M[color]
    }
    colnames(MEs_s) <- new_colnames
    MEs_list[[study]] <- MEs_s
  }
  # Calculate module-trait correlations
  cat("Calculating module-trait correlations...\n")
  all_correlations <- list()
  all_pvalues <- list()
  
  for (study in names(multiExpr)) {
    ME <- MEs_list[[study]]
    traits <- multiTrait[[study]]
    common_samples <- intersect(rownames(ME), rownames(traits))
    ME <- ME[common_samples, , drop = FALSE]
    traits <- traits[common_samples, , drop = FALSE]
    
    moduleTraitCor <- bicor(ME, traits, use = "pairwise.complete.obs")
    moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(ME))
    
    all_correlations[[study]] <- moduleTraitCor
    all_pvalues[[study]] <- moduleTraitPvalue
  }
  
  # Consensus correlations
  all_trait_names <- unique(unlist(lapply(all_correlations, colnames)))
  all_module_names <- rownames(all_correlations[[1]])
  
  consensus_cor <- matrix(NA, nrow = length(all_module_names), ncol = length(all_trait_names))
  rownames(consensus_cor) <- all_module_names
  colnames(consensus_cor) <- all_trait_names
  
  consensus_pval <- matrix(NA, nrow = length(all_module_names), ncol = length(all_trait_names))
  rownames(consensus_pval) <- all_module_names
  colnames(consensus_pval) <- all_trait_names
  
  for (module in all_module_names) {
    for (trait in all_trait_names) {
      cors <- c()
      pvals <- c()
      n_samples <- c()
      
      for (study in names(all_correlations)) {
        if (trait %in% colnames(all_correlations[[study]])) {
          cors <- c(cors, all_correlations[[study]][module, trait])
          pvals <- c(pvals, all_pvalues[[study]][module, trait])
          n_samples <- c(n_samples, nrow(multiExpr[[study]]$data))
        }
      }
      
      if (length(cors) >= 1) {
        total_n <- sum(n_samples)
        weights <- n_samples / total_n
        consensus_cor[module, trait] <- sum(cors * weights, na.rm = TRUE)
        
        # Filter out NA p-values before combining
        valid_pvals <- !is.na(pvals)
        pvals_valid <- pvals[valid_pvals]
        
        if (length(pvals_valid) > 1) {
          # Limit p-values to avoid log(0) or numerical issues
          pvals_safe <- pmax(pvals_valid, 1e-300)
          chi_sq <- -2 * sum(log(pvals_safe))
          consensus_pval[module, trait] <- pchisq(chi_sq, df = 2 * length(pvals_safe), lower.tail = FALSE)
        } else if (length(pvals_valid) == 1) {
          consensus_pval[module, trait] <- pvals_valid[1]
        } else {
          consensus_pval[module, trait] <- NA
        }
      }
    }
  }
  
  # Calculate preservation
  cat("Calculating ORA with consensus modules...\n")
  ora_results <- calculate_ora_with_consensus(original_colors, common_proteins, consensus_module_file)
  
  if (!is.null(ora_results)) {
    # Add M names to ORA results
    ora_results$New_Module_Name <- color_to_M[ora_results$New_Module]
    
    # Reorder columns
    ora_results <- ora_results[, c("New_Module", "New_Module_Name", "Consensus_Module", 
                                   "Overlap", "New_Module_Size", "Consensus_Module_Size", 
                                   "Universe_Size", "P_value", "FDR", "Odds_Ratio")]
    
    write.csv(ora_results, 
              file.path(output_dir, paste0(scenario$name, "_ORA_vs_consensus.csv")),
              row.names = FALSE)
    cat("  Saved ORA results\n")
  }
  
  # Cell type enrichment
  cat("Cell type enrichment analysis...\n")
  tryCatch({
    data("markers_df_human_brain", package = "BRETIGEA")
    cell_types_map <- c(neurons = "neu", astrocytes = "ast", oligodendrocytes = "oli",
                        microglia = "mic", endothelial = "end", OPCs = "opc")
    
    bretigea_markers <- list()
    for (cell_type in names(cell_types_map)) {
      cell_abbr <- cell_types_map[cell_type]
      markers <- markers_df_human_brain[markers_df_human_brain$cell == cell_abbr, ]
      bretigea_markers[[cell_type]] <- head(markers$markers, 200)
    }
    
    all_genes <- sapply(strsplit(common_proteins, "\\|"), function(x) x[1])
    fisher_results <- list()
    
    for (module_color in unique(original_colors)) {
      # Include grey module in cell type enrichment analysis
      module_proteins <- common_proteins[original_colors == module_color]
      module_genes <- sapply(strsplit(module_proteins, "\\|"), function(x) x[1])
      
      for (cell_type in names(bretigea_markers)) {
        markers <- bretigea_markers[[cell_type]]
        in_both <- sum(module_genes %in% markers)
        in_module_only <- length(module_genes) - in_both
        in_markers_only <- sum(markers %in% all_genes) - in_both
        background <- length(all_genes) - length(module_genes) - in_markers_only
        
        contingency <- matrix(c(in_both, in_module_only, in_markers_only, background), nrow = 2)
        fisher_test <- fisher.test(contingency, alternative = "greater")
        
        result_key <- paste(module_color, cell_type, sep = "_")
        fisher_results[[result_key]] <- data.frame(
          Module = module_color,
          Module_Name = color_to_M[module_color],
          Cell_Type = cell_type,
          Overlap = in_both,
          OR = fisher_test$estimate,
          P_value = fisher_test$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
    
    fisher_df <- do.call(rbind, fisher_results)
    fisher_df$FDR <- p.adjust(fisher_df$P_value, method = "fdr")
    
    # Reorder columns
    fisher_df <- fisher_df[, c("Module", "Module_Name", "Cell_Type", "Overlap", "OR", "P_value", "FDR")]
    
    write.csv(fisher_df, 
              file.path(output_dir, paste0(scenario$name, "_celltype.csv")),
              row.names = FALSE)
  }, error = function(e) {
    cat("  Cell type analysis failed\n")
  })
  
  # Generate heatmap
  cat("Generating heatmap...\n")
  valid_traits <- colSums(!is.na(consensus_cor)) > 0
  cor_plot <- consensus_cor[, valid_traits, drop = FALSE]
  pval_plot <- consensus_pval[, valid_traits, drop = FALSE]
  
  desired_order <- c("CN", "EOAD_vs_CN", "LOAD_vs_CN", "EODSD_vs_CN", "EODLB_vs_CN", "LODLB_vs_CN",
                     "EOFTD_vs_CN", "LOFTD_vs_CN", "EOOD_vs_CN", "LOOD_vs_CN",
                     "Age", "Sex", "MoCA", "MMSE", "AB42", "tTau", "pTau", "pTau181",
                     "AB42/pTau", "AB40", "NEFL", "YKL40", "pTau217", "pTau231")
  available_traits <- colnames(cor_plot)
  ordered_traits <- desired_order[desired_order %in% available_traits]
  remaining_traits <- setdiff(available_traits, ordered_traits)
  if (length(remaining_traits) > 0) {
    ordered_traits <- c(ordered_traits, remaining_traits)
  }
  
  cor_plot <- cor_plot[, ordered_traits, drop = FALSE]
  pval_plot <- pval_plot[, ordered_traits, drop = FALSE]
  
  combined_matrix <- cor_plot
  
  # Add cell type enrichment
  if (exists("fisher_df") && nrow(fisher_df) > 0) {
    cell_matrix <- matrix(NA, nrow = nrow(combined_matrix), ncol = length(unique(fisher_df$Cell_Type)))
    rownames(cell_matrix) <- rownames(combined_matrix)
    colnames(cell_matrix) <- paste0("Cell_", unique(fisher_df$Cell_Type))
    
    for (i in 1:nrow(fisher_df)) {
      module_name <- fisher_df$Module_Name[i]
      cell_type <- fisher_df$Cell_Type[i]
      if (module_name %in% rownames(cell_matrix)) {
        log_fdr <- -log10(fisher_df$FDR[i])
        cell_matrix[module_name, paste0("Cell_", cell_type)] <- min(log_fdr / 10, 1)
      }
    }
    
    combined_matrix <- cbind(combined_matrix, cell_matrix)
  }
  
  # Create text matrix
  textMatrix <- matrix("", nrow = nrow(combined_matrix), ncol = ncol(combined_matrix))
  
  for (i in 1:nrow(cor_plot)) {
    for (j in 1:ncol(cor_plot)) {
      if (!is.na(cor_plot[i, j])) {
        pval <- pval_plot[i, j]
        if (!is.na(pval) && pval < 0.01) {
          textMatrix[i, j] <- paste0(signif(cor_plot[i, j], 2), "\n(<0.01)")
        } else if (!is.na(pval)) {
          textMatrix[i, j] <- paste0(signif(cor_plot[i, j], 2), "\n(", signif(pval, 1), ")")
        } else {
          textMatrix[i, j] <- signif(cor_plot[i, j], 2)
        }
      }
    }
  }
  
  png(file.path(output_dir, paste0(scenario$name, "_heatmap.png")), 
      width = 16, height = 10, units = "in", res = 300)
  
  par(mar = c(10, 10, 3, 3))
  labeledHeatmap(Matrix = combined_matrix,
                xLabels = colnames(combined_matrix),
                yLabels = rownames(combined_matrix),
                ySymbols = rownames(combined_matrix),
                colorLabels = FALSE,
                colors = blueWhiteRed(50),
                textMatrix = textMatrix,
                setStdMargins = FALSE,
                cex.text = 0.4,
                zlim = c(-1, 1),
                main = paste0(scenario$name, ": Module-Trait Relationships"))
  
  dev.off()
  
  write.csv(cor_plot, 
            file.path(output_dir, paste0(scenario$name, "_correlations.csv")))
  write.csv(pval_plot, 
            file.path(output_dir, paste0(scenario$name, "_pvalues.csv")))
  
  # Save module assignments with M names
  module_assignments <- data.frame(
    Protein = common_proteins,
    Module = original_colors,
    Module_Name = color_to_M[original_colors],
    stringsAsFactors = FALSE
  )
  
  write.csv(module_assignments,
            file.path(output_dir, paste0(scenario$name, "_module_assignments.csv")),
            row.names = FALSE)
  
  cat("Scenario", scenario_idx, "complete\n\n")
}

cat("================================================================================\n")
cat("ALL ANALYSES COMPLETE\n")
cat("================================================================================\n\n")
cat("Output directory:", output_base, "\n")
cat("Generated files per scenario:\n")
cat("  - _correlations.csv\n")
cat("  - _pvalues.csv\n")
cat("  - _ORA_vs_consensus.csv\n")
cat("  - _celltype.csv\n")
cat("  - _module_assignments.csv\n")
cat("  - _heatmap.png\n\n")
