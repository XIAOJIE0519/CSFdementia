#!/usr/bin/env Rscript
cat("================================================================================\n")
cat("WGCNA MODULE GO/KEGG ENRICHMENT ANALYSIS FOR MULTIPLE DIAGNOSES\n")
cat("================================================================================\n\n")

# Load required libraries
cat("Loading libraries...\n")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(tidyverse)
})

cat("All libraries loaded successfully\n\n")

# 定义输入和输出路径
base_dir <- "F:/1a-EOD-CSF-protein/wgcna_consensus_new"
output_file <- "F:/1a-EOD-CSF-protein/wgcna_consensus_new/GOKEGG.csv"

# 获取所有诊断文件夹
diagnosis_folders <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
cat("Found diagnosis folders:", paste(diagnosis_folders, collapse=", "), "\n\n")

# 存储所有富集结果
all_results <- list()

# 对每个诊断进行分析
cat("================================================================================\n")
cat("PROCESSING EACH DIAGNOSIS\n")
cat("================================================================================\n\n")

for (diagnosis in diagnosis_folders) {
  cat("Processing diagnosis:", diagnosis, "\n")
  
  # 构建module_assignments文件路径
  module_file <- file.path(base_dir, diagnosis, paste0(diagnosis, "_module_assignments.csv"))
  
  if (!file.exists(module_file)) {
    cat("  WARNING: Module assignment file not found, skipping\n\n")
    next
  }
  
  # 读取模块分配
  module_df <- read.csv(module_file, stringsAsFactors = FALSE)
  cat("  Loaded", nrow(module_df), "proteins\n")
  
  # 获取所有模块
  modules <- unique(module_df$Module)
  modules <- modules[modules != "grey"]  # 排除grey模块
  cat("  Modules found:", length(modules), "\n")
  
  # 对每个模块进行富集分析
  for (module_color in modules) {
    cat("    Module:", module_color, "\n")
    
    # 提取该模块的蛋白
    module_proteins <- module_df$Protein[module_df$Module == module_color]
    
    if (length(module_proteins) < 3) {
      cat("      WARNING: Too few proteins (", length(module_proteins), "), skipping\n")
      next
    }
    
    cat("      Proteins:", length(module_proteins), "\n")
    
    # 转换为Entrez ID
    gene_to_entrez <- NULL
    tryCatch({
      gene_to_entrez <- suppressMessages(bitr(
        module_proteins,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
      ))
    }, error = function(e) {
      cat("      ERROR: Failed to map genes to Entrez\n")
    })
    
    if (is.null(gene_to_entrez) || nrow(gene_to_entrez) < 3) {
      cat("      WARNING: Too few genes mapped to Entrez\n")
      next
    }
    
    gene_entrez <- gene_to_entrez$ENTREZID
    cat("      Genes for enrichment:", length(gene_entrez), "\n")
    
    # 执行富集分析
    top_bp <- NA
    top_mf <- NA
    top_cc <- NA
    top_kegg <- NA
    
    # GO Biological Process
    tryCatch({
      ego_bp <- enrichGO(
        gene = gene_entrez,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        readable = TRUE
      )
      
      if (!is.null(ego_bp) && nrow(ego_bp@result) > 0) {
        top_bp <- ego_bp@result$Description[1]
        cat("        GO:BP:", nrow(ego_bp@result), "terms\n")
      }
    }, error = function(e) {})
    
    # GO Molecular Function
    tryCatch({
      ego_mf <- enrichGO(
        gene = gene_entrez,
        OrgDb = org.Hs.eg.db,
        ont = "MF",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        readable = TRUE
      )
      
      if (!is.null(ego_mf) && nrow(ego_mf@result) > 0) {
        top_mf <- ego_mf@result$Description[1]
        cat("        GO:MF:", nrow(ego_mf@result), "terms\n")
      }
    }, error = function(e) {})
    
    # GO Cellular Component
    tryCatch({
      ego_cc <- enrichGO(
        gene = gene_entrez,
        OrgDb = org.Hs.eg.db,
        ont = "CC",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        readable = TRUE
      )
      
      if (!is.null(ego_cc) && nrow(ego_cc@result) > 0) {
        top_cc <- ego_cc@result$Description[1]
        cat("        GO:CC:", nrow(ego_cc@result), "terms\n")
      }
    }, error = function(e) {})
    
    # KEGG Pathway
    tryCatch({
      ekegg <- enrichKEGG(
        gene = gene_entrez,
        organism = 'hsa',
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05
      )
      
      if (!is.null(ekegg) && nrow(ekegg@result) > 0) {
        ekegg_readable <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
        top_kegg <- ekegg_readable@result$Description[1]
        cat("        KEGG:", nrow(ekegg@result), "pathways\n")
      }
    }, error = function(e) {})
    
    # 存储结果
    result_key <- paste(diagnosis, module_color, sep="_")
    all_results[[result_key]] <- data.frame(
      Diagnosis = diagnosis,
      Module = module_color,
      N_Proteins = length(module_proteins),
      Top_GO_BP = ifelse(is.na(top_bp), "None", top_bp),
      Top_GO_MF = ifelse(is.na(top_mf), "None", top_mf),
      Top_GO_CC = ifelse(is.na(top_cc), "None", top_cc),
      Top_KEGG = ifelse(is.na(top_kegg), "None", top_kegg),
      stringsAsFactors = FALSE
    )
  }
  
  cat("\n")
}

# 合并所有结果
if (length(all_results) > 0) {
  enrichment_summary <- do.call(rbind, all_results)
  rownames(enrichment_summary) <- NULL
  
  # 保存CSV结果
  write.csv(enrichment_summary, output_file, row.names = FALSE)
  cat("================================================================================\n")
  cat("ANALYSIS COMPLETE\n")
  cat("================================================================================\n")
  cat("Total results:", nrow(enrichment_summary), "module-diagnosis combinations\n")
  cat("Output file:", output_file, "\n")
  cat("================================================================================\n")
} else {
  cat("ERROR: No enrichment results generated\n")
}
