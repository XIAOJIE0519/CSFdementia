#!/usr/bin/env Rscript
cat("================================================================================\n")
cat("WGCNA MODULE GO/KEGG ENRICHMENT ANALYSIS\n")
cat("================================================================================\n\n")

# Load required libraries
cat("Loading libraries...\n")
library(clusterProfiler, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)
library(tidyverse, quietly = TRUE)

cat("All libraries loaded successfully\n\n")

# 创建输出目录
output_dir <- "F:/1a-EOD-CSF-protein/enrichment_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n\n")
}

# 读取模块分配文件
module_file <- "F:/1a-EOD-CSF-protein/wgcna_consensus_main/consensus_module_assignments.csv"

if (!file.exists(module_file)) {
  stop("Module assignment file not found: ", module_file)
}

module_df <- read.csv(module_file, stringsAsFactors = FALSE)
cat("Loaded module assignments:", nrow(module_df), "proteins\n")
cat("Modules found:", length(unique(module_df$Module)), "\n\n")

# 读取pvalues文件以获取模块顺序
pval_file <- "F:/1a-EOD-CSF-protein/wgcna_consensus_main/consensus_module_trait_pvalues.csv"
if (!file.exists(pval_file)) {
  stop("Pvalue file not found: ", pval_file)
}

pval_df <- read.csv(pval_file, row.names = 1)
module_order_full <- rownames(pval_df)
cat("Module order from pvalues file:", length(module_order_full), "modules\n")

# 提取颜色名（从M5_pink提取pink）
module_order <- sapply(module_order_full, function(x) {
  parts <- strsplit(x, "_")[[1]]
  if (length(parts) >= 2) {
    return(paste(parts[-1], collapse="_"))
  } else {
    return(x)
  }
})
names(module_order) <- NULL
cat("Module color order:", paste(module_order, collapse=", "), "\n\n")

# 存储所有富集结果
all_enrichment_results <- list()

# 对每个模块进行富集分析
cat("================================================================================\n")
cat("PERFORMING ENRICHMENT ANALYSIS FOR EACH MODULE\n")
cat("================================================================================\n\n")

for (i in 1:length(module_order)) {
  module_color <- module_order[i]
  module_full_name <- module_order_full[i]
  
  cat("Processing module:", module_full_name, "(", module_color, ")\n")
  
  # 提取该模块的蛋白（使用颜色名匹配Module列）
  module_proteins <- module_df$Protein[module_df$Module == module_color]
  
  if (length(module_proteins) < 3) {
    cat("  WARNING: Too few proteins (", length(module_proteins), "), skipping\n\n")
    all_enrichment_results[[module_full_name]] <- data.frame(
      Module = module_full_name,
      N_Proteins = length(module_proteins),
      Top_GO_BP = NA,
      Top_GO_MF = NA,
      Top_GO_CC = NA,
      Top_KEGG = NA,
      stringsAsFactors = FALSE
    )
    next
  }
  
  cat("  Proteins in module:", length(module_proteins), "\n")
  
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
    cat("  ERROR: Failed to map genes to Entrez\n\n")
  })
  
  if (is.null(gene_to_entrez) || nrow(gene_to_entrez) < 3) {
    cat("  WARNING: Too few genes mapped to Entrez\n\n")
    all_enrichment_results[[module_full_name]] <- data.frame(
      Module = module_full_name,
      N_Proteins = length(module_proteins),
      Top_GO_BP = NA,
      Top_GO_MF = NA,
      Top_GO_CC = NA,
      Top_KEGG = NA,
      stringsAsFactors = FALSE
    )
    next
  }
  
  gene_entrez <- gene_to_entrez$ENTREZID
  cat("  Genes for enrichment:", length(gene_entrez), "\n")
  
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
      cat("    GO:BP:", nrow(ego_bp@result), "terms, Top:", top_bp, "\n")
    }
  }, error = function(e) {
    cat("    GO:BP: No significant terms\n")
  })
  
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
      cat("    GO:MF:", nrow(ego_mf@result), "terms, Top:", top_mf, "\n")
    }
  }, error = function(e) {
    cat("    GO:MF: No significant terms\n")
  })
  
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
      cat("    GO:CC:", nrow(ego_cc@result), "terms, Top:", top_cc, "\n")
    }
  }, error = function(e) {
    cat("    GO:CC: No significant terms\n")
  })
  
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
      cat("    KEGG:", nrow(ekegg@result), "pathways, Top:", top_kegg, "\n")
    }
  }, error = function(e) {
    cat("    KEGG: No significant pathways\n")
  })
  
  # 存储结果
  all_enrichment_results[[module_full_name]] <- data.frame(
    Module = module_full_name,
    N_Proteins = length(module_proteins),
    Top_GO_BP = ifelse(is.na(top_bp), "None", top_bp),
    Top_GO_MF = ifelse(is.na(top_mf), "None", top_mf),
    Top_GO_CC = ifelse(is.na(top_cc), "None", top_cc),
    Top_KEGG = ifelse(is.na(top_kegg), "None", top_kegg),
    stringsAsFactors = FALSE
  )
  
  cat("\n")
}

# 合并所有结果
enrichment_summary <- do.call(rbind, all_enrichment_results)
rownames(enrichment_summary) <- NULL

# 保存CSV结果
csv_output <- file.path(output_dir, "module_enrichment_summary.csv")
write.csv(enrichment_summary, csv_output, row.names = FALSE)
cat("Saved enrichment summary:", csv_output, "\n\n")

# 根据富集结果为每个模块命名
cat("================================================================================\n")
cat("GENERATING MODULE NAMES BASED ON ENRICHMENT\n")
cat("================================================================================\n\n")

# 定义命名函数
generate_module_name <- function(module, top_bp, top_mf, top_cc, top_kegg) {
  # 优先级：KEGG > GO:BP > GO:MF > GO:CC
  
  # 提取关键词的函数
  extract_keywords <- function(term) {
    if (is.na(term) || term == "None") return(NULL)
    
    # 移除常见的前缀和后缀
    term <- gsub("^(positive |negative |regulation of |cellular |biological |molecular )", "", term, ignore.case = TRUE)
    term <- gsub(" process$| pathway$| activity$| component$", "", term, ignore.case = TRUE)
    
    # 分割并取前2-3个关键词
    words <- strsplit(term, " ")[[1]]
    keywords <- words[1:min(3, length(words))]
    
    # 转换为简短名称
    name <- paste(keywords, collapse = "_")
    name <- gsub("[^A-Za-z0-9_]", "", name)
    
    return(name)
  }
  
  # 按优先级尝试命名
  if (!is.na(top_kegg) && top_kegg != "None") {
    name <- extract_keywords(top_kegg)
    if (!is.null(name)) return(name)
  }
  
  if (!is.na(top_bp) && top_bp != "None") {
    name <- extract_keywords(top_bp)
    if (!is.null(name)) return(name)
  }
  
  if (!is.na(top_mf) && top_mf != "None") {
    name <- extract_keywords(top_mf)
    if (!is.null(name)) return(name)
  }
  
  if (!is.na(top_cc) && top_cc != "None") {
    name <- extract_keywords(top_cc)
    if (!is.null(name)) return(name)
  }
  
  # 如果没有富集结果，使用模块颜色
  return(module)
}

# 为每个模块生成名称
module_names <- character(nrow(enrichment_summary))

for (i in 1:nrow(enrichment_summary)) {
  row <- enrichment_summary[i, ]
  module_names[i] <- generate_module_name(
    row$Module,
    row$Top_GO_BP,
    row$Top_GO_MF,
    row$Top_GO_CC,
    row$Top_KEGG
  )
}

# 创建模块名称映射（使用完整的模块名）
module_mapping <- data.frame(
  Module_Full = module_order_full,
  Module_Color = module_order,
  Module_Name = module_names,
  stringsAsFactors = FALSE
)

# 输出到module_enrich.txt
txt_output <- "F:/1a-EOD-CSF-protein/module_enrich.txt"
writeLines(paste0(module_mapping$Module_Full, "_", module_mapping$Module_Name), txt_output)

cat("Module names generated:\n")
for (i in 1:nrow(module_mapping)) {
  cat(sprintf("  %s_%s\n", 
              module_mapping$Module_Full[i], 
              module_mapping$Module_Name[i]))
}

cat("\nSaved module names to:", txt_output, "\n")

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n")
cat("Output files:\n")
cat("  1. Enrichment summary CSV:", csv_output, "\n")
cat("  2. Module names TXT:", txt_output, "\n")
cat("================================================================================\n")
