#!/usr/bin/env Rscript
cat("================================================================================\n")
cat("GOKEGG ENRICHMENT VISUALIZATION AND GSEA ANALYSIS\n")
cat("================================================================================\n\n")

# Load required libraries
cat("Loading libraries...\n")
library(clusterProfiler, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)
library(enrichplot, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(gground, quietly = TRUE)
library(ggprism, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(GseaVis, quietly = TRUE)

cat("All libraries loaded successfully\n\n")

# 创建输出目录
output_dir <- "F:/1a-EOD-CSF-protein/enrichment_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n\n")
}

# 定义颜色方案
pal <- c('#c3e1e6', '#f3dfb7', '#dcc6dc', '#96c38e')
# 定义更深的基因字颜色
pal_dark <- c('#5a9fb5', '#d4a03a', '#a888a8', '#5a8a5e')

# 定义六个meta结果文件
meta_files <- c(
  'EOAD_vs_CN',
  'LOAD_vs_CN',
  'EOAD_vs_LOAD',
  'EOD_vs_CN',
  'LOD_vs_CN',
  'EOD_vs_LOD'
)

# ========== 第一部分：GO/KEGG富集分析和绘图 ==========
cat("PART 1: GO/KEGG Enrichment Analysis and Visualization\n")
cat("----------------------------------------------------------------\n\n")

# 初始化存储所有组合富集结果的列表（供跨组合对比图使用）
all_enrichment_results <- list()

for (comp_name in meta_files) {
  cat("Processing", comp_name, "...\n")
  
  # 读取meta分析结果
  meta_file <- file.path("F:/1a-EOD-CSF-protein/meta", paste0(comp_name, ".csv"))
  
  if (!file.exists(meta_file)) {
    cat("  WARNING: Meta file not found:", meta_file, "\n\n")
    next
  }
  
  meta_df <- read.csv(meta_file, stringsAsFactors = FALSE)
  
  # 筛选FDR < 0.05的蛋白
  sig_proteins <- meta_df$Protein[meta_df$FDR_BH_Stratified < 0.05 & !is.na(meta_df$FDR_BH_Stratified)]
  
  if (length(sig_proteins) < 3) {
    cat("  WARNING: Too few significant proteins (FDR < 0.05):", length(sig_proteins), "\n\n")
    next
  }
  
  cat("  Significant proteins (FDR < 0.05):", length(sig_proteins), "\n")
  
  # 蛋白名已经是基因名，直接转换为Entrez ID
  gene_to_entrez <- NULL
  tryCatch({
    gene_to_entrez <- suppressMessages(bitr(
      sig_proteins,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    ))
  }, error = function(e) {
    cat("  ERROR: Failed to map genes to Entrez:", e$message, "\n\n")
  })
  
  if (is.null(gene_to_entrez) || nrow(gene_to_entrez) < 3) {
    cat("  WARNING: Too few genes mapped to Entrez\n\n")
    next
  }
  
  gene_entrez <- gene_to_entrez$ENTREZID
  cat("  Genes for enrichment:", length(gene_entrez), "\n")
  
  # 执行富集分析
  results <- list()
  
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
      results$GO_BP <- ego_bp@result
      cat("    GO:BP:", nrow(ego_bp@result), "terms\n")
    }
  }, error = function(e) {
    cat("    GO:BP: Error\n")
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
      results$GO_MF <- ego_mf@result
      cat("    GO:MF:", nrow(ego_mf@result), "terms\n")
    }
  }, error = function(e) {
    cat("    GO:MF: Error\n")
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
      results$GO_CC <- ego_cc@result
      cat("    GO:CC:", nrow(ego_cc@result), "terms\n")
    }
  }, error = function(e) {
    cat("    GO:CC: Error\n")
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
      results$KEGG <- ekegg_readable@result
      cat("    KEGG:", nrow(ekegg@result), "pathways\n")
    }
  }, error = function(e) {
    cat("    KEGG: Error\n")
  })
  
  # 如果没有富集结果，跳过
  if (length(results) == 0) {
    cat("  WARNING: No significant enrichment found\n\n")
    next
  }
  
  # 合并结果并准备绘图数据
  combined_results <- do.call(rbind, lapply(names(results), function(source) {
    df <- results[[source]]
    df$Source <- source
    
    df <- df[, c("Source", "ID", "Description", "GeneRatio", "BgRatio", 
                 "pvalue", "p.adjust", "qvalue", "geneID", "Count")]
    
    colnames(df) <- c("Source", "Term_ID", "Term_Name", "GeneRatio", "BgRatio",
                     "P_value", "FDR", "Q_value", "Genes", "Count")
    
    return(df)
  }))
  
  # 提取Source分类
  data <- combined_results %>%
    mutate(
      ONTOLOGY = case_when(
        grepl('GO_BP', Source) ~ 'BP',
        grepl('GO_MF', Source) ~ 'MF',
        grepl('GO_CC', Source) ~ 'CC',
        grepl('KEGG', Source) ~ 'KEGG',
        TRUE ~ 'Other'
      )
    ) %>%
    filter(ONTOLOGY != 'Other')
  
  # 筛选每个分类的top2（共8个），按FDR升序
  use_pathway <- data %>%
    group_by(ONTOLOGY) %>%
    arrange(FDR) %>%
    slice_head(n = 2) %>%
    ungroup() %>%
    mutate(
      ONTOLOGY = factor(ONTOLOGY, levels = rev(c('BP', 'CC', 'MF', 'KEGG'))),
      Count = as.numeric(Count)
    ) %>%
    arrange(ONTOLOGY, FDR) %>%
    mutate(Term_Name = factor(Term_Name, levels = Term_Name)) %>%
    tibble::rowid_to_column('index')
  
  # 构造左侧标记数据
  width <- 0.5
  xaxis_max <- ceiling(max(-log10(use_pathway$FDR), na.rm = TRUE)) + 1
  
  rect.data <- group_by(use_pathway, ONTOLOGY) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    mutate(
      xmin = -3 * width,
      xmax = -2 * width,
      ymax = cumsum(n),
      ymin = lag(ymax, default = 0) + 0.6,
      ymax = ymax + 0.4
    )
  
  # 绘制富集通路图
  p <- use_pathway %>%
    ggplot(aes(-log10(FDR), y = index, fill = ONTOLOGY)) +
    geom_round_col(
      aes(y = Term_Name), width = 0.6, alpha = 0.6
    ) +
    geom_text(
      aes(x = 0.05, label = Term_Name),
      hjust = 0, size = 5
    ) +
    geom_text(
      aes(x = 0.1, label = Genes, colour = ONTOLOGY),
      hjust = 0, vjust = 4.5, size = 3.5, fontface = 'italic',
      show.legend = FALSE
    ) +
    # 基因数量
    geom_point(
      aes(x = -width, size = Count),
      shape = 21, colour = "black"
    ) +
    geom_text(
      aes(x = -width, label = Count)
    ) +
    scale_size_continuous(name = 'Count', range = c(5, 12), breaks = function(x) {
      breaks <- pretty(x, n = 4)
      breaks <- floor(breaks)
      breaks <- unique(breaks)
      return(breaks[1:min(4, length(breaks))])
    }) +
    # 分类标签
    geom_round_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
          fill = ONTOLOGY),
      data = rect.data,
      radius = unit(2, 'mm'),
      inherit.aes = FALSE
    ) +
    geom_text(
      aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
      data = rect.data,
      inherit.aes = FALSE
    ) +
    geom_segment(
      aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
      linewidth = 1.5,
      inherit.aes = FALSE
    ) +
    labs(y = NULL, x = "-log10(FDR)") +
    scale_fill_manual(name = 'Category', values = pal, breaks = c('BP', 'CC', 'MF', 'KEGG')) +
    scale_colour_manual(values = pal_dark, breaks = c('BP', 'CC', 'MF', 'KEGG')) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(c(0, 0.05))
    ) +
    theme_prism() +
    theme(
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text()
    )
  
  # 保存图形
  output_name <- file.path(output_dir, paste0(comp_name, '_enrichment.png'))
  ggsave(output_name, plot = p, width = 10, height = 9, dpi = 300)
  cat("  [OK] Saved:", basename(output_name), "\n")
  pdf(file.path(output_dir, paste0(comp_name, '_enrichment.pdf')), width = 10, height = 9)
  print(p)
  dev.off()
  cat("  [OK] PDF saved:\n\n")
  
  # 存储每个比较组的富集结果供后续跨组合对比图使用
  all_enrichment_results[[comp_name]] <- combined_results
}

# ========== 跨组合对比点图：四类各top4通路，六组比较 ==========
cat("Generating cross-comparison dot plot...\n")

# 收集所有组合的所有通路 FDR 数据
comp_colors <- c(
  'EOD_vs_CN'   = '#E41A1C',
  'LOD_vs_CN'   = '#377EB8',
  'EOD_vs_LOD'  = '#4DAF4A',
  'EOAD_vs_CN'  = '#FF7F00',
  'LOAD_vs_CN'  = '#984EA3',
  'EOAD_vs_LOAD' = '#A65628'
)

# 合并所有组合的富集结果，标记来源
all_combined <- lapply(names(all_enrichment_results), function(comp) {
  df <- all_enrichment_results[[comp]]
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$Comparison <- comp
  df
})
all_combined <- do.call(rbind, Filter(Negate(is.null), all_combined))

if (!is.null(all_combined) && nrow(all_combined) > 0) {
  # all_combined 已含 Term_Name, FDR, Source, ONTOLOGY 等列
  all_combined <- all_combined %>%
    mutate(ONTOLOGY = case_when(
      grepl('GO_BP', Source) ~ 'BP',
      grepl('GO_MF', Source) ~ 'MF',
      grepl('GO_CC', Source) ~ 'CC',
      grepl('KEGG',  Source) ~ 'KEGG',
      TRUE ~ 'Other'
    )) %>%
    filter(ONTOLOGY != 'Other')
  
  # 每类取在任意组合中 FDR 最小、且至少一个组合显著（FDR<0.05）的 top4 通路
  top_terms <- all_combined %>%
    group_by(ONTOLOGY, Term_Name) %>%
    summarise(
      min_FDR  = min(FDR, na.rm = TRUE),
      any_sig  = any(FDR < 0.05, na.rm = TRUE),
      .groups  = 'drop'
    ) %>%
    filter(any_sig) %>%
    group_by(ONTOLOGY) %>%
    arrange(min_FDR) %>%
    slice_head(n = 4) %>%
    ungroup() %>%
    arrange(ONTOLOGY, min_FDR) %>%
    mutate(row_idx = row_number())
  
  # 筛选这些通路的所有组合数据
  dot_data <- all_combined %>%
    filter(Term_Name %in% top_terms$Term_Name) %>%
    left_join(top_terms %>% select(Term_Name, row_idx), by = 'Term_Name') %>%
    mutate(
      ONTOLOGY   = factor(ONTOLOGY, levels = c('BP', 'MF', 'CC', 'KEGG')),
      Term_Name  = factor(Term_Name, levels = rev(top_terms$Term_Name)),
      Comparison = factor(Comparison, levels = names(comp_colors)),
      neg_logFDR = -log10(FDR + 1e-300),
      is_sig     = FDR < 0.05
    )
  
  # 构造交替背景条纹数据（每行通路交替灰白）
  n_terms <- nlevels(dot_data$Term_Name)
  stripe_data <- data.frame(
    ymin       = seq(0.5, n_terms - 0.5, by = 1)[seq(1, n_terms, 2)],
    ymax       = seq(1.5, n_terms + 0.5, by = 1)[seq(1, n_terms, 2)],
    Term_Name  = factor(levels(dot_data$Term_Name)[seq(1, n_terms, 2)],
                        levels = levels(dot_data$Term_Name))
  )
  
  # 绘制跨组合对比点图
  p_dot <- ggplot(dot_data, aes(x = neg_logFDR, y = Term_Name)) +
    # 交替背景条纹
    geom_rect(
      data = stripe_data,
      aes(xmin = -Inf, xmax = Inf, ymin = as.numeric(Term_Name) - 0.5,
          ymax = as.numeric(Term_Name) + 0.5),
      fill = '#f5f5f5', colour = NA, inherit.aes = FALSE
    ) +
    # 显著阈值线
    geom_vline(xintercept = -log10(0.05), linetype = 'dashed',
               colour = 'grey60', linewidth = 0.5) +
    # 不显著点（空心小点）
    geom_point(
      data = subset(dot_data, !is_sig),
      aes(colour = Comparison),
      size = 2, alpha = 0.35, shape = 1,
      position = position_dodge(width = 0.65)
    ) +
    # 显著点（实心大点）
    geom_point(
      data = subset(dot_data, is_sig),
      aes(colour = Comparison),
      size = 3.5, alpha = 0.9, shape = 19,
      position = position_dodge(width = 0.65)
    ) +
    facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
    scale_colour_manual(name = 'Comparison', values = comp_colors) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5),
                       expand = expansion(c(0.02, 0.05))) +
    labs(x = '-log10(FDR)', y = NULL,
         caption = 'Filled: FDR < 0.05; Open: FDR ≥ 0.05') +
    theme_prism(base_size = 11) +
    theme(
      strip.text.y    = element_text(angle = 0, hjust = 0),
      legend.position = 'right',
      panel.spacing   = unit(0.4, 'lines'),
      plot.caption    = element_text(size = 8, colour = 'grey50'),
      axis.text.y     = element_text(size = 9)
    )
  
  dot_w <- 13
  dot_h <- max(8, nrow(top_terms) * 0.55 + 3)
  dot_output <- file.path(output_dir, 'cross_comparison_dot.png')
  ggsave(dot_output, plot = p_dot, width = dot_w, height = dot_h, dpi = 300)
  cat("[OK] Cross-comparison dot plot saved:", basename(dot_output), "\n")
  pdf(file.path(output_dir, 'cross_comparison_dot.pdf'), width = dot_w, height = dot_h)
  print(p_dot)
  dev.off()
  cat("[OK] Cross-comparison dot PDF saved\n\n")
} else {
  cat("WARNING: No enrichment data for cross-comparison plot\n\n")
}

cat("\n")
cat("================================================================================\n")
cat("PART 2: GSEA Analysis and Visualization\n")
cat("================================================================================\n\n")

# 输出全量 GO/KEGG 富集结果 CSV
if (length(all_enrichment_results) > 0) {
  gokegg_all <- do.call(rbind, lapply(names(all_enrichment_results), function(comp) {
    df <- all_enrichment_results[[comp]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$Comparison <- comp
    df
  }))
  gokegg_all <- Filter(Negate(is.null), list(gokegg_all))
  if (length(gokegg_all) > 0) {
    gokegg_all <- gokegg_all[[1]]
    gokegg_csv <- file.path(output_dir, "GOKEGG_all_results.csv")
    write.csv(gokegg_all, gokegg_csv, row.names = FALSE)
    cat("[OK] Full GO/KEGG results saved:", basename(gokegg_csv), "\n\n")
  }
}

# 初始化 GSEA 全量结果累积列表
all_gsea_results <- list()

# 处理每个组合的GSEA
for (comp_name in meta_files) {
  tryCatch({
    cat("Processing", comp_name, "GSEA...\n")
    
    # 加载meta分析结果
    meta_file <- file.path("F:/1a-EOD-CSF-protein/meta", paste0(comp_name, ".csv"))
    
    if (!file.exists(meta_file)) {
      cat("  WARNING: Meta file not found:", meta_file, "\n\n")
      next
    }
    
    meta_df <- read.csv(meta_file, stringsAsFactors = FALSE)
    
    # 使用所有蛋白（基因名已经是SYMBOL格式）
    valid_proteins <- !is.na(meta_df$Weighted_Effect) & (meta_df$Weighted_Effect != 0)
    
    if (sum(valid_proteins) < 10) {
      cat("  WARNING: Too few valid proteins:", sum(valid_proteins), "\n\n")
      next
    }
    
    meta_df_valid <- meta_df[valid_proteins, ]
    
    # 基因名已经是SYMBOL，直接转换为Entrez ID
    gene_mapping <- suppressMessages(bitr(
      meta_df_valid$Protein,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    ))
    
    if (is.null(gene_mapping) || nrow(gene_mapping) < 10) {
      cat("  WARNING: Gene mapping failed or too few genes\n\n")
      next
    }
    
    # 匹配回原始数据
    matched_idx <- match(gene_mapping$SYMBOL, meta_df_valid$Protein)
    
    # 创建排序的基因列表（使用Weighted_Effect）
    gene_list <- meta_df_valid$Weighted_Effect[matched_idx]
    names(gene_list) <- gene_mapping$ENTREZID
    
    # 按效应值降序排序
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    cat("  Gene list:", length(gene_list), "genes\n")
    
    # 运行GSEA for GO BP
    gsea_go_bp <- gseGO(
      geneList = gene_list,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    # 运行GSEA for GO MF
    gsea_go_mf <- gseGO(
      geneList = gene_list,
      OrgDb = org.Hs.eg.db,
      ont = "MF",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    # 运行GSEA for GO CC
    gsea_go_cc <- gseGO(
      geneList = gene_list,
      OrgDb = org.Hs.eg.db,
      ont = "CC",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    # 运行GSEA for KEGG
    gsea_kegg <- gseKEGG(
      geneList = gene_list,
      organism = 'hsa',
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE
    )

    # 合并 BP/MF/CC/KEGG 结果
    combined_results <- data.frame()
    
    if (!is.null(gsea_go_bp) && nrow(gsea_go_bp@result) > 0) {
      combined_results <- rbind(combined_results, gsea_go_bp@result)
    }
    if (!is.null(gsea_go_mf) && nrow(gsea_go_mf@result) > 0) {
      combined_results <- rbind(combined_results, gsea_go_mf@result)
    }
    if (!is.null(gsea_go_cc) && nrow(gsea_go_cc@result) > 0) {
      combined_results <- rbind(combined_results, gsea_go_cc@result)
    }
    if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
      combined_results <- rbind(combined_results, gsea_kegg@result)
    }
    
    if (nrow(combined_results) == 0) {
      cat("  WARNING: No GSEA pathway results\n\n")
      next
    }
    
    # 筛选显著通路
    significant_data <- combined_results %>%
    filter(p.adjust < 0.05) %>%
    arrange(p.adjust)
  
    if (nrow(significant_data) == 0) {
      cat("  WARNING: No significant GSEA pathways (FDR < 0.05)\n\n")
    next
  }
  
    # 选择top通路（2个正向，2个负向）
    if (nrow(significant_data) > 4) {
      pos_pathways <- significant_data %>% filter(NES > 0) %>% slice_head(n = 2)
      neg_pathways <- significant_data %>% filter(NES < 0) %>% slice_head(n = 2)
      significant_data <- rbind(pos_pathways, neg_pathways) %>% arrange(p.adjust)
    }
    
    cat("  Selected pathways:", nrow(significant_data), "\n")
    
    # 构造TERM2GENE用于重新运行GSEA
    gs2gene <- data.frame()
    for (i in 1:nrow(significant_data)) {
      genes <- unlist(strsplit(as.character(significant_data$core_enrichment[i]), "/"))
      term_name <- as.character(significant_data$Description[i])
      gs2gene <- rbind(gs2gene, data.frame(term = term_name, gene = genes, stringsAsFactors = FALSE))
    }
    
    # 重新运行GSEA以获取gseaResult对象
    gseaRes <- GSEA(gene_list, 
                    TERM2GENE = gs2gene, 
                    pvalueCutoff = 1,
                    minGSSize = 1, 
                    maxGSSize = 10000)
    
    # 获取通路ID
    geneSetID <- as.character(significant_data$Description)
    
    # 定义颜色
    n_pathways <- length(geneSetID)
    if (n_pathways == 1) {
      curve_colors <- c("#009D73")
    } else if (n_pathways == 2) {
      curve_colors <- c("#009D73", "#E59F24")
    } else if (n_pathways == 3) {
      curve_colors <- c("#009D73", "#5BB3E4", "#E59F24")
    } else {
      curve_colors <- c("#009D73", "#5BB3E4", "#E59F24", "#000000")
    }
    
    # 绘制GSEA图
    p <- gseaNb(object = gseaRes,
                geneSetID = geneSetID,
                     subPlot = 2,
                     addPval = TRUE,
                pvalX = 0.15, 
                pvalY = 0.3,
                     rmHt = TRUE,
                     termWidth = 50,
                     base_size = 12,
                     legend.position = c(0.8, 0.85),
                curveCol = curve_colors[1:n_pathways]
    )
    
    # 保存图形
    output_name <- file.path(output_dir, paste0(comp_name, "_GSEA.png"))
    ggsave(output_name, plot = p, width = 10, height = 8, dpi = 300)
    cat("  [OK] GSEA plot saved:", basename(output_name), "\n")
    pdf(file.path(output_dir, paste0(comp_name, "_GSEA.pdf")), width = 10, height = 8)
    print(p)
    dev.off()
    cat("  [OK] GSEA PDF saved\n\n")
    
    # 存储全量 GSEA 结果（不筛选，保留所有通路）
    all_gsea_results[[comp_name]] <- combined_results %>%
      mutate(Comparison = comp_name)
    
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n\n")
  })
}

# 输出全量 GSEA 结果 CSV
if (length(all_gsea_results) > 0) {
  gsea_all <- do.call(rbind, all_gsea_results)
  gsea_csv <- file.path(output_dir, "GSEA_all_results.csv")
  write.csv(gsea_all, gsea_csv, row.names = FALSE)
  cat("[OK] Full GSEA results saved:", basename(gsea_csv), "\n\n")
}

cat("\n================================================================================\n")
cat("ALL ANALYSES COMPLETE\n")
cat("================================================================================\n")
cat("Generated 12 plots (6 enrichment + 6 GSEA)\n")
cat("Output directory:", output_dir, "\n")
cat("================================================================================\n")
