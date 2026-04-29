library(ComplexHeatmap)
library(circlize)
library(grid)
library(openxlsx)

cat("Starting new consensus heatmap generation...\n")

# 设置工作目录
setwd("F:/1a-EOD-CSF-protein/1a-figure")
cat("Working directory set to:", getwd(), "\n")

# 创建输出目录
out_dir <- "F:/1a-EOD-CSF-protein/wgcna_consensus_new/figure"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 读取 GOKEGG 模块名称映射
gokegg <- read.csv("F:/1a-EOD-CSF-protein/wgcna_consensus_new/GOKEGG.csv", stringsAsFactors = FALSE)
cat("GOKEGG mapping loaded:", nrow(gokegg), "rows\n")

# 定义颜色映射
color_map <- c(
  black = "#000000", blue = "#0000FF", brown = "#A52A2A", cyan = "#00FFFF",
  darkgreen = "#006400", darkred = "#8B0000", green = "#00FF00", greenyellow = "#ADFF2F",
  grey = "#BEBEBE", grey60 = "#999999", lightcyan = "#E0FFFF", lightgreen = "#90EE90",
  lightyellow = "#FFFFE0", magenta = "#FF00FF", midnightblue = "#191970", pink = "#FFC0CB",
  purple = "#A020F0", red = "#FF0000", royalblue = "#4169E1", salmon = "#FA8072",
  tan = "#D2B48C", turquoise = "#40E0D0", yellow = "#FFFF00", gold = "#FFD700"
)

# 读取 consensus 模块名称映射（用于 ORA 列名）
consensus_assignments <- read.csv(
  "F:/1a-EOD-CSF-protein/wgcna_consensus_main/consensus_module_assignments.csv",
  stringsAsFactors = FALSE
)
# 建立颜色 -> Module_Name 映射（如 turquoise -> M1_turquoise）
consensus_name_map <- setNames(
  consensus_assignments$Module_Name,
  consensus_assignments$Module
)
consensus_name_map <- consensus_name_map[!duplicated(names(consensus_name_map))]
cat("Consensus module name map loaded:", length(consensus_name_map), "modules\n")

# 定义主模块顺序（consensus模块，颜色名）
main_module_colors <- c("pink", "yellow", "cyan", "red", "grey", "midnightblue",
                        "turquoise", "green", "greenyellow", "purple", "salmon",
                        "lightcyan", "blue", "brown", "magenta", "black", "tan")
# 将主模块颜色转换为 Mx_color 名称作为列名
main_modules_named <- sapply(main_module_colors, function(col) {
  if (col %in% names(consensus_name_map)) consensus_name_map[col] else col
})

# 定义诊断文件夹
diagnoses <- c("EODLB", "EODSD", "EOFTD", "LOAD", "LODLB", "LOFTD")

# 固定表型列顺序：只保留 CN、各诊断表型、Age、Sex
# 各诊断的 *_vs_CN 列在读取后会重命名为诊断名本身（如 EODLB_vs_CN -> EODLB）
trait_keep <- c("CN", "EODLB", "EODSD", "EOFTD", "EOAD", "EOD",
                "LOAD", "LODLB", "LOFTD", "Age", "Sex")

# 为每个诊断生成热图
for (diagnosis in diagnoses) {
  cat("\n=== Processing", diagnosis, "===\n")
  
  # 读取数据
  corr_file <- paste0("F:/1a-EOD-CSF-protein/wgcna_consensus_new/", diagnosis, "/", diagnosis, "_correlations.csv")
  pval_file <- paste0("F:/1a-EOD-CSF-protein/wgcna_consensus_new/", diagnosis, "/", diagnosis, "_pvalues.csv")
  ora_file  <- paste0("F:/1a-EOD-CSF-protein/wgcna_consensus_new/", diagnosis, "/", diagnosis, "_ORA_vs_consensus.csv")
  
  corr_data <- read.csv(corr_file, row.names = 1)
  pval_data <- read.csv(pval_file, row.names = 1)
  ora_data  <- read.csv(ora_file)
  
  cat("Data loaded for", diagnosis, "\n")
  
  # 获取该诊断的模块（从correlations.csv的行名）
  actual_modules <- rownames(corr_data)
  
  # 创建模块信息数据框
  diag_modules <- data.frame(
    Module = actual_modules,
    color  = sapply(strsplit(actual_modules, "_"), function(x) x[2]),
    stringsAsFactors = FALSE
  )
  
  # 从 GOKEGG.csv 中按 Diagnosis + Module 颜色匹配 Name
  diag_modules$Name <- sapply(diag_modules$color, function(col) {
    idx <- which(gokegg$Diagnosis == diagnosis & gokegg$Module == col)
    if (length(idx) > 0) {
      return(as.character(gokegg$Name[idx[1]]))
    } else {
      return(col)  # 找不到时回退到颜色名
    }
  })
  
  # 将 *_vs_CN 列重命名为诊断名（如 EODLB_vs_CN -> EODLB）
  vs_cn_col <- paste0(diagnosis, "_vs_CN")
  if (vs_cn_col %in% colnames(corr_data)) {
    colnames(corr_data)[colnames(corr_data) == vs_cn_col] <- diagnosis
    colnames(pval_data)[colnames(pval_data) == vs_cn_col] <- diagnosis
  }
  
  # 提取相关性和p值矩阵（动态筛选：按 trait_keep 顺序保留实际存在的列）
  available_traits <- trait_keep[trait_keep %in% colnames(corr_data)]
  mat_trait <- as.matrix(corr_data[, available_traits, drop = FALSE])
  pval_trait <- as.matrix(pval_data[, available_traits, drop = FALSE])
  
  # 生成显著性标记
  stars_trait <- ifelse(pval_trait < 0.001, "***", 
                       ifelse(pval_trait < 0.01, "**", 
                             ifelse(pval_trait < 0.05, "*", "")))
  
  # 准备ORA矩阵（主模块的p值，列名使用 Mx_color 格式）
  n_modules <- nrow(diag_modules)
  mat_ora <- matrix(NA, nrow = n_modules, ncol = length(main_module_colors))
  rownames(mat_ora) <- diag_modules$Module
  colnames(mat_ora) <- main_modules_named
  
  for (i in 1:n_modules) {
    new_module_color <- diag_modules$color[i]
    for (j in 1:length(main_module_colors)) {
      consensus_module <- main_module_colors[j]
      idx <- which(ora_data$New_Module == new_module_color & 
                   ora_data$Consensus_Module == consensus_module)
      if (length(idx) > 0) {
        mat_ora[i, j] <- ora_data$P_value[idx]
      }
    }
  }
  
  # 转换为-log10(p)用于可视化
  mat_ora_log <- -log10(mat_ora + 1e-300)
  mat_ora_log[is.infinite(mat_ora_log)] <- 300
  
  # 行名使用 GOKEGG Name
  new_row_names <- diag_modules$Name
  rownames(mat_trait)   <- new_row_names
  rownames(pval_trait)  <- new_row_names
  rownames(stars_trait) <- new_row_names
  rownames(mat_ora_log) <- new_row_names
  
  # 提取颜色用于左侧注释
  row_colors <- sapply(diag_modules$color, function(x) {
    if (x %in% names(color_map)) color_map[x] else "#808080"
  })
  names(row_colors) <- new_row_names
  
  # 左侧颜色条
  left_anno <- rowAnnotation(
    Module = new_row_names,
    col = list(Module = row_colors),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    simple_anno_size = unit(5, "mm")
  )
  
  # 定义颜色函数
  col_fun_corr <- colorRamp2(c(-1, 0, 1), c("#6baed6", "white", "#fc8d59"))
  col_fun_ora  <- colorRamp2(c(0, 1.3, 2, 5), c("white", "#FFF7AC", "#ECB477", "#D2691E"))
  
  # 设置尺寸
  MY_ROW_HEIGHT      <- unit(0.8, "cm")
  MY_COL_WIDTH_TRAIT <- unit(1.2, "cm")
  MY_COL_WIDTH_ORA   <- unit(1.0, "cm")
  
  # 顶部标注
  top_anno_trait <- HeatmapAnnotation(
    group = anno_block(
      gp = gpar(fill = "white", col = "grey40", lwd = 3),
      labels = "Traits", 
      labels_gp = gpar(fontsize = 12)
    )
  )
  
  top_anno_ora <- HeatmapAnnotation(
    group = anno_block(
      gp = gpar(fill = "white", col = "grey40", lwd = 3),
      labels = "ORA vs Main Modules", 
      labels_gp = gpar(fontsize = 12)
    )
  )
  
  # 构建热图1：Traits
  ht1 <- Heatmap(mat_trait, 
                 name = "Correlation", 
                 col = col_fun_corr,
                 width  = ncol(mat_trait) * MY_COL_WIDTH_TRAIT, 
                 height = nrow(mat_trait) * MY_ROW_HEIGHT,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 show_row_names = TRUE,
                 row_names_side = "left",
                 row_names_gp = gpar(fontsize = 10),
                 top_annotation = top_anno_trait, 
                 left_annotation = left_anno,
                 column_title = NULL,
                 na_col = "white",
                 rect_gp = gpar(col = "grey40", lwd = 2),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   if (is.na(mat_trait[i, j])) {
                     grid.rect(x, y, width, height, gp = gpar(col = "grey40", lwd = 2, fill = "white"))
                     grid.text("NA", x, y, gp = gpar(fontsize = 8, col = "grey50"))
                   } else {
                     grid.text(stars_trait[i, j], x, y, gp = gpar(fontsize = 9))
                   }
                 },
                 column_names_rot = 45,
                 column_names_gp = gpar(fontsize = 10),
                 heatmap_legend_param = list(
                   title = "Correlation",
                   border = "black",
                   legend_height = unit(4, "cm")
                 )
  )
  
  # 构建热图2：ORA
  ht2 <- Heatmap(mat_ora_log, 
                 name = "ORA -log10(p)", 
                 col = col_fun_ora,
                 width  = ncol(mat_ora_log) * MY_COL_WIDTH_ORA, 
                 height = nrow(mat_ora_log) * MY_ROW_HEIGHT,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 show_row_names = FALSE,
                 top_annotation = top_anno_ora,
                 column_title = NULL,
                 na_col = "white",
                 rect_gp = gpar(col = "grey40", lwd = 2),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   if (!is.na(mat_ora[i, j])) {
                     if (mat_ora[i, j] < 0.001) {
                       grid.text("***", x, y, gp = gpar(fontsize = 8))
                     } else if (mat_ora[i, j] < 0.01) {
                       grid.text("**", x, y, gp = gpar(fontsize = 8))
                     } else if (mat_ora[i, j] < 0.05) {
                       grid.text("*", x, y, gp = gpar(fontsize = 8))
                     }
                   }
                 },
                 column_names_rot = 45,
                 column_names_gp = gpar(fontsize = 9),
                 heatmap_legend_param = list(
                   title = "ORA\n-log10(p-value)",
                   border = "black",
                   legend_height = unit(4, "cm")
                 )
  )
  
  # 导出数据到Excel
  wb <- createWorkbook()
  
  addWorksheet(wb, "Correlation")
  writeData(wb, "Correlation", cbind(Module = rownames(mat_trait), mat_trait), rowNames = FALSE)
  
  addWorksheet(wb, "Pvalue")
  writeData(wb, "Pvalue", cbind(Module = rownames(pval_trait), pval_trait), rowNames = FALSE)
  
  addWorksheet(wb, "ORA")
  writeData(wb, "ORA", cbind(Module = rownames(mat_ora), mat_ora), rowNames = FALSE)
  
  excel_file <- file.path(out_dir, paste0(diagnosis, "_heatmap_data.xlsx"))
  saveWorkbook(wb, excel_file, overwrite = TRUE)
  cat("Data exported to", excel_file, "\n")
  
  # 绘制并保存热图
  png_file <- file.path(out_dir, paste0(diagnosis, "_heatmap.png"))
  png(png_file, width = 36, height = max(12, nrow(mat_trait) * 0.8), units = "in", res = 300)
  draw(ht1 + ht2, 
       ht_gap = unit(3, "mm"), 
       merge_legend = FALSE, 
       heatmap_legend_side = "bottom")
  dev.off()
  cat("Heatmap saved to", png_file, "\n")
  
  pdf_file <- file.path(out_dir, paste0(diagnosis, "_heatmap.pdf"))
  pdf(pdf_file, width = 36, height = max(12, nrow(mat_trait) * 0.8))
  draw(ht1 + ht2, 
       ht_gap = unit(3, "mm"), 
       merge_legend = FALSE, 
       heatmap_legend_side = "bottom")
  dev.off()
  cat("PDF saved to", pdf_file, "\n")
}

cat("\nAll heatmaps generated successfully!\n")
