# 设置工作目录
setwd("F:/1a-EOD-CSF-protein/1a-figure")

cat("=== PCA分析开始 ===\n")

tryCatch({
  # 加载必要的包
  cat("加载R包...\n")
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
  })
  cat("R包加载完成\n")
  
  # 读取数据
  cat("读取combined_expression_matrices.csv...\n")
  data <- read.csv("F:/1a-EOD-CSF-protein/combine/combined_expression_matrices.csv", 
                   stringsAsFactors = FALSE, check.names = FALSE)
  
  cat("数据维度:", nrow(data), "行 x", ncol(data), "列\n")
  
  # 查找APOE|P02649列的位置
  apoe_col <- which(colnames(data) == "APOE|P02649")
  if (length(apoe_col) == 0) {
    stop("未找到APOE|P02649列")
  }
  
  cat("APOE|P02649列位置:", apoe_col, "\n")
  
  # 提取样本信息和蛋白质数据
  sample_info <- data[, 1:6]
  colnames(sample_info) <- c("Study", "Batch", "Sample_ID", "Age", "Sex", "Diagnosis")
  protein_data <- data[, apoe_col:ncol(data)]
  
  cat("样本信息列:", ncol(sample_info), "\n")
  cat("蛋白质数据列:", ncol(protein_data), "\n")
  cat("总样本数:", nrow(data), "\n")
  
  # 检查并移除缺失值过多的样本和蛋白质
  cat("处理缺失值...\n")
  missing_per_sample <- rowSums(is.na(protein_data)) / ncol(protein_data)
  valid_samples <- missing_per_sample < 0.5
  cat("移除缺失值过多的样本后剩余:", sum(valid_samples), "个样本\n")
  
  sample_info <- sample_info[valid_samples, ]
  protein_data <- protein_data[valid_samples, ]
  
  # 移除缺失值超过50%的蛋白质
  missing_per_protein <- colSums(is.na(protein_data)) / nrow(protein_data)
  valid_proteins <- missing_per_protein < 0.5
  cat("移除缺失值过多的蛋白质后剩余:", sum(valid_proteins), "个蛋白质\n")
  
  protein_data <- protein_data[, valid_proteins]
  
  # 对剩余的缺失值进行填充（使用列均值）
  cat("填充剩余缺失值...\n")
  for (i in 1:ncol(protein_data)) {
    if (any(is.na(protein_data[, i]))) {
      protein_data[is.na(protein_data[, i]), i] <- mean(protein_data[, i], na.rm = TRUE)
    }
  }
  
  # 转换为数值矩阵
  protein_matrix <- as.matrix(protein_data)
  cat("蛋白质矩阵维度:", nrow(protein_matrix), "x", ncol(protein_matrix), "\n")
  
  # 执行PCA
  cat("执行PCA分析...\n")
  pca_result <- prcomp(protein_matrix, center = TRUE, scale. = TRUE)
  cat("PCA分析完成\n")
  
  # 提取PC1和PC2
  pca_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Study = sample_info$Study,
    Sample_ID = sample_info$Sample_ID,
    Age = sample_info$Age,
    Sex = sample_info$Sex,
    Diagnosis = sample_info$Diagnosis
  )
  
  # 计算方差解释比例
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100
  
  cat("PC1方差解释:", round(var_explained[1], 2), "%\n")
  cat("PC2方差解释:", round(var_explained[2], 2), "%\n")
  
  # 定义研究的颜色（12个研究）
  study_colors <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
    "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
    "#999999", "#66C2A5", "#FC8D62", "#8DA0CB"
  )
  
  # 绘制PCA图
  cat("绘制PCA图...\n")
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Study)) +
    geom_point(size = 2, alpha = 0.6) +
    scale_color_manual(values = study_colors) +
    labs(
      title = "PCA Analysis of CSF Protein Expression",
      subtitle = paste0(nrow(pca_data), " samples from 12 studies"),
      x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 2), "%)")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  # 保存图形
  output_file <- "PCA_plot.png"
  ggsave(output_file, plot = p, width = 12, height = 8, dpi = 300)
  cat("PCA图已保存:", output_file, "\n")
  
  # 保存PCA结果数据
  write.csv(pca_data, "PCA_results.csv", row.names = FALSE)
  cat("PCA结果数据已保存: PCA_results.csv\n")
  
  # 打印研究统计
  cat("\n各研究样本数:\n")
  print(table(pca_data$Study))
  
  cat("\n=== PCA分析完成 ===\n")
  
}, error = function(e) {
  cat("\n!!! 错误 !!!\n")
  cat("错误信息:", conditionMessage(e), "\n")
  print(traceback())
})
