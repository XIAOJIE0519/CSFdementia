# ============================================================================
# EOD Single-Cell Integration Pipeline - Correlation Analysis
# ============================================================================
# 
# 模块间相关性热图 + 疾病-模块相关性连线图
# 
# 输入：
#   - EOD/figures/module_heatmaps/module_expression_by_celltype_disease.csv
# 
# 输出：
#   - EOD/figures/correlation_heatmap.png (相关性热图+连线)
# 
# Author: Bioinformatics Analysis
# Date: 2026-03-08
# ============================================================================

# 清空环境变量
rm(list = ls())

# 加载必需的包
suppressPackageStartupMessages({
  library(linkET)
  library(ggplot2)
  library(dplyr)
})

cat("\n")
cat("================================================================================\n")
cat("EOD Pipeline - Correlation Analysis\n")
cat("================================================================================\n\n")

# 设置随机种子
set.seed(123)

# ============================================================================
# 1. 加载模块表达数据（按真实样本）
# ============================================================================

cat("Loading module expression data by sample...\n")

input_file <- "EOD/figures/module_heatmaps/module_expression_by_sample.csv"
if (!file.exists(input_file)) {
  stop("Error: Module expression file not found: ", input_file)
}

# 读取数据
module_expr_with_meta <- read.csv(input_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

cat(sprintf("  Loaded: %d samples x %d columns\n", nrow(module_expr_with_meta), ncol(module_expr_with_meta)))

# 分离元数据和模块表达
meta_cols <- c("disease_group", "data_source", "subtype")
sample_metadata <- module_expr_with_meta[, meta_cols]
module_expr <- module_expr_with_meta[, !colnames(module_expr_with_meta) %in% meta_cols]

# 按指定顺序排列模块
module_order <- c("M5_pink", "M6_yellow", "M14_cyan", "M9_red", "M17_grey",
                  "M12_midnightblue", "M1_turquoise", "M11_green", "M13_greenyellow",
                  "M15_purple", "M10_salmon", "M16_lightcyan", "M2_blue",
                  "M4_brown", "M8_magenta", "M7_black", "M3_tan")

# 只保留存在的模块
module_order <- module_order[module_order %in% colnames(module_expr)]
module_expr <- module_expr[, module_order]

# 保存全名用于行标签
module_full_names <- colnames(module_expr)

# 简化模块名称用于列标签（M1_pink → M1）
colnames(module_expr) <- gsub("(M\\d+)_.*", "\\1", colnames(module_expr))

cat(sprintf("  Samples: %d\n", nrow(module_expr)))
cat(sprintf("  Modules: %d\n", ncol(module_expr)))
cat(sprintf("  Module order: %s\n", paste(colnames(module_expr)[1:5], "...", collapse=", ")))
cat(sprintf("  Disease groups: %s\n", paste(unique(sample_metadata$disease_group), collapse=", ")))

# ============================================================================
# 2. 计算模块-模块相关性（用于热图）
# ============================================================================

cat("\nCalculating module-module correlations...\n")

# 使用correlate函数
module_cor_result <- correlate(module_expr, method = "pearson")

# 修改行名为全名，列名保持简化名
attr(module_cor_result, "row.names") <- module_full_names

cat(sprintf("  Correlation matrix: %d x %d\n", nrow(module_cor_result), ncol(module_cor_result)))

# ============================================================================
# 3. 计算疾病-模块相关性（用于连线）
# ============================================================================

cat("\nCalculating disease-module correlations...\n")

# 创建疾病分组矩阵（虚拟变量）
disease_groups <- unique(sample_metadata$disease_group)
env_matrix <- model.matrix(~ sample_metadata$disease_group - 1)
colnames(env_matrix) <- disease_groups
rownames(env_matrix) <- rownames(module_expr)

# 直接计算相关性和p值（不使用correlate函数）
n_diseases <- ncol(env_matrix)
n_modules <- ncol(module_expr)

disease_module_cor <- matrix(NA, nrow = n_diseases, ncol = n_modules)
rownames(disease_module_cor) <- colnames(env_matrix)
colnames(disease_module_cor) <- colnames(module_expr)

disease_module_pval <- matrix(NA, nrow = n_diseases, ncol = n_modules)
rownames(disease_module_pval) <- colnames(env_matrix)
colnames(disease_module_pval) <- colnames(module_expr)

for (i in 1:n_diseases) {
  for (j in 1:n_modules) {
    test_result <- cor.test(env_matrix[, i], module_expr[, j], method = "pearson")
    disease_module_cor[i, j] <- test_result$estimate
    disease_module_pval[i, j] <- test_result$p.value
  }
}

cat(sprintf("  Correlation matrix: %d diseases x %d modules\n", 
            nrow(disease_module_cor), ncol(disease_module_cor)))

# ============================================================================
# 4. 准备连线数据（疾病→模块）
# ============================================================================

cat("\nPreparing coupling data for visualization...\n")

cat(sprintf("  Correlation matrix: %d x %d\n", nrow(disease_module_cor), ncol(disease_module_cor)))
cat(sprintf("  P-value matrix: %d x %d\n", nrow(disease_module_pval), ncol(disease_module_pval)))

# 转换为长格式
coupling_data <- data.frame()

for (i in 1:nrow(disease_module_cor)) {
  for (j in 1:ncol(disease_module_cor)) {
    coupling_data <- rbind(coupling_data, data.frame(
      spec = rownames(disease_module_cor)[i],
      env = colnames(disease_module_cor)[j],
      r = disease_module_cor[i, j],
      p = disease_module_pval[i, j]
    ))
  }
}

# 添加分类变量（只保留相关系数分类，不用p值分类）
coupling_data <- coupling_data %>%
  mutate(
    rd = cut(abs(r), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
             labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")),
    linetype = cut(r, breaks = c(-Inf, 0, Inf),
                   labels = c("Negative", "Positive"))
  )

cat(sprintf("  Prepared %d coupling relationships\n", nrow(coupling_data)))

# 打印显著结果
sig_results <- coupling_data %>% filter(p < 0.05)
cat(sprintf("  Significant correlations (p < 0.05): %d\n", nrow(sig_results)))

# 检查分类变量的级别
cat(sprintf("  R categories: %s\n", paste(levels(coupling_data$rd), collapse=", ")))
cat(sprintf("  Linetype categories: %s\n", paste(levels(coupling_data$linetype), collapse=", ")))

# ============================================================================
# 5. 可视化：模块间相关性热图 + 疾病-模块连线图
# ============================================================================

cat("\nGenerating correlation heatmap with coupling...\n")

# 创建输出目录
output_dir <- "EOD/figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 绘制组合图
p1 <- qcorrplot(module_cor_result, type = 'upper', diag = TRUE) +
  geom_square(colour = NA) +
  geom_couple(aes(colour = r, size = rd, linetype = linetype),
              data = coupling_data,
              curvature = nice_curvature()) +
  scale_linetype_manual(values = c("Negative" = 2, "Positive" = 1)) +
  scale_fill_gradientn(colors = c("#6baed6", "white", "#fc8d59"),
                       limits = c(-1, 1),
                       name = "Pearson's r") +
  scale_colour_gradientn(colors = c("#6baed6", "white", "#fc8d59"),
                         limits = c(-1, 1),
                         guide = "none") +
  scale_size_manual(values = c("0-0.2" = 0.5, 
                               "0.2-0.4" = 1.0, 
                               "0.4-0.6" = 1.5,
                               "0.6-0.8" = 2.0,
                               "0.8-1.0" = 2.5)) +
  guides(fill = guide_colorbar(title = "Pearson's r", order = 1),
         size = guide_legend(title = "|r|", order = 2),
         linetype = guide_legend(title = "Direction", order = 3)) +
  labs(title = "Module Correlation Heatmap with Disease-Module Coupling") +
  theme(axis.text.x.top = element_text(size = 10, angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 10, angle = 0, face = "bold"),
        text = element_text(size = 10, face = 1),
        title = element_text(size = 12, face = "bold"))

# 保存图形
output_file <- file.path(output_dir, "correlation_heatmap.png")
ggsave(output_file, plot = p1, width = 14, height = 12, dpi = 300)
cat(sprintf("  Saved: %s\n", output_file))

# 保存数据
write.csv(as.matrix(module_cor_result), file.path(output_dir, "module_module_correlation.csv"))
write.csv(disease_module_cor, file.path(output_dir, "disease_module_correlation.csv"))
write.csv(coupling_data, file.path(output_dir, "disease_module_coupling_data.csv"), row.names = FALSE)

cat("\n")
cat("================================================================================\n")
cat("Correlation Analysis Complete\n")
cat("================================================================================\n\n")
