# ============================================================================
# EOD Single-Cell Integration Pipeline - Mantel Test Analysis
# ============================================================================
# 
# 模块间相关性热图 + 疾病-模块Mantel检验连线图
# 
# 输入：
#   - EOD/figures/module_heatmaps/module_expression_by_disease.csv
# 
# 输出：
#   - EOD/figures/mantel_correlation_heatmap.png (相关性热图+Mantel检验)
#   - EOD/figures/mantel_results.csv (Mantel检验结果)
# 
# Author: Bioinformatics Analysis
# Date: 2026-03-07
# ============================================================================

# 清空环境变量
rm(list = ls())

# 加载必需的包
suppressPackageStartupMessages({
  library(linkET)
  library(ggplot2)
  library(dplyr)
  library(vegan)
  library(RColorBrewer)
})

cat("\n")
cat("================================================================================\n")
cat("EOD Pipeline - Mantel Test Analysis\n")
cat("================================================================================\n\n")

# 设置随机种子
set.seed(123)

# ============================================================================
# 1. 加载模块表达数据（细胞类型×疾病组 × 模块）
# ============================================================================

cat("Loading module expression data...\n")

input_file <- "EOD/figures/module_heatmaps/module_expression_by_celltype_disease.csv"
if (!file.exists(input_file)) {
  stop("Error: Module expression file not found: ", input_file)
}

# 读取数据（第一列是CellType_DiseaseGroup，作为行名）
module_expr <- read.csv(input_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

cat(sprintf("  Loaded: %d samples x %d modules\n", nrow(module_expr), ncol(module_expr)))

# 提取细胞类型和疾病分组信息
sample_info <- data.frame(
  sample_name = rownames(module_expr),
  stringsAsFactors = FALSE
)

# 分割样本名称（格式：CellType_DiseaseGroup）
sample_info$cell_type <- sapply(strsplit(sample_info$sample_name, "_"), function(x) {
  paste(x[1:(length(x)-1)], collapse="_")
})
sample_info$disease_group <- sapply(strsplit(sample_info$sample_name, "_"), function(x) {
  x[length(x)]
})

cat(sprintf("  Cell types: %s\n", paste(unique(sample_info$cell_type), collapse=", ")))
cat(sprintf("  Disease groups: %s\n", paste(unique(sample_info$disease_group), collapse=", ")))

# ============================================================================
# 2. 提取模块编号和颜色
# ============================================================================

cat("\nExtracting module information...\n")

# 从列名提取模块编号和颜色（格式：M1_turquoise）
module_info <- data.frame(
  full_name = colnames(module_expr),
  stringsAsFactors = FALSE
)

module_info$module_num <- sapply(strsplit(module_info$full_name, "_"), function(x) x[1])
module_info$color_name <- sapply(strsplit(module_info$full_name, "_"), function(x) x[2])

cat(sprintf("  Extracted %d modules\n", nrow(module_info)))

# 定义WGCNA标准颜色映射
wgcna_colors <- c(
  "turquoise" = "#40E0D0",
  "blue" = "#0000FF",
  "tan" = "#D2B48C",
  "brown" = "#A52A2A",
  "pink" = "#FFC0CB",
  "yellow" = "#FFFF00",
  "black" = "#000000",
  "magenta" = "#FF00FF",
  "red" = "#FF0000",
  "salmon" = "#FA8072",
  "green" = "#00FF00",
  "midnightblue" = "#191970",
  "greenyellow" = "#ADFF2F",
  "cyan" = "#00FFFF",
  "purple" = "#A020F0",
  "lightcyan" = "#E0FFFF",
  "grey" = "#BEBEBE"
)

# 映射颜色
module_info$hex_color <- wgcna_colors[module_info$color_name]

cat("  Module colors mapped\n")

# ============================================================================
# 3. 计算模块间相关性（用于热图）
# ============================================================================

cat("\nCalculating module-module correlations...\n")

# 计算Pearson相关性（行=样本，列=模块）
module_cor <- cor(module_expr, method = "pearson")

cat(sprintf("  Correlation matrix: %d x %d\n", nrow(module_cor), ncol(module_cor)))

# ============================================================================
# 4. 准备Mantel检验数据
# ============================================================================

cat("\nPreparing data for Mantel test...\n")

# 创建环境因子矩阵（疾病分组的虚拟变量）
disease_groups <- unique(sample_info$disease_group)
env_matrix <- model.matrix(~ sample_info$disease_group - 1)
colnames(env_matrix) <- disease_groups
rownames(env_matrix) <- rownames(module_expr)

cat(sprintf("  Module data: %d samples x %d modules\n", nrow(module_expr), ncol(module_expr)))
cat(sprintf("  Environment data: %d samples x %d disease groups\n", nrow(env_matrix), ncol(env_matrix)))

# ============================================================================
# 5. Mantel检验：疾病-模块关联
# ============================================================================

cat("\nPerforming Mantel test (Disease-Module associations)...\n")

# 按疾病分组进行Mantel检验
# spec_select: 每个疾病组作为一个"物种组"
disease_select_list <- lapply(disease_groups, function(dg) {
  which(colnames(env_matrix) == dg)
})
names(disease_select_list) <- disease_groups

# 执行Mantel检验
# env_matrix作为"物种"（疾病分组），module_expr作为"环境"（模块表达）
mantel_result <- mantel_test(
  env_matrix,           # 疾病分组作为"物种"
  module_expr,          # 模块表达作为"环境"
  mantel_fun = 'mantel',
  spec_dist = 'euclidean',    # 疾病分组使用欧氏距离
  env_dist = 'euclidean',     # 模块表达使用欧氏距离
  spec_select = disease_select_list
) %>%
  mutate(
    rd = cut(r, breaks = c(-Inf, -0.4, -0.2, 0.2, 0.4, Inf),
             labels = c("<= -0.4", "-0.4 - -0.2", "-0.2 - 0.2", "0.2 - 0.4", ">= 0.4")),
    linetype = cut(r, breaks = c(-Inf, 0, Inf),
                   labels = c("Negative", "Positive")),
    pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
             labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))
  )

# 调整因子水平顺序
mantel_result$rd <- factor(mantel_result$rd, 
                           levels = c("-0.2 - 0.2", "0.2 - 0.4", ">= 0.4", 
                                      "-0.4 - -0.2", "<= -0.4"))

cat("  Mantel test completed\n")
cat(sprintf("  Total tests: %d\n", nrow(mantel_result)))

# 保存结果
output_dir <- "EOD/figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

result_file <- file.path(output_dir, "mantel_results.csv")
write.csv(mantel_result, result_file, row.names = FALSE)
cat(sprintf("  Saved results to: %s\n", result_file))

# 打印显著结果
cat("\nSignificant Mantel test results (p < 0.05):\n")
sig_results <- mantel_result %>% filter(p < 0.05)
if (nrow(sig_results) > 0) {
  print(sig_results)
} else {
  cat("  No significant results found.\n")
}

# ============================================================================
# 6. 可视化：模块间相关性热图 + 疾病-模块Mantel检验连线图
# ============================================================================

cat("\nGenerating correlation heatmap with Mantel test...\n")

# 创建相关性数据框（用于qcorrplot）
module_cor_df <- as.data.frame(module_cor)

# 绘制热图
p1 <- qcorrplot(module_cor_df, type = 'upper', diag = TRUE) +
  geom_square() +
  geom_mark(sep = '\n', size = 2.5,
            sig_level = c(0.05, 0.01, 0.001),
            mark = c("*", "**", "***"),
            fontface = 1) +
  geom_couple(aes(color = pd, size = rd, linetype = linetype),
              data = mantel_result,
              curvature = nice_curvature(),
              label.size = 3,
              label.fontface = 1) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_fill_gradientn(colors = c("#2ca02c", "white", "#FF8C00"),  # 绿色-白色-橙色
                       limits = c(-1, 1)) +
  scale_size_manual(values = c(0.8, 1.3, 2)) +
  scale_color_manual(values = c("< 0.01" = "#d62728",      # 红色
                                "0.01 - 0.05" = "#ff7f0e",  # 橙色
                                ">= 0.05" = "#E0E0E0")) +   # 灰色
  guides(color = guide_legend(title = "Mantel's p", order = 1),
         size = guide_legend(title = "Mantel's r", order = 2),
         fill = guide_colorbar(title = "Pearson's r", order = 4),
         linetype = guide_legend(title = "Correlation", order = 3)) +
  labs(title = "Module Correlation Heatmap with Mantel Test (Disease-Module)") +
  theme(axis.text.x.top = element_text(size = 10, angle = 45, hjust = 0, vjust = 0, face = "bold"),
        axis.text.y = element_text(size = 10, angle = 0, face = "bold"),
        text = element_text(size = 10, face = 1),
        title = element_text(size = 12, face = "bold"))

# 保存图形
output_file <- file.path(output_dir, "mantel_correlation_heatmap.png")
ggsave(output_file, plot = p1, width = 14, height = 12, dpi = 300)
cat(sprintf("  Saved: %s\n", output_file))

cat("\n")
cat("================================================================================\n")
cat("Mantel Test Analysis Complete\n")
cat("================================================================================\n\n")
