library(ComplexHeatmap)
library(circlize)
library(grid)
library(openxlsx)

cat("Starting heatmap generation...\n")

# 设置工作目录
setwd("F:/1a-EOD-CSF-protein/1a-figure")
cat("Working directory set to:", getwd(), "\n")

# 读取数据
cat("Reading data...\n")
module_enrich <- read.csv("F:/1a-EOD-CSF-protein/module_enrich.txt")
cat("module_enrich loaded\n")
corr_main <- read.csv("F:/1a-EOD-CSF-protein/wgcna_consensus_main/consensus_module_trait_correlations.csv", row.names = 1)
cat("corr_main loaded\n")
pval_main <- read.csv("F:/1a-EOD-CSF-protein/wgcna_consensus_main/consensus_module_trait_pvalues.csv", row.names = 1)
cat("pval_main loaded\n")
celltype_enrich <- read.csv("F:/1a-EOD-CSF-protein/wgcna_consensus_main/module_celltype_fisher_enrichment.csv")
cat("celltype_enrich loaded\n")
preservation <- read.csv("F:/1a-EOD-CSF-protein/wgcna_consensus_test/module_preservation_by_study.csv")
cat("preservation loaded\n")

# 定义颜色映射
color_map <- c(
  black = "#000000", blue = "#0000FF", brown = "#A52A2A", cyan = "#00FFFF",
  darkgreen = "#006400", darkred = "#8B0000", green = "#00FF00", greenyellow = "#ADFF2F",
  grey = "#BEBEBE", grey60 = "#999999", lightcyan = "#E0FFFF", lightgreen = "#90EE90",
  lightyellow = "#FFFFE0", magenta = "#FF00FF", midnightblue = "#191970", pink = "#FFC0CB",
  purple = "#A020F0", red = "#FF0000", royalblue = "#4169E1", salmon = "#FA8072",
  tan = "#D2B48C", turquoise = "#40E0D0", yellow = "#FFFF00"
)

# 按module_enrich.txt的顺序排列
ordered_modules <- module_enrich$Module
module_names_new <- module_enrich$New.Name
new_row_names <- module_names_new

# 提取颜色
row_colors <- sapply(ordered_modules, function(x) {
  color_name <- strsplit(x, "_")[[1]][2]
  color_map[color_name]
})
names(row_colors) <- new_row_names

# 按顺序排列数据
corr_main <- corr_main[ordered_modules, ]
pval_main <- pval_main[ordered_modules, ]
rownames(corr_main) <- new_row_names
rownames(pval_main) <- new_row_names

# ========== Main图数据准备 ==========
# Biomarker顺序: MoCA, MMSE, AB42, AB40, AB42/pTau, pTau, pTau181, pTau217, pTau231, tTau, NEFL, YKL40
biomarker_cols <- c("MoCA", "MMSE", "AB42", "AB40", "AB42.pTau", "pTau", "pTau181", "pTau217", "pTau231", "tTau", "NEFL", "YKL40")
available_biomarkers <- biomarker_cols[biomarker_cols %in% colnames(corr_main)]
cat("Available biomarkers:", paste(available_biomarkers, collapse=", "), "\n")

# 定义study颜色
study_colors <- c(
  "study_1" = "#D62728",   # 红色
  "study_11" = "#1F77B4",  # 蓝色
  "study_4" = "#2CA02C",   # 绿色
  "study_6" = "#FF7F0E",   # 橙色
  "study_7" = "#9467BD",   # 紫色
  "study_9" = "#8C564B"    # 棕色
)

# ========== 分析原始数据中各变量在各study的存在情况 ==========
cat("Analyzing variable presence in each study...\n")
raw_data <- read.csv('F:/1a-EOD-CSF-protein/combine/combined_expression_matrices.csv')

# 定义变量映射关系（热图列名 -> 原始数据列名）
# 注意：CSV读取时会将空格和特殊字符转换为点号
var_mapping <- list(
  "CN" = "Diagnosis_Derived",
  "EOD_vs_CN" = "Diagnosis_Derived",
  "EOAD_vs_CN" = "Diagnosis_Derived",
  "EODSD_vs_CN" = "Diagnosis_Derived",
  "EODLB_vs_CN" = "Diagnosis_Derived",
  "EOFTD_vs_CN" = "Diagnosis_Derived",
  "EOD_vs_Other" = "Diagnosis_Derived",
  "EOAD_vs_Other" = "Diagnosis_Derived",
  "EODSD_vs_Other" = "Diagnosis_Derived",
  "EODLB_vs_Other" = "Diagnosis_Derived",
  "EOFTD_vs_Other" = "Diagnosis_Derived",
  "EOD" = "Diagnosis_Derived",
  "EOAD" = "Diagnosis_Derived",
  "EODSD" = "Diagnosis_Derived",
  "EODLB" = "Diagnosis_Derived",
  "EOFTD" = "Diagnosis_Derived",
  "Age" = "Age",
  "Sex" = "Sex",
  "MoCA" = "Cognitive.Score",      # CSV读取时 "Cognitive Score" -> "Cognitive.Score"
  "MMSE" = "Cognitive.Score",
  "AB42" = "AB42",
  "AB40" = "AB40",
  "AB42.pTau" = "AB42.pTau",       # CSV读取时 "AB42/pTau" -> "AB42.pTau"
  "pTau" = "pTau",
  "pTau181" = "pTau181",
  "pTau217" = "pTau217",
  "pTau231" = "pTau231",
  "tTau" = "tTau",
  "NEFL" = "NEFL",
  "YKL40" = "YKL40"
)

# 检查每个变量在每个study中是否有数据
check_var_in_study <- function(var_name, study_name) {
  raw_var_name <- var_mapping[[var_name]]
  if(is.null(raw_var_name)) return(FALSE)
  
  study_data <- raw_data[raw_data$Study == study_name, ]
  if(nrow(study_data) == 0) return(FALSE)
  
  if(!(raw_var_name %in% colnames(study_data))) return(FALSE)
  
  # 对于诊断相关的变量，需要检查是否有对应的诊断类型
  if(raw_var_name == "Diagnosis_Derived") {
    # 提取诊断类型（如 EOFTD_vs_CN -> EOFTD, EOD_vs_CN -> EOD/EOAD）
    if(var_name == "CN") {
      # CN在所有有Diagnosis_Derived数据的study中都存在
      return(sum(!is.na(study_data[[raw_var_name]])) > 0)
    } else if(grepl("EOFTD", var_name)) {
      # 检查是否有EOFTD诊断
      return(sum(study_data[[raw_var_name]] == "EOFTD", na.rm = TRUE) > 0)
    } else if(grepl("EODSD", var_name)) {
      # 检查是否有EODSD诊断
      return(sum(study_data[[raw_var_name]] == "EODSD", na.rm = TRUE) > 0)
    } else if(grepl("EODLB", var_name)) {
      # 检查是否有EODLB诊断
      return(sum(study_data[[raw_var_name]] == "EODLB", na.rm = TRUE) > 0)
    } else if(grepl("EOAD", var_name)) {
      # 检查是否有EOAD诊断
      return(sum(study_data[[raw_var_name]] == "EOAD", na.rm = TRUE) > 0)
    } else if(grepl("EOD", var_name)) {
      # EOD包括EOAD, EODSD, EODLB, EOFTD等所有早发性痴呆
      eod_diagnoses <- c("EOAD", "EODSD", "EODLB", "EOFTD", "EOOD", "EOND")
      return(sum(study_data[[raw_var_name]] %in% eod_diagnoses, na.rm = TRUE) > 0)
    } else {
      # 其他诊断相关变量
      return(sum(!is.na(study_data[[raw_var_name]])) > 0)
    }
  }
  
  # 对于非诊断变量，检查是否有非NA值
  has_data <- sum(!is.na(study_data[[raw_var_name]])) > 0
  return(has_data)
}

# 为训练集创建变量存在矩阵
train_studies <- c("study_1", "study_11")
test_studies <- c("study_4", "study_6", "study_7", "study_9")

cat("Variable presence analysis complete\n")

mat1_main <- as.matrix(corr_main[, c("CN", "EOD_vs_CN", "EOAD_vs_CN", "EOD_vs_Other", "EOAD_vs_Other")])
mat2_main <- as.matrix(corr_main[, c("Age", "Sex")])
mat3_main <- as.matrix(corr_main[, available_biomarkers])

# 创建训练集study来源矩阵（使用与主热图相同的列数）
cat("Creating study source matrices for main...\n")
study_mat1_main <- matrix("", nrow = 1, ncol = ncol(mat1_main))
colnames(study_mat1_main) <- colnames(mat1_main)

study_mat2_main <- matrix("", nrow = 1, ncol = ncol(mat2_main))
colnames(study_mat2_main) <- colnames(mat2_main)

study_mat3_main <- matrix("", nrow = 1, ncol = ncol(mat3_main))
colnames(study_mat3_main) <- colnames(mat3_main)
cat("Study source matrices created\n")

# 准备细胞类型数据（-log10(P)）
cell_types <- c("neurons", "astrocytes", "oligodendrocytes", "microglia", "endothelial", "OPCs")
mat4_main <- matrix(NA, nrow = nrow(corr_main), ncol = length(cell_types))
rownames(mat4_main) <- rownames(corr_main)
colnames(mat4_main) <- cell_types

for (i in 1:nrow(mat4_main)) {
  module_name <- ordered_modules[i]
  for (j in 1:length(cell_types)) {
    cell_type <- cell_types[j]
    idx <- which(celltype_enrich$Module == gsub("M[0-9]+_", "", module_name) & 
                 celltype_enrich$Cell_Type == cell_type)
    if (length(idx) > 0) {
      mat4_main[i, j] <- -log10(celltype_enrich$P_value[idx])
    }
  }
}

pval1_main <- as.matrix(pval_main[, c("CN", "EOD_vs_CN", "EOAD_vs_CN", "EOD_vs_Other", "EOAD_vs_Other")])
pval2_main <- as.matrix(pval_main[, c("Age", "Sex")])
pval3_main <- as.matrix(pval_main[, available_biomarkers])

stars1_main <- ifelse(pval1_main < 0.001, "***", ifelse(pval1_main < 0.01, "**", ifelse(pval1_main < 0.05, "*", "")))
stars2_main <- ifelse(pval2_main < 0.001, "***", ifelse(pval2_main < 0.01, "**", ifelse(pval2_main < 0.05, "*", "")))
stars3_main <- ifelse(pval3_main < 0.001, "***", ifelse(pval3_main < 0.01, "**", ifelse(pval3_main < 0.05, "*", "")))

# 修改列名
colnames(mat1_main) <- c("CN", "EOD", "EOAD", "EOD", "EOAD")

# 绘图设置
col_fun1 <- colorRamp2(c(-1, 0, 1), c("#6baed6", "white", "#fc8d59"))
col_fun2 <- colorRamp2(c(-1, 0, 1), c("#9ebcda", "white", "#e0ecf4"))

# Biomarker 新的颜色映射（绿-白-黄橙）
color_palette_neg <- colorRampPalette(c("#AAD09D", "#ECF4DD"))(length(seq(-2.1, -0.101, by = 0.2)))
color_palette_zero <- c("white")
color_palette_pos <- colorRampPalette(c("#FFF7AC", "#ECB477"))(length(seq(0.101, 2.1, by = 0.2)))
breaks_bio <- c(seq(-2.1, -0.101, by = 0.2), 0, seq(0.101, 2.1, by = 0.2))
colors_bio <- c(color_palette_neg, color_palette_zero, color_palette_pos)
col_fun3 <- colorRamp2(breaks_bio, colors_bio)

# Cell 颜色映射（白到淡紫色，最大值3）
col_fun4 <- colorRamp2(c(0, 1.5, 3), c("white", "#D8BFD8", "#9370DB"))

MY_ROW_HEIGHT = unit(0.8, "cm")
MY_COL_WIDTH  = unit(1.5, "cm")
STUDY_ROW_HEIGHT = unit(0.3, "cm")  # study来源行高度

# 左侧颜色条
left_anno <- rowAnnotation(
  Module = new_row_names,
  col = list(Module = row_colors),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  simple_anno_size = unit(5, "mm")
)

# 顶部标注
top_anno1_main <- HeatmapAnnotation(
  group = anno_block(
    gp = gpar(fill = "white", col = "grey40", lwd = 3),
    labels = c("CN", "versus CN", "versus Other"), 
    labels_gp = gpar(fontsize = 12)
  ),
  annotation_name_side = "left"
)

top_anno2 <- HeatmapAnnotation(
  group = anno_block(
    gp = gpar(fill = "white", col = "grey40", lwd = 3),
    labels = "Demography", 
    labels_gp = gpar(fontsize = 12)
  )
)

top_anno3 <- HeatmapAnnotation(
  group = anno_block(
    gp = gpar(fill = "white", col = "grey40", lwd = 3),
    labels = "Biomarker", 
    labels_gp = gpar(fontsize = 12)
  )
)

top_anno4 <- HeatmapAnnotation(
  group = anno_block(
    gp = gpar(fill = "white", col = "grey40", lwd = 3),
    labels = "Cell", 
    labels_gp = gpar(fontsize = 12)
  )
)

# 构建Main热图
# 创建自定义annotation函数来绘制分割的study格子（训练集：2个study）
anno_study_split_2 <- function(col_names) {
  AnnotationFunction(
    fun = function(index, k, n) {
      n_col <- length(index)
      for (i in 1:n_col) {
        var_name <- col_names[index[i]]
        
        # 每个格子分成两半
        x_left <- (i - 1) / n_col
        x_right <- i / n_col
        x_mid <- (x_left + x_right) / 2
        width_half <- (x_right - x_left) / 2
        
        # 检查study_1是否有该变量
        has_study1 <- check_var_in_study(var_name, "study_1")
        # 检查study_11是否有该变量
        has_study11 <- check_var_in_study(var_name, "study_11")
        
        # 左半部分 - study_1（只有有数据时才填充颜色）
        if(has_study1) {
          grid.rect(x = x_mid - width_half/2, y = 0.5,
                   width = width_half, height = 1,
                   gp = gpar(fill = study_colors["study_1"], col = "grey40", lwd = 1),
                   default.units = "npc")
        } else {
          grid.rect(x = x_mid - width_half/2, y = 0.5,
                   width = width_half, height = 1,
                   gp = gpar(fill = "white", col = "grey40", lwd = 1),
                   default.units = "npc")
        }
        
        # 右半部分 - study_11（只有有数据时才填充颜色）
        if(has_study11) {
          grid.rect(x = x_mid + width_half/2, y = 0.5,
                   width = width_half, height = 1,
                   gp = gpar(fill = study_colors["study_11"], col = "grey40", lwd = 1),
                   default.units = "npc")
        } else {
          grid.rect(x = x_mid + width_half/2, y = 0.5,
                   width = width_half, height = 1,
                   gp = gpar(fill = "white", col = "grey40", lwd = 1),
                   default.units = "npc")
        }
      }
    },
    var_import = list(study_colors = study_colors, check_var_in_study = check_var_in_study, col_names = col_names),
    n = length(col_names),
    height = STUDY_ROW_HEIGHT
  )
}

bottom_anno1_main <- HeatmapAnnotation(
  study = anno_study_split_2(colnames(mat1_main)),
  annotation_name_side = "left",
  show_legend = FALSE
)

ht1_main <- Heatmap(mat1_main, 
               name = "Module\nCorrelation", col = col_fun1,
               width = ncol(mat1_main) * MY_COL_WIDTH, 
               height = nrow(mat1_main) * MY_ROW_HEIGHT,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_row_names = TRUE,
               row_names_side = "left",
               row_names_gp = gpar(fontsize = 10),
               top_annotation = top_anno1_main, 
               bottom_annotation = bottom_anno1_main,
               left_annotation = left_anno,
               column_split = factor(c("CN", "vs_CN", "vs_CN", "vs_Other", "vs_Other"), 
                                    levels = c("CN", "vs_CN", "vs_Other")),
               column_gap = unit(3, "mm"),
               column_title = NULL,
               na_col = "white",
               rect_gp = gpar(col = "grey40", lwd = 3),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (is.na(mat1_main[i, j])) {
                   grid.rect(x, y, width, height, gp = gpar(col = "grey40", lwd = 3, fill = "white"))
                   grid.text("NA", x, y, gp = gpar(fontsize = 10, col = "grey50"))
                 } else {
                 grid.text(stars1_main[i, j], x, y, gp = gpar(fontsize = 10))
                 }
               },
               column_names_centered = TRUE, column_names_rot = 0,
               column_names_gp = gpar(fontsize = 10),
               heatmap_legend_param = list(
                 title = "Module\nCorrelation",
                 border = "black",
                 legend_height = unit(4, "cm")
               )
)

bottom_anno2_main <- HeatmapAnnotation(
  study = anno_study_split_2(colnames(mat2_main)),
  show_legend = FALSE
)

ht2_main <- Heatmap(mat2_main, 
               name = "Module\nCorrelation", col = col_fun1,
               width = ncol(mat2_main) * MY_COL_WIDTH, 
               height = nrow(mat2_main) * MY_ROW_HEIGHT,
               cluster_rows = FALSE,
               cluster_columns = FALSE, show_row_names = FALSE,
               top_annotation = top_anno2,
               bottom_annotation = bottom_anno2_main,
               na_col = "white",
               rect_gp = gpar(col = "grey40", lwd = 3),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (is.na(mat2_main[i, j])) {
                   grid.rect(x, y, width, height, gp = gpar(col = "grey40", lwd = 3, fill = "white"))
                   grid.text("NA", x, y, gp = gpar(fontsize = 10, col = "grey50"))
                 } else {
                 grid.text(stars2_main[i, j], x, y, gp = gpar(fontsize = 10))
                 }
               },
               column_names_rot = 0,
               column_names_gp = gpar(fontsize = 10),
               show_heatmap_legend = FALSE
)

bottom_anno3_main <- HeatmapAnnotation(
  study = anno_study_split_2(colnames(mat3_main)),
  show_legend = FALSE
)

ht3_main <- Heatmap(mat3_main, 
               name = "Biomarker", col = col_fun3,
               width = ncol(mat3_main) * MY_COL_WIDTH, 
               height = nrow(mat3_main) * MY_ROW_HEIGHT,
               cluster_rows = FALSE,
               cluster_columns = FALSE, show_row_names = FALSE,
               top_annotation = top_anno3,
               bottom_annotation = bottom_anno3_main,
               na_col = "white",
               rect_gp = gpar(col = "grey40", lwd = 3),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (is.na(mat3_main[i, j])) {
                   grid.rect(x, y, width, height, gp = gpar(col = "grey40", lwd = 3, fill = "white"))
                   grid.text("NA", x, y, gp = gpar(fontsize = 10, col = "grey50"))
                 } else {
                 grid.text(stars3_main[i, j], x, y, gp = gpar(fontsize = 10))
                 }
               },
               column_names_rot = 0,
               column_names_gp = gpar(fontsize = 10),
               heatmap_legend_param = list(
                 title = "Biomarker\nCorrelation",
                 border = "black",
                 legend_height = unit(4, "cm"),
                 at = c(-1, -0.5, 0, 0.5, 1),
                 labels = c("-1", "-0.5", "0", "0.5", "1")
               )
)

ht4_main <- Heatmap(mat4_main, 
               name = "Cell Type\n-log10(p)", col = col_fun4,
               width = ncol(mat4_main) * MY_COL_WIDTH, 
               height = nrow(mat4_main) * MY_ROW_HEIGHT,
               cluster_rows = FALSE,
               cluster_columns = FALSE, show_row_names = FALSE,
               top_annotation = top_anno4,
               na_col = "white",
               rect_gp = gpar(col = "grey40", lwd = 3),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (!is.na(mat4_main[i, j]) && mat4_main[i, j] >= -log10(0.05)) {
                   grid.text("*", x, y, gp = gpar(fontsize = 10))
                 }
               },
               column_names_rot = 0,
               column_names_gp = gpar(fontsize = 10),
               heatmap_legend_param = list(
                 title = "Cell Type\n-log10(p-value)",
                 border = "black",
                 legend_height = unit(4, "cm")
               )
)

# 导出Main数据到Excel
wb_main <- createWorkbook()

# 训练集相关度
addWorksheet(wb_main, "Main_Correlation")
main_corr_export <- cbind(
  Module = rownames(mat1_main),
  mat1_main,
  mat2_main,
  mat3_main,
  mat4_main
)
writeData(wb_main, "Main_Correlation", main_corr_export, rowNames = FALSE)

# 训练集P值
addWorksheet(wb_main, "Main_Pvalue")
main_pval_export <- cbind(
  Module = rownames(pval1_main),
  pval1_main,
  pval2_main,
  pval3_main,
  matrix(NA, nrow = nrow(mat4_main), ncol = ncol(mat4_main), 
         dimnames = list(rownames(mat4_main), paste0(colnames(mat4_main), "_pval")))
)
writeData(wb_main, "Main_Pvalue", main_pval_export, rowNames = FALSE)

saveWorkbook(wb_main, "consensus_main_heatmap_data.xlsx", overwrite = TRUE)
print("Main data exported to Excel")

png("consensus_main_heatmap.png", width = 32, height = 22, units = "in", res = 300)
ht_list_main <- ht1_main + ht2_main + ht3_main + ht4_main
draw(ht_list_main, ht_gap = unit(c(3, 3, 3), "mm"), merge_legend = FALSE, heatmap_legend_side = "bottom")

# 添加study legend
lgd <- Legend(labels = c("study_1", "study_11"), 
             title = "Study Source",
             legend_gp = gpar(fill = c(study_colors["study_1"], study_colors["study_11"])),
             border = "black")
draw(lgd, x = unit(0.85, "npc"), y = unit(0.1, "npc"), just = c("left", "bottom"))

dev.off()
print("Main heatmap saved")

pdf("consensus_main_heatmap.pdf", width = 32, height = 22)
ht_list_main <- ht1_main + ht2_main + ht3_main + ht4_main
draw(ht_list_main, ht_gap = unit(c(3, 3, 3), "mm"), merge_legend = FALSE, heatmap_legend_side = "bottom")
draw(lgd, x = unit(0.85, "npc"), y = unit(0.1, "npc"), just = c("left", "bottom"))
dev.off()
print("Main heatmap PDF saved")

# ========== Test图数据准备 ==========
corr_test <- read.csv("F:/1a-EOD-CSF-protein/wgcna_consensus_test/test_studies_consensus_correlations.csv", row.names = 1)
pval_test <- read.csv("F:/1a-EOD-CSF-protein/wgcna_consensus_test/test_studies_consensus_pvalues.csv", row.names = 1)

corr_test <- corr_test[ordered_modules, ]
pval_test <- pval_test[ordered_modules, ]
rownames(corr_test) <- new_row_names
rownames(pval_test) <- new_row_names

# Test biomarkers
available_biomarkers_test <- biomarker_cols[biomarker_cols %in% colnames(corr_test)]
cat("Available biomarkers (test):", paste(available_biomarkers_test, collapse=", "), "\n")

mat1_test <- as.matrix(corr_test[, c("CN", "EOD_vs_CN", "EOAD_vs_CN", "EODSD_vs_CN", 
                                      "EODLB_vs_CN", "EOFTD_vs_CN",
                                      "EOD_vs_Other", "EOAD_vs_Other", "EODSD_vs_Other",
                                      "EODLB_vs_Other", "EOFTD_vs_Other")])
mat2_test <- as.matrix(corr_test[, c("Age", "Sex")])
mat3_test <- as.matrix(corr_test[, available_biomarkers_test])

# 创建测试集study来源矩阵（使用与主热图相同的列数）
study_mat1_test <- matrix("", nrow = 1, ncol = ncol(mat1_test))
colnames(study_mat1_test) <- colnames(mat1_test)

study_mat2_test <- matrix("", nrow = 1, ncol = ncol(mat2_test))
colnames(study_mat2_test) <- colnames(mat2_test)

study_mat3_test <- matrix("", nrow = 1, ncol = ncol(mat3_test))
colnames(study_mat3_test) <- colnames(mat3_test)

# 准备保守度数据（Zsummary）
preservation_studies <- c("study_4", "study_6", "study_7", "study_9", "study_12")
mat4_test <- matrix(NA, nrow = nrow(corr_test), ncol = length(preservation_studies))
rownames(mat4_test) <- rownames(corr_test)
colnames(mat4_test) <- preservation_studies

for (i in 1:nrow(mat4_test)) {
  module_name <- ordered_modules[i]
  # 提取颜色名（如 M14_blue -> blue）
  color_name <- strsplit(module_name, "_")[[1]][2]
  idx <- which(preservation$Module == color_name)
  if (length(idx) > 0) {
    for (j in 1:length(preservation_studies)) {
      study <- preservation_studies[j]
      val <- preservation[idx, study]
      mat4_test[i, j] <- ifelse(is.na(val) || val < 0, 0, val)
    }
  }
}

pval1_test <- as.matrix(pval_test[, c("CN", "EOD_vs_CN", "EOAD_vs_CN", "EODSD_vs_CN", 
                                       "EODLB_vs_CN", "EOFTD_vs_CN",
                                       "EOD_vs_Other", "EOAD_vs_Other", "EODSD_vs_Other",
                                       "EODLB_vs_Other", "EOFTD_vs_Other")])
pval2_test <- as.matrix(pval_test[, c("Age", "Sex")])
pval3_test <- as.matrix(pval_test[, available_biomarkers_test])

stars1_test <- ifelse(pval1_test < 0.001, "***", ifelse(pval1_test < 0.01, "**", ifelse(pval1_test < 0.05, "*", "")))
stars2_test <- ifelse(pval2_test < 0.001, "***", ifelse(pval2_test < 0.01, "**", ifelse(pval2_test < 0.05, "*", "")))
stars3_test <- ifelse(pval3_test < 0.001, "***", ifelse(pval3_test < 0.01, "**", ifelse(pval3_test < 0.05, "*", "")))

colnames(mat1_test) <- c("CN", "EOD", "EOAD", "EODSD", "EODLB", "EOFTD",
                         "EOD", "EOAD", "EODSD", "EODLB", "EOFTD")

top_anno1_test <- HeatmapAnnotation(
  group = anno_block(
    gp = gpar(fill = "white", col = "grey40", lwd = 3),
    labels = c("CN", "versus CN", "versus Other"), 
    labels_gp = gpar(fontsize = 12)
  ),
  annotation_name_side = "left"
)

top_anno4_test <- HeatmapAnnotation(
  group = anno_block(
    gp = gpar(fill = "white", col = "grey40", lwd = 3),
    labels = "Preservation", 
    labels_gp = gpar(fontsize = 12)
  )
)

# Preservation颜色映射（0-10）
col_fun4_test <- colorRamp2(c(0, 1, 2, 10), c("#9EBCC8", "#9DCFC2", "#B5DBBA", "#D3E4A6"))

# 创建自定义annotation函数来绘制分割成4份的study格子（测试集：4个study）
anno_study_split_4 <- function(col_names) {
  AnnotationFunction(
    fun = function(index, k, n) {
      n_col <- length(index)
      for (i in 1:n_col) {
        var_name <- col_names[index[i]]
        
        # 每个格子分成四份
        x_left <- (i - 1) / n_col
        x_right <- i / n_col
        x_mid <- (x_left + x_right) / 2
        width_quarter <- (x_right - x_left) / 4
        
        # 检查各个study是否有该变量
        has_study4 <- check_var_in_study(var_name, "study_4")
        has_study6 <- check_var_in_study(var_name, "study_6")
        has_study7 <- check_var_in_study(var_name, "study_7")
        has_study9 <- check_var_in_study(var_name, "study_9")
        
        # 第1份 - study_4（只有有数据时才填充颜色）
        if(has_study4) {
          grid.rect(x = x_mid - 3*width_quarter/2, y = 0.5,
                   width = width_quarter, height = 1,
                   gp = gpar(fill = study_colors["study_4"], col = "grey40", lwd = 1),
                   default.units = "npc")
        } else {
          grid.rect(x = x_mid - 3*width_quarter/2, y = 0.5,
                   width = width_quarter, height = 1,
                   gp = gpar(fill = "white", col = "grey40", lwd = 1),
                   default.units = "npc")
        }
        
        # 第2份 - study_6（只有有数据时才填充颜色）
        if(has_study6) {
          grid.rect(x = x_mid - width_quarter/2, y = 0.5,
                   width = width_quarter, height = 1,
                   gp = gpar(fill = study_colors["study_6"], col = "grey40", lwd = 1),
                   default.units = "npc")
        } else {
          grid.rect(x = x_mid - width_quarter/2, y = 0.5,
                   width = width_quarter, height = 1,
                   gp = gpar(fill = "white", col = "grey40", lwd = 1),
                   default.units = "npc")
        }
        
        # 第3份 - study_7（只有有数据时才填充颜色）
        if(has_study7) {
          grid.rect(x = x_mid + width_quarter/2, y = 0.5,
                   width = width_quarter, height = 1,
                   gp = gpar(fill = study_colors["study_7"], col = "grey40", lwd = 1),
                   default.units = "npc")
        } else {
          grid.rect(x = x_mid + width_quarter/2, y = 0.5,
                   width = width_quarter, height = 1,
                   gp = gpar(fill = "white", col = "grey40", lwd = 1),
                   default.units = "npc")
        }
        
        # 第4份 - study_9（只有有数据时才填充颜色）
        if(has_study9) {
          grid.rect(x = x_mid + 3*width_quarter/2, y = 0.5,
                   width = width_quarter, height = 1,
                   gp = gpar(fill = study_colors["study_9"], col = "grey40", lwd = 1),
                   default.units = "npc")
        } else {
          grid.rect(x = x_mid + 3*width_quarter/2, y = 0.5,
                   width = width_quarter, height = 1,
                   gp = gpar(fill = "white", col = "grey40", lwd = 1),
                   default.units = "npc")
        }
      }
    },
    var_import = list(study_colors = study_colors, check_var_in_study = check_var_in_study, col_names = col_names),
    n = length(col_names),
    height = STUDY_ROW_HEIGHT
  )
}

bottom_anno1_test <- HeatmapAnnotation(
  study = anno_study_split_4(colnames(mat1_test)),
  show_legend = FALSE
)

ht1_test <- Heatmap(mat1_test, 
               name = "Module\nCorrelation", col = col_fun1,
               width = ncol(mat1_test) * MY_COL_WIDTH, 
               height = nrow(mat1_test) * MY_ROW_HEIGHT,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_row_names = TRUE,
               row_names_side = "left",
               row_names_gp = gpar(fontsize = 10),
               top_annotation = top_anno1_test, 
               bottom_annotation = bottom_anno1_test,
               left_annotation = left_anno,
               column_split = factor(c("CN", rep("vs_CN", 5), rep("vs_Other", 5)), 
                                    levels = c("CN", "vs_CN", "vs_Other")),
               column_gap = unit(3, "mm"),
               column_title = NULL,
               na_col = "white",
               rect_gp = gpar(col = "grey40", lwd = 3),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (is.na(mat1_test[i, j])) {
                   grid.rect(x, y, width, height, gp = gpar(col = "grey40", lwd = 3, fill = "white"))
                   grid.text("NA", x, y, gp = gpar(fontsize = 10, col = "grey50"))
                 } else {
                 grid.text(stars1_test[i, j], x, y, gp = gpar(fontsize = 10))
                 }
               },
               column_names_centered = TRUE, column_names_rot = 0,
               column_names_gp = gpar(fontsize = 10),
               heatmap_legend_param = list(
                 title = "Module\nCorrelation",
                 border = "black",
                 legend_height = unit(4, "cm")
               )
)

bottom_anno2_test <- HeatmapAnnotation(
  study = anno_study_split_4(colnames(mat2_test)),
  show_legend = FALSE
)

ht2_test <- Heatmap(mat2_test, 
               name = "Module\nCorrelation", col = col_fun1,
               width = ncol(mat2_test) * MY_COL_WIDTH, 
               height = nrow(mat2_test) * MY_ROW_HEIGHT,
               cluster_rows = FALSE,
               cluster_columns = FALSE, show_row_names = FALSE,
               top_annotation = top_anno2,
               bottom_annotation = bottom_anno2_test,
               na_col = "white",
               rect_gp = gpar(col = "grey40", lwd = 3),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (is.na(mat2_test[i, j])) {
                   grid.rect(x, y, width, height, gp = gpar(col = "grey40", lwd = 3, fill = "white"))
                   grid.text("NA", x, y, gp = gpar(fontsize = 10, col = "grey50"))
                 } else {
                 grid.text(stars2_test[i, j], x, y, gp = gpar(fontsize = 10))
                 }
               },
               column_names_rot = 0,
               column_names_gp = gpar(fontsize = 10),
               show_heatmap_legend = FALSE
)

bottom_anno3_test <- HeatmapAnnotation(
  study = anno_study_split_4(colnames(mat3_test)),
  show_legend = FALSE
)

ht3_test <- Heatmap(mat3_test, 
               name = "Biomarker", col = col_fun3,
               width = ncol(mat3_test) * MY_COL_WIDTH, 
               height = nrow(mat3_test) * MY_ROW_HEIGHT,
               cluster_rows = FALSE,
               cluster_columns = FALSE, show_row_names = FALSE,
               top_annotation = top_anno3,
               bottom_annotation = bottom_anno3_test,
               na_col = "white",
               rect_gp = gpar(col = "grey40", lwd = 3),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (is.na(mat3_test[i, j])) {
                   grid.rect(x, y, width, height, gp = gpar(col = "grey40", lwd = 3, fill = "white"))
                   grid.text("NA", x, y, gp = gpar(fontsize = 10, col = "grey50"))
                 } else {
                 grid.text(stars3_test[i, j], x, y, gp = gpar(fontsize = 10))
                 }
               },
               column_names_rot = 0,
               column_names_gp = gpar(fontsize = 10),
               heatmap_legend_param = list(
                 title = "Biomarker\nCorrelation",
                 border = "black",
                 legend_height = unit(4, "cm"),
                 at = c(-1, -0.5, 0, 0.5, 1),
                 labels = c("-1", "-0.5", "0", "0.5", "1")
               )
)

ht4_test <- Heatmap(mat4_test, 
               name = "Preservation\nZsummary", col = col_fun4_test,
               width = ncol(mat4_test) * MY_COL_WIDTH, 
               height = nrow(mat4_test) * MY_ROW_HEIGHT,
               cluster_rows = FALSE,
               cluster_columns = FALSE, show_row_names = FALSE,
               top_annotation = top_anno4_test,
               na_col = "white",
               rect_gp = gpar(col = "grey40", lwd = 3),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (!is.na(mat4_test[i, j])) {
                   if (mat4_test[i, j] >= 10) {
                     grid.text("##", x, y, gp = gpar(fontsize = 10))
                   } else if (mat4_test[i, j] >= 1.96) {
                     grid.text("#", x, y, gp = gpar(fontsize = 10))
                   }
                 }
               },
               column_names_rot = 0,
               column_names_gp = gpar(fontsize = 10),
               heatmap_legend_param = list(
                 title = "Preservation\nZsummary",
                 border = "black",
                 legend_height = unit(4, "cm"),
                 at = c(0, 2, 10),
                 labels = c("0", "2", "10")
               )
)

# 导出Test数据到Excel
wb_test <- createWorkbook()

# 测试集相关度
addWorksheet(wb_test, "Test_Correlation")
test_corr_export <- cbind(
  Module = rownames(mat1_test),
  mat1_test,
  mat2_test,
  mat3_test,
  mat4_test
)
writeData(wb_test, "Test_Correlation", test_corr_export, rowNames = FALSE)

# 测试集P值
addWorksheet(wb_test, "Test_Pvalue")
test_pval_export <- cbind(
  Module = rownames(pval1_test),
  pval1_test,
  pval2_test,
  pval3_test,
  matrix(NA, nrow = nrow(mat4_test), ncol = ncol(mat4_test), 
         dimnames = list(rownames(mat4_test), paste0(colnames(mat4_test), "_pval")))
)
writeData(wb_test, "Test_Pvalue", test_pval_export, rowNames = FALSE)

saveWorkbook(wb_test, "consensus_test_heatmap_data.xlsx", overwrite = TRUE)
print("Test data exported to Excel")

png("consensus_test_heatmap.png", width = 42, height = 22, units = "in", res = 300)
ht_list_test <- ht1_test + ht2_test + ht3_test + ht4_test
draw(ht_list_test, ht_gap = unit(c(3, 3, 3), "mm"), merge_legend = FALSE, heatmap_legend_side = "bottom")

# 添加study legend
lgd_test <- Legend(labels = c("study_4", "study_6", "study_7", "study_9"), 
                  title = "Study Source",
                  legend_gp = gpar(fill = c(study_colors["study_4"], study_colors["study_6"], 
                                           study_colors["study_7"], study_colors["study_9"])),
                  border = "black")
draw(lgd_test, x = unit(0.85, "npc"), y = unit(0.1, "npc"), just = c("left", "bottom"))

dev.off()
print("Test heatmap saved")

pdf("consensus_test_heatmap.pdf", width = 42, height = 22)
ht_list_test <- ht1_test + ht2_test + ht3_test + ht4_test
draw(ht_list_test, ht_gap = unit(c(3, 3, 3), "mm"), merge_legend = FALSE, heatmap_legend_side = "bottom")
draw(lgd_test, x = unit(0.85, "npc"), y = unit(0.1, "npc"), just = c("left", "bottom"))
dev.off()
print("Test heatmap PDF saved")
