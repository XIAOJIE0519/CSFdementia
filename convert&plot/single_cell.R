library(ComplexHeatmap)
library(circlize)
library(grid)

# 输出目录
out_dir <- "F:/1a-痴呆亚型整合/EOD/figures/module_heatmaps"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 读取数据
# ============================================================
mat_raw <- read.csv(
  "F:/1a-痴呆亚型整合/EOD/figures/module_heatmaps/module_expression_by_celltype_disease.csv",
  row.names = 1, check.names = FALSE
)
cat("Data loaded:", nrow(mat_raw), "rows x", ncol(mat_raw), "cols\n")

# ============================================================
# 模块顺序（与热图.R一致，来自 module_enrich.txt）
# ============================================================
module_order <- c(
  "M5_pink", "M6_yellow", "M14_cyan", "M9_red", "M17_grey",
  "M12_midnightblue", "M1_turquoise", "M11_green", "M13_greenyellow",
  "M15_purple", "M10_salmon", "M16_lightcyan", "M2_blue",
  "M4_brown", "M8_magenta", "M7_black", "M3_tan"
)

# 简短行标签（仅 M1、M4 等数字编号）
module_short_names <- gsub("_.*", "", module_order)  # M5, M6, ...

# 模块颜色映射
color_map <- c(
  black = "#000000", blue = "#0000FF", brown = "#A52A2A", cyan = "#00FFFF",
  darkgreen = "#006400", darkred = "#8B0000", green = "#00FF00", greenyellow = "#ADFF2F",
  grey = "#BEBEBE", grey60 = "#999999", lightcyan = "#E0FFFF", lightgreen = "#90EE90",
  lightyellow = "#FFFFE0", magenta = "#FF00FF", midnightblue = "#191970", pink = "#FFC0CB",
  purple = "#A020F0", red = "#FF0000", royalblue = "#4169E1", salmon = "#FA8072",
  tan = "#D2B48C", turquoise = "#40E0D0", yellow = "#FFFF00"
)

module_colors <- sapply(module_order, function(m) {
  col_name <- sub("M[0-9]+_", "", m)
  if (col_name %in% names(color_map)) color_map[[col_name]] else "#808080"
})
names(module_colors) <- module_short_names

# ============================================================
# 列顺序：按细胞类型分组，每组内 CN/EOAD/LOAD/EOFTD/LOFTD
# ============================================================
cell_types_order <- c(
  "Astrocytes", "EndothelialCells", "ExcitatoryNeurons", "InhibitoryNeurons",
  "Macrophages", "Microglia", "Oligodendrocytes", "OPCs", "Pericytes"
)
disease_order <- c("CN", "EOAD", "LOAD", "EOFTD", "LOFTD")

col_order <- c()
for (ct in cell_types_order) {
  for (dg in disease_order) {
    cname <- paste0(ct, "_", dg)
    if (cname %in% rownames(mat_raw)) col_order <- c(col_order, cname)
  }
}

# ============================================================
# 构建矩阵：行=模块，列=CellType_DiseaseGroup
# ============================================================
avail_modules <- module_order[module_order %in% colnames(mat_raw)]
mat <- t(as.matrix(mat_raw[col_order, avail_modules]))
rownames(mat) <- module_short_names[module_order %in% avail_modules]
cat("Matrix dim:", nrow(mat), "rows x", ncol(mat), "cols\n")

# ============================================================
# 解析列的细胞类型和疾病分组（归一化前需要）
# ============================================================
col_parts     <- strsplit(col_order, "_")
col_celltypes <- sapply(col_parts, function(x) paste(x[-length(x)], collapse = "_"))
col_diseases  <- sapply(col_parts, function(x) x[length(x)])

# ============================================================
# 归一化：每个细胞类型 × 模块内，对5个疾病表型做 z-score
# ============================================================
mat_norm <- mat
for (ct in cell_types_order) {
  ct_idx <- which(col_celltypes == ct)
  if (length(ct_idx) < 2) next
  for (ri in seq_len(nrow(mat))) {
    vals <- mat[ri, ct_idx]
    m <- mean(vals, na.rm = TRUE)
    s <- sd(vals, na.rm = TRUE)
    if (!is.na(s) && s > 0) {
      mat_norm[ri, ct_idx] <- (vals - m) / s
    } else {
      mat_norm[ri, ct_idx] <- 0
    }
  }
}
mat <- mat_norm
cat("Normalization done (z-score within each cell type x module)\n")

# ============================================================
# 配色（来自 eod_06_visualization.py）
# ============================================================
cell_type_colors <- c(
  ExcitatoryNeurons  = "#1f77b4",
  InhibitoryNeurons  = "#ff7f0e",
  Astrocytes         = "#2ca02c",
  Oligodendrocytes   = "#d62728",
  OPCs               = "#9467bd",
  Microglia          = "#8c564b",
  Macrophages        = "#e377c2",
  EndothelialCells   = "#7f7f7f",
  Pericytes          = "#bcbd22"
)

disease_colors <- c(
  CN    = "#2ca02c",
  EOAD  = "#d62728",
  LOAD  = "#ff7f0e",
  EOFTD = "#9467bd",
  LOFTD = "#8c564b"
)

# ============================================================
# 顶部注释：疾病分类（第一排）+ 细胞类型（第二排）
# ============================================================
top_anno <- HeatmapAnnotation(
  Disease  = col_diseases,
  CellType = col_celltypes,
  col = list(
    Disease  = disease_colors,
    CellType = cell_type_colors
  ),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 9),
  simple_anno_size     = unit(4, "mm"),
  show_legend          = TRUE,
  annotation_legend_param = list(
    Disease  = list(title = "Disease",   title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8)),
    CellType = list(title = "Cell Type", title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8))
  )
)

# ============================================================
# 左侧注释：模块颜色块
# ============================================================
row_module_colors <- module_colors[rownames(mat)]
left_anno <- rowAnnotation(
  Module = rownames(mat),
  col    = list(Module = row_module_colors),
  show_legend          = FALSE,
  show_annotation_name = FALSE,
  simple_anno_size     = unit(5, "mm")
)

# ============================================================
# 颜色函数（归一化后 z-score，范围约 -2 到 2）
# ============================================================
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#313695", "white", "#A50026")
)

# ============================================================
# 行分组（每个模块一间隔，row_split 按行名分组）
# ============================================================
row_split <- factor(rownames(mat), levels = rownames(mat))

# ============================================================
# 列分组（细胞类型间隔）
# ============================================================
col_split <- factor(col_celltypes, levels = cell_types_order)

# ============================================================
# 构建热图
# ============================================================
CELL_SIZE <- unit(0.5, "cm")  # 正方形格子边长

ht <- Heatmap(
  mat,
  name = "Z-score",
  col  = col_fun,

  # 行
  cluster_rows         = FALSE,
  row_split            = row_split,
  row_gap              = unit(1.5, "mm"),
  row_title            = NULL,
  show_row_names       = TRUE,
  row_names_side       = "left",
  row_names_gp         = gpar(fontsize = 9),
  left_annotation      = left_anno,
  # 列
  cluster_columns      = FALSE,
  column_split         = col_split,
  column_gap           = unit(2, "mm"),
  show_column_names    = FALSE,
  column_title_rot     = 45,
  column_title_gp      = gpar(fontsize = 9),
  top_annotation       = top_anno,

  # 正方形格子
  width  = CELL_SIZE * ncol(mat),
  height = CELL_SIZE * nrow(mat),

  # 单元格
  border               = FALSE,
  rect_gp              = gpar(col = NA),

  # 图例
  heatmap_legend_param = list(
    title             = "Z-score",
    legend_height     = unit(4, "cm"),
    title_gp          = gpar(fontsize = 9),
    labels_gp         = gpar(fontsize = 8),
    border            = "black"
  )
)

# ============================================================
# 输出 PNG
# ============================================================
# 输出尺寸：格子0.5cm + 行间隔 + 注释 + 图例留白
fig_w <- ncol(mat) * 0.5 / 2.54 + 8
fig_h <- nrow(mat) * 0.5 / 2.54 + nrow(mat) * 0.06 + 5

# 绘图函数（含每5列分组边框，覆盖所有行分组）
draw_with_border <- function() {
  draw(ht, merge_legend = FALSE, heatmap_legend_side = "right")
  # row_split 将每行分为单独一组，需对每个 row_slice x column_slice 组合画框
  n_row_slices <- nrow(mat)
  n_col_slices <- length(cell_types_order)
  for (ri in seq_len(n_row_slices)) {
    for (ci in seq_len(n_col_slices)) {
      decorate_heatmap_body("Z-score", row_slice = ri, column_slice = ci, {
        grid.rect(gp = gpar(col = "grey40", lwd = 1.5, fill = NA))
      })
    }
  }
}

png_file <- file.path(out_dir, "module_expression_heatmap.png")
png(png_file, width = fig_w, height = fig_h, units = "in", res = 300)
draw_with_border()
dev.off()
cat("PNG saved to", png_file, "\n")

pdf_file <- file.path(out_dir, "module_expression_heatmap.pdf")
pdf(pdf_file, width = fig_w, height = fig_h)
draw_with_border()
dev.off()
cat("PDF saved to", pdf_file, "\n")

cat("Done!\n")
