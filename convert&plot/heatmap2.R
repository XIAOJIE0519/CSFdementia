# 1. 安装并加载包
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

# 2. 准备模拟数据 (实际操作中请替换为你的矩阵)
# 假设 mat 是一个基因在不同细胞类型中的平均表达量矩阵
set.seed(123)
mat <- matrix(rnorm(500), nrow = 50)
rownames(mat) <- paste0("Gene", 1:50)
colnames(mat) <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")

# 3. 数据标准化 (Z-score)，这是热图显示差异的关键
mat_scaled <- t(scale(t(mat)))

# 4. 设置颜色映射 (紫色 -> 白色 -> 橙色)
col_fun = colorRamp2(c(-2, 0, 4), c("#440154", "white", "#E66101"))

# 5. 定义列注释
cell_colors <- c("Naive CD4 T" = "#4E79A7", 
  "Memory CD4 T" = "#F28E2B", 
  "CD14+ Mono" = "#E15759",
  "B" = "#76B7B2",
  "CD8 T" = "#59A14F",
  "FCGR3A+ Mono" = "#EDC948",
  "NK" = "#B07AA1",
  "DC" = "#FF9DA7",
                 "Platelet" = "#9C755F")

column_ha = HeatmapAnnotation(
  cell_anno = colnames(mat),
  col = list(cell_anno = cell_colors),
  show_legend = TRUE
)

# 6. 定义行注释 - 修正：确保长度与行数一致
# 假设每个细胞类型有5-6个marker基因
gene_groups <- rep(colnames(mat), length.out = nrow(mat))  # 修正：使用length.out确保长度匹配

row_ha = rowAnnotation(
  gene_anno = gene_groups,
  col = list(gene_anno = cell_colors),
  show_legend = FALSE
)

# 7. 绘图
Heatmap(mat_scaled, 
  name = "Expression", 
  col = col_fun,
        cluster_rows = FALSE,      # 通常 Marker 基因按顺序排，不聚类
  cluster_columns = FALSE, 
  show_row_names = TRUE,
  row_names_side = "right",
  top_annotation = column_ha,
  left_annotation = row_ha,
        column_names_rot = 90,
        rect_gp = gpar(col = "white", lwd = 0.5)) # 格子之间的白线