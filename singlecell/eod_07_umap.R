# ============================================================================
# EOD Single-Cell Integration Pipeline - Circular UMAP Visualization
# ============================================================================
# 
# 生成环形UMAP图，展示细胞类型、疾病分组和数据来源的三层注释
# 
# 输入：
#   - EOD/umap_coordinates_metadata.csv (UMAP坐标和元数据)
#   - EOD/pca_harmony_coordinates.csv (可选：Harmony校正的PCA)
# 
# 输出：
#   - EOD/figures/umap_circular_plot.png (环形UMAP图)
# 
# 说明：
#   - UMAP坐标已经基于Harmony校正的PCA计算，去除了批次效应
#   - 细胞按生物学特征（细胞类型、疾病状态）聚集，而非数据源
# 
# Author: Bioinformatics Analysis
# Date: 2026-03-07
# ============================================================================

# 清空环境变量
rm(list = ls())

# 加载必需的包
suppressPackageStartupMessages({
  library(Seurat)
  library(RColorBrewer)
  library(plot1cell)
})

cat("\n")
cat("================================================================================\n")
cat("EOD Pipeline - Circular UMAP Visualization\n")
cat("================================================================================\n\n")

# ============================================================================
# 1. 加载数据
# ============================================================================

cat("Loading data from CSV...\n")

# 从CSV文件加载UMAP坐标和元数据
csv_file <- "EOD/umap_coordinates_metadata.csv"

if (!file.exists(csv_file)) {
  stop("Error: Input file not found: ", csv_file)
}

# 读取CSV文件
umap_data <- read.csv(csv_file, row.names = 1, stringsAsFactors = FALSE)

cat(sprintf("  Loaded: %d cells\n", nrow(umap_data)))
cat("  Note: UMAP coordinates are computed from Harmony-corrected PCA\n")
cat("        Batch effects have been removed during integration\n")

# 提取UMAP坐标
umap_coords <- as.matrix(umap_data[, c("UMAP_1", "UMAP_2")])

# 提取元数据
metadata <- umap_data[, c("cell_type", "disease_group", "data_source", "sample_id")]

# 创建一个最小的表达矩阵（只用于初始化Seurat对象）
n_cells <- nrow(umap_data)
dummy_matrix <- Matrix::Matrix(0, nrow = 100, ncol = n_cells, sparse = TRUE)
rownames(dummy_matrix) <- paste0("Gene_", 1:100)
colnames(dummy_matrix) <- rownames(umap_data)

# 创建Seurat对象
sco <- CreateSeuratObject(counts = dummy_matrix, meta.data = metadata)

# 添加UMAP坐标
sco[["umap"]] <- CreateDimReducObject(embeddings = umap_coords, key = "UMAP_", assay = "RNA")

cat(sprintf("  Created Seurat object: %d cells\n", ncol(sco)))

# 检查是否有Harmony PCA数据（可选，用于高级分析）
pca_harmony_file <- "EOD/pca_harmony_coordinates.csv"
if (file.exists(pca_harmony_file)) {
  cat("\n  Found Harmony-corrected PCA coordinates\n")
  pca_harmony <- read.csv(pca_harmony_file, row.names = 1)
  sco[["pca_harmony"]] <- CreateDimReducObject(
    embeddings = as.matrix(pca_harmony), 
    key = "PC_", 
    assay = "RNA"
  )
  cat(sprintf("    Added %d PCs to Seurat object\n", ncol(pca_harmony)))
}

# ============================================================================
# 2. 检查和准备元数据
# ============================================================================

cat("\nPreparing metadata...\n")

# 检查必需的列
required_cols <- c("cell_type", "disease_group", "data_source")
missing_cols <- setdiff(required_cols, colnames(sco@meta.data))

if (length(missing_cols) > 0) {
  stop("Error: Missing required columns in metadata: ", paste(missing_cols, collapse = ", "))
}

# 确保元数据列是factor类型
sco$cell_type <- factor(sco$cell_type)
sco$disease_group <- factor(sco$disease_group, levels = c("CN", "EOAD", "LOAD", "EOFTD", "LOFTD"))
sco$data_source <- factor(sco$data_source)

# 设置细胞类型为默认分组（必须是factor）
Idents(sco) <- sco$cell_type

# 计算每个细胞类型的细胞数量
cell_type_counts <- table(sco$cell_type)

# 创建带细胞数量的标签
cell_type_labels <- paste0(names(cell_type_counts), " (", cell_type_counts, ")")
names(cell_type_labels) <- names(cell_type_counts)

# 验证Idents设置
cat("  Checking Idents...\n")
cat(sprintf("    Class: %s\n", class(Idents(sco))))
cat(sprintf("    Has levels: %s\n", !is.null(levels(Idents(sco)))))
cat(sprintf("    Number of levels: %d\n", length(levels(Idents(sco)))))

cat("\n  Cell types with counts:\n")
for (ct in names(cell_type_counts)) {
  cat(sprintf("    %s: %d cells\n", ct, cell_type_counts[ct]))
}

cat("\n  Cell types:\n")
print(table(sco$cell_type))

cat("\n  Disease groups:\n")
print(table(sco$disease_group))

cat("\n  Data sources:\n")
print(table(sco$data_source))

# ============================================================================
# 3. 定义颜色方案
# ============================================================================

cat("\nDefining color schemes...\n")

# 细胞类型颜色（与Python保持一致）
cell_type_colors <- c(
  "ExcitatoryNeurons" = "#1f77b4",
  "InhibitoryNeurons" = "#ff7f0e",
  "Astrocytes" = "#2ca02c",
  "Oligodendrocytes" = "#d62728",
  "OPCs" = "#9467bd",
  "Microglia" = "#8c564b",
  "Macrophages" = "#e377c2",
  "EndothelialCells" = "#7f7f7f",
  "Pericytes" = "#bcbd22",
  "Fibroblasts" = "#17becf"
)

# 疾病分组颜色
disease_colors <- c(
  "CN" = "#2ca02c",
  "EOAD" = "#d62728",
  "LOAD" = "#ff7f0e",
  "EOFTD" = "#9467bd",
  "LOFTD" = "#8c564b"
)

# 数据来源颜色
data_sources <- unique(sco$data_source)
n_sources <- length(data_sources)
source_colors <- colorRampPalette(brewer.pal(n = min(12, n_sources), name = "Set3"))(n_sources)
names(source_colors) <- data_sources

# 提取对应的颜色向量
cell_types_present <- levels(sco$cell_type)
cluster_colors <- cell_type_colors[cell_types_present]

# 移除NA值
cluster_colors <- cluster_colors[!is.na(cluster_colors)]

cat(sprintf("  Cell types with colors: %d\n", length(cluster_colors)))

# ============================================================================
# 4. 准备环形图数据
# ============================================================================

cat("\nPreparing circlize data...\n")

# 检查可用的reductions
cat("  Available reductions:\n")
print(names(sco@reductions))

# 确保UMAP reduction存在
if (!"umap" %in% names(sco@reductions)) {
  cat("  Warning: 'umap' reduction not found, checking alternatives...\n")
  
  # 检查是否有其他UMAP相关的reduction
  umap_reductions <- grep("umap", names(sco@reductions), ignore.case = TRUE, value = TRUE)
  
  if (length(umap_reductions) > 0) {
    cat(sprintf("  Found alternative UMAP reduction: %s\n", umap_reductions[1]))
    # 重命名为标准的"umap"
    sco@reductions[["umap"]] <- sco@reductions[[umap_reductions[1]]]
  } else {
    stop("Error: No UMAP reduction found in Seurat object")
  }
}

# 转换Seurat对象为环状图所需数据格式
# prepare_circlize_data需要Idents是factor且有levels
cat("  Calling prepare_circlize_data...\n")
circ_data <- prepare_circlize_data(sco, scale = 0.8)

cat("  Data prepared successfully\n")

# ============================================================================
# 5. 绘制环形UMAP图
# ============================================================================

cat("\nGenerating circular UMAP plot...\n")

# 设置随机种子以保证可重复性
set.seed(1234)

# 创建输出目录
output_dir <- "EOD/figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 输出文件路径
output_file <- file.path(output_dir, "umap_circular_plot.png")
output_pdf  <- file.path(output_dir, "umap_circular_plot.pdf")

# 打开PNG设备
png(filename = output_file, width = 12, height = 12, units = 'in', res = 300)

# 绘制环形UMAP主图
plot_circlize(circ_data,
              do.label = TRUE,           # 显示细胞类型标签
              pt.size = 0.01,            # 点大小
              col.use = cluster_colors,  # 细胞类型颜色
              bg.color = 'white',        # 背景颜色
              kde2d.n = 200,             # KDE密度估计分辨率
              repel = TRUE,              # 标签避免重叠
              label.cex = 0.6)           # 标签字体大小

# 添加外围信息轨道：第一层 - 疾病分组
add_track(circ_data, 
          group = "disease_group", 
          colors = disease_colors, 
          track_num = 2)

# 添加外围信息轨道：第二层 - 数据来源
add_track(circ_data, 
          group = "data_source", 
          colors = source_colors, 
          track_num = 3)

# 添加图例：细胞类型（带细胞数量）
legend("topright", 
       legend = cell_type_labels[names(cluster_colors)],
       fill = cluster_colors,
       title = "Cell Type (n cells)",
       cex = 0.7,
       bty = "n")

# 添加图例：疾病分组
legend("topleft",
       legend = names(disease_colors),
       fill = disease_colors,
       title = "Disease Group",
       cex = 0.7,
       bty = "n")

# 添加图例：数据来源
legend("bottomleft",
       legend = names(source_colors),
       fill = source_colors,
       title = "Data Source",
       cex = 0.7,
       bty = "n")

# 关闭设备
dev.off()

cat(sprintf("\n  Saved: %s\n", output_file))

# ---- 同样内容输出为 PDF ----
pdf(file = output_pdf, width = 12, height = 12)

plot_circlize(circ_data,
              do.label = TRUE,
              pt.size = 0.01,
              col.use = cluster_colors,
              bg.color = 'white',
              kde2d.n = 200,
              repel = TRUE,
              label.cex = 0.6)

add_track(circ_data,
          group = "disease_group",
          colors = disease_colors,
          track_num = 2)

add_track(circ_data,
          group = "data_source",
          colors = source_colors,
          track_num = 3)

legend("topright",
       legend = cell_type_labels[names(cluster_colors)],
       fill = cluster_colors,
       title = "Cell Type (n cells)",
       cex = 0.7,
       bty = "n")

legend("topleft",
       legend = names(disease_colors),
       fill = disease_colors,
       title = "Disease Group",
       cex = 0.7,
       bty = "n")

legend("bottomleft",
       legend = names(source_colors),
       fill = source_colors,
       title = "Data Source",
       cex = 0.7,
       bty = "n")

dev.off()

cat(sprintf("  Saved: %s\n", output_pdf))

# ============================================================================
# 6. 输出统计信息
# ============================================================================

cat("\n")
cat("================================================================================\n")
cat("Summary Statistics\n")
cat("================================================================================\n\n")

cat(sprintf("Total cells: %d\n", ncol(sco)))
cat(sprintf("Cell types: %d\n", length(unique(sco$cell_type))))
cat(sprintf("Disease groups: %d\n", length(unique(sco$disease_group))))
cat(sprintf("Data sources: %d\n", length(unique(sco$data_source))))

cat("\n")
cat("================================================================================\n")
cat("Circular UMAP visualization complete!\n")
cat("================================================================================\n\n")

cat("Output files:\n")
cat(sprintf("  - %s\n", output_file))
cat(sprintf("  - %s\n", output_pdf))
cat("\nThe circular UMAP shows:\n")
cat("  - Inner scatter: UMAP coordinates colored by cell type\n")
cat("  - Outer ring 1 (innermost): Cell type labels with KDE contours\n")
cat("  - Outer ring 2 (middle): Disease group composition\n")
cat("  - Outer ring 3 (outermost): Data source composition\n")

cat("\n")
cat("================================================================================\n")
cat("Done!\n")
cat("================================================================================\n\n")
