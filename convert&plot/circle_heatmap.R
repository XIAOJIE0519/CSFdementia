library(circlize)
library(ComplexHeatmap)
library(grid)
library(openxlsx)

cat("Starting dual-semicircle heatmap generation...\n")
setwd("F:/1a-EOD-CSF-protein/1a-figure")

# ==================== 数据读取 ====================
cat("Reading data...\n")
module_order <- read.csv("F:/1a-EOD-CSF-protein/wgcna_consensus_main/module_dendrogram_order.csv")
corr_main <- read.csv("F:/1a-EOD-CSF-protein/wgcna_consensus_main/consensus_module_trait_correlations.csv", row.names = 1)
pval_main <- read.csv("F:/1a-EOD-CSF-protein/wgcna_consensus_main/consensus_module_trait_pvalues.csv", row.names = 1)
celltype_enrich <- read.csv("F:/1a-EOD-CSF-protein/wgcna_consensus_main/module_celltype_fisher_enrichment.csv")
preservation <- read.csv("F:/1a-EOD-CSF-protein/wgcna_consensus_test/module_preservation_by_study.csv")
corr_test <- read.csv("F:/1a-EOD-CSF-protein/wgcna_consensus_test/test_studies_consensus_correlations.csv", row.names = 1)
pval_test <- read.csv("F:/1a-EOD-CSF-protein/wgcna_consensus_test/test_studies_consensus_pvalues.csv", row.names = 1)

# 颜色映射
color_map <- c(
  black = "#000000", blue = "#0000FF", brown = "#A52A2A", cyan = "#00FFFF",
  darkgreen = "#006400", darkred = "#8B0000", green = "#00FF00", greenyellow = "#ADFF2F",
  grey = "#BEBEBE", grey60 = "#999999", lightcyan = "#E0FFFF", lightgreen = "#90EE90",
  lightyellow = "#FFFFE0", magenta = "#FF00FF", midnightblue = "#191970", pink = "#FFC0CB",
  purple = "#A020F0", red = "#FF0000", royalblue = "#4169E1", salmon = "#FA8072",
  tan = "#D2B48C", turquoise = "#40E0D0", yellow = "#FFFF00"
)

module_names_full <- c(
  "M1 Proteasome Degradation", "M2 Carbon Metabolism", "M3 HPV-Related Synapse Organization",
  "M4 Axon Guidance", "M5 Extracellular Matrix Adhesion", "M6 Synaptic Adhesion",
  "M7 IgSF CAM Neuronal Projection", "M8 Synaptic Membrane Adhesion", "M9 MAPK Signaling",
  "M10 Actin Cytoskeleton & Focal Adhesion", "M11 Angiogenesis", "M12 Lysosomal Glycan Metabolism",
  "M13 ECM Organization & Protein Digestion", "M14 Complement & Coagulation",
  "M15 Immunoglobulin & Glycosphingolipid", "M16 Cornified Envelope", "M17 Hemoglobin & Nitrogen Metabolism"
)

ordered_modules <- module_order$Module
n_modules <- length(ordered_modules)
module_colors <- sapply(ordered_modules, function(x) color_map[strsplit(x, "_")[[1]][2]])

# ==================== 准备数据 ====================
corr_main <- corr_main[ordered_modules, ]
pval_main <- pval_main[ordered_modules, ]
corr_test <- corr_test[ordered_modules, ]
pval_test <- pval_test[ordered_modules, ]

biomarker_cols <- c("MoCA", "MMSE", "AB42", "tTau", "pTau", "pTau181", "AB42.pTau", "AB40", "NEFL", "YKL40", "pTau217", "pTau231")
available_biomarkers <- biomarker_cols[biomarker_cols %in% colnames(corr_main)]
available_biomarkers_test <- biomarker_cols[biomarker_cols %in% colnames(corr_test)]
cell_types <- c("neurons", "astrocytes", "oligodendrocytes", "microglia", "endothelial", "OPCs")

# Main数据 - 左半圆（Biomarker+细胞类型）
mat_bio <- as.matrix(corr_main[, available_biomarkers])
pval_bio <- as.matrix(pval_main[, available_biomarkers])

mat_cell <- matrix(NA, nrow = n_modules, ncol = length(cell_types))
colnames(mat_cell) <- cell_types
for (i in 1:n_modules) {
  color_name <- gsub("M[0-9]+_", "", ordered_modules[i])
  for (j in 1:length(cell_types)) {
    idx <- which(celltype_enrich$Module == color_name & celltype_enrich$Cell_Type == cell_types[j])
    if (length(idx) > 0) mat_cell[i, j] <- -log10(celltype_enrich$P_value[idx])
  }
}

data_left_main <- cbind(mat_bio, mat_cell)
pval_left_main <- cbind(pval_bio, matrix(NA, nrow = n_modules, ncol = ncol(mat_cell)))

# Main数据 - 右半圆（表型+人口学）
mat_cn <- as.matrix(corr_main[, "CN", drop = FALSE])
mat_vs_cn <- as.matrix(corr_main[, c("EOD_vs_CN", "EOAD_vs_CN")])
mat_vs_other <- as.matrix(corr_main[, c("EOD_vs_Other", "EOAD_vs_Other")])
mat_demo <- as.matrix(corr_main[, c("Age", "Sex")])
pval_cn <- as.matrix(pval_main[, "CN", drop = FALSE])
pval_vs_cn <- as.matrix(pval_main[, c("EOD_vs_CN", "EOAD_vs_CN")])
pval_vs_other <- as.matrix(pval_main[, c("EOD_vs_Other", "EOAD_vs_Other")])
pval_demo <- as.matrix(pval_main[, c("Age", "Sex")])

data_right_main <- cbind(mat_cn, mat_vs_cn, mat_vs_other, mat_demo)
pval_right_main <- cbind(pval_cn, pval_vs_cn, pval_vs_other, pval_demo)

# Test数据 - 左半圆（Biomarker+保守度）
mat_bio_t <- as.matrix(corr_test[, available_biomarkers_test])
pval_bio_t <- as.matrix(pval_test[, available_biomarkers_test])

preservation_studies <- c("study_4", "study_6", "study_7", "study_9")
mat_pres <- matrix(NA, nrow = n_modules, ncol = length(preservation_studies))
colnames(mat_pres) <- preservation_studies
for (i in 1:n_modules) {
  color_name <- strsplit(ordered_modules[i], "_")[[1]][2]
  idx <- which(preservation$Module == color_name)
  if (length(idx) > 0) {
    for (j in 1:length(preservation_studies)) {
      val <- preservation[idx, preservation_studies[j]]
      mat_pres[i, j] <- ifelse(is.na(val) || val < 0, 0, val)
    }
  }
}

data_left_test <- cbind(mat_bio_t, mat_pres)
pval_left_test <- cbind(pval_bio_t, matrix(NA, nrow = n_modules, ncol = ncol(mat_pres)))

# Test数据 - 右半圆（表型+人口学）
mat_cn_t <- as.matrix(corr_test[, "CN", drop = FALSE])
mat_vs_cn_t <- as.matrix(corr_test[, c("EOD_vs_CN", "EOAD_vs_CN", "EODSD_vs_CN", "EODLB_vs_CN", "EOFTD_vs_CN")])
mat_vs_other_t <- as.matrix(corr_test[, c("EOD_vs_Other", "EOAD_vs_Other", "EODSD_vs_Other", "EODLB_vs_Other", "EOFTD_vs_Other")])
mat_demo_t <- as.matrix(corr_test[, c("Age", "Sex")])
pval_cn_t <- as.matrix(pval_test[, "CN", drop = FALSE])
pval_vs_cn_t <- as.matrix(pval_test[, c("EOD_vs_CN", "EOAD_vs_CN", "EODSD_vs_CN", "EODLB_vs_CN", "EOFTD_vs_CN")])
pval_vs_other_t <- as.matrix(pval_test[, c("EOD_vs_Other", "EOAD_vs_Other", "EODSD_vs_Other", "EODLB_vs_Other", "EOFTD_vs_Other")])
pval_demo_t <- as.matrix(pval_test[, c("Age", "Sex")])

data_right_test <- cbind(mat_cn_t, mat_vs_cn_t, mat_vs_other_t, mat_demo_t)
pval_right_test <- cbind(pval_cn_t, pval_vs_cn_t, pval_vs_other_t, pval_demo_t)

# ==================== 双半圆绘图函数 ====================
draw_dual_semicircle <- function(data_left, pval_left, data_right, pval_right, 
                                 module_colors, module_names, title, filename, is_test = FALSE) {
  
  n_modules <- nrow(data_left)
  n_left <- ncol(data_left)
  n_right <- ncol(data_right)
  
  cat("Drawing", n_modules, "modules: left", n_left, "tracks, right", n_right, "tracks\n")
  
  # 创建扇区：每个模块在左右半圆各一个
  # 右半圆反转顺序以实现水平对称
  sectors_left <- paste0("M", 1:n_modules, "_L")
  sectors_right <- paste0("M", n_modules:1, "_R")
  all_sectors <- c(sectors_left, sectors_right)
  
  # 右半圆的数据也需要反转
  data_right <- data_right[n_modules:1, , drop = FALSE]
  pval_right <- pval_right[n_modules:1, , drop = FALSE]
  module_colors_right <- rev(module_colors)
  module_names_right <- rev(module_names)
  
  # 颜色函数
  col_fun1 <- colorRamp2(c(-1, 0, 1), c("#6baed6", "white", "#fc8d59"))
  
  # Biomarker颜色
  breaks_bio <- c(seq(-2.1, -0.101, by = 0.2), 0, seq(0.101, 2.1, by = 0.2))
  colors_bio <- c(colorRampPalette(c("#AAD09D", "#ECF4DD"))(length(seq(-2.1, -0.101, by = 0.2))), 
                  "white", 
                  colorRampPalette(c("#FFF7AC", "#ECB477"))(length(seq(0.101, 2.1, by = 0.2))))
  col_fun3 <- colorRamp2(breaks_bio, colors_bio)
  
  col_fun4 <- colorRamp2(c(0, 1.5, 3), c("white", "#D8BFD8", "#9370DB"))
  col_fun_pres <- colorRamp2(c(0, 1, 2, 10), c("#9EBCC8", "#9DCFC2", "#B5DBBA", "#D3E4A6"))
  
  # 计算track高度
  base_h <- 0.025
  gap_h <- 0.008
  
  # 初始化circos - 左半圆270°-90°，右半圆90°-270°，中间各留40°
  circos.clear()
  gaps <- rep(0.3, length(all_sectors))
  gaps[n_modules] <- 40  # 左半圆到右半圆
  gaps[length(all_sectors)] <- 40  # 右半圆到左半圆
  
  circos.par(start.degree = 90, gap.after = gaps, cell.padding = c(0,0,0,0), track.margin = c(0.001, 0.001))
  
  png(filename, width = 28, height = 28, units = "in", res = 300)
  par(mar = c(2, 2, 2, 2))
  
  circos.initialize(factors = all_sectors, xlim = c(0, 1))
  
  # 计算需要的填充使两边最内层对齐
  max_tracks <- max(n_left, n_right)
  padding_left <- (max_tracks - n_left) * base_h
  padding_right <- (max_tracks - n_right) * base_h
  
  if (padding_left > 0) {
    circos.track(factors = sectors_left, ylim = c(0, 1), track.height = padding_left, 
                 bg.border = NA, panel.fun = function(x, y) {})
  }
  if (padding_right > 0) {
    circos.track(factors = sectors_right, ylim = c(0, 1), track.height = padding_right, 
                 bg.border = NA, panel.fun = function(x, y) {})
  }
  
  # 绘制左半圆tracks
  track_idx_left <- 1
  for (i in 1:n_left) {
    col_name <- colnames(data_left)[i]
    values_left <- data_left[, i]
    pvals_left <- pval_left[, i]
    
    # 判断是Biomarker还是Cell/Preservation
    is_biomarker <- col_name %in% biomarker_cols
    is_cell <- col_name %in% cell_types
    is_pres <- grepl("study_", col_name)
    
    if (is_biomarker) {
      col_fun <- col_fun3
    } else if (is_cell) {
      col_fun <- col_fun4
    } else if (is_pres) {
      col_fun <- col_fun_pres
    } else {
      col_fun <- col_fun1
    }
    
    circos.track(
      factors = sectors_left,
      ylim = c(0, 1),
      track.height = base_h,
      bg.border = NA,
      panel.fun = function(x, y) {
        sector_idx <- get.cell.meta.data("sector.index")
        module_idx <- as.integer(gsub("M|_L", "", sector_idx))
        value <- values_left[module_idx]
        pval <- pvals_left[module_idx]
        
        if (!is.na(value)) {
          xlim <- get.cell.meta.data("xlim")
          ylim <- get.cell.meta.data("ylim")
          circos.rect(xlim[1], ylim[1], xlim[2], ylim[2], col = col_fun(value), border = "grey40", lwd = 0.8)
          
          if (is_pres) {
            if (value >= 10) circos.text(0.5, 0.5, "##", cex = 0.45, col = "black", font = 1)
            else if (value >= 1.96) circos.text(0.5, 0.5, "#", cex = 0.45, col = "black", font = 1)
          } else if (is_cell && !is.na(value) && value >= -log10(0.05)) {
            circos.text(0.5, 0.5, "*", cex = 0.45, col = "black", font = 1)
          } else if (!is.na(pval)) {
            star <- ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "")))
            if (star != "") circos.text(0.5, 0.5, star, cex = 0.45, col = "black", font = 1)
          }
        }
      }
    )
    
    track_idx_left <- track_idx_left + 1
  }
  
  # 绘制右半圆tracks
  track_idx_right <- 1
  for (i in 1:n_right) {
    col_name <- colnames(data_right)[i]
    values_right <- data_right[, i]
    pvals_right <- pval_right[, i]
    
    # 右半圆是表型+人口学，使用col_fun1
    col_fun <- col_fun1
    
    circos.track(
      factors = sectors_right,
      ylim = c(0, 1),
      track.height = base_h,
      bg.border = NA,
      panel.fun = function(x, y) {
        sector_idx <- get.cell.meta.data("sector.index")
        module_idx <- as.integer(gsub("M|_R", "", sector_idx))
        # 注意：右半圆数据已经反转，所以直接用module_idx
        value <- values_right[module_idx]
        pval <- pvals_right[module_idx]
        
        if (!is.na(value)) {
          xlim <- get.cell.meta.data("xlim")
          ylim <- get.cell.meta.data("ylim")
          circos.rect(xlim[1], ylim[1], xlim[2], ylim[2], col = col_fun(value), border = "grey40", lwd = 0.8)
          
          if (!is.na(pval)) {
            star <- ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "")))
            if (star != "") circos.text(0.5, 0.5, star, cex = 0.45, col = "black", font = 1)
          }
        }
      }
    )
    
    track_idx_right <- track_idx_right + 1
  }
  
  # 添加模块颜色track
  all_module_colors <- c(module_colors, module_colors_right)
  circos.track(factors = all_sectors, ylim = c(0, 1), track.height = 0.025, 
               bg.col = all_module_colors, bg.border = "white", bg.lwd = 2,
               panel.fun = function(x, y) {})
  
  # 添加模块名称（与半径平行，垂直于弧）
  # 左半圆
  for (i in 1:n_modules) {
    sector_idx <- paste0("M", i, "_L")
    theta <- mean(circlize:::get.sector.data(sector_idx)[c("start.degree", "end.degree")])
    theta_rad <- theta * pi / 180
    
    radius <- 1.02
    x_pos <- radius * cos(theta_rad)
    y_pos <- radius * sin(theta_rad)
    
    # 与半径平行（垂直于弧）
    # 左半圆：270°-90°，文字方向与半径一致
    if (theta >= 90 && theta <= 270) {
      text_angle <- theta - 90
      text_adj <- c(0.5, 1)
    } else {
      text_angle <- theta + 90
      text_adj <- c(0.5, 0)
    }
    
    text(x_pos, y_pos, module_names[i], srt = text_angle, adj = text_adj, 
         cex = 0.6, col = "grey20", font = 2)
  }
  
  # 右半圆
  for (i in 1:n_modules) {
    sector_idx <- paste0("M", n_modules - i + 1, "_R")
    theta <- mean(circlize:::get.sector.data(sector_idx)[c("start.degree", "end.degree")])
    theta_rad <- theta * pi / 180
    
    radius <- 1.02
    x_pos <- radius * cos(theta_rad)
    y_pos <- radius * sin(theta_rad)
    
    # 与半径平行（垂直于弧）
    # 右半圆：90°-270°，文字方向与半径一致
    if (theta >= 90 && theta <= 270) {
      text_angle <- theta - 90
      text_adj <- c(0.5, 1)
    } else {
      text_angle <- theta + 90
      text_adj <- c(0.5, 0)
    }
    
    text(x_pos, y_pos, module_names_right[i], srt = text_angle, adj = text_adj, 
         cex = 0.6, col = "grey20", font = 2)
  }
  
  # 添加标题
  text(0, 0, title, cex = 2.5, font = 2, col = "grey20")
  
  # 添加图例（底部中间位置）
  lgd1 <- Legend(title = "Module\nCorrelation", col_fun = col_fun1, at = c(-1, 0, 1),
                 direction = "horizontal", title_position = "topcenter", legend_width = unit(3.5, "cm"),
                 title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 9), border = "black")
  
  lgd2 <- Legend(title = "Biomarker\nCorrelation", col_fun = col_fun3, at = c(-1, -0.5, 0, 0.5, 1),
                 direction = "horizontal", title_position = "topcenter", legend_width = unit(3.5, "cm"),
                 title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 9), border = "black")
  
  if (is_test) {
    lgd3 <- Legend(title = "Preservation\nZsummary", col_fun = col_fun_pres, at = c(0, 2, 10),
                   direction = "horizontal", title_position = "topcenter", legend_width = unit(3.5, "cm"),
                   title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 9), border = "black")
  } else {
    lgd3 <- Legend(title = "Cell Type\n-log10(p-value)", col_fun = col_fun4, at = c(0, 1.5, 3),
                   direction = "horizontal", title_position = "topcenter", legend_width = unit(3.5, "cm"),
                   title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 9), border = "black")
  }
  
  draw(lgd1, x = unit(0.5, "npc"), y = unit(0.18, "npc"), just = c("center", "center"))
  draw(lgd2, x = unit(0.5, "npc"), y = unit(0.10, "npc"), just = c("center", "center"))
  draw(lgd3, x = unit(0.5, "npc"), y = unit(0.02, "npc"), just = c("center", "center"))
  
  dev.off()
  circos.clear()
  cat("Saved:", filename, "\n")
}

# ==================== 生成图形 ====================
cat("\nGenerating Training Set...\n")
draw_dual_semicircle(data_left_main, pval_left_main, data_right_main, pval_right_main,
                     module_colors, module_names_full, "Training Set", 
                     "circular_heatmap_main.png", is_test = FALSE)

cat("\nGenerating Test Set...\n")
draw_dual_semicircle(data_left_test, pval_left_test, data_right_test, pval_right_test,
                     module_colors, module_names_full, "Test Set", 
                     "circular_heatmap_test.png", is_test = TRUE)

cat("\nAll done!\n")
