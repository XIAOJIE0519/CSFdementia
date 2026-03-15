# ============================================================================
# 蛋白质互作网络可视化
# 基于STRING数据库的模块内/模块间蛋白质互作关系
# 使用 ggraph + geom_mark_hull 实现参考图风格
# ============================================================================

# 加载必需的R包
suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(ggraph)
  library(ggforce)
  library(tidygraph)
  library(data.table)
})

# ============================================================================
# 1. 设置工作目录和读取数据
# ============================================================================

work_dir <- "F:/1a-EOD-CSF-protein"
output_dir <- "F:/1a-EOD-CSF-protein/1a-figure"
setwd(work_dir)

cat("读取数据文件...\n")

# 读取STRING互作数据（已过滤）
ppi_data <- fread("module_correlations/module_protein_correlations_STRING_filtered.csv", 
                  data.table = FALSE, stringsAsFactors = FALSE)

# 读取模块-蛋白-基因对应关系
module_protein_gene <- fread("module_correlations/module_protein_gene_summary.csv", 
                             data.table = FALSE, stringsAsFactors = FALSE)

cat("数据读取完成！\n")
cat("  互作关系数:", nrow(ppi_data), "\n")
cat("  蛋白质数:", length(unique(c(ppi_data$Protein1, ppi_data$Protein2))), "\n\n")

# ============================================================================
# 2. 定义模块颜色映射（与热图保持一致）
# ============================================================================

module_color_map <- c(
  M1_black = "#000000", 
  M2_blue = "#0000FF", 
  M3_brown = "#A52A2A", 
  M4_cyan = "#00FFFF",
  M5_darkgreen = "#006400", 
  M6_darkred = "#8B0000", 
  M7_green = "#00FF00", 
  M8_greenyellow = "#ADFF2F",
  M9_grey = "#BEBEBE", 
  M9_grey60 = "#999999", 
  M10_lightcyan = "#E0FFFF", 
  M11_lightgreen = "#90EE90",
  M12_lightyellow = "#FFFFE0", 
  M13_magenta = "#FF00FF", 
  M14_midnightblue = "#191970", 
  M15_pink = "#FFC0CB",
  M16_purple = "#A020F0", 
  M17_red = "#FF0000", 
  M18_royalblue = "#4169E1", 
  M19_salmon = "#FA8072",
  M20_tan = "#D2B48C", 
  M21_turquoise = "#40E0D0", 
  M22_yellow = "#FFFF00"
)

# ============================================================================
# 3. 准备网络数据
# ============================================================================

cat("准备网络数据...\n")

# 创建蛋白质ID到基因名和模块的映射
protein_to_gene <- module_protein_gene %>%
  select(Protein, Gene, Module) %>%
  distinct()

# ============================================================================
# 4. 构建蛋白质ID网络
# ============================================================================

cat("构建蛋白质ID网络...\n")

# 边表
edges_protein <- data.frame(
  from = ppi_data$Protein1,
  to = ppi_data$Protein2,
  weight = ppi_data$STRING_Weight,
  stringsAsFactors = FALSE
)

# 节点表
nodes_protein <- data.frame(
  name = unique(c(edges_protein$from, edges_protein$to)),
  stringsAsFactors = FALSE
) %>%
  left_join(protein_to_gene %>% select(Protein, Module), 
            by = c("name" = "Protein")) %>%
  mutate(Module = ifelse(is.na(Module), "Unknown", Module)) %>%
  filter(Module != "Unknown") %>%  # 移除未知模块的节点
  distinct(name, .keep_all = TRUE)

# 过滤边表，只保留有模块信息的节点
edges_protein <- edges_protein %>%
  filter(from %in% nodes_protein$name & to %in% nodes_protein$name)

# 创建igraph对象
net_protein <- graph_from_data_frame(
  d = edges_protein,
  directed = FALSE,
  vertices = nodes_protein
)

# 转换为tidygraph对象
tg_protein <- as_tbl_graph(net_protein)

cat("  蛋白质网络节点数:", vcount(net_protein), "\n")
cat("  蛋白质网络边数:", ecount(net_protein), "\n")

# ============================================================================
# 5. 构建基因名网络
# ============================================================================

cat("构建基因名网络...\n")

# 为边表添加基因名
edges_gene <- edges_protein %>%
  left_join(protein_to_gene %>% select(Protein, Gene), 
            by = c("from" = "Protein")) %>%
  rename(from_gene = Gene) %>%
  left_join(protein_to_gene %>% select(Protein, Gene), 
            by = c("to" = "Protein")) %>%
  rename(to_gene = Gene) %>%
  filter(!is.na(from_gene) & !is.na(to_gene) & 
         from_gene != "" & to_gene != "") %>%
  select(from = from_gene, to = to_gene, weight)

# 节点表
nodes_gene <- data.frame(
  name = unique(c(edges_gene$from, edges_gene$to)),
  stringsAsFactors = FALSE
) %>%
  left_join(protein_to_gene %>% select(Gene, Module), 
            by = c("name" = "Gene")) %>%
  mutate(Module = ifelse(is.na(Module), "Unknown", Module)) %>%
  filter(Module != "Unknown") %>%
  distinct(name, .keep_all = TRUE)

# 过滤边表
edges_gene <- edges_gene %>%
  filter(from %in% nodes_gene$name & to %in% nodes_gene$name)

# 创建igraph对象
net_gene <- graph_from_data_frame(
  d = edges_gene,
  directed = FALSE,
  vertices = nodes_gene
)

# 转换为tidygraph对象
tg_gene <- as_tbl_graph(net_gene)

cat("  基因网络节点数:", vcount(net_gene), "\n")
cat("  基因网络边数:", ecount(net_gene), "\n\n")

# ============================================================================
# 6. 优化布局：调整边权重实现模块聚集
# ============================================================================

cat("重新计算布局引力权重...\n")

# 提取节点模块信息
node_info_protein <- nodes_protein %>% select(name, Module)

# 为边表打上模块标签，并计算新的布局权重
edges_for_layout_protein <- edges_protein %>%
  left_join(node_info_protein, by = c("from" = "name")) %>%
  rename(from_module = Module) %>%
  left_join(node_info_protein, by = c("to" = "name")) %>%
  rename(to_module = Module) %>%
  mutate(
    # 核心魔法：同模块边权重放大20倍，跨模块边权重缩小到0.01倍（更强的聚集效果）
    layout_weight = ifelse(from_module == to_module, weight * 20, weight * 0.01),
    # 标记边类型
    edge_type = ifelse(from_module == to_module, "Intra", "Inter")
  )

cat("  模块内边数:", sum(edges_for_layout_protein$edge_type == "Intra"), "\n")
cat("  跨模块边数:", sum(edges_for_layout_protein$edge_type == "Inter"), "\n")

# 重建带有layout_weight的网络对象
net_protein_opt <- graph_from_data_frame(
  d = edges_for_layout_protein,
  directed = FALSE,
  vertices = nodes_protein
)
tg_protein_opt <- as_tbl_graph(net_protein_opt)

# ============================================================================
# 7. 使用ggraph绘制蛋白质ID网络（优化版）
# ============================================================================

cat("绘制蛋白质ID网络图（优化布局）...\n")

set.seed(123)

# 准备颜色向量
modules_in_net <- unique(nodes_protein$Module)
color_vector <- setNames(module_color_map[modules_in_net], modules_in_net)

# 使用FR布局 + 自定义权重 + 增加迭代次数
layout_protein <- create_layout(
  tg_protein_opt, 
  layout = 'fr',
  weights = E(tg_protein_opt)$layout_weight,
  niter = 3000  # 增加到3000次迭代，让模块更紧密
)

# 过滤有效模块（>=3个节点）
module_counts <- table(layout_protein$Module)
valid_modules <- names(module_counts[module_counts >= 3])
layout_protein_filtered <- layout_protein %>%
  filter(Module %in% valid_modules)

cat("  有效模块数（>=3个节点）:", length(valid_modules), "\n")

# 绘制网络
tryCatch({
  p_protein <- ggraph(layout_protein) + 
    
    # 1. 绘制背景色块（只有边框，没有填充）
    geom_mark_hull(
      data = layout_protein_filtered,
      aes(x = x, y = y, group = Module, color = Module, label = Module),
      fill = NA,               # 关键：不填充背景色
      concavity = 2,           
      expand = unit(8, "mm"),  
      label.fontsize = 11,     
      label.fill = "white",    
      label.colour = "black",  
      label.fontface = "bold", 
      con.cap = 0,             
      radius = unit(3, "mm")   
    ) +
    
    # 2. 绘制网络连线（跨模块边用深色）
    geom_edge_link(aes(alpha = edge_type, width = edge_type, color = edge_type)) +
    scale_edge_alpha_manual(values = c("Intra" = 0.4, "Inter" = 0.6)) +
    scale_edge_width_manual(values = c("Intra" = 0.6, "Inter" = 1.2)) +
    scale_edge_color_manual(values = c("Intra" = "#CCCCCC", "Inter" = "#333333")) +
    
    # 3. 绘制节点
    geom_node_point(aes(color = Module), size = 5, show.legend = FALSE) +
    
    # 4. 添加节点文字标签
    geom_node_text(aes(label = name), repel = TRUE, size = 2.5, color = "black", 
                   fontface = "bold", bg.color = "white", bg.r = 0.1, 
                   max.overlaps = 50, force = 1.5, point.padding = 0.1) +
    
    # 5. 设置配色方案
    scale_color_manual(values = color_vector) +
    
    # 6. 主题设置
    theme_graph() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(title = "Protein-Protein Interaction Network (Protein ID)")
  
  # 保存图片
  ggsave(
    filename = file.path(output_dir, "protein_interaction_network_proteinID.png"),
    plot = p_protein,
    width = 20,
    height = 14,
    dpi = 300,
    bg = "white"
  )
  
  cat("  蛋白质ID网络图已保存\n")
  
}, error = function(e) {
  cat("  警告: geom_mark_hull出错，尝试不使用背景色块绘图...\n")
  cat("  错误信息:", conditionMessage(e), "\n")
  
  # 备用方案
  p_protein_simple <- ggraph(layout_protein) + 
    geom_edge_link(aes(alpha = edge_type, width = edge_type, color = edge_type)) +
    scale_edge_alpha_manual(values = c("Intra" = 0.4, "Inter" = 0.6)) +
    scale_edge_width_manual(values = c("Intra" = 0.6, "Inter" = 1.2)) +
    scale_edge_color_manual(values = c("Intra" = "#CCCCCC", "Inter" = "#333333")) +
    geom_node_point(aes(color = Module), size = 5, show.legend = TRUE) +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.5, color = "black", 
                   fontface = "bold", bg.color = "white", bg.r = 0.1, max.overlaps = 50) +
    scale_color_manual(values = color_vector) +
    theme_graph() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
    ) +
    labs(title = "Protein-Protein Interaction Network (Protein ID)", color = "Module")
  
  ggsave(
    filename = file.path(output_dir, "protein_interaction_network_proteinID.png"),
    plot = p_protein_simple,
    width = 20,
    height = 14,
    dpi = 300,
    bg = "white"
  )
  
  cat("  蛋白质ID网络图已保存（无背景色块版本）\n")
})

# ============================================================================
# 8. 优化布局：调整边权重实现模块聚集（基因网络）
# ============================================================================

cat("重新计算基因网络布局引力权重...\n")

# 提取节点模块信息
node_info_gene <- nodes_gene %>% select(name, Module)

# 为边表打上模块标签，并计算新的布局权重
edges_for_layout_gene <- edges_gene %>%
  left_join(node_info_gene, by = c("from" = "name")) %>%
  rename(from_module = Module) %>%
  left_join(node_info_gene, by = c("to" = "name")) %>%
  rename(to_module = Module) %>%
  mutate(
    # 同模块边权重放大20倍，跨模块边权重缩小到0.01倍（更强的聚集效果）
    layout_weight = ifelse(from_module == to_module, weight * 20, weight * 0.01),
    # 标记边类型
    edge_type = ifelse(from_module == to_module, "Intra", "Inter")
  )

cat("  模块内边数:", sum(edges_for_layout_gene$edge_type == "Intra"), "\n")
cat("  跨模块边数:", sum(edges_for_layout_gene$edge_type == "Inter"), "\n")

# 重建带有layout_weight的网络对象
net_gene_opt <- graph_from_data_frame(
  d = edges_for_layout_gene,
  directed = FALSE,
  vertices = nodes_gene
)
tg_gene_opt <- as_tbl_graph(net_gene_opt)

# ============================================================================
# 9. 使用ggraph绘制基因名网络（优化版）
# ============================================================================

cat("绘制基因名网络图（优化布局）...\n")

set.seed(123)

# 准备颜色向量
modules_in_gene_net <- unique(nodes_gene$Module)
color_vector_gene <- setNames(module_color_map[modules_in_gene_net], modules_in_gene_net)

# 使用FR布局 + 自定义权重 + 增加迭代次数
layout_gene <- create_layout(
  tg_gene_opt, 
  layout = 'fr',
  weights = E(tg_gene_opt)$layout_weight,
  niter = 3000  # 增加到3000次迭代
)

# 过滤有效模块（>=3个节点）
module_counts_gene <- table(layout_gene$Module)
valid_modules_gene <- names(module_counts_gene[module_counts_gene >= 3])
layout_gene_filtered <- layout_gene %>%
  filter(Module %in% valid_modules_gene)

cat("  有效模块数（>=3个节点）:", length(valid_modules_gene), "\n")

# 绘制网络
tryCatch({
  p_gene <- ggraph(layout_gene) + 
    
    # 1. 绘制背景色块（只有边框，没有填充）
    geom_mark_hull(
      data = layout_gene_filtered,
      aes(x = x, y = y, group = Module, color = Module, label = Module),
      fill = NA,               # 关键：不填充背景色
      concavity = 2,           
      expand = unit(8, "mm"),  
      label.fontsize = 11,     
      label.fill = "white",    
      label.colour = "black",  
      label.fontface = "bold", 
      con.cap = 0,             
      radius = unit(3, "mm")   
    ) +
    
    # 2. 绘制网络连线（跨模块边用深色）
    geom_edge_link(aes(alpha = edge_type, width = edge_type, color = edge_type)) +
    scale_edge_alpha_manual(values = c("Intra" = 0.4, "Inter" = 0.6)) +
    scale_edge_width_manual(values = c("Intra" = 0.6, "Inter" = 1.2)) +
    scale_edge_color_manual(values = c("Intra" = "#CCCCCC", "Inter" = "#333333")) +
    
    # 3. 绘制节点
    geom_node_point(aes(color = Module), size = 5, show.legend = FALSE) +
    
    # 4. 添加节点文字标签
    geom_node_text(aes(label = name), repel = TRUE, size = 2.5, color = "black", 
                   fontface = "bold", bg.color = "white", bg.r = 0.1, 
                   max.overlaps = 50, force = 1.5, point.padding = 0.1) +
    
    # 5. 设置配色方案
    scale_color_manual(values = color_vector_gene) +
    
    # 6. 主题设置
    theme_graph() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(title = "Protein-Protein Interaction Network (Gene Name)")
  
  # 保存图片
  ggsave(
    filename = file.path(output_dir, "protein_interaction_network_geneName.png"),
    plot = p_gene,
    width = 20,
    height = 14,
    dpi = 300,
    bg = "white"
  )
  
  cat("  基因名网络图已保存\n")
  
}, error = function(e) {
  cat("  警告: geom_mark_hull出错，尝试不使用背景色块绘图...\n")
  cat("  错误信息:", conditionMessage(e), "\n")
  
  # 备用方案
  p_gene_simple <- ggraph(layout_gene) + 
    geom_edge_link(aes(alpha = edge_type, width = edge_type, color = edge_type)) +
    scale_edge_alpha_manual(values = c("Intra" = 0.4, "Inter" = 0.6)) +
    scale_edge_width_manual(values = c("Intra" = 0.6, "Inter" = 1.2)) +
    scale_edge_color_manual(values = c("Intra" = "#CCCCCC", "Inter" = "#333333")) +
    geom_node_point(aes(color = Module), size = 5, show.legend = TRUE) +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.5, color = "black", 
                   fontface = "bold", bg.color = "white", bg.r = 0.1, max.overlaps = 50) +
    scale_color_manual(values = color_vector_gene) +
    theme_graph() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
    ) +
    labs(title = "Protein-Protein Interaction Network (Gene Name)", color = "Module")
  
  ggsave(
    filename = file.path(output_dir, "protein_interaction_network_geneName.png"),
    plot = p_gene_simple,
    width = 20,
    height = 14,
    dpi = 300,
    bg = "white"
  )
  
  cat("  基因名网络图已保存（无背景色块版本）\n")
})

# ============================================================================
# 10. 输出统计摘要
# ============================================================================

cat("\n=== 绘图完成 ===\n")
cat("输出目录:", output_dir, "\n")
cat("已生成两张网络图:\n")
cat("  1. protein_interaction_network_proteinID.png (蛋白质ID标注)\n")
cat("  2. protein_interaction_network_geneName.png (基因名标注)\n\n")

cat("网络特征:\n")
cat("  - 图片尺寸: 20x14英寸\n")
cat("  - 分辨率: 300 DPI\n")
cat("  - 布局算法: Fruchterman-Reingold (FR) + 权重优化\n")
cat("  - 节点颜色: 模块颜色实心填充\n")
cat("  - 模块背景: geom_mark_hull，透明度15%\n")
cat("  - 边权重策略: 模块内x10，跨模块x0.05\n")
cat("  - 边视觉效果: 模块内边明显，跨模块边弱化\n\n")

cat("蛋白质网络统计:\n")
cat("  - 节点数:", vcount(net_protein), "\n")
cat("  - 边数:", ecount(net_protein), "\n")
cat("  - 模块数:", length(unique(nodes_protein$Module)), "\n")
cat("  - 模块内边:", sum(edges_for_layout_protein$edge_type == "Intra"), "\n")
cat("  - 跨模块边:", sum(edges_for_layout_protein$edge_type == "Inter"), "\n\n")

cat("基因网络统计:\n")
cat("  - 节点数:", vcount(net_gene), "\n")
cat("  - 边数:", ecount(net_gene), "\n")
cat("  - 模块数:", length(unique(nodes_gene$Module)), "\n")
cat("  - 模块内边:", sum(edges_for_layout_gene$edge_type == "Intra"), "\n")
cat("  - 跨模块边:", sum(edges_for_layout_gene$edge_type == "Inter"), "\n\n")

cat("优化策略:\n")
cat("  ✓ 引力操控: 同模块节点强吸引(x20)，跨模块节点弱吸引(x0.01)\n")
cat("  ✓ 边权重调整: 模块内边x20，跨模块边x0.01\n")
cat("  ✓ 视觉强化: 跨模块边深色(#333333)、更粗(1.2)、更明显\n")
cat("  ✓ 模块内边: 浅色(#CCCCCC)、较细(0.6)、较淡\n")
cat("  ✓ 迭代优化: FR布局迭代3000次，确保模块紧密聚集\n")
cat("  ✓ 色块优化: 只显示边框，无背景填充\n\n")

cat("所有任务已完成！\n")
