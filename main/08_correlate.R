#!/usr/bin/env Rscript
# PPI Network Visualization for WGCNA Hub Genes
# Filters hub genes (kME > 0.7), identifies significant correlations in both EOD and CN,
# validates with STRING database, and creates publication-quality network visualization

# Install required packages if not available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_packages <- c("igraph", "ggraph", "shadowtext")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(RColorBrewer)
  library(shadowtext)
})

cat("================================================================================\n")
cat("PPI NETWORK VISUALIZATION FOR WGCNA HUB GENES\n")
cat("================================================================================\n\n")

# ================================================================================
# 1. Load and filter hub genes (kME > 0.7), select top 10 per module
# ================================================================================
cat("Step 1: Loading hub genes with kME > 0.7...\n")

hub_genes <- read_csv("wgcna_consensus_main/module_hub_genes.csv", 
                      show_col_types = FALSE)

# Filter kME > 0.7
hub_genes_filtered <- hub_genes %>%
  filter(kME > 0.7) %>%
  select(Protein, Module, Module_Name, kME)

cat(sprintf("  Total hub genes: %d\n", nrow(hub_genes)))
cat(sprintf("  Hub genes with kME > 0.7: %d\n", nrow(hub_genes_filtered)))

# Extract gene symbols from Protein IDs (format: SYMBOL|UNIPROT or just SYMBOL)
hub_genes_filtered <- hub_genes_filtered %>%
  mutate(Gene_Symbol = str_extract(Protein, "^[^|]+"))

# Select top 10 genes per module by kME
hub_genes_filtered <- hub_genes_filtered %>%
  group_by(Module_Name) %>%
  arrange(desc(kME)) %>%
  slice_head(n = 10) %>%
  ungroup()

cat(sprintf("  After selecting top 10 per module: %d\n", nrow(hub_genes_filtered)))
cat(sprintf("  Unique gene symbols: %d\n\n", 
            length(unique(hub_genes_filtered$Gene_Symbol))))

# ================================================================================
# 2. Load correlation data from EOD and CN
# ================================================================================
cat("Step 2: Loading correlation meta-analysis data...\n")

# Function to load and filter correlation data
load_correlations <- function(file_path, condition_name) {
  cat(sprintf("  Loading %s correlations...\n", condition_name))
  
  # Read in chunks due to large file size
  chunk_size <- 100000
  filtered_chunks <- list()
  total_rows <- 0
  first_chunk <- TRUE
  
  # Get hub gene symbols for filtering
  hub_symbols <- unique(hub_genes_filtered$Gene_Symbol)
  
  # Process chunks
  repeat {
    chunk <- tryCatch({
      if (first_chunk) {
        read_csv(file_path, 
                 n_max = chunk_size,
                 show_col_types = FALSE)
      } else {
        read_csv(file_path, 
                 skip = total_rows + 1,
                 n_max = chunk_size,
                 col_names = FALSE,
                 show_col_types = FALSE)
      }
    }, error = function(e) {
      cat(sprintf("    Error reading chunk: %s\n", e$message))
      return(NULL)
    })
    
    if (is.null(chunk) || nrow(chunk) == 0) break
    
    # For subsequent chunks, add column names from first chunk
    if (!first_chunk && ncol(chunk) >= 2) {
      # Get column names from first read
      first_read <- read_csv(file_path, n_max = 1, show_col_types = FALSE)
      if (ncol(chunk) == ncol(first_read)) {
        colnames(chunk) <- colnames(first_read)
      }
    }
    
    first_chunk <- FALSE
    total_rows <- total_rows + nrow(chunk)
    
    # Check if required columns exist
    if (!all(c("Protein1", "Protein2", "Pearson_meta", "P_value", "N_Studies") %in% colnames(chunk))) {
      cat(sprintf("    WARNING: Missing required columns in chunk\n"))
      break
    }
    
    # Convert numeric columns
    chunk <- chunk %>%
      mutate(
        Pearson_meta = as.numeric(Pearson_meta),
        P_value = as.numeric(P_value),
        N_Studies = as.integer(N_Studies)
      )
    
    # Extract gene symbols
    chunk <- chunk %>%
      mutate(
        Gene1 = str_extract(Protein1, "^[^|]+"),
        Gene2 = str_extract(Protein2, "^[^|]+")
      )
    
    # Filter to hub genes only
    chunk_filtered <- chunk %>%
      filter(Gene1 %in% hub_symbols & Gene2 %in% hub_symbols)
    
    if (nrow(chunk_filtered) > 0) {
      filtered_chunks[[length(filtered_chunks) + 1]] <- chunk_filtered
    }
    
    if (total_rows %% 1000000 == 0) {
      cat(sprintf("    Processed %s rows...\n", format(total_rows, big.mark = ",")))
    }
    
    # Limit to reasonable number of rows to avoid memory issues
    if (total_rows >= 10000000) {
      cat(sprintf("    Reached 10M rows limit, stopping...\n"))
      break
    }
  }
  
  if (length(filtered_chunks) == 0) {
    cat(sprintf("    WARNING: No correlations found for hub genes in %s\n", condition_name))
    return(tibble())
  }
  
  df_all <- bind_rows(filtered_chunks)
  
  cat(sprintf("    Total rows processed: %s\n", format(total_rows, big.mark = ",")))
  cat(sprintf("    Correlations for hub genes: %s\n\n", format(nrow(df_all), big.mark = ",")))
  
  return(df_all)
}

# Load EOD correlations
eod_corr <- load_correlations(
  "correlation_meta/EOD/correlation_meta_analysis.csv",
  "EOD"
)

# Load CN correlations
cn_corr <- load_correlations(
  "correlation_meta/CN/correlation_meta_analysis.csv",
  "CN"
)

# ================================================================================
# 3. Filter significant correlations in both EOD and CN
# ================================================================================
cat("Step 3: Filtering significant correlations (FDR < 0.01, |r| > 0.5)...\n")

# Filter EOD: FDR < 0.01, |r| > 0.5
eod_sig <- eod_corr %>%
  filter(FDR_BH_Stratified < 0.01, abs(Pearson_meta) > 0.5) %>%
  mutate(Pair = paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "_")) %>%
  select(Pair, Gene1, Gene2, Protein1, Protein2, 
         EOD_r = Pearson_meta, EOD_FDR = FDR_BH_Stratified, EOD_N = N_Studies)

cat(sprintf("  EOD significant pairs: %d\n", nrow(eod_sig)))

# Filter CN: FDR < 0.01, |r| > 0.5
cn_sig <- cn_corr %>%
  filter(FDR_BH_Stratified < 0.01, abs(Pearson_meta) > 0.5) %>%
  mutate(Pair = paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "_")) %>%
  select(Pair, Gene1, Gene2, Protein1, Protein2,
         CN_r = Pearson_meta, CN_FDR = FDR_BH_Stratified, CN_N = N_Studies)

cat(sprintf("  CN significant pairs: %d\n", nrow(cn_sig)))

# Find pairs significant in BOTH EOD and CN
both_sig <- inner_join(eod_sig, cn_sig, by = "Pair") %>%
  select(Pair, 
         Gene1 = Gene1.x, Gene2 = Gene2.x,
         Protein1 = Protein1.x, Protein2 = Protein2.x,
         EOD_r, EOD_FDR, EOD_N,
         CN_r, CN_FDR, CN_N)

cat(sprintf("  Pairs significant in BOTH EOD and CN: %d\n\n", nrow(both_sig)))

if (nrow(both_sig) == 0) {
  cat("ERROR: No pairs found significant in both conditions\n")
  cat("Exiting...\n")
  quit(status = 1)
}

# ================================================================================
# 4. Load local STRING database for PPI validation
# ================================================================================
cat("Step 4: Loading local STRING database for PPI validation...\n")

# Load local STRING file
string_file <- "STRING.csv"
if (!file.exists(string_file)) {
  cat(sprintf("ERROR: STRING file not found: %s\n", string_file))
  cat("Exiting...\n")
  quit(status = 1)
}

string_data <- read_csv(string_file, show_col_types = FALSE)
cat(sprintf("  Loaded STRING database: %s interactions\n", 
            format(nrow(string_data), big.mark = ",")))

# Get unique genes
unique_genes <- unique(c(both_sig$Gene1, both_sig$Gene2))
cat(sprintf("  Unique genes to validate: %d\n", length(unique_genes)))

# Filter STRING to our genes and confidence > 0.5
ppi_edges <- string_data %>%
  filter(
    source_gene %in% unique_genes & target_gene %in% unique_genes,
    weight > 0.5
  ) %>%
  mutate(
    Gene1 = source_gene,
    Gene2 = target_gene,
    Pair = paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "_"),
    STRING_confidence = weight
  )

cat(sprintf("  PPI edges with confidence > 0.5: %d\n\n", nrow(ppi_edges)))

# ================================================================================
# 5. Filter by STRING confidence (> 0.5)
# ================================================================================
cat("Step 5: Filtering by STRING confidence > 0.5...\n")

cat(sprintf("  PPI edges already filtered: %d\n\n", nrow(ppi_edges)))

# ================================================================================
# 6. Combine correlation and PPI data
# ================================================================================
cat("Step 6: Combining correlation and PPI data...\n")

# Merge with correlation data (keep only one set of Gene1/Gene2)
network_edges <- both_sig %>%
  inner_join(ppi_edges %>% select(Pair, STRING_confidence), 
             by = "Pair")

cat(sprintf("  Edges validated by STRING: %d\n", nrow(network_edges)))

# Calculate correlation difference (EOD - CN)
network_edges <- network_edges %>%
  mutate(r_diff = EOD_r - CN_r)

# ================================================================================
# 7. Filter proteins with >= 3 connections
# ================================================================================
cat("Step 7: Filtering proteins with >= 3 connections...\n")

# Count connections per protein
protein_connections <- bind_rows(
  network_edges %>% select(Gene = Gene1),
  network_edges %>% select(Gene = Gene2)
) %>%
  count(Gene, name = "n_connections")

# Keep proteins with >= 3 connections
valid_proteins <- protein_connections %>%
  filter(n_connections >= 3) %>%
  pull(Gene)

cat(sprintf("  Proteins with >= 3 connections: %d\n", length(valid_proteins)))

# Filter edges
network_edges_filtered <- network_edges %>%
  filter(Gene1 %in% valid_proteins & Gene2 %in% valid_proteins)

cat(sprintf("  Final edges: %d\n\n", nrow(network_edges_filtered)))

if (nrow(network_edges_filtered) == 0) {
  cat("ERROR: No edges remain after filtering\n")
  cat("Exiting...\n")
  quit(status = 1)
}

# ================================================================================
# 8. Build igraph network
# ================================================================================
cat("Step 8: Building network graph...\n")

# Create igraph object
g <- graph_from_data_frame(
  network_edges_filtered %>% select(Gene1, Gene2, STRING_confidence, r_diff, EOD_r, CN_r),
  directed = FALSE
)

# Gene symbols are already in vertex names
V(g)$symbol <- V(g)$name

# Add module information - use only M number (M1, M2, M3, etc.)
gene_to_module <- hub_genes_filtered %>%
  select(Gene_Symbol, Module_Name) %>%
  distinct() %>%
  deframe()

V(g)$module <- gene_to_module[V(g)$symbol]

# Extract M number only (M1, M2, M3, etc.)
V(g)$module_label <- str_extract(V(g)$module, "^M\\d+")

# Extract color name for color mapping
V(g)$module_color <- str_extract(V(g)$module, "[^_]+$")

cat(sprintf("  Network nodes: %d\n", vcount(g)))
cat(sprintf("  Network edges: %d\n\n", ecount(g)))

# ================================================================================
# 9. Prepare visualization
# ================================================================================
cat("Step 9: Preparing visualization...\n")

# Define WGCNA module colors (exact mapping)
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

# Get unique module labels in network (sorted numerically)
modules_in_network <- unique(V(g)$module_label)
modules_in_network <- modules_in_network[!is.na(modules_in_network)]
# Sort by numeric value (M1, M2, M3... not M1, M10, M11...)
modules_in_network <- modules_in_network[order(as.numeric(str_extract(modules_in_network, "\\d+")))]

cat(sprintf("  Modules in network: %d\n", length(modules_in_network)))
cat(sprintf("    %s\n", paste(modules_in_network, collapse = ", ")))

# ================================================================================
# 10. Create publication-quality network plot
# ================================================================================
cat("\nStep 10: Creating network visualization...\n")

# Set seed for reproducible layout
set.seed(42)

# Create plot
p <- ggraph(g, layout = "stress") +
  # Edges: color by r_diff (EOD - CN), width by STRING confidence
  geom_edge_link(
    aes(
      width = STRING_confidence,
      color = r_diff
    ),
    alpha = 0.9
  ) +
  scale_edge_width_continuous(
    range = c(0.5, 3),
    breaks = c(0.5, 0.6, 0.7, 0.8, 0.9),
    labels = c("0.5", "0.6", "0.7", "0.8", "0.9"),
    name = "STRING\nconfidence"
  ) +
  scale_edge_colour_gradient2(
    low = "#0571B0",
    mid = "#F7F7F7", 
    high = "#CA0020",
    midpoint = 0,
    limits = c(-0.5, 0.5),
    breaks = seq(-0.5, 0.5, 0.25),
    name = "Correlation\ndifference\n(EOD - CN)",
    guide = guide_edge_colourbar(
      barwidth = 1.5,
      barheight = 12,
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  # Nodes: color by module color
  geom_node_point(
    aes(colour = module_color),
    size = 10,
    alpha = 0.95
  ) +
  scale_colour_manual(
    values = wgcna_colors,
    name = "WGCNA\nModule",
    breaks = {
      # Get unique colors and their corresponding labels
      color_label_df <- data.frame(
        color = V(g)$module_color,
        label = V(g)$module_label
      ) %>%
        distinct() %>%
        filter(!is.na(color) & !is.na(label)) %>%
        mutate(num = as.numeric(str_extract(label, "\\d+"))) %>%
        arrange(num)
      color_label_df$color
    },
    labels = {
      # Get unique colors and their corresponding labels
      color_label_df <- data.frame(
        color = V(g)$module_color,
        label = V(g)$module_label
      ) %>%
        distinct() %>%
        filter(!is.na(color) & !is.na(label)) %>%
        mutate(num = as.numeric(str_extract(label, "\\d+"))) %>%
        arrange(num)
      color_label_df$label
    },
    na.value = "grey50",
    guide = guide_legend(override.aes = list(size = 7))
  ) +
  # Node labels - black text with white shadow
  geom_shadowtext(
    aes(x = x, y = y, label = symbol),
    size = 3.8,
    fontface = "bold",
    color = "black",
    bg.colour = "white",
    bg.r = 0.15
  ) +
  # Theme
  theme_void() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing = unit(0.5, "cm"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(1, "cm"),
    plot.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Save plot to correlation_meta folder
output_dir <- "correlation_meta"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_file <- file.path(output_dir, "08_PPI_network_hub_genes.pdf")
ggsave(output_file, p, width = 16, height = 12, dpi = 300)

cat(sprintf("\n✓ Network plot saved: %s\n", output_file))

# Also save as PNG
output_png <- file.path(output_dir, "08_PPI_network_hub_genes.png")
ggsave(output_png, p, width = 16, height = 12, dpi = 300)

cat(sprintf("✓ Network plot saved: %s\n", output_png))

# ================================================================================
# 11. Save network data
# ================================================================================
cat("\nStep 11: Saving network data...\n")

# Save edge list with all information
edge_data <- network_edges_filtered %>%
  left_join(
    hub_genes_filtered %>% select(Gene_Symbol, Module1 = Module_Name),
    by = c("Gene1" = "Gene_Symbol")
  ) %>%
  left_join(
    hub_genes_filtered %>% select(Gene_Symbol, Module2 = Module_Name),
    by = c("Gene2" = "Gene_Symbol")
  ) %>%
  select(Gene1, Gene2, Module1, Module2,
         EOD_r, EOD_FDR, EOD_N,
         CN_r, CN_FDR, CN_N,
         r_diff, STRING_confidence)

write_csv(edge_data, file.path(output_dir, "08_PPI_network_edges.csv"))
cat(sprintf("✓ Edge data saved: %s\n", file.path(output_dir, "08_PPI_network_edges.csv")))

# Save node list
node_data <- tibble(
  Gene = V(g)$symbol,
  Module = V(g)$module,
  Module_Label = V(g)$module_label,
  Module_Color = V(g)$module_color,
  Degree = degree(g)
) %>%
  left_join(
    hub_genes_filtered %>% select(Gene_Symbol, kME),
    by = c("Gene" = "Gene_Symbol")
  ) %>%
  arrange(desc(Degree))

write_csv(node_data, file.path(output_dir, "08_PPI_network_nodes.csv"))
cat(sprintf("✓ Node data saved: %s\n", file.path(output_dir, "08_PPI_network_nodes.csv")))

# ================================================================================
# 12. Summary statistics
# ================================================================================
cat("\n================================================================================\n")
cat("NETWORK SUMMARY\n")
cat("================================================================================\n\n")

cat(sprintf("Network Statistics:\n"))
cat(sprintf("  Nodes: %d\n", vcount(g)))
cat(sprintf("  Edges: %d\n", ecount(g)))
cat(sprintf("  Density: %.4f\n", edge_density(g)))
cat(sprintf("  Average degree: %.2f\n", mean(degree(g))))
cat(sprintf("  Clustering coefficient: %.4f\n", transitivity(g, type = "global")))

cat(sprintf("\nModules represented:\n"))
module_counts <- node_data %>%
  count(Module, Module_Label, Module_Color) %>%
  arrange(as.numeric(str_extract(Module_Label, "\\d+")))

for (i in 1:nrow(module_counts)) {
  cat(sprintf("  %s (%s): %d proteins\n", 
              module_counts$Module_Label[i],
              module_counts$Module[i],
              module_counts$n[i]))
}

cat(sprintf("\nTop 10 hub proteins (by degree):\n"))
top_hubs <- node_data %>%
  slice_head(n = 10)

for (i in 1:nrow(top_hubs)) {
  cat(sprintf("  %d. %s (%s): degree=%d, kME=%.3f\n",
              i,
              top_hubs$Gene[i],
              top_hubs$Module[i],
              top_hubs$Degree[i],
              top_hubs$kME[i]))
}

cat(sprintf("\nCorrelation statistics:\n"))
cat(sprintf("  EOD r range: [%.3f, %.3f]\n", 
            min(edge_data$EOD_r), max(edge_data$EOD_r)))
cat(sprintf("  CN r range: [%.3f, %.3f]\n",
            min(edge_data$CN_r), max(edge_data$CN_r)))
cat(sprintf("  r_diff range: [%.3f, %.3f]\n",
            min(edge_data$r_diff), max(edge_data$r_diff)))
cat(sprintf("  STRING confidence range: [%.3f, %.3f]\n",
            min(edge_data$STRING_confidence), max(edge_data$STRING_confidence)))

cat("\n================================================================================\n")
cat("✓ PPI NETWORK ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")

cat("Output files (in correlation_meta/):\n")
cat("  - 08_PPI_network_hub_genes.pdf: Network visualization (PDF)\n")
cat("  - 08_PPI_network_hub_genes.png: Network visualization (PNG)\n")
cat("  - 08_PPI_network_edges.csv: Edge list with correlations\n")
cat("  - 08_PPI_network_nodes.csv: Node list with module info\n")
cat("\n")
