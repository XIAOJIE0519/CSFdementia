"""
EOD Single-Cell Integration Pipeline - Visualization Module
============================================================

Generate core visualizations for single-cell analysis.

This module creates:
1. QC plots (quality control metrics)
2. PCA plots (principal component analysis)
3. UMAP plots (uniform manifold approximation and projection)
4. Cell type composition plots

Scientific Standards:
- Publication-quality figures
- Comprehensive quality control visualization
- Clear and informative plots

Author: Bioinformatics Analysis
Date: 2026-02-14
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path
import scanpy as sc

# Set style
sns.set_style('whitegrid')
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

print("\n" + "="*80)
print("EOD Pipeline - Visualization Module")
print("="*80)

# ============================================================================
# QC Visualization
# ============================================================================

def plot_qc_metrics(adata, output_dir, prefix="qc"):
    """Plot QC metrics"""
    print("\nGenerating QC plots...")
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Quality Control Metrics', fontsize=16, fontweight='bold')
    
    # Plot 1: n_genes distribution
    axes[0, 0].hist(adata.obs['n_genes_by_counts'], bins=100, edgecolor='black')
    axes[0, 0].set_xlabel('Number of genes')
    axes[0, 0].set_ylabel('Number of cells')
    axes[0, 0].set_title('Genes per cell')
    
    # Plot 2: total_counts distribution
    axes[0, 1].hist(adata.obs['total_counts'], bins=100, edgecolor='black')
    axes[0, 1].set_xlabel('Total counts')
    axes[0, 1].set_ylabel('Number of cells')
    axes[0, 1].set_title('UMI counts per cell')
    axes[0, 1].set_xscale('log')
    
    # Plot 3: pct_counts_mt distribution
    axes[0, 2].hist(adata.obs['pct_counts_mt'], bins=100, edgecolor='black')
    axes[0, 2].set_xlabel('Mitochondrial %')
    axes[0, 2].set_ylabel('Number of cells')
    axes[0, 2].set_title('Mitochondrial percentage')
    
    # Plot 4: Scatter n_genes vs total_counts
    axes[1, 0].scatter(adata.obs['total_counts'], adata.obs['n_genes_by_counts'], 
                      s=1, alpha=0.5)
    axes[1, 0].set_xlabel('Total counts')
    axes[1, 0].set_ylabel('Number of genes')
    axes[1, 0].set_title('Genes vs UMI counts')
    axes[1, 0].set_xscale('log')
    
    # Plot 5: Scatter n_genes vs pct_counts_mt
    axes[1, 1].scatter(adata.obs['n_genes_by_counts'], adata.obs['pct_counts_mt'],
                      s=1, alpha=0.5)
    axes[1, 1].set_xlabel('Number of genes')
    axes[1, 1].set_ylabel('Mitochondrial %')
    axes[1, 1].set_title('Genes vs Mitochondrial %')
    
    # Plot 6: Cells per sample
    sample_counts = adata.obs['sample_id'].value_counts().head(20)
    axes[1, 2].barh(range(len(sample_counts)), sample_counts.values)
    axes[1, 2].set_yticks(range(len(sample_counts)))
    axes[1, 2].set_yticklabels(sample_counts.index, fontsize=8)
    axes[1, 2].set_xlabel('Number of cells')
    axes[1, 2].set_title('Cells per sample (top 20)')
    
    plt.tight_layout()
    output_file = Path(output_dir) / f"{prefix}_metrics.png"
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_file}")

# ============================================================================
# PCA Visualization
# ============================================================================

def plot_pca(adata, output_dir, color_by=['subtype', 'data_source', 'cell_type']):
    """Plot PCA colored by different variables"""
    print("\nGenerating PCA plots...")
    
    for var in color_by:
        if var not in adata.obs.columns:
            continue
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        unique_vals = adata.obs[var].unique()
        colors = plt.cm.tab20(np.linspace(0, 1, len(unique_vals)))
        
        for i, val in enumerate(unique_vals):
            mask = adata.obs[var] == val
            ax.scatter(adata.obsm['X_pca'][mask, 0], 
                      adata.obsm['X_pca'][mask, 1],
                      s=5, alpha=0.6, label=val, c=[colors[i]])
        
        ax.set_xlabel(f'PC1 ({adata.uns["pca"]["variance_ratio"][0]:.1%})')
        ax.set_ylabel(f'PC2 ({adata.uns["pca"]["variance_ratio"][1]:.1%})')
        ax.set_title(f'PCA colored by {var}')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        
        plt.tight_layout()
        output_file = Path(output_dir) / f"pca_by_{var}.png"
        plt.savefig(output_file, bbox_inches='tight')
        plt.close()
        
        print(f"  Saved: {output_file}")
    
    # Variance explained plot
    fig, ax = plt.subplots(figsize=(10, 6))
    variance_ratio = adata.uns['pca']['variance_ratio'][:30]
    ax.bar(range(len(variance_ratio)), variance_ratio)
    ax.set_xlabel('Principal Component')
    ax.set_ylabel('Variance Explained')
    ax.set_title('PCA Variance Explained')
    
    plt.tight_layout()
    output_file = Path(output_dir) / "pca_variance.png"
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_file}")

# ============================================================================
# UMAP Visualization
# ============================================================================

def plot_umap(adata, output_dir):
    """Plot UMAP colored by cell type, disease group, and data source"""
    print("\nGenerating UMAP plots...")
    
    if 'X_umap' not in adata.obsm:
        print("  Computing UMAP coordinates...")
        
        # 检查是否有Harmony校正的PCA（批次校正后的结果）
        if 'X_pca_harmony' in adata.obsm:
            print("  Using Harmony-corrected PCA (batch-corrected)...")
            use_rep = 'X_pca_harmony'
        else:
            print("  Warning: Harmony-corrected PCA not found, using original PCA...")
            print("  This may result in batch effects in UMAP!")
            use_rep = 'X_pca'
        
        # 使用eod_01配置的参数：n_neighbors=15, n_pcs=30
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep=use_rep)
        sc.tl.umap(adata, min_dist=0.3)
    
    # Create disease_group column if not exists
    if 'disease_group' not in adata.obs.columns:
        print("  Creating disease_group column...")
        
        def get_disease_group(row):
            subtype = row['subtype']
            age_class = row['age_class'] if 'age_class' in row else 'HC'
            
            if subtype == 'HC':
                return 'CN'
            elif subtype == 'AD' and age_class == 'EOD':
                return 'EOAD'
            elif subtype == 'AD' and age_class == 'LOD':
                return 'LOAD'
            elif subtype == 'FTD' and age_class == 'EOD':
                return 'EOFTD'
            elif subtype == 'FTD' and age_class == 'LOD':
                return 'LOFTD'
            else:
                return 'Other'
        
        adata.obs['disease_group'] = adata.obs.apply(get_disease_group, axis=1)
    
    # Define colors
    cell_type_colors = {
        'ExcitatoryNeurons': '#1f77b4', 'InhibitoryNeurons': '#ff7f0e',
        'Astrocytes': '#2ca02c', 'Oligodendrocytes': '#d62728',
        'OPCs': '#9467bd', 'Microglia': '#8c564b',
        'Macrophages': '#e377c2', 'EndothelialCells': '#7f7f7f',
        'Pericytes': '#bcbd22', 'Fibroblasts': '#17becf'
    }
    
    disease_colors = {
        'CN': '#2ca02c', 'EOAD': '#d62728', 'LOAD': '#ff7f0e',
        'EOFTD': '#9467bd', 'LOFTD': '#8c564b'
    }
    
    # Plot UMAP by cell type
    fig, ax = plt.subplots(figsize=(10, 8))
    for ct in sorted(adata.obs['cell_type'].unique()):
        mask = adata.obs['cell_type'] == ct
        ax.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
                  c=cell_type_colors.get(ct, '#CCCCCC'), s=5, alpha=0.6, label=ct)
    ax.set_xlabel('UMAP 1', fontsize=12, fontweight='bold')
    ax.set_ylabel('UMAP 2', fontsize=12, fontweight='bold')
    ax.set_title('UMAP by Cell Type', fontsize=14, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    output_file = Path(output_dir) / "umap_by_celltype.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"  Saved: {output_file}")
    
    # Plot UMAP by disease group
    fig, ax = plt.subplots(figsize=(10, 8))
    for dg in ['CN', 'EOAD', 'LOAD', 'EOFTD', 'LOFTD']:
        mask = adata.obs['disease_group'] == dg
        if mask.sum() > 0:
            ax.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
                      c=disease_colors.get(dg, '#CCCCCC'), s=5, alpha=0.6, label=dg)
    ax.set_xlabel('UMAP 1', fontsize=12, fontweight='bold')
    ax.set_ylabel('UMAP 2', fontsize=12, fontweight='bold')
    ax.set_title('UMAP by Disease Group', fontsize=14, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    output_file = Path(output_dir) / "umap_by_disease.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"  Saved: {output_file}")
    
    # Plot UMAP by data source
    fig, ax = plt.subplots(figsize=(10, 8))
    data_sources = sorted(adata.obs['data_source'].unique())
    source_palette = plt.cm.Set3(np.linspace(0, 1, len(data_sources)))
    source_colors = dict(zip(data_sources, source_palette))
    for ds in data_sources:
        mask = adata.obs['data_source'] == ds
        ax.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
                  c=[source_colors[ds]], s=5, alpha=0.6, label=ds)
    ax.set_xlabel('UMAP 1', fontsize=12, fontweight='bold')
    ax.set_ylabel('UMAP 2', fontsize=12, fontweight='bold')
    ax.set_title('UMAP by Data Source', fontsize=14, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    output_file = Path(output_dir) / "umap_by_source.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"  Saved: {output_file}")

# ============================================================================
# Circular UMAP with Three-Ring Annotations (Python Version)
# ============================================================================

def plot_circular_umap_annotations(adata, output_dir):
    """
    Plot circular UMAP with three-ring annotations
    - Inner ring: Cell type
    - Middle ring: Disease group
    - Outer ring: Data source
    """
    print("\nGenerating circular UMAP with three-ring annotations...")
    
    import matplotlib.patches as mpatches
    from matplotlib.collections import LineCollection
    
    # Get UMAP coordinates
    umap_coords = adata.obsm['X_umap']
    
    # Define colors
    cell_type_colors = {
        'ExcitatoryNeurons': '#1f77b4', 'InhibitoryNeurons': '#ff7f0e',
        'Astrocytes': '#2ca02c', 'Oligodendrocytes': '#d62728',
        'OPCs': '#9467bd', 'Microglia': '#8c564b',
        'Macrophages': '#e377c2', 'EndothelialCells': '#7f7f7f',
        'Pericytes': '#bcbd22', 'Fibroblasts': '#17becf'
    }
    
    disease_colors = {
        'CN': '#2ca02c', 'EOAD': '#d62728', 'LOAD': '#ff7f0e',
        'EOFTD': '#9467bd', 'LOFTD': '#8c564b'
    }
    
    # Create figure with polar projection
    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111, projection='polar')
    
    # Convert UMAP coordinates to polar
    x, y = umap_coords[:, 0], umap_coords[:, 1]
    
    # Center coordinates
    x_centered = x - x.mean()
    y_centered = y - y.mean()
    
    # Convert to polar
    theta = np.arctan2(y_centered, x_centered)
    r = np.sqrt(x_centered**2 + y_centered**2)
    
    # Normalize radius
    r_norm = (r - r.min()) / (r.max() - r.min()) * 0.6 + 0.1
    
    # Plot cells colored by cell type
    for ct in sorted(adata.obs['cell_type'].unique()):
        mask = adata.obs['cell_type'] == ct
        ax.scatter(theta[mask], r_norm[mask], 
                  c=cell_type_colors.get(ct, '#CCCCCC'),
                  s=1, alpha=0.5, label=ct)
    
    # Add three annotation rings
    ring_positions = [0.75, 0.85, 0.95]
    ring_labels = ['Cell Type', 'Disease Group', 'Data Source']
    
    # Ring 1: Cell type (inner)
    cell_types = sorted(adata.obs['cell_type'].unique())
    n_ct = len(cell_types)
    for i, ct in enumerate(cell_types):
        theta_start = 2 * np.pi * i / n_ct
        theta_end = 2 * np.pi * (i + 1) / n_ct
        theta_mid = (theta_start + theta_end) / 2
        
        # Draw arc
        theta_arc = np.linspace(theta_start, theta_end, 100)
        r_arc = np.full_like(theta_arc, ring_positions[0])
        ax.plot(theta_arc, r_arc, color=cell_type_colors.get(ct, '#CCCCCC'), linewidth=8)
    
    # Ring 2: Disease group (middle)
    disease_groups = ['CN', 'EOAD', 'LOAD', 'EOFTD', 'LOFTD']
    n_dg = len(disease_groups)
    for i, dg in enumerate(disease_groups):
        theta_start = 2 * np.pi * i / n_dg
        theta_end = 2 * np.pi * (i + 1) / n_dg
        
        # Draw arc
        theta_arc = np.linspace(theta_start, theta_end, 100)
        r_arc = np.full_like(theta_arc, ring_positions[1])
        ax.plot(theta_arc, r_arc, color=disease_colors.get(dg, '#CCCCCC'), linewidth=8)
    
    # Ring 3: Data source (outer)
    data_sources = sorted(adata.obs['data_source'].unique())
    n_ds = len(data_sources)
    source_palette = plt.cm.Set3(np.linspace(0, 1, n_ds))
    for i, ds in enumerate(data_sources):
        theta_start = 2 * np.pi * i / n_ds
        theta_end = 2 * np.pi * (i + 1) / n_ds
        
        # Draw arc
        theta_arc = np.linspace(theta_start, theta_end, 100)
        r_arc = np.full_like(theta_arc, ring_positions[2])
        ax.plot(theta_arc, r_arc, color=source_palette[i], linewidth=8)
    
    # Remove grid and labels
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines['polar'].set_visible(False)
    
    # Add title
    ax.set_title('UMAP with Three-Ring Annotations\n(Inner: Cell Type | Middle: Disease Group | Outer: Data Source)',
                fontsize=16, fontweight='bold', pad=20)
    
    # Add legend for cell types
    legend_elements = [mpatches.Patch(facecolor=cell_type_colors.get(ct, '#CCCCCC'), 
                                     label=ct) for ct in cell_types]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.1, 1.0),
             fontsize=10, title='Cell Type', title_fontsize=12)
    
    plt.tight_layout()
    output_file = Path(output_dir) / "umap_circular_annotations.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"  Saved: {output_file}")

# ============================================================================
# Cell Type Composition
# ============================================================================

def plot_cell_type_composition(adata, output_dir):
    """Plot cell type composition by disease group"""
    print("\nGenerating cell type composition plots...")
    
    # Overall cell type distribution
    fig, ax = plt.subplots(figsize=(12, 6))
    celltype_counts = adata.obs['cell_type'].value_counts()
    colors = plt.cm.tab10(np.linspace(0, 1, len(celltype_counts)))
    ax.barh(range(len(celltype_counts)), celltype_counts.values, color=colors)
    ax.set_yticks(range(len(celltype_counts)))
    ax.set_yticklabels(celltype_counts.index, fontsize=10)
    ax.set_xlabel('Number of Cells', fontsize=12, fontweight='bold')
    ax.set_title('Cell Type Distribution', fontsize=14, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    output_file = Path(output_dir) / "celltype_distribution.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"  Saved: {output_file}")
    
    # Cell type by disease group
    if 'disease_group' in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(14, 8))
        
        ct_disease = pd.crosstab(adata.obs['cell_type'], adata.obs['disease_group'], normalize='columns')
        ct_disease = ct_disease[['CN', 'EOAD', 'LOAD', 'EOFTD', 'LOFTD']]
        ct_disease.plot(kind='bar', stacked=True, ax=ax, colormap='Set3')
        
        ax.set_xlabel('Cell Type', fontsize=12, fontweight='bold')
        ax.set_ylabel('Proportion', fontsize=12, fontweight='bold')
        ax.set_title('Cell Type Composition by Disease Group', fontsize=14, fontweight='bold')
        ax.legend(title='Disease Group', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.xticks(rotation=45, ha='right')
        
        plt.tight_layout()
        output_file = Path(output_dir) / "celltype_by_disease.png"
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        plt.close()
        
        print(f"  Saved: {output_file}")


# ============================================================================
# Main Visualization Function
# ============================================================================

def generate_all_visualizations(adata, output_dir):
    """Generate core visualizations"""
    print("="*80)
    print("Generating Core Visualizations")
    print("="*80)
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # QC plots
    if 'n_genes_by_counts' in adata.obs.columns:
        plot_qc_metrics(adata, output_dir, prefix="final_qc")
    
    # PCA plots
    if 'X_pca' in adata.obsm.keys():
        plot_pca(adata, output_dir)
    
    # UMAP plots (compute if not exists)
    plot_umap(adata, output_dir)
    
    # Circular UMAP with three-ring annotations (Python version)
    if 'X_umap' in adata.obsm.keys() and 'cell_type' in adata.obs.columns:
        plot_circular_umap_annotations(adata, output_dir)
    
    # Cell type composition plots
    if 'cell_type' in adata.obs.columns:
        plot_cell_type_composition(adata, output_dir)
    
    print("\n" + "="*80)
    print("Core Visualizations Complete")
    print("="*80)

# ============================================================================
# Execute if run directly
# ============================================================================

if __name__ == '__main__':
    import sys
    from pathlib import Path
    
    print("\n" + "="*80)
    print("EOD Visualization - Core Pipeline")
    print("="*80)
    
    # Default paths
    default_input = Path("./EOD/eod_annotated.h5ad")
    default_output = Path("./EOD/figures")
    
    # Check input file
    if not default_input.exists():
        print(f"\n✗ ERROR: Input file not found: {default_input}")
        sys.exit(1)
    
    # Create output directory
    default_output.mkdir(exist_ok=True)
    
    # Load data
    print(f"\nLoading data from: {default_input}")
    adata = sc.read_h5ad(default_input)
    print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    
    # Generate core visualizations (QC, PCA, UMAP, cell type composition)
    generate_all_visualizations(
        adata=adata,
        output_dir=default_output
    )
    
    # Don't save h5ad (too large), only save necessary CSV files
    print("\nSkipping h5ad save (file too large, CSV exports are sufficient)")
    
    # Export UMAP coordinates and metadata to CSV for R scripts
    print("\nExporting UMAP coordinates and metadata to CSV...")
    umap_df = pd.DataFrame(
        adata.obsm['X_umap'],
        columns=['UMAP_1', 'UMAP_2'],
        index=adata.obs_names
    )
    
    # Combine UMAP coordinates with metadata
    metadata_cols = ['cell_type', 'disease_group', 'data_source', 'sample_id']
    export_df = pd.concat([umap_df, adata.obs[metadata_cols]], axis=1)
    
    # Save to CSV
    csv_file = Path("./EOD/umap_coordinates_metadata.csv")
    export_df.to_csv(csv_file, index=True)
    print(f"  Saved: {csv_file}")
    print(f"  Exported: {len(export_df):,} cells with UMAP coordinates and metadata")
    
    # Also export Harmony-corrected PCA if available (for advanced R analysis)
    if 'X_pca_harmony' in adata.obsm:
        print("\nExporting Harmony-corrected PCA coordinates...")
        pca_harmony_df = pd.DataFrame(
            adata.obsm['X_pca_harmony'][:, :30],  # First 30 PCs
            columns=[f'PC{i+1}' for i in range(30)],
            index=adata.obs_names
        )
        pca_file = Path("./EOD/pca_harmony_coordinates.csv")
        pca_harmony_df.to_csv(pca_file, index=True)
        print(f"  Saved: {pca_file}")
        print(f"  Exported: 30 Harmony-corrected PCs for {len(pca_harmony_df):,} cells")
    
    # Calculate and export module expression data for Mantel test
    print("\n" + "="*80)
    print("Calculating Module Expression Data for Mantel Test")
    print("="*80)
    
    module_file = Path("./module_protein_gene_summary.csv")
    if module_file.exists():
        print(f"\nLoading module information from: {module_file}")
        module_df = pd.read_csv(module_file)
        print(f"  Loaded {len(module_df)} protein-module mappings")
        
        # Get unique modules
        modules = sorted(module_df['Module'].unique())
        print(f"  Found {len(modules)} modules")
        
        # Create output directory for module data
        module_output_dir = Path("./EOD/figures/module_heatmaps")
        module_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Calculate module expression by cell type × disease group
        print(f"\nCalculating module expression by cell type and disease group...")
        
        disease_groups = ['CN', 'EOAD', 'LOAD', 'EOFTD', 'LOFTD']
        cell_types = sorted(adata.obs['cell_type'].unique())
        
        # Create index: CellType_DiseaseGroup
        sample_names = []
        for ct in cell_types:
            for dg in disease_groups:
                sample_names.append(f"{ct}_{dg}")
        
        module_expr = pd.DataFrame(index=sample_names, columns=modules, dtype=float)
        
        # Calculate mean expression for each module in each cell type × disease group
        for module in modules:
            # Get genes in this module
            module_genes = module_df[module_df['Module'] == module]['Protein'].tolist()
            
            # Filter to genes present in data
            available_genes = [g for g in module_genes if g in adata.var_names]
            
            if len(available_genes) == 0:
                print(f"  Warning: No genes found for {module}")
                module_expr[module] = 0.0
                continue
            
            # Calculate mean expression for each cell type × disease group
            for ct in cell_types:
                for dg in disease_groups:
                    sample_name = f"{ct}_{dg}"
                    
                    # Get cells matching both cell type and disease group
                    mask = (adata.obs['cell_type'] == ct) & (adata.obs['disease_group'] == dg)
                    n_cells = mask.sum()
                    
                    if n_cells == 0:
                        module_expr.loc[sample_name, module] = 0.0
                        continue
                    
                    # Get expression matrix for these genes and cells
                    gene_indices = [i for i, g in enumerate(adata.var_names) if g in available_genes]
                    cell_indices = np.where(mask)[0]
                    
                    # Extract expression data
                    expr_data = adata.X[np.ix_(cell_indices, gene_indices)]
                    
                    # Convert to dense if sparse
                    if hasattr(expr_data, 'toarray'):
                        expr_data = expr_data.toarray()
                    
                    # Calculate mean expression across genes and cells
                    mean_expr = np.mean(expr_data)
                    module_expr.loc[sample_name, module] = mean_expr
            
            print(f"  Processed: {module}")
        
        # Save module expression data (cell type × disease group)
        output_file = module_output_dir / "module_expression_by_celltype_disease.csv"
        module_expr.to_csv(output_file)
        print(f"\n  Saved: {output_file}")
        print(f"  Samples: {len(sample_names)} (cell types × disease groups), Modules: {len(modules)}")
        
        # Calculate normalized expression (within each cell type × disease group)
        print(f"\nCalculating normalized module expression...")
        module_expr_normalized = module_expr.copy()
        
        # Z-score normalization within each cell type × disease group
        for ct in cell_types:
            for dg in disease_groups:
                sample_name = f"{ct}_{dg}"
                
                # Get expression values for this sample
                expr_values = module_expr.loc[sample_name, :].values
                
                # Z-score normalization
                mean_val = np.mean(expr_values)
                std_val = np.std(expr_values)
                
                if std_val > 0:
                    module_expr_normalized.loc[sample_name, :] = (expr_values - mean_val) / std_val
                else:
                    module_expr_normalized.loc[sample_name, :] = 0.0
        
        # Save normalized data
        output_file_norm = module_output_dir / "module_expression_normalized.csv"
        module_expr_normalized.to_csv(output_file_norm)
        print(f"  Saved: {output_file_norm}")
        print(f"  Normalized within each cell type × disease group combination")
        
        # Also save overall disease group summary (for reference)
        disease_summary = pd.DataFrame(index=disease_groups, columns=modules, dtype=float)
        for dg in disease_groups:
            for module in modules:
                # Average across all cell types for this disease group
                dg_samples = [s for s in sample_names if s.endswith(f"_{dg}")]
                disease_summary.loc[dg, module] = module_expr.loc[dg_samples, module].mean()
        
        output_file2 = module_output_dir / "module_expression_by_disease.csv"
        disease_summary.to_csv(output_file2)
        print(f"  Saved: {output_file2}")
        print(f"  Disease groups: {len(disease_groups)}, Modules: {len(modules)}")
        
        # Calculate module expression by real sample (for correlation analysis)
        print(f"\nCalculating module expression by real sample...")
        
        # Get unique samples
        unique_samples = sorted(adata.obs['sample_id'].unique())
        print(f"  Found {len(unique_samples)} unique samples")
        
        module_expr_by_sample = pd.DataFrame(index=unique_samples, columns=modules, dtype=float)
        
        # Calculate mean expression for each module in each sample
        for module in modules:
            # Get genes in this module
            module_genes = module_df[module_df['Module'] == module]['Protein'].tolist()
            
            # Filter to genes present in data
            available_genes = [g for g in module_genes if g in adata.var_names]
            
            if len(available_genes) == 0:
                module_expr_by_sample[module] = 0.0
                continue
            
            # Calculate mean expression for each sample
            for sample in unique_samples:
                # Get cells from this sample
                mask = adata.obs['sample_id'] == sample
                n_cells = mask.sum()
                
                if n_cells == 0:
                    module_expr_by_sample.loc[sample, module] = 0.0
                    continue
                
                # Get expression matrix for these genes and cells
                gene_indices = [i for i, g in enumerate(adata.var_names) if g in available_genes]
                cell_indices = np.where(mask)[0]
                
                # Extract expression data
                expr_data = adata.X[np.ix_(cell_indices, gene_indices)]
                
                # Convert to dense if sparse
                if hasattr(expr_data, 'toarray'):
                    expr_data = expr_data.toarray()
                
                # Calculate mean expression across genes and cells
                mean_expr = np.mean(expr_data)
                module_expr_by_sample.loc[sample, module] = mean_expr
        
        # Add sample metadata
        sample_metadata = adata.obs.groupby('sample_id').first()[['disease_group', 'data_source', 'subtype']]
        module_expr_by_sample_with_meta = pd.concat([sample_metadata, module_expr_by_sample], axis=1)
        
        # Save sample-level data
        output_file_sample = module_output_dir / "module_expression_by_sample.csv"
        module_expr_by_sample_with_meta.to_csv(output_file_sample)
        print(f"  Saved: {output_file_sample}")
        print(f"  Samples: {len(unique_samples)}, Modules: {len(modules)}")
        
        print("\n✓ Module expression data calculation complete!")
    else:
        print(f"\n⚠ Warning: Module file not found: {module_file}")
        print("  Skipping module expression calculation")
    
    print("\n" + "="*80)
    print("All processing complete!")
    print("="*80)
    print(f"\nOutput directory: {default_output}/")
    print("\nGenerated files:")
    print("  - final_qc_metrics.png (Quality control)")
    print("  - pca_by_*.png (PCA plots)")
    print("  - pca_variance.png (PCA variance)")
    print("  - umap_by_celltype.png (UMAP by cell type)")
    print("  - umap_by_disease.png (UMAP by disease group)")
    print("  - umap_by_source.png (UMAP by data source)")
    print("  - celltype_distribution.png (Cell type counts)")
    print("  - celltype_by_disease.png (Cell type composition)")
    
    print("\n" + "="*80)
    print("Done!")
    print("="*80)
