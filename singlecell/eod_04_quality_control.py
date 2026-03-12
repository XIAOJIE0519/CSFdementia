"""
EOD Single-Cell Integration Pipeline - Quality Control Module
==============================================================

Perform rigorous quality control following Nature-level standards.

This module:
1. Calculate QC metrics (mitochondrial %, ribosomal %, etc.)
2. Filter low-quality cells
3. Filter low-expression genes
4. Normalize and log-transform data
5. Select highly variable genes

Scientific Standards:
- Stringent QC thresholds
- Document all filtering steps
- Preserve raw counts
- Maintain reproducibility

Author: Bioinformatics Analysis
Date: 2026-02-14
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse
from datetime import datetime
import gc

print("\n" + "="*80)
print("EOD Pipeline - Step 4: Quality Control")
print("="*80)
print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# ============================================================================
# QC Metric Calculation
# ============================================================================

def calculate_qc_metrics(adata):
    """
    Calculate quality control metrics
    
    Parameters:
        adata: AnnData object
    
    Returns:
        adata: AnnData with QC metrics
    """
    
    print(f"\n{'='*80}")
    print("Calculating QC Metrics")
    print(f"{'='*80}")
    
    # Identify mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    n_mt_genes = adata.var['mt'].sum()
    print(f"  Mitochondrial genes: {n_mt_genes}")
    
    # Identify ribosomal genes
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    n_ribo_genes = adata.var['ribo'].sum()
    print(f"  Ribosomal genes: {n_ribo_genes}")
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt', 'ribo'],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    
    print(f"\nQC metrics calculated:")
    print(f"  n_genes_by_counts: median = {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"  total_counts: median = {adata.obs['total_counts'].median():.0f}")
    print(f"  pct_counts_mt: median = {adata.obs['pct_counts_mt'].median():.2f}%")
    print(f"  pct_counts_ribo: median = {adata.obs['pct_counts_ribo'].median():.2f}%")
    
    return adata

# ============================================================================
# Cell Filtering
# ============================================================================

def filter_cells_by_qc(adata, min_genes=200, max_genes=10000, max_pct_mito=20, 
                       min_counts=500, max_counts=100000):
    """
    Filter cells based on QC metrics
    
    Parameters:
        adata: AnnData object
        min_genes: Minimum genes per cell
        max_genes: Maximum genes per cell (doublet filter)
        max_pct_mito: Maximum mitochondrial percentage
        min_counts: Minimum UMI counts
        max_counts: Maximum UMI counts
    
    Returns:
        adata: Filtered AnnData object
    """
    
    print(f"\n{'='*80}")
    print("Cell Quality Control Filtering")
    print(f"{'='*80}")
    
    print(f"Before filtering: {adata.n_obs:,} cells")
    
    print(f"\nFiltering criteria:")
    print(f"  min_genes: {min_genes}")
    print(f"  max_genes: {max_genes}")
    print(f"  max_pct_mito: {max_pct_mito}%")
    print(f"  min_counts: {min_counts}")
    print(f"  max_counts: {max_counts}")
    
    # Apply filters
    n_before = adata.n_obs
    
    # Filter by gene count
    adata = adata[adata.obs['n_genes_by_counts'] >= min_genes, :].copy()
    n_after_min_genes = adata.n_obs
    print(f"\nAfter min_genes filter: {n_after_min_genes:,} cells ({n_before - n_after_min_genes:,} removed)")
    
    adata = adata[adata.obs['n_genes_by_counts'] <= max_genes, :].copy()
    n_after_max_genes = adata.n_obs
    print(f"After max_genes filter: {n_after_max_genes:,} cells ({n_after_min_genes - n_after_max_genes:,} removed)")
    
    # Filter by mitochondrial percentage
    adata = adata[adata.obs['pct_counts_mt'] <= max_pct_mito, :].copy()
    n_after_mito = adata.n_obs
    print(f"After mito filter: {n_after_mito:,} cells ({n_after_max_genes - n_after_mito:,} removed)")
    
    # Filter by total counts
    adata = adata[adata.obs['total_counts'] >= min_counts, :].copy()
    n_after_min_counts = adata.n_obs
    print(f"After min_counts filter: {n_after_min_counts:,} cells ({n_after_mito - n_after_min_counts:,} removed)")
    
    adata = adata[adata.obs['total_counts'] <= max_counts, :].copy()
    n_after_max_counts = adata.n_obs
    print(f"After max_counts filter: {n_after_max_counts:,} cells ({n_after_min_counts - n_after_max_counts:,} removed)")
    
    print(f"\nFinal: {adata.n_obs:,} cells ({n_before - adata.n_obs:,} removed, {adata.n_obs/n_before*100:.1f}% retained)")
    
    return adata

# ============================================================================
# Gene Filtering
# ============================================================================

def filter_genes_by_expression(adata, min_cells=3):
    """
    Filter genes by minimum cell expression
    
    Parameters:
        adata: AnnData object
        min_cells: Minimum cells expressing gene
    
    Returns:
        adata: Filtered AnnData object
    """
    
    print(f"\n{'='*80}")
    print("Gene Expression Filtering")
    print(f"{'='*80}")
    
    print(f"Before filtering: {adata.n_vars:,} genes")
    print(f"Filtering criterion: expressed in ≥{min_cells} cells")
    
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    print(f"After filtering: {adata.n_vars:,} genes")
    print(f"Removed: {adata.n_vars:,} genes")
    
    return adata

# ============================================================================
# Normalization
# ============================================================================

def normalize_and_log_transform(adata, target_sum=1e4):
    """
    Normalize and log-transform data (manual implementation to avoid numba)
    
    Parameters:
        adata: AnnData object
        target_sum: Target sum for normalization
    
    Returns:
        adata: Normalized AnnData object
    """
    
    print(f"\n{'='*80}")
    print("Normalization and Log Transformation")
    print(f"{'='*80}")
    
    # Save raw counts
    print("  Saving raw counts to layer...")
    if issparse(adata.X):
        adata.layers['counts'] = adata.X.copy()
    else:
        from scipy.sparse import csr_matrix
        adata.layers['counts'] = csr_matrix(adata.X)
    
    # Manual normalization (avoid numba bug)
    print(f"  Normalizing to target sum = {target_sum}...")
    
    if issparse(adata.X):
        counts_per_cell = np.array(adata.X.sum(axis=1)).flatten()
    else:
        counts_per_cell = adata.X.sum(axis=1)
    
    # Avoid division by zero
    counts_per_cell[counts_per_cell == 0] = 1
    
    # Normalize
    if issparse(adata.X):
        from scipy.sparse import diags
        scaling_factors = target_sum / counts_per_cell
        adata.X = diags(scaling_factors) @ adata.X
    else:
        adata.X = adata.X / counts_per_cell[:, None] * target_sum
    
    # Log transform
    print("  Log1p transformation...")
    if issparse(adata.X):
        adata.X.data = np.log1p(adata.X.data)
    else:
        adata.X = np.log1p(adata.X)
    
    print("  Normalization complete")
    
    del counts_per_cell
    gc.collect()
    
    return adata

# ============================================================================
# Highly Variable Genes
# ============================================================================

def select_highly_variable_genes(adata, n_top_genes=2000):
    """
    Select highly variable genes (manual implementation)
    
    Parameters:
        adata: AnnData object
        n_top_genes: Number of top genes to select
    
    Returns:
        adata: AnnData with HVG annotation
    """
    
    print(f"\n{'='*80}")
    print("Selecting Highly Variable Genes")
    print(f"{'='*80}")
    
    print(f"  Target: {n_top_genes} genes")
    print(f"  Computing gene statistics...")
    
    # Calculate gene statistics
    if issparse(adata.X):
        gene_means = np.array(adata.X.mean(axis=0)).flatten()
        gene_vars = np.array(adata.X.power(2).mean(axis=0)).flatten() - gene_means**2
    else:
        gene_means = adata.X.mean(axis=0)
        gene_vars = adata.X.var(axis=0)
    
    # Calculate dispersion
    gene_dispersions = gene_vars / (gene_means + 1e-12)
    
    # Select top genes
    top_gene_indices = np.argsort(gene_dispersions)[-n_top_genes:]
    
    adata.var['highly_variable'] = False
    adata.var.iloc[top_gene_indices, adata.var.columns.get_loc('highly_variable')] = True
    
    n_hvg = adata.var['highly_variable'].sum()
    print(f"  Selected {n_hvg:,} highly variable genes")
    
    del gene_means, gene_vars, gene_dispersions
    gc.collect()
    
    return adata

# ============================================================================
# Complete QC Pipeline
# ============================================================================

def run_complete_qc_pipeline(adata, qc_params, n_top_genes=2000):
    """
    Run complete quality control pipeline
    
    Parameters:
        adata: AnnData object
        qc_params: Dictionary of QC parameters
        n_top_genes: Number of highly variable genes
    
    Returns:
        adata: QC'd AnnData object
    """
    
    print("="*80)
    print("Complete Quality Control Pipeline")
    print("="*80)
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    print(f"\nInitial state:")
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")
    
    # Step 1: Calculate QC metrics
    adata = calculate_qc_metrics(adata)
    
    # Step 2: Filter cells
    adata = filter_cells_by_qc(
        adata,
        min_genes=qc_params['min_genes'],
        max_genes=qc_params['max_genes'],
        max_pct_mito=qc_params['max_mito_pct'],
        min_counts=qc_params['min_counts'],
        max_counts=qc_params['max_counts']
    )
    
    # Step 3: Filter genes
    adata = filter_genes_by_expression(adata, min_cells=qc_params['min_cells'])
    
    # Step 4: Normalize and log-transform
    adata = normalize_and_log_transform(adata, target_sum=1e4)
    
    # Step 5: Select highly variable genes
    adata = select_highly_variable_genes(adata, n_top_genes=n_top_genes)
    
    # Step 6: Compute PCA (CRITICAL for downstream analysis)
    print(f"\n{'='*80}")
    print("Computing PCA")
    print(f"{'='*80}")
    
    n_hvg = adata.var['highly_variable'].sum()
    print(f"  Using {n_hvg:,} highly variable genes")
    
    # Check for NaN/Inf values and replace with 0
    if issparse(adata.X):
        # For sparse matrices, check data array
        if np.any(~np.isfinite(adata.X.data)):
            print("  Warning: Found NaN/Inf in sparse matrix, replacing with 0...")
            adata.X.data = np.nan_to_num(adata.X.data, nan=0.0, posinf=0.0, neginf=0.0)
    else:
        # For dense matrices
        if np.any(~np.isfinite(adata.X)):
            print("  Warning: Found NaN/Inf in matrix, replacing with 0...")
            adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Determine n_comps: must be <= min(n_cells, n_hvg)
    max_comps = min(adata.n_obs, n_hvg)
    n_comps = min(50, max_comps)
    
    print(f"  Computing {n_comps} principal components (max possible: {max_comps})...")
    
    sc.tl.pca(
        adata,
        n_comps=n_comps,
        use_highly_variable=True,
        svd_solver='arpack',
        random_state=42
    )
    
    print(f"  ✓ PCA computed: {adata.obsm['X_pca'].shape}")
    n_show = min(5, n_comps)
    print(f"  Variance explained (first {n_show} PCs): {adata.uns['pca']['variance_ratio'][:n_show]}")
    
    print("\n" + "="*80)
    print("Quality Control Complete")
    print("="*80)
    print(f"Final state:")
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")
    print(f"  Highly variable genes: {adata.var['highly_variable'].sum():,}")
    print(f"  PCA: {adata.obsm['X_pca'].shape}")
    
    return adata

# ============================================================================
# Execute if run directly
# ============================================================================

if __name__ == '__main__':
    print("\n" + "="*80)
    print("Testing QC Module")
    print("="*80)
    
    print("Module loaded successfully")
