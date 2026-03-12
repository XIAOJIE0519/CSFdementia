"""
EOD Single-Cell Integration Pipeline - Batch Integration Module
================================================================

Perform batch effect correction using Harmony integration.

This module addresses technical variation between datasets while preserving
biological heterogeneity. Critical for multi-cohort studies.

Scientific Rationale:
- Three data sources (GSE250280, GSE272082, Rexach) have different:
  * Sequencing platforms
  * Library preparation protocols
  * Tissue processing methods
- Batch effects can confound biological signals
- Harmony: Fast, scalable, preserves cell type structure

Reference:
Korsunsky et al. (2019) Nature Methods
"Fast, sensitive and accurate integration of single-cell data with Harmony"

Author: Bioinformatics Analysis
Date: 2026-02-15
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse
import gc

print("\n" + "="*80)
print("EOD Pipeline - Step 4.5: Batch Integration")
print("="*80)

def run_harmony_integration(adata, batch_key='data_source', n_pcs=30, verbose=True):
    """
    Perform Harmony batch integration on PCA space
    
    Parameters:
        adata: AnnData object with PCA computed
        batch_key: Column in adata.obs for batch information
        n_pcs: Number of PCs to use
        verbose: Print progress
    
    Returns:
        adata: AnnData with corrected PCA in adata.obsm['X_pca_harmony']
    
    Scientific Notes:
    - Harmony iteratively corrects PCA embeddings
    - Preserves local structure within batches
    - Does not over-correct biological variation
    """
    
    if verbose:
        print(f"\n{'='*80}")
        print("Harmony Batch Integration")
        print(f"{'='*80}")
        print(f"Batch key: {batch_key}")
        print(f"Number of PCs: {n_pcs}")
        
        # Check batch distribution
        batch_counts = adata.obs[batch_key].value_counts()
        print(f"\nBatch distribution:")
        for batch, count in batch_counts.items():
            print(f"  {batch}: {count:,} cells")
    
    # Import harmonypy
    import harmonypy as hm
    
    if verbose:
        print("\nRunning Harmony integration with enhanced parameters...")
        print("  This may take 5-10 minutes for large datasets...")
    
    # Extract PCA matrix
    if 'X_pca' not in adata.obsm.keys():
        raise ValueError("PCA not computed. Run PCA first.")
    
    pca_matrix = adata.obsm['X_pca'][:, :n_pcs]
    
    # Prepare batch metadata
    meta_data = adata.obs[[batch_key]].copy()
    
    # Run Harmony with stronger integration parameters
    # theta: 增大theta值可以更强力地去除批次效应（默认2，增加到4）
    # max_iter_harmony: 增加迭代次数确保收敛（默认10，增加到20）
    # sigma: 调整聚类的平滑度
    if verbose:
        print("\n  Enhanced Harmony parameters:")
        print("    theta=4.0 (stronger batch correction, default=2.0)")
        print("    max_iter_harmony=20 (more iterations, default=10)")
        print("    epsilon_harmony=1e-5 (tighter convergence)")
    
    ho = hm.run_harmony(
        pca_matrix,
        meta_data,
        batch_key,
        theta=4.0,              # 增强批次校正强度（默认2.0）
        max_iter_harmony=20,    # 增加迭代次数（默认10）
        epsilon_harmony=1e-5,   # 更严格的收敛标准（默认1e-4）
        verbose=verbose
    )
    
    # Store corrected PCA
    # Harmony returns Z_corr with shape (n_pcs, n_cells)
    # We need (n_cells, n_pcs) for scanpy
    
    # Get corrected PCA from Harmony object
    corrected_pca = ho.Z_corr
    
    # Convert to numpy array if needed
    if not isinstance(corrected_pca, np.ndarray):
        corrected_pca = np.array(corrected_pca)
        if verbose:
            print(f"\n  Converting to numpy array: {corrected_pca.shape}")
    
    if verbose:
        print(f"\n  Harmony output shape: {corrected_pca.shape}")
        print(f"  Expected shape: ({adata.n_obs}, {n_pcs})")
    
    # Check if transpose is needed
    if corrected_pca.shape[0] == n_pcs and corrected_pca.shape[1] == adata.n_obs:
        # Z_corr is (n_pcs, n_cells), need to transpose
        corrected_pca = corrected_pca.T
        if verbose:
            print(f"  Transposing: ({n_pcs}, {adata.n_obs}) -> {corrected_pca.shape}")
    elif corrected_pca.shape[0] == adata.n_obs and corrected_pca.shape[1] == n_pcs:
        # Already in correct shape
        if verbose:
            print(f"  Already correct shape: {corrected_pca.shape}")
    else:
        raise ValueError(f"Unexpected Harmony output shape: {corrected_pca.shape}")
    
    # Verify final shape
    if corrected_pca.shape != (adata.n_obs, n_pcs):
        raise ValueError(f"Corrected PCA shape mismatch: got {corrected_pca.shape}, expected ({adata.n_obs}, {n_pcs})")
    
    adata.obsm['X_pca_harmony'] = corrected_pca
    
    if verbose:
        print("\n✓ Harmony integration completed")
        print(f"  Corrected PCA stored in adata.obsm['X_pca_harmony']")
        # Get number of iterations safely
        try:
            if hasattr(ho, 'objective_harmony'):
                if isinstance(ho.objective_harmony, np.ndarray):
                    n_iters = ho.objective_harmony.shape[0]
                elif isinstance(ho.objective_harmony, list):
                    n_iters = len(ho.objective_harmony)
                else:
                    n_iters = "unknown"
                print(f"  Converged after {n_iters} iterations")
        except:
            pass
    
    del pca_matrix, meta_data, ho
    gc.collect()
    
    return adata

def evaluate_integration_quality(adata, batch_key='data_source', verbose=True):
    """
    Evaluate integration quality using mixing metrics
    
    Parameters:
        adata: AnnData with integrated embeddings
        batch_key: Batch column
        verbose: Print results
    
    Returns:
        metrics: Dictionary of quality metrics
    """
    
    if verbose:
        print(f"\n{'='*80}")
        print("Integration Quality Assessment")
        print(f"{'='*80}")
    
    metrics = {}
    
    # 1. Local batch mixing (kNN-based)
    if 'neighbors' in adata.uns.keys():
        # Calculate batch mixing entropy for each cell
        from scipy.stats import entropy
        
        batch_labels = adata.obs[batch_key].values
        n_batches = len(np.unique(batch_labels))
        
        # Get kNN graph
        knn_indices = adata.uns['neighbors']['indices']
        
        mixing_scores = []
        for i in range(adata.n_obs):
            # Get batch distribution in neighborhood
            neighbor_batches = batch_labels[knn_indices[i]]
            batch_counts = pd.Series(neighbor_batches).value_counts()
            batch_probs = batch_counts / len(neighbor_batches)
            
            # Calculate entropy (higher = better mixing)
            mix_score = entropy(batch_probs, base=n_batches)
            mixing_scores.append(mix_score)
        
        metrics['mean_mixing_entropy'] = np.mean(mixing_scores)
        metrics['median_mixing_entropy'] = np.median(mixing_scores)
        
        if verbose:
            print(f"\nBatch mixing entropy:")
            print(f"  Mean: {metrics['mean_mixing_entropy']:.3f}")
            print(f"  Median: {metrics['median_mixing_entropy']:.3f}")
            print(f"  (Range: 0-1, higher = better mixing)")
    
    # 2. Silhouette score (batch vs cell type)
    if verbose:
        print(f"\n✓ Integration quality assessed")
    
    return metrics

# ============================================================================
# Execute if run directly
# ============================================================================

if __name__ == '__main__':
    print("\n" + "="*80)
    print("Testing Batch Integration Module")
    print("="*80)
    
    print("Module loaded successfully")
