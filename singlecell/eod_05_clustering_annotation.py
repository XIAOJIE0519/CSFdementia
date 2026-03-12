"""
EOD Single-Cell Integration Pipeline - Clustering and Annotation Module
========================================================================

Perform dimensionality reduction, clustering, and cell type annotation.

This module:
1. PCA (incremental for large datasets)
2. Neighbor graph construction
3. UMAP (skipped for speed, can add later)
4. Leiden clustering
5. Cell type annotation using marker genes

Scientific Standards:
- Memory-efficient algorithms
- Reproducible random seeds
- Avoid numba dependencies
- Document all parameters

Author: Bioinformatics Analysis
Date: 2026-02-14
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse
from sklearn.decomposition import IncrementalPCA
from sklearn.neighbors import NearestNeighbors
from datetime import datetime
import gc

print("\n" + "="*80)
print("EOD Pipeline - Step 5: Clustering and Annotation")
print("="*80)
print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# ============================================================================
# PCA
# ============================================================================

def compute_pca_incremental(adata, n_comps=30, use_hvg=True):
    """
    Compute PCA using IncrementalPCA for memory efficiency
    
    Parameters:
        adata: AnnData object
        n_comps: Number of principal components
        use_hvg: Use only highly variable genes
    
    Returns:
        adata: AnnData with PCA results
    """
    
    print(f"\n{'='*80}")
    print("Computing PCA")
    print(f"{'='*80}")
    
    # Subset to HVG if requested
    if use_hvg and 'highly_variable' in adata.var.columns:
        adata_pca = adata[:, adata.var['highly_variable']].copy()
        print(f"  Using {adata_pca.n_vars:,} highly variable genes")
    else:
        adata_pca = adata
        print(f"  Using all {adata_pca.n_vars:,} genes")
    
    n_cells = adata_pca.n_obs
    
    # Use IncrementalPCA for large datasets
    if n_cells > 100000:
        print(f"  Using IncrementalPCA for {n_cells:,} cells...")
        
        batch_size = 10000
        ipca = IncrementalPCA(n_components=n_comps, batch_size=batch_size)
        
        # Fit in batches
        n_batches = (n_cells + batch_size - 1) // batch_size
        print(f"  Fitting {n_batches} batches...")
        
        for i in range(0, n_cells, batch_size):
            end_idx = min(i + batch_size, n_cells)
            batch = adata_pca.X[i:end_idx]
            
            if issparse(batch):
                batch = batch.toarray()
            
            ipca.partial_fit(batch)
            
            if (i // batch_size + 1) % 10 == 0:
                print(f"    Batch {i//batch_size + 1}/{n_batches}")
            
            gc.collect()
        
        # Transform in batches
        print("  Transforming data...")
        X_pca = np.zeros((n_cells, n_comps))
        
        for i in range(0, n_cells, batch_size):
            end_idx = min(i + batch_size, n_cells)
            batch = adata_pca.X[i:end_idx]
            
            if issparse(batch):
                batch = batch.toarray()
            
            X_pca[i:end_idx] = ipca.transform(batch)
            gc.collect()
        
        # Store results
        adata.obsm['X_pca'] = X_pca
        adata.varm['PCs'] = np.zeros((adata.n_vars, n_comps))
        if use_hvg:
            adata.varm['PCs'][adata.var['highly_variable'], :] = ipca.components_.T
        
        adata.uns['pca'] = {
            'variance': ipca.explained_variance_,
            'variance_ratio': ipca.explained_variance_ratio_
        }
        
        print(f"  Explained variance (first 5 PCs): {ipca.explained_variance_ratio_[:5]}")
        
        del ipca, X_pca
        
    else:
        print(f"  Using standard PCA for {n_cells:,} cells...")
        from sklearn.decomposition import PCA
        
        X = adata_pca.X
        if issparse(X):
            X = X.toarray()
        
        pca = PCA(n_components=n_comps, random_state=42)
        X_pca = pca.fit_transform(X)
        
        adata.obsm['X_pca'] = X_pca
        adata.varm['PCs'] = np.zeros((adata.n_vars, n_comps))
        if use_hvg:
            adata.varm['PCs'][adata.var['highly_variable'], :] = pca.components_.T
        
        adata.uns['pca'] = {
            'variance': pca.explained_variance_,
            'variance_ratio': pca.explained_variance_ratio_
        }
        
        print(f"  Explained variance (first 5 PCs): {pca.explained_variance_ratio_[:5]}")
        
        del pca, X_pca
    
    gc.collect()
    
    return adata

# ============================================================================
# Neighbor Graph
# ============================================================================

def compute_neighbors(adata, n_neighbors=15, n_pcs=30):
    """
    Compute neighbor graph using sklearn
    
    Parameters:
        adata: AnnData object
        n_neighbors: Number of neighbors
        n_pcs: Number of PCs to use
    
    Returns:
        adata: AnnData with neighbor graph
    """
    
    print(f"\n{'='*80}")
    print("Computing Neighbor Graph")
    print(f"{'='*80}")
    
    print(f"  Using {n_pcs} PCs")
    print(f"  Finding {n_neighbors} nearest neighbors...")
    
    X_pca = adata.obsm['X_pca'][:, :n_pcs]
    
    nn = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto', metric='euclidean', n_jobs=1)
    nn.fit(X_pca)
    
    distances, indices = nn.kneighbors(X_pca)
    
    # Create sparse matrices
    from scipy.sparse import csr_matrix
    n_obs = adata.n_obs
    
    rows = np.repeat(np.arange(n_obs), n_neighbors)
    cols = indices.flatten()
    data = np.ones(len(rows))
    
    connectivities = csr_matrix((data, (rows, cols)), shape=(n_obs, n_obs))
    
    distances_data = distances.flatten()
    distances_matrix = csr_matrix((distances_data, (rows, cols)), shape=(n_obs, n_obs))
    
    adata.obsp['distances'] = distances_matrix
    adata.obsp['connectivities'] = connectivities
    adata.uns['neighbors'] = {
        'connectivities_key': 'connectivities',
        'distances_key': 'distances',
        'params': {'n_neighbors': n_neighbors, 'method': 'sklearn'}
    }
    
    print("  Neighbor graph computed")
    
    del nn, distances, indices, X_pca
    gc.collect()
    
    return adata

# ============================================================================
# UMAP (Optional - Skipped for Speed)
# ============================================================================

def compute_umap_real(adata, skip=False):
    """
    Compute UMAP (now enabled by default)
    
    Parameters:
        adata: AnnData object
        skip: Skip UMAP computation
    
    Returns:
        adata: AnnData with UMAP
    """
    
    print(f"\n{'='*80}")
    print("UMAP Computation")
    print(f"{'='*80}")
    
    if skip:
        print("  Skipping UMAP (only for visualization, not needed for analysis)")
        print("  Creating placeholder...")
        adata.obsm['X_umap'] = np.zeros((adata.n_obs, 2))
    else:
        print("  Computing UMAP...")
        print("  This may take 10-20 minutes for large datasets...")
        
        try:
            import umap
            
            reducer = umap.UMAP(
                n_components=2,
                min_dist=0.3,
                n_neighbors=15,
                metric='euclidean',
                random_state=42,
                n_jobs=1,
                low_memory=True,
                verbose=True
            )
            
            X_umap = reducer.fit_transform(adata.obsm['X_pca'][:, :30])
            adata.obsm['X_umap'] = X_umap
            
            print("  UMAP computed successfully")
            
            del reducer, X_umap
            gc.collect()
            
        except Exception as e:
            print(f"  ERROR computing UMAP: {e}")
            print("  Creating placeholder...")
            adata.obsm['X_umap'] = np.zeros((adata.n_obs, 2))
    
    return adata

# ============================================================================
# Leiden Clustering
# ============================================================================

def compute_leiden_clustering(adata, resolution=1.0):
    """
    Compute Leiden clustering
    
    Parameters:
        adata: AnnData object
        resolution: Clustering resolution
    
    Returns:
        adata: AnnData with cluster labels
    """
    
    print(f"\n{'='*80}")
    print("Leiden Clustering")
    print(f"{'='*80}")
    
    print(f"  Resolution: {resolution}")
    
    try:
        import leidenalg
        import igraph as ig
        
        print("  Building graph...")
        sources, targets = adata.obsp['connectivities'].nonzero()
        weights = adata.obsp['connectivities'].data
        
        g = ig.Graph(directed=False)
        g.add_vertices(adata.n_obs)
        edges = list(zip(sources.tolist(), targets.tolist()))
        g.add_edges(edges)
        g.es['weight'] = weights
        
        print("  Running Leiden algorithm...")
        partition = leidenalg.find_partition(
            g,
            leidenalg.RBConfigurationVertexPartition,
            weights='weight',
            resolution_parameter=resolution,
            seed=42
        )
        
        adata.obs['leiden'] = [str(x) for x in partition.membership]
        
        del g, partition
        
    except ImportError:
        print("  WARNING: leidenalg not available, using KMeans...")
        from sklearn.cluster import KMeans
        
        n_clusters = 20
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        clusters = kmeans.fit_predict(adata.obsm['X_pca'][:, :30])
        adata.obs['leiden'] = [str(x) for x in clusters]
        
        del kmeans
    
    n_clusters = adata.obs['leiden'].nunique()
    print(f"  Identified {n_clusters} clusters")
    
    gc.collect()
    
    return adata

# ============================================================================
# Cell Type Annotation
# ============================================================================

def annotate_cell_types(adata, marker_file):
    """
    Annotate cell types using marker genes
    
    Parameters:
        adata: AnnData object
        marker_file: Path to marker gene file
    
    Returns:
        adata: AnnData with cell type annotations
    """
    
    print(f"\n{'='*80}")
    print("Cell Type Annotation")
    print(f"{'='*80}")
    
    if not marker_file.exists():
        print(f"  WARNING: Marker file not found: {marker_file}")
        print("  Skipping cell type annotation")
        adata.obs['cell_type'] = 'Unknown'
        adata.obs['cell_type_score'] = 0.0
        return adata
    
    # Load marker genes
    print(f"  Loading markers from: {marker_file}")
    marker_df = pd.read_excel(marker_file)
    
    print(f"  Loaded {len(marker_df)} cell types")
    
    # Create marker dictionary
    marker_dict = {}
    
    for idx, row in marker_df.iterrows():
        cell_type = row['cellTypes']
        markers = [str(m).strip() for m in row[1:] if pd.notna(m) and str(m).strip() != '']
        
        # Filter to genes in dataset
        markers_in_data = [m for m in markers if m in adata.var_names]
        
        if len(markers_in_data) > 0:
            marker_dict[cell_type] = markers_in_data
            print(f"  {cell_type}: {len(markers_in_data)}/{len(markers)} markers found")
    
    # Calculate scores (manual implementation)
    print("\n  Computing cell type scores...")
    
    for cell_type, markers in marker_dict.items():
        score_name = f'score_{cell_type.replace(" ", "_").replace("/", "_")}'
        
        # Get marker gene indices
        marker_indices = [i for i, gene in enumerate(adata.var_names) if gene in markers]
        
        if len(marker_indices) > 0:
            # Calculate mean expression
            if issparse(adata.X):
                marker_expr = adata.X[:, marker_indices].toarray()
            else:
                marker_expr = adata.X[:, marker_indices]
            
            scores = np.mean(marker_expr, axis=1)
            adata.obs[score_name] = scores
            
            print(f"    {cell_type}: OK")
        
        gc.collect()
    
    # Assign cell types
    print("\n  Assigning cell types...")
    
    score_columns = [col for col in adata.obs.columns if col.startswith('score_')]
    
    if len(score_columns) == 0:
        print("  WARNING: No scores calculated")
        adata.obs['cell_type'] = 'Unknown'
        adata.obs['cell_type_score'] = 0.0
    else:
        score_matrix = adata.obs[score_columns].values
        max_score_idx = np.argmax(score_matrix, axis=1)
        
        cell_type_names = [col.replace('score_', '').replace('_', ' ') for col in score_columns]
        adata.obs['cell_type'] = [cell_type_names[i] for i in max_score_idx]
        adata.obs['cell_type_score'] = np.max(score_matrix, axis=1)
        
        print("\n  Cell type distribution:")
        print(adata.obs['cell_type'].value_counts())
    
    return adata

# ============================================================================
# Complete Pipeline
# ============================================================================

def run_clustering_and_annotation(adata, n_pcs=30, n_neighbors=15, resolution=1.0, 
                                  marker_file=None, skip_umap=True):
    """
    Run complete clustering and annotation pipeline
    
    Parameters:
        adata: AnnData object
        n_pcs: Number of PCs
        n_neighbors: Number of neighbors
        resolution: Clustering resolution
        marker_file: Path to marker gene file
        skip_umap: Skip UMAP computation
    
    Returns:
        adata: Processed AnnData object
    """
    
    print("="*80)
    print("Clustering and Annotation Pipeline")
    print("="*80)
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Step 1: PCA
    adata = compute_pca_incremental(adata, n_comps=n_pcs, use_hvg=True)
    
    # Step 2: Neighbors
    adata = compute_neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    # Step 3: UMAP (optional)
    adata = compute_umap_real(adata, skip=skip_umap)
    
    # Step 4: Leiden clustering
    adata = compute_leiden_clustering(adata, resolution=resolution)
    
    # Step 5: Cell type annotation
    if marker_file is not None:
        adata = annotate_cell_types(adata, marker_file)
    
    print("\n" + "="*80)
    print("Clustering and Annotation Complete")
    print("="*80)
    
    return adata

# ============================================================================
# Execute if run directly
# ============================================================================

if __name__ == '__main__':
    print("\n" + "="*80)
    print("Testing Clustering Module")
    print("="*80)
    
    print("Module loaded successfully")
