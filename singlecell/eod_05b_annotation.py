"""
EOD Single-Cell Integration Pipeline - Cell Type Annotation Module (Revised)
============================================================================

Scientific cell type annotation using:
1. sc.tl.score_genes for marker expression scoring
2. Leiden cluster-level voting mechanism
3. Remove small clusters (<50 cells) after annotation

This approach follows Nature/Cell standards:
- Cluster-level annotation (not single-cell argmax)
- Statistical robustness through voting
- Biological interpretability

Reference:
Hao et al. (2021) Cell - Seurat v4
Traag et al. (2019) Scientific Reports - Leiden algorithm

Author: Bioinformatics Analysis
Date: 2026-02-17
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("\n" + "="*80)
print("EOD Pipeline - Cell Type Annotation (Cluster-based)")
print("="*80)

# ============================================================================
# Marker Gene Definitions (Updated Format)
# ============================================================================

MARKER_GENES = {
    'ExcitatoryNeurons': ['SLC17A7', 'CAMK2A', 'SATB2'],
    'InhibitoryNeurons': ['GAD1', 'GAD2', 'SLC32A1'],
    'Astrocytes': ['AQP4', 'SLC1A2', 'ALDH1L1', 'GFAP'],
    'Oligodendrocytes': ['MOBP', 'MOG', 'PLP1', 'MBP'],
    'OPCs': ['PDGFRA', 'CSPG4', 'VCAN'],
    'Microglia': ['C3', 'CSF1R', 'P2RY12', 'C1QA'],
    'Macrophages': ['MRC1', 'CD163', 'LYVE1', 'CD206'],
    'EndothelialCells': ['CLDN5', 'FLT1', 'PECAM1'],
    'Pericytes': ['PDGFRB', 'DCN'],
    'Fibroblasts': ['COL1A1', 'COL1A2', 'LUM']
}

def load_marker_genes_from_dict(marker_dict, adata, verbose=True):
    """
    Load and validate marker genes from dictionary
    
    Parameters:
        marker_dict: Dictionary of {celltype: [genes]}
        adata: AnnData object
        verbose: Print validation info
    
    Returns:
        validated_markers: Dictionary with only genes present in data
    """
    
    if verbose:
        print(f"\n{'='*80}")
        print("Loading and Validating Marker Genes")
        print(f"{'='*80}")
        print(f"Total cell types: {len(marker_dict)}")
    
    validated_markers = {}
    available_genes = set(adata.var_names)
    
    for celltype, genes in marker_dict.items():
        # Filter to available genes
        valid_genes = [g for g in genes if g in available_genes]
        
        if len(valid_genes) > 0:
            validated_markers[celltype] = valid_genes
            
            if verbose:
                missing = len(genes) - len(valid_genes)
                status = "✓" if missing == 0 else "⚠"
                print(f"  {status} {celltype}: {len(valid_genes)}/{len(genes)} markers found")
                if missing > 0:
                    missing_genes = [g for g in genes if g not in available_genes]
                    print(f"      Missing: {missing_genes}")
        else:
            if verbose:
                print(f"  ✗ {celltype}: No markers found in data")
    
    if verbose:
        print(f"\nValidated cell types: {len(validated_markers)}")
    
    return validated_markers

def score_cell_types(adata, marker_dict, verbose=True):
    """
    Score each cell for each cell type using sc.tl.score_genes
    
    Parameters:
        adata: AnnData object
        marker_dict: Dictionary of validated markers
        verbose: Print progress
    
    Returns:
        adata: AnnData with scores in adata.obs
    
    Scientific Notes:
    - sc.tl.score_genes computes mean expression of marker genes
    - Subtracts background (random gene set)
    - Robust to library size differences
    """
    
    if verbose:
        print(f"\n{'='*80}")
        print("Scoring Cell Types")
        print(f"{'='*80}")
    
    for celltype, genes in marker_dict.items():
        score_name = f'score_{celltype}'
        
        try:
            sc.tl.score_genes(
                adata,
                gene_list=genes,
                score_name=score_name,
                use_raw=False
            )
            
            if verbose:
                mean_score = adata.obs[score_name].mean()
                print(f"  ✓ {celltype}: mean score = {mean_score:.3f}")
        
        except Exception as e:
            if verbose:
                print(f"  ✗ {celltype}: Error - {e}")
    
    return adata

def annotate_clusters_by_voting(adata, marker_dict, cluster_key='leiden', 
                                min_cluster_size=50, verbose=True):
    """
    Annotate clusters using voting mechanism, then remove small clusters
    
    For each cluster:
    1. Calculate mean score for each cell type
    2. Assign cell type with highest mean score
    3. Mark small clusters as "Other cells"
    4. Remove cells marked as "Other cells"
    
    Parameters:
        adata: AnnData with scores
        marker_dict: Marker dictionary
        cluster_key: Cluster column
        min_cluster_size: Minimum cells per cluster
        verbose: Print results
    
    Returns:
        adata: AnnData with 'cell_type' annotation (Other cells removed)
    
    Scientific Rationale:
    - Cluster-level annotation is more robust than single-cell
    - Voting reduces noise from dropout events
    - Small clusters likely represent doublets or low-quality cells
    """
    
    if verbose:
        print(f"\n{'='*80}")
        print("Cluster-Level Cell Type Annotation")
        print(f"{'='*80}")
        print(f"Cluster key: {cluster_key}")
        print(f"Minimum cluster size: {min_cluster_size} cells")
    
    # Get cluster sizes
    cluster_sizes = adata.obs[cluster_key].value_counts()
    
    if verbose:
        print(f"\nTotal clusters: {len(cluster_sizes)}")
        print(f"Clusters >= {min_cluster_size} cells: {(cluster_sizes >= min_cluster_size).sum()}")
        print(f"Small clusters (< {min_cluster_size}): {(cluster_sizes < min_cluster_size).sum()}")
    
    # Initialize annotation
    adata.obs['cell_type'] = 'Unassigned'
    
    # Get score columns
    score_cols = [f'score_{ct}' for ct in marker_dict.keys()]
    
    # Annotate each cluster
    cluster_annotations = {}
    
    for cluster in adata.obs[cluster_key].unique():
        cluster_mask = adata.obs[cluster_key] == cluster
        cluster_size = cluster_mask.sum()
        
        # Mark small clusters
        if cluster_size < min_cluster_size:
            cluster_annotations[cluster] = 'Other cells'
            if verbose:
                print(f"  Cluster {cluster} ({cluster_size:,} cells): Other cells (too small)")
            continue
        
        # Calculate mean scores for this cluster
        cluster_scores = {}
        for celltype in marker_dict.keys():
            score_col = f'score_{celltype}'
            if score_col in adata.obs.columns:
                mean_score = adata.obs.loc[cluster_mask, score_col].mean()
                cluster_scores[celltype] = mean_score
        
        # Assign cell type with highest score
        if cluster_scores:
            best_celltype = max(cluster_scores, key=cluster_scores.get)
            best_score = cluster_scores[best_celltype]
            
            cluster_annotations[cluster] = best_celltype
            
            if verbose:
                print(f"  Cluster {cluster} ({cluster_size:,} cells): {best_celltype} (score={best_score:.3f})")
        else:
            cluster_annotations[cluster] = 'Other cells'
    
    # Apply annotations
    adata.obs['cell_type'] = adata.obs[cluster_key].map(cluster_annotations)
    
    # Summary before filtering
    if verbose:
        print(f"\n{'='*80}")
        print("Annotation Summary (Before Filtering)")
        print(f"{'='*80}")
        
        celltype_counts = adata.obs['cell_type'].value_counts()
        print(f"\nCell type distribution:")
        for celltype, count in celltype_counts.items():
            pct = count / adata.n_obs * 100
            print(f"  {celltype}: {count:,} cells ({pct:.1f}%)")
    
    # Remove "Other cells"
    n_before = adata.n_obs
    other_cells_mask = adata.obs['cell_type'] == 'Other cells'
    n_other = other_cells_mask.sum()
    
    if n_other > 0:
        if verbose:
            print(f"\n{'='*80}")
            print("Filtering Small Clusters")
            print(f"{'='*80}")
            print(f"Removing {n_other:,} cells marked as 'Other cells'...")
        
        # CRITICAL: Remove obsp matrices before subsetting to avoid memory explosion
        # obsp contains neighbor graphs (n_cells x n_cells) which cause huge memory allocation during indexing
        if len(adata.obsp.keys()) > 0:
            if verbose:
                print(f"  Clearing obsp matrices: {list(adata.obsp.keys())}")
            adata.obsp.clear()
        
        # Now safe to subset
        adata = adata[~other_cells_mask].copy()
        
        if verbose:
            print(f"  ✓ Removed {n_other:,} cells")
            print(f"  Remaining cells: {adata.n_obs:,} ({adata.n_obs/n_before*100:.1f}%)")
            
            print(f"\nFinal cell type distribution:")
            celltype_counts_final = adata.obs['cell_type'].value_counts()
            for celltype, count in celltype_counts_final.items():
                pct = count / adata.n_obs * 100
                print(f"  {celltype}: {count:,} cells ({pct:.1f}%)")
    
    return adata

def validate_annotations(adata, marker_dict, verbose=True):
    """
    Validate annotations by checking marker expression
    
    Parameters:
        adata: Annotated AnnData
        marker_dict: Marker dictionary
        verbose: Print validation
    
    Returns:
        validation_df: DataFrame with validation metrics
    """
    
    if verbose:
        print(f"\n{'='*80}")
        print("Annotation Validation")
        print(f"{'='*80}")
    
    validation_results = []
    
    for celltype in adata.obs['cell_type'].unique():
        if celltype == 'Other cells' or celltype == 'Unassigned':
            continue
        
        if celltype not in marker_dict:
            continue
        
        # Get cells of this type
        celltype_mask = adata.obs['cell_type'] == celltype
        n_cells = celltype_mask.sum()
        
        # Check marker expression
        markers = marker_dict[celltype]
        available_markers = [m for m in markers if m in adata.var_names]
        
        if len(available_markers) == 0:
            continue
        
        # Calculate mean expression
        marker_expr = adata[:, available_markers].X
        if hasattr(marker_expr, 'toarray'):
            marker_expr = marker_expr.toarray()
        
        celltype_expr = marker_expr[celltype_mask, :].mean(axis=0)
        other_expr = marker_expr[~celltype_mask, :].mean(axis=0)
        
        fold_change = np.mean(celltype_expr) / (np.mean(other_expr) + 1e-10)
        
        validation_results.append({
            'cell_type': celltype,
            'n_cells': n_cells,
            'n_markers': len(available_markers),
            'mean_expr_in_type': np.mean(celltype_expr),
            'mean_expr_other': np.mean(other_expr),
            'fold_enrichment': fold_change
        })
        
        if verbose:
            print(f"  {celltype}: {n_cells:,} cells, {fold_change:.2f}x enrichment")
    
    validation_df = pd.DataFrame(validation_results)
    
    return validation_df

# ============================================================================
# Main Annotation Pipeline
# ============================================================================

def run_cell_type_annotation(adata, marker_dict=None, cluster_key='leiden',
                            min_cluster_size=50, verbose=True):
    """
    Complete cell type annotation pipeline
    
    Parameters:
        adata: AnnData with clustering
        marker_dict: Marker gene dictionary (uses default if None)
        cluster_key: Cluster column
        min_cluster_size: Minimum cluster size
        verbose: Print progress
    
    Returns:
        adata: Annotated AnnData (with Other cells removed)
        validation_df: Validation metrics
    """
    
    if marker_dict is None:
        marker_dict = MARKER_GENES
    
    print("="*80)
    print("Cell Type Annotation Pipeline")
    print("="*80)
    
    # Step 1: Load and validate markers
    validated_markers = load_marker_genes_from_dict(marker_dict, adata, verbose)
    
    if len(validated_markers) == 0:
        raise ValueError("No valid markers found in data")
    
    # Step 2: Score cell types
    adata = score_cell_types(adata, validated_markers, verbose)
    
    # Step 3: Annotate clusters by voting and remove Other cells
    adata = annotate_clusters_by_voting(
        adata,
        validated_markers,
        cluster_key=cluster_key,
        min_cluster_size=min_cluster_size,
        verbose=verbose
    )
    
    # Step 4: Validate annotations
    validation_df = validate_annotations(adata, validated_markers, verbose)
    
    print("\n" + "="*80)
    print("✓ Cell Type Annotation Complete")
    print("="*80)
    
    return adata, validation_df

# ============================================================================
# Execute if run directly
# ============================================================================

if __name__ == '__main__':
    print("\n" + "="*80)
    print("Testing Cell Type Annotation Module")
    print("="*80)
    
    print("Module loaded successfully")
    print(f"\nDefault marker genes:")
    for ct, genes in MARKER_GENES.items():
        print(f"  {ct}: {len(genes)} markers")
