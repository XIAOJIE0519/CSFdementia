"""
EOD Single-Cell Integration Pipeline - Gene ID Conversion Module
=================================================================

Convert ENSG IDs to Gene Symbols and process gene names.

This module:
1. Converts ENSG IDs to Gene Symbols using mygene
2. Merges duplicate gene names (sum counts)
3. Removes LINC genes
4. Maintains detailed conversion logs

Scientific Standards:
- Document all gene ID conversions
- Handle conversion failures appropriately
- Preserve biological information
- Ensure reproducibility

Author: Bioinformatics Analysis
Date: 2026-02-14
"""

import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import issparse, csr_matrix, lil_matrix
from datetime import datetime
import gc

print("\n" + "="*80)
print("EOD Pipeline - Step 3: Gene ID Conversion")
print("="*80)
print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# ============================================================================
# Gene ID Conversion Functions
# ============================================================================

def convert_ensg_to_symbol_batch(gene_list, verbose=True):
    """
    Convert ENSG IDs to Gene Symbols using mygene
    
    Parameters:
        gene_list: List of gene IDs (mix of ENSG and symbols)
        verbose: Print progress
    
    Returns:
        conversion_dict: Dictionary mapping old ID to new symbol
        failed_genes: List of genes that failed conversion
    """
    
    try:
        import mygene
    except ImportError:
        print("ERROR: mygene not installed")
        print("Please run: pip install mygene")
        return None, None
    
    if verbose:
        print(f"\n{'='*80}")
        print("Gene ID Conversion (ENSG → Gene Symbol)")
        print(f"{'='*80}")
        print(f"Total genes: {len(gene_list):,}")
    
    # Separate ENSG and non-ENSG genes
    ensg_genes = [g for g in gene_list if str(g).startswith('ENSG')]
    non_ensg_genes = [g for g in gene_list if not str(g).startswith('ENSG')]
    
    if verbose:
        print(f"  ENSG format: {len(ensg_genes):,}")
        print(f"  Non-ENSG format: {len(non_ensg_genes):,}")
    
    conversion_dict = {}
    failed_genes = []
    
    # Non-ENSG genes keep their names
    for gene in non_ensg_genes:
        conversion_dict[gene] = gene
    
    if len(ensg_genes) == 0:
        if verbose:
            print("  No ENSG genes to convert")
        return conversion_dict, failed_genes
    
    # Convert ENSG genes using mygene
    if verbose:
        print(f"\nConverting {len(ensg_genes):,} ENSG IDs...")
        print("  This may take a few minutes...")
    
    mg = mygene.MyGeneInfo()
    
    # Query in batches
    batch_size = 1000
    all_results = []
    
    for i in range(0, len(ensg_genes), batch_size):
        batch = ensg_genes[i:i+batch_size]
        
        if verbose and len(ensg_genes) > batch_size:
            print(f"  Batch {i//batch_size + 1}/{(len(ensg_genes)-1)//batch_size + 1}...")
        
        try:
            results = mg.querymany(
                batch,
                scopes='ensembl.gene',
                fields='symbol',
                species='human',
                returnall=True
            )
            all_results.extend(results['out'])
        except Exception as e:
            print(f"  ERROR in batch {i//batch_size + 1}: {e}")
            # Add failed genes
            failed_genes.extend(batch)
    
    # Parse results
    for result in all_results:
        query = result['query']
        
        if 'symbol' in result:
            conversion_dict[query] = result['symbol']
        else:
            failed_genes.append(query)
    
    if verbose:
        n_converted = len([g for g in ensg_genes if g in conversion_dict])
        print(f"\nConversion results:")
        print(f"  Successfully converted: {n_converted:,} / {len(ensg_genes):,}")
        print(f"  Failed: {len(failed_genes):,}")
        
        if len(failed_genes) > 0 and len(failed_genes) <= 10:
            print(f"  Failed genes: {failed_genes}")
        elif len(failed_genes) > 10:
            print(f"  Failed genes (first 10): {failed_genes[:10]}")
    
    return conversion_dict, failed_genes

def apply_gene_conversion_to_adata(adata, conversion_dict, remove_failed=True):
    """
    Apply gene ID conversion to AnnData object
    
    Parameters:
        adata: AnnData object
        conversion_dict: Conversion dictionary
        remove_failed: Remove genes that failed conversion
    
    Returns:
        adata: Updated AnnData object
    """
    
    print(f"\nApplying gene ID conversion...")
    print(f"  Original genes: {adata.n_vars:,}")
    
    # Save original gene IDs
    adata.var['original_gene_id'] = adata.var_names.copy()
    
    # Convert gene names
    new_gene_names = []
    failed_indices = []
    
    for i, gene in enumerate(adata.var_names):
        if gene in conversion_dict:
            new_gene_names.append(conversion_dict[gene])
        else:
            new_gene_names.append(gene)
            if str(gene).startswith('ENSG'):
                failed_indices.append(i)
    
    # Update gene names
    adata.var_names = new_gene_names
    
    # Remove failed genes if requested
    if remove_failed and len(failed_indices) > 0:
        print(f"  Removing {len(failed_indices)} failed conversion genes")
        keep_indices = [i for i in range(adata.n_vars) if i not in failed_indices]
        adata = adata[:, keep_indices].copy()
    
    print(f"  Genes after conversion: {adata.n_vars:,}")
    
    return adata

def merge_duplicate_genes_efficient(adata, method='sum'):
    """
    Merge duplicate gene names efficiently using sparse matrix operations
    
    Parameters:
        adata: AnnData object
        method: 'sum' or 'mean'
    
    Returns:
        adata: AnnData with merged genes
    """
    
    print(f"\n{'='*80}")
    print("Merging Duplicate Genes")
    print(f"{'='*80}")
    
    # Check for duplicates
    gene_counts = pd.Series(adata.var_names).value_counts()
    duplicates = gene_counts[gene_counts > 1]
    
    print(f"Original genes: {adata.n_vars:,}")
    print(f"Unique genes: {len(gene_counts):,}")
    print(f"Duplicate genes: {len(duplicates)}")
    
    if len(duplicates) == 0:
        print("  No duplicate genes found")
        return adata
    
    print(f"\nDuplicate genes (showing first 10):")
    for gene, count in duplicates.head(10).items():
        print(f"  {gene}: {count} occurrences")
    
    # Fast merge using sparse matrix operations
    print(f"\nMerging using '{method}' method (sparse matrix)...")
    
    # Create mapping from old indices to new indices
    unique_genes = gene_counts.index.tolist()
    gene_to_new_idx = {gene: i for i, gene in enumerate(unique_genes)}
    
    # Build aggregation matrix (old_genes x new_genes)
    # Each column sums the corresponding old gene columns
    from scipy.sparse import csr_matrix, lil_matrix
    
    n_old = adata.n_vars
    n_new = len(unique_genes)
    
    # Use lil_matrix for efficient construction
    agg_matrix = lil_matrix((n_old, n_new))
    
    for old_idx, gene in enumerate(adata.var_names):
        new_idx = gene_to_new_idx[gene]
        agg_matrix[old_idx, new_idx] = 1
    
    agg_matrix = agg_matrix.tocsr()
    
    # Multiply: (cells x old_genes) @ (old_genes x new_genes) = (cells x new_genes)
    if issparse(adata.X):
        X_merged = adata.X @ agg_matrix
    else:
        X_merged = csr_matrix(adata.X) @ agg_matrix
    
    # Handle mean
    if method == 'mean':
        # Divide by gene counts
        gene_count_array = np.array([gene_counts[g] for g in unique_genes])
        X_merged = X_merged / gene_count_array[np.newaxis, :]
    
    # Create new AnnData
    adata_merged = sc.AnnData(
        X=X_merged,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=unique_genes)
    )
    
    print(f"\nMerging complete:")
    print(f"  Final genes: {adata_merged.n_vars:,}")
    print(f"  Removed: {adata.n_vars - adata_merged.n_vars:,} duplicate entries")
    
    del agg_matrix, X_merged
    gc.collect()
    
    return adata_merged

def remove_linc_genes(adata):
    """
    Remove LINC genes from dataset
    
    Parameters:
        adata: AnnData object
    
    Returns:
        adata: AnnData without LINC genes
    """
    
    print(f"\n{'='*80}")
    print("Removing LINC Genes")
    print(f"{'='*80}")
    
    linc_genes = [g for g in adata.var_names if str(g).startswith('LINC')]
    
    print(f"LINC genes found: {len(linc_genes):,} / {adata.n_vars:,}")
    
    if len(linc_genes) == 0:
        print("  No LINC genes to remove")
        return adata
    
    if len(linc_genes) <= 20:
        print(f"LINC genes: {linc_genes}")
    else:
        print(f"LINC genes (first 20): {linc_genes[:20]}")
    
    # Remove LINC genes
    keep_genes = [g for g in adata.var_names if not str(g).startswith('LINC')]
    adata = adata[:, keep_genes].copy()
    
    print(f"\nRemaining genes: {adata.n_vars:,}")
    
    return adata

# ============================================================================
# Main Processing Function
# ============================================================================

def process_gene_ids_complete(adata, remove_linc=True, merge_duplicates=True):
    """
    Complete gene ID processing pipeline
    
    Steps:
    1. Convert ENSG → Gene Symbol
    2. Merge duplicate genes (MUST be before removing LINC)
    3. Remove LINC genes
    
    Parameters:
        adata: AnnData object
        remove_linc: Remove LINC genes
        merge_duplicates: Merge duplicate gene names
    
    Returns:
        adata: Processed AnnData object
    """
    
    print("="*80)
    print("Gene ID Processing Pipeline")
    print("="*80)
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"\nInitial state:")
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")
    
    # Step 1: Convert ENSG IDs
    gene_list = adata.var_names.tolist()
    conversion_dict, failed_genes = convert_ensg_to_symbol_batch(gene_list, verbose=True)
    
    if conversion_dict is not None:
        adata = apply_gene_conversion_to_adata(adata, conversion_dict, remove_failed=True)
    
    # Step 2: Merge duplicates (MUST be before removing LINC)
    if merge_duplicates:
        adata = merge_duplicate_genes_efficient(adata, method='sum')
    
    # Step 3: Remove LINC genes
    if remove_linc:
        adata = remove_linc_genes(adata)
    
    print("\n" + "="*80)
    print("Gene ID Processing Complete")
    print("="*80)
    print(f"Final state:")
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")
    
    return adata

# ============================================================================
# Execute if run directly
# ============================================================================

if __name__ == '__main__':
    print("\n" + "="*80)
    print("Testing Gene ID Conversion Module")
    print("="*80)
    
    # Test with sample genes
    test_genes = [
        'ENSG00000130203',  # APOE
        'ENSG00000142192',  # APP
        'APOE',
        'APP',
        'LINC00996'
    ]
    
    print(f"Test genes: {test_genes}")
    
    conversion_dict, failed = convert_ensg_to_symbol_batch(test_genes, verbose=True)
    
    if conversion_dict:
        print(f"\nConversion results:")
        for gene, symbol in conversion_dict.items():
            print(f"  {gene} → {symbol}")
