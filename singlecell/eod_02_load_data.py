"""
EOD Single-Cell Integration Pipeline - Data Loading Module
===========================================================

Load and preprocess single-cell data from multiple sources for EOD analysis.

This module:
1. Loads data from GSE250280, GSE272082, and Rexach et al.
2. Filters samples based on EODsample.csv
3. Adds comprehensive metadata
4. Performs initial quality checks

Scientific Standards:
- Preserve raw counts for downstream analysis
- Maintain sample provenance and batch information
- Document all filtering steps
- Ensure reproducibility

Author: Bioinformatics Analysis
Date: 2026-02-14
"""

import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import os
import gc
from glob import glob
from pathlib import Path
from datetime import datetime

print("\n" + "="*80)
print("EOD Pipeline - Step 2: Data Loading")
print("="*80)
print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# ============================================================================
# Helper Functions
# ============================================================================

def load_10x_h5_data(h5_file, sample_id, dataset_name):
    """
    Load 10x Genomics h5 file with error handling
    
    Parameters:
        h5_file: Path to h5 file
        sample_id: Sample identifier
        dataset_name: Dataset name
    
    Returns:
        adata: AnnData object or None if loading fails
    """
    try:
        adata = sc.read_10x_h5(h5_file)
        adata.obs['sample_id'] = sample_id
        adata.obs['data_source'] = dataset_name
        adata.var_names_make_unique()
        return adata
    except Exception as e:
        print(f"  ERROR loading {h5_file}: {e}")
        return None

def load_h5ad_data(h5ad_file, dataset_name):
    """
    Load h5ad file (for Rexach data)
    
    Parameters:
        h5ad_file: Path to h5ad file
        dataset_name: Dataset name
    
    Returns:
        adata: AnnData object or None if loading fails
    """
    try:
        adata = sc.read_h5ad(h5ad_file)
        adata.obs['data_source'] = dataset_name
        return adata
    except Exception as e:
        print(f"  ERROR loading {h5ad_file}: {e}")
        return None

# ============================================================================
# Main Loading Function
# ============================================================================

def load_eod_datasets(sample_metadata):
    """
    Load all datasets specified in EODsample.csv
    
    Parameters:
        sample_metadata: DataFrame from EODsample.csv
    
    Returns:
        all_datasets: List of AnnData objects
        loading_log: DataFrame with loading statistics
    """
    
    print("\n" + "="*80)
    print("Loading EOD Datasets")
    print("="*80)
    
    all_datasets = []
    loading_log = []
    
    # Get unique data sources from sample metadata
    data_sources = sample_metadata['data_source'].unique()
    print(f"\nData sources to load: {list(data_sources)}")
    
    # ========================================================================
    # 1. Load GSE250280
    # ========================================================================
    
    if 'GSE250280' in data_sources:
        print(f"\n{'='*80}")
        print("Loading GSE250280")
        print(f"{'='*80}")
        
        # Get samples from this dataset
        gse250280_samples = sample_metadata[sample_metadata['data_source'] == 'GSE250280']
        print(f"  Expected samples: {len(gse250280_samples)}")
        
        # Find h5 files
        h5_files = glob("./GSE250280/*.h5")
        print(f"  Found {len(h5_files)} h5 files")
        
        dataset_adatas = []
        
        for h5_file in sorted(h5_files):
            filename = os.path.basename(h5_file)
            sample_gsm = filename.split('.')[0]
            
            # Check if this sample is in our EOD list
            if sample_gsm not in gse250280_samples['sample_id'].values:
                continue
            
            # Load data
            adata = load_10x_h5_data(h5_file, sample_gsm, 'GSE250280')
            
            if adata is None:
                continue
            
            # Get metadata for this sample
            sample_meta = gse250280_samples[gse250280_samples['sample_id'] == sample_gsm].iloc[0]
            
            # Add metadata - direct assignment, not fillna
            adata.obs['main_type'] = str(sample_meta['main_type'])
            adata.obs['subtype'] = str(sample_meta['subtype'])
            adata.obs['age'] = sample_meta.get('age', np.nan)
            adata.obs['age_class'] = str(sample_meta.get('age_class', 'NA'))
            adata.obs['sex'] = str(sample_meta.get('sex', 'NA'))
            adata.obs['position'] = str(sample_meta.get('position', 'NA'))
            
            print(f"  Loaded {sample_gsm}: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
            print(f"    Type: {sample_meta['subtype']}, Age: {sample_meta.get('age', 'NA')}, Sex: {sample_meta.get('sex', 'NA')}")
            
            dataset_adatas.append(adata)
            
            loading_log.append({
                'dataset': 'GSE250280',
                'sample_id': sample_gsm,
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'subtype': sample_meta['subtype'],
                'age': sample_meta.get('age', np.nan),
                'sex': sample_meta.get('sex', 'NA')
            })
        
        if len(dataset_adatas) > 0:
            # Concatenate samples
            gse250280_combined = ad.concat(
                dataset_adatas, 
                join='outer',
                label='sample_id',
                keys=[a.obs['sample_id'].iloc[0] for a in dataset_adatas],
                index_unique='-'
            )
            all_datasets.append(gse250280_combined)
            
            print(f"\n  GSE250280 Summary:")
            print(f"    Loaded samples: {len(dataset_adatas)}")
            print(f"    Total cells: {gse250280_combined.n_obs:,}")
            print(f"    Total genes: {gse250280_combined.n_vars:,}")
    
    # ========================================================================
    # 2. Load GSE272082
    # ========================================================================
    
    if 'GSE272082' in data_sources:
        print(f"\n{'='*80}")
        print("Loading GSE272082")
        print(f"{'='*80}")
        
        # Get samples from this dataset
        gse272082_samples = sample_metadata[sample_metadata['data_source'] == 'GSE272082']
        print(f"  Expected samples: {len(gse272082_samples)}")
        
        # Find h5 files
        h5_files = glob("./GSE272082/*_filtered_feature_bc_matrix.h5")
        print(f"  Found {len(h5_files)} filtered h5 files")
        
        dataset_adatas = []
        
        for h5_file in sorted(h5_files):
            filename = os.path.basename(h5_file)
            sample_gsm = filename.split('_')[0]
            
            # Check if this sample is in our EOD list
            if sample_gsm not in gse272082_samples['sample_id'].values:
                continue
            
            # Load data
            adata = load_10x_h5_data(h5_file, sample_gsm, 'GSE272082')
            
            if adata is None:
                continue
            
            # Get metadata for this sample
            sample_meta = gse272082_samples[gse272082_samples['sample_id'] == sample_gsm].iloc[0]
            
            # Add metadata - direct assignment, not fillna
            adata.obs['main_type'] = str(sample_meta['main_type'])
            adata.obs['subtype'] = str(sample_meta['subtype'])
            adata.obs['age'] = sample_meta.get('age', np.nan)
            adata.obs['age_class'] = str(sample_meta.get('age_class', 'NA'))
            adata.obs['sex'] = str(sample_meta.get('sex', 'NA'))
            adata.obs['position'] = str(sample_meta.get('position', 'NA'))
            
            print(f"  Loaded {sample_gsm}: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
            print(f"    Type: {sample_meta['subtype']}, Age: {sample_meta.get('age', 'NA')}, Sex: {sample_meta.get('sex', 'NA')}")
            
            dataset_adatas.append(adata)
            
            loading_log.append({
                'dataset': 'GSE272082',
                'sample_id': sample_gsm,
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'subtype': sample_meta['subtype'],
                'age': sample_meta.get('age', np.nan),
                'sex': sample_meta.get('sex', 'NA')
            })
        
        if len(dataset_adatas) > 0:
            # Concatenate samples
            gse272082_combined = ad.concat(
                dataset_adatas,
                join='outer',
                label='sample_id',
                keys=[a.obs['sample_id'].iloc[0] for a in dataset_adatas],
                index_unique='-'
            )
            all_datasets.append(gse272082_combined)
            
            print(f"\n  GSE272082 Summary:")
            print(f"    Loaded samples: {len(dataset_adatas)}")
            print(f"    Total cells: {gse272082_combined.n_obs:,}")
            print(f"    Total genes: {gse272082_combined.n_vars:,}")
    
    # ========================================================================
    # 3. Load GSE174367
    # ========================================================================
    
    if 'GSE174367' in data_sources:
        print(f"\n{'='*80}")
        print("Loading GSE174367")
        print(f"{'='*80}")
        
        # Get samples from this dataset
        gse174367_samples = sample_metadata[sample_metadata['data_source'] == 'GSE174367']
        print(f"  Expected samples: {len(gse174367_samples)}")
        
        # Find h5 files
        h5_files = glob("./GSE174367/*.h5")
        print(f"  Found {len(h5_files)} h5 files")
        
        dataset_adatas = []
        
        for h5_file in sorted(h5_files):
            filename = os.path.basename(h5_file)
            sample_gsm = filename.split('.')[0]
            
            # Check if this sample is in our EOD list
            if sample_gsm not in gse174367_samples['sample_id'].values:
                continue
            
            # Load data
            adata = load_10x_h5_data(h5_file, sample_gsm, 'GSE174367')
            
            if adata is None:
                continue
            
            # Get metadata for this sample
            sample_meta = gse174367_samples[gse174367_samples['sample_id'] == sample_gsm].iloc[0]
            
            # Add metadata - direct assignment, not fillna
            adata.obs['main_type'] = str(sample_meta['main_type'])
            adata.obs['subtype'] = str(sample_meta['subtype'])
            adata.obs['age'] = sample_meta.get('age', np.nan)
            adata.obs['age_class'] = str(sample_meta.get('age_class', 'NA'))
            adata.obs['sex'] = str(sample_meta.get('sex', 'NA'))
            adata.obs['position'] = str(sample_meta.get('position', 'NA'))
            
            print(f"  Loaded {sample_gsm}: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
            print(f"    Type: {sample_meta['subtype']}, Age: {sample_meta.get('age', 'NA')}, Sex: {sample_meta.get('sex', 'NA')}")
            
            dataset_adatas.append(adata)
            
            loading_log.append({
                'dataset': 'GSE174367',
                'sample_id': sample_gsm,
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'subtype': sample_meta['subtype'],
                'age': sample_meta.get('age', np.nan),
                'sex': sample_meta.get('sex', 'NA')
            })
        
        if len(dataset_adatas) > 0:
            # Concatenate samples
            gse174367_combined = ad.concat(
                dataset_adatas, 
                join='outer',
                label='sample_id',
                keys=[a.obs['sample_id'].iloc[0] for a in dataset_adatas],
                index_unique='-'
            )
            all_datasets.append(gse174367_combined)
            
            print(f"\n  GSE174367 Summary:")
            print(f"    Loaded samples: {len(dataset_adatas)}")
            print(f"    Total cells: {gse174367_combined.n_obs:,}")
            print(f"    Total genes: {gse174367_combined.n_vars:,}")
    
    # ========================================================================
    # 4. Load Rexach et al. - Load each sample×position individually
    # ========================================================================
    
    if 'Rexach et al.' in data_sources:
        print(f"\n{'='*80}")
        print("Loading Rexach et al.")
        print(f"{'='*80}")
        
        # Get samples from this dataset
        rexach_samples = sample_metadata[sample_metadata['data_source'] == 'Rexach et al.']
        print(f"  Expected samples: {len(rexach_samples)}")
        
        # Load Rexach h5ad file once
        rexach_file = "./Rexach et al 2024/Rexach.h5ad"
        
        if not os.path.exists(rexach_file):
            print(f"  ERROR: Rexach file not found: {rexach_file}")
        else:
            print(f"  Loading Rexach data file...")
            adata_rexach_full = load_h5ad_data(rexach_file, 'Rexach et al.')
            
            if adata_rexach_full is not None:
                print(f"  Total cells in file: {adata_rexach_full.n_obs:,}")
                
                # Check for donor_id column
                if 'donor_id' not in adata_rexach_full.obs.columns:
                    print(f"  ERROR: donor_id column not found in Rexach data")
                else:
                    # Load each sample×position combination separately
                    rexach_adatas = []
                    
                    for idx, sample_row in rexach_samples.iterrows():
                        sample_id = sample_row['sample_id']
                        position = sample_row['position']
                        
                        print(f"  Loading {sample_id} - {position}...")
                        
                        # Filter by donor_id
                        mask = adata_rexach_full.obs['donor_id'] == sample_id
                        
                        # Filter by tissue/position
                        if 'tissue' in adata_rexach_full.obs.columns:
                            mask = mask & (adata_rexach_full.obs['tissue'] == position)
                        
                        # Extract this sample
                        adata_sample = adata_rexach_full[mask].copy()
                        
                        if adata_sample.n_obs == 0:
                            print(f"    WARNING: No cells found for {sample_id} - {position}")
                            loading_log.append({
                                'dataset': 'Rexach et al.',
                                'sample_id': f"{sample_id}_{position}",
                                'n_cells': 0,
                                'n_genes': 0,
                                'subtype': sample_row['subtype'],
                                'age': sample_row.get('age', np.nan),
                                'sex': sample_row.get('sex', 'NA'),
                                'position': position
                            })
                            continue
                        
                        # Add metadata - direct assignment like original code
                        adata_sample.obs['sample_id'] = sample_id
                        adata_sample.obs['data_source'] = 'Rexach et al.'
                        adata_sample.obs['main_type'] = str(sample_row['main_type'])
                        adata_sample.obs['subtype'] = str(sample_row['subtype'])
                        adata_sample.obs['age'] = sample_row.get('age', np.nan)
                        adata_sample.obs['age_class'] = str(sample_row.get('age_class', 'NA'))
                        adata_sample.obs['position'] = str(position)
                        adata_sample.obs['sex'] = str(sample_row.get('sex', 'NA'))
                        
                        print(f"    Loaded: {adata_sample.n_obs:,} cells")
                        
                        rexach_adatas.append(adata_sample)
                        
                        loading_log.append({
                            'dataset': 'Rexach et al.',
                            'sample_id': f"{sample_id}_{position}",
                            'n_cells': adata_sample.n_obs,
                            'n_genes': adata_sample.n_vars,
                            'subtype': sample_row['subtype'],
                            'age': sample_row.get('age', np.nan),
                            'sex': sample_row.get('sex', 'NA'),
                            'position': position
                        })
                    
                    # Concatenate all Rexach samples
                    if len(rexach_adatas) > 0:
                        rexach_combined = ad.concat(
                            rexach_adatas,
                            join='outer',
                            label='sample_id',
                            keys=[f"{s.obs['sample_id'].iloc[0]}_{s.obs['position'].iloc[0]}" for s in rexach_adatas],
                            index_unique='-'
                        )
                        all_datasets.append(rexach_combined)
                        
                        print(f"\n  Rexach et al. Summary:")
                        print(f"    Loaded samples: {len(rexach_adatas)}")
                        print(f"    Total cells: {rexach_combined.n_obs:,}")
                        print(f"    Total genes: {rexach_combined.n_vars:,}")
                        
                        del rexach_adatas
                    
                    # Clean up
                    del adata_rexach_full
                    gc.collect()
    
    # ========================================================================
    # Summary
    # ========================================================================
    
    print("\n" + "="*80)
    print("Data Loading Summary")
    print("="*80)
    
    loading_log_df = pd.DataFrame(loading_log)
    
    print(f"\nTotal datasets loaded: {len(all_datasets)}")
    print(f"Total samples loaded: {len(loading_log_df)}")
    print(f"Total cells (before QC): {sum([d.n_obs for d in all_datasets]):,}")
    
    print(f"\nSamples by dataset:")
    print(loading_log_df.groupby('dataset').size())
    
    print(f"\nSamples by subtype:")
    print(loading_log_df.groupby('subtype').size())
    
    print(f"\nCells by dataset:")
    dataset_cells = loading_log_df.groupby('dataset')['n_cells'].sum()
    for dataset, n_cells in dataset_cells.items():
        print(f"  {dataset}: {n_cells:,} cells")
    
    # Save loading log
    log_file = OUTPUT_DIR / "eod_loading_log.csv"
    loading_log_df.to_csv(log_file, index=False)
    print(f"\nLoading log saved to: {log_file}")
    
    return all_datasets, loading_log_df

# ============================================================================
# Execute if run directly
# ============================================================================

if __name__ == '__main__':
    print("\n" + "="*80)
    print("Testing EOD Data Loading")
    print("="*80)
    
    # This will be called from main pipeline
    print("Module loaded successfully")
