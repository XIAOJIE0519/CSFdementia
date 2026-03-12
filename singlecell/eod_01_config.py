"""
EOD (Early-Onset Dementia) Single-Cell RNA-seq Integration Pipeline
====================================================================

Configuration and Setup Module

This pipeline integrates 63 samples from EODsample.csv following Nature-level
standards for multi-cohort single-cell analysis.

Scientific Rationale:
- Focus on early-onset dementia (EOD) samples for age-matched comparison
- Include matched healthy controls for differential analysis
- Preserve biological heterogeneity across brain regions
- Enable cross-dataset batch effect correction

Author: Bioinformatics Analysis
Date: 2026-02-14
Version: 1.0
"""

import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import os
import sys
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Disable numba to avoid NumPy 2.4 compatibility issues
os.environ['NUMBA_DISABLE_JIT'] = '1'
os.environ['NUMBA_DISABLE_PERFORMANCE_WARNINGS'] = '1'

# Set scanpy settings for reproducibility
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white', frameon=False)
sc.settings.n_jobs = 1  # Avoid parallel processing issues

print("="*80)
print("EOD SINGLE-CELL INTEGRATION PIPELINE v1.0")
print("Nature-Level Quality Standards")
print("="*80)
print("\nStep 1: Configuration and Setup")
print("-"*80)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Data directories
DATA_DIR = Path("./")
OUTPUT_DIR = Path("./EOD")
OUTPUT_DIR.mkdir(exist_ok=True)

# Sample metadata file (63 EOD-related samples)
SAMPLE_CSV = DATA_DIR / "EODsample.csv"

# Dataset file patterns (same as main pipeline)
DATASET_PATTERNS = {
    'GSE250280': {'pattern': '*.h5', 'format': '10x_h5'},
    'GSE272082': {'pattern': '*_filtered_feature_bc_matrix.h5', 'format': '10x_h5'},
    'Rexach et al.': {'pattern': 'Rexach.h5ad', 'format': 'h5ad'}
}

# Quality control thresholds (Nature-level stringent)
QC_PARAMS = {
    'min_genes': 200,        # Minimum genes per cell
    'min_cells': 3,          # Minimum cells per gene  
    'max_genes': 10000,      # Maximum genes per cell (doublet filter)
    'max_mito_pct': 20,      # Maximum mitochondrial percentage
    'min_counts': 500,       # Minimum UMI counts per cell
    'max_counts': 100000     # Maximum UMI counts per cell
}

# Gene ID conversion settings
GENE_ID_PARAMS = {
    'remove_linc': True,      # Remove LINC genes
    'merge_duplicates': True,  # Merge duplicate gene names (sum)
    'species': 'human'
}

# Highly variable genes and dimensionality reduction
N_TOP_GENES = 2000          # Highly variable genes
N_PCS = 30                  # Principal components
N_NEIGHBORS = 15            # Neighbors for graph construction

# Clustering parameters
CLUSTERING_PARAMS = {
    'resolution': 1.0,       # Leiden clustering resolution
    'n_clusters_kmeans': 20  # Fallback KMeans clusters
}

# Memory optimization
MEMORY_EFFICIENT = True
CHUNK_SIZE = 10000          # Process cells in chunks for large operations

# Parallel processing for neighbors computation
N_JOBS = 8                  # Number of parallel jobs for neighbors (set to CPU cores)

print(f"\nConfiguration:")
print(f"  Sample metadata: {SAMPLE_CSV}")
print(f"  Output directory: {OUTPUT_DIR}")
print(f"  Expected samples: 63")
print(f"\nQuality control parameters:")
for key, value in QC_PARAMS.items():
    print(f"  {key}: {value}")
print(f"\nGene ID conversion:")
for key, value in GENE_ID_PARAMS.items():
    print(f"  {key}: {value}")
print(f"\nDimensionality reduction:")
print(f"  Highly variable genes: {N_TOP_GENES}")
print(f"  Principal components: {N_PCS}")
print(f"  Neighbors: {N_NEIGHBORS}")

# ============================================================================
# LOAD AND VALIDATE SAMPLE METADATA
# ============================================================================

print(f"\n{'='*80}")
print("Loading EOD Sample Metadata")
print(f"{'='*80}")

if not SAMPLE_CSV.exists():
    print(f"ERROR: Sample metadata file not found: {SAMPLE_CSV}")
    sys.exit(1)

# Load sample metadata
sample_df = pd.read_csv(SAMPLE_CSV)

# Standardize column names
sample_df.columns = sample_df.columns.str.strip().str.lower().str.replace(' ', '_')

print(f"\nLoaded {len(sample_df)} samples")
print(f"Columns: {sample_df.columns.tolist()}")

# Validate required columns
required_cols = ['sample_id', 'data_source', 'main_type', 'subtype']
missing_cols = [col for col in required_cols if col not in sample_df.columns]
if missing_cols:
    print(f"ERROR: Missing required columns: {missing_cols}")
    sys.exit(1)

# Summary by data source
print(f"\nSamples by data source:")
source_counts = sample_df['data_source'].value_counts()
for source, count in source_counts.items():
    print(f"  {source}: {count} samples")

# Summary by disease type
print(f"\nSamples by disease type:")
disease_counts = sample_df['subtype'].value_counts()
for disease, count in disease_counts.items():
    print(f"  {disease}: {count} samples")

# Summary by main type
print(f"\nSamples by main type:")
main_type_counts = sample_df['main_type'].value_counts()
for main_type, count in main_type_counts.items():
    print(f"  {main_type}: {count} samples")

# Check for Rexach duplicates (same donor_id, different positions)
rexach_samples = sample_df[sample_df['data_source'] == 'Rexach et al.']
if len(rexach_samples) > 0:
    print(f"\nRexach et al. samples: {len(rexach_samples)}")
    
    # Check for duplicate sample_ids with different positions
    rexach_dup_ids = rexach_samples['sample_id'].value_counts()
    rexach_dup_ids = rexach_dup_ids[rexach_dup_ids > 1]
    
    if len(rexach_dup_ids) > 0:
        print(f"\nRexach samples with multiple brain regions:")
        for sample_id, count in rexach_dup_ids.items():
            positions = rexach_samples[rexach_samples['sample_id'] == sample_id]['position'].tolist()
            subtype = rexach_samples[rexach_samples['sample_id'] == sample_id]['subtype'].iloc[0]
            print(f"  {sample_id} ({subtype}): {count} regions - {positions}")

# Age distribution
if 'age' in sample_df.columns:
    print(f"\nAge distribution:")
    print(f"  Mean: {sample_df['age'].mean():.1f} years")
    print(f"  Median: {sample_df['age'].median():.1f} years")
    print(f"  Range: {sample_df['age'].min():.0f} - {sample_df['age'].max():.0f} years")
    
    # Age by disease
    print(f"\nAge by disease subtype:")
    age_by_disease = sample_df.groupby('subtype')['age'].agg(['mean', 'median', 'count'])
    print(age_by_disease)

# Sex distribution
if 'sex' in sample_df.columns:
    print(f"\nSex distribution:")
    sex_counts = sample_df['sex'].value_counts()
    for sex, count in sex_counts.items():
        print(f"  {sex}: {count} samples ({count/len(sample_df)*100:.1f}%)")

# Brain region distribution
if 'position' in sample_df.columns:
    print(f"\nBrain region distribution:")
    position_counts = sample_df['position'].value_counts()
    for position, count in position_counts.items():
        print(f"  {position}: {count} samples")

print("\n" + "="*80)
print("Configuration completed successfully")
print("="*80)
print(f"\nReady to process {len(sample_df)} EOD-related samples")
print(f"  Disease samples: {len(sample_df[sample_df['main_type'] == 'Disease'])}")
print(f"  Control samples: {len(sample_df[sample_df['main_type'] == 'HC'])}")
