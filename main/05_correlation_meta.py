#!/usr/bin/env python3
"""
Protein-Protein Correlation Meta-Analysis
Calculates pairwise protein correlations within each study using ALL samples
Performs Fisher's Z-transformation and sample-size weighted z-score meta-analysis
Following the same methodology as meta2.R
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import Dict, List
import warnings
from scipy import stats
from sklearn.impute import KNNImputer
warnings.filterwarnings('ignore')

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
log = logging.getLogger(__name__)

# Study priority rules (same as meta2.R)
STUDY_PRIORITY = {
    "study_1": 1,
    "study_2": 2,
    "study_5": 3,
    "study_3_Soma": 4,
    "study_3_Olink": 5,
    "study_10": 6,
    "study_8": 10,
    "study_9": 11,
    "study_4": 20,
    "study_6": 21,
    "study_7": 22,
    "study_11": 23
}

# Duplicate groups (same as meta2.R)
DUPLICATE_GROUP1 = ["study_1", "study_2", "study_5", "study_3_Soma", "study_3_Olink", "study_10"]
DUPLICATE_GROUP2 = ["study_8", "study_9"]

N_IMPUTATIONS = 12  # Number of KNN imputations


def get_study_priority(study_name: str) -> int:
    """Get study priority (lower = higher priority)."""
    return STUDY_PRIORITY.get(study_name, 100)


def select_studies_for_pair(pair_data: pd.DataFrame) -> pd.DataFrame:
    """
    Apply study selection rules based on duplicate groups.
    Same logic as meta2.R select_studies_for_protein().
    """
    # Group 1: if multiple studies present, keep only highest priority
    in_group1 = pair_data['Study'].isin(DUPLICATE_GROUP1)
    if in_group1.sum() >= 2:
        group1_data = pair_data[in_group1]
        priorities = group1_data['Study'].apply(get_study_priority)
        keep_study = priorities.idxmin()
        pair_data = pair_data[~in_group1 | (pair_data.index == keep_study)]
    
    # Group 2: if multiple studies present, keep only highest priority
    in_group2 = pair_data['Study'].isin(DUPLICATE_GROUP2)
    if in_group2.sum() >= 2:
        group2_data = pair_data[in_group2]
        priorities = group2_data['Study'].apply(get_study_priority)
        keep_study = priorities.idxmin()
        pair_data = pair_data[~in_group2 | (pair_data.index == keep_study)]
    
    return pair_data


def count_unique_studies(study_names: List[str]) -> int:
    """
    Count unique non-duplicate studies.
    Same logic as meta2.R count_unique_studies().
    """
    unique_count = 0
    
    # Check if any study from group1 is present
    if any(s in DUPLICATE_GROUP1 for s in study_names):
        unique_count += 1
    
    # Check if any study from group2 is present
    if any(s in DUPLICATE_GROUP2 for s in study_names):
        unique_count += 1
    
    # Count other studies (not in any duplicate group)
    other_studies = [s for s in study_names if s not in DUPLICATE_GROUP1 + DUPLICATE_GROUP2]
    unique_count += len(set(other_studies))
    
    return unique_count


def calculate_study_correlations(expression_df: pd.DataFrame, study_name: str) -> pd.DataFrame:
    """
    Calculate pairwise protein correlations using ALL samples.
    Uses KNN imputation (12 times with different random seeds) and Pearson correlation.
    """
    log.info(f"  Processing {study_name}...")
    
    # Quality control: remove proteins with >50% missing
    missing_pct = expression_df.isna().sum() / len(expression_df)
    valid_proteins = missing_pct[missing_pct < 0.5].index
    expression_df = expression_df[valid_proteins]
    
    n_proteins = len(expression_df.columns)
    n_samples = len(expression_df)
    
    log.info(f"    Samples: {n_samples}, Proteins: {n_proteins}")
    
    if n_proteins < 2:
        log.warning(f"    Too few proteins, skipping")
        return pd.DataFrame()
    
    proteins = expression_df.columns.tolist()
    
    # Perform KNN imputation N_IMPUTATIONS times with different random seeds
    log.info(f"    Performing KNN imputation ({N_IMPUTATIONS} times)...")
    
    # Convert to numpy for faster computation
    expr_array = expression_df.values
    
    # Adjust k for small sample sizes (k must be < n_samples)
    knn_k = min(5, max(1, n_samples - 1))
    
    all_corr_matrices = []
    
    import sys
    for imp_idx in range(N_IMPUTATIONS):
        # KNN imputation with k=5 and different random seed
        np.random.seed(imp_idx)
        imputer = KNNImputer(n_neighbors=knn_k, weights='uniform')
        expr_imputed = imputer.fit_transform(expr_array)
    
        # Calculate correlation matrix
        corr_matrix = np.corrcoef(expr_imputed.T)
        all_corr_matrices.append(corr_matrix)
        
        if (imp_idx + 1) % 4 == 0:
            log.info(f"      Completed {imp_idx + 1}/{N_IMPUTATIONS} imputations")
            sys.stdout.flush()
    
    # Average correlation matrices across imputations
    avg_corr_matrix = np.mean(all_corr_matrices, axis=0)
    
    log.info(f"    Extracting correlations...")
    sys.stdout.flush()
    results = []
    
    for i in range(n_proteins):
        if (i + 1) % 500 == 0:
            log.info(f"      Progress: {i+1}/{n_proteins} proteins")
            sys.stdout.flush()
        
        for j in range(i+1, n_proteins):
            r = avg_corr_matrix[i, j]
            
            # Skip invalid correlations
            if np.isnan(r) or np.isinf(r):
                continue
            
            # Clip to valid range
            r = np.clip(r, -0.9999, 0.9999)
            
            # Calculate p-value using t-distribution
            if abs(r) >= 0.9999:
                pval = 0.0
            else:
                t_stat = r * np.sqrt((n_samples - 2) / (1 - r**2))
                pval = 2 * (1 - stats.t.cdf(abs(t_stat), df=n_samples-2))
            
            results.append({
                'Protein1': proteins[i],
                'Protein2': proteins[j],
                'Pearson': r,
                'P_value': pval,
                'N': n_samples,
                'Study': study_name
            })
    
    df_results = pd.DataFrame(results)
    log.info(f"    Found {len(df_results)} protein pairs")
    sys.stdout.flush()
    
    return df_results


def fisher_z_transform(r: float) -> float:
    """Fisher's Z transformation for correlation coefficients."""
    r = np.clip(r, -0.9999, 0.9999)
    return 0.5 * np.log((1 + r) / (1 - r))


def inverse_fisher_z(z: float) -> float:
    """Inverse Fisher's Z transformation."""
    return (np.exp(2 * z) - 1) / (np.exp(2 * z) + 1)


def perform_weighted_meta(pearson_values: np.ndarray, n_samples: np.ndarray, p_values: np.ndarray) -> Dict:
    """
    Perform sample-size weighted average and weighted z-score meta-analysis.
    Same methodology as meta2.R perform_weighted_meta().
    
    Returns meta-analyzed correlation (back-transformed from Fisher's Z).
    """
    # Remove NA values
    valid = ~(np.isnan(pearson_values) | np.isnan(n_samples) | np.isnan(p_values))
    pearson_values = pearson_values[valid]
    n_samples = n_samples[valid]
    p_values = p_values[valid]
    
    if len(pearson_values) < 2:
        return {
            'pearson_meta': np.nan,
            'ci_lower': np.nan,
            'ci_upper': np.nan,
            'z_score': np.nan,
            'p_value': np.nan,
            'i_squared': np.nan
        }
    
    # Fisher Z-transform correlations
    z_values = np.array([fisher_z_transform(r) for r in pearson_values])
    
    # Standard errors for Fisher's Z (1/sqrt(n-3))
    se_values = np.array([1.0 / np.sqrt(max(n - 3, 1)) for n in n_samples])
    
    # Sample-size weighted average of Fisher's Z
    weights_n = n_samples / np.sum(n_samples)
    weighted_z = np.sum(weights_n * z_values)
    
    # Standard error of weighted effect
    se_weighted = np.sqrt(np.sum((weights_n * se_values) ** 2))
    
    # 95% confidence interval (in Z scale)
    ci_lower_z = weighted_z - 1.96 * se_weighted
    ci_upper_z = weighted_z + 1.96 * se_weighted
    
    # Weighted z-score method
    # Convert p-values to z-scores with direction
    z_scores = stats.norm.ppf(p_values / 2, loc=0, scale=1)
    z_scores = np.abs(z_scores) * np.sign(pearson_values)
    
    # Weight by square root of sample size
    weights_z = np.sqrt(n_samples)
    
    # Weighted z-score
    z_meta = np.sum(weights_z * z_scores) / np.sqrt(np.sum(weights_z ** 2))
    
    # Meta-analytic p-value
    p_meta = 2 * (1 - stats.norm.cdf(abs(z_meta)))
    
    # Calculate I-squared for heterogeneity
    # Using inverse variance weights
    weights_iv = 1 / (se_values ** 2)
    q_stat = np.sum(weights_iv * (z_values - weighted_z) ** 2)
    df = len(z_values) - 1
    i_squared = max(0, (q_stat - df) / q_stat) if q_stat > 0 else 0
    
    # Back-transform to correlation scale
    pearson_meta = inverse_fisher_z(weighted_z)
    ci_lower_r = inverse_fisher_z(ci_lower_z)
    ci_upper_r = inverse_fisher_z(ci_upper_z)
        
    return {
        'pearson_meta': pearson_meta,
        'ci_lower': ci_lower_r,
        'ci_upper': ci_upper_r,
        'z_score': z_meta,
        'p_value': p_meta,
        'i_squared': i_squared
    }


def perform_loo_analysis(pearson_values: np.ndarray, n_samples: np.ndarray) -> Dict:
    """
    Perform leave-one-out sensitivity analysis.
    Same methodology as meta2.R perform_loo_analysis().
    """
    valid = ~(np.isnan(pearson_values) | np.isnan(n_samples))
    pearson_values = pearson_values[valid]
    n_samples = n_samples[valid]
    
    if len(pearson_values) < 3:
        return {'mean': np.nan, 'ci_lower': np.nan, 'ci_upper': np.nan}
    
    # Fisher Z-transform
    z_values = np.array([fisher_z_transform(r) for r in pearson_values])
    
    loo_effects = []
    
    for i in range(len(z_values)):
        loo_z = np.delete(z_values, i)
        loo_n = np.delete(n_samples, i)
        weights = loo_n / np.sum(loo_n)
        loo_z_weighted = np.sum(weights * loo_z)
        # Back-transform to correlation
        loo_effects.append(inverse_fisher_z(loo_z_weighted))
    
    return {
        'mean': np.mean(loo_effects),
        'ci_lower': np.min(loo_effects),
        'ci_upper': np.max(loo_effects)
    }


def perform_correlation_meta_analysis(df_corr: pd.DataFrame) -> pd.DataFrame:
    """
    Perform meta-analysis on protein-protein correlations.
    Groups by protein pair and applies sample-size weighted z-score method.
    Same methodology as meta2.R.
    """
    log.info(f"\nPerforming meta-analysis...")
    log.info(f"  Total correlation records: {len(df_corr)}")
    
    # Create protein pair identifier (sorted to ensure consistency)
    df_corr['Pair'] = df_corr.apply(
        lambda row: tuple(sorted([row['Protein1'], row['Protein2']])), 
        axis=1
    )
    
    # Group by protein pair
    protein_pairs = df_corr.groupby('Pair')
    log.info(f"  Unique protein pairs: {len(protein_pairs)}")
    
    results = []
    
    for i, (pair, df_pair_original) in enumerate(protein_pairs):
        if (i + 1) % 100000 == 0:
            log.info(f"    Progress: {i+1}/{len(protein_pairs)} ({100*(i+1)/len(protein_pairs):.1f}%)")
        
        # Save original study information BEFORE deduplication
        n_studies_original = len(df_pair_original)
        total_n_original = df_pair_original['N'].sum()
        studies_list_original = '; '.join(df_pair_original['Study'].tolist())
        
        # Apply study selection rules (deduplication)
        df_pair = select_studies_for_pair(df_pair_original)
        
        # Skip if less than 2 studies after deduplication
        if len(df_pair) < 2:
            continue
        
        # Calculate unique non-duplicate study count
        unique_study_count = count_unique_studies(df_pair['Study'].tolist())
        
        # Perform weighted meta-analysis (using deduplicated data)
        meta = perform_weighted_meta(
            df_pair['Pearson'].values,
            df_pair['N'].values,
            df_pair['P_value'].values
        )
        
        # Skip if meta result is NA
        if np.isnan(meta['pearson_meta']):
            continue
        
        # Perform leave-one-out analysis
        loo = perform_loo_analysis(
            df_pair['Pearson'].values,
            df_pair['N'].values
        )
        
        # Create unique study names list (after deduplication)
        unique_study_names = '; '.join(df_pair['Study'].unique())
        
        # Store results
        results.append({
            'Protein1': pair[0],
            'Protein2': pair[1],
            'Pearson_meta': meta['pearson_meta'],
            'CI_Lower': meta['ci_lower'],
            'CI_Upper': meta['ci_upper'],
            'Z_score': meta['z_score'],
            'P_value': meta['p_value'],
            'I_squared': meta['i_squared'],
            'LOO_Mean': loo['mean'],
            'LOO_CI_Lower': loo['ci_lower'],
            'LOO_CI_Upper': loo['ci_upper'],
            'N_Studies': n_studies_original,
            'Total_N': total_n_original,
            'Studies': studies_list_original,
            'N_Unique_Studies': unique_study_count,
            'Unique_Study_Names': unique_study_names
        })
    
    df_results = pd.DataFrame(results)
    
    log.info(f"  Protein pairs with >=2 studies: {len(df_results)}")
    
    if len(df_results) == 0:
        return df_results
    
    # Filter proteins with N_Unique_Studies >= 2
    df_results = df_results[df_results['N_Unique_Studies'] >= 2].copy()
    
    log.info(f"  Protein pairs with >=2 unique studies: {len(df_results)}")
    
    if len(df_results) == 0:
        return df_results
    
    # Stratified FDR correction by number of unique studies
    df_results['FDR_BH_Stratified'] = np.nan
    
    for n_unique in df_results['N_Unique_Studies'].unique():
        idx = df_results['N_Unique_Studies'] == n_unique
        if idx.sum() > 0:
            from scipy.stats import false_discovery_control
            try:
                df_results.loc[idx, 'FDR_BH_Stratified'] = false_discovery_control(
                    df_results.loc[idx, 'P_value'].values,
                    method='bh'
                )
            except:
                # Fallback to manual BH
                pvals = df_results.loc[idx, 'P_value'].values
                n = len(pvals)
                sorted_indices = np.argsort(pvals)
                sorted_pvals = pvals[sorted_indices]
                
                fdr_manual = np.zeros(n)
                for i in range(n):
                    fdr_manual[sorted_indices[i]] = min(1.0, sorted_pvals[i] * n / (i + 1))
                
                # Enforce monotonicity
                for i in range(n - 2, -1, -1):
                    if fdr_manual[sorted_indices[i]] > fdr_manual[sorted_indices[i + 1]]:
                        fdr_manual[sorted_indices[i]] = fdr_manual[sorted_indices[i + 1]]
                
                df_results.loc[idx, 'FDR_BH_Stratified'] = fdr_manual
    
    # Sort by p-value
    df_results = df_results.sort_values('P_value')
    
    log.info(f"  Completed: {len(df_results)} protein pairs analyzed")
    
    return df_results


def process_diagnosis_group(df_all: pd.DataFrame, group_name: str, diagnosis_list: List[str], output_base_dir: Path):
    """
    Process a specific diagnosis group (EOD, LOD, or CN).
    Calculate within-study correlations and perform meta-analysis.
    """
    log.info("\n" + "="*80)
    log.info(f"PROCESSING GROUP: {group_name}")
    log.info(f"Diagnoses: {', '.join(diagnosis_list)}")
    log.info("="*80)
    
    # Create output directory for this group
    output_dir = output_base_dir / group_name
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Filter samples by diagnosis
    df_group = df_all[df_all['Diagnosis_Derived'].isin(diagnosis_list)].copy()
    log.info(f"  Total samples in {group_name}: {len(df_group)}")
    
    if len(df_group) < 10:
        log.warning(f"  Insufficient samples for {group_name}, skipping")
        return
    
    # Get studies present in this group
    studies = sorted(df_group['Study'].unique())
    log.info(f"  Studies with {group_name} samples: {', '.join(studies)}")
    
    # Calculate within-study correlations
    log.info(f"\nCalculating within-study correlations for {group_name}...")
    
    all_correlations = []
    
    for study in studies:
        df_study = df_group[df_group['Study'] == study].copy()
        
        log.info(f"  {study}: {len(df_study)} {group_name} samples")
        
        # Find protein columns (start from first column without clinical keywords)
        clinical_keywords = ['Study', 'Batch', 'GUID', 'Age', 'Sex', 'Cognitive', 'AB', 'Tau', 'NEFL', 'YKL40', 'Diagnosis']
        protein_start_idx = None
        for i, col in enumerate(df_study.columns):
            if not any(keyword.lower() in col.lower() for keyword in clinical_keywords):
                protein_start_idx = i
                break
        
        if protein_start_idx is None:
            log.warning(f"    No protein columns found, skipping")
            continue
        
        # Extract expression data (protein columns only)
        expression_df = df_study.iloc[:, protein_start_idx:].copy()
        
        # Keep only numeric columns (exclude any remaining meta/string columns)
        expression_df = expression_df.select_dtypes(include=[np.number])
        
        # Remove all-NA proteins
        expression_df = expression_df.dropna(axis=1, how='all')
        
        if len(expression_df) < 3 or len(expression_df.columns) < 2:
            log.warning(f"    Insufficient data (n={len(expression_df)}), skipping")
            continue
        
        # Calculate correlations with KNN imputation
        # Skip if already computed
        output_path = output_dir / f"{study}_correlations.csv"
        if output_path.exists():
            log.info(f"    Already computed, loading: {output_path.name}")
            df_corr = pd.read_csv(output_path)
        else:
            df_corr = calculate_study_correlations(expression_df, study)
        
        if len(df_corr) > 0:
            all_correlations.append(df_corr)
            
            # Save study-specific results
            if not output_path.exists():
                df_corr.to_csv(output_path, index=False)
            log.info(f"    Saved: {output_path.name}")
    
    if not all_correlations:
        log.warning(f"  No correlations calculated for {group_name}")
        return
    
    # Check if meta-analysis already done
    meta_path = output_dir / "correlation_meta_analysis.csv"
    if meta_path.exists():
        log.info(f"  Meta-analysis already exists, skipping: {meta_path.name}")
        return
    
    # Combine all correlations
    log.info(f"\nPerforming meta-analysis for {group_name}...")
    
    df_all_corr = pd.concat(all_correlations, ignore_index=True)
    log.info(f"  Total correlation records: {len(df_all_corr)}")
    
    # Perform meta-analysis
    df_meta = perform_correlation_meta_analysis(df_all_corr)
    
    if df_meta is not None and len(df_meta) > 0:
        # Save meta-analysis results
        meta_path = output_dir / "correlation_meta_analysis.csv"
        df_meta.to_csv(meta_path, index=False)
        log.info(f"  Saved: {meta_path.name}")
        
        # Summary statistics
        log.info(f"\n  Summary for {group_name}:")
        log.info(f"    Total protein pairs: {len(df_meta)}")
        log.info(f"    Significant (P<0.05): {(df_meta['P_value'] < 0.05).sum()}")
        log.info(f"    Significant (FDR<0.05): {(df_meta['FDR_BH_Stratified'] < 0.05).sum()}")
        log.info(f"    Mean N_Studies: {df_meta['N_Studies'].mean():.2f}")
        log.info(f"    Mean N_Unique_Studies: {df_meta['N_Unique_Studies'].mean():.2f}")
        log.info(f"    Mean Total_N: {df_meta['Total_N'].mean():.0f}")
        log.info(f"    Mean I-squared: {df_meta['I_squared'].mean():.1f}%")
        
        # Top correlations
        if len(df_meta) >= 5:
            log.info(f"\n  Top 5 positive correlations in {group_name}:")
            top_pos = df_meta.nlargest(5, 'Pearson_meta')
            for _, row in top_pos.iterrows():
                log.info(f"    {row['Protein1']} - {row['Protein2']}: "
                        f"r={row['Pearson_meta']:.3f}, P={row['P_value']:.2e}, "
                        f"N_Unique={row['N_Unique_Studies']}")
            
            log.info(f"\n  Top 5 negative correlations in {group_name}:")
            top_neg = df_meta.nsmallest(5, 'Pearson_meta')
            for _, row in top_neg.iterrows():
                log.info(f"    {row['Protein1']} - {row['Protein2']}: "
                        f"r={row['Pearson_meta']:.3f}, P={row['P_value']:.2e}, "
                        f"N_Unique={row['N_Unique_Studies']}")


def main():
    """Main analysis pipeline."""
    log.info("="*80)
    log.info("PROTEIN-PROTEIN CORRELATION META-ANALYSIS BY DIAGNOSIS GROUP")
    log.info("Method: Sample-size weighted average + Weighted z-score")
    log.info("KNN imputation: 12 times")
    log.info("Groups: EOD, LOD, CN")
    log.info("="*80)
    
    # Setup paths
    combine_dir = Path("F:/1a-EOD-CSF-protein/combine")
    output_base_dir = Path("F:/1a-EOD-CSF-protein/correlation_meta")
    output_base_dir.mkdir(exist_ok=True)
    
    expr_path = combine_dir / "combined_expression_matrices.csv"
    
    if not expr_path.exists():
        log.error(f"Expression data not found: {expr_path}")
        return
    
    # Load data
    log.info(f"\nLoading expression data...")
    df_all = pd.read_csv(expr_path, low_memory=False)
    log.info(f"  Loaded {len(df_all)} samples")
    
    # Check diagnosis distribution
    log.info(f"\nDiagnosis distribution:")
    for diag, count in df_all['Diagnosis_Derived'].value_counts().items():
        log.info(f"  {diag}: {count}")
    
    # Define diagnosis groups
    diagnosis_groups = {
        'EOD': ['EOAD', 'EOOD', 'EOFTD', 'EODSD', 'EODLB'],
        'LOD': ['LOAD', 'LOOD', 'LOFTD', 'LODLB'],
        'CN': ['CN']
    }
    
    # Process each diagnosis group
    for group_name, diagnosis_list in diagnosis_groups.items():
        process_diagnosis_group(df_all, group_name, diagnosis_list, output_base_dir)
    
    log.info("\n" + "="*80)
    log.info("ANALYSIS COMPLETE")
    log.info("="*80)
    log.info(f"Output directory: {output_base_dir}")
    log.info(f"Method: Sample-size weighted average + Weighted z-score")
    log.info(f"Groups processed: {', '.join(diagnosis_groups.keys())}")
    log.info("="*80)


if __name__ == "__main__":
    import time
    start_time = time.time()
    main()
    elapsed = time.time() - start_time
    log.info(f"\nTotal execution time: {elapsed:.1f}s ({elapsed/60:.1f}min)")
