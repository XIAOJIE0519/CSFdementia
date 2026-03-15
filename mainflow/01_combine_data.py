#!/usr/bin/env python3
"""
Data Integration Script for CSF Protein Studies
Combines expression matrices and differential expression results with demographics
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import re
from typing import Dict, List, Tuple

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
log = logging.getLogger(__name__)


def load_demographics(demo_path: Path) -> pd.DataFrame:
    """Load demographics data."""
    log.info("Loading demographics data...")
    df = pd.read_csv(demo_path)
    log.info(f"  Loaded {len(df)} demographic records")
    log.info(f"  First few rows:")
    log.info(f"\n{df.head()}")
    return df


def load_biomarkers(biomarker_path: Path) -> pd.DataFrame:
    """Load biomarker data."""
    log.info("Loading biomarker data...")
    df = pd.read_csv(biomarker_path)
    log.info(f"  Loaded {len(df)} biomarker records")
    log.info(f"  First few rows:")
    log.info(f"\n{df.head()}")
    return df


def extract_study_name(filename: str) -> str:
    """Extract study name from filename."""
    # study_3_Olink_cleaned.csv -> study_3
    # study_1_cleaned.csv -> study_1
    match = re.match(r'(study_\d+)', filename)
    if match:
        return match.group(1)
    return None


def extract_comparison(filename: str) -> str:
    """Extract comparison from differential expression filename (not used in current version)."""
    match = re.search(r'diff_+res_(.+)\.csv$', filename)
    if match:
        return match.group(1)
    return None


def get_study_name_from_file(filename: str) -> str:
    """
    Extract study name from filename, keeping Olink/Soma distinction.
    
    Examples:
    - study_1_cleaned.csv -> study_1
    - study_3_Olink_cleaned.csv -> study_3_Olink
    - study_3_Soma_cleaned.csv -> study_3_Soma
    """
    if 'study_3_Olink' in filename:
        return 'study_3_Olink'
    elif 'study_3_Soma' in filename:
        return 'study_3_Soma'
    else:
        match = re.match(r'(study_\d+)', filename)
        if match:
            return match.group(1)
    return None


def process_expression_matrix(file_path: Path, demographics: pd.DataFrame, 
                              biomarkers: pd.DataFrame, study_name: str) -> pd.DataFrame:
    """
    Process expression matrix file from transform-final directory.
    
    Args:
        file_path: Path to cleaned expression matrix CSV
        demographics: Demographics dataframe
        biomarkers: Biomarkers dataframe
        study_name: Study identifier (e.g., 'study_1', 'study_3_Olink')
    
    Returns:
        Combined dataframe with demographics, biomarkers and protein expression
    """
    log.info(f"  Processing expression matrix: {file_path.name}")
    
    # Load expression data - first column is sample ID, rest are proteins
    df = pd.read_csv(file_path)
    id_col = df.columns[0]
    protein_cols = df.columns[1:]
    
    log.info(f"    Loaded {len(df)} samples with {len(protein_cols)} proteins")
    
    # Extract study number for matching
    study_num_match = re.match(r'study_(\d+)', study_name)
    if not study_num_match:
        log.warning(f"    Cannot extract study number from {study_name}")
        return None
    
    study_num = f"Study{study_num_match.group(1)}"
    study_num_int = study_num_match.group(1)
    
    # Filter demographics for this study
    demo_study = demographics[demographics['Study'] == study_num].copy()
    
    if len(demo_study) == 0:
        log.warning(f"    No demographics found for {study_num}")
        return None
    
    # Convert IDs to string for matching
    df[id_col] = df[id_col].astype(str)
    demo_study['Batch'] = demo_study['Batch'].astype(str)
    demo_study['GUID'] = demo_study['GUID'].astype(str)
    
    # Try Batch matching first, then GUID
    batch_matches = df[id_col].isin(demo_study['Batch'])
    if batch_matches.sum() > 0:
        log.info(f"    Matched {batch_matches.sum()} samples by Batch")
        merge_key = 'Batch'
        df_matched = df[batch_matches].copy()
        df_matched = df_matched.rename(columns={id_col: merge_key})
    else:
        guid_matches = df[id_col].isin(demo_study['GUID'])
        if guid_matches.sum() > 0:
            log.info(f"    Matched {guid_matches.sum()} samples by GUID")
            merge_key = 'GUID'
            df_matched = df[guid_matches].copy()
            df_matched = df_matched.rename(columns={id_col: merge_key})
        else:
            log.warning(f"    No matches found for {study_name}")
            return None
    
    # Merge with demographics
    df_combined = pd.merge(
        df_matched,
        demo_study[['Batch', 'GUID', 'Age', 'Sex', 'Diagnosis_Derived']],
        on=merge_key,
        how='left'
    )
    
    # Now match biomarkers for this study
    biomarker_cols = ['Cognitive Score', 'AB42', 'tTau', 'pTau', 'pTau181', 
                      'AB42/pTau', 'AB40', 'NEFL', 'YKL40', 'pTau217', 'pTau231']
    
    # Filter biomarkers for this study
    bio_study = biomarkers[biomarkers['Study'] == f"study_{study_num_int}"].copy()
    
    if len(bio_study) > 0:
        log.info(f"    Found {len(bio_study)} biomarker records for study_{study_num_int}")
        
        # Convert biomarker IDs to string
        bio_study['ID'] = bio_study['ID'].astype(str)
        
        # Try to match biomarkers by Batch or GUID
        df_combined['Batch_str'] = df_combined['Batch'].astype(str)
        df_combined['GUID_str'] = df_combined['GUID'].astype(str)
        
        # Match by Batch first
        batch_bio_matches = df_combined['Batch_str'].isin(bio_study['ID'])
        guid_bio_matches = df_combined['GUID_str'].isin(bio_study['ID'])
        
        log.info(f"    Biomarker matches: {batch_bio_matches.sum()} by Batch, {guid_bio_matches.sum()} by GUID")
        
        # Create biomarker dataframe for merging - only include columns that exist
        available_bio_cols = [col for col in biomarker_cols if col in bio_study.columns]
        bio_merge = bio_study[['ID'] + available_bio_cols].copy()
        
        # Try Batch matching first
        if batch_bio_matches.sum() > 0:
            bio_merge_batch = bio_merge.rename(columns={'ID': 'Batch_str'})
            df_combined = pd.merge(df_combined, bio_merge_batch, on='Batch_str', how='left')
            log.info(f"    Merged biomarkers by Batch")
        else:
            # Try GUID matching
            bio_merge_guid = bio_merge.rename(columns={'ID': 'GUID_str'})
            df_combined = pd.merge(df_combined, bio_merge_guid, on='GUID_str', how='left')
            log.info(f"    Merged biomarkers by GUID")
        
        # Drop temporary columns
        df_combined = df_combined.drop(columns=['Batch_str', 'GUID_str'])
        
        # Add missing biomarker columns with NaN
        for col in biomarker_cols:
            if col not in df_combined.columns:
                df_combined[col] = np.nan
    else:
        log.info(f"    No biomarker data found for study_{study_num_int}")
        # Add empty biomarker columns
        for col in biomarker_cols:
            df_combined[col] = np.nan
    
    # Add study identifier
    df_combined.insert(0, 'Study', study_name)
    
    # Ensure all biomarker columns exist in the correct order
    for col in biomarker_cols:
        if col not in df_combined.columns:
            df_combined[col] = np.nan
    
    # Reorder columns: Study, Batch, GUID, Age, Sex, 11 Biomarkers, Diagnosis, Proteins
    # Get protein columns (exclude metadata and biomarker columns)
    metadata_cols = ['Study', 'Batch', 'GUID', 'Age', 'Sex', 'Diagnosis_Derived']
    all_meta_bio_cols = metadata_cols + biomarker_cols
    protein_cols_final = [col for col in df_combined.columns if col not in all_meta_bio_cols]
    
    cols_order = ['Study', 'Batch', 'GUID', 'Age', 'Sex'] + biomarker_cols + ['Diagnosis_Derived'] + protein_cols_final
    df_combined = df_combined[cols_order]
    
    log.info(f"    Final combined data: {len(df_combined)} samples")
    
    return df_combined





def main():
    """Main function to combine all data."""
    log.info("="*80)
    log.info("DATA INTEGRATION PIPELINE")
    log.info("="*80)
    
    # Setup paths
    transform_final_dir = Path("transform-final")
    s3_dir             = Path("S3")
    output_dir         = Path("combine")
    output_dir.mkdir(exist_ok=True)

    demographics_path = Path("demographics.csv")
    biomarker_path    = Path("biomarker.csv")

    # Load demographics and biomarkers
    demographics = load_demographics(demographics_path)
    biomarkers   = load_biomarkers(biomarker_path)

    # Collect expression files:
    #   - transform-final/: study_1,2,4,5,6,7,8,9,10,11 (exclude study_3_* and summary)
    #   - S3/: study_3_Olink_cleaned.csv, study_3_Soma_cleaned.csv
    tf_files = sorted([
        f for f in transform_final_dir.glob("*.csv")
        if f.name != "cleaning_summary_report.csv"
        and "study_3" not in f.name
    ])
    s3_files = sorted([
        s3_dir / "study_3_Olink_cleaned.csv",
        s3_dir / "study_3_Soma_cleaned.csv",
    ])
    # Verify S3 files exist
    s3_files = [f for f in s3_files if f.exists()]
    expression_files = tf_files + s3_files

    log.info(f"\nExpression files to process ({len(expression_files)} total):")
    for f in expression_files:
        log.info(f"  {f}")
    
    # Process expression matrices
    log.info("\n" + "="*80)
    log.info("PROCESSING EXPRESSION MATRICES")
    log.info("="*80)
    
    all_expression_data = []
    
    for file_path in expression_files:
        # Extract study name from filename
        study_name = get_study_name_from_file(file_path.name)
        if not study_name:
            log.warning(f"Cannot extract study name from {file_path.name}, skipping")
            continue
        
        df_combined = process_expression_matrix(file_path, demographics, biomarkers, study_name)
        
        if df_combined is not None:
            all_expression_data.append(df_combined)
    
    # Combine all expression data
    if all_expression_data:
        log.info("\n" + "="*80)
        log.info("COMBINING ALL EXPRESSION MATRICES")
        log.info("="*80)
        df_all_expression = pd.concat(all_expression_data, axis=0, ignore_index=True)
        
        output_path = output_dir / "combined_expression_matrices.csv"
        df_all_expression.to_csv(output_path, index=False)
        log.info(f"✓ Saved combined expression data: {output_path}")
        log.info(f"  Total samples: {len(df_all_expression)}")
        log.info(f"  Total columns: {len(df_all_expression.columns)}")
        log.info(f"  Studies: {sorted(df_all_expression['Study'].unique())}")
        
        # Count non-null biomarker values
        biomarker_cols = ['Cognitive Score', 'AB42', 'tTau', 'pTau', 'pTau181', 
                         'AB42/pTau', 'AB40', 'NEFL', 'YKL40', 'pTau217', 'pTau231']
        for col in biomarker_cols:
            non_null = df_all_expression[col].notna().sum()
            log.info(f"  {col}: {non_null} non-null values")
    
    # Summary
    log.info("\n" + "="*80)
    log.info("INTEGRATION SUMMARY")
    log.info("="*80)
    log.info(f"Expression matrices processed: {len(all_expression_data)}")
    log.info(f"Output directory: {output_dir.absolute()}")
    log.info("="*80)
    log.info("✓ DATA INTEGRATION COMPLETE")
    log.info("="*80)


if __name__ == "__main__":
    import time
    start_time = time.time()
    main()
    elapsed = time.time() - start_time
    log.info(f"\nTotal execution time: {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
