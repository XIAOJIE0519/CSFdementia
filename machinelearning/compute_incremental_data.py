#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compute incremental AUC data for lollipop plot visualization.
Output: incremental_auc_data.csv
"""

import os
import warnings
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import StandardScaler
import joblib

warnings.filterwarnings('ignore')

def direction_aware_roc(y_true, scores):
    fpr, tpr, _ = roc_curve(y_true, scores)
    roc_auc = auc(fpr, tpr)
    if roc_auc < 0.5:
        fpr, tpr, _ = roc_curve(y_true, -scores)
        roc_auc = auc(fpr, tpr)
    return fpr, tpr, roc_auc

def logistic_combine_auc(y, base_vals, add_vals):
    """Fit a 2-feature LogisticRegression (CV) to get optimal combination AUC."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_predict, StratifiedKFold
    try:
        valid = ~(np.isnan(base_vals) | np.isnan(add_vals))
        if valid.sum() < 20 or len(np.unique(y[valid])) < 2:
            return None
        bv = base_vals[valid].astype(float)
        av = add_vals[valid].astype(float)
        yv = y[valid]
        bv = (bv - bv.mean()) / (bv.std() + 1e-9)
        av = (av - av.mean()) / (av.std() + 1e-9)
        X2 = np.column_stack([bv, av])
        lr = LogisticRegression(C=1.0, max_iter=1000, random_state=42)
        # Choose cv folds based on minority class size
        n_minority = int(np.bincount(yv.astype(int)).min())
        n_cv = min(5, n_minority) if n_minority >= 3 else None
        if n_cv is not None:
            cv = StratifiedKFold(n_splits=n_cv, shuffle=True, random_state=42)
            y_prob = cross_val_predict(lr, X2, yv, cv=cv,
                                       method='predict_proba', n_jobs=1)[:, 1]
        else:
            # Too few minority samples for CV: fit on all, report train AUC
            lr.fit(X2, yv)
            y_prob = lr.predict_proba(X2)[:, 1]
        fpr, tpr, _ = roc_curve(yv, y_prob)
        roc_auc = auc(fpr, tpr)
        return roc_auc if roc_auc >= 0.5 else 1.0 - roc_auc
    except Exception:
        return None

# Load & normalise
df_combined = pd.read_csv('../combine/combined_expression_matrices.csv',
                          low_memory=False)
eod_diagnoses = ['EOAD', 'EODSD', 'EOFTD', 'EOOD', 'EODLB']
df_combined['Label'] = df_combined['Diagnosis_Derived'].apply(
    lambda x: 1 if x in eod_diagnoses else 0)
df_combined = df_combined[
    df_combined['Diagnosis_Derived'].isin(eod_diagnoses + ['CN'])].copy()

non_protein_cols = [
    'Study', 'Batch', 'GUID', 'Age', 'Sex', 'Cognitive Score',
    'AB42', 'tTau', 'pTau', 'pTau181', 'AB42/pTau', 'AB40',
    'NEFL', 'YKL40', 'pTau217', 'pTau231', 'Diagnosis_Derived', 'Label',
]
meta_cols = [c for c in non_protein_cols if c in df_combined.columns]

norm_list = []
for study in df_combined['Study'].unique():
    ds = df_combined[df_combined['Study'] == study].copy()
    pcols = [c for c in ds.columns if c not in meta_cols]
    pcols_ok = [c for c in pcols if ds[c].notna().sum() > 0
             and pd.api.types.is_numeric_dtype(ds[c])]
    sc = StandardScaler()
    ds[pcols_ok] = sc.fit_transform(ds[pcols_ok].fillna(0))
    norm_list.append(ds)
df_normalized = pd.concat(norm_list, ignore_index=True)

# Load model
ee_model = joblib.load('./models/unified_LogisticRegression.pkl')
with open('./results/core_features.txt', 'r') as f:
    lines = f.readlines()
core_features = [l.strip().split('. ', 1)[1] for l in lines[2:] if '. ' in l]

study_biomarkers = {
    'study_1':  ['AB42', 'pTau', 'tTau'],
    'study_4':  ['AB42', 'AB40', 'tTau', 'NEFL', 'YKL40',
                 'pTau181', 'pTau217', 'pTau231'],
    'study_9':  ['AB42', 'tTau', 'pTau181'],
    'study_11': ['AB42', 'AB40', 'pTau', 'tTau'],
}
studies_to_analyze = list(study_biomarkers.keys())
INCR_BIOMARKERS = ['pTau', 'tTau', 'pTau181', 'pTau217', 'AB42', 'AB40']

# Collect data
results = []

for study in studies_to_analyze:
    df_s = df_normalized[df_normalized['Study'] == study].copy()
    if len(df_s) == 0 or df_s['Label'].nunique() < 2:
        continue
    y = df_s['Label'].values
    
    # EODstage score
    model_score = None
    try:
        X_s = df_s[core_features].fillna(0).values
        model_score = ee_model.predict_proba(X_s)[:, 1]
    except Exception:
        pass
    
    if model_score is None:
        continue
    
    for bio in INCR_BIOMARKERS:
        if bio not in study_biomarkers.get(study, []):
            continue
        if bio not in df_s.columns:
            continue
        
        bio_vals = df_s[bio].values
        valid = ~np.isnan(bio_vals)
        if valid.sum() < 20:
            continue
        
        # Base AUC
        _, _, auc_base = direction_aware_roc(y[valid], bio_vals[valid])
        
        # Base + EODstage
        auc_comb = logistic_combine_auc(y, bio_vals, model_score)
        
        if auc_comb is not None:
            results.append({
                'Biomarker': bio,
                'Study': study,
                'AUC_Base': auc_base,
                'AUC_Plus_EODstage': auc_comb,
                'Delta': auc_comb - auc_base,
            })

df_results = pd.DataFrame(results)
os.makedirs('./compare', exist_ok=True)
df_results.to_csv('./compare/incremental_auc_data.csv', index=False)
print(f'Saved incremental_auc_data.csv ({len(df_results)} rows)')
print(df_results.head(10))
