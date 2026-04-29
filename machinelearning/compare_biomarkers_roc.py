#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compare ROC curves: clinical biomarkers alone vs biomarker + EODstage (EasyEnsemble-SVM)
One figure per study, each figure shows all available biomarkers and their +EODstage combined curves.
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import StandardScaler
import joblib

warnings.filterwarnings('ignore')

rcParams.update({
    'font.family':      'serif',
    'font.serif':       ['Times New Roman'],
    'mathtext.fontset': 'stix',
    'font.size':        11,
    'axes.linewidth':   0.8,
    'figure.dpi':       300,
    'pdf.fonttype':     42,
})

OUTPUT_DIR = './compare'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
FIXED_COLORS = {
    'Age':    '#1f77b4', 'Sex':    '#ff7f0e',
    'AB42':   '#17becf', 'AB40':   '#bcbd22',
    'pTau':   '#aec7e8', 'tTau':   '#ffbb78',
    'pTau181':'#98df8a', 'pTau217':'#ff9896',
    'pTau231':'#c5b0d5', 'NEFL':   '#c49c94',
    'YKL40':  '#f7b6d2',
}
MODEL_COLOR = '#d62728'

# Darker/saturated paired colours for +EODstage versions
COMBO_COLORS = {
    'AB42':    '#0a5e6b', 'AB40':    '#6b6b00',
    'pTau':    '#2c5f85', 'tTau':    '#a06000',
    'pTau181': '#2d7a2d', 'pTau217': '#8b0000',
    'pTau231': '#4a3070', 'NEFL':    '#6b3a2a',
    'YKL40':   '#8b3060',
}

def direction_aware_roc(y_true, scores):
    fpr, tpr, _ = roc_curve(y_true, scores)
    roc_auc = auc(fpr, tpr)
    if roc_auc < 0.5:
        fpr, tpr, _ = roc_curve(y_true, -scores)
        roc_auc = auc(fpr, tpr)
    return fpr, tpr, roc_auc

def z_combine_roc(y, vals1, vals2):
    """Combine two predictors via LogisticRegression (CV-adaptive), return (fpr, tpr, auc)."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_predict, StratifiedKFold
    try:
        valid = ~(np.isnan(vals1) | np.isnan(vals2))
        if valid.sum() < 20 or len(np.unique(y[valid])) < 2:
            return None
        v1 = vals1[valid].astype(float)
        v2 = vals2[valid].astype(float)
        yv = y[valid]
        v1 = (v1 - v1.mean()) / (v1.std() + 1e-9)
        v2 = (v2 - v2.mean()) / (v2.std() + 1e-9)
        X2 = np.column_stack([v1, v2])
        lr = LogisticRegression(C=1.0, max_iter=1000, random_state=42)
        n_minority = int(np.bincount(yv.astype(int)).min())
        n_cv = min(5, n_minority) if n_minority >= 3 else None
        if n_cv is not None:
            cv = StratifiedKFold(n_splits=n_cv, shuffle=True, random_state=42)
            y_prob = cross_val_predict(lr, X2, yv, cv=cv,
                                       method='predict_proba', n_jobs=1)[:, 1]
        else:
            lr.fit(X2, yv)
            y_prob = lr.predict_proba(X2)[:, 1]
        fpr, tpr, _ = roc_curve(yv, y_prob)
        roc_auc = auc(fpr, tpr)
        if roc_auc < 0.5:
            fpr, tpr, _ = roc_curve(yv, -y_prob)
            roc_auc = auc(fpr, tpr)
        return fpr, tpr, roc_auc
    except Exception:
        return None

# ============================================================================
# Load & normalise
# ============================================================================
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

# ============================================================================
# Load model
# ============================================================================
ee_model = joblib.load('./models/unified_LogisticRegression.pkl')
with open('./results/core_features.txt', 'r') as f:
    lines = f.readlines()
core_features = [l.strip().split('. ', 1)[1] for l in lines[2:] if '. ' in l]
n_features = len(core_features)

# ============================================================================
# Study config  (study_7 excluded)
# ============================================================================
study_biomarkers = {
    'study_1':  ['AB42', 'pTau', 'tTau'],
    'study_4':  ['AB42', 'AB40', 'tTau', 'NEFL', 'YKL40',
                 'pTau181', 'pTau217', 'pTau231'],
    'study_9':  ['AB42', 'tTau', 'pTau181'],
    'study_11': ['AB42', 'AB40', 'pTau', 'tTau'],
}
studies_to_analyze = list(study_biomarkers.keys())

# ============================================================================
# Per-study ROC figures
# ============================================================================
for study in studies_to_analyze:
    df_s = df_normalized[df_normalized['Study'] == study].copy()
    if len(df_s) == 0 or df_s['Label'].nunique() < 2:
        continue
    
    y_true = df_s['Label'].values
    n_eod  = int(y_true.sum())
    n_cn   = int((y_true == 0).sum())

    # EODstage score (LogisticRegression)
    model_score = None
    try:
        X_s = df_s[core_features].fillna(0).values
        model_score = ee_model.predict_proba(X_s)[:, 1]
    except Exception:
        pass

    roc_data = {}   # name -> {fpr, tpr, auc}

    # EODstage alone
    if model_score is not None:
        fpr, tpr, r = direction_aware_roc(y_true, model_score)
        roc_data[f'EODstage ({n_features} proteins)'] = {
            'fpr': fpr, 'tpr': tpr, 'auc': r,
            'color': MODEL_COLOR, 'lw': 2.5, 'ls': '-', 'zorder': 10
        }

    # Biomarker alone + biomarker + EODstage
    for bio in study_biomarkers.get(study, []):
        if bio not in df_s.columns:
            continue
        bio_vals = df_s[bio].values
        valid = ~np.isnan(bio_vals)
        if valid.sum() < 20 or len(np.unique(y_true[valid])) < 2:
            continue

        fpr, tpr, r = direction_aware_roc(y_true[valid], bio_vals[valid])
        roc_data[bio] = {
            'fpr': fpr, 'tpr': tpr, 'auc': r,
            'color': FIXED_COLORS.get(bio, '#888888'),
            'lw': 1.6, 'ls': '-', 'zorder': 5
        }

        # Combined with EODstage
        if model_score is not None:
            res = z_combine_roc(y_true, bio_vals, model_score)
            if res is not None:
                fpr_c, tpr_c, r_c = res
                roc_data[f'{bio}+EODstage'] = {
                    'fpr': fpr_c, 'tpr': tpr_c, 'auc': r_c,
                    'color': COMBO_COLORS.get(bio, '#333333'),
                    'lw': 1.6, 'ls': '--', 'zorder': 6
                }

    if not roc_data:
        continue

    # Sort by AUC descending
    sorted_items = sorted(roc_data.items(),
                          key=lambda x: x[1]['auc'], reverse=True)
    n_curves = len(sorted_items)

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    ax.plot([0, 1], [0, 1], 'k--', lw=0.7, alpha=0.35)

    for name, data in sorted_items:
            ax.plot(data['fpr'], data['tpr'], 
                color=data['color'], lw=data['lw'],
                linestyle=data['ls'], zorder=data['zorder'],
                label=f"{name} (AUC\u2009=\u2009{data['auc']:.3f})")

    ax.set_xlabel('1 \u2212 Specificity (False Positive Rate)',
                  fontsize=11, fontweight='bold')
    ax.set_ylabel('Sensitivity (True Positive Rate)',
                  fontsize=11, fontweight='bold')
    ax.set_title(
        f'{study}: ROC Curves \u2014 EOD vs CN\n'
        f'n\u2009=\u2009{len(df_s)}  (EOD\u2009=\u2009{n_eod}, CN\u2009=\u2009{n_cn})\n'
        f'Solid\u2009=\u2009biomarker alone \u2009|\u2009 Dashed\u2009=\u2009biomarker\u2009+\u2009EODstage',
        fontsize=11, fontweight='bold', pad=10)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
        ax.set_aspect('equal')
    ax.grid(True, alpha=0.2, linestyle='--', lw=0.5)
    for sp in ax.spines.values():
        sp.set_linewidth(0.7)

    leg_loc  = 'lower right' if n_curves <= 8 else 'best'
    leg_size = 7.5 if n_curves > 9 else 8.5
    leg = ax.legend(loc=leg_loc, fontsize=leg_size, frameon=True,
                    fancybox=False, edgecolor='#333333',
                    framealpha=0.95, handlelength=2.0, labelspacing=0.3)
    leg.get_frame().set_linewidth(0.5)
        
        plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, f'{study}_roc_comparison.png'),
                dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(OUTPUT_DIR, f'{study}_roc_comparison.pdf'),
                bbox_inches='tight')
    plt.close(fig)
