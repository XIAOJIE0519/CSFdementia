#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
EOD Classification Pipeline — Evaluation Phase (Steps 9-10)
Depends on easyensemble.py outputs in ./results/ and ./models/
Outputs: lodo_cv_results.csv, external_test_results.csv, summary_comparison.csv
"""

import os
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, AdaBoostClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from catboost import CatBoostClassifier
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    confusion_matrix, roc_curve, auc, precision_recall_curve, roc_auc_score
)
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from imblearn.ensemble import EasyEnsembleClassifier
import joblib
import warnings

warnings.filterwarnings("ignore")
RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

RESULTS_DIR   = './results'
MODEL_SAVE_DIR = './models'
os.makedirs(RESULTS_DIR, exist_ok=True)

print("=" * 80)
print("EOD Classification Pipeline — Evaluation Phase")
print("=" * 80)


# ============================================================================
# CatBoostWrapper (local copy, identical to easyensemble.py)
# ============================================================================
class CatBoostWrapper(BaseEstimator, ClassifierMixin):
    def __init__(self, iterations=100, learning_rate=0.1, depth=5,
                 random_state=42, verbose=0, l2_leaf_reg=3):
        self.iterations = iterations; self.learning_rate = learning_rate
        self.depth = depth; self.random_state = random_state
        self.verbose = verbose; self.l2_leaf_reg = l2_leaf_reg

    def __sklearn_tags__(self):
        from sklearn.utils._tags import ClassifierTags
        tags = super().__sklearn_tags__()
        tags.classifier_tags = ClassifierTags()
        tags.estimator_type = 'classifier'
        return tags

    def fit(self, X, y):
        self.model_ = CatBoostClassifier(
            iterations=self.iterations, learning_rate=self.learning_rate,
            depth=self.depth, random_seed=self.random_state,
            verbose=self.verbose, l2_leaf_reg=self.l2_leaf_reg)
        self.model_.fit(X, y)
        self.classes_ = np.unique(y)
        return self

    def predict(self, X):       return self.model_.predict(X).ravel()
    def predict_proba(self, X): return self.model_.predict_proba(X)


# ============================================================================
# Load artefacts from training phase
# ============================================================================
print("\nLoading training artefacts...")
unified_models     = joblib.load(os.path.join(RESULTS_DIR, 'unified_models.pkl'))
optimal_thresholds = joblib.load(os.path.join(RESULTS_DIR, 'optimal_thresholds.pkl'))
core_features      = joblib.load(os.path.join(RESULTS_DIR, 'core_features.pkl'))
df_test_external   = joblib.load(os.path.join(RESULTS_DIR, 'df_test_external.pkl'))
test_studies       = joblib.load(os.path.join(RESULTS_DIR, 'test_studies.pkl'))
print(f"  Core features ({len(core_features)}): {core_features}")
print(f"  Test studies  ({len(test_studies)}): {test_studies}")

# Re-normalise training data for LODO CV
train_studies = ['study_1', 'study_4', 'study_6', 'study_7', 'study_9', 'study_11']
df_combined = pd.read_csv("../combine/combined_expression_matrices.csv",
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
df_norm_list = []
for study in df_combined['Study'].unique():
    ds = df_combined[df_combined['Study'] == study].copy()
    pcols = [c for c in ds.columns
             if c not in meta_cols
             and ds[c].notna().sum() > 0
             and pd.api.types.is_numeric_dtype(ds[c])]
    ds[pcols] = StandardScaler().fit_transform(ds[pcols].fillna(0))
    df_norm_list.append(ds)
df_normalized = pd.concat(df_norm_list, ignore_index=True)
df_train_all  = df_normalized[
    df_normalized['Study'].isin(train_studies)].copy()
print(f"Training data reloaded: {len(df_train_all)} samples")


# ============================================================================
# Helper functions
# ============================================================================
def find_optimal_threshold_youden(y_true, y_prob):
    fpr, tpr, thresholds = roc_curve(y_true, y_prob)
    return thresholds[np.argmax(tpr - fpr)]

def calculate_metrics_with_threshold(y_true, y_prob, threshold):
    y_pred = (y_prob >= threshold).astype(int)
    cm = confusion_matrix(y_true, y_pred)
    tn, fp, fn, tp = cm.ravel() if cm.shape == (2, 2) else (0, 0, 0, 0)
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    try:
        roc_auc = roc_auc_score(y_true, y_prob) if len(np.unique(y_true)) > 1 else 0.5
    except Exception:
        roc_auc = 0.5
    try:
        prec_c, rec_c, _ = precision_recall_curve(y_true, y_prob)
        auprc = auc(rec_c, prec_c)
    except Exception:
        auprc = 0.5
    return {'ROC_AUC': roc_auc, 'AUPRC': auprc,
            'Accuracy':    accuracy_score(y_true, y_pred),
            'Precision':   precision_score(y_true, y_pred, zero_division=0),
            'Recall':      recall_score(y_true, y_pred, zero_division=0),
            'Specificity': specificity,
            'F1':          f1_score(y_true, y_pred, zero_division=0)}

def train_easyensemble_model(X_train, y_train, base_clf, param_grid):
    gs = GridSearchCV(base_clf, param_grid, cv=5, n_jobs=1,
                      scoring='roc_auc', verbose=0)
    gs.fit(X_train, y_train)
    ee = EasyEnsembleClassifier(estimator=gs.best_estimator_,
                                n_estimators=10,
                                random_state=RANDOM_STATE, n_jobs=1)
    ee.fit(X_train, y_train)
    return ee, gs.best_params_


# Model definitions for LODO CV
base_classifiers = {
    'RandomForest':      RandomForestClassifier(random_state=RANDOM_STATE, n_jobs=1),
    'ExtraTrees':        ExtraTreesClassifier(random_state=RANDOM_STATE, n_jobs=1),
    'XGBoost':           XGBClassifier(random_state=RANDOM_STATE,
                                       eval_metric='logloss', n_jobs=1),
    'LightGBM':          LGBMClassifier(random_state=RANDOM_STATE,
                                        n_jobs=1, verbose=-1),
    'CatBoost':          CatBoostWrapper(random_state=RANDOM_STATE, verbose=0),
    'AdaBoost':          AdaBoostClassifier(random_state=RANDOM_STATE),
    'SVM':               SVC(probability=True, random_state=RANDOM_STATE),
    'MLP':               MLPClassifier(random_state=RANDOM_STATE, max_iter=1000),
    'LogisticRegression':LogisticRegression(max_iter=2000,
                                            random_state=RANDOM_STATE, n_jobs=1),
    'ElasticNet':        LogisticRegression(penalty='elasticnet', solver='saga',
                                            max_iter=2000, random_state=RANDOM_STATE,
                                            n_jobs=1),
    'LDA':               LinearDiscriminantAnalysis(),
}

param_grids = {
    'RandomForest':  {'n_estimators': [200, 300, 500], 'max_depth': [8, 12, 20, None],
                      'min_samples_leaf': [1, 2, 4], 'max_features': ['sqrt', 'log2']},
    'ExtraTrees':    {'n_estimators': [200, 300], 'max_depth': [8, 15, None],
                      'min_samples_leaf': [1, 2]},
    'XGBoost':       {'n_estimators': [100, 200, 300], 'learning_rate': [0.01, 0.05, 0.1],
                      'max_depth': [3, 5, 7], 'subsample': [0.8, 1.0],
                      'colsample_bytree': [0.8, 1.0]},
    'LightGBM':      {'n_estimators': [100, 200, 300], 'learning_rate': [0.01, 0.05, 0.1],
                      'max_depth': [3, 5, 7], 'num_leaves': [15, 31, 63],
                      'min_child_samples': [5, 10, 20]},
    'CatBoost':      {'iterations': [100, 200, 300], 'learning_rate': [0.01, 0.05, 0.1],
                      'depth': [4, 6, 8], 'l2_leaf_reg': [1, 3, 5]},
    'AdaBoost':      {'n_estimators': [100, 200, 300],
                      'learning_rate': [0.01, 0.1, 0.5, 1.0]},
    'SVM':           {'C': [0.01, 0.1, 1, 10, 100], 'kernel': ['rbf', 'linear'],
                      'gamma': ['scale', 'auto']},
    'MLP':           {'hidden_layer_sizes': [(64,), (128, 64), (256, 128, 64), (128,)],
                      'alpha': [0.0001, 0.001, 0.01],
                      'learning_rate_init': [0.001, 0.01]},
    'LogisticRegression': {'C': [0.01, 0.1, 1, 10, 100]},
    'ElasticNet':    {'C': [0.01, 0.1, 1, 10, 100],
                      'l1_ratio': [0.1, 0.3, 0.5, 0.7, 0.9]},
    'LDA':           {'solver': ['lsqr'], 'shrinkage': [None, 'auto', 0.1, 0.5]},
}

# ============================================================================
# Step 9: LODO cross-validation
# ============================================================================
print("\n" + "=" * 80)
print("Step 9: LODO Cross-Validation")
print("=" * 80)

lodo_results = []
for test_study in train_studies:
    tr_studies = [s for s in train_studies if s != test_study]
    df_tr = df_train_all[df_train_all['Study'].isin(tr_studies)].copy()
    df_te = df_train_all[df_train_all['Study'] == test_study].copy()
    X_tr  = df_tr[core_features].fillna(0).values
    y_tr  = df_tr['Label'].values
    X_te  = df_te[core_features].fillna(0).values
    y_te  = df_te['Label'].values
    print(f"\n  LODO test: {test_study} "
          f"(EOD={y_te.sum()}, CN={(y_te==0).sum()})")
    for model_name in base_classifiers:
        try:
            ee, _ = train_easyensemble_model(
                X_tr, y_tr,
                base_classifiers[model_name], param_grids[model_name])
            thr  = find_optimal_threshold_youden(
                y_tr, ee.predict_proba(X_tr)[:, 1])
            mets = calculate_metrics_with_threshold(
                y_te, ee.predict_proba(X_te)[:, 1], thr)
        except Exception as e:
            print(f"    {model_name}: skipped — {e}"); continue
        lodo_results.append({'Test_Study': test_study, 'Model': model_name,
                             'Threshold': thr, **mets})
        print(f"    {model_name}: AUC={mets['ROC_AUC']:.4f}")

pd.DataFrame(lodo_results).to_csv(
    os.path.join(RESULTS_DIR, 'lodo_cv_results.csv'), index=False)
print(f"\nSaved: lodo_cv_results.csv ({len(lodo_results)} rows)")

# ============================================================================
# Step 10: External test set evaluation
# ============================================================================
print("\n" + "=" * 80)
print("Step 10: External Test Set Evaluation")
print("=" * 80)

external_results = []
for ts in test_studies:
    df_ts = df_test_external[df_test_external['Study'] == ts].copy()
    if len(df_ts) == 0:
        continue
    y_ts = df_ts['Label'].values
    print(f"\n  {ts}: {len(df_ts)} samples "
          f"(EOD={y_ts.sum()}, CN={(y_ts==0).sum()})")
    # Align features: zero-pad missing columns
    X_ts = np.zeros((len(df_ts), len(core_features)))
    for j, f in enumerate(core_features):
        if f in df_ts.columns:
            X_ts[:, j] = df_ts[f].fillna(0).values
    for model_name, ee_model in unified_models.items():
        try:
            y_prob = ee_model.predict_proba(X_ts)[:, 1]
            thr    = optimal_thresholds[model_name]
            mets   = calculate_metrics_with_threshold(y_ts, y_prob, thr)
        except Exception as e:
            print(f"    {model_name}: skipped — {e}"); continue
        external_results.append({'Test_Study': ts, 'Model': model_name,
                                 'Threshold': thr, **mets})
        print(f"    {model_name}: AUC={mets['ROC_AUC']:.4f}")

pd.DataFrame(external_results).to_csv(
    os.path.join(RESULTS_DIR, 'external_test_results.csv'), index=False)
print(f"\nSaved: external_test_results.csv ({len(external_results)} rows)")

# ============================================================================
# Summary comparison
# ============================================================================
print("\n" + "=" * 80)
print("Summary Comparison")
print("=" * 80)

unified_results = pd.read_csv(
    os.path.join(RESULTS_DIR, 'unified_training_results.csv'))

summary = []
for model_name in base_classifiers:
    tr  = unified_results[unified_results['Model'] == model_name]
    lo  = [r for r in lodo_results     if r['Model'] == model_name]
    ex  = [r for r in external_results if r['Model'] == model_name]
    summary.append({
        'Model':                   model_name,
        'Unified_Training_AUC':    tr['ROC_AUC'].values[0]  if len(tr) else np.nan,
        'Unified_Training_F1':     tr['F1'].values[0]       if len(tr) else np.nan,
        'LODO_CV_Avg_AUC':         np.mean([r['ROC_AUC'] for r in lo]) if lo else np.nan,
        'LODO_CV_Avg_F1':          np.mean([r['F1']      for r in lo]) if lo else np.nan,
        'External_Test_Avg_AUC':   np.mean([r['ROC_AUC'] for r in ex]) if ex else np.nan,
        'External_Test_Avg_F1':    np.mean([r['F1']      for r in ex]) if ex else np.nan,
    })

pd.DataFrame(summary).to_csv(
    os.path.join(RESULTS_DIR, 'summary_comparison.csv'), index=False)

print("\n" + "=" * 80)
print("Evaluation Phase Complete!")
print("=" * 80)
print(f"  lodo_cv_results.csv       — {len(lodo_results)} rows")
print(f"  external_test_results.csv — {len(external_results)} rows")
print(f"  summary_comparison.csv    — {len(summary)} models")
