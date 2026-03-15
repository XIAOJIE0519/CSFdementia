#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
EOD Classification Pipeline — Training Phase (Steps 1-8)
Feature Selection: FDR<0.01 & Present in all 6 training studies
  -> SVM-RFE / Lasso / Boruta (Votes=3)
Outputs: trained models, feature info, per-dataset CV, between-dataset CV
"""

import os
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, AdaBoostClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.linear_model import LogisticRegression, LassoCV
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.feature_selection import RFE
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from catboost import CatBoostClassifier
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    confusion_matrix, roc_curve, auc, precision_recall_curve, roc_auc_score
)
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from imblearn.ensemble import EasyEnsembleClassifier
import joblib
import warnings

warnings.filterwarnings("ignore")
RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

MODEL_SAVE_DIR = './models'
RESULTS_DIR    = './results'
PLOTS_DIR      = './plots'
os.makedirs(MODEL_SAVE_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR,    exist_ok=True)
os.makedirs(PLOTS_DIR,      exist_ok=True)


class CatBoostWrapper(BaseEstimator, ClassifierMixin):
    """Sklearn-compatible wrapper for CatBoostClassifier."""
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


print("=" * 80)
print("EOD Classification Pipeline — Training Phase")
print("=" * 80)

# ============================================================================
# Step 1: Load data
# ============================================================================
print("\n" + "=" * 80)
print("Step 1: Loading Data")
print("=" * 80)

df_combined = pd.read_csv("../combine/combined_expression_matrices.csv",
                          low_memory=False)
print(f"Combined data shape: {df_combined.shape}")

eod_diagnoses = ['EOAD', 'EODSD', 'EOFTD', 'EOOD', 'EODLB']
df_combined['Label'] = df_combined['Diagnosis_Derived'].apply(
    lambda x: 1 if x in eod_diagnoses else 0)
df_combined = df_combined[
    df_combined['Diagnosis_Derived'].isin(eod_diagnoses + ['CN'])].copy()
print(f"After filtering: {df_combined.shape}, "
      f"EOD={df_combined['Label'].sum()}, "
      f"CN={(df_combined['Label']==0).sum()}")

train_studies = ['study_1', 'study_4', 'study_6', 'study_7', 'study_9', 'study_11']

non_protein_cols = [
    'Study', 'Batch', 'GUID', 'Age', 'Sex', 'Cognitive Score',
    'AB42', 'tTau', 'pTau', 'pTau181', 'AB42/pTau', 'AB40',
    'NEFL', 'YKL40', 'pTau217', 'pTau231', 'Diagnosis_Derived', 'Label',
]
meta_cols = [c for c in non_protein_cols if c in df_combined.columns]

# Per-study normalisation
print("\nNormalising each study independently...")
df_normalized_list = []
for study in df_combined['Study'].unique():
    ds = df_combined[df_combined['Study'] == study].copy()
    pcols = [c for c in ds.columns
             if c not in meta_cols
             and ds[c].notna().sum() > 0
             and pd.api.types.is_numeric_dtype(ds[c])]
    print(f"  {study}: {len(ds)} samples, {len(pcols)} proteins")
    ds[pcols] = StandardScaler().fit_transform(ds[pcols].fillna(0))
    df_normalized_list.append(ds)
df_combined_normalized = pd.concat(df_normalized_list, ignore_index=True)
print(f"Total samples after normalisation: {len(df_combined_normalized)}")

# study_12: independent CSV, CN + EOND only, normalised separately
print("\nLoading study_12 (CN and EOND only) ...")
df_s12_raw = pd.read_csv("F:/1a-EOD-CSF-protein/S12/study_12_final.csv")
df_s12 = df_s12_raw[df_s12_raw['Group_New'].isin(['CN', 'EOND'])].copy().reset_index(drop=True)
df_s12['Study'] = 'study_12'
df_s12['Age']   = df_s12_raw.iloc[:, 1][df_s12.index].values
df_s12['Label'] = (df_s12['Group_New'] == 'EOND').astype(int)
non_prot_s12 = {df_s12_raw.columns[0], 'Group_New', 'Study', 'Age', 'Label'}
protein_cols_s12 = [c for c in df_s12.columns
                    if c not in non_prot_s12 and df_s12[c].notna().sum() > 0]
df_s12[protein_cols_s12] = StandardScaler().fit_transform(
    df_s12[protein_cols_s12].fillna(0))
print(f"  study_12: {len(df_s12)} samples "
      f"(EOND={df_s12['Label'].sum()}, CN={(df_s12['Label']==0).sum()}), "
      f"{len(protein_cols_s12)} proteins")

# Train / test split
df_train_all = df_combined_normalized[
    df_combined_normalized['Study'].isin(train_studies)].copy()
df_test_external = pd.concat([
    df_combined_normalized[
        ~df_combined_normalized['Study'].isin(train_studies)],
    df_s12
], ignore_index=True)
test_studies = sorted(df_test_external['Study'].unique().tolist())
print(f"\nTraining studies ({len(train_studies)}): {train_studies}")
print(f"Test studies     ({len(test_studies)}): {test_studies}")
print(f"Training : {len(df_train_all)} samples "
      f"(EOD={df_train_all['Label'].sum()}, CN={(df_train_all['Label']==0).sum()})")
print(f"Test     : {len(df_test_external)} samples "
      f"(EOD={df_test_external['Label'].sum()}, CN={(df_test_external['Label']==0).sum()})")

# Persist test set for evaluate.py
joblib.dump(df_test_external, os.path.join(RESULTS_DIR, 'df_test_external.pkl'))
joblib.dump(test_studies,     os.path.join(RESULTS_DIR, 'test_studies.pkl'))

# ============================================================================
# Step 2: Proteins common to all 6 training studies
# ============================================================================
print("\n" + "=" * 80)
print("Step 2: Finding common proteins in all 6 training studies")
print("=" * 80)

protein_cols_all = [c for c in df_train_all.columns if c not in meta_cols]
proteins_in_all_6 = [
    p for p in protein_cols_all
    if all(df_train_all[df_train_all['Study'] == s][p].notna().sum() > 0
           for s in train_studies)
]
print(f"Proteins present in all 6 training studies: {len(proteins_in_all_6)}")

# ============================================================================
# Step 3: Feature selection — FDR<0.01, SVM-RFE / Lasso / Boruta
# ============================================================================
print("\n" + "=" * 80)
print("Step 3: Feature Selection (FDR<0.01 & SVM-RFE/Lasso/Boruta, Votes=3)")
print("=" * 80)

df_diff = pd.read_csv("../meta/EOD_vs_CN.csv")
sig_proteins = df_diff[df_diff['FDR_BH_Stratified'] < 0.01]['Protein'].tolist()
print(f"Significant proteins (FDR<0.01): {len(sig_proteins)}")
initial_features = [p for p in sig_proteins if p in proteins_in_all_6]
print(f"Candidate features (in all 6 & FDR<0.01): {len(initial_features)}")

X_fs = df_train_all[initial_features].fillna(0).values
y_fs = df_train_all['Label'].values
X_fs_scaled = StandardScaler().fit_transform(X_fs)

print("Running SVM-RFE...")
rfe = RFE(estimator=SVC(kernel='linear', random_state=RANDOM_STATE),
          n_features_to_select=max(10, len(initial_features) // 3), step=0.1)
rfe.fit(X_fs_scaled, y_fs)

print("Running Lasso...")
lasso = LassoCV(cv=5, random_state=RANDOM_STATE, max_iter=10000, n_jobs=-1)
lasso.fit(X_fs_scaled, y_fs)
lasso_coef = np.abs(lasso.coef_)

print("Running Boruta (RF importance)...")
rf_bor = RandomForestClassifier(n_estimators=100, random_state=RANDOM_STATE,
                                n_jobs=-1)
rf_bor.fit(X_fs_scaled, y_fs)
rf_imp = rf_bor.feature_importances_

lasso_thr = (np.percentile(lasso_coef[lasso_coef > 0], 33)
             if np.any(lasso_coef > 0) else 0.0)
rf_thr = np.percentile(rf_imp, 67)

feature_details = []
for i, feat in enumerate(initial_features):
    row = df_diff[df_diff['Protein'] == feat]
    feature_details.append({
        'Feature':           feat,
        'FDR':               row['FDR_BH_Stratified'].values[0] if len(row) else np.nan,
        'Weighted_Effect':   row['Weighted_Effect'].values[0]    if len(row) else np.nan,
        'SVM_RFE_Ranking':   int(rfe.ranking_[i]),
        'SVM_RFE_Selected':  bool(rfe.support_[i]),
        'Lasso_Coefficient': float(lasso_coef[i]),
        'Boruta_Importance': float(rf_imp[i]),
        'Lasso_Selected':    bool(lasso_coef[i] > lasso_thr),
        'Boruta_Selected':   bool(rf_imp[i] > rf_thr),
    })

df_fd = pd.DataFrame(feature_details)
df_fd['Votes'] = (df_fd['SVM_RFE_Selected'].astype(int)
                  + df_fd['Lasso_Selected'].astype(int)
                  + df_fd['Boruta_Selected'].astype(int))
df_fd['Core_Feature'] = df_fd['Votes'] == 3
df_fd.to_csv(os.path.join(RESULTS_DIR, 'feature_selection_detailed.csv'), index=False)

core_features = df_fd[df_fd['Core_Feature']]['Feature'].tolist()
print(f"\nCore features (Votes=3): {len(core_features)}")
for i, f in enumerate(core_features, 1):
    print(f"  {i}. {f}")

with open(os.path.join(RESULTS_DIR, 'core_features.txt'), 'w') as fh:
    fh.write(f"Total core features: {len(core_features)}\n")
    fh.write("=" * 50 + "\n")
    for i, f in enumerate(core_features, 1):
        fh.write(f"{i}. {f}\n")

joblib.dump(core_features, os.path.join(RESULTS_DIR, 'core_features.pkl'))

# ============================================================================
# Step 4: Model definitions
# ============================================================================
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
    'RandomForest':  {'n_estimators': [200, 300, 500],
                      'max_depth': [8, 12, 20, None],
                      'min_samples_leaf': [1, 2, 4],
                      'max_features': ['sqrt', 'log2']},
    'ExtraTrees':    {'n_estimators': [200, 300],
                      'max_depth': [8, 15, None],
                      'min_samples_leaf': [1, 2]},
    'XGBoost':       {'n_estimators': [100, 200, 300],
                      'learning_rate': [0.01, 0.05, 0.1],
                      'max_depth': [3, 5, 7],
                      'subsample': [0.8, 1.0],
                      'colsample_bytree': [0.8, 1.0]},
    'LightGBM':      {'n_estimators': [100, 200, 300],
                      'learning_rate': [0.01, 0.05, 0.1],
                      'max_depth': [3, 5, 7],
                      'num_leaves': [15, 31, 63],
                      'min_child_samples': [5, 10, 20]},
    'CatBoost':      {'iterations': [100, 200, 300],
                      'learning_rate': [0.01, 0.05, 0.1],
                      'depth': [4, 6, 8],
                      'l2_leaf_reg': [1, 3, 5]},
    'AdaBoost':      {'n_estimators': [100, 200, 300],
                      'learning_rate': [0.01, 0.1, 0.5, 1.0]},
    'SVM':           {'C': [0.01, 0.1, 1, 10, 100],
                      'kernel': ['rbf', 'linear'],
                      'gamma': ['scale', 'auto']},
    'MLP':           {'hidden_layer_sizes': [(64,), (128, 64),
                                             (256, 128, 64), (128,)],
                      'alpha': [0.0001, 0.001, 0.01],
                      'learning_rate_init': [0.001, 0.01]},
    'LogisticRegression': {'C': [0.01, 0.1, 1, 10, 100]},
    'ElasticNet':    {'C': [0.01, 0.1, 1, 10, 100],
                      'l1_ratio': [0.1, 0.3, 0.5, 0.7, 0.9]},
    'LDA':           {'solver': ['lsqr'],
                      'shrinkage': [None, 'auto', 0.1, 0.5]},
}

# ============================================================================
# Step 5: Helper functions
# ============================================================================
def find_optimal_threshold_youden(y_true, y_prob):
    fpr, tpr, thresholds = roc_curve(y_true, y_prob)
    idx = np.argmax(tpr - fpr)
    return thresholds[idx]

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
            'Accuracy': accuracy_score(y_true, y_pred),
            'Precision': precision_score(y_true, y_pred, zero_division=0),
            'Recall': recall_score(y_true, y_pred, zero_division=0),
            'Specificity': specificity,
            'F1': f1_score(y_true, y_pred, zero_division=0)}

def train_easyensemble_model(X_train, y_train, base_clf, param_grid):
    gs = GridSearchCV(base_clf, param_grid, cv=5, n_jobs=1,
                      scoring='roc_auc', verbose=0)
    gs.fit(X_train, y_train)
    ee = EasyEnsembleClassifier(estimator=gs.best_estimator_,
                                n_estimators=10,
                                random_state=RANDOM_STATE, n_jobs=1)
    ee.fit(X_train, y_train)
    return ee, gs.best_params_

# ============================================================================
# Step 6: Train unified models on all 6 training studies
# ============================================================================
print("\n" + "=" * 80)
print("Step 6: Training Unified Models (studies 1,4,6,7,9,11)")
print("=" * 80)

X_train_all = df_train_all[core_features].fillna(0).values
y_train_all  = df_train_all['Label'].values

unified_models    = {}
optimal_thresholds = {}
unified_results   = []

for model_name in base_classifiers:
    print(f"\nTraining {model_name}...")
    ee, best_params = train_easyensemble_model(
        X_train_all, y_train_all,
        base_classifiers[model_name], param_grids[model_name])
    y_prob = ee.predict_proba(X_train_all)[:, 1]
    thr    = find_optimal_threshold_youden(y_train_all, y_prob)
    mets   = calculate_metrics_with_threshold(y_train_all, y_prob, thr)
    print(f"  Optimal threshold (Youden): {thr:.4f}  Train AUC: {mets['ROC_AUC']:.4f}")
    unified_models[model_name]     = ee
    optimal_thresholds[model_name] = thr
    unified_results.append({'Model': model_name, 'Optimal_Threshold': thr,
                             **mets, 'Best_Params': str(best_params)})
    joblib.dump(ee, os.path.join(MODEL_SAVE_DIR, f'unified_{model_name}.pkl'))

pd.DataFrame(unified_results).to_csv(
    os.path.join(RESULTS_DIR, 'unified_training_results.csv'), index=False)
pd.DataFrame([{'Model': k, 'Optimal_Threshold': v}
              for k, v in optimal_thresholds.items()]).to_csv(
    os.path.join(RESULTS_DIR, 'optimal_thresholds.csv'), index=False)
joblib.dump(unified_models,     os.path.join(RESULTS_DIR, 'unified_models.pkl'))
joblib.dump(optimal_thresholds, os.path.join(RESULTS_DIR, 'optimal_thresholds.pkl'))

# ============================================================================
# Step 7: Per-dataset cross-validation (within each training study)
# ============================================================================
print("\n" + "=" * 80)
print("Step 7: Per-Dataset Cross-Validation")
print("=" * 80)

per_dataset_results = []
for study in train_studies:
    df_s = df_train_all[df_train_all['Study'] == study].copy()
    if len(df_s) < 20:
        continue
    X_s = df_s[core_features].fillna(0).values
    y_s = df_s['Label'].values
    min_cls = np.bincount(y_s).min()
    n_splits = min(5, min_cls)
    if n_splits < 2:
        print(f"  {study}: skipped (minority class too small: {min_cls})")
        continue
    print(f"\n  {study}: {len(df_s)} samples, {n_splits}-fold CV")
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True,
                          random_state=RANDOM_STATE)
    for model_name in base_classifiers:
        fold_aucs = []
        for tr_idx, val_idx in skf.split(X_s, y_s):
            if (len(np.unique(y_s[tr_idx])) < 2
                    or len(np.unique(y_s[val_idx])) < 2):
                continue
            try:
                ee, _ = train_easyensemble_model(
                    X_s[tr_idx], y_s[tr_idx],
                    base_classifiers[model_name], param_grids[model_name])
                fold_aucs.append(
                    roc_auc_score(y_s[val_idx],
                                  ee.predict_proba(X_s[val_idx])[:, 1]))
            except Exception as e:
                print(f"    [{model_name}] fold skipped: {e}")
        if not fold_aucs:
            continue
        per_dataset_results.append({
            'Study': study, 'Model': model_name,
            'Mean_AUC': np.mean(fold_aucs),
            'Std_AUC':  np.std(fold_aucs)})
        print(f"    {model_name}: AUC = {np.mean(fold_aucs):.4f} "
              f"± {np.std(fold_aucs):.4f}")

pd.DataFrame(per_dataset_results).to_csv(
    os.path.join(RESULTS_DIR, 'per_dataset_cv_results.csv'), index=False)

# ============================================================================
# Step 8: Between-dataset cross-validation
# ============================================================================
print("\n" + "=" * 80)
print("Step 8: Between-Dataset Cross-Validation")
print("=" * 80)

between_dataset_results = []
for test_study in train_studies:
    tr_studies = [s for s in train_studies if s != test_study]
    df_tr = df_train_all[df_train_all['Study'].isin(tr_studies)].copy()
    df_te = df_train_all[df_train_all['Study'] == test_study].copy()
    X_tr = df_tr[core_features].fillna(0).values
    y_tr = df_tr['Label'].values
    X_te = df_te[core_features].fillna(0).values
    y_te = df_te['Label'].values
    print(f"\n  Train: {tr_studies}  |  Test: {test_study}")
    for model_name in base_classifiers:
        try:
            ee, _ = train_easyensemble_model(
                X_tr, y_tr,
                base_classifiers[model_name], param_grids[model_name])
            test_auc = roc_auc_score(y_te,
                                     ee.predict_proba(X_te)[:, 1])
        except Exception as e:
            print(f"    {model_name}: skipped — {e}")
            continue
        between_dataset_results.append({
            'Train_Studies': ','.join(tr_studies),
            'Test_Study':    test_study,
            'Model':         model_name,
            'Test_AUC':      test_auc})
        print(f"    {model_name}: AUC = {test_auc:.4f}")

pd.DataFrame(between_dataset_results).to_csv(
    os.path.join(RESULTS_DIR, 'between_dataset_cv_results.csv'), index=False)

print("\n" + "=" * 80)
print("Training Phase Complete — run evaluate.py for LODO CV & external test")
print("=" * 80)
