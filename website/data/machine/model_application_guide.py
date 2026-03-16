#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
EODstage Model Application Guide
Model: EasyEnsemble-LogisticRegression
  - EasyEnsembleClassifier(n_estimators=10)
  - Base: Pipeline(RandomUnderSampler → LogisticRegression(C=0.1, solver='lbfgs',
                   max_iter=2000, random_state=42))
Core features: 15 CSF proteins (see results/core_features.txt)
"""

import os
import pandas as pd
import numpy as np
import joblib
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve, auc
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("EODstage — EasyEnsemble-LogisticRegression Application Guide")
print("=" * 80)

# ============================================================================
# Model specification
# ============================================================================
print("""
Model: EasyEnsembleClassifier
  n_estimators : 10
  Base estimator pipeline:
    1. RandomUnderSampler (auto strategy, no replacement)
    2. LogisticRegression
         C            = 0.1
         solver       = lbfgs
         max_iter     = 2000
         random_state = 42
         n_jobs       = 1

Core features (15 CSF proteins, selected by SVM-RFE ∩ Lasso ∩ Boruta, FDR<0.01):
  MIF, DDAH1, ENO2, PEBP1, PAM, SPON1, SOD1, RTN4R,
  VASN, SOD2, GFRA2, CA4, CANT1, GLRX, GAS6

Optimal threshold (Youden index on training set): 0.5500
Training LODO CV AUC : 0.877 ± (per study 0.691–0.992)
External test AUC    : 0.860 (average across 7 independent cohorts)
""")

# ============================================================================
# Step-by-step application code
# ============================================================================
print("=" * 80)
print("Step-by-Step: Applying EODstage to a New Cohort")
print("=" * 80)

print("""
STEP 1 — Prepare data
----------------------
Required columns: the 15 core protein columns listed above.
Each sample is one CSF specimen. Missing values are imputed with 0
after per-cohort StandardScaler normalisation.

STEP 2 — Normalise within your cohort (critical)
-------------------------------------------------
The model was trained with per-study StandardScaler (fit on each study
independently). You MUST apply a fresh StandardScaler to your new cohort.
Do NOT reuse the training scaler.

STEP 3 — Load model and predict
---------------------------------
See complete example below.

STEP 4 — Interpret probability
--------------------------------
  ≥ 0.70  → High EOD risk
  0.40–0.70 → Borderline
  < 0.40  → Low risk (likely CN)

Use the Youden-optimal threshold (0.5500) for binary classification.
""")

# ============================================================================
# Complete application example
# ============================================================================
example = '''
import pandas as pd
import numpy as np
import joblib
from sklearn.preprocessing import StandardScaler

# --- 1. Load your new cohort ---
df_new = pd.read_csv('your_new_cohort.csv')

CORE_FEATURES = [
    'MIF', 'DDAH1', 'ENO2', 'PEBP1', 'PAM',
    'SPON1', 'SOD1', 'RTN4R', 'VASN', 'SOD2',
    'GFRA2', 'CA4', 'CANT1', 'GLRX', 'GAS6',
]

# --- 2. Normalise within your cohort ---
sc = StandardScaler()
X_new = sc.fit_transform(df_new[CORE_FEATURES].fillna(0).values)

# --- 3. Load model and threshold ---
model     = joblib.load('./models/unified_LogisticRegression.pkl')
threshold = 0.5500  # Youden-index optimal on training set

# --- 4. Predict ---
y_prob = model.predict_proba(X_new)[:, 1]   # EOD probability
y_pred = (y_prob >= threshold).astype(int)  # 0=CN, 1=EOD

# --- 5. Save results ---
df_new['EODstage_Score']      = y_prob
df_new['EODstage_Prediction'] = y_pred
df_new['Risk_Group'] = pd.cut(
    y_prob,
    bins=[0, 0.40, 0.70, 1.0],
    labels=['Low', 'Borderline', 'High'],
    include_lowest=True
)
df_new.to_csv('eodstage_predictions.csv', index=False)
print(df_new[['EODstage_Score', 'EODstage_Prediction', 'Risk_Group']].head())
'''
print("Complete Application Example:")
print(example)

# ============================================================================
# Important notes
# ============================================================================
print("=" * 80)
print("Important Notes")
print("=" * 80)
print("""
1. NORMALISATION
   Always fit a fresh StandardScaler on your cohort.
   Mixing scalers across cohorts will degrade performance.

2. MISSING PROTEINS
   If a core protein is absent in your platform, fillna(0) is applied.
   Performance may degrade proportionally to the number of missing proteins.
   Validate on a labelled subset before clinical use.

3. SAMPLE SIZE
   Recommend ≥ 30 samples for stable normalisation.
   With fewer samples, per-sample z-scores become unreliable.

4. COHORT SHIFT
   Model was trained on SomaScan CSF proteomics (studies 1,4,6,7,9,11).
   Platform differences (e.g. Olink, MS-based) may reduce performance.
   study_3_Olink achieved AUC=0.986 suggesting good cross-platform transfer
   for most proteins; study_12 (MS-based, EOND label) achieved AUC=0.642.

5. SUBTYPE COVERAGE
   Training set includes EOAD, EODLB, EOFTD, EOOD, EODSD.
   The model is an EOD-vs-CN binary classifier, not a subtype classifier.

6. THRESHOLD
   Threshold 0.5500 was optimised by Youden index on the pooled training set.
   For clinical use, recalibrate on a local labelled validation set.
""")

print("=" * 80)
print("Model files:")
print("  ./models/unified_LogisticRegression.pkl  — EasyEnsemble model")
print("  ./results/core_features.pkl              — list of 15 feature names")
print("  ./results/optimal_thresholds.csv         — Youden thresholds per model")
print("=" * 80)
