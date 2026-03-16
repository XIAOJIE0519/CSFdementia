#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SHAP Analysis — EasyEnsemble-LogisticRegression

Strategy:
  1. Load EasyEnsemble-LogisticRegression
  2. Extract averaged logistic regression coefficients across 10 sub-estimators
  3. Build a single LogisticRegression surrogate with averaged coefs
  4. Use shap.LinearExplainer for exact, instant SHAP values

Outputs to ./shap/:
  summary_beeswarm.png/pdf
  summary_bar.png/pdf
  shap_importance.csv
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import joblib
import shap
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV

warnings.filterwarnings('ignore')
np.random.seed(42)

OUTPUT_DIR = './shap'
os.makedirs(OUTPUT_DIR, exist_ok=True)

rcParams.update({
    'font.family':      'serif',
    'font.serif':       ['Times New Roman'],
    'mathtext.fontset': 'stix',
    'font.size':        11,
    'axes.linewidth':   0.8,
    'figure.dpi':       300,
    'pdf.fonttype':     42,
    'ps.fonttype':      42,
})

print('=' * 70)
print('SHAP Analysis — EasyEnsemble-LogisticRegression')
print('=' * 70)

# ============================================================================
# 1. Load & normalise training data  (mirrors easyensemble.py)
# ============================================================================
print('\n[1] Loading and normalising training data ...')

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

train_studies = ['study_1', 'study_4', 'study_6', 'study_7', 'study_9', 'study_11']
norm_list = []
for study in df_combined['Study'].unique():
    ds = df_combined[df_combined['Study'] == study].copy()
    pcols = [c for c in ds.columns
             if c not in meta_cols
             and ds[c].notna().sum() > 0
             and pd.api.types.is_numeric_dtype(ds[c])]
    sc = StandardScaler()
    ds[pcols] = sc.fit_transform(ds[pcols].fillna(0))
    norm_list.append(ds)
df_normalized = pd.concat(norm_list, ignore_index=True)
df_train = df_normalized[df_normalized['Study'].isin(train_studies)].copy()
print(f'   Training samples: {len(df_train)}  '
      f'(EOD={df_train["Label"].sum()}, '
      f'CN={(df_train["Label"]==0).sum()})')

# ============================================================================
# 2. Load core features & EasyEnsemble-LogisticRegression
# ============================================================================
print('\n[2] Loading model and core features ...')

core_features = joblib.load('./results/core_features.pkl')
n_feat = len(core_features)
print(f'   Core features ({n_feat}): {core_features}')

ee_model = joblib.load('./models/unified_LogisticRegression.pkl')
print(f'   Model: {type(ee_model).__name__}, n_estimators={ee_model.n_estimators}')

X_train = df_train[core_features].fillna(0).values
y_train = df_train['Label'].values

# ============================================================================
# 3. Average LogisticRegression coefficients across 10 sub-estimators
#    Each sub-estimator is a Pipeline(RandomUnderSampler, LogisticRegression)
# ============================================================================
print('\n[3] Averaging LogisticRegression coefficients across sub-estimators ...')

coef_list = []
intercept_list = []
for pipe in ee_model.estimators_:
    lr = pipe.named_steps['classifier']
    coef_list.append(lr.coef_[0])
    intercept_list.append(lr.intercept_[0])

coef_mean      = np.mean(coef_list, axis=0)   # shape (n_features,)
intercept_mean = np.mean(intercept_list)
print(f'   Averaged {len(coef_list)} sub-estimator coefficients.')

# Build surrogate LogisticRegression with averaged parameters
from sklearn.linear_model import LogisticRegression as LR
surrogate = LR(C=0.1, max_iter=2000, random_state=42)
# Fit on full training data to set classes_, then override coef/intercept
surrogate.fit(X_train, y_train)
surrogate.coef_[0]    = coef_mean
surrogate.intercept_[0] = intercept_mean

# Verify surrogate quality
from sklearn.metrics import roc_auc_score
y_prob_ee  = ee_model.predict_proba(X_train)[:, 1]
y_prob_sur = surrogate.predict_proba(X_train)[:, 1]
auc_ee  = roc_auc_score(y_train, y_prob_ee)
auc_sur = roc_auc_score(y_train, y_prob_sur)
print(f'   EasyEnsemble-LogReg AUC on train: {auc_ee:.4f}')
print(f'   Surrogate (avg coef) AUC on train: {auc_sur:.4f}')

# ============================================================================
# 4. SHAP via LinearExplainer (exact, instant)
# ============================================================================
print('\n[4] Computing SHAP values via LinearExplainer ...')

explainer   = shap.LinearExplainer(surrogate, X_train,
                                    feature_perturbation='correlation_dependent')
shap_values = explainer.shap_values(X_train)

if isinstance(shap_values, list):
    shap_values = shap_values[1] if len(shap_values) > 1 else shap_values[0]

base_value = float(np.atleast_1d(explainer.expected_value)[0])
print(f'   SHAP matrix shape: {shap_values.shape}')
print(f'   Base value:        {base_value:.4f}')

explanation = shap.Explanation(
    values=shap_values,
    base_values=np.full(len(X_train), base_value),
    data=X_train,
    feature_names=core_features,
)

mean_abs_shap  = np.abs(shap_values).mean(axis=0)
feat_order     = np.argsort(mean_abs_shap)[::-1]
ordered_features = [core_features[i] for i in feat_order]

# ============================================================================
# 5. NC-style SHAP bar + LOWESS + donut figure
# ============================================================================
print('\n[5] Plotting NC-style SHAP bar + LOWESS + donut ...')

from matplotlib.patches import FancyArrowPatch, Wedge, Patch
from matplotlib.lines import Line2D
from scipy.stats import pearsonr

# ---- Feature biological categories (15 features) ----
# 0=Neuroinflammation, 1=Oxidative stress, 2=Neuronal injury,
# 3=Synaptic/axonal signaling, 4=Vascular/extracellular
# feat_sorted order (high→low |SHAP|) will be applied at runtime
# Fixed mapping by feature name:
_cat_map = {
    'MIF':   0, 'DDAH1': 0, 'GAS6':  0,
    'SOD1':  1, 'SOD2':  1, 'GLRX':  1,
    'ENO2':  2, 'PEBP1': 2,
    'SPON1': 3, 'RTN4R': 3, 'GFRA2': 3, 'PAM': 3,
    'VASN':  4, 'CA4':   4, 'CANT1': 4,
}
Cluster_index = [_cat_map[f] for f in core_features]
group_name    = ['Neuroinflammation', 'Oxidative stress',
                 'Neuronal injury',   'Synaptic/axonal signaling',
                 'Vascular/extracellular']
Color_Index   = ['#f67c7f', '#fcdf8e', '#90dbcd', '#7fc1ce', '#b5a0d4']

# Sort features by mean |SHAP| descending
feat_sorted   = [core_features[i] for i in feat_order]   # high→low
shap_sorted   = shap_values[:, feat_order]                # (n, 15) high→low
mean_abs_sorted = mean_abs_shap[feat_order]
cluster_sorted  = [Cluster_index[i] for i in feat_order]

top_k   = len(core_features)   # show all 15
curve_span_ratio = 1 / 3

# ---- Layout: left=bar, right=lowess, top=donut ----
fig = plt.figure(figsize=(5.5, 7.5))
gs  = fig.add_gridspec(
    2, 2,
    width_ratios=[3, 1.6],
    height_ratios=[1, 5],
    hspace=0.04, wspace=0.06,
    left=0.22, right=0.97, top=0.93, bottom=0.08
)
ax_donut  = fig.add_subplot(gs[0, 0])   # top-left: donut
ax_legend = fig.add_subplot(gs[0, 1])   # top-right: legend
ax_bar    = fig.add_subplot(gs[1, 0])   # bottom-left: bar
ax_curve  = fig.add_subplot(gs[1, 1])   # bottom-right: lowess

y_pos = np.arange(top_k)[::-1]          # 0 at bottom, top_k-1 at top

# ---- (A) Bar chart ----
for idx, (feat, yi, ci) in enumerate(
        zip(feat_sorted[:top_k], y_pos, cluster_sorted[:top_k])):
    val   = mean_abs_sorted[idx]
    col   = Color_Index[ci]
    # grey background bar
    ax_bar.barh(yi, mean_abs_sorted[0] * 1.05,
                color='#eeeeee', height=0.72, zorder=1)
    # coloured foreground bar
    ax_bar.barh(yi, val, color=col, height=0.72,
                alpha=0.85, zorder=2)
    # dot on bar tip
    ax_bar.scatter(val, yi, color=col, s=38, zorder=3,
                   edgecolors='white', linewidths=0.5)

ax_bar.set_yticks(y_pos)
ax_bar.set_yticklabels(feat_sorted[:top_k], fontsize=8.5,
                        fontfamily='serif')
ax_bar.set_xlabel('Mean |SHAP value|', fontsize=8.5, fontfamily='serif')
ax_bar.set_xlim(0, mean_abs_sorted[0] * 1.18)
ax_bar.set_ylim(-0.6, top_k - 0.4)
ax_bar.spines['top'].set_visible(False)
ax_bar.spines['right'].set_visible(False)
ax_bar.tick_params(axis='x', labelsize=7.5)
ax_bar.tick_params(axis='y', length=0)

# vertical dashed zero line
ax_bar.axvline(0, color='#aaaaaa', lw=0.6, ls='--')

# ---- (B) LOWESS scatter + curve ----
from statsmodels.nonparametric.smoothers_lowess import lowess as sm_lowess

ax_curve.set_ylim(-0.6, top_k - 0.4)
ax_curve.set_yticks([])
ax_curve.spines['top'].set_visible(False)
ax_curve.spines['right'].set_visible(False)
ax_curve.spines['left'].set_visible(False)

# x-axis: SHAP value range
shap_all = shap_sorted[:, :top_k]
xmin = shap_all.min() * 1.1
xmax = shap_all.max() * 1.1
ax_curve.set_xlim(xmin, xmax)
ax_curve.axvline(0, color='#aaaaaa', lw=0.6, ls='--')
ax_curve.set_xlabel('SHAP value', fontsize=8.5, fontfamily='serif')

# shared y-axis link
ax_curve.set_ylim(ax_bar.get_ylim())

for idx in range(top_k):
    yi   = y_pos[idx]
    sv   = shap_sorted[:, idx]           # shap values for this feature
    fv   = X_train[:, feat_order[idx]]   # feature values
    ci   = cluster_sorted[idx]
    col  = Color_Index[ci]

    # scatter (small, semi-transparent)
    sc = ax_curve.scatter(
        sv, np.full(len(sv), yi) + np.random.uniform(-0.22, 0.22, len(sv)),
        c=fv, cmap='RdYlBu_r', s=4, alpha=0.35, linewidths=0,
        vmin=np.percentile(fv, 5), vmax=np.percentile(fv, 95)
    )

    # LOWESS curve
    if len(sv) > 10:
        try:
            sorted_idx = np.argsort(sv)
            lw_out = sm_lowess(yi + np.zeros(len(sv)), sv,
                               frac=0.4, return_sorted=True)
            ax_curve.plot(lw_out[:, 0], lw_out[:, 1],
                          color=col, lw=1.2, alpha=0.85, zorder=4)
        except Exception:
            pass

# black connector line from bar tip to curve
for idx in range(top_k):
    yi    = y_pos[idx]
    x_bar = mean_abs_sorted[idx]
    x_end = mean_abs_sorted[0] * curve_span_ratio
    # draw in bar axes as extension
    ax_bar.annotate('', xy=(mean_abs_sorted[0] * 1.05, yi),
                    xytext=(x_bar, yi),
                    arrowprops=dict(arrowstyle='-', color='#888888',
                                   lw=0.5, linestyle='dotted'),
                    annotation_clip=False)

# ---- (C) Donut chart ----
cluster_counts = [cluster_sorted[:top_k].count(i)
                  for i in range(len(group_name))]
wedge_props = dict(width=0.42, edgecolor='white', linewidth=1.5)
ax_donut.pie(
    cluster_counts,
    colors=Color_Index,
    wedgeprops=wedge_props,
    startangle=90,
    counterclock=False,
)
# centre text
ax_donut.text(0, 0, f'n={top_k}\nfeatures',
              ha='center', va='center', fontsize=7,
              fontfamily='serif', color='#444444')
ax_donut.set_aspect('equal')
ax_donut.axis('off')

# ---- (D) Legend ----
ax_legend.axis('off')
legend_patches = [
    Patch(facecolor=Color_Index[i], label=group_name[i],
          edgecolor='white', linewidth=0.5)
    for i in range(len(group_name))
]
ax_legend.legend(
    handles=legend_patches,
    loc='center left',
    fontsize=7.5,
    frameon=False,
    handlelength=1.2,
    handletextpad=0.5,
    labelspacing=0.6,
)

# ---- Title ----
fig.suptitle(
    'EODstage — CSF Protein SHAP Importance\n'
    '(EasyEnsemble-LogisticRegression, training cohort)',
    fontsize=9.5, fontfamily='serif', fontweight='bold', y=0.975
)

plt.savefig(os.path.join(OUTPUT_DIR, 'summary_bar.png'),
            dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(OUTPUT_DIR, 'summary_bar.pdf'),
            bbox_inches='tight')
plt.close(fig)
print('   Saved: summary_bar.png/pdf  (NC-style bar+lowess+donut)')

# ---- keep beeswarm as backup ----
print('\n[6] Plotting beeswarm (backup) ...')
fig_bsw, ax_bsw = plt.subplots(figsize=(7, 5))
plt.sca(ax_bsw)
shap.summary_plot(
    shap_values, X_train,
    feature_names=core_features,
    max_display=top_k,
    show=False, plot_type='dot',
    color_bar=True, alpha=0.5, plot_size=None,
)
ax_bsw.set_xlabel('SHAP value', fontsize=10)
ax_bsw.axvline(0, color='#555', lw=0.6, ls='--', alpha=0.5)
for sp in ax_bsw.spines.values(): sp.set_linewidth(0.6)
plt.tight_layout()
fig_bsw.savefig(os.path.join(OUTPUT_DIR, 'summary_beeswarm.png'),
                dpi=300, bbox_inches='tight')
fig_bsw.savefig(os.path.join(OUTPUT_DIR, 'summary_beeswarm.pdf'),
                bbox_inches='tight')
plt.close(fig_bsw)
print('   Saved: summary_beeswarm.png/pdf')

# ============================================================================
# 7. Save importance CSV
# ============================================================================
mean_shap_signed = shap_values.mean(axis=0)
df_imp = pd.DataFrame({
    'Feature':        core_features,
    'Mean_Abs_SHAP':  mean_abs_shap,
    'Mean_SHAP':      mean_shap_signed,
    'Rank':           pd.Series(mean_abs_shap).rank(ascending=False).values.astype(int),
})
df_imp = df_imp.sort_values('Mean_Abs_SHAP', ascending=False)
df_imp.to_csv(os.path.join(OUTPUT_DIR, 'shap_importance.csv'), index=False)
print('   Saved: shap_importance.csv')

print('\n' + '=' * 70)
print('SHAP analysis complete!')
print(f'Outputs: {OUTPUT_DIR}/')
for f in ['summary_bar.png', 'summary_bar.pdf',
          'summary_beeswarm.png', 'summary_beeswarm.pdf',
          'shap_importance.csv']:
    print(f'  {f}')
print('=' * 70)
