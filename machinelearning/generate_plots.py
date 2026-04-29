#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generate radar plots for training and test results (Nature style)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os

# =========================
# 0. Nature Style Settings
# =========================
config = {
    "font.family": 'serif',
    "font.serif": ['Times New Roman'], 
    "mathtext.fontset": 'stix',
    "font.size": 12,
    "axes.linewidth": 1.0,
    "figure.dpi": 300,
}
rcParams.update(config)

PLOTS_DIR = './plots'
os.makedirs(PLOTS_DIR, exist_ok=True)

# =========================
# 1. Helper Functions
# =========================

def _draw_ring_label(ax, text, angle_rad, radius, fontsize=10):
    """Draw labels on the outer ring of radar chart"""
    display_theta = ax.get_theta_direction() * angle_rad + ax.get_theta_offset()
    angle_deg = (np.degrees(display_theta) + 360) % 360
    rotation = angle_deg - 90
    rotation = (rotation + 180) % 360 - 180
    if rotation < -90 or rotation > 90:
        rotation += 180
        rotation = (rotation + 180) % 360 - 180
    ax.text(
        angle_rad, radius, text,
        ha="center", va="center",
        rotation=rotation, rotation_mode="anchor",
        fontsize=fontsize, fontweight="bold",
        color="black", zorder=5,
    )

def plot_single_radar(ax, model_names, metric_names, values, line_colors, ring_colors,
                      title="", legend=True, legend_fontsize=9, label_fontsize=9):
    """
    Plot a single radar chart
    
    Parameters:
    -----------
    ax : matplotlib axis with polar projection
    model_names : list of str
        Names of models to plot
    metric_names : list of str
        Names of metrics (will be displayed on outer ring)
    values : numpy array
        Shape (n_models, n_metrics), values should be in range [0, 1]
    line_colors : list of str
        Colors for each model line
    ring_colors : list of str
        Colors for outer ring segments
    title : str
        Title for this radar chart
    legend : bool
        Whether to show legend
    legend_fontsize : float
        Font size for legend
    label_fontsize : float
        Font size for metric labels on outer ring
    """
    from matplotlib.lines import Line2D
    
    n_models = len(model_names)
    n_metrics = len(metric_names)
    angles = np.linspace(0, 2 * np.pi, n_metrics, endpoint=False)
    angles_closed = np.concatenate([angles, [angles[0]]])
    
    # Fixed range: 0 to 1
    r_grid_max = 1.0
    r_ring_inner = r_grid_max
    r_ring_width = r_grid_max * 0.15
    r_outer = r_ring_inner + r_ring_width
    r_ring_mid = r_ring_inner + r_ring_width / 2
    
    ax.set_facecolor("white")
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_ylim(0, r_outer + 0.05)
    ax.set_axisbelow(True)
    
    # Draw outer ring with metric labels
    seg_width = 2 * np.pi / n_metrics
    for i, (ang, label) in enumerate(zip(angles, metric_names)):
        ax.bar(
            ang - seg_width / 2, r_ring_width,
            width=seg_width, bottom=r_ring_inner,
            align="edge",
            color=ring_colors[i % len(ring_colors)],
            edgecolor="white", linewidth=1.5, zorder=1,
        )
        _draw_ring_label(ax, label, ang, radius=r_ring_mid, fontsize=label_fontsize)
    
    # Grid settings - Fixed ticks at 0, 0.25, 0.5, 0.75, 1.0
    ax.set_xticks(angles)
    ax.set_xticklabels([""] * n_metrics)
    ax.set_rlabel_position(0)
    
    yticks = [0.25, 0.5, 0.75, 1.0]
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{t:.2f}" for t in yticks], fontsize=8, fontweight="bold", color="black")
    
    grid_color = "#7F7F7F"
    ax.yaxis.grid(True, linestyle="--", linewidth=0.8, color=grid_color, alpha=0.6)
    ax.xaxis.grid(True, linestyle="--", linewidth=0.7, color=grid_color, alpha=0.5)
    
    ax.spines["polar"].set_color("#333333")
    ax.spines["polar"].set_linewidth(1.5)
    theta = np.linspace(0, 2 * np.pi, 720)
    ax.plot(theta, np.full_like(theta, r_outer), color="#333333", linewidth=1.0, alpha=0.9, zorder=2)
    
    # Plot data for each model
    for i in range(n_models):
        row = values[i, :]
        row_closed = np.concatenate([row, [row[0]]])
        c = line_colors[i % len(line_colors)]
        ax.fill(angles_closed, row_closed, color=c, alpha=0.05, zorder=3)
        ax.plot(
            angles_closed, row_closed,
            color=c, linewidth=1.4,
            marker="o", markersize=4,
            markerfacecolor=c,
            markeredgecolor="white",
            markeredgewidth=0.8,
            zorder=4,
        )
    
    # Title
    if title:
        ax.set_title(title, fontsize=11, fontweight='bold', pad=10)
    
    # Legend
    if legend:
        legend_elements = [
            Line2D([0], [0], color=line_colors[i % len(line_colors)],
                   marker="o", markersize=4,
                   markerfacecolor=line_colors[i % len(line_colors)],
                   markeredgecolor="white", markeredgewidth=0.8,
                   linewidth=1.4, label=model_names[i])
            for i in range(n_models)
        ]
        ax.legend(
            handles=legend_elements,
            loc="upper left",
            bbox_to_anchor=(1.05, 1.0),
            fontsize=legend_fontsize,
            frameon=True,
            fancybox=False,
            shadow=False,
            edgecolor="black",
            framealpha=1.0,
        )

# =========================
# 2. Color Schemes
# =========================

line_colors = ["#C00000", "#2E75B6", "#E69F00", "#2CA02C", "#9467BD", "#8C564B", "#D55E00"]
ring_colors = [
    "#B7DEE8", "#C5D9F1", "#F4CCCC", "#F9CB9C",
    "#FFE599", "#D9EAD3", "#D0E0E3", "#D9D2E9",
]

# =========================
# 3. Read Data
# =========================

df_lodo = pd.read_csv('./results/lodo_cv_results.csv')
df_test = pd.read_csv('./results/external_test_results.csv')

# Rename ROC_AUC to AUROC for display
df_lodo = df_lodo.rename(columns={'ROC_AUC': 'AUROC'})
df_test = df_test.rename(columns={'ROC_AUC': 'AUROC'})

# Handle NaN values in AUROC - replace with 0.5 (random classifier baseline)
df_lodo['AUROC'] = df_lodo['AUROC'].fillna(0.5)
df_test['AUROC'] = df_test['AUROC'].fillna(0.5)
df_lodo['AUPRC'] = df_lodo['AUPRC'].fillna(0.5)
df_test['AUPRC'] = df_test['AUPRC'].fillna(0.5)

train_studies = ['study_1', 'study_4', 'study_6', 'study_7', 'study_9', 'study_11']
test_studies = df_test['Test_Study'].unique().tolist()
models = df_lodo['Model'].unique().tolist()
metrics = ['AUROC', 'AUPRC', 'Accuracy', 'Precision', 'Recall', 'F1']

# =========================
# 4. Plot 1: Training + Test Summary (2 radars)
# =========================

print("Generating Plot 1: Training + Test Summary...")

fig = plt.figure(figsize=(16, 7))
fig.patch.set_facecolor('white')

# Left: Training Summary (average across 6 studies)
ax1 = fig.add_subplot(121, projection='polar')
train_summary = []
for model in models:
    model_data = df_lodo[df_lodo['Model'] == model]
    row = [model_data[metric].mean() for metric in metrics]
    train_summary.append(row)
train_summary = np.array(train_summary)

plot_single_radar(
    ax1, models, metrics, train_summary,
    line_colors, ring_colors,
    title="Training Set (LODO CV Average)",
    legend=True, legend_fontsize=9, label_fontsize=9
)

# Right: Test Summary (average across all test studies)
ax2 = fig.add_subplot(122, projection='polar')
test_summary = []
for model in models:
    model_data = df_test[df_test['Model'] == model]
    row = [model_data[metric].mean() for metric in metrics]
    test_summary.append(row)
test_summary = np.array(test_summary)

plot_single_radar(
    ax2, models, metrics, test_summary,
    line_colors, ring_colors,
    title="External Test Set (Average)",
    legend=True, legend_fontsize=9, label_fontsize=9
)

plt.suptitle("Training vs Test Performance Summary", fontsize=16, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(os.path.join(PLOTS_DIR, 'summary_train_test.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(PLOTS_DIR, 'summary_train_test.pdf'), bbox_inches='tight')
plt.close()

print(f"  Saved: summary_train_test.png/pdf")

# =========================
# 5. Plot 2: Training Studies (6 radars)
# =========================

print("Generating Plot 2: Training Studies (6 radars)...")

fig = plt.figure(figsize=(18, 12))
fig.patch.set_facecolor('white')

for idx, study in enumerate(train_studies):
    ax = fig.add_subplot(2, 3, idx+1, projection='polar')
    
    study_data = []
    for model in models:
        model_row = df_lodo[(df_lodo['Test_Study'] == study) & (df_lodo['Model'] == model)]
        row = [model_row[metric].values[0] if len(model_row) > 0 else 0 for metric in metrics]
        study_data.append(row)
    study_data = np.array(study_data)
    
    plot_single_radar(
        ax, models, metrics, study_data,
        line_colors, ring_colors,
        title=f"{study}",
        legend=(idx == 0),  # Only show legend on first subplot
        legend_fontsize=8, label_fontsize=8
    )

plt.suptitle("Training Studies Performance (LODO CV)", fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.99])
plt.savefig(os.path.join(PLOTS_DIR, 'training_studies_6radars.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(PLOTS_DIR, 'training_studies_6radars.pdf'), bbox_inches='tight')
plt.close()

print(f"  Saved: training_studies_6radars.png/pdf")

# =========================
# 6. Plot 3: Test Studies (n radars)
# =========================

print(f"Generating Plot 3: Test Studies ({len(test_studies)} radars)...")

# Calculate grid layout
n_test = len(test_studies)
ncols = 3
nrows = (n_test + ncols - 1) // ncols

fig = plt.figure(figsize=(18, 6*nrows))
fig.patch.set_facecolor('white')

for idx, study in enumerate(test_studies):
    ax = fig.add_subplot(nrows, ncols, idx+1, projection='polar')
    
    study_data = []
    for model in models:
        model_row = df_test[(df_test['Test_Study'] == study) & (df_test['Model'] == model)]
        row = [model_row[metric].values[0] if len(model_row) > 0 else 0 for metric in metrics]
        study_data.append(row)
    study_data = np.array(study_data)
    
    plot_single_radar(
        ax, models, metrics, study_data,
        line_colors, ring_colors,
        title=f"{study}",
        legend=(idx == 0),  # Only show legend on first subplot
        legend_fontsize=8, label_fontsize=8
    )

plt.suptitle("External Test Studies Performance", fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.99])
plt.savefig(os.path.join(PLOTS_DIR, 'test_studies_nradars.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(PLOTS_DIR, 'test_studies_nradars.pdf'), bbox_inches='tight')
plt.close()

print(f"  Saved: test_studies_nradars.png/pdf")

print("\n" + "="*80)
print("All plots generated successfully!")
print("="*80)
print(f"Plots saved to: {PLOTS_DIR}")
print("\nGenerated 3 plots:")
print("  1. summary_train_test: Training + Test summary (2 radars)")
print("  2. training_studies_6radars: 6 training studies (6 radars)")
print(f"  3. test_studies_nradars: {len(test_studies)} test studies ({len(test_studies)} radars)")
