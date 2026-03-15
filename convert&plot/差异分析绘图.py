import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# ─── 路径配置 ───
meta_path = 'F:/1a-EOD-CSF-protein/meta/'
out_dir   = 'F:/1a-EOD-CSF-protein/1a-figure/差异分析+富集/'

import os
os.makedirs(out_dir, exist_ok=True)

# ─── 两张图的配置 ───
# 每组：(group1_name, group2_name, file_g1_vs_CN, file_g2_vs_CN, file_inter, output_stem)
COMPARISONS = [
    (
        'CN', 'EOD', 'LOD',
        'EOD_vs_CN.csv',
        'LOD_vs_CN.csv',
        'EOD_vs_LOD.csv',
        'CN_EOD_LOD'
    ),
    (
        'CN', 'EOAD', 'LOAD',
        'EOAD_vs_CN.csv',
        'LOAD_vs_CN.csv',
        'EOAD_vs_LOAD.csv',
        'CN_EOAD_LOAD'
    ),
]

N_STUDIES_MIN = 2   # N_Studies > 2
THRESH_P     = 0.05

# ─── 颜色 ───
COLORS = {
    'Both':    '#5da5da',
    'G1_only': '#faa43a',
    'G2_only': '#60bd68',
    'NS':      '#d3d3d3',
}


def load_meta(fname, effect_col='Weighted_Effect', p_col='FDR_BH_Stratified'):
    df = pd.read_csv(meta_path + fname)
    df = df[df['N_Studies'] > N_STUDIES_MIN][['Protein', effect_col, p_col]].copy()
    df.columns = ['Protein', 'Effect', 'P']
    return df


def build_df(g1_name, g2_name, f_g1, f_g2, f_inter):
    d1 = load_meta(f_g1)
    d2 = load_meta(f_g2)
    di = load_meta(f_inter)  # inter: g1 vs g2 直接比较

    d1.columns = ['Protein', f'{g1_name}_Effect', f'{g1_name}_P']
    d2.columns = ['Protein', f'{g2_name}_Effect', f'{g2_name}_P']
    di.columns = ['Protein', 'Inter_Effect', 'Inter_P']

    df = d1.merge(d2, on='Protein', how='inner')
    df = df.merge(di[['Protein', 'Inter_P']], on='Protein', how='left')

    def sig_group(row):
        s1 = row[f'{g1_name}_P'] < THRESH_P
        s2 = row[f'{g2_name}_P'] < THRESH_P
        if s1 and s2:   return 'Both'
        elif s1:        return 'G1_only'
        elif s2:        return 'G2_only'
        else:           return 'NS'

    df['Sig']       = df.apply(sig_group, axis=1)
    df['Inter_sig'] = df['Inter_P'].notna() & (df['Inter_P'] < THRESH_P)
    return df


def plot_one(df, cn_name, g1_name, g2_name, out_stem):
    """
    横轴：g1 vs CN 的效应值
    纵轴：g2 vs CN 的效应值
    颜色：Both / G1_only / G2_only / NS
    形状：圆圈=互作不显著，三角=互作显著
    """
    fig, ax = plt.subplots(figsize=(8, 7))
    
    color_map = {
        'Both':   COLORS['Both'],
        'G1_only': COLORS['G1_only'],
        'G2_only': COLORS['G2_only'],
        'NS':     COLORS['NS'],
    }
    label_map = {
        'Both':   'Both',
        'G1_only': f'{g1_name} only',
        'G2_only': f'{g2_name} only',
        'NS':     'NS',
    }

    draw_order = ['NS', 'Both', 'G1_only', 'G2_only']

    for grp in draw_order:
        sub = df[df['Sig'] == grp]
        # 圆圈（互作不显著）
        circ = sub[~sub['Inter_sig']]
        if len(circ):
            ax.scatter(circ[f'{g1_name}_Effect'], circ[f'{g2_name}_Effect'],
                       c=color_map[grp], marker='o', alpha=0.8, s=30,
                       edgecolors='none',
                       zorder=2 if grp == 'NS' else 3)
        # 三角（互作显著）
        tri = sub[sub['Inter_sig']]
        if len(tri):
            ax.scatter(tri[f'{g1_name}_Effect'], tri[f'{g2_name}_Effect'],
                       c=color_map[grp], marker='^', alpha=0.9, s=55,
                       edgecolors='black', linewidths=0.5,
                       zorder=4 if grp == 'NS' else 5)

    # 辅助线
    lim = 0.55
    ax.axhline(0, color='black', linewidth=0.8, zorder=1)
    ax.axvline(0, color='black', linewidth=0.8, zorder=1)
    ax.plot([-lim, lim], [-lim, lim], 'k--', alpha=0.3, zorder=1)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    # 图例 — 颜色
    color_handles = [
        Patch(facecolor=color_map[g], label=label_map[g])
        for g in ['Both', 'G1_only', 'G2_only', 'NS']
    ]
    # 图例 — 形状（互作）
    shape_handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',
               markeredgecolor='black', markersize=8, label='Interaction NS'),
        Line2D([0], [0], marker='^', color='w', markerfacecolor='grey',
               markeredgecolor='black', markersize=8, label='Interaction sig.'),
    ]

    leg1 = ax.legend(handles=color_handles, title='Significance vs. CN',
                     loc='upper left', bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=9)
    ax.add_artist(leg1)
    ax.legend(handles=shape_handles, title=f'{g1_name} vs. {g2_name}',
              loc='upper left', bbox_to_anchor=(1.02, 0.55), frameon=False, fontsize=9)

    ax.set_title(f"{cn_name} / {g1_name} / {g2_name}", fontsize=13, fontweight='bold', pad=12)
    ax.set_xlabel(f"Effect size: {g1_name} vs. {cn_name}", fontsize=11)
    ax.set_ylabel(f"Effect size: {g2_name} vs. {cn_name}", fontsize=11)
    
    plt.tight_layout()

    png_path = out_dir + out_stem + '.png'
    pdf_path = out_dir + out_stem + '.pdf'
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close()
    print(f"已保存: {png_path}")
    print(f"已保存: {pdf_path}")


# ─── 主流程 ───
for cn_name, g1_name, g2_name, f_g1, f_g2, f_inter, stem in COMPARISONS:
    df = build_df(g1_name, g2_name, f_g1, f_g2, f_inter)
    print(f"\n{g1_name} vs {g2_name}: {len(df)} 个蛋白")
    print(df['Sig'].value_counts().to_string())
    plot_one(df, cn_name, g1_name, g2_name, stem)

print("\n全部完成！")
