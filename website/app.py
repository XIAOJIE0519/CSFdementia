from flask import Flask, render_template, jsonify, request, send_from_directory
import pandas as pd
import numpy as np
import os, io, warnings, logging, joblib
warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__)
app.config['JSON_AS_ASCII'] = False
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
data_cache = {}

CORE_FEATURES = ['MIF','DDAH1','ENO2','PEBP1','PAM','SPON1','SOD1','RTN4R','VASN','SOD2','GFRA2','CA4','CANT1','GLRX','GAS6']
ML_THRESHOLD = 0.55
CELL_TYPE_COLS = ['neurons','astrocytes','oligodendrocytes','microglia','endothelial','OPCs']


def _safe(val):
    if val is None: return None
    try:
        if np.isnan(val) or np.isinf(val): return None
    except Exception: pass
    return val


def _read_xlsx_safe(path):
    try:
        return pd.read_excel(path, engine='openpyxl')
    except Exception:
        try:
            import openpyxl
            wb = openpyxl.load_workbook(path, data_only=True, keep_links=False)
            ws = wb.active
            rows = list(ws.iter_rows(values_only=True))
            if not rows: return pd.DataFrame()
            return pd.DataFrame(rows[1:], columns=rows[0])
        except Exception as e:
            logger.error(f"xlsx failed {path}: {e}")
            return pd.DataFrame()


def _read_xlsx_safe_sheet(path, sheet):
    try:
        return pd.read_excel(path, sheet_name=sheet, engine='openpyxl')
    except Exception as e:
        logger.error(f"xlsx sheet {sheet} failed {path}: {e}")
        return pd.DataFrame()


def load_data():
    global data_cache
    try:
        logger.info("Loading data...")
        meta_comps = ['EOAD_vs_CN','LOAD_vs_CN','EOD_vs_CN','LOD_vs_CN','EOAD_vs_LOAD','EOD_vs_LOD']
        data_cache['meta'] = {}
        data_cache['meta_significant'] = {}
        for comp in meta_comps:
            fp = os.path.join(DATA_DIR, 'DE', f'{comp}.csv')
            if os.path.exists(fp):
                df = pd.read_csv(fp)
                data_cache['meta'][comp] = df
                sig = df[df['FDR_BH_Stratified'] < 0.05].sort_values('FDR_BH_Stratified')
                data_cache['meta_significant'][comp] = sig
                logger.info(f"  {comp}: {len(df)} total, {len(sig)} sig")

        study_dfs = []
        de_dir = os.path.join(DATA_DIR, 'DE')
        for fname in os.listdir(de_dir):
            if fname.startswith('study_') and fname.endswith('.csv'):
                try: study_dfs.append(pd.read_csv(os.path.join(de_dir, fname), low_memory=False))
                except Exception as e: logger.warning(f"  skip {fname}: {e}")
        data_cache['study_de'] = pd.concat(study_dfs, ignore_index=True) if study_dfs else pd.DataFrame()
        logger.info(f"  study DE rows: {len(data_cache['study_de'])}")

        data_cache['enrichment_gokegg'] = {}
        data_cache['enrichment_gsea'] = {}
        gk_path = os.path.join(DATA_DIR, 'enrichment', 'GOKEGG_all_results.csv')
        gs_path = os.path.join(DATA_DIR, 'enrichment', 'GSEA_all_results.csv')
        if os.path.exists(gk_path):
            gk = pd.read_csv(gk_path, low_memory=False)
            for c in gk['Comparison'].unique():
                data_cache['enrichment_gokegg'][c] = gk[gk['Comparison']==c].sort_values('FDR').head(50)
        if os.path.exists(gs_path):
            gs = pd.read_csv(gs_path, low_memory=False)
            for c in gs['Comparison'].unique():
                data_cache['enrichment_gsea'][c] = gs[gs['Comparison']==c].sort_values('p.adjust').head(10)
        logger.info(f"  GOKEGG comps: {list(data_cache['enrichment_gokegg'].keys())}")

        data_cache['wgcna'] = {
            'module_assignments': pd.read_csv(os.path.join(DATA_DIR,'wgcna','consensus_module_assignments.csv')),
            'main_corr':   _read_xlsx_safe_sheet(os.path.join(DATA_DIR,'wgcna','consensus_main_heatmap_data.xlsx'), 0),
            'main_pval':   _read_xlsx_safe_sheet(os.path.join(DATA_DIR,'wgcna','consensus_main_heatmap_data.xlsx'), 1),
            'test_corr':   _read_xlsx_safe_sheet(os.path.join(DATA_DIR,'wgcna','consensus_test_heatmap_data.xlsx'), 0),
            'test_pval':   _read_xlsx_safe_sheet(os.path.join(DATA_DIR,'wgcna','consensus_test_heatmap_data.xlsx'), 1),
        }
        logger.info(f"  WGCNA main: {data_cache['wgcna']['main_corr'].shape}")

        data_cache['corr_index'] = {}
        for grp in ['CN', 'EOD', 'LOD']:
            fpath = os.path.join(DATA_DIR, 'correlation', f'{grp}.csv')
            if not os.path.exists(fpath):
                data_cache['corr_index'][grp] = {}
                continue
            logger.info(f'  Loading correlation {grp}...')
            df = pd.read_csv(fpath, low_memory=False,
                             usecols=['Protein1','Protein2','Pearson_meta','P_value','FDR_BH_Stratified','N_Studies'])
            df['_abs'] = df['Pearson_meta'].abs()
            df.sort_values('_abs', ascending=False, inplace=True)
            # Build compact tuple index: {protein: [(partner, corr, pval, fdr, n_studies), ...]}
            idx = {}
            for row in df.itertuples(index=False):
                p1, p2 = row.Protein1, row.Protein2
                corr = float(row.Pearson_meta) if row.Pearson_meta==row.Pearson_meta else None
                pval = float(row.P_value)      if row.P_value==row.P_value           else None
                fdr  = float(row.FDR_BH_Stratified) if row.FDR_BH_Stratified==row.FDR_BH_Stratified else None
                ns   = int(row.N_Studies)      if row.N_Studies==row.N_Studies       else None
                if p1 not in idx: idx[p1] = []
                if p2 not in idx: idx[p2] = []
                if len(idx[p1]) < 10: idx[p1].append((p2, corr, pval, fdr, ns))
                if len(idx[p2]) < 10: idx[p2].append((p1, corr, pval, fdr, ns))
            data_cache['corr_index'][grp] = idx
            del df
            logger.info(f'  Correlation {grp}: {len(idx)} proteins indexed')

        try:
            mp = os.path.join(DATA_DIR,'machine','unified_LogisticRegression.pkl')
            data_cache['ml_model'] = joblib.load(mp) if os.path.exists(mp) else None
            logger.info(f"  ML model: {'loaded' if data_cache['ml_model'] else 'NOT FOUND'}")
        except Exception as e:
            logger.warning(f"  ML model load failed (non-fatal): {e}")
            data_cache['ml_model'] = None
        logger.info("All data loaded.")
        return True
    except Exception as e:
        logger.error(f"load_data failed: {e}")
        import traceback; traceback.print_exc()
        return False


def _query_correlation(protein_id, group, top_n=10):
    idx = data_cache.get('corr_index', {}).get(group, {})
    tuples = idx.get(protein_id, [])
    return [{
        'partner':     t[0],
        'correlation': _safe(t[1]),
        'p_value':     _safe(t[2]),
        'fdr':         _safe(t[3]),
        'n_studies':   t[4],
    } for t in tuples[:top_n]]


def _parse_heatmap_row(corr_row, pval_row, hm_cols):
    """Parse heatmap row using separate correlation and p-value rows.
    - hm_cols: columns from the correlation sheet (excluding 'Module')
    - pval_row: corresponding row from p-value sheet (same Module, different col names)
    Cell-type cols have _pval suffix in pval sheet; study_ cols have no p-value.
    Correlation sheet may have duplicate-renamed cols (EOD.1, EOAD.1) which are ignored.
    """
    traits = []
    # Build p-value lookup from pval_row by position (same order as corr sheet)
    # pval sheet cols (excl Module) align positionally with corr sheet cols (excl Module)
    pval_cols = [c for c in (pval_row.index.tolist() if pval_row is not None else []) if c != 'Module']
    corr_cols_clean = [c for c in hm_cols if c != 'Module']
    # Build positional map: corr_col -> pval value
    pval_by_pos = {}
    for i, cc in enumerate(corr_cols_clean):
        if i < len(pval_cols) and pval_row is not None:
            pval_by_pos[cc] = pval_row.get(pval_cols[i])
        else:
            pval_by_pos[cc] = None
    # Skip duplicate-renamed cols (pandas adds .1, .2 suffix for duplicate names)
    skip = {c for c in corr_cols_clean if '.' in c and c.split('.')[-1].isdigit()}
    for col in corr_cols_clean:
        if col in skip:
            continue
        corr_val = corr_row.get(col)
        if col in CELL_TYPE_COLS:
            # Cell-type ORA: corr sheet stores the p-value directly
            pval = corr_row.get(col)
            traits.append({'trait': col, 'type': 'cell_type_ora',
                           'p_value': _safe(float(pval)) if pval is not None and pd.notna(pval) else None,
                           'correlation': None, 'n_studies': None})
        elif col.startswith('study_'):
            # Conservancy: corr sheet stores the value directly
            val = corr_row.get(col)
            traits.append({'trait': col, 'type': 'conservancy',
                           'n_studies': _safe(float(val)) if val is not None and pd.notna(val) else None,
                           'correlation': None, 'p_value': None})
        else:
            pval = pval_by_pos.get(col)
            traits.append({'trait': col, 'type': 'trait_corr',
                           'correlation': _safe(float(corr_val)) if corr_val is not None and pd.notna(corr_val) else None,
                           'p_value': _safe(float(pval)) if pval is not None and pd.notna(pval) else None,
                           'n_studies': None})
    return traits


def _build_shap_force_png(protein_values):
    import shap, matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from sklearn.linear_model import LogisticRegression as LR
    model = data_cache['ml_model']
    if model is None:
        raise ValueError("ML model not loaded")
    x = np.array([[protein_values.get(f, 0.0) for f in CORE_FEATURES]], dtype=float)
    coef_mean      = np.mean([p.named_steps['classifier'].coef_[0]      for p in model.estimators_], axis=0)
    intercept_mean = np.mean([p.named_steps['classifier'].intercept_[0] for p in model.estimators_])
    surrogate = LR(C=0.1, max_iter=2000, random_state=42)
    surrogate.fit(np.zeros((2, len(CORE_FEATURES))), [0, 1])
    surrogate.coef_[0]      = coef_mean
    surrogate.intercept_[0] = intercept_mean
    explainer = shap.LinearExplainer(surrogate, np.zeros((1, len(CORE_FEATURES))),
                                     feature_perturbation='interventional')
    sv = explainer.shap_values(x)
    if isinstance(sv, list):
        sv = sv[1][0] if len(sv) > 1 else sv[0][0]
    else:
        sv = sv[0]
    base_val = float(np.atleast_1d(explainer.expected_value)[0])
    prob = float(surrogate.predict_proba(x)[0, 1])
    pos_f = sorted([(CORE_FEATURES[i], sv[i], x[0,i]) for i in range(15) if sv[i] >= 0], key=lambda t: t[1], reverse=True)
    neg_f = sorted([(CORE_FEATURES[i], sv[i], x[0,i]) for i in range(15) if sv[i] <  0], key=lambda t: t[1])
    fig, ax = plt.subplots(figsize=(13, 3.0))
    fig.patch.set_facecolor('#f8f9fa'); ax.set_facecolor('#f8f9fa')
    y0, bar_h = 0.5, 0.55
    cursor = base_val
    for fname, val, fval in pos_f:
        ax.barh(y0, val, left=cursor, height=bar_h, color='#d73027', alpha=0.88, edgecolor='white', linewidth=0.5)
        if abs(val) > 0.012:
            ax.text(cursor+val/2, y0, f'{fname}\n={fval:.2f}', ha='center', va='center',
                    fontsize=6.5, color='white', fontweight='bold', fontfamily='monospace')
        cursor += val
    cursor = base_val
    for fname, val, fval in neg_f:
        ax.barh(y0, val, left=cursor, height=bar_h, color='#4575b4', alpha=0.88, edgecolor='white', linewidth=0.5)
        if abs(val) > 0.012:
            ax.text(cursor+val/2, y0, f'{fname}\n={fval:.2f}', ha='center', va='center',
                    fontsize=6.5, color='white', fontweight='bold', fontfamily='monospace')
        cursor += val
    ax.axvline(base_val, color='#555', lw=1.2, ls='--', alpha=0.7)
    ax.axvline(prob, color='#222', lw=2.0, alpha=0.9)
    ax.set_yticks([]); ax.set_ylim(0, 1)
    for sp in ['top','left','right']: ax.spines[sp].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.8)
    ax.set_xlabel('Model output value', fontsize=9)
    ax.set_title(f'SHAP Force Plot  |  EOD Risk Score: {prob:.3f}  (threshold={ML_THRESHOLD})',
                 fontsize=10, fontweight='bold', pad=10)
    ax.legend(handles=[
        Patch(facecolor='#d73027', alpha=0.88, label='Increases EOD risk'),
        Patch(facecolor='#4575b4', alpha=0.88, label='Decreases EOD risk'),
    ], loc='upper right', fontsize=8, frameon=True, framealpha=0.8)
    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=180, bbox_inches='tight', facecolor='#f8f9fa')
    plt.close(fig); buf.seek(0)
    return buf.getvalue()


# ───────────────────────── routes ─────────────────────────

@app.route('/data/figure/<path:filename>')
def serve_figure(filename):
    return send_from_directory(os.path.join(DATA_DIR, 'figure'), filename)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/protein/<protein_id>')
def protein_detail_page(protein_id):
    return render_template('protein_detail.html', protein_id=protein_id)

@app.route('/api/search_protein')
def search_protein():
    query = request.args.get('q', '').strip().upper()
    if not query: return jsonify({'proteins': []})
    results = set()
    for df in data_cache['meta'].values():
        hits = df[df['Protein'].str.upper().str.contains(query, na=False)]['Protein']
        results.update(hits.tolist())
    return jsonify({'proteins': sorted(results)[:10]})

@app.route('/api/comparisons')
def get_comparisons():
    return jsonify({'comparisons': [
        {'id': 'EOAD_vs_CN',   'name': 'EOAD vs CN'},
        {'id': 'LOAD_vs_CN',   'name': 'LOAD vs CN'},
        {'id': 'EOD_vs_CN',    'name': 'EOD vs CN'},
        {'id': 'LOD_vs_CN',    'name': 'LOD vs CN'},
        {'id': 'EOAD_vs_LOAD', 'name': 'EOAD vs LOAD'},
        {'id': 'EOD_vs_LOD',   'name': 'EOD vs LOD'},
    ]})

@app.route('/api/differential_expression/<comparison>')
def get_differential_expression(comparison):
    try:
        if comparison not in data_cache['meta_significant']:
            return jsonify({'error': 'Comparison not found'}), 404
        sig_df  = data_cache['meta_significant'][comparison]
        meta_df = data_cache['meta'][comparison]
        volcano_data = []
        for _, row in sig_df.iterrows():
            fdr = row['FDR_BH_Stratified']
            if pd.notna(row.get('Weighted_Effect')) and pd.notna(fdr) and fdr > 0:
                volcano_data.append({
                    'protein': row['Protein'],
                    'log2fc':  _safe(float(row['Weighted_Effect'])),
                    'neg_log10_fdr': _safe(float(-np.log10(fdr))),
                    'fdr':     _safe(float(fdr)),
                    'n_studies': int(row['N_Studies']) if pd.notna(row.get('N_Studies')) else None,
                })
        sig_proteins = []
        for _, row in sig_df.head(50).iterrows():
            sig_proteins.append({
                'protein':   row['Protein'],
                'log2fc':    _safe(float(row['Weighted_Effect']))   if pd.notna(row.get('Weighted_Effect'))   else None,
                'fdr':       _safe(float(row['FDR_BH_Stratified'])) if pd.notna(row.get('FDR_BH_Stratified')) else None,
                'n_studies': int(row['N_Studies'])                  if pd.notna(row.get('N_Studies'))         else None,
            })
        enrichment_gokegg = []
        enrichment_gsea = []
        if comparison in data_cache['enrichment_gokegg']:
            for _, row in data_cache['enrichment_gokegg'][comparison].iterrows():
                enrichment_gokegg.append({
                    'source': row.get('Source','ORA'), 'term_id': row.get('Term_ID',''),
                    'term_name': row.get('Term_Name',''), 'gene_ratio': row.get('GeneRatio',''),
                    'p_value': _safe(float(row['P_value'])) if pd.notna(row.get('P_value')) else None,
                    'fdr':     _safe(float(row['FDR']))     if pd.notna(row.get('FDR'))     else None,
                    'count':   int(row['Count'])            if pd.notna(row.get('Count'))   else 0,
                    'genes': row.get('Genes',''), 'type': 'ORA',
                })
        if comparison in data_cache['enrichment_gsea']:
            for _, row in data_cache['enrichment_gsea'][comparison].iterrows():
                enrichment_gsea.append({
                    'source': 'GSEA', 'term_id': row.get('ID',''),
                    'term_name': row.get('Description',''), 'gene_ratio': str(row.get('setSize','')),
                    'p_value': _safe(float(row['pvalue']))   if pd.notna(row.get('pvalue'))   else None,
                    'fdr':     _safe(float(row['p.adjust'])) if pd.notna(row.get('p.adjust')) else None,
                    'count':   int(row['setSize'])           if pd.notna(row.get('setSize'))  else 0,
                    'genes': row.get('core_enrichment',''), 'type': 'GSEA',
                    'nes': _safe(float(row['NES'])) if pd.notna(row.get('NES')) else None,
                })
        return jsonify({
            'comparison': comparison, 'volcano_data': volcano_data,
            'significant_proteins': sig_proteins, 'total_proteins': len(meta_df),
            'significant_count': len(sig_df),
            'enrichment_gokegg': enrichment_gokegg,
            'enrichment_gsea': enrichment_gsea,
            'enrichment': enrichment_gokegg + enrichment_gsea,
        })
    except Exception as e:
        logger.error(f"DE error: {e}"); import traceback; traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@app.route('/api/protein/<protein_id>')
def get_protein_detail(protein_id):
    try:
        result = {'protein_id': protein_id, 'differential_expression': {}, 'meta_analysis': {}, 'external_links': {}}
        gene_name  = protein_id.split('|')[0]
        uniprot_id = protein_id.split('|')[1] if '|' in protein_id else protein_id
        result['external_links'] = {
            'uniprot': f"https://www.uniprot.org/uniprotkb/{uniprot_id}",
            'genecard': f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene_name}",
        }
        study_df = data_cache.get('study_de', pd.DataFrame())
        if not study_df.empty:
            pdata = study_df[study_df['Protein'] == protein_id]
            for _, row in pdata.iterrows():
                comp = row.get('Comparison',''); study = row.get('Study','')
                if comp not in result['differential_expression']:
                    result['differential_expression'][comp] = []
                result['differential_expression'][comp].append({
                    'study':   study,
                    'log2fc':  _safe(float(row['Log2FC']))  if pd.notna(row.get('Log2FC'))  else None,
                    'p_value': _safe(float(row['P_value'])) if pd.notna(row.get('P_value')) else None,
                    'fdr':     _safe(float(row['FDR']))     if pd.notna(row.get('FDR'))     else None,
                })
        if not result['differential_expression']:
            return jsonify({'error': 'Protein not found'}), 404
        for comp, mdf in data_cache['meta'].items():
            pmeta = mdf[mdf['Protein'] == protein_id]
            if len(pmeta) > 0:
                row = pmeta.iloc[0]
                result['meta_analysis'][comp] = {
                    'weighted_effect': _safe(float(row['Weighted_Effect']))    if pd.notna(row.get('Weighted_Effect'))    else None,
                    'ci_lower':        _safe(float(row['CI_Lower']))           if pd.notna(row.get('CI_Lower'))           else None,
                    'ci_upper':        _safe(float(row['CI_Upper']))           if pd.notna(row.get('CI_Upper'))           else None,
                    'p_value':         _safe(float(row['P_value']))            if pd.notna(row.get('P_value'))            else None,
                    'fdr':             _safe(float(row['FDR_BH_Stratified'])) if pd.notna(row.get('FDR_BH_Stratified')) else None,
                    'n_studies':       int(row['N_Studies'])                   if pd.notna(row.get('N_Studies'))         else None,
                    'studies':         row['Studies']                          if pd.notna(row.get('Studies'))           else None,
                }
        mod_df = data_cache['wgcna']['module_assignments']
        pmod   = mod_df[mod_df['Protein'] == protein_id]
        if len(pmod) > 0:
            result['module'] = {'module': pmod.iloc[0]['Module'], 'module_name': pmod.iloc[0]['Module_Name']}
        corr_by_group = {}
        for grp in ['CN','EOD','LOD']:
            corr_by_group[grp] = _query_correlation(protein_id, grp, top_n=10)
        result['correlations_by_group'] = corr_by_group
        result['correlations'] = corr_by_group.get('EOD', [])
        return jsonify(result)
    except Exception as e:
        logger.error(f"Protein detail error: {e}"); import traceback; traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@app.route('/api/modules')
def get_modules():
    try:
        main_corr = data_cache['wgcna']['main_corr']
        mod_asgn  = data_cache['wgcna']['module_assignments']
        prefix_to_asgn = {}
        for mn in mod_asgn['Module_Name'].unique():
            prefix = mn.split('_')[0]
            prefix_to_asgn[prefix] = mn
        modules = []
        for _, row in main_corr.iterrows():
            mname  = row['Module']
            prefix = mname.split('_')[0]
            asgn_name = prefix_to_asgn.get(prefix)
            color_src  = asgn_name.split('_')[1] if asgn_name and '_' in asgn_name else mname.split('_')[1] if '_' in mname else mname
            cnt = len(mod_asgn[mod_asgn['Module_Name'] == asgn_name]) if asgn_name else 0
            modules.append({'module_name': mname, 'color': color_src,
                            'protein_count': cnt, 'asgn_name': asgn_name})
        return jsonify({'modules': modules})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/module/<module_name>')
def get_module_detail(module_name):
    try:
        main_corr = data_cache['wgcna']['main_corr']
        main_pval = data_cache['wgcna']['main_pval']
        test_corr = data_cache['wgcna']['test_corr']
        test_pval = data_cache['wgcna']['test_pval']
        mod_asgn  = data_cache['wgcna']['module_assignments']
        main_row  = main_corr[main_corr['Module'] == module_name]
        test_row  = test_corr[test_corr['Module'] == module_name] if not test_corr.empty else pd.DataFrame()
        if main_row.empty:
            return jsonify({'error': 'Module not found'}), 404
        main_corr_row = main_row.iloc[0]
        main_pval_row = main_pval[main_pval['Module'] == module_name].iloc[0] if not main_pval.empty and (main_pval['Module'] == module_name).any() else None
        test_corr_row = test_row.iloc[0] if not test_row.empty else None
        test_pval_row = test_pval[test_pval['Module'] == module_name].iloc[0] if test_corr_row is not None and not test_pval.empty and (test_pval['Module'] == module_name).any() else None
        main_cols = [c for c in main_corr.columns if c != 'Module']
        test_cols = [c for c in test_corr.columns  if c != 'Module'] if test_corr_row is not None else []
        train_traits = _parse_heatmap_row(main_corr_row, main_pval_row, main_cols)
        test_traits  = _parse_heatmap_row(test_corr_row, test_pval_row, test_cols) if test_corr_row is not None else []
        prefix = module_name.split('_')[0]
        asgn_name = None
        for mn in mod_asgn['Module_Name'].unique():
            if mn.split('_')[0] == prefix:
                asgn_name = mn
                break
        proteins = mod_asgn[mod_asgn['Module_Name'] == asgn_name]['Protein'].tolist() if asgn_name else []
        color_src = asgn_name.split('_')[1] if asgn_name and '_' in asgn_name else module_name.split('_')[1] if '_' in module_name else module_name
        return jsonify({
            'module_name': module_name, 'color': color_src,
            'protein_count': len(proteins), 'proteins': proteins,
            'train_traits': train_traits, 'test_traits': test_traits,
        })
    except Exception as e:
        logger.error(f"Module detail error: {e}"); import traceback; traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@app.route('/api/ml_predict', methods=['POST'])
def ml_predict():
    try:
        data = request.get_json(force=True)
        protein_values = data.get('proteins', {})
        missing = [f for f in CORE_FEATURES if f not in protein_values]
        if missing:
            return jsonify({'error': f'Missing proteins: {missing}'}), 400
        model = data_cache['ml_model']
        if model is None:
            return jsonify({'error': 'ML model not loaded'}), 500
        x     = np.array([[protein_values[f] for f in CORE_FEATURES]], dtype=float)
        score = float(model.predict_proba(x)[0, 1])
        prediction = 'EOD' if score >= ML_THRESHOLD else 'CN'
        png_bytes = _build_shap_force_png(protein_values)
        import base64
        shap_b64 = base64.b64encode(png_bytes).decode('utf-8')
        return jsonify({'score': round(score,4), 'threshold': ML_THRESHOLD,
                        'prediction': prediction, 'shap_png': shap_b64})
    except Exception as e:
        logger.error(f"ML predict error: {e}"); import traceback; traceback.print_exc()
        return jsonify({'error': str(e)}), 500


# Load data once at module level (works for both direct run and gunicorn)
logger.info("Loading data...")
if not load_data():
    logger.error("Data load failed!")
else:
    logger.info("Data loaded, app ready.")

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000, use_reloader=False)
