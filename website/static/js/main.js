// CSF Proteomics Dementia Database - Main JavaScript

let currentComparison = null;
let modulesData = [];
let searchTimeout = null;

document.addEventListener('DOMContentLoaded', function() {
    initializeTabs();
    loadVolcanoPlots();
    loadModules();
    setupSearch();
    initMLInputs();
    setTimeout(() => { document.body.classList.add('loaded'); }, 1500);
});

function initializeTabs() {
    document.querySelectorAll('.tab').forEach(tab => {
        tab.addEventListener('click', function() { switchTab(this.getAttribute('data-tab')); });
    });
}
function switchTab(tabName) {
    document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
    document.querySelector(`[data-tab="${tabName}"]`).classList.add('active');
    document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
    document.getElementById(`${tabName}-tab`).classList.add('active');
}
function setupSearch() {
    const input = document.getElementById('protein-search');
    input.addEventListener('input', function(e) {
        clearTimeout(searchTimeout);
        const q = e.target.value.trim();
        if (q.length < 2) { hideSuggestions(); return; }
        searchTimeout = setTimeout(() => fetchSuggestions(q), 300);
    });
    input.addEventListener('keypress', function(e) {
        if (e.key === 'Enter') { hideSuggestions(); searchProtein(); }
    });
    document.addEventListener('click', function(e) {
        if (!e.target.closest('.search-box')) hideSuggestions();
    });
}
async function fetchSuggestions(query) {
    try {
        const r = await fetch(`/api/search_protein?q=${encodeURIComponent(query)}`);
        const d = await r.json();
        if (d.proteins && d.proteins.length > 0) showSuggestions(d.proteins);
        else hideSuggestions();
    } catch(e) { console.error('Suggestions:', e); }
}
function showSuggestions(proteins) {
    const div = document.getElementById('search-suggestions');
    div.innerHTML = '';
    proteins.forEach(p => {
        const item = document.createElement('div');
        item.className = 'suggestion-item';
        item.textContent = p;
        item.onclick = () => { document.getElementById('protein-search').value = p; hideSuggestions(); searchProtein(); };
        div.appendChild(item);
    });
    div.style.display = 'block';
}
function hideSuggestions() { document.getElementById('search-suggestions').style.display = 'none'; }
function searchProtein() {
    const q = document.getElementById('protein-search').value.trim();
    if (!q) { alert('请输入蛋白质名称'); return; }
    window.location.href = `/protein/${encodeURIComponent(q)}`;
}
function searchProteinById(id) { window.location.href = `/protein/${encodeURIComponent(id)}`; }

async function loadVolcanoPlots() {
    try {
        const r = await fetch('/api/comparisons');
        const d = await r.json();
        const grid = document.getElementById('volcano-grid');
        grid.innerHTML = '';
        d.comparisons.forEach(comp => {
            const card = document.createElement('div');
            card.className = 'volcano-card';
            card.onclick = () => showDifferentialExpression(comp.id);
            card.innerHTML = `<h3>${comp.name}</h3><p style="color:var(--text-secondary);font-size:0.9rem;margin-bottom:1rem;">${t('comp-'+comp.id)}</p><div class="volcano-placeholder">${t('click-to-view')}</div>`;
            grid.appendChild(card);
        });
    } catch(e) { console.error('loadVolcanoPlots:', e); }
}
window.addEventListener('languageChanged', function() {
    loadVolcanoPlots();
    if (modulesData.length > 0) loadModules();
});

async function showDifferentialExpression(comparison) {
    currentComparison = comparison;
    try {
        const r = await fetch(`/api/differential_expression/${comparison}`);
        const data = await r.json();
        document.getElementById('de-modal').style.display = 'block';
        const names = { zh: {
            'EOAD_vs_CN':'EOAD vs CN - 早发性阿尔茨海默病 vs 认知正常',
            'LOAD_vs_CN':'LOAD vs CN - 晚发性阿尔茨海默病 vs 认知正常',
            'EOD_vs_CN':'EOD vs CN - 早发性痴呆 vs 认知正常',
            'LOD_vs_CN':'LOD vs CN - 晚发性痴呆 vs 认知正常',
            'EOAD_vs_LOAD':'EOAD vs LOAD - 早发性 vs 晚发性阿尔茨海默病',
            'EOD_vs_LOD':'EOD vs LOD - 早发性 vs 晚发性痴呆'
        }, en: {
            'EOAD_vs_CN':'EOAD vs CN - Early-Onset AD vs CN',
            'LOAD_vs_CN':'LOAD vs CN - Late-Onset AD vs CN',
            'EOD_vs_CN':'EOD vs CN - Early-Onset Dementia vs CN',
            'LOD_vs_CN':'LOD vs CN - Late-Onset Dementia vs CN',
            'EOAD_vs_LOAD':'EOAD vs LOAD - Early vs Late-Onset AD',
            'EOD_vs_LOD':'EOD vs LOD - Early vs Late-Onset Dementia'
        }};
        document.getElementById('de-modal-title').textContent = (names[currentLanguage]||names.en)[comparison]||comparison;
        document.getElementById('de-total-proteins').textContent = data.total_proteins;
        document.getElementById('de-sig-proteins').textContent  = data.significant_count;
        plotVolcano(data.volcano_data);
        const fmt = n => (n===null||n===undefined)?'N/A':Number(n.toPrecision(3));
        const fdrBadge = v => v!==null&&v<0.05?`<span class="badge" style="background:#10b981;">${fmt(v)}</span>`:fmt(v);
        const sigBody = document.getElementById('sig-protein-table-body');
        sigBody.innerHTML = '';
        (data.significant_proteins||[]).forEach(p => {
            const tr = document.createElement('tr');
            tr.innerHTML = `<td style="font-family:'JetBrains Mono',monospace;">${p.protein}</td><td>${fmt(p.log2fc)}</td><td>${fdrBadge(p.fdr)}</td><td>${p.n_studies||'N/A'}</td><td><a href="/protein/${encodeURIComponent(p.protein)}" class="external-link" style="font-size:0.8rem;">${t('view')}</a></td>`;
            sigBody.appendChild(tr);
        });
        const gokeggBody = document.getElementById('gokegg-table-body');
        const gseaBody   = document.getElementById('gsea-table-body');
        gokeggBody.innerHTML = '';
        gseaBody.innerHTML   = '';
        const fdrBadgeE = v => v!==null&&v<0.05?`<span class="badge" style="background:#10b981;">${fmt(v)}</span>`:fmt(v);
        if (data.enrichment_gokegg && data.enrichment_gokegg.length > 0) {
            data.enrichment_gokegg.forEach(term => {
                const tr = document.createElement('tr');
                const badge = `<span class="badge badge-info" style="font-size:0.7rem;">${term.source}</span>`;
                tr.innerHTML = `<td>${badge}</td><td><strong>${term.term_name}</strong><br><small style="color:var(--text-secondary);">${term.term_id}</small></td><td>${term.gene_ratio}</td><td>${fdrBadgeE(term.fdr)}</td><td>${term.count}</td>`;
                gokeggBody.appendChild(tr);
            });
        } else {
            gokeggBody.innerHTML = '<tr><td colspan="5" style="text-align:center;color:var(--text-secondary);">暂无GO/KEGG富集结果</td></tr>';
        }
        if (data.enrichment_gsea && data.enrichment_gsea.length > 0) {
            data.enrichment_gsea.forEach(term => {
                const tr = document.createElement('tr');
                const nes = term.nes !== null && term.nes !== undefined ? Number(term.nes.toPrecision(3)) : 'N/A';
                const nesColor = (typeof term.nes === 'number' && term.nes > 0) ? '#ef4444' : '#3b82f6';
                tr.innerHTML = `<td><strong>${term.term_name}</strong><br><small style="color:var(--text-secondary);">${term.term_id}</small></td><td style="font-weight:700;color:${nesColor};">${nes}</td><td>${fdrBadgeE(term.fdr)}</td><td>${term.count}</td>`;
                gseaBody.appendChild(tr);
            });
        } else {
            gseaBody.innerHTML = '<tr><td colspan="4" style="text-align:center;color:var(--text-secondary);">暂无GSEA富集结果</td></tr>';
        }
    } catch(e) { console.error('showDE:', e); alert('加载失败'); }
}
function plotVolcano(volcanoData) {
    const colors = getChartColors();
    let d = volcanoData;
    if (d.length > 600) d = d.slice().sort((a,b)=>a.fdr-b.fdr).slice(0,600);
    Plotly.newPlot('volcano-chart', [{
        x:d.map(p=>p.log2fc), y:d.map(p=>p.neg_log10_fdr),
        mode:'markers', type:'scattergl', text:d.map(p=>p.protein),
        marker:{size:7, color:d.map(p=>p.log2fc>0?'#ef4444':'#3b82f6'), opacity:0.75, line:{width:0.5,color:'white'}},
        hovertemplate:'<b>%{text}</b><br>Log2FC:%{x:.3f}<br>-log10(FDR):%{y:.3f}<extra></extra>'
    }], {
        title:{text:t('volcano-chart-title'),font:{size:16,color:colors.textColor,family:'Poppins'}},
        xaxis:{title:{text:'Log2 Fold Change'},zeroline:true,zerolinewidth:2,zerolinecolor:colors.gridColor,gridcolor:colors.gridColor,color:colors.textColor},
        yaxis:{title:{text:'-Log10(FDR)'},gridcolor:colors.gridColor,color:colors.textColor},
        hovermode:'closest',plot_bgcolor:colors.background,paper_bgcolor:colors.background,
        font:{color:colors.textColor,family:'Poppins'},margin:{t:50,r:30,b:50,l:60},autosize:true
    }, {responsive:true,displaylogo:false,modeBarButtonsToRemove:['lasso2d','select2d']});
}
window.addEventListener('themeChanged', function() {
    const m = document.getElementById('de-modal');
    if (m && m.style.display==='block' && currentComparison) showDifferentialExpression(currentComparison);
});
function closeDeModal() { document.getElementById('de-modal').style.display='none'; }
document.addEventListener('click', function(e) {
    const m = document.getElementById('de-modal');
    if (e.target===m) closeDeModal();
});

// ─── Modules ───
async function loadModules() {
    try {
        const r = await fetch('/api/modules');
        const data = await r.json();
        modulesData = data.modules;
        const grid = document.getElementById('module-grid');
        grid.innerHTML = '';
        const colorMap = {
            'black':'#000000','red':'#ef4444','lightcyan':'#22d3ee','turquoise':'#14b8a6',
            'midnightblue':'#1e3a8a','brown':'#92400e','magenta':'#db2777','blue':'#3b82f6',
            'green':'#10b981','lightgreen':'#84cc16','greenyellow':'#a3e635','cyan':'#06b6d4',
            'salmon':'#fb7185','pink':'#f472b6','purple':'#a855f7','grey60':'#6b7280',
            'tan':'#d97706','yellow':'#eab308'
        };
        data.modules.forEach(mod => {
            const card = document.createElement('div');
            card.className = 'module-card';
            card.style.setProperty('--module-color', colorMap[mod.color]||'#0ea5e9');
            card.onclick = () => showModuleDetail(mod.module_name);
            card.innerHTML = `<h4>${mod.module_name}</h4><p class="protein-count">${mod.protein_count}${t('proteins')}</p>`;
            grid.appendChild(card);
        });
    } catch(e) { console.error('loadModules:', e); }
}

async function showModuleDetail(moduleName) {
    try {
        const r = await fetch(`/api/module/${encodeURIComponent(moduleName)}`);
        const data = await r.json();
        const div = document.getElementById('module-detail');
        div.style.display = 'block';
        document.getElementById('module-name').textContent = moduleName;
        document.getElementById('module-protein-count').textContent = data.protein_count;
        plotTraitCorrelation(data.train_traits, data.test_traits);
        const tbody = document.getElementById('module-protein-table-body');
        tbody.innerHTML = '';
        (data.proteins||[]).forEach(protein => {
            const tr = document.createElement('tr');
            tr.innerHTML = `<td style="font-family:'JetBrains Mono',monospace;">${protein}</td><td><a href="/protein/${encodeURIComponent(protein)}" class="external-link" style="font-size:0.8rem;">${t('view')}</a></td>`;
            tbody.appendChild(tr);
        });
        div.scrollIntoView({behavior:'smooth'});
    } catch(e) { console.error('showModuleDetail:', e); alert('加载失败'); }
}

function plotTraitCorrelation(trainTraits, testTraits) {
    const container = document.getElementById('trait-correlation-container');
    container.innerHTML = '';
    const fmt = n => (n===null||n===undefined)?'N/A':Number(n.toPrecision(3));
    const getColor = corr => {
        if (corr===null||corr===undefined) return '#cbd5e1';
        const c = Math.max(-0.5, Math.min(0.5, corr)), tt = c+0.5;
        if (tt < 0.5) {
            const f=tt/0.5;
            return `rgb(${Math.round(71+f*176)},${Math.round(117+f*130)},${Math.round(174+f*73)})`;
        } else {
            const f=(tt-0.5)/0.5;
            return `rgb(${Math.round(247-f*34)},${Math.round(247-f*185)},${Math.round(247-f*168)})`;
        }
    };
    const getTC = bg => {
        const m=bg.match(/(\d+),(\d+),(\d+)/);
        if (!m) return '#0f172a';
        return (parseInt(m[1])*299+parseInt(m[2])*587+parseInt(m[3])*114)/1000<140?'#f1f5f9':'#0f172a';
    };
    const buildSection = (title, traits) => {
        const sec = document.createElement('div');
        sec.className = 'trait-section';
        sec.innerHTML = `<h4 class="trait-section-title">${title}</h4>`;
        const traitItems = (traits||[]).filter(t=>t.type==='trait_corr');
        const cellItems  = (traits||[]).filter(t=>t.type==='cell_type_ora');
        const consItems  = (traits||[]).filter(t=>t.type==='conservancy');
        if (traitItems.length > 0) {
            const grid = document.createElement('div');
            grid.className = 'trait-grid';
            traitItems.forEach(trait => {
                const bg=getColor(trait.correlation), tc=getTC(bg);
                const card = document.createElement('div');
                card.className = 'trait-card';
                card.style.cssText = `background:${bg};border-color:${bg};color:${tc};`;
                card.innerHTML = `
                    <div class="trait-card-header" title="${trait.trait}" style="color:${tc};font-weight:700;">${trait.trait}</div>
                    <div class="trait-card-correlation" style="color:${tc};font-size:1.8rem;font-weight:800;">${trait.correlation!==null?fmt(trait.correlation):'N/A'}</div>
                    <div class="trait-card-stats" style="border-top-color:${tc};">
                        <div class="trait-card-stat"><div class="trait-card-stat-label" style="color:${tc};">p</div><div style="color:${tc};">${fmt(trait.p_value)}</div></div>
                    </div>`;
                grid.appendChild(card);
            });
            sec.appendChild(grid);
        }
        if (cellItems.length > 0) {
            const h5 = document.createElement('h5');
            h5.style.cssText = 'margin:1rem 0 0.5rem;color:var(--text-secondary);font-size:0.85rem;font-weight:600;';
            h5.textContent = 'Cell-type ORA p-value';
            sec.appendChild(h5);
            const wrap = document.createElement('div');
            wrap.style.cssText = 'display:flex;flex-wrap:wrap;gap:0.4rem;';
            cellItems.forEach(ct => {
                const sig = ct.p_value!==null && ct.p_value<0.05;
                const chip = document.createElement('div');
                chip.className = 'cell-type-card';
                chip.style.borderColor = sig?'var(--primary-color)':'';
                chip.innerHTML = `<span class="ct-name">${ct.trait}</span><span class="ct-p" style="color:${sig?'var(--primary-color)':'var(--text-secondary)'};">p=${ct.p_value!==null?ct.p_value.toExponential(2):'N/A'}</span>`;
                wrap.appendChild(chip);
            });
            sec.appendChild(wrap);
        }
        if (consItems.length > 0) {
            const h5c = document.createElement('h5');
            h5c.style.cssText = 'margin:1rem 0 0.5rem;color:var(--text-secondary);font-size:0.85rem;font-weight:600;';
            h5c.textContent = 'Module Preservation (Zsummary)';
            sec.appendChild(h5c);
            const wrap2 = document.createElement('div');
            wrap2.style.cssText = 'display:flex;flex-wrap:wrap;gap:0.4rem;';
            consItems.forEach(cs => {
                const high = cs.n_studies!==null && cs.n_studies>=2;
                const chip2 = document.createElement('div');
                chip2.className = 'cell-type-card';
                chip2.style.borderColor = high?'var(--primary-color)':'';
                chip2.innerHTML = `<span class="ct-name">${cs.trait}</span><span class="ct-p" style="color:${high?'var(--primary-color)':'var(--text-secondary)'};">${cs.n_studies!==null?cs.n_studies.toFixed(2):'N/A'}</span>`;
                wrap2.appendChild(chip2);
            });
            sec.appendChild(wrap2);
        }
        return sec;
    };
    container.appendChild(buildSection('Main Dataset (Training)', trainTraits));
    if (testTraits && testTraits.length > 0)
        container.appendChild(buildSection('Apply Dataset (Test)', testTraits));
}

// ─── ML Prediction ───
const ML_FEATURES = ['MIF','DDAH1','ENO2','PEBP1','PAM','SPON1','SOD1','RTN4R','VASN','SOD2','GFRA2','CA4','CANT1','GLRX','GAS6'];
const ML_DEMO = {MIF:0.82,DDAH1:1.14,ENO2:0.67,PEBP1:-0.45,PAM:1.03,SPON1:0.91,SOD1:0.58,RTN4R:-0.32,VASN:0.74,SOD2:0.61,GFRA2:0.88,CA4:-0.21,CANT1:0.53,GLRX:0.79,GAS6:0.44};

function initMLInputs() {
    const container = document.getElementById('ml-protein-inputs');
    if (!container) return;
    container.innerHTML = '';
    ML_FEATURES.forEach(f => {
        const item = document.createElement('div');
        item.className = 'ml-protein-input-item';
        item.innerHTML = `<label for="ml-input-${f}">${f}</label><input type="number" id="ml-input-${f}" step="0.01" placeholder="z-score">`;
        container.appendChild(item);
    });
}
function clearMLInputs() {
    ML_FEATURES.forEach(f => { const el=document.getElementById(`ml-input-${f}`); if(el) el.value=''; });
    document.getElementById('ml-result-card').style.display='none';
}
function fillMLDemo() {
    ML_FEATURES.forEach(f => { const el=document.getElementById(`ml-input-${f}`); if(el) el.value=ML_DEMO[f]||0; });
}
async function runMLPrediction() {
    const proteins={};
    ML_FEATURES.forEach(f => { const v=parseFloat(document.getElementById(`ml-input-${f}`)?.value); proteins[f]=isNaN(v)?0:v; });
    const btn=document.getElementById('ml-predict-btn');
    btn.disabled=true;
    btn.innerHTML='<span style="display:flex;align-items:center;gap:0.5rem;"><span class="ml-loading-spinner"></span>'+t('ml-calculating')+'</span>';
    try {
        const r=await fetch('/api/ml_predict',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({proteins})});
        const data=await r.json();
        if (data.error){alert(t('ml-predict-failed')+': '+data.error);return;}
        const score=data.score, thr=data.threshold, isEOD=data.prediction==='EOD';
        const pct=Math.round(score*100), col=isEOD?'#ef4444':'#10b981';
        document.getElementById('ml-score-display').innerHTML=`
            <div><div style="font-size:0.8rem;color:var(--text-secondary);margin-bottom:0.3rem;">${t('ml-score-label')}</div><div class="ml-score-value" style="color:${col};">${score.toFixed(3)}</div></div>
            <div class="ml-score-bar-wrap"><div style="font-size:0.75rem;color:var(--text-secondary);margin-bottom:0.4rem;">0 ──── ${t('ml-threshold')} ${thr} ──── 1</div>
            <div class="ml-score-bar-bg"><div class="ml-score-bar-fill" style="width:${pct}%;background:${col};"></div><div class="ml-score-bar-threshold" style="left:${thr*100}%;"></div></div></div>`;
        const vEl=document.getElementById('ml-verdict');
        vEl.className='ml-verdict '+(isEOD?'positive':'negative');
        vEl.innerHTML=isEOD
            ?`⚠️ ${t('ml-verdict-eod-pre')}<strong>${t('ml-verdict-eod-label')}</strong>${t('ml-verdict-eod-post').replace('{score}',score.toFixed(3)).replace('{thr}',thr)}`
            :`✓ ${t('ml-verdict-low-pre')}<strong>${t('ml-verdict-low-label')}</strong>${t('ml-verdict-low-post').replace('{score}',score.toFixed(3)).replace('{thr}',thr)}`;
        document.getElementById('ml-shap-img').src=`data:image/png;base64,${data.shap_png}`;
        const rc=document.getElementById('ml-result-card');
        rc.style.display='block';
        rc.scrollIntoView({behavior:'smooth',block:'start'});
    } catch(e){console.error('ML predict:',e);alert(t('ml-request-failed'));}
    finally{btn.disabled=false;btn.innerHTML='<span data-i18n="ml-predict">'+t('ml-predict')+'</span>';}
}
