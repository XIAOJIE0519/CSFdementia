// 国际化翻译
const translations = {
    zh: {
        'title': 'CSF蛋白组学痴呆数据库',
        'subtitle': '脑脊液蛋白组学痴呆数据库 - 多队列整合分析平台',
        'tab-protein': 'Protein',
        'tab-module': 'Module',
        'tab-ml': 'Machine Learning',
        'tab-source': 'Source Data',
        'search-title': '蛋白质查询',
        'search-placeholder': 'APOE|P02649',
        'search-button': '搜索',
        'volcano-title': '差异表达火山图',
        'volcano-desc': '点击火山图查看该比较组的详细差异表达结果和富集分析（仅显示FDR < 0.05的显著结果）',
        'module-title': 'WGCNA 共识模块分析',
        'module-desc': '点击模块查看详细的性状相关性分析和模块内蛋白列表',
        'module-detail-title': '模块详情',
        'protein-count': '蛋白数量',
        'trait-correlation': '性状相关性分析',
        'module-proteins': '模块内蛋白列表（显示前5个，可滚动查看全部）',
        'protein-id': '蛋白质ID',
        'action': '操作',
        'ml-calculating': '计算中...',
        'ml-predict-failed': '预测失败',
        'ml-request-failed': '预测请求失败',
        'ml-score-label': 'EOD 风险评分',
        'ml-threshold': '阈值',
        'ml-verdict-eod-pre': '模型预测：',
        'ml-verdict-eod-label': '可能为早发性痴呆（EOD）',
        'ml-verdict-eod-post': '（评分 {score} ≥ 阈值 {thr}）。请结合临床综合判断，本模型仅供研究参考。',
        'ml-verdict-low-pre': '模型预测：',
        'ml-verdict-low-label': '低 EOD 风险',
        'ml-verdict-low-post': '（评分 {score} < 阈值 {thr}）。请结合临床综合判断，本模型仅供研究参考。',
        'ml-title': '早发性痴呆 CSF 蛋白风险评分',
        'ml-intro': '输入15个核心蛋白的标准化z-score值，输出EOD风险评分及SHAP解释图。阈值=0.55，LODO AUC=0.877，外部验证AUC=0.860。',
        'ml-lodo-auc': 'LODO AUC',
        'ml-ext-auc': '外部验证 AUC',
        'ml-core-proteins': '核心蛋白',
        'ml-threshold-label': '判断阈值',
        'ml-input-desc': '请输入各蛋白经 <strong>per-cohort StandardScaler</strong> 标准化的 z-score 值。',
        'ml-clear': '清空',
        'ml-demo': '示例数据',
        'ml-shap-title': 'SHAP 力图',
        'ml-shap-desc': '红色：增加风险；蓝色：降低风险。',
        'ml-desc': '输入15个核心CSF蛋白的标准化表达值（z-score），模型将输出早发性痴呆（EOD）风险评分及SHAP解释图。阈值=0.55（Youden指数优化）。模型：EasyEnsemble-LogisticRegression，训练LODO AUC=0.877，外部验证AUC=0.860。',
        'ml-accuracy': '模型准确率',
        'ml-auc': '训练集 LODO AUC',
        'ml-features': '核心蛋白特征',
        'ml-input-title': '输入蛋白表达值',
        'ml-predict': '计算风险评分',
        'ml-result-title': '预测结果',
        'ml-shap': 'SHAP特征重要性',
        'ml-developing': '功能开发中',
        'ml-coming': '机器学习模型和SHAP分析结果即将上线',
        'total-proteins': '总蛋白数',
        'sig-proteins': '显著差异蛋白 (FDR<0.05)',
        'volcano-chart': '火山图',
        'volcano-chart-title': '火山图',
        'sig-proteins-table': '显著差异蛋白 (FDR < 0.05, 显示前5个，可滚动查看全部)',
        'protein': '蛋白质',
        'p-value': 'P值',
        'n-studies': '研究数量',
        'enrichment-table': '富集分析结果（显示前5个，可滚动查看全部）',
        'source': '来源',
        'pathway': '通路/功能',
        'gene-ratio': '基因比例',
        'gene-count': '基因数量',
        'view': '查看',
        'proteins': '个蛋白',
        'click-to-view': '点击查看火山图和差异表达结果',
        'back-home': '返回首页',
        'protein-info-title': '蛋白质信息',
        'expression-trend-title': '差异表达趋势（各研究）',
        'meta-analysis-title': 'Meta分析结果',
        'study-expression-title': '各研究中的表达情况',
        'correlation-title': '蛋白质相关性分析（CN / EOD / LOD 各前10个）',
        'module-label': '所属模块',
        'comparison': '比较组',
        'study': '研究',
        'weighted-effect': '加权效应值',
        'correlated-protein': '相关蛋白',
        'correlation-coef': '相关系数',
        'comp-EOAD_vs_CN': '早发性阿尔茨海默病 vs 认知正常',
        'comp-LOAD_vs_CN': '晚发性阿尔茨海默病 vs 认知正常',
        'comp-EOD_vs_CN': '早发性痴呆 vs 认知正常',
        'comp-LOD_vs_CN': '晚发性痴呆 vs 认知正常',
        'comp-EOAD_vs_LOAD': '早发性 vs 晚发性阿尔茨海默病',
        'comp-EOD_vs_LOD': '早发性 vs 晚发性痴呆',
        'source-title': '数据与代码开放获取',
        'source-warning-title': '重要声明',
        'source-warning-1': '请勿爬取本网站数据！本网站不提供API接口，禁止任何形式的数据爬取和网站攻击行为。',
        'source-warning-2': '所有数据和代码已在GitHub开源，请通过正规渠道获取。',
        'source-github-title': 'GitHub 仓库',
        'source-github-desc': '本项目的所有资源均已在GitHub开源，包括：',
        'source-item-1': '网站源代码：完整的前端和后端代码',
        'source-item-2': '分析脚本：数据处理、统计分析、可视化的完整代码',
        'source-item-3': '原始数据：最大限度可公开的原始数据文件',
        'source-item-4': '文档说明：详细的使用文档和数据说明',
        'source-visit-github': '访问 GitHub 仓库',
        'source-contact-title': '联系我们',
        'source-contact-desc': '如有任何疑问、建议或合作意向，欢迎联系：',
        'source-contact-person': '联系人：',
        'source-contact-email': '邮箱：',
        'source-terms-title': '使用条款',
        'source-terms-desc': '使用本网站和相关数据时，请遵守以下条款：',
        'source-term-1': '引用数据时请注明来源',
        'source-term-2': '禁止用于商业用途（如需商业使用请联系作者）',
        'source-term-3': '禁止爬取网站数据，请从GitHub获取',
        'source-term-4': '禁止对网站进行任何形式的攻击',
        'source-term-5': '遵守学术诚信和数据使用规范'
    },
    en: {
        'title': 'CSF Proteomics Dementia Database',
        'subtitle': 'CSF Proteomics Dementia Database - Multi-cohort Integrated Analysis Platform',
        'tab-protein': 'Protein',
        'tab-module': 'Module',
        'tab-ml': 'Machine Learning',
        'search-title': 'Protein Query',
        'search-placeholder': 'APOE|P02649',
        'search-button': 'Search',
        'volcano-title': 'Differential Expression Volcano Plots',
        'volcano-desc': 'Click volcano plot to view detailed differential expression results and enrichment analysis (only showing FDR < 0.05)',
        'module-title': 'WGCNA Consensus Module Analysis',
        'module-desc': 'Click module to view detailed trait correlation analysis and protein list',
        'module-detail-title': 'Module Details',
        'protein-count': 'Protein Count',
        'trait-correlation': 'Trait Correlation Analysis',
        'module-proteins': 'Module Proteins (showing first 5, scroll to view all)',
        'protein-id': 'Protein ID',
        'action': 'Action',
        'ml-calculating': 'Calculating...',
        'ml-predict-failed': 'Prediction failed',
        'ml-request-failed': 'Prediction request failed',
        'ml-score-label': 'EOD Risk Score',
        'ml-threshold': 'Threshold',
        'ml-verdict-eod-pre': 'Model predicts: ',
        'ml-verdict-eod-label': 'Likely Early-Onset Dementia (EOD)',
        'ml-verdict-eod-post': ' (score {score} >= threshold {thr}). Please combine with clinical judgment. This model is for research purposes only.',
        'ml-verdict-low-pre': 'Model predicts: ',
        'ml-verdict-low-label': 'Low EOD Risk',
        'ml-verdict-low-post': ' (score {score} < threshold {thr}). Please combine with clinical judgment. This model is for research purposes only.',
        'ml-title': 'Early-Onset Dementia CSF Protein Risk Score',
        'ml-intro': 'Enter standardized z-score values for 15 core CSF proteins. The model outputs an EOD risk score and SHAP force plot. Threshold=0.55, LODO AUC=0.877, external validation AUC=0.860.',
        'ml-lodo-auc': 'LODO AUC',
        'ml-ext-auc': 'External Validation AUC',
        'ml-core-proteins': 'Core Proteins',
        'ml-threshold-label': 'Decision Threshold',
        'ml-input-desc': 'Enter per-cohort <strong>StandardScaler</strong> normalized z-score values for each protein.',
        'ml-clear': 'Clear',
        'ml-demo': 'Demo Data',
        'ml-shap-title': 'SHAP Force Plot',
        'ml-shap-desc': 'Red: increases risk; Blue: decreases risk.',
        'ml-desc': 'Enter standardized z-score values for 15 core CSF proteins. The model outputs an Early-Onset Dementia (EOD) risk score and SHAP force plot. Threshold=0.55 (Youden optimized). Model: EasyEnsemble-LR, LODO AUC=0.877, external validation AUC=0.860.',
        'ml-accuracy': 'Model Accuracy',
        'ml-auc': 'Training LODO AUC',
        'ml-features': 'Core Protein Features',
        'ml-input-title': 'Enter Protein Expression Values',
        'ml-predict': 'Calculate Risk Score',
        'ml-result-title': 'Prediction Result',
        'ml-shap': 'SHAP Feature Importance',
        'ml-developing': 'Under Development',
        'ml-coming': 'Machine learning models and SHAP analysis results coming soon',
        'total-proteins': 'Total Proteins',
        'sig-proteins': 'Significant Proteins (FDR<0.05)',
        'volcano-chart': 'Volcano Plot',
        'volcano-chart-title': 'Volcano Plot',
        'sig-proteins-table': 'Significant Proteins (FDR < 0.05, showing first 5, scroll to view all)',
        'protein': 'Protein',
        'p-value': 'P-value',
        'n-studies': 'N Studies',
        'enrichment-table': 'Enrichment Analysis (showing first 5, scroll to view all)',
        'source': 'Source',
        'pathway': 'Pathway/Function',
        'gene-ratio': 'Gene Ratio',
        'gene-count': 'Gene Count',
        'view': 'View',
        'proteins': ' Proteins',
        'click-to-view': 'Click to view volcano plot and differential expression results',
        'back-home': 'Back to Home',
        'protein-info-title': 'Protein Information',
        'expression-trend-title': 'Differential Expression Trends (All Studies)',
        'meta-analysis-title': 'Meta-Analysis Results',
        'study-expression-title': 'Expression in Individual Studies',
        'correlation-title': 'Protein Correlation Analysis (CN / EOD / LOD, Top 10 Each)',
        'module-label': 'Module',
        'comparison': 'Comparison',
        'study': 'Study',
        'weighted-effect': 'Weighted Effect',
        'correlated-protein': 'Correlated Protein',
        'correlation-coef': 'Correlation Coefficient',
        'comp-EOAD_vs_CN': 'Early-Onset AD vs Cognitively Normal',
        'comp-LOAD_vs_CN': 'Late-Onset AD vs Cognitively Normal',
        'comp-EOD_vs_CN': 'Early-Onset Dementia vs Cognitively Normal',
        'comp-LOD_vs_CN': 'Late-Onset Dementia vs Cognitively Normal',
        'comp-EOAD_vs_LOAD': 'Early-Onset vs Late-Onset AD',
        'comp-EOD_vs_LOD': 'Early-Onset vs Late-Onset Dementia',
        'source-title': 'Data & Code Open Access',
        'source-warning-title': 'Important Notice',
        'source-warning-1': 'Do NOT scrape this website! This website does not provide API. Any form of data scraping or website attacks is strictly prohibited.',
        'source-warning-2': 'All data and code are open-sourced on GitHub. Please obtain them through official channels.',
        'source-github-title': 'GitHub Repository',
        'source-github-desc': 'All resources of this project are open-sourced on GitHub, including:',
        'source-item-1': 'Website Source Code: Complete frontend and backend code',
        'source-item-2': 'Analysis Scripts: Complete code for data processing, statistical analysis, and visualization',
        'source-item-3': 'Raw Data: Maximum publicly available raw data files',
        'source-item-4': 'Documentation: Detailed usage documentation and data description',
        'source-visit-github': 'Visit GitHub Repository',
        'source-contact-title': 'Contact Us',
        'source-contact-desc': 'For any questions, suggestions, or collaboration inquiries, please contact:',
        'source-contact-person': 'Contact:',
        'source-contact-email': 'Email:',
        'source-terms-title': 'Terms of Use',
        'source-terms-desc': 'When using this website and related data, please comply with the following terms:',
        'source-term-1': 'Cite the source when using the data',
        'source-term-2': 'Commercial use is prohibited (contact the author for commercial use)',
        'source-term-3': 'Do not scrape website data, obtain from GitHub instead',
        'source-term-4': 'Do not attack the website in any form',
        'source-term-5': 'Follow academic integrity and data usage standards'
    }
};

// 当前语言
let currentLanguage = 'zh';

// 获取翻译文本
function t(key) {
    return translations[currentLanguage][key] || key;
}

// 切换语言
function switchLanguage(lang) {
    currentLanguage = lang;
    
    // 更新按钮状态
    document.querySelectorAll('.language-switcher button').forEach(btn => {
        btn.classList.remove('active');
    });
    document.getElementById(`lang-${lang}`).classList.add('active');
    
    // 更新页面文本
    document.querySelectorAll('[data-i18n]').forEach(element => {
        const key = element.getAttribute('data-i18n');
        if (translations[lang][key]) {
            // Use innerHTML to support embedded HTML tags in translations
            element.innerHTML = translations[lang][key];
        }
    });
    
    // 更新placeholder
    document.querySelectorAll('[data-i18n-placeholder]').forEach(element => {
        const key = element.getAttribute('data-i18n-placeholder');
        if (translations[lang][key]) {
            element.placeholder = translations[lang][key];
        }
    });
    
    // 更新HTML lang属性
    document.documentElement.lang = lang === 'zh' ? 'zh-CN' : 'en';
    
    // 保存语言偏好
    localStorage.setItem('language', lang);
    
    // 触发自定义事件，通知其他组件语言已更改
    window.dispatchEvent(new CustomEvent('languageChanged', { detail: { language: lang } }));
}

// 页面加载时恢复语言设置
document.addEventListener('DOMContentLoaded', function() {
    const savedLang = localStorage.getItem('language') || 'zh';
    if (savedLang !== 'zh') {
        switchLanguage(savedLang);
    }
});
