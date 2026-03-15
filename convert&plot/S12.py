#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Study 12 单独分析
1. 使用ID_mapping_final.csv进行列名映射
2. 未找到的通过UniProt API查询
3. API也找不到的删除
4. 合并重复测量（多列映射到同一基因取平均）
5. 进行差异分析：CN vs EOND, CN vs LOND, EOND vs LOND
"""

import requests
import pandas as pd
import time
import os
import subprocess
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

def query_uniprot(query_id, is_human=True):
    """
    查询UniProt数据库获取蛋白信息
    
    参数:
        query_id: UniProt ID
        is_human: 是否限制为人类
    
    返回:
        dict: 包含蛋白信息的字典
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    query_parts = [query_id]
    if is_human:
        query_parts.append("(organism_id:9606)")
    
    params = {
        'query': ' AND '.join(query_parts),
        'format': 'json',
        'size': 1
    }
    
    try:
        response = requests.get(base_url, params=params, timeout=10)
        response.raise_for_status()
        
        data = response.json()
        
        if data.get('results') and len(data['results']) > 0:
            result = data['results'][0]
            
            uniprot_id = result.get('primaryAccession', 'NA')
            
            genes = result.get('genes', [])
            gene_name = 'NA'
            if genes and len(genes) > 0:
                gene_name = genes[0].get('geneName', {}).get('value', 'NA')
            
            protein_name = result.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'NA')
            if protein_name == 'NA':
                submitted_names = result.get('proteinDescription', {}).get('submissionNames', [])
                if submitted_names:
                    protein_name = submitted_names[0].get('fullName', {}).get('value', 'NA')
            
            entry_type = result.get('entryType', 'Unknown')
            is_reviewed_status = 'Yes' if entry_type == 'UniProtKB reviewed (Swiss-Prot)' else 'No'
            
            return {
                'Query': query_id,
                'UniProt_ID': uniprot_id,
                'Gene_Name': gene_name,
                'Protein_Name': protein_name,
                'UniProtKB_Reviewed': is_reviewed_status
            }
        else:
            return {
                'Query': query_id,
                'UniProt_ID': 'NA',
                'Gene_Name': 'NA',
                'Protein_Name': 'NA',
                'UniProtKB_Reviewed': 'NA'
            }
            
    except Exception as e:
        return {
            'Query': query_id,
            'UniProt_ID': 'NA',
            'Gene_Name': 'NA',
            'Protein_Name': 'NA',
            'UniProtKB_Reviewed': f'Error: {str(e)}'
        }

def main():
    print("="*80)
    print("Study 12 分析")
    print("="*80)
    
    # 文件路径
    input_file = 'F:/1a-EOD-CSF-protein/S12/study_12.csv'
    output_dir = 'F:/1a-EOD-CSF-protein/S12'
    mapping_file = 'F:/1a-EOD-CSF-protein/uniport/ID_mapping_final.csv'
    
    # 读取数据
    print(f"\n读取数据: {input_file}")
    df = pd.read_csv(input_file)
    print(f"数据维度: {df.shape}")
    print(f"分组信息:\n{df['Group_New'].value_counts()}")
    
    # 获取蛋白列
    protein_columns = [col for col in df.columns if col not in ['SampleID', 'Group_New']]
    print(f"\n蛋白列数: {len(protein_columns)}")
    
    # 步骤1: 从ID_mapping_final.csv加载映射
    print("\n" + "="*80)
    print("步骤1: 从ID_mapping_final.csv加载映射")
    print("="*80)
    
    mapping_df = pd.read_csv(mapping_file, encoding='utf-8-sig')
    print(f"映射文件记录数: {len(mapping_df)}")
    
    # 构建映射字典（所有研究的映射）
    existing_mapping = {}
    for _, row in mapping_df.iterrows():
        original_col = row['Original_Column']
        gene_name = row['Gene_Name']
        if pd.notna(gene_name) and gene_name != 'NA':
            existing_mapping[original_col] = gene_name
    
    print(f"已有映射数: {len(existing_mapping)}")
    
    # 步骤2: 检查哪些列已有映射，哪些需要查询
    print("\n" + "="*80)
    print("步骤2: 检查映射覆盖情况")
    print("="*80)
    
    mapped_cols = []
    unmapped_cols = []
    
    for col in protein_columns:
        if col in existing_mapping:
            mapped_cols.append(col)
        else:
            unmapped_cols.append(col)
    
    print(f"已有映射的列: {len(mapped_cols)}")
    print(f"需要查询的列: {len(unmapped_cols)}")
    
    # 步骤3: 对未映射的列进行UniProt查询
    if unmapped_cols:
        print("\n" + "="*80)
        print("步骤3: 查询未映射的蛋白")
        print("="*80)
        
        new_mappings = []
        for i, col in enumerate(unmapped_cols, 1):
            print(f"[{i}/{len(unmapped_cols)}] 查询: {col}...", end='')
            result = query_uniprot(col, is_human=True)
            
            if result['Gene_Name'] != 'NA':
                print(f" 找到: {result['Gene_Name']}")
                existing_mapping[col] = result['Gene_Name']
                new_mappings.append({
                    'Original_Column': col,
                    'Gene_Name': result['Gene_Name'],
                    'UniProt_ID': result['UniProt_ID'],
                    'Protein_Name': result['Protein_Name'],
                    'Source': 'API_Query'
                })
            else:
                print(f" 未找到")
            
            time.sleep(0.2)
        
        # 保存新查询的映射
        if new_mappings:
            new_mapping_df = pd.DataFrame(new_mappings)
            new_mapping_file = os.path.join(output_dir, 'study_12_new_mappings.csv')
            new_mapping_df.to_csv(new_mapping_file, index=False, encoding='utf-8-sig')
            print(f"\n新查询的映射已保存: {new_mapping_file}")
            print(f"新增映射数: {len(new_mappings)}")
    
    # 步骤4: 构建基因名到原始列的映射（处理重复）
    print("\n" + "="*80)
    print("步骤4: 构建清洗后的数据")
    print("="*80)
    
    gene_to_cols = defaultdict(list)
    deleted_cols = []
    
    for col in protein_columns:
        if col in existing_mapping:
            gene_name = existing_mapping[col]
            gene_to_cols[gene_name].append(col)
        else:
            deleted_cols.append(col)
    
    print(f"成功映射的列: {len(protein_columns) - len(deleted_cols)}")
    print(f"删除的列（无法映射）: {len(deleted_cols)}")
    print(f"唯一基因数: {len(gene_to_cols)}")
    
    # 检查重复映射
    merged_genes = [(gene, cols) for gene, cols in gene_to_cols.items() if len(cols) > 1]
    print(f"需要合并的基因（多列取平均）: {len(merged_genes)}")
    
    if merged_genes:
        print("\n前5个需要合并的基因:")
        for gene, cols in merged_genes[:5]:
            print(f"  {gene} ← {len(cols)}列: {', '.join(cols[:3])}{' ...' if len(cols) > 3 else ''}")
    
    # 构建新数据框
    new_df = pd.DataFrame()
    new_df['SampleID'] = df['SampleID']
    
    for gene_name, orig_cols in gene_to_cols.items():
        if len(orig_cols) == 1:
            new_df[gene_name] = df[orig_cols[0]]
        else:
            # 多列取平均
            cols_data = [df[col] for col in orig_cols]
            new_df[gene_name] = pd.concat(cols_data, axis=1).mean(axis=1)
    
    new_df['Group_New'] = df['Group_New']
    
    print(f"\n清洗后数据维度: {new_df.shape}")
    
    # 保存清洗后的数据
    cleaned_file = os.path.join(output_dir, 'study_12_cleaned.csv')
    new_df.to_csv(cleaned_file, index=False, encoding='utf-8-sig')
    print(f"清洗后数据已保存: {cleaned_file}")
    
    # 步骤5: 创建R脚本进行差异分析
    print("\n" + "="*80)
    print("步骤5: 创建R脚本进行差异分析")
    print("="*80)
    
    r_script = '''
# Study 12 差异分析
library(limma)

cat("\\n", "="*80, "\\n")
cat("Study 12 差异分析\\n")
cat("="*80, "\\n")

# 读取数据
data <- read.csv("F:/1a-EOD-CSF-protein/S12/study_12_cleaned.csv", row.names = 1)

# 分离表达矩阵和分组信息
groups <- data$Group_New
expr_matrix <- data[, -ncol(data)]
expr_matrix <- t(expr_matrix)  # 转置：行为基因，列为样本

# 检查数据
cat("\\n表达矩阵维度:", dim(expr_matrix), "\\n")
cat("分组信息:\\n")
print(table(groups))

# 定义三个比较
comparisons <- list(
  list(name = "CN_vs_EOND", group1 = "CN", group2 = "EOND"),
  list(name = "CN_vs_LOND", group1 = "CN", group2 = "LOND"),
  list(name = "EOND_vs_LOND", group1 = "EOND", group2 = "LOND")
)

# 对每个比较进行差异分析
for (comp in comparisons) {
  cat("\\n", "="*80, "\\n")
  cat("分析:", comp$name, "\\n")
  cat("="*80, "\\n")
  
  # 选择两组样本
  idx <- groups %in% c(comp$group1, comp$group2)
  expr_sub <- expr_matrix[, idx]
  groups_sub <- factor(groups[idx], levels = c(comp$group1, comp$group2))
  
  cat("样本数:", sum(idx), "\\n")
  cat("分组:", table(groups_sub), "\\n")
  
  # 构建设计矩阵
  design <- model.matrix(~ 0 + groups_sub)
  colnames(design) <- levels(groups_sub)
  
  # 拟合线性模型
  fit <- lmFit(expr_sub, design)
  
  # 构建对比矩阵
  contrast_str <- paste0(comp$group2, "-", comp$group1)
  contrast.matrix <- makeContrasts(contrasts = contrast_str, levels = design)
  
  # 对比分析
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # 提取结果
  results <- topTable(fit2, number = Inf, sort.by = "P")
  
  # 添加统计信息
  results$Gene <- rownames(results)
  results <- results[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  colnames(results) <- c("Gene", "logFC", "AveExpr", "t_statistic", "P_value", "FDR", "B_statistic")
  
  # 添加显著性标记
  results$Significant <- ifelse(results$FDR < 0.05, "Yes", "No")
  
  # 保存结果
  output_file <- paste0("F:/1a-EOD-CSF-protein/S12/DE_", comp$name, ".csv")
  write.csv(results, output_file, row.names = FALSE)
  cat("结果已保存:", output_file, "\\n")
  
  # 统计显著差异基因
  sig_up <- sum(results$FDR < 0.05 & results$logFC > 0)
  sig_down <- sum(results$FDR < 0.05 & results$logFC < 0)
  cat("显著上调基因:", sig_up, "\\n")
  cat("显著下调基因:", sig_down, "\\n")
}

cat("\\n", "="*80, "\\n")
cat("差异分析完成！\\n")
cat("="*80, "\\n")
'''
    
    # 保存R脚本
    r_script_file = os.path.join(output_dir, 'differential_analysis.R')
    with open(r_script_file, 'w', encoding='utf-8') as f:
        f.write(r_script)
    print(f"R脚本已保存: {r_script_file}")
    
    # 运行R脚本
    print("\n" + "="*80)
    print("步骤6: 运行R脚本进行差异分析")
    print("="*80)
    
    r_path = 'E:/R-4.4.2/bin/x64/Rscript.exe'
    try:
        result = subprocess.run([r_path, r_script_file], 
                              capture_output=True, 
                              text=True, 
                              encoding='utf-8',
                              timeout=600)
        print(result.stdout)
        if result.stderr:
            print("警告/错误信息:")
            print(result.stderr)
        
        if result.returncode == 0:
            print("\n差异分析完成！")
        else:
            print(f"\nR脚本执行失败，返回码: {result.returncode}")
    except Exception as e:
        print(f"运行R脚本时出错: {str(e)}")
    
    # 生成汇总报告
    print("\n" + "="*80)
    print("汇总报告")
    print("="*80)
    print(f"原始蛋白列数: {len(protein_columns)}")
    print(f"从ID_mapping_final.csv映射: {len(mapped_cols)}")
    print(f"通过API新查询: {len(new_mappings) if unmapped_cols else 0}")
    print(f"无法映射（已删除）: {len(deleted_cols)}")
    print(f"最终基因数: {len(gene_to_cols)}")
    print(f"合并的基因数: {len(merged_genes)}")
    print(f"映射成功率: {(len(protein_columns) - len(deleted_cols)) / len(protein_columns) * 100:.2f}%")
    
    print("\n" + "="*80)
    print("所有分析完成！")
    print("="*80)

if __name__ == '__main__':
    main()
