#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
重新查询study 7, 8, 9的UniProt数据
使用size=5获取更多结果，观察一致性
"""

import requests
import pandas as pd
import time
import os

def query_uniprot_multiple(query_id, size=5):
    """
    查询UniProt数据库获取多个结果
    
    参数:
        query_id: UniProt ID或基因名
        size: 返回结果数量
    
    返回:
        list: 包含多个蛋白信息的列表
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    # 只强制人类物种
    params = {
        'query': f"{query_id} AND (organism_id:9606)",
        'format': 'json',
        'size': size
    }
    
    try:
        print(f"  查询: {query_id}...", end='')
        response = requests.get(base_url, params=params, timeout=10)
        response.raise_for_status()
        
        data = response.json()
        results = []
        
        if data.get('results') and len(data['results']) > 0:
            print(f" 找到 {len(data['results'])} 个结果")
            
            for idx, result in enumerate(data['results']):
                # 提取信息
                uniprot_id = result.get('primaryAccession', 'NA')
                
                # 获取基因名
                genes = result.get('genes', [])
                gene_name = 'NA'
                all_gene_names = []
                if genes and len(genes) > 0:
                    gene_name = genes[0].get('geneName', {}).get('value', 'NA')
                    all_gene_names.append(gene_name)
                    synonyms = genes[0].get('synonyms', [])
                    for syn in synonyms:
                        if 'value' in syn:
                            all_gene_names.append(syn['value'])
                
                # 获取蛋白名
                protein_name = result.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'NA')
                if protein_name == 'NA':
                    submitted_names = result.get('proteinDescription', {}).get('submissionNames', [])
                    if submitted_names:
                        protein_name = submitted_names[0].get('fullName', {}).get('value', 'NA')
                
                # 检查是否reviewed
                entry_type = result.get('entryType', 'Unknown')
                is_reviewed = 'Yes' if entry_type == 'UniProtKB reviewed (Swiss-Prot)' else 'No'
                
                results.append({
                    'Result_Rank': idx + 1,
                    'UniProt_ID': uniprot_id,
                    'Gene_Name': gene_name,
                    'All_Gene_Names': ';'.join(all_gene_names) if all_gene_names else 'NA',
                    'Protein_Name': protein_name,
                    'UniProtKB_Reviewed': is_reviewed
                })
        else:
            print(f" 未找到")
            results.append({
                'Result_Rank': 0,
                'UniProt_ID': 'NA',
                'Gene_Name': 'NA',
                'All_Gene_Names': 'NA',
                'Protein_Name': 'NA',
                'UniProtKB_Reviewed': 'NA'
            })
            
    except Exception as e:
        print(f" 错误: {str(e)}")
        results.append({
            'Result_Rank': 0,
            'UniProt_ID': 'NA',
            'Gene_Name': 'NA',
            'All_Gene_Names': 'NA',
            'Protein_Name': 'NA',
            'UniProtKB_Reviewed': f'Error: {str(e)}'
        })
    
    return results

def process_mapping_file(mapping_file, study_id, output_dir, test_mode=True):
    """
    处理单个mapping文件
    
    参数:
        mapping_file: mapping文件路径
        study_id: 研究编号
        output_dir: 输出目录
        test_mode: 是否测试模式（只处理前10个）
    """
    print(f"\n{'='*80}")
    print(f"处理 Study {study_id}: {mapping_file}")
    print(f"{'='*80}")
    
    # 读取原始mapping文件
    df = pd.read_csv(mapping_file, encoding='utf-8-sig')
    print(f"文件包含 {len(df)} 条记录")
    
    if test_mode:
        df = df.head(10)
        print(f"测试模式：只处理前 {len(df)} 条记录")
    
    # 存储新的查询结果
    all_results = []
    
    # 处理每条记录
    for idx, row in df.iterrows():
        query_id = row['Query']
        original_column = row['Original_Column']
        original_uniprot = row['UniProt_ID']
        original_gene = row['Gene_Name']
        
        print(f"\n[{idx+1}/{len(df)}] 原始Query: {query_id}")
        print(f"  原始结果: UniProt={original_uniprot}, Gene={original_gene}")
        
        # 查询获取多个结果
        results = query_uniprot_multiple(query_id, size=5)
        
        # 检查一致性
        first_result = results[0] if results else None
        is_consistent = False
        if first_result and first_result['UniProt_ID'] != 'NA':
            # 检查第一个结果是否与原始结果一致
            if first_result['UniProt_ID'] == original_uniprot:
                is_consistent = True
                print(f"  [OK] 一致: 第一个结果与原始结果相同")
            else:
                print(f"  [DIFF] 不一致: 第一个结果 {first_result['UniProt_ID']} != 原始 {original_uniprot}")
        
        # 添加到结果列表
        for result in results:
            result['Query'] = query_id
            result['Study'] = study_id
            result['Original_Column'] = original_column
            result['Original_UniProt_ID'] = original_uniprot
            result['Original_Gene_Name'] = original_gene
            result['Is_Consistent'] = 'Yes' if (is_consistent and result['Result_Rank'] == 1) else 'No'
            all_results.append(result)
        
        time.sleep(0.3)  # 避免请求过快
    
    # 保存结果
    results_df = pd.DataFrame(all_results)
    
    # 重新排列列顺序
    columns_order = [
        'Query', 'UniProt_ID', 'Gene_Name', 'All_Gene_Names', 
        'Protein_Name', 'UniProtKB_Reviewed', 'Study', 'Original_Column',
        'Result_Rank', 'Original_UniProt_ID', 'Original_Gene_Name', 'Is_Consistent'
    ]
    results_df = results_df[columns_order]
    
    output_file = os.path.join(output_dir, f'new_study_{study_id}_uniprot_mapping.csv')
    results_df.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"\n新查询结果已保存: {output_file}")
    
    # 统计一致性
    consistent_count = len(results_df[results_df['Is_Consistent'] == 'Yes'])
    total_queries = len(df)
    print(f"\n一致性统计: {consistent_count}/{total_queries} ({consistent_count/total_queries*100:.1f}%)")
    
    return results_df

def main():
    # 定义输入输出目录
    input_dir = 'F:/1a-EOD-CSF-protein/uniport'
    output_dir = 'F:/1a-EOD-CSF-protein/uniport'
    
    # 要处理的研究
    studies = ['7', '8', '9']
    
    # 测试模式：只处理前10个
    test_mode = False
    
    all_results = []
    
    for study_id in studies:
        mapping_file = os.path.join(input_dir, f'study_{study_id}_uniprot_mapping.csv')
        
        if not os.path.exists(mapping_file):
            print(f"\n警告: 文件不存在 - {mapping_file}")
            continue
        
        try:
            results = process_mapping_file(mapping_file, study_id, output_dir, test_mode=test_mode)
            all_results.append(results)
        except Exception as e:
            print(f"\n错误: 处理 Study {study_id} 时出错: {str(e)}")
            import traceback
            traceback.print_exc()
            continue
    
    print(f"\n{'='*80}")
    print("处理完成！")
    print(f"{'='*80}")

if __name__ == '__main__':
    main()
