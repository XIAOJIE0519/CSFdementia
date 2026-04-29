#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
批量查询UniProt数据库获取蛋白名和基因名
处理多个研究的差异表达数据
"""

import requests
import pandas as pd
import time
import os
import re

def extract_protein_id(column_name, study_id):
    """
    根据不同研究提取蛋白ID
    
    参数:
        column_name: 列名
        study_id: 研究编号 (1-11)
    
    返回:
        list: 提取的UniProt ID列表
    """
    # 跳过ID列和Group列
    if column_name in ['ID', 'SampleID', 'Group_New', 'GUID']:
        return None
    
    # Study 1, 2, 3_Olink, 4, 10, 11: 格式为 UBA52|M0R1V7 或 M0R1V7
    if study_id in ['1', '2', '3_Olink', '4', '10', '11']:
        if '|' in column_name:
            # 提取|后面的部分
            parts = column_name.split('|')[1]
            # 如果有逗号，分割成多个
            if ',' in parts:
                return [p.strip() for p in parts.split(',')]
            return [parts.strip()]
        else:
            # 直接使用列名
            return [column_name.strip()]
    
    # Study 3_Soma: 格式为 PNP|P00491^SL005262@PNP.10039.32
    elif study_id == '3_Soma':
        if '|' in column_name:
            after_pipe = column_name.split('|')[1]
            # 如果有^，提取|后^前的部分
            if '^' in after_pipe:
                uniprot_id = after_pipe.split('^')[0]
                return [uniprot_id.strip()]
            else:
                return [after_pipe.strip()]
        else:
            return [column_name.strip()]
    
    # Study 5: 格式为 C3_P01024
    elif study_id == '5':
        if '_' in column_name:
            return [column_name.split('_')[1].strip()]
        else:
            return [column_name.strip()]
    
    # Study 6: 格式为 P02768
    elif study_id == '6':
        return [column_name.strip()]
    
    # Study 7: 格式为 CHGB...96 或 NOMO1;NOMO2;NOMO3
    elif study_id == '7':
        # 删除...后面的内容
        if '...' in column_name:
            column_name = column_name.split('...')[0]
        
        # 如果有分号，拆分
        if ';' in column_name:
            return [name.strip() for name in column_name.split(';')]
        else:
            return [column_name.strip()]
    
    # Study 8, 9: 直接使用列名
    elif study_id in ['8', '9']:
        return [column_name.strip()]
    
    return None

def query_uniprot(query_id, is_human=True, is_reviewed=False):
    """
    查询UniProt数据库获取蛋白信息
    
    参数:
        query_id: UniProt ID或基因名
        is_human: 是否限制为人类（强制）
        is_reviewed: 是否限制为reviewed条目（不强制，仅记录状态）
    
    返回:
        dict: 包含蛋白信息的字典
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    # 构建查询：只强制人类物种，不强制reviewed状态
    query_parts = [query_id]
    if is_human:
        query_parts.append("(organism_id:9606)")
    
    params = {
        'query': ' AND '.join(query_parts),
        'format': 'json',
        'size': 1
    }
    
    try:
        print(f"  查询: {query_id}...", end='')
        response = requests.get(base_url, params=params, timeout=10)
        response.raise_for_status()
        
        data = response.json()
        
        if data.get('results') and len(data['results']) > 0:
            result = data['results'][0]
            
            # 提取信息
            uniprot_id = result.get('primaryAccession', 'NA')
            
            # 获取基因名
            genes = result.get('genes', [])
            gene_name = 'NA'
            all_gene_names = []
            if genes and len(genes) > 0:
                gene_name = genes[0].get('geneName', {}).get('value', 'NA')
                # 获取所有基因名（包括同义词）
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
            is_reviewed_status = 'Yes' if entry_type == 'UniProtKB reviewed (Swiss-Prot)' else 'No'
            
            print(f" 找到: {gene_name}")
            
            return {
                'Query': query_id,
                'UniProt_ID': uniprot_id,
                'Gene_Name': gene_name,
                'All_Gene_Names': ';'.join(all_gene_names) if all_gene_names else 'NA',
                'Protein_Name': protein_name,
                'UniProtKB_Reviewed': is_reviewed_status
            }
        else:
            print(f" 未找到")
            return {
                'Query': query_id,
                'UniProt_ID': 'NA',
                'Gene_Name': 'NA',
                'All_Gene_Names': 'NA',
                'Protein_Name': 'NA',
                'UniProtKB_Reviewed': 'NA'
            }
            
    except Exception as e:
        print(f" 错误: {str(e)}")
        return {
            'Query': query_id,
            'UniProt_ID': 'NA',
            'Gene_Name': 'NA',
            'All_Gene_Names': 'NA',
            'Protein_Name': 'NA',
            'UniProtKB_Reviewed': f'Error: {str(e)}'
        }

def process_study(study_file, study_id, output_dir, test_mode=True):
    """
    处理单个研究文件
    
    参数:
        study_file: 研究文件路径
        study_id: 研究编号
        output_dir: 输出目录
        test_mode: 是否测试模式（只处理前5个）
    """
    print(f"\n{'='*80}")
    print(f"处理 Study {study_id}: {study_file}")
    print(f"{'='*80}")
    
    # 读取文件
    df = pd.read_csv(study_file)
    print(f"文件包含 {len(df.columns)} 列")
    
    # 获取所有列名（排除ID列和Group列）
    protein_columns = [col for col in df.columns if col not in ['ID', 'SampleID', 'Group_New', 'GUID']]
    
    if test_mode:
        protein_columns = protein_columns[:5]
        print(f"测试模式：只处理前 {len(protein_columns)} 个蛋白")
    
    # 存储查询结果和映射关系
    query_results = []
    column_mapping = {}  # 原始列名 -> 新基因名
    
    # 处理每个蛋白列
    for col in protein_columns:
        protein_ids = extract_protein_id(col, study_id)
        
        if protein_ids is None:
            continue
        
        print(f"\n原始列名: {col}")
        print(f"  提取的ID: {protein_ids}")
        
        # 对于study 7，如果有多个ID（分号分隔），需要查询每个并合并
        if len(protein_ids) > 1:
            # 查询每个ID
            sub_results = []
            for pid in protein_ids:
                result = query_uniprot(pid, is_human=True, is_reviewed=False)
                result['Study'] = study_id
                result['Original_Column'] = col
                sub_results.append(result)
                query_results.append(result)
                time.sleep(0.3)
            
            # 使用第一个成功的基因名作为列名
            for res in sub_results:
                if res['Gene_Name'] != 'NA':
                    column_mapping[col] = res['Gene_Name']
                    break
            if col not in column_mapping:
                column_mapping[col] = col  # 如果都失败，保持原名
        else:
            # 单个ID
            result = query_uniprot(protein_ids[0], is_human=True, is_reviewed=False)
            result['Study'] = study_id
            result['Original_Column'] = col
            query_results.append(result)
            
            # 映射列名
            if result['Gene_Name'] != 'NA':
                column_mapping[col] = result['Gene_Name']
            else:
                column_mapping[col] = col
            
            time.sleep(0.3)
    
    # 保存查询结果
    results_df = pd.DataFrame(query_results)
    results_file = os.path.join(output_dir, f'study_{study_id}_uniprot_mapping.csv')
    results_df.to_csv(results_file, index=False, encoding='utf-8-sig')
    print(f"\n查询结果已保存: {results_file}")
    
    # 重命名数据文件的列
    df_renamed = df.copy()
    df_renamed.rename(columns=column_mapping, inplace=True)
    
    # 保存重命名后的文件
    renamed_file = os.path.join(output_dir, f'study_{study_id}_renamed.csv')
    df_renamed.to_csv(renamed_file, index=False, encoding='utf-8-sig')
    print(f"重命名后的数据已保存: {renamed_file}")
    
    return results_df

def main():
    # 定义研究文件
    base_dir = 'F:/1a-EOD-CSF-protein/DE-results'
    output_dir = 'F:/1a-EOD-CSF-protein/uniport'
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    studies = {
        '2': 'study_2.csv',
        '1': 'study_1.csv',
        '3_Olink': 'study_3_Olink.csv',
        '3_Soma': 'study_3_Soma.csv',
        '4': 'study_4.csv',
        '5': 'study_5.csv',
        '6': 'study_6.csv',
        '7': 'study_7.csv',
        '8': 'study_8.csv',
        '9': 'study_9.csv',
        '10': 'study_10.csv',
        '11': 'study_11.csv'
    }
    
    # 测试模式：只处理前5个蛋白
    test_mode = False
    
    all_results = []
    
    for study_id, filename in studies.items():
        study_file = os.path.join(base_dir, filename)
        
        if not os.path.exists(study_file):
            print(f"\n警告: 文件不存在 - {study_file}")
            continue
        
        try:
            results = process_study(study_file, study_id, output_dir, test_mode=test_mode)
            all_results.append(results)
        except Exception as e:
            print(f"\n错误: 处理 Study {study_id} 时出错: {str(e)}")
            import traceback
            traceback.print_exc()
            continue
    
    # 合并所有结果
    if all_results:
        combined_results = pd.concat(all_results, ignore_index=True)
        combined_file = os.path.join(output_dir, 'all_studies_uniprot_mapping.csv')
        combined_results.to_csv(combined_file, index=False, encoding='utf-8-sig')
        print(f"\n{'='*80}")
        print(f"所有研究的查询结果已合并保存: {combined_file}")
        print(f"{'='*80}")
    
    print("\n处理完成！")

if __name__ == '__main__':
    main()
