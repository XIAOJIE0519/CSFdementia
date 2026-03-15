#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
批量查询UniProt数据库获取蛋白名和基因名
"""

import requests
import pandas as pd
import time
from urllib.parse import quote

def query_uniprot(gene_name):
    """
    查询UniProt数据库获取蛋白信息
    
    参数:
        gene_name: 基因名或蛋白名
    
    返回:
        dict: 包含蛋白信息的字典
    """
    # UniProt REST API endpoint
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    # 构建查询参数 - 限制为人类蛋白，reviewed条目
    params = {
        'query': f'{gene_name} AND (organism_id:9606) AND (reviewed:true)',
        'format': 'json',
        'size': 1  # 只获取第一个结果
    }
    
    try:
        print(f"正在查询: {gene_name}...")
        response = requests.get(base_url, params=params, timeout=10)
        response.raise_for_status()
        
        data = response.json()
        
        if data.get('results') and len(data['results']) > 0:
            result = data['results'][0]
            
            # 提取信息
            uniprot_id = result.get('primaryAccession', 'N/A')
            protein_name = result.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'N/A')
            
            # 如果没有recommendedName，尝试获取submittedName
            if protein_name == 'N/A':
                submitted_names = result.get('proteinDescription', {}).get('submissionNames', [])
                if submitted_names:
                    protein_name = submitted_names[0].get('fullName', {}).get('value', 'N/A')
            
            # 获取基因名
            genes = result.get('genes', [])
            gene_symbol = 'N/A'
            if genes and len(genes) > 0:
                gene_symbol = genes[0].get('geneName', {}).get('value', 'N/A')
            
            return {
                'Query': gene_name,
                'UniProt_ID': uniprot_id,
                'Gene_Symbol': gene_symbol,
                'Protein_Name': protein_name,
                'Status': 'Found'
            }
        else:
            print(f"  未找到结果: {gene_name}")
            return {
                'Query': gene_name,
                'UniProt_ID': 'N/A',
                'Gene_Symbol': 'N/A',
                'Protein_Name': 'N/A',
                'Status': 'Not Found'
            }
            
    except requests.exceptions.RequestException as e:
        print(f"  查询出错 {gene_name}: {str(e)}")
        return {
            'Query': gene_name,
            'UniProt_ID': 'N/A',
            'Gene_Symbol': 'N/A',
            'Protein_Name': 'N/A',
            'Status': f'Error: {str(e)}'
        }

def main():
    # 要查询的基因/蛋白列表
    query_list = [
        'NRP2',
        'CADM3',
        'UNC5C',
        'VWC2',
        'Siglec.9',
        'CLM.6',
        'EZR',
        'SMOC2',
        'NBL1',
        'EF4',
        'SCARB2'
    ]
    
    # 清理查询名称（去除特殊字符）
    cleaned_queries = []
    for q in query_list:
        # 处理特殊格式
        cleaned = q.replace('.', '-')  # Siglec.9 -> Siglec-9
        if cleaned.startswith('CLM-'):
            cleaned = 'CXADR'  # CLM-6 通常指CXADR
        cleaned_queries.append(cleaned)
    
    results = []
    
    # 批量查询
    for i, (original, cleaned) in enumerate(zip(query_list, cleaned_queries)):
        result = query_uniprot(cleaned)
        result['Original_Query'] = original
        results.append(result)
        
        # 避免请求过快
        if i < len(query_list) - 1:
            time.sleep(0.5)
    
    # 转换为DataFrame
    df = pd.DataFrame(results)
    
    # 重新排列列顺序
    df = df[['Original_Query', 'Query', 'UniProt_ID', 'Gene_Symbol', 'Protein_Name', 'Status']]
    
    # 保存结果
    output_file = 'uniprot_query_results.csv'
    df.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"\n查询完成！结果已保存到: {output_file}")
    
    # 打印结果
    print("\n" + "="*80)
    print("查询结果:")
    print("="*80)
    for _, row in df.iterrows():
        print(f"\n原始查询: {row['Original_Query']}")
        print(f"  清理后: {row['Query']}")
        print(f"  UniProt ID: {row['UniProt_ID']}")
        print(f"  基因名: {row['Gene_Symbol']}")
        print(f"  蛋白名: {row['Protein_Name']}")
        print(f"  状态: {row['Status']}")

if __name__ == '__main__':
    main()
