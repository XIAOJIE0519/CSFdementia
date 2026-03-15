#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
去除study 7, 8, 9中的重复Query记录
当同一个study中有相同的Query时，只保留一条记录
"""

import pandas as pd
import os

def deduplicate_studies_789(input_file, output_file):
    """
    去除study 7, 8, 9中的重复Query
    
    参数:
        input_file: 输入文件路径
        output_file: 输出文件路径
    """
    print(f"读取文件: {input_file}")
    df = pd.read_csv(input_file, encoding='utf-8-sig')
    print(f"原始记录数: {len(df)}")
    
    # 分离study 7, 8, 9和其他研究
    df_789 = df[df['Study'].isin(['7', '8', '9'])].copy()
    df_other = df[~df['Study'].isin(['7', '8', '9'])].copy()
    
    print(f"\nStudy 7, 8, 9 记录数: {len(df_789)}")
    print(f"其他研究记录数: {len(df_other)}")
    
    # 统计去重前的情况
    print("\n去重前统计:")
    for study in ['7', '8', '9']:
        study_df = df_789[df_789['Study'] == study]
        unique_queries = study_df['Query'].nunique()
        total_rows = len(study_df)
        duplicates = total_rows - unique_queries
        print(f"  Study {study}: {total_rows} 条记录, {unique_queries} 个唯一Query, {duplicates} 条重复")
    
    # 对study 7, 8, 9按照Study和Query去重，保留第一条
    df_789_dedup = df_789.drop_duplicates(subset=['Study', 'Query'], keep='first')
    
    print("\n去重后统计:")
    for study in ['7', '8', '9']:
        study_df = df_789_dedup[df_789_dedup['Study'] == study]
        unique_queries = study_df['Query'].nunique()
        total_rows = len(study_df)
        print(f"  Study {study}: {total_rows} 条记录, {unique_queries} 个唯一Query")
    
    # 合并所有数据
    df_final = pd.concat([df_other, df_789_dedup], ignore_index=True)
    
    # 按照原始顺序排序（先按Study，再按原始索引）
    # 为了保持原始顺序，我们按Study排序
    study_order = ['2', '1', '3_Olink', '3_Soma', '4', '5', '6', '7', '8', '9', '10', '11']
    df_final['Study'] = pd.Categorical(df_final['Study'], categories=study_order, ordered=True)
    df_final = df_final.sort_values('Study').reset_index(drop=True)
    
    print(f"\n最终记录数: {len(df_final)}")
    print(f"减少了 {len(df) - len(df_final)} 条重复记录")
    
    # 保存结果
    df_final.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"\n结果已保存到: {output_file}")
    
    return df_final

def main():
    input_file = 'F:/1a-EOD-CSF-protein/uniport/all_studies_uniprot_mapping - final.csv'
    output_file = 'F:/1a-EOD-CSF-protein/uniport/ID_mapping_final.csv'
    
    if not os.path.exists(input_file):
        print(f"错误: 输入文件不存在 - {input_file}")
        return
    
    deduplicate_studies_789(input_file, output_file)
    
    print("\n处理完成！")

if __name__ == '__main__':
    main()
