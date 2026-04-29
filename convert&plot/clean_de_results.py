#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
数据清洗脚本：根据UniProt映射关系转换DE-results数据集的列名

功能：
1. 读取ID_mapping_final.csv获取Original_Column到Gene_Name的映射关系
2. 对每个研究数据集，保留第一列（SampleID），转换其他列名为Gene_Name
3. 处理一对多映射：如果一个Original_Column映射到多个Gene_Name，复制该列
4. 删除无法映射（Gene_Name为NA）的列
5. 生成清洗报告

作者：AI Assistant
日期：2026-02-24
标准：Nature级别的科学严谨性
"""

import pandas as pd
import numpy as np
import os
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

class DataCleaner:
    """数据清洗类"""
    
    def __init__(self, mapping_file, input_dir, output_dir):
        """
        初始化数据清洗器
        
        参数:
            mapping_file: ID映射文件路径
            input_dir: 输入数据目录
            output_dir: 输出数据目录
        """
        self.mapping_file = mapping_file
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.mapping_dict = {}
        self.cleaning_report = []
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        
    def load_mapping(self):
        """
        加载映射关系
        
        返回:
            dict: {study: {original_column: [gene_names]}}
        """
        print("=" * 80)
        print("步骤 1: 加载映射关系")
        print("=" * 80)
        
        df = pd.read_csv(self.mapping_file, encoding='utf-8-sig')
        print(f"映射文件记录数: {len(df)}")
        
        # 按研究分组构建映射字典
        for study in df['Study'].unique():
            study_df = df[df['Study'] == study]
            self.mapping_dict[study] = defaultdict(list)
            
            for _, row in study_df.iterrows():
                original_col = row['Original_Column']
                gene_name = row['Gene_Name']
                
                # 只添加有效的基因名（非NA）
                if pd.notna(gene_name) and gene_name != 'NA':
                    self.mapping_dict[study][original_col].append(gene_name)
            
            print(f"  Study {study}: {len(self.mapping_dict[study])} 个唯一的Original_Column")
        
        print(f"\n总共加载了 {len(self.mapping_dict)} 个研究的映射关系")
        return self.mapping_dict
    
    def clean_single_study(self, study_id, study_file):
        """
        清洗单个研究数据集
        
        参数:
            study_id: 研究ID（如'1', '2', '3_Olink'等）
            study_file: 研究数据文件路径
        
        返回:
            tuple: (cleaned_df, report_dict)
        """
        print(f"\n{'='*80}")
        print(f"处理 Study {study_id}: {os.path.basename(study_file)}")
        print(f"{'='*80}")
        
        # 读取原始数据
        df = pd.read_csv(study_file, encoding='utf-8-sig')
        original_cols = df.columns.tolist()
        n_original_cols = len(original_cols)
        
        print(f"原始数据: {df.shape[0]} 行 × {df.shape[1]} 列")
        
        # 第一列是SampleID，保留
        id_col = original_cols[0]
        data_cols = original_cols[1:]
        
        print(f"ID列: {id_col}")
        print(f"数据列数: {len(data_cols)}")
        
        # 获取该研究的映射关系
        if study_id not in self.mapping_dict:
            print(f"警告: 未找到Study {study_id}的映射关系，跳过")
            return None, None
        
        study_mapping = self.mapping_dict[study_id]
        
        # 第一步：构建基因名到原始列的映射（多对一）
        gene_to_cols = defaultdict(list)
        for col in data_cols:
            if col in study_mapping:
                gene_names = study_mapping[col]
                for gene_name in gene_names:
                    gene_to_cols[gene_name].append(col)
        
        # 构建新的数据框
        new_df = pd.DataFrame()
        new_df[id_col] = df[id_col]
        
        # 统计信息
        mapped_cols = []
        unmapped_cols = []
        expanded_cols = []  # 一对多映射的列（一个原始列映射到多个基因）
        merged_genes = []  # 多对一映射的基因（多个原始列映射到一个基因）
        
        # 处理每一列
        for col in data_cols:
            if col in study_mapping:
                gene_names = study_mapping[col]
                
                if len(gene_names) == 1:
                    # 一对一映射
                    mapped_cols.append((col, gene_names[0]))
                else:
                    # 一对多映射：复制列
                    for gene_name in gene_names:
                        expanded_cols.append((col, gene_name))
                    mapped_cols.append((col, f"{len(gene_names)} genes"))
            else:
                # 未映射的列
                unmapped_cols.append(col)
        
        # 第二步：对每个基因名，如果有多个原始列映射到它，取平均值
        for gene_name, orig_cols in gene_to_cols.items():
            if len(orig_cols) == 1:
                # 只有一个原始列映射到这个基因
                new_df[gene_name] = df[orig_cols[0]]
            else:
                # 多个原始列映射到同一个基因，取平均值
                cols_data = [df[col] for col in orig_cols]
                new_df[gene_name] = pd.concat(cols_data, axis=1).mean(axis=1)
                merged_genes.append((gene_name, orig_cols))
        
        # 计算最终列数（去重后的基因数）
        final_gene_count = len(gene_to_cols)
        
        # 生成报告
        report = {
            'Study': study_id,
            'Original_Columns': n_original_cols - 1,  # 不计ID列
            'Mapped_Columns': len(mapped_cols),
            'Unmapped_Columns': len(unmapped_cols),
            'Final_Columns': final_gene_count,
            'Expanded_Columns': len(expanded_cols),
            'Merged_Genes': len(merged_genes),
            'Unmapped_List': '; '.join(unmapped_cols[:10]) + ('...' if len(unmapped_cols) > 10 else '')
        }
        
        print(f"\n清洗结果:")
        print(f"  原始数据列: {report['Original_Columns']}")
        print(f"  成功映射列: {report['Mapped_Columns']}")
        print(f"  未映射列（删除）: {report['Unmapped_Columns']}")
        print(f"  一对多扩展: {report['Expanded_Columns']}")
        print(f"  合并基因数（多列→一基因，取平均）: {report['Merged_Genes']}")
        print(f"  最终数据列: {report['Final_Columns']}")
        
        if unmapped_cols:
            print(f"\n  未映射的列（前10个）:")
            for col in unmapped_cols[:10]:
                print(f"    - {col}")
            if len(unmapped_cols) > 10:
                print(f"    ... 还有 {len(unmapped_cols) - 10} 个")
        
        if expanded_cols:
            print(f"\n  一对多扩展的列（前5个）:")
            expansion_summary = defaultdict(list)
            for orig_col, gene in expanded_cols:
                expansion_summary[orig_col].append(gene)
            
            for i, (orig_col, genes) in enumerate(list(expansion_summary.items())[:5]):
                print(f"    - {orig_col} → {', '.join(genes)}")
            if len(expansion_summary) > 5:
                print(f"    ... 还有 {len(expansion_summary) - 5} 个")
        
        if merged_genes:
            print(f"\n  合并的基因（多列取平均，前5个）:")
            for i, (gene_name, orig_cols) in enumerate(merged_genes[:5]):
                print(f"    - {gene_name} ← 平均({len(orig_cols)}列): {', '.join(orig_cols[:3])}{' ...' if len(orig_cols) > 3 else ''}")
            if len(merged_genes) > 5:
                print(f"    ... 还有 {len(merged_genes) - 5} 个")
        
        # 保存清洗后的数据
        output_file = os.path.join(self.output_dir, f'study_{study_id}_cleaned.csv')
        new_df.to_csv(output_file, index=False, encoding='utf-8-sig')
        print(f"\n[OK] 已保存: {output_file}")
        print(f"  最终数据: {new_df.shape[0]} 行 × {new_df.shape[1]} 列")
        
        return new_df, report
    
    def process_all_studies(self):
        """
        处理所有研究数据集
        """
        print("\n" + "=" * 80)
        print("步骤 2: 处理所有研究数据集")
        print("=" * 80)
        
        # 定义要处理的研究及其文件
        studies = [
            ('1', 'study_1.csv'),
            ('2', 'study_2.csv'),
            ('3_Olink', 'study_3_Olink.csv'),
            ('3_Soma', 'study_3_Soma.csv'),
            ('4', 'study_4.csv'),
            ('5', 'study_5.csv'),
            ('6', 'study_6.csv'),
            ('7', 'study_7.csv'),
            ('8', 'study_8.csv'),
            ('9', 'study_9.csv'),
            ('10', 'study_10.csv'),
            ('11', 'study_11.csv'),
        ]
        
        for study_id, filename in studies:
            study_file = os.path.join(self.input_dir, filename)
            
            if not os.path.exists(study_file):
                print(f"\n警告: 文件不存在 - {study_file}")
                continue
            
            try:
                cleaned_df, report = self.clean_single_study(study_id, study_file)
                if report:
                    self.cleaning_report.append(report)
            except Exception as e:
                print(f"\n错误: 处理Study {study_id}时出错: {str(e)}")
                import traceback
                traceback.print_exc()
                continue
    
    def generate_summary_report(self):
        """
        生成汇总报告
        """
        print("\n" + "=" * 80)
        print("步骤 3: 生成汇总报告")
        print("=" * 80)
        
        if not self.cleaning_report:
            print("没有生成任何报告")
            return
        
        # 创建报告DataFrame
        report_df = pd.DataFrame(self.cleaning_report)
        
        # 保存报告
        report_file = os.path.join(self.output_dir, 'cleaning_summary_report.csv')
        report_df.to_csv(report_file, index=False, encoding='utf-8-sig')
        
        print(f"\n汇总报告:")
        print(report_df.to_string(index=False))
        
        print(f"\n[OK] 报告已保存: {report_file}")
        
        # 打印总体统计
        print("\n" + "=" * 80)
        print("总体统计")
        print("=" * 80)
        print(f"处理的研究数: {len(report_df)}")
        print(f"原始总列数: {report_df['Original_Columns'].sum()}")
        print(f"成功映射列数: {report_df['Mapped_Columns'].sum()}")
        print(f"删除列数: {report_df['Unmapped_Columns'].sum()}")
        print(f"最终总列数: {report_df['Final_Columns'].sum()}")
        print(f"扩展列数: {report_df['Expanded_Columns'].sum()}")
        print(f"合并基因数: {report_df['Merged_Genes'].sum()}")
        print(f"映射成功率: {report_df['Mapped_Columns'].sum() / report_df['Original_Columns'].sum() * 100:.2f}%")

def main():
    """主函数"""
    print("\n" + "=" * 80)
    print("数据清洗脚本 - Nature级别标准")
    print("=" * 80)
    
    # 配置路径
    mapping_file = 'F:/1a-EOD-CSF-protein/uniport/ID_mapping_final.csv'
    input_dir = 'F:/1a-EOD-CSF-protein/DE-results'
    output_dir = 'F:/1a-EOD-CSF-protein/transform-final'
    
    # 验证输入文件
    if not os.path.exists(mapping_file):
        print(f"错误: 映射文件不存在 - {mapping_file}")
        return
    
    if not os.path.exists(input_dir):
        print(f"错误: 输入目录不存在 - {input_dir}")
        return
    
    # 创建清洗器并执行
    cleaner = DataCleaner(mapping_file, input_dir, output_dir)
    
    # 步骤1: 加载映射关系
    cleaner.load_mapping()
    
    # 步骤2: 处理所有研究
    cleaner.process_all_studies()
    
    # 步骤3: 生成汇总报告
    cleaner.generate_summary_report()
    
    print("\n" + "=" * 80)
    print("[OK] 所有处理完成！")
    print("=" * 80)

if __name__ == '__main__':
    main()
