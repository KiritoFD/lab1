#!/usr/bin/env python3
"""
比对结果可视化工具：将TSV格式的比对结果以图形方式展示
需要安装matplotlib: pip install matplotlib
"""

import os
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

def parse_tsv(tsv_file):
    """解析TSV格式的比对结果"""
    alignments = []
    with open(tsv_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue  # 跳过注释行
            
            fields = line.strip().split('\t')
            if len(fields) < 12:
                print(f"警告: 行格式不正确: {line}", file=sys.stderr)
                continue
                
            alignment = {
                'q_name': fields[0],
                'q_len': int(fields[1]),
                'q_st': int(fields[2]),
                'q_en': int(fields[3]),
                'r_name': fields[4],
                'r_len': int(fields[5]),
                'r_st': int(fields[6]),
                'r_en': int(fields[7]),
                'strand': fields[8],
                'score': int(fields[9]),
                'edit_distance': fields[10],
                'cigar': fields[11]
            }
            alignments.append(alignment)
    
    return alignments

def visualize_alignments(alignments, output_file=None):
    """可视化比对结果"""
    if not alignments:
        print("错误: 没有找到有效的比对结果", file=sys.stderr)
        return False
        
    # 按查询名称分组
    queries = {}
    for aln in alignments:
        q_name = aln['q_name']
        if q_name not in queries:
            queries[q_name] = []
        queries[q_name].append(aln)
    
    # 为每个查询创建一个子图
    n_queries = len(queries)
    fig, axes = plt.subplots(n_queries, 1, figsize=(12, 3 * n_queries), squeeze=False)
    
    # 创建颜色映射 - 根据比对得分着色
    colors = ['#ffcccc', '#ff9999', '#ff6666', '#ff3333', '#ff0000']
    cmap = LinearSegmentedColormap.from_list("score_cmap", colors)
    
    # 最高得分（用于归一化颜色）
    max_score = max(aln['score'] for aln in alignments)
    
    for i, (q_name, alns) in enumerate(queries.items()):
        ax = axes[i, 0]
        
        # 获取查询序列长度
        q_len = alns[0]['q_len']
        
        # 画查询序列条带
        ax.add_patch(patches.Rectangle((0, 0.45), q_len, 0.1, 
                                       facecolor='lightgray', edgecolor='black'))
        
        # 画每个比对段
        for aln in alns:
            # 归一化得分以获取颜色
            norm_score = aln['score'] / max_score if max_score > 0 else 0
            color = cmap(norm_score)
            
            # 在查询序列上绘制比对段
            ax.add_patch(patches.Rectangle((aln['q_st'], 0.45), 
                                          aln['q_en'] - aln['q_st'], 0.1, 
                                          facecolor=color, edgecolor='black'))
            
            # 添加标签
            ax.text(aln['q_st'] + (aln['q_en'] - aln['q_st'])/2, 0.6, 
                   f"{aln['r_name']}:{aln['r_st']}-{aln['r_en']} ({aln['strand']})",
                   ha='center', rotation=45, fontsize=8)
            
            # 如果是负链比对，添加特殊标记
            if aln['strand'] == '-':
                ax.add_patch(patches.ConnectionPatch(
                    xyA=(aln['q_st'], 0.45), xyB=(aln['q_en'], 0.55),
                    coordsA="data", coordsB="data", color='blue', linewidth=0.5
                ))
                ax.add_patch(patches.ConnectionPatch(
                    xyA=(aln['q_en'], 0.45), xyB=(aln['q_st'], 0.55),
                    coordsA="data", coordsB="data", color='blue', linewidth=0.5
                ))
        
        # 设置坐标轴
        ax.set_xlim(0, q_len)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xlabel('Query Position')
        ax.set_title(f'Query: {q_name} (Length: {q_len}bp)')
        
        # 添加刻度线
        ticks = list(range(0, q_len+1, max(1, q_len//10)))
        ax.set_xticks(ticks)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"可视化结果已保存至: {output_file}")
    else:
        plt.show()
    
    plt.close()
    return True

def main():
    parser = argparse.ArgumentParser(description='可视化SV-Aware序列比对结果')
    parser.add_argument('tsv_file', help='包含比对结果的TSV文件')
    parser.add_argument('-o', '--output', default=None, help='输出图像文件名 (例如 output.png)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.tsv_file):
        print(f"错误: 文件 '{args.tsv_file}' 不存在", file=sys.stderr)
        return 1
    
    alignments = parse_tsv(args.tsv_file)
    if not alignments:
        print("错误: 未找到有效的比对结果", file=sys.stderr)
        return 1
    
    success = visualize_alignments(alignments, args.output)
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
