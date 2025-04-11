#!/usr/bin/env python3
"""
测试拷贝数变异(重复)检测
"""
import os
import sys
import argparse
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sv_aligner.data_types import AlignmentSegment
from sv_aligner.sv_detector import detect_duplication, format_duplication_report

def create_duplication_example():
    """创建一个包含拷贝数变异的比对结果"""
    # 示例：一个参考序列上的区域在查询序列中重复出现
    ref_name = "reference"
    ref_len = 1000
    query_name = "query"
    query_len = 1200
    
    # 第一个片段：第一个副本
    seg1 = AlignmentSegment(
        q_name=query_name, 
        q_len=query_len, 
        q_st=100, 
        q_en=200,  # 100bp长的片段
        r_name=ref_name, 
        r_len=ref_len, 
        r_st=500, 
        r_en=600,  # 对应参考序列的500-600
        strand="+", 
        score=100, 
        edit_distance=0, 
        cigar="100M"
    )
    
    # 第二个片段：第二个副本，几乎相同的序列但在查询序列中的不同位置
    seg2 = AlignmentSegment(
        q_name=query_name, 
        q_len=query_len, 
        q_st=300, 
        q_en=400,  # 另一个100bp长的片段
        r_name=ref_name, 
        r_len=ref_len, 
        r_st=500, 
        r_en=600,  # 映射到相同的参考区域
        strand="+", 
        score=98,  # 略低的分数表示有少量变异
        edit_distance=2,  # 2个碱基的差异
        cigar="100M"
    )
    
    # 第三个片段：随机匹配，不是重复的一部分
    seg3 = AlignmentSegment(
        q_name=query_name, 
        q_len=query_len, 
        q_st=500, 
        q_en=700,
        r_name=ref_name, 
        r_len=ref_len, 
        r_st=700, 
        r_en=900,
        strand="+", 
        score=200, 
        edit_distance=0, 
        cigar="200M"
    )
    
    # 第四个片段：反向互补的重复
    seg4 = AlignmentSegment(
        q_name=query_name, 
        q_len=query_len, 
        q_st=800, 
        q_en=900,
        r_name=ref_name, 
        r_len=ref_len, 
        r_st=500, 
        r_en=600,
        strand="-",  # 注意这里是负链
        score=95, 
        edit_distance=5, 
        cigar="100M"
    )
    
    return [seg1, seg2, seg3, seg4]

def main():
    parser = argparse.ArgumentParser(description="测试拷贝数变异检测")
    parser.add_argument("--output", "-o", default=None, help="输出文件路径")
    parser.add_argument("--min-length", type=int, default=50, help="最小重复长度")
    parser.add_argument("--min-similarity", type=float, default=0.85, help="最小相似度阈值")
    parser.add_argument("--use-english", action="store_true", help="使用英文输出")
    
    args = parser.parse_args()
    
    # 创建测试数据
    alignments = create_duplication_example()
    
    # 使用中文或英文
    use_chinese = not args.use_english
    
    if use_chinese:
        print("检测拷贝数变异(重复)...")
    else:
        print("Detecting copy number variations (duplications)...")
    
    # 检测重复
    duplications = detect_duplication(
        alignments, 
        min_length=args.min_length,
        min_similarity=args.min_similarity
    )
    
    # 输出结果
    output_file = args.output
    output_handle = open(output_file, "w", encoding="utf-8") if output_file else sys.stdout
    
    try:
        if use_chinese:
            output_handle.write("## 拷贝数变异检测结果 ##\n\n")
        else:
            output_handle.write("## COPY NUMBER VARIATION DETECTION RESULTS ##\n\n")
        
        if duplications:
            for i, dup in enumerate(duplications, 1):
                output_handle.write(f"[{i}] {format_duplication_report(dup, use_chinese)}\n\n")
        else:
            if use_chinese:
                output_handle.write("未检测到拷贝数变异。\n")
            else:
                output_handle.write("No copy number variations detected.\n")
    finally:
        if output_file and output_handle != sys.stdout:
            output_handle.close()
            if use_chinese:
                print(f"结果已保存到 {output_file}")
            else:
                print(f"Results saved to {output_file}")
    
    if use_chinese:
        print(f"检测到 {len(duplications)} 个拷贝数变异")
    else:
        print(f"Detected {len(duplications)} copy number variations")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
