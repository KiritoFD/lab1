#!/usr/bin/env python3
"""
拷贝数变异(重复)检测示例脚本
"""
import os
import sys
import argparse
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sv_aligner.data_types import AlignmentSegment
from sv_aligner.duplication_detector import detect_and_report_duplications

def create_example_alignments():
    """创建测试用的比对结果"""
    # 创建一个包含重复的例子：参考序列上的同一个区域在查询序列上出现两次
    ref_name = "reference"
    ref_len = 1000
    query_name = "query"
    query_len = 1200
    
    alignments = [
        # 第一个片段：查询序列中的第一个重复
        AlignmentSegment(
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
        ),
        # 第二个片段：查询序列中的第二个重复，映射到同一参考位置
        AlignmentSegment(
            q_name=query_name,
            q_len=query_len,
            q_st=300,
            q_en=400,  # 另一个100bp的片段
            r_name=ref_name,
            r_len=ref_len,
            r_st=500,
            r_en=600,  # 映射到相同的参考区域
            strand="+",
            score=95,
            edit_distance=5,
            cigar="100M"
        ),
        # 第三个片段：无关的片段
        AlignmentSegment(
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
        ),
        # 第四个片段：第三个重复，但是反向互补的
        AlignmentSegment(
            q_name=query_name,
            q_len=query_len,
            q_st=800,
            q_en=900,
            r_name=ref_name,
            r_len=ref_len,
            r_st=500,
            r_en=600,
            strand="-",  # 反向链
            score=90,
            edit_distance=10,
            cigar="100M"
        )
    ]
    
    return alignments

def main():
    parser = argparse.ArgumentParser(description="拷贝数变异(重复)检测示例")
    parser.add_argument("--output", "-o", default=None, help="输出文件路径")
    parser.add_argument("--min-length", type=int, default=50, help="最小重复长度")
    parser.add_argument("--min-similarity", type=float, default=0.85, help="最小相似度阈值")
    parser.add_argument("--use-english", action="store_true", help="使用英文输出")
    
    args = parser.parse_args()
    
    # 使用中文或英文
    use_chinese = not args.use_english
    
    if use_chinese:
        print("创建示例比对结果...")
    else:
        print("Creating example alignments...")
    
    # 创建示例比对结果
    alignments = create_example_alignments()
    
    if use_chinese:
        print(f"创建了 {len(alignments)} 个比对片段")
        print("检测拷贝数变异(重复)...")
    else:
        print(f"Created {len(alignments)} alignment segments")
        print("Detecting duplications...")
    
    # 检测重复
    duplications = detect_and_report_duplications(
        alignments,
        min_length=args.min_length,
        min_similarity=args.min_similarity,
        output_file=args.output,
        use_chinese=use_chinese
    )
    
    if use_chinese:
        print(f"检测完成。发现 {len(duplications)} 个重复。")
    else:
        print(f"Detection completed. Found {len(duplications)} duplications.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
