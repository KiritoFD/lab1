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

def create_example_duplication():
    """创建一个包含重复的示例比对结果"""
    # 示例：查询序列中的两个不同区域映射到参考序列的相同区域
    alignments = [
        # 第一个重复片段
        AlignmentSegment(
            q_name="query1", q_len=600, q_st=0, q_en=100,
            r_name="ref1", r_len=400, r_st=200, r_en=300,
            strand="+", score=100, edit_distance=0, cigar="100M"
        ),
        # 第二个重复片段（映射到相同的参考区域）
        AlignmentSegment(
            q_name="query1", q_len=600, q_st=200, q_en=300,
            r_name="ref1", r_len=400, r_st=200, r_en=300,
            strand="+", score=95, edit_distance=2, cigar="100M"
        ),
        # 不相关的比对片段
        AlignmentSegment(
            q_name="query1", q_len=600, q_st=350, q_en=500,
            r_name="ref1", r_len=400, r_st=0, r_en=150,
            strand="+", score=150, edit_distance=0, cigar="150M"
        ),
        # 另一个重复，但反向互补
        AlignmentSegment(
            q_name="query2", q_len=500, q_st=0, q_en=100,
            r_name="ref1", r_len=400, r_st=100, r_en=200,
            strand="+", score=100, edit_distance=0, cigar="100M"
        ),
        AlignmentSegment(
            q_name="query2", q_len=500, q_st=150, q_en=250,
            r_name="ref1", r_len=400, r_st=100, r_en=200,
            strand="-", score=95, edit_distance=3, cigar="100M"
        )
    ]
    return alignments

def main():
    parser = argparse.ArgumentParser(description="拷贝数变异(重复)检测示例")
    parser.add_argument("--output", "-o", default=None,
                        help="输出文件路径")
    parser.add_argument("--min-length", type=int, default=50,
                        help="最小重复长度 (默认: 50bp)")
    parser.add_argument("--min-similarity", type=float, default=0.85,
                        help="最小相似度阈值 (默认: 0.85)")
    parser.add_argument("--use-english", action="store_true", 
                        help="使用英语输出 (默认: 中文)")
    args = parser.parse_args()
    
    # 设置语言
    use_chinese = not args.use_english
    
    if use_chinese:
        print("创建示例比对结果...")
    else:
        print("Creating example alignments...")
    
    alignments = create_example_duplication()
    
    if use_chinese:
        print("检测拷贝数变异(重复)...")
    else:
        print("Detecting duplications...")
    
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
