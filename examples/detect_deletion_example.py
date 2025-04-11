#!/usr/bin/env python3
"""
大型缺失检测示例脚本
"""
import os
import sys
import argparse
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sv_aligner.data_types import AlignmentSegment
from sv_aligner.sv_detector import detect_and_report_deletions

def create_example_deletion():
    """创建一个包含缺失的示例比对结果"""
    # 示例查询序列中的缺失：缺少参考序列中100-200位置的片段
    alignments = [
        # 缺失左侧片段
        AlignmentSegment(
            q_name="query1", q_len=300, q_st=0, q_en=100,
            r_name="ref1", r_len=400, r_st=0, r_en=100,
            strand="+", score=100, edit_distance=0, cigar="100M"
        ),
        # 缺失右侧片段
        AlignmentSegment(
            q_name="query1", q_len=300, q_st=100, q_en=300,
            r_name="ref1", r_len=400, r_st=200, r_en=400,
            strand="+", score=200, edit_distance=0, cigar="200M"
        ),
        # 另一个查询序列，没有缺失
        AlignmentSegment(
            q_name="query2", q_len=400, q_st=0, q_en=400,
            r_name="ref1", r_len=400, r_st=0, r_en=400,
            strand="+", score=400, edit_distance=0, cigar="400M"
        )
    ]
    return alignments

def main():
    parser = argparse.ArgumentParser(description="Deletion Detection Example")
    parser.add_argument("--output", "-o", default=None,
                        help="Output file path")
    parser.add_argument("--min-size", type=int, default=50,
                        help="Minimum deletion size (default: 50bp)")
    parser.add_argument("--use-chinese", action="store_true", 
                        help="Use Chinese for output messages (default: English)")
    args = parser.parse_args()
    
    # Set default language
    use_chinese = args.use_chinese
    
    if use_chinese:
        print("创建示例比对结果...")
    else:
        print("Creating example alignments...")
    
    alignments = create_example_deletion()
    
    if use_chinese:
        print("检测大型缺失...")
    else:
        print("Detecting large deletions...")
    
    deletions = detect_and_report_deletions(
        alignments, 
        min_size=args.min_size,
        output_file=args.output,
        use_chinese=use_chinese
    )
    
    if use_chinese:
        print(f"检测完成。发现 {len(deletions)} 个缺失。")
    else:
        print(f"Detection completed. Found {len(deletions)} deletions.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
