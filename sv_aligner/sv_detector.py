"""
结构变异检测模块
提供用于分析比对结果并检测各类结构变异的功能
"""
from typing import List, Tuple, Dict, Optional
from sv_aligner.data_types import AlignmentSegment
from collections import defaultdict
import time
import sys

# 导入拷贝数变异检测器
from sv_aligner.duplication_detector import (
    detect_duplication, 
    format_duplication_report, 
    DuplicationPattern
)

def detect_deletion(alignments: List[AlignmentSegment], min_size: int = 50, max_query_gap: int = 5) -> List[Dict]:
    """
    检测查询序列中的大型缺失
    
    大型缺失定义为: 查询序列中连续或几乎连续的区域，在参考序列中映射位置有较大间隔
    
    Args:
        alignments: 比对结果列表，需按查询序列位置排序
        min_size: 最小缺失大小，小于此阈值的缺失不会被报告
        max_query_gap: 查询序列中允许的最大间隙，超过此值不被视为连续
    
    Returns:
        检测到的缺失列表，每个缺失为一个字典，包含以下字段:
        - query_name: 查询序列名称
        - query_left_pos: 缺失左侧的查询位置
        - query_right_pos: 缺失右侧的查询位置
        - ref_name: 参考序列名称
        - ref_left_pos: 缺失左侧的参考位置
        - ref_right_pos: 缺失右侧的参考位置
        - deletion_size: 缺失大小（参考序列中）
        - left_segment: 缺失左侧的比对片段
        - right_segment: 缺失右侧的比对片段
    """
    import sys
    
    if not alignments or len(alignments) < 2:
        print("Debug: Not enough alignment segments to detect deletions (need at least 2)", file=sys.stderr)
        return []
    
    # 按查询序列位置排序
    sorted_alignments = sorted(alignments, key=lambda a: (a.q_name, a.q_st))
    print(f"Debug: Analyzing {len(sorted_alignments)} alignment segments", file=sys.stderr)
    deletions = []
    
    for i in range(len(sorted_alignments)-1):
        curr = sorted_alignments[i]
        next_seg = sorted_alignments[i+1]
        
        # 调试输出
        print(f"Debug: Comparing segments {i} and {i+1}:", file=sys.stderr)
        print(f"  Segment {i}: q:{curr.q_st}-{curr.q_en}, r:{curr.r_st}-{curr.r_en}, {curr.strand}", file=sys.stderr)
        print(f"  Segment {i+1}: q:{next_seg.q_st}-{next_seg.q_en}, r:{next_seg.r_st}-{next_seg.r_en}, {next_seg.strand}", file=sys.stderr)
        
        # 检查是否为同一查询序列
        if curr.q_name != next_seg.q_name:
            print(f"  Skip: Different query sequences", file=sys.stderr)
            continue
        
        # 检查是否为同一参考序列
        if curr.r_name != next_seg.r_name:
            print(f"  Skip: Different reference sequences", file=sys.stderr)
            continue
            
        # 检查是否为同一链方向
        if curr.strand != next_seg.strand:
            print(f"  Skip: Different strands", file=sys.stderr)
            continue
        
        # 计算查询和参考序列中的间隔
        query_gap = next_seg.q_st - curr.q_en
        ref_gap = next_seg.r_st - curr.r_en if curr.strand == '+' else curr.r_st - next_seg.r_en
        
        print(f"  Query gap: {query_gap}", file=sys.stderr)
        print(f"  Reference gap: {ref_gap}", file=sys.stderr)
        
        # 对于负链，需要考虑坐标方向
        if curr.strand == '-':
            # 确保参考坐标按正确顺序
            if ref_gap < 0:
                ref_gap = -ref_gap
        
        # 检查查询序列间隔是否小于阈值，且参考序列间隔大于最小缺失大小
        if 0 <= query_gap <= max_query_gap and ref_gap >= min_size:
            print(f"  Deletion detected! Size: {ref_gap}bp", file=sys.stderr)
            deletion = {
                "query_name": curr.q_name,
                "query_left_pos": curr.q_en,
                "query_right_pos": next_seg.q_st,
                "ref_name": curr.r_name,
                "ref_left_pos": curr.r_en if curr.strand == '+' else curr.r_st,
                "ref_right_pos": next_seg.r_st if curr.strand == '+' else next_seg.r_en,
                "deletion_size": ref_gap,
                "left_segment": curr,
                "right_segment": next_seg
            }
            deletions.append(deletion)
        else:
            if query_gap > max_query_gap:
                print(f"  Skip: Query gap too large ({query_gap} > {max_query_gap})", file=sys.stderr)
            if ref_gap < min_size:
                print(f"  Skip: Reference gap too small ({ref_gap} < {min_size})", file=sys.stderr)
    
    return deletions

def format_deletion_report(deletion: Dict, use_chinese: bool = False) -> str:
    """
    Format report for a single deletion
    
    Args:
        deletion: Deletion information dictionary
        use_chinese: Whether to use Chinese for output messages
        
    Returns:
        Formatted deletion report text
    """
    strand = deletion["left_segment"].strand
    
    if use_chinese:
        report = [
            f"缺失检测 (DELETION):",
            f"  查询序列: {deletion['query_name']}",
            f"  查询位置: {deletion['query_left_pos']}-{deletion['query_right_pos']} (间隔: {deletion['query_right_pos']-deletion['query_left_pos']}bp)",
            f"  参考序列: {deletion['ref_name']}",
            f"  参考位置: {deletion['ref_left_pos']}-{deletion['ref_right_pos']} (链方向: {strand})",
            f"  缺失大小: {deletion['deletion_size']}bp",
            f"  左侧片段: q:{deletion['left_segment'].q_st}-{deletion['left_segment'].q_en} → r:{deletion['left_segment'].r_st}-{deletion['left_segment'].r_en}",
            f"  右侧片段: q:{deletion['right_segment'].q_st}-{deletion['right_segment'].q_en} → r:{deletion['right_segment'].r_st}-{deletion['right_segment'].r_en}",
        ]
    else:
        report = [
            f"DELETION DETECTED:",
            f"  Query sequence: {deletion['query_name']}",
            f"  Query position: {deletion['query_left_pos']}-{deletion['query_right_pos']} (gap: {deletion['query_right_pos']-deletion['query_left_pos']}bp)",
            f"  Reference sequence: {deletion['ref_name']}",
            f"  Reference position: {deletion['ref_left_pos']}-{deletion['ref_right_pos']} (strand: {strand})",
            f"  Deletion size: {deletion['deletion_size']}bp",
            f"  Left segment: q:{deletion['left_segment'].q_st}-{deletion['left_segment'].q_en} → r:{deletion['left_segment'].r_st}-{deletion['left_segment'].r_en}",
            f"  Right segment: q:{deletion['right_segment'].q_st}-{deletion['right_segment'].q_en} → r:{deletion['right_segment'].r_st}-{deletion['right_segment'].r_en}",
        ]
    
    return "\n".join(report)

def detect_and_report_deletions(alignments: List[AlignmentSegment], min_size: int = 50, 
                                max_query_gap: int = 5, output_file: Optional[str] = None,
                                use_chinese: bool = False) -> List[Dict]:
    """
    Detect and report large deletions
    
    Args:
        alignments: List of alignment results
        min_size: Minimum deletion size
        max_query_gap: Maximum gap in query sequence to consider segments adjacent
        output_file: Path to output file, if None output goes to standard output
        use_chinese: Whether to use Chinese for output messages
        
    Returns:
        List of detected deletions
    """
    import sys
    
    # Detect deletions
    deletions = detect_deletion(alignments, min_size, max_query_gap)
    
    # Prepare output
    output_handle = open(output_file, 'w', encoding='utf-8') if output_file else sys.stdout
    try:
        # Write summary
        if use_chinese:
            output_handle.write(f"## 大型缺失检测结果 ##\n")
            output_handle.write(f"总检测到 {len(deletions)} 个大型缺失 (最小大小: {min_size}bp)\n\n")
        else:
            output_handle.write(f"## LARGE DELETION DETECTION RESULTS ##\n")
            output_handle.write(f"Total detected: {len(deletions)} large deletions (minimum size: {min_size}bp)\n\n")
        
        # Write detailed information for each deletion
        for i, deletion in enumerate(deletions, 1):
            output_handle.write(f"[{i}] {format_deletion_report(deletion, use_chinese)}\n\n")
    
    finally:
        if output_file and output_handle != sys.stdout:
            output_handle.close()
            if use_chinese:
                print(f"缺失检测结果已写入 {output_file}")
            else:
                print(f"Deletion detection results written to {output_file}")
    
    return deletions

def detect_structural_variations(alignments: List[AlignmentSegment], 
                                params: Dict,
                                output_file: Optional[str] = None,
                                use_chinese: bool = True) -> Dict:
    """
    检测并报告所有类型的结构变异
    
    Args:
        alignments: 比对结果列表
        params: 参数字典，包含各种检测阈值
        output_file: 输出文件路径，如果为None则输出到标准输出
        use_chinese: 是否使用中文输出
        
    Returns:
        包含各类结构变异的字典
    """
    import sys
    import os
    
    results = {}
    
    # 创建输出目录（如果需要）
    if output_file:
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
    # 准备输出
    output_handle = open(output_file, 'w', encoding='utf-8') if output_file else sys.stdout
    try:
        # 输出标题
        if use_chinese:
            output_handle.write("# 结构变异检测结果 #\n\n")
        else:
            output_handle.write("# STRUCTURAL VARIATION DETECTION RESULTS #\n\n")
        
        # 1. 检测缺失
        if params.get('detect_deletion', True):
            if use_chinese:
                output_handle.write("## 1. 大型缺失检测\n\n")
            else:
                output_handle.write("## 1. Large Deletion Detection\n\n")
            
            deletions = detect_deletion(alignments, 
                                       params.get('min_deletion_size', 50), 
                                       params.get('max_query_gap', 5))
            
            results['deletions'] = deletions
            
            if deletions:
                for i, deletion in enumerate(deletions, 1):
                    output_handle.write(f"[1.{i}] {format_deletion_report(deletion, use_chinese)}\n\n")
            else:
                if use_chinese:
                    output_handle.write("未检测到大型缺失。\n\n")
                else:
                    output_handle.write("No large deletions detected.\n\n")
        
        # 2. 检测重复
        if params.get('detect_duplication', True):
            if use_chinese:
                output_handle.write("## 2. 拷贝数变异(重复)检测\n\n")
            else:
                output_handle.write("## 2. Copy Number Variation (Duplication) Detection\n\n")
            
            duplications = detect_duplication(alignments, 
                                            params.get('min_duplication_length', 50),
                                            params.get('min_duplication_similarity', 0.85))
            
            results['duplications'] = duplications
            
            if duplications:
                for i, duplication in enumerate(duplications, 1):
                    output_handle.write(f"[2.{i}] {format_duplication_report(duplication, use_chinese)}\n\n")
            else:
                if use_chinese:
                    output_handle.write("未检测到拷贝数变异(重复)。\n\n")
                else:
                    output_handle.write("No duplications detected.\n\n")
        
        # TODO: 添加更多结构变异类型的检测
                
        # 输出摘要
        if use_chinese:
            output_handle.write("# 检测结果摘要 #\n")
            output_handle.write(f"缺失: {len(results.get('deletions', []))} 个\n")
            output_handle.write(f"重复: {len(results.get('duplications', []))} 个\n")
        else:
            output_handle.write("# Detection Summary #\n")
            output_handle.write(f"Deletions: {len(results.get('deletions', []))}\n")
            output_handle.write(f"Duplications: {len(results.get('duplications', []))}\n")
    
    finally:
        if output_file and output_handle != sys.stdout:
            output_handle.close()
            if use_chinese:
                print(f"结构变异检测结果已写入 {output_file}", file=sys.stderr)
            else:
                print(f"Structural variation detection results written to {output_file}", file=sys.stderr)
    
    return results
