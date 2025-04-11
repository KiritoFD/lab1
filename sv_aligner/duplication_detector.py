"""
拷贝数变异(重复)检测模块
基于哈希表和序列相似性方法检测拷贝数变异
直接从C代码翻译而来
"""
import time
import sys
from typing import List, Dict, Tuple, Optional, Set
from collections import defaultdict
from sv_aligner.data_types import AlignmentSegment
from sv_aligner.utils import reverse_complement

# 重复片段的数据结构 - 对应C代码中的RepeatPattern
class DuplicationPattern:
    """表示检测到的重复片段"""
    def __init__(self, query_name: str, query_len: int,
                 q_st1: int, q_en1: int, q_st2: int, q_en2: int,
                 r_st1: int, r_en1: int, r_st2: int, r_en2: int,
                 r_name: str, r_len: int, 
                 strand1: str, strand2: str, 
                 length: int, similarity: float,
                 seg1: AlignmentSegment, seg2: AlignmentSegment):
        """初始化重复片段对象"""
        self.query_name = query_name    # 查询序列名称
        self.query_len = query_len      # 查询序列长度
        self.q_st1 = q_st1              # 第一个重复片段在查询序列上的起始位置
        self.q_en1 = q_en1              # 第一个重复片段在查询序列上的结束位置
        self.q_st2 = q_st2              # 第二个重复片段在查询序列上的起始位置
        self.q_en2 = q_en2              # 第二个重复片段在查询序列上的结束位置
        self.r_st1 = r_st1              # 第一个重复片段在参考序列上的起始位置
        self.r_en1 = r_en1              # 第一个重复片段在参考序列上的结束位置
        self.r_st2 = r_st2              # 第二个重复片段在参考序列上的起始位置
        self.r_en2 = r_en2              # 第二个重复片段在参考序列上的结束位置
        self.r_name = r_name            # 参考序列名称
        self.r_len = r_len              # 参考序列长度
        self.strand1 = strand1          # 第一个片段的链方向
        self.strand2 = strand2          # 第二个片段的链方向
        self.length = length            # 重复长度
        self.similarity = similarity    # 相似度
        self.segment1 = seg1            # 第一个比对片段
        self.segment2 = seg2            # 第二个比对片段
        
        # 以下属性对应C版本特有的属性
        self.position = r_st1           # 参考序列上的位置
        self.repeat_count = 2           # 重复次数，默认为2
        self.is_reverse = strand1 != strand2  # 是否为反向重复

    def __str__(self):
        """返回重复的字符串表示"""
        return (f"重复: {self.query_name} 位置1:{self.q_st1}-{self.q_en1} ({self.strand1}) "
                f"位置2:{self.q_st2}-{self.q_en2} ({self.strand2}), "
                f"长度:{self.length}bp, 相似度:{self.similarity:.2f}")


# 计算哈希值，对应C代码中的hash_function
def hash_function(s: str, size: int) -> int:
    """
    计算字符串的哈希值
    与C代码中的哈希函数保持一致
    """
    hash_val = 5381
    for c in s:
        hash_val = ((hash_val << 5) + hash_val) + ord(c)
    return hash_val % size


# DNA序列相似度计算，对应C代码中的calculate_similarity
def calculate_similarity(seq1: str, seq2: str) -> float:
    """
    计算两个DNA序列的相似度
    返回匹配碱基占总长度的比例
    """
    # 确保两个序列长度相同
    min_len = min(len(seq1), len(seq2))
    seq1 = seq1[:min_len]
    seq2 = seq2[:min_len]
    
    if min_len == 0:
        return 0.0
        
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return matches / min_len


# 模糊哈希函数，对应C代码中的fuzzy_hash_function
def fuzzy_hash_function(s: str, size: int) -> int:
    """
    生成对小变异不敏感的哈希值
    使用滑动窗口对局部特征进行哈希
    """
    hash_val = 5381
    window_size = 3  # k-mer大小，与C代码保持一致
    
    for i in range(max(0, len(s) - window_size + 1)):
        local_hash = 0
        for j in range(min(window_size, len(s) - i)):
            local_hash = ((local_hash << 5) + local_hash) + ord(s[i + j])
        hash_val = ((hash_val << 5) + hash_val) + local_hash
    
    return hash_val % size


# 计算两个比对片段之间的相似度
def calculate_segment_similarity(seg1: AlignmentSegment, seg2: AlignmentSegment) -> float:
    """
    估算两个比对片段之间的序列相似度
    使用编辑距离或比对分数进行估算
    """
    # 使用比对得分和编辑距离估算相似度
    if seg1.edit_distance is not None and seg2.edit_distance is not None:
        len1 = seg1.q_en - seg1.q_st
        len2 = seg2.q_en - seg2.q_st
        avg_len = (len1 + len2) / 2
        
        if avg_len <= 0:
            return 0.0
            
        # 归一化编辑距离
        edit_dist = (seg1.edit_distance + seg2.edit_distance) / 2
        similarity = 1.0 - (edit_dist / avg_len)
        return max(0.0, min(1.0, similarity))
    
    # 如果没有编辑距离，使用比对得分
    score1 = seg1.score if seg1.score is not None else 0
    score2 = seg2.score if seg2.score is not None else 0
    len1 = seg1.q_en - seg1.q_st
    len2 = seg2.q_en - seg2.q_st
    
    # 简单估计：假设完全匹配得分为2分/碱基
    max_possible_score = 2 * min(len1, len2)
    if max_possible_score <= 0:
        return 0.0
        
    return min(1.0, (score1 + score2) / (2 * max_possible_score))


# 对应C代码中的build_sequence_hashmap函数
def build_sequence_hashmap(sequence_dict: Dict[str, str], length: int) -> Dict[int, List[Tuple[str, int, str]]]:
    """
    构建序列片段的哈希映射
    
    Args:
        sequence_dict: 序列字典，键为序列名，值为序列内容
        length: 序列片段长度
        
    Returns:
        哈希值到位置信息的映射
    """
    hashmap = {}
    
    for seq_name, sequence in sequence_dict.items():
        seq_len = len(sequence)
        
        for i in range(max(0, seq_len - length + 1)):
            segment = sequence[i:i+length]
            h = fuzzy_hash_function(segment, 16384)  # 使用与C代码相同的哈希表大小
            
            if h not in hashmap:
                hashmap[h] = []
            hashmap[h].append((seq_name, i, segment))
    
    return hashmap


# 核心函数：检测重复，对应C代码中的find_repeats函数
def detect_duplication(alignments: List[AlignmentSegment], 
                       min_length: int = 50, 
                       min_similarity: float = 0.85,
                       max_length: int = 101) -> List[DuplicationPattern]:
    """
    检测查询序列中的拷贝数变异(重复)
    
    Args:
        alignments: 比对结果列表
        min_length: 最小重复长度，对应C代码中的MIN_LENGTH
        min_similarity: 最小相似度阈值
        max_length: 最大重复长度，对应C代码中的MAX_LENGTH
        
    Returns:
        重复片段列表
    """
    if len(alignments) < 2:
        print("DEBUG: 比对片段数量不足，无法检测重复", file=sys.stderr)
        return []
    
    print(f"DEBUG: 开始检测拷贝数变异，共有 {len(alignments)} 个比对片段...", file=sys.stderr)
    start_time = time.time()
    
    # 按参考序列位置排序
    ref_sorted = sorted(alignments, key=lambda a: (a.r_name, a.r_st, a.r_en))
    
    # 创建参考区域索引，使用离散化的bins
    bin_size = 50  # 对应C代码中的窗口大小
    ref_regions = defaultdict(list)
    
    for seg in ref_sorted:
        # 跳过过短的片段
        if seg.r_en - seg.r_st < min_length or seg.q_en - seg.q_st < min_length:
            continue
        
        # 计算当前片段覆盖的所有bins
        start_bin = seg.r_st // bin_size
        end_bin = (seg.r_en - 1) // bin_size + 1
        
        for bin_idx in range(start_bin, end_bin):
            bin_key = (seg.r_name, bin_idx, seg.strand)
            ref_regions[bin_key].append(seg)
    
    print(f"DEBUG: 创建了 {len(ref_regions)} 个参考位置bins", file=sys.stderr)
    
    # 存储检测到的重复
    duplications = []
    
    # 检查每个区域中的片段对，寻找潜在重复
    for bin_key, segments in ref_regions.items():
        if len(segments) < 2:
            continue
            
        # 按查询序列名称分组
        by_query = defaultdict(list)
        for seg in segments:
            by_query[seg.q_name].append(seg)
        
        # 对每个查询序列单独处理
        for q_name, q_segs in by_query.items():
            if len(q_segs) < 2:
                continue
                
            # 分析此查询序列中的所有segment对
            for i in range(len(q_segs)):
                for j in range(i+1, len(q_segs)):
                    seg1 = q_segs[i]
                    seg2 = q_segs[j]
                    
                    # 1. 计算参考序列重叠
                    r_overlap_start = max(seg1.r_st, seg2.r_st)
                    r_overlap_end = min(seg1.r_en, seg2.r_en)
                    r_overlap = r_overlap_end - r_overlap_start
                    
                    # 2. 检查查询序列区域是否不重叠（不同位置的拷贝）
                    q_overlap_start = max(seg1.q_st, seg2.q_st)
                    q_overlap_end = min(seg1.q_en, seg2.q_en)
                    q_overlap = max(0, q_overlap_end - q_overlap_start)
                    
                    # 3. 检查重复条件
                    # (a) 参考区域有足够重叠
                    # (b) 查询区域几乎不重叠
                    if r_overlap >= min_length and q_overlap <= 0.1 * min_length:
                        # 计算相似度
                        similarity = calculate_segment_similarity(seg1, seg2)
                        
                        if similarity >= min_similarity:
                            dup_length = min(seg1.q_en - seg1.q_st, seg2.q_en - seg2.q_st)
                            
                            # 限制重复长度范围
                            if min_length <= dup_length <= max_length:
                                # 避免重复添加
                                is_duplicate = False
                                for existing_dup in duplications:
                                    if (existing_dup.q_st1 == seg1.q_st and 
                                        existing_dup.q_en1 == seg1.q_en and
                                        existing_dup.q_st2 == seg2.q_st and
                                        existing_dup.q_en2 == seg2.q_en):
                                        is_duplicate = True
                                        break
                                
                                if not is_duplicate:
                                    print(f"DEBUG: 找到重复: q1:{seg1.q_st}-{seg1.q_en} 和 q2:{seg2.q_st}-{seg2.q_en}, "
                                          f"相似度: {similarity:.2f}", file=sys.stderr)
                                    
                                    duplication = DuplicationPattern(
                                        query_name=q_name,
                                        query_len=seg1.q_len,
                                        q_st1=seg1.q_st,
                                        q_en1=seg1.q_en,
                                        q_st2=seg2.q_st,
                                        q_en2=seg2.q_en,
                                        r_st1=seg1.r_st,
                                        r_en1=seg1.r_en,
                                        r_st2=seg2.r_st,
                                        r_en2=seg2.r_en,
                                        r_name=seg1.r_name,
                                        r_len=seg1.r_len,
                                        strand1=seg1.strand,
                                        strand2=seg2.strand,
                                        length=dup_length,
                                        similarity=similarity,
                                        seg1=seg1,
                                        seg2=seg2
                                    )
                                    duplications.append(duplication)
    
    # 按照重复长度和相似度排序，按照C代码中的quick_sort_repeats逻辑
    duplications.sort(key=lambda d: (d.length, d.similarity), reverse=True)
    
    end_time = time.time()
    print(f"DEBUG: 共检测到 {len(duplications)} 个重复，耗时: {(end_time - start_time)*1000:.2f}ms", file=sys.stderr)
    return duplications


# 格式化重复报告
def format_duplication_report(duplication: DuplicationPattern, use_chinese: bool = True) -> str:
    """
    格式化重复报告
    
    Args:
        duplication: 重复信息
        use_chinese: 是否使用中文输出
        
    Returns:
        格式化的重复报告文本
    """
    if use_chinese:
        report = [
            f"重复检测 (DUPLICATION):",
            f"  查询序列: {duplication.query_name}",
            f"  第一片段: 查询位置 {duplication.q_st1}-{duplication.q_en1} -> 参考位置 {duplication.r_name}:{duplication.r_st1}-{duplication.r_en1} ({duplication.strand1}链)",
            f"  第二片段: 查询位置 {duplication.q_st2}-{duplication.q_en2} -> 参考位置 {duplication.r_name}:{duplication.r_st2}-{duplication.r_en2} ({duplication.strand2}链)",
            f"  重复长度: {duplication.length}bp",
            f"  序列相似度: {duplication.similarity:.2f}",
        ]
    else:
        report = [
            f"DUPLICATION DETECTED:",
            f"  Query sequence: {duplication.query_name}",
            f"  First segment: Query {duplication.q_st1}-{duplication.q_en1} -> Reference {duplication.r_name}:{duplication.r_st1}-{duplication.r_en1} ({duplication.strand1} strand)",
            f"  Second segment: Query {duplication.q_st2}-{duplication.q_en2} -> Reference {duplication.r_name}:{duplication.r_st2}-{duplication.r_en2} ({duplication.strand2} strand)",
            f"  Duplication length: {duplication.length}bp",
            f"  Sequence similarity: {duplication.similarity:.2f}",
        ]
    
    return "\n".join(report)


# 检测并报告重复的完整流程
def detect_and_report_duplications(alignments: List[AlignmentSegment], 
                                  min_length: int = 50, 
                                  min_similarity: float = 0.85,
                                  output_file: Optional[str] = None,
                                  use_chinese: bool = True) -> List[DuplicationPattern]:
    """
    检测并报告重复
    
    Args:
        alignments: 比对结果列表
        min_length: 最小重复长度
        min_similarity: 最小相似度阈值
        output_file: 输出文件路径，如果为None则输出到标准输出
        use_chinese: 是否使用中文输出
        
    Returns:
        检测到的重复列表
    """
    # 检测重复
    start_time = time.time()
    duplications = detect_duplication(alignments, min_length, min_similarity)
    end_time = time.time()
    
    # 准备输出
    output_handle = open(output_file, 'w', encoding='utf-8') if output_file else sys.stdout
    try:
        # 输出摘要信息
        if use_chinese:
            output_handle.write(f"## 拷贝数变异(重复)检测结果 ##\n")
            output_handle.write(f"总检测到 {len(duplications)} 个重复 (最小长度: {min_length}bp, 最小相似度: {min_similarity:.2f})\n")
            output_handle.write(f"检测用时: {(end_time - start_time)*1000:.2f}毫秒\n\n")
        else:
            output_handle.write(f"## DUPLICATION DETECTION RESULTS ##\n")
            output_handle.write(f"Total detected: {len(duplications)} duplications (minimum length: {min_length}bp, minimum similarity: {min_similarity:.2f})\n")
            output_handle.write(f"Detection time: {(end_time - start_time)*1000:.2f}ms\n\n")
        
        # 输出每个重复的详细信息
        for i, duplication in enumerate(duplications, 1):
            output_handle.write(f"[{i}] {format_duplication_report(duplication, use_chinese)}\n\n")
    
    finally:
        if output_file and output_handle != sys.stdout:
            output_handle.close()
            if use_chinese:
                print(f"重复检测结果已写入 {output_file}", file=sys.stderr)
            else:
                print(f"Duplication detection results written to {output_file}", file=sys.stderr)
    
    return duplications


# 演示用函数，用于创建包含拷贝数变异的测试数据
def create_test_duplication_data() -> Tuple[List[AlignmentSegment], List[DuplicationPattern]]:
    """
    创建测试用的重复数据
    返回包含已知重复的比对片段列表
    """
    query_name = "test_query"
    query_len = 1000
    ref_name = "test_ref"
    ref_len = 2000
    
    # 创建两个比对片段，它们映射到参考序列的相同区域
    seg1 = AlignmentSegment(
        q_name=query_name,
        q_len=query_len,
        q_st=100,
        q_en=200,
        r_name=ref_name,
        r_len=ref_len,
        r_st=500,
        r_en=600,
        strand="+",
        score=100,
        edit_distance=0,
        cigar="100M"
    )
    
    seg2 = AlignmentSegment(
        q_name=query_name,
        q_len=query_len,
        q_st=300,
        q_en=400,
        r_name=ref_name,
        r_len=ref_len,
        r_st=500,
        r_en=600,
        strand="+",
        score=90,
        edit_distance=5,
        cigar="100M"
    )
    
    # 创建已知的重复信息
    duplication = DuplicationPattern(
        query_name=query_name,
        query_len=query_len,
        q_st1=100,
        q_en1=200,
        q_st2=300,
        q_en2=400,
        r_st1=500,
        r_en1=600,
        r_st2=500,
        r_en2=600,
        r_name=ref_name,
        r_len=ref_len,
        strand1="+",
        strand2="+",
        length=100,
        similarity=0.95,
        seg1=seg1,
        seg2=seg2
    )
    
    return [seg1, seg2], [duplication]
