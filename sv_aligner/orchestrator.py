import os
import sys
import time
from typing import Dict, List, Optional, Tuple

from sv_aligner.indexer import IndexBuilder
from sv_aligner.data_types import AlignmentSegment
from sv_aligner.fasta_parser import read_fasta

class Orchestrator:
    """协调比对工作流程的主类"""
    
    def __init__(self, params: Dict):
        """初始化协调器"""
        self.params = params
        self.reference_index = None
        self.reference_info = None
        self.index_builder = IndexBuilder()
    
    def load_reference_and_index(self, ref_path: str, index_path: Optional[str] = None) -> None:
        """加载参考序列和索引"""
        print(f"加载参考序列 {ref_path}...", file=sys.stderr)
        
        # 检查参考文件是否存在
        if not os.path.exists(ref_path):
            raise FileNotFoundError(f"参考序列文件不存在: {ref_path}")
        
        # 尝试加载现有索引
        if index_path and os.path.exists(index_path):
            try:
                self.reference_index, self.reference_info = IndexBuilder.load_index(index_path)
                # 验证索引是否有效
                if not self.reference_info:
                    print("警告: 加载的索引不包含任何参考序列信息。重新构建索引...", file=sys.stderr)
                    self._build_new_index(ref_path, index_path)
                else:
                    return
            except Exception as e:
                print(f"加载索引失败: {e}，将重新构建索引", file=sys.stderr)
                self._build_new_index(ref_path, index_path)
        else:
            # 构建新索引
            self._build_new_index(ref_path, index_path)
    
    def _build_new_index(self, ref_path: str, index_path: Optional[str] = None) -> None:
        """构建新的索引"""
        try:
            print(f"检测到序列文件 {ref_path}，准备构建索引...", file=sys.stderr)
            self.reference_index, self.reference_info = self.index_builder.build(
                ref_path, self.params['k'], self.params['w']
            )
            
            # 检查索引是否有效
            if not self.reference_info:
                raise RuntimeError(f"无法从 {ref_path} 创建有效的索引。请检查参考序列文件。")
                
            # 如果提供了索引路径，保存索引
            if index_path:
                self.index_builder.save_index(index_path)
        except Exception as e:
            raise RuntimeError(f"构建参考索引时出错: {e}")

    def align_query(self, query_path: str) -> List[AlignmentSegment]:
        """比对查询序列到参考序列"""
        if not self.reference_index or not self.reference_info:
            raise RuntimeError("必须先加载参考序列和索引")
            
        # 检查查询文件是否存在
        if not os.path.exists(query_path):
            raise FileNotFoundError(f"查询序列文件不存在: {query_path}")
            
        print(f"开始比对查询序列 {query_path}...", file=sys.stderr)
        
        start_time = time.time()
        final_alignments = []
        
        # 读取查询序列
        query_count = 0
        
        for q_name, q_seq in read_fasta(query_path):
            query_count += 1
            
            if query_count % 10 == 0:
                print(f"处理第 {query_count} 条查询序列 ({q_name})...", file=sys.stderr)
            
            # TODO: 实现完整的比对工作流程:
            # 1. 找到查询序列的minimizers
            # 2. 与参考序列索引匹配，生成种子
            # 3. 链接种子
            # 4. 执行局部比对
            # 5. 选择非重叠片段
            
            # 临时: 创建一个占位符比对结果
            if len(q_seq) >= self.params['min_segment_length']:
                # 确保reference_info至少包含一个有效的参考序列
                if not self.reference_info:
                    raise RuntimeError("参考序列索引无效，无法进行比对。")
                
                # 获取第一个参考序列的信息
                first_ref_id = next(iter(self.reference_info))
                ref_name, ref_len = self.reference_info[first_ref_id]
                
                # 修改为只映射查询序列前三分之一的碱基，以便留空间检测结构变异
                alignment_length = min(len(q_seq) // 3, 100)
                
                dummy_alignment = AlignmentSegment(
                    q_name=q_name,
                    q_len=len(q_seq),
                    q_st=0,
                    q_en=alignment_length,
                    r_name=ref_name,
                    r_len=ref_len,
                    r_st=0,
                    r_en=alignment_length,
                    strand='+',
                    score=50,
                    edit_distance=0,
                    cigar=f"{alignment_length}M"
                )
                final_alignments.append(dummy_alignment)
        
        if query_count == 0:
            print(f"警告: 在 {query_path} 中未找到任何查询序列。", file=sys.stderr)
        
        end_time = time.time()
        print(f"比对完成。处理了 {query_count} 条查询序列。用时: {end_time - start_time:.2f} 秒。", 
              file=sys.stderr)
        
        return final_alignments
