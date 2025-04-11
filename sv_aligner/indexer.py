import pickle
import sys
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

from sv_aligner.utils import get_minimizers
from sv_aligner.fasta_parser import read_fasta
from sv_aligner.data_types import ReferenceLocation

class IndexBuilder:
    """构建参考序列的索引"""
    
    def __init__(self):
        self.reference_index: Dict[int, List[ReferenceLocation]] = defaultdict(list)
        self.reference_info: Dict[int, Tuple[str, int]] = {}  # ref_id -> (name, length)

    def build(self, ref_fasta_path: str, k: int, w: int) -> Tuple[Dict[int, List[ReferenceLocation]], Dict[int, Tuple[str, int]]]:
        """从参考FASTA或纯序列文件构建minimizer索引"""
        print(f"构建索引 {ref_fasta_path}，k={k}, w={w}...", file=sys.stderr)
        ref_id_counter = 0
        
        # 添加文件存在检查
        try:
            sequence_count = 0
            for header, seq in read_fasta(ref_fasta_path):
                sequence_count += 1
                if not seq:
                    print(f"警告: 跳过空序列 '{header}'", file=sys.stderr)
                    continue
                
                # 检查序列是否包含有效的DNA字符
                valid_chars = set("ACGTN")
                upper_seq = seq.upper()
                invalid_chars = set(upper_seq) - valid_chars
                if invalid_chars:
                    print(f"警告: 序列 '{header}' 包含非DNA字符: {', '.join(invalid_chars)}", 
                          file=sys.stderr)
                    # 过滤非标准字符
                    seq = ''.join(c if c.upper() in valid_chars else 'N' for c in seq)
                    print(f"已将非标准字符替换为'N'", file=sys.stderr)

                ref_name = header  # 已在read_fasta中提取头部第一部分
                self.reference_info[ref_id_counter] = (ref_name, len(seq))
                print(f"  索引 {ref_name} (ID: {ref_id_counter}, 长度: {len(seq)})...", file=sys.stderr)

                minimizers = get_minimizers(seq, k, w)
                
                print(f"  找到 {len(minimizers)} 个minimizers", file=sys.stderr)

                for mini_hash, pos, is_forward_strand in minimizers:
                    # 存储位置
                    location = ReferenceLocation(ref_id=ref_id_counter, pos=pos, strand=is_forward_strand)
                    self.reference_index[mini_hash].append(location)

                ref_id_counter += 1
                
            if sequence_count == 0:
                raise ValueError(f"未在 {ref_fasta_path} 中找到任何序列。请检查文件格式。")
                
            print(f"索引构建完成。已索引 {ref_id_counter} 个序列。", file=sys.stderr)
            print(f"  索引中唯一minimizers总数: {len(self.reference_index)}", file=sys.stderr)
            
            return self.reference_index, self.reference_info
            
        except FileNotFoundError:
            raise FileNotFoundError(f"参考序列文件未找到: {ref_fasta_path}")
        except Exception as e:
            raise RuntimeError(f"构建索引时出错: {e}")

    def save_index(self, index_path: str):
        """保存索引到文件"""
        try:
            with open(index_path, 'wb') as outfile:
                pickle.dump((self.reference_index, self.reference_info), outfile)
            print(f"索引已保存到 {index_path}", file=sys.stderr)
        except Exception as e:
            print(f"保存索引到 {index_path} 时出错: {e}", file=sys.stderr)

    @staticmethod
    def load_index(index_path: str) -> Tuple[Dict[int, List[ReferenceLocation]], Dict[int, Tuple[str, int]]]:
        """从文件加载索引"""
        print(f"从 {index_path} 加载索引...", file=sys.stderr)
        try:
            with open(index_path, 'rb') as infile:
                reference_index, reference_info = pickle.load(infile)
            print("索引加载成功。", file=sys.stderr)
            return reference_index, reference_info
        except FileNotFoundError:
            print(f"错误: 索引文件未找到 {index_path}", file=sys.stderr)
            raise
        except Exception as e:
            print(f"从 {index_path} 加载索引时出错: {e}", file=sys.stderr)
            raise
