from dataclasses import dataclass, field
from typing import Optional, List, Tuple

@dataclass(frozen=True)
class ReferenceLocation:
    """参考序列中的位置"""
    ref_id: int    # 对应参考序列的索引
    pos: int       # 0-based起始位置
    strand: bool   # True表示正向链，False表示反向链

@dataclass(frozen=True)
class Seed:
    """查询序列和参考序列之间的种子匹配"""
    q_pos: int      # 查询序列位置
    ref_id: int     # 参考序列ID
    r_pos: int      # 参考序列位置
    strand: bool    # True表示同向，False表示反向
    ref_name: str = ""  # 可选的参考序列名称

    def __lt__(self, other):
        """定义排序顺序：先按查询位置，再按参考ID，最后按参考位置"""
        if not isinstance(other, Seed):
            return NotImplemented
        return (self.q_pos, self.ref_id, self.r_pos) < \
               (other.q_pos, other.ref_id, other.r_pos)

@dataclass
class Chain:
    """种子链"""
    seeds: List[Seed]  # 链中的种子
    q_start: int       # 查询序列起始位置
    q_end: int         # 查询序列结束位置
    r_start: int       # 参考序列起始位置
    r_end: int         # 参考序列结束位置
    strand: bool       # 链方向
    score: float       # 链得分
    ref_id: int = -1   # 参考序列ID (从seeds推导)

    def __post_init__(self):
        """初始化后处理：如果没有提供ref_id，从seeds中推导"""
        if self.ref_id == -1 and self.seeds:
            self.ref_id = self.seeds[0].ref_id

@dataclass
class AlignmentSegment:
    """代表最终的比对片段"""
    q_name: str    # 查询序列名称
    q_len: int     # 查询序列长度
    q_st: int      # 查询序列起始位置
    q_en: int      # 查询序列结束位置
    r_name: str    # 参考序列名称
    r_len: int     # 参考序列长度
    r_st: int      # 参考序列起始位置
    r_en: int      # 参考序列结束位置
    strand: str    # '+' 或 '-'
    score: int     # 比对得分
    edit_distance: Optional[int] = None  # 编辑距离
    cigar: Optional[str] = None          # CIGAR字符串
    chain_score: Optional[float] = field(default=None, compare=False)  # 链得分(内部使用)

    def __lt__(self, other):
        """定义排序顺序：按查询起始位置"""
        if not isinstance(other, AlignmentSegment):
            return NotImplemented
        return self.q_st < other.q_st

    def to_tsv_line(self) -> str:
        """转换为TSV格式的行"""
        cigar_str = self.cigar if self.cigar else '*'
        ed_str = str(self.edit_distance) if self.edit_distance is not None else '*'

        return "\t".join(map(str, [
            self.q_name, self.q_len, self.q_st, self.q_en,
            self.r_name, self.r_len, self.r_st, self.r_en,
            self.strand, self.score, ed_str, cigar_str
        ]))
