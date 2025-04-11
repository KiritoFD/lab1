import pytest
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sv_aligner.data_types import ReferenceLocation, Seed, AlignmentSegment

class TestDataTypes:
    """测试数据类型的基本功能"""
    
    def test_reference_location(self):
        """测试ReferenceLocation类"""
        loc = ReferenceLocation(ref_id=0, pos=10, strand=True)
        assert loc.ref_id == 0
        assert loc.pos == 10
        assert loc.strand is True
        
        # 测试不可变性
        with pytest.raises(AttributeError):
            loc.ref_id = 1
    
    def test_seed(self):
        """测试Seed类"""
        seed = Seed(q_pos=5, ref_id=0, r_pos=15, strand=True)
        assert seed.q_pos == 5
        assert seed.ref_id == 0
        assert seed.r_pos == 15
        assert seed.strand is True
        
        # 测试比较操作
        seed1 = Seed(q_pos=5, ref_id=0, r_pos=15, strand=True)
        seed2 = Seed(q_pos=10, ref_id=0, r_pos=20, strand=True)
        assert seed1 < seed2  # 基于q_pos比较
        
        seed3 = Seed(q_pos=5, ref_id=1, r_pos=15, strand=True)
        assert seed1 < seed3  # 相同q_pos, 基于ref_id比较
        
        seed4 = Seed(q_pos=5, ref_id=0, r_pos=20, strand=True)
        assert seed1 < seed4  # 相同q_pos和ref_id, 基于r_pos比较
    
    def test_alignment_segment(self):
        """测试AlignmentSegment类"""
        segment = AlignmentSegment(
            q_name="query1", q_len=100, q_st=10, q_en=20,
            r_name="ref1", r_len=200, r_st=50, r_en=60,
            strand="+", score=40, edit_distance=2, cigar="10M"
        )
        
        # 测试基本属性
        assert segment.q_name == "query1"
        assert segment.q_len == 100
        assert segment.q_st == 10
        assert segment.q_en == 20
        assert segment.r_name == "ref1"
        assert segment.r_len == 200
        assert segment.r_st == 50
        assert segment.r_en == 60
        assert segment.strand == "+"
        assert segment.score == 40
        assert segment.edit_distance == 2
        assert segment.cigar == "10M"
        
        # 测试比较操作
        segment1 = AlignmentSegment(q_name="query1", q_len=100, q_st=10, q_en=20,
                                    r_name="ref1", r_len=200, r_st=50, r_en=60,
                                    strand="+", score=40)
        segment2 = AlignmentSegment(q_name="query1", q_len=100, q_st=30, q_en=40,
                                    r_name="ref1", r_len=200, r_st=70, r_en=80,
                                    strand="+", score=50)
        
        assert segment1 < segment2  # 基于q_st比较
        
        # 测试TSV格式化
        tsv_line = segment.to_tsv_line()
        assert "query1" in tsv_line
        assert "10" in tsv_line
        assert "20" in tsv_line
        assert "ref1" in tsv_line
        assert "+" in tsv_line
        assert "40" in tsv_line
        assert "2" in tsv_line
        assert "10M" in tsv_line
