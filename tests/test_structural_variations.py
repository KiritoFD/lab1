import os
import pytest
import tempfile
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sv_aligner.utils import reverse_complement

class TestStructuralVariations:
    """测试各类结构变异的检测能力"""
    
    def setup_method(self):
        """在每个测试方法前设置环境"""
        self.temp_dir = tempfile.TemporaryDirectory()
        
        # 创建参考序列
        self.ref_file = os.path.join(self.temp_dir.name, "reference.fa")
        with open(self.ref_file, "w") as f:
            f.write(">reference\n")
            f.write("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")
    
    def teardown_method(self):
        """在每个测试方法后清理环境"""
        self.temp_dir.cleanup()
    
    def test_deletion(self):
        """测试检测大型缺失的能力"""
        # 创建一个包含大型缺失的查询序列
        query_file = os.path.join(self.temp_dir.name, "deletion.fa")
        with open(query_file, "w") as f:
            f.write(">deletion\n")
            # 从参考序列中删除中间的10个字符
            f.write("ACGTACGTACGTACGTACGTTACGTACGTACGT\n")
        
        try:
            from sv_aligner.orchestrator import Orchestrator
            
            params = {
                'k': 3,
                'w': 1,
                'match_score': 2,
                'mismatch_penalty': -3,
                'gap_open_penalty': -5,
                'gap_extend_penalty': -2,
                'min_segment_length': 4,
                'min_alignment_score': 5,
                'min_chain_score': 5
            }
            
            orchestrator = Orchestrator(params)
            orchestrator.load_reference_and_index(self.ref_file)
            results = orchestrator.align_query(query_file)
            
            # 按查询序列位置排序
            results.sort(key=lambda x: x.q_st)
            
            # 收集查询片段和对应参考片段
            query_segments = [(segment.q_st, segment.q_en) for segment in results]
            ref_segments = [(segment.r_st, segment.r_en) for segment in results]
            
            # 验证存在两个查询位置相连但参考位置有间隔的片段
            found_deletion = False
            for i in range(len(results) - 1):
                query_gap = results[i+1].q_st - results[i].q_en
                ref_gap = results[i+1].r_st - results[i].r_en
                
                if query_gap <= 1 and ref_gap > 5:  # 查询连续但参考有间隔
                    found_deletion = True
                    break
            
            assert found_deletion, "没有检测到预期的缺失"
            
        except ImportError:
            pytest.skip("Orchestrator未实现，跳过测试")
    
    def test_insertion(self):
        """测试检测大型插入的能力"""
        # 创建一个包含大型插入的查询序列
        query_file = os.path.join(self.temp_dir.name, "insertion.fa")
        with open(query_file, "w") as f:
            f.write(">insertion\n")
            # 在中间插入10个字符
            f.write("ACGTACGTACGTAAAAAAAAAATACGTACGTACGT\n")
        
        try:
            from sv_aligner.orchestrator import Orchestrator
            
            params = {
                'k': 3,
                'w': 1,
                'match_score': 2,
                'mismatch_penalty': -3,
                'gap_open_penalty': -5,
                'gap_extend_penalty': -2,
                'min_segment_length': 4,
                'min_alignment_score': 5,
                'min_chain_score': 5
            }
            
            orchestrator = Orchestrator(params)
            orchestrator.load_reference_and_index(self.ref_file)
            results = orchestrator.align_query(query_file)
            
            # 按查询序列位置排序
            results.sort(key=lambda x: x.q_st)
            
            # 验证存在两个参考位置相连但查询位置有间隔的片段
            found_insertion = False
            for i in range(len(results) - 1):
                query_gap = results[i+1].q_st - results[i].q_en
                ref_gap = results[i+1].r_st - results[i].r_en
                
                if ref_gap <= 1 and query_gap > 5:  # 参考连续但查询有间隔
                    found_insertion = True
                    break
            
            assert found_insertion, "没有检测到预期的插入"
            
        except ImportError:
            pytest.skip("Orchestrator未实现，跳过测试")
    
    def test_inversion(self):
        """测试检测倒置的能力"""
        # 创建一个包含倒置的查询序列
        query_file = os.path.join(self.temp_dir.name, "inversion.fa")
        
        # 从参考序列中取一段，反向互补，然后放回
        reference_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        prefix = reference_seq[:10]
        inverted_region = reverse_complement(reference_seq[10:30])
        suffix = reference_seq[30:]
        
        query_seq = prefix + inverted_region + suffix
        
        with open(query_file, "w") as f:
            f.write(">inversion\n")
            f.write(query_seq + "\n")
        
        try:
            from sv_aligner.orchestrator import Orchestrator
            
            params = {
                'k': 3,
                'w': 1,
                'match_score': 2,
                'mismatch_penalty': -3,
                'gap_open_penalty': -5,
                'gap_extend_penalty': -2,
                'min_segment_length': 4,
                'min_alignment_score': 5,
                'min_chain_score': 5
            }
            
            orchestrator = Orchestrator(params)
            orchestrator.load_reference_and_index(self.ref_file)
            results = orchestrator.align_query(query_file)
            
            # 检查是否有反向链的片段
            assert any(segment.strand == "-" for segment in results), "没有检测到预期的倒置"
            
        except ImportError:
            pytest.skip("Orchestrator未实现，跳过测试")
    
    def test_translocation(self):
        """测试检测易位的能力"""
        query_file = os.path.join(self.temp_dir.name, "translocation.fa")
        
        # 创建一个包含易位的查询序列（位置顺序颠倒）
        reference_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        part1 = reference_seq[:10]
        part2 = reference_seq[30:40]
        part3 = reference_seq[10:30]
        
        query_seq = part1 + part2 + part3
        
        with open(query_file, "w") as f:
            f.write(">translocation\n")
            f.write(query_seq + "\n")
        
        try:
            from sv_aligner.orchestrator import Orchestrator
            
            params = {
                'k': 3,
                'w': 1,
                'match_score': 2,
                'mismatch_penalty': -3,
                'gap_open_penalty': -5,
                'gap_extend_penalty': -2,
                'min_segment_length': 4,
                'min_alignment_score': 5,
                'min_chain_score': 5
            }
            
            orchestrator = Orchestrator(params)
            orchestrator.load_reference_and_index(self.ref_file)
            results = orchestrator.align_query(query_file)
            
            # 按查询序列位置排序
            results.sort(key=lambda x: x.q_st)
            
            # 检查参考序列上的位置顺序是否被打乱
            if len(results) >= 3:
                r_positions = [segment.r_st for segment in results]
                # 如果不是单调递增的，说明检测到了易位
                assert r_positions != sorted(r_positions), "没有检测到预期的易位"
            
        except ImportError:
            pytest.skip("Orchestrator未实现，跳过测试")
    
    def test_duplication(self):
        """测试检测重复的能力"""
        query_file = os.path.join(self.temp_dir.name, "duplication.fa")
        
        # 创建一个包含重复的查询序列
        reference_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        part1 = reference_seq[:10]
        part2 = reference_seq[10:20]
        
        # 重复part2
        query_seq = part1 + part2 + part2
        
        with open(query_file, "w") as f:
            f.write(">duplication\n")
            f.write(query_seq + "\n")
        
        try:
            from sv_aligner.orchestrator import Orchestrator
            
            params = {
                'k': 3,
                'w': 1,
                'match_score': 2,
                'mismatch_penalty': -3,
                'gap_open_penalty': -5,
                'gap_extend_penalty': -2,
                'min_segment_length': 4,
                'min_alignment_score': 5,
                'min_chain_score': 5
            }
            
            orchestrator = Orchestrator(params)
            orchestrator.load_reference_and_index(self.ref_file)
            results = orchestrator.align_query(query_file)
            
            # 检查是否有多个查询片段映射到同一参考片段
            ref_regions = [(segment.r_st, segment.r_en) for segment in results]
            
            # 如果有重复映射，应该存在至少一个参考区域在列表中出现多次
            has_duplicates = len(ref_regions) != len(set(ref_regions))
            
            if not has_duplicates:
                # 另一种检测方法：两个查询片段映射到重叠的参考区域
                found_overlap = False
                for i in range(len(results)):
                    for j in range(i+1, len(results)):
                        r1_st, r1_en = results[i].r_st, results[i].r_en
                        r2_st, r2_en = results[j].r_st, results[j].r_en
                        
                        # 检查参考区域是否有重叠
                        if max(r1_st, r2_st) < min(r1_en, r2_en):
                            found_overlap = True
                            break
                    if found_overlap:
                        break
                        
                has_duplicates = found_overlap
            
            assert has_duplicates, "没有检测到预期的重复"
            
        except ImportError:
            pytest.skip("Orchestrator未实现，跳过测试")
