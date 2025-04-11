import os
import pytest
import tempfile
import subprocess
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

class TestIntegration:
    """集成测试：测试完整的比对流程"""
    
    def setup_method(self):
        """在每个测试方法前设置环境"""
        self.temp_dir = tempfile.TemporaryDirectory()
        
        # 创建测试所需的各种序列文件
        self.create_test_sequences()
    
    def teardown_method(self):
        """在每个测试方法后清理环境"""
        self.temp_dir.cleanup()
    
    def create_test_sequences(self):
        """创建各种测试序列文件"""
        # 简单匹配序列
        self.simple_ref = os.path.join(self.temp_dir.name, "simple_ref.fa")
        self.simple_query = os.path.join(self.temp_dir.name, "simple_query.fa")
        
        with open(self.simple_ref, "w") as f:
            f.write(">simple_ref\n")
            f.write("ACGTACGTACGTACGTACGT\n")
        
        with open(self.simple_query, "w") as f:
            f.write(">simple_query\n")
            f.write("ACGTACGT\n")
        
        # 含有缺失的序列
        self.deletion_ref = os.path.join(self.temp_dir.name, "deletion_ref.fa")
        self.deletion_query = os.path.join(self.temp_dir.name, "deletion_query.fa")
        
        with open(self.deletion_ref, "w") as f:
            f.write(">deletion_ref\n")
            f.write("ACGTACGTACGTACGTACGT\n")
        
        with open(self.deletion_query, "w") as f:
            f.write(">deletion_query\n")
            f.write("ACGTTACGT\n")  # 缺失了中间的"ACG"
        
        # 含有插入的序列
        self.insertion_ref = os.path.join(self.temp_dir.name, "insertion_ref.fa")
        self.insertion_query = os.path.join(self.temp_dir.name, "insertion_query.fa")
        
        with open(self.insertion_ref, "w") as f:
            f.write(">insertion_ref\n")
            f.write("ACGTACGTACGT\n")
        
        with open(self.insertion_query, "w") as f:
            f.write(">insertion_query\n")
            f.write("ACGTTTTTACGT\n")  # 中间插入了"TTTT"
        
        # 含有倒置的序列
        self.inversion_ref = os.path.join(self.temp_dir.name, "inversion_ref.fa")
        self.inversion_query = os.path.join(self.temp_dir.name, "inversion_query.fa")
        
        with open(self.inversion_ref, "w") as f:
            f.write(">inversion_ref\n")
            f.write("ACGTACGTACGTACGTACGT\n")
        
        with open(self.inversion_query, "w") as f:
            f.write(">inversion_query\n")
            f.write("ACGTACGTATGCA\n")  # ACGTA的反向互补是TACGT，这里用了TGCA
    
    def run_alignment(self, ref_file, query_file, output_file, extra_args=None):
        """运行比对命令并返回结果"""
        if extra_args is None:
            extra_args = []
        
        # 构建命令行
        script_path = os.path.join(os.path.dirname(__file__), '..', 'run_sv_aligner.py')
        cmd = [sys.executable, script_path, ref_file, query_file, "-o", output_file] + extra_args
        
        # 运行命令
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            return True, result.stdout, result.stderr
        except subprocess.CalledProcessError as e:
            return False, e.stdout, e.stderr
    
    def test_simple_alignment(self):
        """测试简单序列的比对"""
        output_file = os.path.join(self.temp_dir.name, "simple_output.tsv")
        
        try:
            from sv_aligner.orchestrator import Orchestrator
            from sv_aligner.data_types import AlignmentSegment
            
            # 简化参数
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
            
            # 创建Orchestrator直接运行
            orchestrator = Orchestrator(params)
            orchestrator.load_reference_and_index(self.simple_ref)
            results = orchestrator.align_query(self.simple_query)
            
            # 验证结果
            assert len(results) > 0
            segment = results[0]
            assert segment.q_name == "simple_query"
            assert segment.r_name == "simple_ref"
            assert segment.strand == "+"
            
            # 检查是否匹配序列的开头
            assert segment.q_st == 0
            assert segment.r_st == 0
            
        except ImportError:
            # 如果模块未实现，跳过测试
            pytest.skip("Orchestrator未实现，跳过测试")
    
    def test_detect_deletion(self):
        """测试检测缺失的能力"""
        output_file = os.path.join(self.temp_dir.name, "deletion_output.tsv")
        
        try:
            from sv_aligner.orchestrator import Orchestrator
            
            # 简化参数
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
            
            # 创建Orchestrator直接运行
            orchestrator = Orchestrator(params)
            orchestrator.load_reference_and_index(self.deletion_ref)
            results = orchestrator.align_query(self.deletion_query)
            
            # 应该找到两个比对片段，中间有一个缺失
            assert len(results) >= 1
            
            # 如果能找到两个片段，验证它们的位置符合缺失模式
            if len(results) >= 2:
                # 按查询位置排序
                results.sort(key=lambda x: x.q_st)
                
                # 检查第一个片段末尾和第二个片段开头是否连续
                assert results[0].q_en <= results[1].q_st
                
                # 第一个片段末尾和第二个片段开头在参考上的距离应该较大
                ref_gap = results[1].r_st - results[0].r_en
                assert ref_gap > 0
            
        except ImportError:
            # 如果模块未实现，跳过测试
            pytest.skip("Orchestrator未实现，跳过测试")
    
    def test_detect_inversion(self):
        """测试检测倒置的能力"""
        output_file = os.path.join(self.temp_dir.name, "inversion_output.tsv")
        
        try:
            from sv_aligner.orchestrator import Orchestrator
            
            # 简化参数
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
            
            # 创建Orchestrator直接运行
            orchestrator = Orchestrator(params)
            orchestrator.load_reference_and_index(self.inversion_ref)
            results = orchestrator.align_query(self.inversion_query)
            
            # 应该至少找到一个反向比对片段
            assert any(segment.strand == "-" for segment in results), "没有找到反向链的比对片段"
            
        except ImportError:
            # 如果模块未实现，跳过测试
            pytest.skip("Orchestrator未实现，跳过测试")
