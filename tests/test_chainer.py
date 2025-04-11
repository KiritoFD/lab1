import pytest
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sv_aligner.chainer import Chainer
from sv_aligner.data_types import Seed

class TestChainer:
    def setup_method(self):
        """在每个测试方法前设置环境"""
        # 创建一个基本的参数字典
        self.params = {
            'k': 4,
            'gap_penalty_factor': 0.01,
            'dist_diff_penalty_factor': 0.05,
            'max_gap': 1000,
            'max_dist_diff': 200,
            'min_chain_score': 10
        }
        
        # 创建一个简单的Chainer实例
        self.chainer = Chainer(self.params)
    
    def test_calculate_gap_penalty(self):
        """测试间隙惩罚计算功能"""
        # 创建两个种子
        seed_i = Seed(q_pos=10, ref_id=0, r_pos=100, strand=True)
        seed_j = Seed(q_pos=20, ref_id=0, r_pos=110, strand=True)
        
        # 计算间隙惩罚
        penalty = self.chainer._calculate_gap_penalty(seed_i, seed_j)
        
        # 期望的惩罚 = gap_penalty_factor * (q_dist + r_dist) + dist_diff_penalty_factor * dist_diff
        # q_dist = 10, r_dist = 10, dist_diff = 0
        expected_penalty = 0.01 * (10 + 10) + 0.05 * 0
        assert penalty == expected_penalty
        
        # 测试不同方向的种子
        seed_k = Seed(q_pos=30, ref_id=0, r_pos=120, strand=False)  # 反向链
        penalty = self.chainer._calculate_gap_penalty(seed_i, seed_k)
        assert penalty == float('inf')  # 不同方向应该得到无穷大的惩罚
    
    def test_chain_seeds_simple(self):
        """测试简单情况下的种子链接"""
        # 创建一组共线性的种子
        seeds = [
            Seed(q_pos=0, ref_id=0, r_pos=100, strand=True),
            Seed(q_pos=10, ref_id=0, r_pos=110, strand=True),
            Seed(q_pos=20, ref_id=0, r_pos=120, strand=True),
            Seed(q_pos=30, ref_id=0, r_pos=130, strand=True)
        ]
        
        # 链接种子
        chains = self.chainer.chain_seeds(seeds)
        
        # 应该只有一条链，包含所有种子
        assert len(chains) == 1
        assert len(chains[0].seeds) == 4
        
        # 验证链的属性
        chain = chains[0]
        assert chain.q_start == 0
        assert chain.q_end == 34  # 30 + 4 (假设种子长度为参数k值)
        assert chain.r_start == 100
        assert chain.r_end == 134
        assert chain.strand is True
        assert chain.score > self.params['min_chain_score']
    
    def test_chain_seeds_multiple_chains(self):
        """测试识别多条独立链的功能"""
        # 创建两组独立的共线性种子
        seeds = [
            # 第一组
            Seed(q_pos=0, ref_id=0, r_pos=100, strand=True),
            Seed(q_pos=10, ref_id=0, r_pos=110, strand=True),
            
            # 第二组（与第一组不共线，因为在不同的参考序列上）
            Seed(q_pos=50, ref_id=1, r_pos=200, strand=True),
            Seed(q_pos=60, ref_id=1, r_pos=210, strand=True)
        ]
        
        # 链接种子
        chains = self.chainer.chain_seeds(seeds)
        
        # 应该有两条独立的链
        assert len(chains) == 2
        
        # 验证第一条链
        assert chains[0].ref_id == seeds[0].ref_id
        assert len(chains[0].seeds) == 2
        
        # 验证第二条链
        assert chains[1].ref_id == seeds[2].ref_id
        assert len(chains[1].seeds) == 2
    
    def test_chain_seeds_with_reversal(self):
        """测试处理倒置的能力"""
        # 创建一组包含倒置的种子
        seeds = [
            # 正向链种子
            Seed(q_pos=0, ref_id=0, r_pos=100, strand=True),
            Seed(q_pos=10, ref_id=0, r_pos=110, strand=True),
            
            # 反向链种子
            Seed(q_pos=30, ref_id=0, r_pos=200, strand=False),
            Seed(q_pos=40, ref_id=0, r_pos=190, strand=False)  # r_pos递减表示反向链
        ]
        
        # 链接种子
        chains = self.chainer.chain_seeds(seeds)
        
        # 应该有两条独立的链（一条正向，一条反向）
        assert len(chains) == 2
        
        # 找到正向链和反向链
        forward_chain = [c for c in chains if c.strand is True][0]
        reverse_chain = [c for c in chains if c.strand is False][0]
        
        # 验证正向链
        assert len(forward_chain.seeds) == 2
        assert forward_chain.seeds[0].q_pos == 0
        assert forward_chain.seeds[1].q_pos == 10
        
        # 验证反向链
        assert len(reverse_chain.seeds) == 2
        assert reverse_chain.seeds[0].q_pos == 30
        assert reverse_chain.seeds[1].q_pos == 40
    
    def test_chain_seeds_with_gap(self):
        """测试处理有间隙的种子链"""
        # 创建一组种子，中间有一个大间隙
        seeds = [
            Seed(q_pos=0, ref_id=0, r_pos=100, strand=True),
            Seed(q_pos=10, ref_id=0, r_pos=110, strand=True),
            # 大间隙
            Seed(q_pos=500, ref_id=0, r_pos=600, strand=True),  # 间隙内但太远
            Seed(q_pos=520, ref_id=0, r_pos=620, strand=True)
        ]
        
        # 设置较小的最大间隙参数
        small_gap_params = dict(self.params)
        small_gap_params['max_gap'] = 100  # 只允许小于100的间隙
        chainer_small_gap = Chainer(small_gap_params)
        
        # 链接种子
        chains = chainer_small_gap.chain_seeds(seeds)
        
        # 应该有两条独立的链
        assert len(chains) == 2
        
        # 验证链内种子数量
        assert len(chains[0].seeds) == 2
        assert len(chains[1].seeds) == 2
        
        # 现在测试大间隙参数
        large_gap_params = dict(self.params)
        large_gap_params['max_gap'] = 1000  # 允许大于500的间隙
        chainer_large_gap = Chainer(large_gap_params)
        
        # 链接种子
        chains = chainer_large_gap.chain_seeds(seeds)
        
        # 如果间隙参数足够大，可能会形成单链，取决于gap_penalty_factor的具体实现
        # 这里我们不对链的数量做硬性断言，而是检查是否所有种子都被链接
        total_seeds = sum(len(chain.seeds) for chain in chains)
        assert total_seeds == 4
    
    def test_empty_seeds(self):
        """测试处理空种子列表"""
        chains = self.chainer.chain_seeds([])
        assert len(chains) == 0
