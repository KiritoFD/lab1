import pytest
from unittest.mock import patch
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sv_aligner.seeder import Seeder
from sv_aligner.data_types import ReferenceLocation, Seed

class TestSeeder:
    def setup_method(self):
        """在每个测试方法前设置环境"""
        # 创建一个测试用的参考索引
        self.ref_index = {
            123: [ReferenceLocation(ref_id=0, pos=0, strand=True),    # 假设哈希123对应ACGT
                  ReferenceLocation(ref_id=0, pos=6, strand=True)],
            456: [ReferenceLocation(ref_id=0, pos=3, strand=True),    # 假设哈希456对应TACG
                  ReferenceLocation(ref_id=1, pos=0, strand=False)],
            789: [ReferenceLocation(ref_id=1, pos=4, strand=True)]    # 假设哈希789对应GCTA
        }
        
        # 创建参考序列信息
        self.ref_info = {
            0: ("ref1", 12),  # ref1长度为12
            1: ("ref2", 8)    # ref2长度为8
        }
    
    def test_find_seeds_basic(self):
        """测试基本的种子识别功能"""
        seeder = Seeder(self.ref_index, self.ref_info)
        
        # 使用猴子补丁替换get_minimizers函数
        with patch('sv_aligner.seeder.get_minimizers') as mock_get_minimizers:
            # 设置mock返回值
            mock_get_minimizers.return_value = [(123, 0, True), (456, 3, True), (123, 6, True)]
            
            seeds = seeder.find_seeds("ACGTACGTACGT", k=4, w=1)
            
            # 验证找到的种子数量
            assert len(seeds) == 4  # 应该找到4个种子（2个来自哈希123，2个来自哈希456）
            
            # 验证种子内容
            seed_tuples = [(seed.q_pos, seed.ref_id, seed.r_pos, seed.strand) for seed in seeds]
            assert (0, 0, 0, True) in seed_tuples  # 查询位置0 -> 参考序列0位置0，正向链
            assert (0, 0, 6, True) in seed_tuples  # 查询位置0 -> 参考序列0位置6，正向链
            assert (3, 0, 3, True) in seed_tuples  # 查询位置3 -> 参考序列0位置3，正向链
            assert (3, 1, 0, False) in seed_tuples  # 查询位置3 -> 参考序列1位置0，反向链
    
    def test_find_seeds_empty(self):
        """测试空序列的种子识别"""
        seeder = Seeder(self.ref_index, self.ref_info)
        
        seeds = seeder.find_seeds("", k=4, w=1)
        assert len(seeds) == 0  # 空序列应该没有种子
    
    def test_find_seeds_no_match(self):
        """测试没有匹配的序列"""
        seeder = Seeder(self.ref_index, self.ref_info)
        
        # 使用猴子补丁替换get_minimizers函数
        with patch('sv_aligner.seeder.get_minimizers') as mock_get_minimizers:
            # 设置mock返回值，返回不在索引中的哈希
            mock_get_minimizers.return_value = [(999, 0, True), (888, 3, True)]
            
            seeds = seeder.find_seeds("NNNNNNNN", k=4, w=1)
            assert len(seeds) == 0  # 应该没有找到种子
    
    def test_seeds_sorting(self):
        """测试种子排序功能"""
        seeder = Seeder(self.ref_index, self.ref_info)
        
        # 使用猴子补丁替换get_minimizers函数
        with patch('sv_aligner.seeder.get_minimizers') as mock_get_minimizers:
            # 设置mock返回值，返回无序的位置
            mock_get_minimizers.return_value = [(123, 6, True), (456, 3, True), (123, 0, True)]
            
            seeds = seeder.find_seeds("ACGTACGTACGT", k=4, w=1)
            
            # 验证种子已按查询位置排序
            assert all(seeds[i].q_pos <= seeds[i+1].q_pos for i in range(len(seeds)-1))
