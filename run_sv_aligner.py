#!/usr/bin/env python3
"""
SV-Aware Sequence Aligner 入口点脚本
"""
import os
import sys

# 添加项目根目录到Python搜索路径
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

# 导入主模块并运行
from sv_aligner.main import main

if __name__ == "__main__":
    main()
