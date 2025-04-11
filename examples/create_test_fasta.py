#!/usr/bin/env python3
"""
创建简单的测试FASTA文件
"""
import os
import sys

def create_test_reference():
    """创建一个简单的参考序列FASTA文件"""
    with open("ref.txt", "w") as f:
        f.write(">reference\n")
        f.write("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")
        f.write("TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\n")
    print(f"已创建参考序列文件: {os.path.abspath('ref.txt')}")

def create_test_query():
    """创建一个简单的查询序列FASTA文件"""
    with open("que.txt", "w") as f:
        f.write(">query1\n")
        f.write("ACGTACGTACGTACGTACGT\n")
        f.write(">query2\n")
        f.write("TGCATGCATGCATGCATGCA\n")
    print(f"已创建查询序列文件: {os.path.abspath('que.txt')}")

if __name__ == "__main__":
    create_test_reference()
    create_test_query()
    print("测试文件创建完成。现在您可以运行:")
    print("python -m sv_aligner.main ref.txt que.txt")
