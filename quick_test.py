#!/usr/bin/env python3
"""
快速测试脚本：自动检查环境并运行SV-Aware比对器测试
"""

import os
import sys
import subprocess
import platform
from pathlib import Path

def check_environment():
    """检查Python环境和必要的依赖项"""
    print(f"Python版本: {platform.python_version()}")
    print(f"操作系统: {platform.system()} {platform.version()}")
    
    # 检查依赖项
    try:
        import parasail
        print("parasail库已安装")
    except ImportError:
        print("警告: parasail库未安装，将使用Python原生Smith-Waterman实现（速度较慢）")
        print("可以通过运行 'pip install parasail' 安装parasail")
    
    # 检查文件存在
    project_dir = Path(__file__).parent
    ref_txt = project_dir / "python" / "ref.txt"
    que_txt = project_dir / "python" / "que.txt"
    
    if not ref_txt.exists():
        print(f"错误: 参考文件 '{ref_txt}' 不存在")
        return False
    
    print(f"参考文件: {ref_txt} (存在)")
    
    if not que_txt.exists():
        print(f"查询文件 '{que_txt}' 不存在，将由test_alignment.py自动生成")
    else:
        print(f"查询文件: {que_txt} (存在)")
    
    return True

def main():
    """快速测试主函数"""
    print("=== SV-Aware序列比对器快速测试 ===")
    
    # 检查环境
    if not check_environment():
        return 1
    
    # 运行测试脚本
    test_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_alignment.py")
    
    print("\n运行测试脚本...")
    try:
        result = subprocess.run([sys.executable, test_script], check=True)
        if result.returncode == 0:
            print("\n✅ 测试成功完成")
            return 0
        else:
            print(f"\n❌ 测试失败，返回码: {result.returncode}")
            return result.returncode
    except subprocess.CalledProcessError as e:
        print(f"\n❌ 测试脚本执行失败: {e}")
        return e.returncode
    except Exception as e:
        print(f"\n❌ 执行测试时出错: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
