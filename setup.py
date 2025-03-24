from setuptools import setup
from Cython.Build import cythonize
from setuptools import Extension

ext_modules = [
    Extension(
        "core",
        sources=["core.pyx", "core_cuda.cu"],
        libraries=["cudart"],
        include_dirs=["/usr/local/cuda/include"],  # 根据CUDA安装路径修改
        library_dirs=["/usr/local/cuda/lib64"],   # 根据CUDA安装路径修改
        extra_compile_args=["-std=c++11", "-Xcompiler", "-fPIC", "-lcudart", "-arch=sm_35"], # 根据GPU架构修改
        extra_link_args=["-lcudart"],
    )
]

setup(
    name="core",
    ext_modules=cythonize(ext_modules),
)
