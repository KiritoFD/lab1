from setuptools import setup, find_packages

setup(
    name="sv_aligner",
    version="0.1.0",
    description="SV-Aware Sequence Aligner",
    author="",
    author_email="",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'sv_aligner=sv_aligner.main:main',
        ],
    },
    python_requires='>=3.6',
    install_requires=[
        'pytest',
    ],
)
