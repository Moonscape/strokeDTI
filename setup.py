from setuptools import setup, find_packages

setup(
    name="strokeDTI",
    version="0.1.0",
    description="Implementation of strokeDTI",
    author="PJJ",
    author_email="ninelightyearsaway@gmail.com",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "matplotlib",
        "numpy",
        "scipy",
        "bio",
        "reportlab",
        "networkx",
        "kaleido",
        "nbformat>=4.2.0",
    ],
    entry_points={
        "console_scripts": [
            "identify_targets=strokeDTI.scripts.identify_targets_cli:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)
