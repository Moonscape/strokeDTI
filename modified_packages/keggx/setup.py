from setuptools import setup, find_packages

setup(
    name="keggx",
    packages=find_packages(),
    package_dir={"keggx": "keggx"},
    package_data={"keggx": ["KEGG_compound_ids.txt"]},
    include_package_data=True,
    version="0.1.0",
    url="https://github.com/iamjli/keggx",
    python_requires=">=3.5",
    classifiers=[
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    license="MIT",
    author="iamjli",
    author_email="iamjli@mit.edu",
    description="Python package for manipulation and visualization of KEGG pathways. [Modified because of compatibility issues]",
    install_requires=[
        "requests",
    ],
)
