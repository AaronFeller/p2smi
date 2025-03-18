from setuptools import setup, find_packages

setup(
    name="p2smi",  # Replace with your package name
    version="0.1.0",  # Version of your package
    author="Aaron Feller",
    author_email="aaronleefeller@gmail.com",
    license="MIT",
    description="A toolkit for conversion of unnatural and natural single-letter amino acid strings (peptides/proteins) to SMILES strings",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/AaronFeller/p2smi",  # Repo URL
    packages=find_packages(),  # Automatically finds all packages in your project
    install_requires=[  # List your dependencies
        "rdkit",
    ],
    classifiers=[  # Metadata about your package
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
    ],
    python_requires=">=3.9",  # Minimum Python version
)