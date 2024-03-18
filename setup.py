#!/usr/bin/env python

from distutils.core import setup


setup(
    name="PepLibGen",
    version="1.2",
    url="https://github.com/fergaljd/cyclops",
    author="Fergal Duffy",
    author_email="fergaljd@gmail.com",
    packages=["PepLibGen", "PepLibGen.StructGen", "PepLibGen.Analysis"],
    scripts=["CycloPs/CycloPs.py", "CycloPs/aa_converter.py"],
)
