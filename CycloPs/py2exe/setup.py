#!/usr/bin/env python
"""
Py2exe script to make pepgui.py a win32 executable
"""
from distutils.core import setup
import py2exe
import glob
import os
import sys
import PepLibGen

package_loc = r"../../Peptide Library Package"
print(("Path exists: ", os.path.isdir(package_loc)))
# Make sure the locally stored PepLibGen is used, not the one on Pythonpath
# This is to allow old CycloPs tags in SVN use old version of PepLibGen
sys.path.insert(0, package_loc)

data_files = [
    (".", [os.path.join(os.path.dirname(os.getcwd()), "zinc_output")]),
    (".", [os.path.join(os.path.dirname(os.getcwd()), "example_peptide_file.txt")]),
]

opts = {"py2exe": {"dll_excludes": ["MSVCR80.dll", "MSVCP80.dll"]}}

print("Imported PepLibGen is...")
print((PepLibGen.__file__))

setup(
    options=opts,
    author="Fergal Duffy",
    author_email="fergaljd@gmail.com",
    windows=["../CycloPs.py"],
    data_files=data_files,
)
