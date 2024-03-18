#!/usr/bin/env python
"""
cx_Freeze script to make pepgui.py an executable
"""
from cx_Freeze import setup, Executable
import shutil
import os
import sys

package_loc = r"../../Peptide Library Package"
print(("Path exists: ", os.path.isdir(package_loc)))
# Make sure the locally stored PepLibGen is used, not the one on Pythonpath
# This is to allow old CycloPs tags in SVN use old version of PepLibGen
sys.path.insert(0, package_loc)

build_dir = "linuxdist"
include_files = [("../zinc_output", build_dir)]
bin_include = [
    "libBLT.2.4.so.8.5",
    "libtk8.5.so.0",
    "libtcl8.5.so.0",
    "liblapack.so.3gf",
    "libblas.so.3gf",
]
bin_path_include = []
bin_excludes = []
includes = ["rdkit", "numpy"]
if sys.platform == "darwin":
    bin_include = []

    #    bin_excludes = ["libicuuc.50.dylib", "libicudata.50.dylib",
    #                 "libDepictor.1.dylib", "libFileParsers.1.dylib",
    #                  "libSubstructMatch.1.dylib", "libChemTransforms.1.dylib",
    #                  "libSmilesParse.1.dylib", "libGraphMol.1.dylib",
    #                  "libRDGeometryLib.1.dylib", "libDataStructs.1.dylib",
    #                  "libRDGeneral.1.dylib", "libForceFieldHelpers.1.dylib",
    #                  "libDistGeometry.1.dylib", "libAlignment.1.dylib",
    #                  "libForceField.1.dylib", "libOptimizer.1.dylib",
    #                  "libEigenSolvers.1.dylib", "libPartialCharges.1.dylib",
    #                  "libSubgraphs.1.dylib"]
    bin_path_includes = ["/usr/local/Cellar/rdkit/2012.12.1/lib/"]


setup(
    name="pepgui",
    version="1b",
    description="test",
    options={
        "build_exe": {
            "includes": includes,
            "excludes": ["matplotlib", "scipy"],
            "copy_dependent_files": True,
            "bin_includes": bin_include,
            "include_files": include_files,
            "bin_excludes": bin_excludes,
            "build_exe": build_dir,
            "bin_path_excludes": ["@loader_path"],
            "bin_path_includes": bin_path_includes,
        }
    },
    executables=[Executable("../CycloPs.py", base=None)],
)

# HACK HACK HACK
# If files are not found by cx_Freeze, add them in.
for lib in bin_include:
    if lib not in os.listdir(build_dir):
        shutil.copy(os.path.join("/usr/lib", lib), build_dir)
