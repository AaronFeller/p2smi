#!/usr/bin/env python
"""
Convert an amino-acid smiles into CycloPs format, with the N-terminal at the end
ad the C-terminal at the start.
"""
import pybel


def convert_smiles(smiles):
    # Find the locations of N and C-termini in the amino acid
    aa = pybel.readstring("smi", smiles)
    n_term_patt = pybel.Smarts("[$([NH2]CC=O)]")
    c_term_patt = pybel.Smarts("[$([OH]C(C)=O)]")
    n_term_idx = n_term_patt.findall(aa)[0][0]
    c_term_idx = c_term_patt.findall(aa)[0][0]
    mol_size = len(aa.atoms)
    # Make external call to babel to move N and C termini
    # obabel -:'CC(C(=O)O)N' -osmi -xl 4 -xf 6
