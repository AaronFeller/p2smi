#!/usr/bin/env python
"""
Functions to determine whether a natural peptide can be chemically synthesised.

Functions include tests for forbidden motifs, required number of charges, 
hydrophobicity.
"""
import sys
import re

from rdkit import Chem
from rdkit.Chem import Crippen

# Hydrophobicity library
# from  Bio.SeqUtils.ProtParam import hydrophobicity
# Don't use kyte and doolittle, use eisenberg
# Author(s): Eisenberg D., Schwarz E., Komarony M., Wall R.
# Reference: J. Mol. Biol. 179:125-142(1984).

hydrophobicity = {
    "A": 0.620,
    "R": -2.530,
    "N": -0.780,
    "D": -0.900,
    "C": 0.290,
    "Q": -0.850,
    "E": -0.740,
    "G": 0.480,
    "H": -0.400,
    "I": 1.380,
    "L": 1.060,
    "K": -1.500,
    "M": 0.640,
    "F": 1.190,
    "P": 0.120,
    "S": -0.180,
    "T": -0.050,
    "W": 0.810,
    "Y": 0.260,
    "V": 1.080,
}


# Motifs that prevent chemical synthesis.
# These will be described as regular expressions
forbidden_motifs = {
    "Over 2 prolines in a row are difficult to synthesise": r"[P]{3,}",
    "DG and DP are difficult to synthesis": r"D[GP]",
    "N or Q  at N-terminus are difficult to synthesise": r"^[NQ]",
}

# When residues have side chains that are bonded to something,
# their hydrophobicity changes
# To work around this, swap for a residue with similar properties to the bonded
# side-chain
hydrophobicity_swaps = {
    "S": "G",
    "D": "N",
    "E": "Q",
    "C": "M",
    "K": "L",
    "Y": "F",
    "T": "V",
}

# Charged residues
charged = ["H", "R", "K", "E", "D"]


class SmilesError(Exception):
    pass


def log_partition_coefficient(smiles):
    """
    Returns the octanol-water partition coefficient given a molecule SMILES
    string
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception as e:
        raise SmilesError("Could not parse SMILES: %s" % smiles)
    return Crippen.MolLogP(mol)


def test_motifs(seq, motifs=forbidden_motifs):
    """
    Determines whether a peptide sequence contains motifs that prevent synthesis

    Returns a list of forbidden motifs found
    """
    # seq = seq.upper() #Allow equivalent D-amino acids
    return [test for test in list(motifs.keys()) if re.search(motifs[test], seq)]


def mod_seq(seq, bond_def):
    """
    Replaces bonded residues in a peptide with non-bonded residues with similar
    properties
    """
    # seq = seq.upper() #Allow equivalent D-amino acids
    mod_seq = ""
    for resi, is_bonded in zip(seq, bond_def):
        if is_bonded == "X":
            mod_seq += resi
        else:
            mod_seq += hydrophobicity_swaps[resi]
    return mod_seq


def test_hydrophobicity(seq):
    """
    Returns peptide hydrophobicity in the form of logP(o/w).
    Positive values are hydrophobic, negative hydrophilic.
    Will get better values using log_partition_coefficient from analysis package
    """
    seq = seq.upper()  # Allow equivalent D-amino acids
    return sum([hydrophobicity[resi] for resi in seq])


def test_charge(seq):
    """
    Returns True if peptide has a charged residue at least every 5 residues
    """
    # seq = seq.upper() #Allow equivalent D-amino acids
    count = 0
    for resi in seq:
        count += 1
        if resi in charged:
            count = 0
        if count >= 5:
            return False
    return True


def run_all_tests(smiles, seq, bond_def=""):
    """
    Run each test function on a seq. Return which ones failed.

    Returns True if passed
    """
    seq = seq.upper()  # Allow equivalent D-amino acids
    fail_strings = []
    if not bond_def:
        bond_def = "X" * len(seq)
    fail_strings.extend(["Failed motifs: " + message for message in test_motifs(seq)])
    hydrophobicity = log_partition_coefficient(smiles)
    hydro_msg = "Failed hydrophobicity: logP %s" % hydrophobicity
    if hydrophobicity > 0:
        fail_strings.append(hydro_msg)

    if not test_charge(mod_seq(seq, bond_def)):
        fail_strings.append("Failed charge: need 1 charged resi every 5 resis")

    if fail_strings:
        return fail_strings
    else:
        return True


def synth_pass(args, strict=False):
    """
    Wraps around run_all_tests but returns True or False only
    """
    seq, bond_def, smiles = args
    if not seq and not bond_def and not smiles:
        return False
    if bond_def[:2] in ["SS", "HT", "SC"]:
        bond_def = bond_def[2:]
    try:
        if run_all_tests(smiles, seq, bond_def) == True:
            return True
        else:
            return False
    except AttributeError:  # Tuple of names, not string
        # try:
        from p2smi.StructGen import aminoacids as aa

        try:
            seq = "".join([aa.aminos[resi]["Letter"] for resi in seq])
        except KeyError:
            # Unknown amino
            if strict:
                return False
            else:
                return True
        if run_all_tests(smiles, seq, bond_def) == True:
            return True
        else:
            return False


# test...
if __name__ == "__main__":
    test = [("ASDF", ""), ("KKSD", "NXXZ"), ("SDFG", "EXXX")]
    for seq, bond_def in test:
        print((run_all_tests(seq, bond_def)))
