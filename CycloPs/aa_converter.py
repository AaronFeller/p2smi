#!/usr/bin/env python
"""Rearrange amino-acid SMILES into a form usable by CycloPs"""
import sys
import openbabel, pybel


def rearrange_smiles(aa_smiles):
    """Rewrite an amino-acid smiles to start with the N-term and end with
    the C-term."""
    mol = pybel.readstring("smi", aa_smiles)
    n_term_pat = pybel.Smarts("[$(NCC(O)=O)]")
    c_term_pat = pybel.Smarts("[$(OC(=O)CN)]")
    # Find location of start and end atoms
    n_term_idx = n_term_pat.findall(mol)[0][0]
    c_term_idx = c_term_pat.findall(mol)[0][0]
    # Rewrite smiles N-term first, then C-term
    rearranger = openbabel.OBConversion()
    rearranger.SetInAndOutFormats("smi", "smi")
    rearranger.AddOption("f", openbabel.OBConversion.OUTOPTIONS, str(n_term_idx))
    rearranger.AddOption("l", openbabel.OBConversion.OUTOPTIONS, str(c_term_idx))
    outmol = openbabel.OBMol()
    rearranger.ReadString(outmol, aa_smiles)
    return rearranger.WriteString(outmol).strip()


def test():
    natural_amino_acid_smiles = [
        "O=C(O)[C@@H](N)C",
        "SC[C@H](N)C(O)=O",
        "OC(C[C@H](N)C(O)=O)=O",
        "O=C(O)[C@@H](N)CCC(O)=O",
        "N[C@H](C(O)=O)CC1=CC=CC=C1",
        "O=C(O)CN",
        "N[C@@H](CC1=CN=CN1)C(O)=O",
        "O=C(O)[C@@H](N)[C@@H](C)CC",
        "N[C@H](C(O)=O)CCCCN",
        "N[C@@H](CC(C)C)C(O)=O",
        "OC([C@@H](N)CCSC)=O",
        "NC(C[C@H](N)C(O)=O)=O",
        "O=C([C@@H]1CCCN1)O",
        "OC([C@@H](N)CCC(N)=O)=O",
        "O=C(O)[C@@H](N)CCCNC(N)=N",
        "OC([C@@H](N)CO)=O",
        "N[C@H](C(O)=O)[C@H](O)C",
        "N[C@H](C(O)=O)C(C)C",
        "O=C(O)[C@@H](N)CC1=CNC2=C1C=CC=C2",
        "N[C@@H](CC1=CC=C(O)C=C1)C(O)=O",
    ]
    for aa_smiles in natural_amino_acid_smiles:
        print((aa_smiles, "-->", rearrange_smiles(aa_smiles)))


def main():
    if len(sys.argv) > 1 and sys.argv[1] != "-t" and sys.argv[1] != "-h":
        for aa_smiles in sys.argv[1:]:
            print((aa_smiles, "-->", rearrange_smiles(aa_smiles)))
    elif len(sys.argv) > 1 and sys.argv[1] == "-t":
        test()
    else:
        print(
            """Rearrange amino-acid smiles to start at the N-terminus 
        and end at the C-terminus>
        Usage: aa_converter.py [-ht] smiles1 [smiles2] ...
        Options: -h.  Print this message and exit.
                 -t.  Do test"""
        )


if __name__ == "__main__":
    main()
