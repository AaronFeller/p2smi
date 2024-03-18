#!/usr/bin/env python
"""Process text files of non-natural aminos from Zinc"""

from rdkit import Chem


def extract_peptides(zinc_file):
    """Extract peptide names and SMILES from output zinc_file"""
    handle = open(zinc_file, "r")
    pep_dict = {}
    with open(zinc_file) as handle:
        for count, line in enumerate(handle):
            smiles, name = line.strip().split()
            # Remove Fmoc placeholder
            if smiles.startswith("[*]"):
                smiles = smiles.replace("[*]", "", 1)
            # Deionise C and N-term if necessary
            if smiles.endswith("[O-]"):
                smiles = smiles.rpartition("[O-]")[0] + "O"
            for pos_ion in ["[NH3+]", "[NH2+]", "[NH+]"]:
                if smiles.startswith(pos_ion):
                    smiles = smiles.replace(pos_ion, "N", 1)
            # If SMILES is parseable, add to peptide dictionary
            if Chem.MolFromSmiles(smiles):
                pep_dict[name] = {
                    "SMILES": smiles,
                    "Code": "",
                    "Formula": "",
                    "Letter": "",
                    "MolWeight": "",
                    "cterm": False,
                    "disulphide": False,
                    "ester": False,
                    "nterm": False,
                }
            else:
                print(("%s not valid, line %d" % (smiles, count + 1)))
    return pep_dict


if __name__ == "__main__":
    import sys

    zinc_file = sys.argv[1]
    outfile = zinc_file.split(".")[0] + "_extracted"
    with open(outfile, "w") as out:
        peps = extract_peptides(zinc_file)
        for peptide in peps:
            out.write(peps[peptide]["SMILES"])
            out.write("\n")
