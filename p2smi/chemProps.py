#!/usr/bin/env python
"""
Chemical data extraction for a molecule.

Molecules are defined by SMILES strings.
This script calculates:
- logP (partition coefficient)
- Lipinski's rule evaluations (with pass/fail reasons)
- Molecular formula
- Molecular weight
- Topological Polar Surface Area (TPSA)
- Number of H-bond donors and acceptors
- Rotatable bonds count
- Ring count
- Fraction of sp3 carbons (Csp3)
- Heavy atom count
- Formal charge
- Overall Lipinski rule pass/fail

Uses RDKit for all chemical computations.
"""

from rdkit import Chem
from rdkit.Chem import (
    Crippen,  # logP calculation
    Descriptors,  # molecular weight, TPSA, etc.
    Lipinski,  # H-bond donors/acceptors, rotatable bonds
    rdMolDescriptors,  # molecular formula, fraction Csp3
    rdmolops,  # formal charge
)


class SmilesError(Exception):
    # Custom exception for invalid SMILES strings
    pass


def log_partition_coefficient(smiles):
    # Calculate logP for a given SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SmilesError(f"{smiles} is not a valid SMILES string")
    return Crippen.MolLogP(mol)


def lipinski_trial(smiles):
    # Evaluate Lipinski's rules, returning lists of passes and failures
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SmilesError(f"{smiles} is not a valid SMILES string")

    passed, failed = [], []
    num_hdonors = Lipinski.NumHDonors(mol)
    num_hacceptors = Lipinski.NumHAcceptors(mol)
    mol_weight = Descriptors.MolWt(mol)
    mol_logp = Crippen.MolLogP(mol)

    if num_hdonors > 5:
        failed.append(f"Over 5 H-bond donors (found {num_hdonors})")
    else:
        passed.append(f"{num_hdonors} H-bond donors")

    if num_hacceptors > 10:
        failed.append(f"Over 10 H-bond acceptors (found {num_hacceptors})")
    else:
        passed.append(f"{num_hacceptors} H-bond acceptors")

    if mol_weight >= 500:
        failed.append(f"Molecular weight over 500 (calculated {mol_weight:.2f})")
    else:
        passed.append(f"Molecular weight: {mol_weight:.2f}")

    if mol_logp >= 5:
        failed.append(f"logP over 5 (calculated {mol_logp:.2f})")
    else:
        passed.append(f"logP: {mol_logp:.2f}")

    return passed, failed


def lipinski_pass(smiles):
    # Return True if Lipinski's rules are passed, otherwise False
    _, failed = lipinski_trial(smiles)
    return not failed


def molecular_formula(smiles):
    # Return molecular formula from a SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SmilesError(f"{smiles} is not a valid SMILES string")
    return rdMolDescriptors.CalcMolFormula(mol)


def tpsa(smiles):
    # Calculate topological polar surface area (TPSA)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SmilesError(f"{smiles} is not a valid SMILES string")
    return Descriptors.TPSA(mol)


def molecule_summary(smiles):
    # Return a dictionary with all computed molecular properties
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SmilesError(f"{smiles} is not a valid SMILES string")

    summary = {
        "Formula": molecular_formula(smiles),
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "logP": round(Crippen.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "H-bond donors": Lipinski.NumHDonors(mol),
        "H-bond acceptors": Lipinski.NumHAcceptors(mol),
        "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),
        "Rings": mol.GetRingInfo().NumRings(),
        "Fraction Csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        "Heavy Atoms": mol.GetNumHeavyAtoms(),
        "Formal Charge": rdmolops.GetFormalCharge(mol),
        "Lipinski pass": lipinski_pass(smiles),
    }
    return summary


if __name__ == "__main__":
    import argparse
    import json

    # Command-line interface for computing and printing molecule summary
    parser = argparse.ArgumentParser(
        description="Analyze SMILES strings for molecular properties."
    )
    parser.add_argument("smiles", type=str, help="SMILES string of the molecule")
    args = parser.parse_args()

    result = molecule_summary(args.smiles)
    print(json.dumps(result, indent=2))
