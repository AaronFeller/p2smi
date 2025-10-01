#!/usr/bin/env python
"""
Chemical data extraction for a molecule (fast path).

Improvements:
- Parse SMILES once per record; reuse a single RDKit Mol.
- No recursive calls that re-parse (e.g., Lipinski pass computed from the same Mol).
- Streamed batch I/O with minimal per-iteration overhead.
- Optional multiprocessing (--procs) for large inputs.
"""

import argparse
import json
from typing import Tuple, Optional

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import (
    Crippen,  # logP
    Descriptors,  # MW, TPSA
    Lipinski,  # donors/acceptors, rotatable bonds
    rdMolDescriptors,  # formula, frac Csp3
    rdmolops,  # formal charge
)

RDLogger.DisableLog("rdApp.*")  # quieter batch runs


class SmilesError(Exception):
    pass


# ---------- Core helpers (Mol-first; no re-parsing) ----------


def make_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SmilesError(f"{smiles} is not a valid SMILES string")
    return mol


def lipinski_trial_mol(mol: Chem.Mol) -> Tuple[list, list]:
    passed, failed = [], []
    hdon = Lipinski.NumHDonors(mol)
    hacc = Lipinski.NumHAcceptors(mol)
    mw = Descriptors.MolWt(mol)
    clogp = Crippen.MolLogP(mol)

    (failed if hdon > 5 else passed).append(
        f"{'Over ' if hdon>5 else ''}5 H-bond donors"
        + (f" (found {hdon})" if hdon > 5 else f" ({hdon})")
    )
    (failed if hacc > 10 else passed).append(
        f"{'Over ' if hacc>10 else ''}10 H-bond acceptors"
        + (f" (found {hacc})" if hacc > 10 else f" ({hacc})")
    )
    (failed if mw >= 500 else passed).append(
        f"{'Molecular weight over 500' if mw>=500 else 'Molecular weight'} (calculated {mw:.2f})"
    )
    (failed if clogp >= 5 else passed).append(
        f"{'logP over 5' if clogp>=5 else 'logP'} (calculated {clogp:.2f})"
    )
    return passed, failed


def molecule_summary_from_mol(smiles: str, mol: Chem.Mol) -> dict:
    # Bind locals (tiny speedup in tight loops)
    _MolWt = Descriptors.MolWt
    _TPSA = Descriptors.TPSA
    _MolLogP = Crippen.MolLogP
    _HDon = Lipinski.NumHDonors
    _HAcc = Lipinski.NumHAcceptors
    _RotB = Lipinski.NumRotatableBonds
    _Rings = mol.GetRingInfo().NumRings
    _FracC = rdMolDescriptors.CalcFractionCSP3
    _Heavy = mol.GetNumHeavyAtoms
    _Charge = rdmolops.GetFormalCharge
    _Formula = rdMolDescriptors.CalcMolFormula

    passed, failed = lipinski_trial_mol(mol)

    return {
        "SMILES": smiles,
        "Formula": _Formula(mol),
        "Molecular weight": round(_MolWt(mol), 2),
        "logP": round(_MolLogP(mol), 2),
        "TPSA": round(_TPSA(mol), 2),
        "H-bond donors": _HDon(mol),
        "H-bond acceptors": _HAcc(mol),
        "Rotatable bonds": _RotB(mol),
        "Rings": _Rings(),
        "Fraction Csp3": round(_FracC(mol), 3),
        "Heavy atoms": _Heavy(),
        "Formal charge": _Charge(mol),
        "Lipinski pass": not failed,
    }


# Backward-compatible single-call helper
def molecule_summary(smiles: str) -> dict:
    return molecule_summary_from_mol(smiles, make_mol(smiles))


# ---------- Batch processing ----------


def parse_smiles_line(line: str) -> Optional[str]:
    s = line.strip()
    if not s:
        return None
    # support lines like "id: SMILES" or just "SMILES"
    return s.split(": ", 1)[-1] if ": " in s else s


def process_line(line: str) -> str:
    s = parse_smiles_line(line)
    if s is None:
        return ""  # skip empties
    try:
        mol = make_mol(s)
        res = molecule_summary_from_mol(s, mol)
        return json.dumps(res)
    except SmilesError as e:
        return json.dumps({"error": f"{e}", "SMILES": s})


# ---------- CLI ----------


def main():
    ap = argparse.ArgumentParser(
        description="Analyze SMILES strings for molecular properties."
    )
    ap.add_argument("-s", "--smiles", help="Single SMILES; prints JSON to stdout.")
    ap.add_argument(
        "-i", "--input_file", help="Text file: one SMILES (or 'id: SMILES') per line."
    )
    ap.add_argument(
        "-o", "--output_file", help="If given with -i, writes JSONL here; else prints."
    )
    ap.add_argument(
        "--procs",
        type=int,
        default=0,
        help="Use N processes for batch mode (0=disable).",
    )
    args = ap.parse_args()

    if not args.smiles and not args.input_file:
        ap.error("At least one of --smiles or --input_file must be provided.")

    # Single SMILES path
    if args.smiles:
        res = molecule_summary(args.smiles)
        print(json.dumps(res, indent=2))
        return

    # Batch path
    if args.input_file:
        if args.procs and args.procs > 1:
            # Multiprocessing for large files
            from multiprocessing import Pool

            with open(args.input_file, "r") as inf:
                if args.output_file:
                    with (
                        open(args.output_file, "w") as outf,
                        Pool(processes=args.procs) as pool,
                    ):
                        for out in pool.imap_unordered(
                            process_line, inf, chunksize=1024
                        ):
                            if out:
                                outf.write(out + "\n")
                else:
                    with Pool(processes=args.procs) as pool:
                        for out in pool.imap_unordered(
                            process_line, inf, chunksize=1024
                        ):
                            if out:
                                print(out)
        else:
            # Single-process fast stream
            if args.output_file:
                with (
                    open(args.input_file, "r") as inf,
                    open(args.output_file, "w") as outf,
                ):
                    for line in inf:
                        out = process_line(line)
                        if out:
                            outf.write(out + "\n")
            else:
                with open(args.input_file, "r") as inf:
                    for line in inf:
                        out = process_line(line)
                        if out:
                            print(out)


if __name__ == "__main__":
    main()
