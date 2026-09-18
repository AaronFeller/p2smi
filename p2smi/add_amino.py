"""CLI for adding validated noncanonical amino acids to a JSON registry."""

import argparse
from pathlib import Path

from rdkit import Chem

from p2smi.registry import RegistryError, Residue, ResidueRegistry


class AminoLibraryError(Exception):
    pass


def _prompt(value, message):
    if value is not None:
        return value
    return input(message).strip()


def _validate_smiles(smiles, label, *, require_dummy=False):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise AminoLibraryError(f"Invalid {label} SMILES: {smiles}")

    dummy_count = sum(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms())
    if require_dummy and dummy_count != 1:
        raise AminoLibraryError(
            f"{label} SMILES must contain exactly one '*' placeholder: {smiles}"
        )
    if not require_dummy and dummy_count != 0:
        raise AminoLibraryError(
            f"Base residue SMILES cannot contain '*' placeholders: {smiles}"
        )
    return mol


def smiles_from_mol_file(mol_file):
    try:
        mol = Chem.MolFromMolFile(str(mol_file), sanitize=True, removeHs=True)
    except OSError as exc:
        raise AminoLibraryError(f"Could not read MOL file: {mol_file}") from exc
    if mol is None:
        raise AminoLibraryError(f"Invalid MOL file: {mol_file}")
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def build_registry_residue(identifier, name, code, letter, smiles):
    mol = _validate_smiles(smiles, "base residue")

    residue = Residue(
        id=identifier,
        name=name,
        code=code or None,
        legacy_letter=letter or None,
        smiles=Chem.MolToSmiles(mol, isomericSmiles=True),
    )
    try:
        ResidueRegistry.validate_residue(residue)
    except RegistryError as exc:
        raise AminoLibraryError(str(exc)) from exc
    return residue


def append_entry_to_registry(registry_file, residue):
    registry = (
        ResidueRegistry.read_json(registry_file)
        if registry_file.exists()
        else ResidueRegistry()
    )
    registry.add(residue)
    registry.write_json(registry_file)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Add a noncanonical amino acid to a v2 JSON residue registry.",
    )
    parser.add_argument(
        "--id", dest="residue_id", help="Stable ASCII v2 residue identifier."
    )
    parser.add_argument("--name", help="Full residue name.")
    parser.add_argument("--code", help="Optional unique residue code alias.")
    parser.add_argument("--letter", help="Optional one-character legacy alias.")
    structure_source = parser.add_mutually_exclusive_group()
    structure_source.add_argument(
        "--smiles", help="Base amino-acid SMILES without '*' placeholders."
    )
    structure_source.add_argument(
        "--mol-file",
        type=Path,
        help=(
            "MDL MOL file for the base residue; its canonical isomeric "
            "SMILES is stored."
        ),
    )
    parser.add_argument(
        "--registry-file",
        type=Path,
        required=True,
        help="JSON residue registry to create or extend.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    name = _prompt(args.name, "Residue name: ")
    smiles = (
        smiles_from_mol_file(args.mol_file)
        if args.mol_file
        else _prompt(args.smiles, "Base residue SMILES: ")
    )

    identifier = _prompt(args.residue_id, "Stable residue ID: ")
    try:
        residue = build_registry_residue(
            identifier,
            name,
            args.code,
            args.letter,
            smiles,
        )
        append_entry_to_registry(args.registry_file, residue)
    except RegistryError as exc:
        raise AminoLibraryError(str(exc)) from exc
    print(f"Added {residue.id} to {args.registry_file}")


if __name__ == "__main__":
    main()
