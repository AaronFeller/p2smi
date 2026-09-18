"""Import fatty-acid structures as peptide acylation modifiers."""

import argparse
import re
from pathlib import Path

from rdkit import Chem

from p2smi.modifiers import ModificationError, ModifierDefinition, ModifierRegistry

FREE_ACID_PATTERN = Chem.MolFromSmarts("[CX3:1](=[OX1])[OX2H1:2]")


def fatty_acid_to_acyl_fragment(mol):
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ModificationError(
            "Lipid structure must contain exactly one connected component"
        )
    matches = mol.GetSubstructMatches(FREE_ACID_PATTERN)
    if len(matches) == 1:
        selected_match = matches[0]
    elif len(matches) == 2:
        symmetry_classes = Chem.CanonicalRankAtoms(
            mol,
            breakTies=False,
            includeChirality=True,
        )
        acid_classes = {symmetry_classes[match[0]] for match in matches}
        if len(acid_classes) != 1:
            raise ModificationError(
                "Lipid diacid termini are not symmetry-equivalent; "
                "the attachment site is ambiguous"
            )
        selected_match = matches[0]
    else:
        raise ModificationError(
            "Lipid must contain one free carboxylic acid or two "
            f"symmetry-equivalent acids; found {len(matches)}"
        )
    hydroxyl_idx = selected_match[2]
    fragment = Chem.RWMol(mol)
    hydroxyl = fragment.GetAtomWithIdx(hydroxyl_idx)
    hydroxyl.SetAtomicNum(0)
    hydroxyl.SetFormalCharge(0)
    hydroxyl.SetNoImplicit(True)
    product = fragment.GetMol()
    Chem.SanitizeMol(product)
    return Chem.MolToSmiles(product, isomericSmiles=True)


def _modifier_id(value):
    identifier = re.sub(r"[^A-Za-z0-9_]+", "_", value.strip()).strip("_")
    if not identifier:
        raise ModificationError("Lipid record has no usable identifier")
    if identifier[0].isdigit():
        identifier = f"lipid_{identifier}"
    return identifier


def _load_lipid_records(input_file, id_property, name_property):
    suffix = input_file.suffix.lower()
    if suffix == ".sdf":
        supplier = Chem.SDMolSupplier(str(input_file), removeHs=True)
        records = (mol for mol in supplier if mol is not None)
    elif suffix == ".mol":
        mol = Chem.MolFromMolFile(str(input_file), removeHs=True)
        records = [] if mol is None else [mol]
    else:
        raise ModificationError("Lipid input must be an SDF or MOL file")
    for index, mol in enumerate(records, start=1):
        raw_id = mol.GetProp(id_property) if mol.HasProp(id_property) else ""
        if not raw_id and mol.HasProp("_Name"):
            raw_id = mol.GetProp("_Name")
        if not raw_id:
            raw_id = f"lipid_{index}"
        name = mol.GetProp(name_property).strip() if mol.HasProp(name_property) else ""
        if not name:
            name = raw_id
        yield _modifier_id(raw_id), name, mol


def import_lipid_modifiers(
    input_file,
    output_file,
    id_property="LM_ID",
    name_property="NAME",
    strict=False,
    max_carbons=20,
    acid_mode="both",
):
    if (
        not isinstance(max_carbons, int)
        or isinstance(max_carbons, bool)
        or max_carbons < 1
    ):
        raise ModificationError("max_carbons must be a positive integer")
    if acid_mode not in {"both", "monoacid", "diacid"}:
        raise ModificationError("acid_mode must be 'both', 'monoacid', or 'diacid'")
    input_file = Path(input_file)
    output_file = Path(output_file)
    existing = (
        ModifierRegistry.read_json(output_file)
        if output_file.exists()
        else ModifierRegistry()
    )
    imported = []
    for identifier, name, mol in _load_lipid_records(
        input_file, id_property, name_property
    ):
        try:
            carbon_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
            if carbon_count > max_carbons:
                raise ModificationError(
                    f"Lipid has {carbon_count} carbons; maximum is {max_carbons}"
                )
            acid_count = len(mol.GetSubstructMatches(FREE_ACID_PATTERN))
            if acid_mode == "monoacid" and acid_count != 1:
                raise ModificationError("Lipid is not a monoacid")
            if acid_mode == "diacid" and acid_count != 2:
                raise ModificationError("Lipid is not a diacid")
            definition = ModifierDefinition(
                id=identifier,
                name=name,
                operation="graft",
                allowed_roles=("n_terminal_amine", "sidechain_amine"),
                fragment_smiles=fatty_acid_to_acyl_fragment(mol),
            )
            definition.validate()
        except ModificationError:
            if strict:
                raise
            continue
        imported.append(definition)
    if not imported:
        raise ModificationError(f"No importable fatty acids found in {input_file}")
    combined = existing.merge(ModifierRegistry(imported))
    combined.write_json(output_file)
    return len(imported)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Import LIPID MAPS-style SDF/MOL fatty acids as acyl modifiers."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        help="LMSD-style SDF or MOL input.",
    )
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument("--id-property", default="LM_ID")
    parser.add_argument("--name-property", default="NAME")
    parser.add_argument(
        "--max-carbons",
        type=int,
        default=20,
        help="Maximum lipid carbon count (default: 20).",
    )
    parser.add_argument(
        "--acid-mode",
        choices=("both", "monoacid", "diacid"),
        default="both",
        help="Import monoacids, symmetric diacids, or both (default: both).",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail on the first structure without one unambiguous acid attachment site.",
    )
    args = parser.parse_args(argv)
    try:
        count = import_lipid_modifiers(
            args.input,
            args.output,
            id_property=args.id_property,
            name_property=args.name_property,
            strict=args.strict,
            max_carbons=args.max_carbons,
            acid_mode=args.acid_mode,
        )
    except ModificationError as exc:
        parser.error(str(exc))
    print(f"Imported {count} lipid modifiers into {args.output}")


if __name__ == "__main__":
    main()
