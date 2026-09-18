import pytest
from rdkit import Chem

from p2smi.add_amino import (
    AminoLibraryError,
    build_registry_residue,
    main,
    smiles_from_mol_file,
)
from p2smi.registry import ResidueRegistry


def test_registry_mode_creates_growing_json_corpus(tmp_path, capsys):
    registry_file = tmp_path / "lab-residues.json"

    main(
        [
            "--registry-file",
            str(registry_file),
            "--id",
            "LAB_001",
            "--name",
            "Lab residue",
            "--smiles",
            "N[C@@H](CC)C(=O)O",
        ]
    )

    main(
        [
            "--registry-file",
            str(registry_file),
            "--id",
            "LAB_002",
            "--name",
            "Second lab residue",
            "--smiles",
            "N[C@@H](CCC)C(=O)O",
        ]
    )

    registry = ResidueRegistry.read_json(registry_file)
    assert [residue.id for residue in registry] == ["LAB_001", "LAB_002"]
    assert "Added LAB_002" in capsys.readouterr().out


def test_registry_entry_rejects_invalid_v2_identifier():
    with pytest.raises(AminoLibraryError):
        build_registry_residue(
            "LAB 001",
            "Lab residue",
            None,
            None,
            "N[C@@H](CC)C(=O)O",
        )


def test_mol_file_input_creates_v2_residue_without_smiles_argument(tmp_path):
    mol_file = tmp_path / "lab-residue.mol"
    source_smiles = "N[C@@H](CC)C(=O)O"
    Chem.MolToMolFile(Chem.MolFromSmiles(source_smiles), str(mol_file))
    registry_file = tmp_path / "lab-residues.json"

    main(
        [
            "--registry-file",
            str(registry_file),
            "--id",
            "LAB_MOL_001",
            "--name",
            "MOL file residue",
            "--mol-file",
            str(mol_file),
        ]
    )

    residue = ResidueRegistry.read_json(registry_file).resolve("LAB_MOL_001")
    assert residue.smiles == smiles_from_mol_file(mol_file)
    assert "constraints" not in registry_file.read_text()
    assert Chem.MolToSmiles(
        Chem.MolFromSmiles(residue.smiles), isomericSmiles=True
    ) == (Chem.MolToSmiles(Chem.MolFromSmiles(source_smiles), isomericSmiles=True))
