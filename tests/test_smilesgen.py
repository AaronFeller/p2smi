from rdkit import Chem
import pytest

from p2smi.utilities.smilesgen import (
    BondSpecError,
    can_scctbond,
    can_scntbond,
    constrained_peptide_smiles,
    get_constraint_type,
    linear_peptide_smiles,
    write_library,
    write_molecule,
)
from p2smi.registry import Residue, ResidueRegistry


def _assert_valid_smiles(smiles):
    assert smiles
    assert "*" not in smiles
    assert Chem.MolFromSmiles(smiles) is not None


def _count_free_carboxylic_acids(smiles):
    mol = Chem.MolFromSmiles(smiles)
    pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    return len(mol.GetSubstructMatches(pattern))


def test_linear_peptide_smiles_uses_valid_reaction_product():
    smiles = linear_peptide_smiles("GA")
    _assert_valid_smiles(smiles)


def test_head_to_tail_constraint_returns_valid_smiles():
    _, pattern, smiles = constrained_peptide_smiles("GA", "HT")
    assert pattern == "HT"
    _assert_valid_smiles(smiles)


def test_sidechain_constraints_handle_terminal_and_carbonyl_placeholders():
    cases = [
        ("AKAAA", "SCXNXXX"),
        ("AAEAA", "SCXXZXX"),
        ("AKEAA", "SCXNZXX"),
    ]

    for sequence, pattern in cases:
        _, out_pattern, smiles = constrained_peptide_smiles(sequence, pattern)
        assert out_pattern == pattern
        _assert_valid_smiles(smiles)


def test_disulphide_constraint_returns_valid_smiles():
    _, pattern, smiles = constrained_peptide_smiles("CACAC", "SSCXXXC")
    assert pattern == "SSCXXXC"
    _assert_valid_smiles(smiles)


def test_scnt_and_scct_patterns_match_documented_constraint_types():
    _, scct_pattern = can_scctbond("AKAAA")
    _, scnt_pattern = can_scntbond("AAAEA")

    assert get_constraint_type(scct_pattern) == "SCCT"
    assert get_constraint_type(scnt_pattern) == "SCNT"


def test_write_molecule_structure_creates_text_sdf(tmp_path):
    assert write_molecule("CCO", "ethanol", "", str(tmp_path), write="structure")

    out_file = tmp_path / "3D-Files" / "ethanol_linear.sdf"
    assert out_file.exists()
    assert "V2000" in out_file.read_text()


def test_write_library_text_writes_expected_records(tmp_path):
    out_file = tmp_path / "library.p2smi"
    count = write_library(
        [("GA", "HT", "CCO"), ("AA", "", "NCC(=O)O")],
        str(out_file),
        write="text",
    )

    assert count == 2
    lines = out_file.read_text().splitlines()
    assert lines == ["GA-HT: CCO", "AA-linear: NCC(=O)O"]


def test_write_library_structure_to_file_writes_sdf_blocks(tmp_path):
    out_file = tmp_path / "library.sdf"
    count = write_library(
        [("GA", "", linear_peptide_smiles("GA"))],
        str(out_file),
        write="structure",
        write_to_file=True,
    )

    assert count == 1
    content = out_file.read_text()
    assert "$$$$" in content
    assert "_linear" in content


def test_linear_peptide_smiles_accepts_comma_delimited_v2_registry_ids():
    registry = ResidueRegistry(
        [
            Residue(
                id="LAB_001",
                name="Lab residue",
                smiles="N[C@@H](CC)C(=O)O",
            ),
            Residue(id="ALA", name="Alanine", smiles="N[C@@H](C)C(=O)O"),
        ]
    )

    smiles = linear_peptide_smiles("LAB_001,ALA", registry=registry)

    _assert_valid_smiles(smiles)


def test_constrained_peptide_smiles_derives_v2_disulphide_sites_from_base_smiles():
    registry = ResidueRegistry(
        [
            Residue(id="CYS", name="Cysteine", smiles="N[C@@H](CS)C(=O)O"),
            Residue(id="ALA", name="Alanine", smiles="N[C@@H](C)C(=O)O"),
        ]
    )

    _, pattern, smiles = constrained_peptide_smiles(
        "CYS,ALA,ALA,CYS", "SSCXXC", registry=registry
    )

    assert pattern == "SSCXXC"
    _assert_valid_smiles(smiles)
    mol = Chem.MolFromSmiles(smiles)
    assert mol.HasSubstructMatch(Chem.MolFromSmarts("[S;X2]-[S;X2]"))
    assert not any(atom.GetAtomMapNum() for atom in mol.GetAtoms())


def test_constrained_peptide_smiles_derives_v2_sidechain_amine_from_base_smiles():
    registry = ResidueRegistry(
        [
            Residue(id="LYS", name="Lysine", smiles="N[C@@H](CCCCN)C(=O)O"),
            Residue(id="ALA", name="Alanine", smiles="N[C@@H](C)C(=O)O"),
        ]
    )

    _, pattern, smiles = constrained_peptide_smiles(
        "LYS,ALA,ALA,ALA", "SCNXXX", registry=registry
    )

    assert pattern == "SCNXXX"
    _assert_valid_smiles(smiles)
    assert _count_free_carboxylic_acids(smiles) == 0


def test_constrained_peptide_smiles_derives_v2_sidechain_acid_from_base_smiles():
    registry = ResidueRegistry(
        [
            Residue(id="GLU", name="Glutamic acid", smiles="N[C@@H](CCC(=O)O)C(=O)O"),
            Residue(id="ALA", name="Alanine", smiles="N[C@@H](C)C(=O)O"),
        ]
    )

    _, pattern, smiles = constrained_peptide_smiles(
        "ALA,ALA,GLU,ALA", "SCXXZX", registry=registry
    )

    assert pattern == "SCXXZX"
    _assert_valid_smiles(smiles)
    assert _count_free_carboxylic_acids(smiles) == 1


def test_constrained_peptide_smiles_rejects_ambiguous_v2_reactive_sites():
    registry = ResidueRegistry(
        [
            Residue(id="DIOH", name="Diol", smiles="N[C@@H](C(O)CO)C(=O)O"),
            Residue(id="ALA", name="Alanine", smiles="N[C@@H](C)C(=O)O"),
        ]
    )

    with pytest.raises(BondSpecError, match="unambiguous ester reaction sites"):
        constrained_peptide_smiles("DIOH,ALA", "SCEX", registry=registry)