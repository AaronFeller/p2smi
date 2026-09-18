import pytest
from rdkit import Chem

from p2smi.chemMods import (
    add_n_methylation,
    add_pegylation,
    is_valid_smiles,
    modify_sequence,
    parse_input_lines,
    process_recipe_sequences,
    process_sequences,
)
from p2smi.modifiers import (
    ModificationError,
    ModifierRegistry,
    builtin_modifier_registry,
)


def test_is_valid_smiles():
    assert is_valid_smiles("C1CCCCC1")  # cyclohexane
    assert not is_valid_smiles("INVALID_SMILES")


def test_add_n_methylation_counts():
    test_seq = "NCC(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O"  # 2 peptide amides
    modified_seq, count = add_n_methylation(test_seq, 1.0)
    assert count == 2
    assert Chem.MolFromSmiles(modified_seq) is not None


def test_add_n_methylation_handles_canonicalized_smiles():
    test_seq = Chem.MolToSmiles(
        Chem.MolFromSmiles("NCC(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O"),
        isomericSmiles=True,
    )
    modified_seq, count = add_n_methylation(test_seq, 1.0)
    assert count == 2
    assert Chem.MolFromSmiles(modified_seq) is not None


def test_add_pegylation_adds_peg():
    test_seq = "N[C@@H](CCCCN)C(=O)O"  # L-Lysine
    peg_seq, peg = add_pegylation(test_seq)
    assert peg_seq != test_seq
    assert peg and "CCO" in peg
    assert Chem.MolFromSmiles(peg_seq) is not None


def test_add_pegylation_handles_canonicalized_smiles():
    test_seq = Chem.MolToSmiles(
        Chem.MolFromSmiles("N[C@@H](CCCCN)C(=O)O"), isomericSmiles=True
    )
    peg_seq, peg = add_pegylation(test_seq)
    assert peg and "CCO" in peg
    assert peg_seq != test_seq
    assert Chem.MolFromSmiles(peg_seq) is not None


def test_parse_input_lines_handles_good_and_bad_lines():
    lines = [
        "peptide1: C1CCCCC1",
        "malformed line without colon",
        "peptide2: C1CCNCC1",
    ]
    parsed = list(parse_input_lines(lines))
    assert parsed[0] == ("peptide1", "C1CCCCC1")
    assert parsed[1] == ("[Malformed line]", None)
    assert parsed[2] == ("peptide2", "C1CCNCC1")


def test_modify_sequence_modifications_included():
    test_seq = "C(=O)N[C@H]C(=O)N[C@H]C"
    mod_seq, mods = modify_sequence(
        test_seq, do_methylate=True, do_pegylate=True, nmeth_residues=0.5
    )
    assert any("N-methylation" in m for m in mods)
    assert any("PEGylation" in m for m in mods)


def test_process_sequences_runs_and_marks_modifications():
    input_lines = [
        "pep1: C(=O)N[C@H]C(=O)N[C@H]C",
        "pep2: N[C@@H](CCCCN)C(=O)O",
        "bad_line",
    ]
    results = list(
        process_sequences(input_lines, nmeth_rate=1.0, peg_rate=1.0, nmeth_residues=0.5)
    )
    assert any("pep1" in r and "N-methylation" in r for r in results)
    assert any("pep2" in r and "PEGylation" in r for r in results)
    assert any("Skipped malformed line" in r for r in results)


def test_process_recipe_sequences_applies_recipe_and_can_skip_errors(tmp_path):
    registry_file = tmp_path / "modifiers.json"
    registry_file.write_text(
        '{"schema_version":1,"kind":"p2smi_modifier_registry",'
        '"modifiers":[],"recipes":[{"id":"cap","steps":['
        '{"modifier":"n_acetyl","site":{"terminus":"N"}}]}]}'
    )
    registry = ModifierRegistry.read_json(
        registry_file, base=builtin_modifier_registry()
    )

    output = list(process_recipe_sequences(["pep: NCC(=O)O"], "cap", registry))
    assert output[0].startswith("pep[n_acetyl]: ")
    assert Chem.MolFromSmiles(output[0].split(": ", 1)[1]) is not None

    skipped = list(process_recipe_sequences(["bad"], "cap", registry, on_error="skip"))
    assert "Skipped malformed line" in skipped[0]
    with pytest.raises(ModificationError):
        list(process_recipe_sequences(["bad"], "cap", registry))
