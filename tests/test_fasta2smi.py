import pytest
from rdkit import Chem

# Mock smilesgen
import p2smi.utilities.smilesgen as smilesgen
from p2smi.fasta2smi import (
    InvalidConstraintError,
    constraint_resolver,
    generate_smiles_strings,
    parse_fasta,
    process_constraints,
    validate_constraint_pattern,
)


def _mock_constraint_functions(monkeypatch):
    monkeypatch.setattr(smilesgen, "can_ssbond", lambda seq: (seq, "SS"))
    monkeypatch.setattr(smilesgen, "can_htbond", lambda seq: (seq, "HT"))
    monkeypatch.setattr(smilesgen, "can_scntbond", lambda seq: (seq, "SCNT"))
    monkeypatch.setattr(smilesgen, "can_scctbond", lambda seq: (seq, "SCCT"))
    monkeypatch.setattr(smilesgen, "can_scscbond", lambda seq: (seq, "SCSC"))
    monkeypatch.setattr(smilesgen, "what_constraints", lambda seq: ["SS", "HT"])


def test_parse_fasta(tmp_path):
    fasta_content = ">seq1|SS\nACDE\n>seq2\nFGHI"
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)

    results = list(parse_fasta(fasta_file))
    assert results == [("ACDE", "SS"), ("FGHI", "")]


def test_constraint_resolver_valid_constraint(monkeypatch):
    _mock_constraint_functions(monkeypatch)
    seq, constr = constraint_resolver("ACDE", "SS")
    assert seq == "ACDE"
    assert constr == "SS"


def test_constraint_resolver_partial_constraint_found(monkeypatch):
    _mock_constraint_functions(monkeypatch)
    seq, constr = constraint_resolver("ACDE", "SCNT")
    assert seq == "ACDE"
    assert constr == "SCNT"


def test_constraint_resolver_partial_constraint_fallback(monkeypatch):
    _mock_constraint_functions(monkeypatch)
    seq, constr = constraint_resolver("ACDE", "SC")
    assert seq == "ACDE"
    assert constr in {"SCNT", "SCCT", "SCSC"}  # any of the SC* matches


def test_constraint_resolver_empty_constraint(monkeypatch):
    _mock_constraint_functions(monkeypatch)
    seq, constr = constraint_resolver("ACDE", "")
    assert seq == "ACDE"
    assert constr == ""


def test_constraint_resolver_invalid_constraint_raises(monkeypatch):
    _mock_constraint_functions(monkeypatch)
    with pytest.raises(InvalidConstraintError):
        constraint_resolver("ACDE", "INVALID")


def test_process_constraints_yields_expected(tmp_path, monkeypatch):
    _mock_constraint_functions(monkeypatch)
    fasta_content = ">seq1|SS\nACDE\n>seq2|HT\nFGHI"
    fasta_file = tmp_path / "test2.fasta"
    fasta_file.write_text(fasta_content)

    results = list(process_constraints(fasta_file))
    assert results == [("ACDE", "SS"), ("FGHI", "HT")]


def test_validate_constraint_pattern_rejects_mismatched_mask_length():
    is_valid, message = validate_constraint_pattern("AAAA", "SCZX")
    assert not is_valid
    assert "Mask length" in message


def test_validate_constraint_pattern_accepts_ester_sidechain_code():
    is_valid, message = validate_constraint_pattern("SAAA", "SCEXXX")
    assert is_valid
    assert "SCCT" in message


def test_validate_constraint_pattern_rejects_unknown_sidechain_code():
    is_valid, message = validate_constraint_pattern("AAAA", "SCQXXX")
    assert not is_valid
    assert "unknown constraint code" in message


def test_generate_smiles_strings_skips_invalid_sequences(tmp_path, monkeypatch, capsys):
    fasta_content = ">valid|HT\nGA\n>invalid|SCZX\nAAAA\n"
    fasta_file = tmp_path / "library.fasta"
    output_file = tmp_path / "out.p2smi"
    fasta_file.write_text(fasta_content)

    written = []

    monkeypatch.setattr(
        smilesgen,
        "constrained_peptide_smiles",
        lambda seq, constr: (seq, constr, f"{seq}-{constr or 'linear'}"),
    )

    def fake_write_library(entries, out_file, write, write_to_file):
        written.extend(list(entries))
        assert out_file == output_file
        assert write == "text"
        assert write_to_file is True

    monkeypatch.setattr(smilesgen, "write_library", fake_write_library)

    generate_smiles_strings(fasta_file, output_file, verbose=True)

    captured = capsys.readouterr()
    assert "Warning:" in captured.out
    assert "[DEBUG]" in captured.out
    assert written == [("GA", "HT", "GA-HT")]


def test_generate_smiles_strings_uses_custom_registry_for_constrained_sequence(
    tmp_path,
):
    registry_file = tmp_path / "residues.json"
    registry_file.write_text(
        '{"schema_version": 1, "residues": ['
        '{"id": "LAB_LYS", "name": "Lab lysine", '
        '"smiles": "N[C@@H](CCCCN)C(=O)O"}]}'
    )
    fasta_file = tmp_path / "custom.fasta"
    fasta_file.write_text(">custom|SCCT\nLAB_LYS,A,A,A\n")
    output_file = tmp_path / "custom.p2smi"

    generate_smiles_strings(fasta_file, output_file, registry_file=str(registry_file))

    output_smiles = output_file.read_text().strip().split(": ", 1)[1]
    assert Chem.MolFromSmiles(output_smiles) is not None
    assert "*" not in output_smiles


def test_generate_smiles_strings_applies_explicit_residue_recipe(tmp_path):
    modifier_file = tmp_path / "modifiers.json"
    modifier_file.write_text(
        '{"schema_version":1,"kind":"p2smi_modifier_registry",'
        '"modifiers":[],"recipes":[{"id":"designed","steps":['
        '{"modifier":"palmitoyl","site":{"residue":2,"role":"sidechain_amine"}},'
        '{"modifier":"n_methyl","site":{"residue":3,"role":"backbone_n"}},'
        '{"modifier":"c_amide","site":{"terminus":"C"}}]}]}'
    )
    fasta_file = tmp_path / "designed.fasta"
    fasta_file.write_text(">designed\nAKAA\n")
    output_file = tmp_path / "designed.p2smi"

    generate_smiles_strings(
        fasta_file,
        output_file,
        modifier_registry_file=str(modifier_file),
        recipe_id="designed",
    )

    label, output_smiles = output_file.read_text().strip().split(": ", 1)
    mol = Chem.MolFromSmiles(output_smiles)
    assert mol is not None
    assert "palmitoyl@2,n_methyl@3,c_amide" in label
    assert mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)N(C)[C@@H]"))
