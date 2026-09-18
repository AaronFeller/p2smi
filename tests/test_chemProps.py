import json

import pytest

from p2smi.chemProps import (
    SmilesError,
    lipinski_trial_mol,
    main,
    molecule_summary,
    make_mol,
    parse_smiles_line,
    process_line,
)


def test_log_partition_coefficient_valid():
    assert (
        round(
            molecule_summary("CCN(CC)C(=O)[C@H]1CN([C@@H]2CC3=CNC4=CC=CC(=C34)C2=C1)C")[
                "logP"
            ]
        )
        == 3
    )  # ethanol approx


def test_log_partition_coefficient_invalid():
    with pytest.raises(SmilesError):
        molecule_summary("INVALID_SMILES")


def lipinski_pass(mol):
    passed, failed = lipinski_trial_mol(mol)
    return len(failed) == 0


# --- Tests ---
def test_lipinski_pass_true():
    mol = make_mol("CCO")  # ethanol
    assert lipinski_pass(mol) is True


def test_lipinski_pass_false():
    big_hydrocarbon = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
    mol = make_mol(big_hydrocarbon)
    assert lipinski_pass(mol) is False


def test_molecular_formula():
    small_hydrocarbon = "CCO"  # ethanol
    assert molecule_summary(small_hydrocarbon)["Formula"] == "C2H6O"


def test_tpsa_known_value():
    assert round(molecule_summary("CCO")["TPSA"]) == 20  # ethanol TPSA


def test_molecule_summary_keys():
    summary = molecule_summary("CCO")
    expected_keys = {
        "SMILES",
        "Formula",
        "Molecular weight",
        "logP",
        "TPSA",
        "H-bond donors",
        "H-bond acceptors",
        "Rotatable bonds",
        "Rings",
        "Fraction Csp3",
        "Heavy atoms",
        "Formal charge",
        "Lipinski pass",
    }
    assert expected_keys.issubset(summary.keys())


def test_molecule_summary_invalid_smiles():
    with pytest.raises(SmilesError):
        molecule_summary("INVALID_SMILES")


def test_parse_smiles_line_preserves_optional_id():
    assert parse_smiles_line("mol1: CCO") == ("mol1", "CCO")
    assert parse_smiles_line("CCO") == (None, "CCO")


def test_process_line_can_include_ids_in_json():
    payload = json.loads(process_line("mol1: CCO", include_id=True))
    assert payload["ID"] == "mol1"
    assert payload["SMILES"] == "CCO"


def test_process_line_strict_raises_for_invalid_smiles():
    with pytest.raises(SmilesError):
        process_line("mol1: INVALID_SMILES", strict=True)


def test_main_batch_include_id_writes_jsonl(tmp_path):
    input_file = tmp_path / "mols.txt"
    output_file = tmp_path / "props.jsonl"
    input_file.write_text("mol1: CCO\nmol2: CCN\n")

    main(["-i", str(input_file), "-o", str(output_file), "--include_id"])

    rows = [json.loads(line) for line in output_file.read_text().splitlines()]
    assert [row["ID"] for row in rows] == ["mol1", "mol2"]
