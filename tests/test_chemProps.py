import pytest

from p2smi.chemProps import (
    SmilesError,
    lipinski_trial_mol,
    molecule_summary,
    make_mol,
    molecule_summary_from_mol,
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
