import pytest

from p2smi.chemProps import (
    SmilesError,
    lipinski_pass,
    lipinski_trial,
    log_partition_coefficient,
    molecular_formula,
    molecule_summary,
    tpsa,
)


def test_log_partition_coefficient_valid():
    assert (
        round(
            log_partition_coefficient(
                "CCN(CC)C(=O)[C@H]1CN([C@@H]2CC3=CNC4=CC=CC(=C34)C2=C1)C"
            ),
            2,
        )
        == 2.91
    )  # ethanol approx


def test_log_partition_coefficient_invalid():
    with pytest.raises(SmilesError):
        log_partition_coefficient("INVALID_SMILES")


def test_lipinski_trial_pass_and_fail():
    passed, failed = lipinski_trial("CCO")  # simple ethanol
    assert "logP: " in next(p for p in passed if p.startswith("logP"))
    passed, failed = lipinski_trial(
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
    )  # very large
    assert any("Molecular weight" in f or "logP" in f for f in failed)


def test_lipinski_pass_true():
    assert lipinski_pass("CCO") is True


def test_lipinski_pass_false():
    big_hydrocarbon = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"  # fails MW and logP
    assert lipinski_pass(big_hydrocarbon) is False


def test_molecular_formula():
    assert molecular_formula("CCO") == "C2H6O"


def test_tpsa_known_value():
    assert round(tpsa("CCO"), 1) == 20.2  # ethanol TPSA


def test_molecule_summary_keys():
    summary = molecule_summary("CCO")
    expected_keys = {
        "Formula",
        "Molecular Weight",
        "logP",
        "TPSA",
        "H-bond donors",
        "H-bond acceptors",
        "Rotatable Bonds",
        "Rings",
        "Fraction Csp3",
        "Heavy Atoms",
        "Formal Charge",
        "Lipinski pass",
    }
    assert expected_keys.issubset(summary.keys())


def test_molecule_summary_invalid_smiles():
    with pytest.raises(SmilesError):
        molecule_summary("INVALID_SMILES")
