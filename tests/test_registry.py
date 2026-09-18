import json

import pytest

from p2smi.registry import (
    RegistryError,
    Residue,
    ResidueRegistry,
    parse_sequence_tokens,
)
from p2smi.utilities.aminoacids import all_aminos


def test_legacy_registry_accepts_letters_and_variable_length_codes():
    registry = ResidueRegistry.from_legacy_mapping(all_aminos)

    assert parse_sequence_tokens("GA", registry) == ["Gly", "Ala"]
    assert parse_sequence_tokens("PG,Ala", registry) == ["PG", "Ala"]


def test_undelimited_codes_are_not_guessed_as_variable_length_tokens():
    registry = ResidueRegistry.from_legacy_mapping(all_aminos)

    assert parse_sequence_tokens("PGA", registry) == ["Pro", "Gly", "Ala"]
    assert parse_sequence_tokens("PG,A", registry) == ["PG", "Ala"]


def test_single_variable_length_id_is_resolved_as_one_token():
    registry = ResidueRegistry(
        [Residue(id="LAB_001", name="Lab residue", smiles="N[C@@H](CC)C(=O)O")]
    )

    assert parse_sequence_tokens("LAB_001", registry) == ["LAB_001"]


def test_registry_json_round_trip_and_overlay(tmp_path):
    base = ResidueRegistry(
        [Residue(id="A", code="A", name="Alanine", smiles="N[C@@H](C)C(=O)O")]
    )
    overlay = ResidueRegistry(
        [
            Residue(
                id="LAB_001",
                name="Lab residue",
                smiles="N[C@@H](CC)C(=O)O",
            )
        ]
    )
    registry_file = tmp_path / "residues.json"
    overlay.write_json(registry_file)

    loaded = ResidueRegistry.read_json(registry_file)
    combined = base.merge(loaded)

    assert parse_sequence_tokens("A,LAB_001", combined) == ["A", "LAB_001"]
    assert json.loads(registry_file.read_text())["schema_version"] == 1


def test_registry_rejects_duplicate_aliases_and_invalid_ids():
    with pytest.raises(RegistryError, match="Duplicate residue alias"):
        ResidueRegistry(
            [
                Residue(id="ONE", code="DUP", name="One", smiles="CC"),
                Residue(id="TWO", code="DUP", name="Two", smiles="CCC"),
            ]
        )

    with pytest.raises(RegistryError, match="ASCII"):
        ResidueRegistry([Residue(id="LAB 001", name="Bad", smiles="CC")])

    with pytest.raises(RegistryError, match="invalid SMILES"):
        ResidueRegistry([Residue(id="BAD", name="Bad", smiles="not smiles")])

    with pytest.raises(RegistryError, match="empty name"):
        ResidueRegistry([Residue(id="BAD", name=3, smiles="CC")])

    with pytest.raises(RegistryError, match="commas"):
        ResidueRegistry([Residue(id="BAD", code="BAD,CODE", name="Bad", smiles="CC")])

    with pytest.raises(RegistryError, match="unknown fields"):
        ResidueRegistry.from_document(
            {
                "schema_version": 1,
                "residues": [
                    {"id": "BAD", "name": "Bad", "smiles": "CC", "formula": "C2H6"}
                ],
            }
        )

    with pytest.raises(RegistryError, match="must be an object"):
        ResidueRegistry.from_document([])
