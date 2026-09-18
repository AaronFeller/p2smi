"""Versioned residue registries and unambiguous peptide sequence tokens."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional

from rdkit import Chem


SCHEMA_VERSION = 1
IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_-]*$")
class RegistryError(ValueError):
    """Raised when residue registry data or sequence tokens are invalid."""


@dataclass(frozen=True)
class Residue:
    """A residue definition with a stable, ASCII identifier."""

    id: str
    name: str
    smiles: str
    code: Optional[str] = None
    legacy_letter: Optional[str] = None

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "Residue":
        if not isinstance(data, Mapping):
            raise RegistryError("Each residue entry must be an object")
        allowed_fields = {"id", "name", "smiles", "code", "legacy_letter"}
        unknown_fields = set(data) - allowed_fields
        if unknown_fields:
            raise RegistryError(f"Residue entry has unknown fields: {sorted(unknown_fields)}")
        try:
            identifier = data["id"]
            name = data["name"]
            smiles = data["smiles"]
        except KeyError as exc:
            raise RegistryError(f"Residue entry is missing required field {exc.args[0]!r}") from exc

        return cls(
            id=identifier,
            name=name,
            smiles=smiles,
            code=data.get("code"),
            legacy_letter=data.get("legacy_letter"),
        )

    def to_dict(self) -> dict[str, Any]:
        entry = {"id": self.id, "name": self.name, "smiles": self.smiles}
        optional_fields = {
            "code": self.code,
            "legacy_letter": self.legacy_letter,
        }
        entry.update({key: value for key, value in optional_fields.items() if value is not None})
        return entry


class ResidueRegistry:
    """Validated registry with stable IDs and optional legacy aliases."""

    def __init__(self, residues: Iterable[Residue] = ()):
        self._by_id: dict[str, Residue] = {}
        self._aliases: dict[str, str] = {}
        for residue in residues:
            self.add(residue)

    def __len__(self) -> int:
        return len(self._by_id)

    def __iter__(self):
        return iter(self._by_id.values())

    def add(self, residue: Residue) -> None:
        self.validate_residue(residue)
        if residue.id in self._by_id or residue.id in self._aliases:
            raise RegistryError(f"Duplicate residue identifier: {residue.id}")

        aliases = [
            alias
            for alias in (residue.code, residue.legacy_letter)
            if alias and alias != residue.id
        ]
        for alias in aliases:
            if alias in self._by_id or alias in self._aliases:
                raise RegistryError(f"Duplicate residue alias: {alias}")

        self._by_id[residue.id] = residue
        for alias in aliases:
            self._aliases[alias] = residue.id

    def resolve(self, token: str) -> Residue:
        try:
            return self._by_id[token]
        except KeyError:
            try:
                return self._by_id[self._aliases[token]]
            except KeyError as exc:
                raise RegistryError(f"Unknown residue token: {token}") from exc

    def merge(self, overlay: "ResidueRegistry") -> "ResidueRegistry":
        return ResidueRegistry([*self, *overlay])

    def to_document(self) -> dict[str, Any]:
        return {
            "schema_version": SCHEMA_VERSION,
            "residues": [residue.to_dict() for residue in self],
        }

    def write_json(self, file_path: Path) -> None:
        file_path.parent.mkdir(parents=True, exist_ok=True)
        file_path.write_text(
            json.dumps(self.to_document(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @classmethod
    def from_document(cls, document: Mapping[str, Any]) -> "ResidueRegistry":
        if not isinstance(document, Mapping):
            raise RegistryError("Registry document must be an object")
        unknown_fields = set(document) - {"schema_version", "residues"}
        if unknown_fields:
            raise RegistryError(f"Registry document has unknown fields: {sorted(unknown_fields)}")
        if document.get("schema_version") != SCHEMA_VERSION:
            raise RegistryError(
                f"Unsupported registry schema version: {document.get('schema_version')!r}"
            )
        residues = document.get("residues")
        if not isinstance(residues, list):
            raise RegistryError("Registry document must contain a residues list")
        return cls(Residue.from_dict(entry) for entry in residues)

    @classmethod
    def read_json(cls, file_path: Path) -> "ResidueRegistry":
        try:
            document = json.loads(file_path.read_text(encoding="utf-8"))
        except FileNotFoundError as exc:
            raise RegistryError(f"Registry file does not exist: {file_path}") from exc
        except json.JSONDecodeError as exc:
            raise RegistryError(f"Invalid JSON registry {file_path}: {exc}") from exc
        return cls.from_document(document)

    @classmethod
    def from_legacy_mapping(cls, mapping: Mapping[str, Mapping[str, Any]]) -> "ResidueRegistry":
        residues = []
        for name, props in mapping.items():
            code = props["Code"]
            identifier = code if IDENTIFIER_PATTERN.fullmatch(code) else f"legacy_{code}"
            residues.append(
                Residue(
                    id=identifier,
                    code=code,
                    name=name,
                    smiles=props["SMILES"],
                    legacy_letter=props.get("Letter"),
                )
            )
        return cls(residues)

    @staticmethod
    def validate_residue(residue: Residue) -> None:
        if not isinstance(residue.id, str):
            raise RegistryError("Residue id must be a string")
        if not IDENTIFIER_PATTERN.fullmatch(residue.id):
            raise RegistryError(
                "Residue id must be ASCII and use letters, digits, underscores, or hyphens"
            )
        if not isinstance(residue.name, str) or not residue.name.strip():
            raise RegistryError(f"Residue {residue.id} has an empty name")
        if not isinstance(residue.smiles, str) or not residue.smiles.strip():
            raise RegistryError(f"Residue {residue.id} has an empty SMILES string")
        if Chem.MolFromSmiles(residue.smiles) is None:
            raise RegistryError(f"Residue {residue.id} has invalid SMILES")
        for field_name, value in (("code", residue.code), ("legacy_letter", residue.legacy_letter)):
            if value is not None and not isinstance(value, str):
                raise RegistryError(f"{field_name} must be a string")
            if value and (value.strip() != value or "," in value):
                raise RegistryError(f"{field_name} cannot contain commas or surrounding whitespace")
        if residue.legacy_letter and len(residue.legacy_letter) != 1:
            raise RegistryError("legacy_letter must be exactly one character")


def parse_sequence_tokens(sequence: str | Iterable[str], registry: ResidueRegistry) -> list[str]:
    """Return canonical residue IDs from a delimited v2 or legacy v1 sequence.

    V2 strings must be comma-delimited. Undelimited strings are interpreted only
    as legacy one-character aliases; variable-length IDs are never guessed.
    """

    if isinstance(sequence, str):
        if not sequence:
            return []
        if "," in sequence:
            raw_tokens = sequence.split(",")
        else:
            try:
                return [registry.resolve(sequence).id]
            except RegistryError:
                raw_tokens = list(sequence)
    else:
        raw_tokens = list(sequence)

    if not raw_tokens or any(not token or token.strip() != token for token in raw_tokens):
        raise RegistryError("Sequence tokens must be non-empty and contain no surrounding whitespace")

    return [registry.resolve(token).id for token in raw_tokens]