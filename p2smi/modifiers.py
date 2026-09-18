"""Deterministic peptide modification recipes built on RDKit molecular graphs."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional

from rdkit import Chem
from rdkit.Chem import rdChemReactions


SCHEMA_VERSION = 1
SITE_ROLES = {
    "n_terminal_amine",
    "c_terminal_carboxyl",
    "backbone_n",
    "backbone_amide_n",
    "sidechain_amine",
    "sidechain_hydroxyl",
    "sidechain_thiol",
}

LIPID_LINKERS = {
    "aeea": "*C(=O)COCCOCCN[1*]",
    "aeea2": "*C(=O)COCCOCCNC(=O)COCCOCCN[1*]",
    "gamma_glu": "*C(=O)CC[C@H](N[1*])C(=O)O",
    "gamma_glu_aeea": ("*C(=O)CC[C@H](NC(=O)COCCOCCN[1*])C(=O)O"),
    "gamma_glu_aeea2": ("*C(=O)CC[C@H](NC(=O)COCCOCCNC(=O)COCCOCCN[1*])C(=O)O"),
}


class ModificationError(ValueError):
    """Raised when a modifier definition, selector, or product is invalid."""


@dataclass(frozen=True)
class SiteSelector:
    residue: Optional[int] = None
    role: Optional[str] = None
    terminus: Optional[str] = None
    atom_index: Optional[int] = None
    smarts: Optional[str] = None
    match_index: Optional[int] = None

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "SiteSelector":
        if not isinstance(data, Mapping):
            raise ModificationError("A site selector must be an object")
        allowed = {"residue", "role", "terminus", "atom_index", "smarts", "match_index"}
        unknown = set(data) - allowed
        if unknown:
            raise ModificationError(f"Unknown site selector fields: {sorted(unknown)}")
        selector = cls(**data)
        selector.validate()
        return selector

    def validate(self) -> None:
        modes = sum(
            value is not None
            for value in (self.role, self.terminus, self.atom_index, self.smarts)
        )
        if modes != 1:
            raise ModificationError(
                "A site selector must specify exactly one selection mode"
            )
        if self.residue is not None and (
            not isinstance(self.residue, int)
            or isinstance(self.residue, bool)
            or self.residue < 1
            or self.terminus is not None
        ):
            raise ModificationError(
                "residue must be a positive integer combined with role, "
                "atom_index, or SMARTS"
            )
        if self.role is not None and self.role not in SITE_ROLES:
            raise ModificationError(f"Unknown site role: {self.role}")
        if self.terminus not in {None, "N", "C"}:
            raise ModificationError("terminus must be 'N' or 'C'")
        if self.atom_index is not None and (
            not isinstance(self.atom_index, int)
            or isinstance(self.atom_index, bool)
            or self.atom_index < 0
        ):
            raise ModificationError("atom_index must be a non-negative integer")
        if self.smarts is not None and Chem.MolFromSmarts(self.smarts) is None:
            raise ModificationError(f"Invalid selector SMARTS: {self.smarts}")
        if self.match_index is not None and (
            self.smarts is None
            or not isinstance(self.match_index, int)
            or isinstance(self.match_index, bool)
            or self.match_index < 0
        ):
            raise ModificationError(
                "match_index requires SMARTS and must be non-negative"
            )

    def to_dict(self) -> dict[str, Any]:
        fields = {
            "residue": self.residue,
            "role": self.role,
            "terminus": self.terminus,
            "atom_index": self.atom_index,
            "smarts": self.smarts,
            "match_index": self.match_index,
        }
        return {key: value for key, value in fields.items() if value is not None}


@dataclass(frozen=True)
class ModifierDefinition:
    id: str
    name: str
    operation: str
    allowed_roles: tuple[str, ...]
    fragment_smiles: Optional[str] = None
    reaction_smarts: Optional[str] = None

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ModifierDefinition":
        if not isinstance(data, Mapping):
            raise ModificationError("Each modifier must be an object")
        allowed = {
            "id",
            "name",
            "operation",
            "allowed_roles",
            "fragment_smiles",
            "reaction_smarts",
        }
        unknown = set(data) - allowed
        if unknown:
            raise ModificationError(f"Unknown modifier fields: {sorted(unknown)}")
        try:
            definition = cls(
                id=data["id"],
                name=data["name"],
                operation=data["operation"],
                allowed_roles=tuple(data["allowed_roles"]),
                fragment_smiles=data.get("fragment_smiles"),
                reaction_smarts=data.get("reaction_smarts"),
            )
        except (KeyError, TypeError) as exc:
            raise ModificationError(f"Malformed modifier definition: {exc}") from exc
        definition.validate()
        return definition

    def validate(self) -> None:
        if (
            not isinstance(self.id, str)
            or not self.id
            or not self.id.replace("_", "").isalnum()
        ):
            raise ModificationError(
                "Modifier id must contain only letters, digits, and underscores"
            )
        if not isinstance(self.name, str) or not self.name.strip():
            raise ModificationError(f"Modifier {self.id} has an empty name")
        if self.operation not in {
            "graft",
            "n_methyl",
            "c_amide",
            "mpeg_acyl",
            "reaction",
        }:
            raise ModificationError(f"Unsupported modifier operation: {self.operation}")
        if not self.allowed_roles or any(
            role not in SITE_ROLES for role in self.allowed_roles
        ):
            raise ModificationError(f"Modifier {self.id} has invalid allowed_roles")
        if self.operation == "graft":
            _validate_graft_fragment(self.fragment_smiles, self.id)
        elif self.fragment_smiles is not None:
            raise ModificationError(
                f"Modifier {self.id} does not accept fragment_smiles"
            )
        if self.operation == "reaction":
            _validate_reaction_smarts(self.reaction_smarts, self.id)
        elif self.reaction_smarts is not None:
            raise ModificationError(
                f"Modifier {self.id} does not accept reaction_smarts"
            )

    def to_dict(self) -> dict[str, Any]:
        data = {
            "id": self.id,
            "name": self.name,
            "operation": self.operation,
            "allowed_roles": list(self.allowed_roles),
        }
        if self.fragment_smiles is not None:
            data["fragment_smiles"] = self.fragment_smiles
        if self.reaction_smarts is not None:
            data["reaction_smarts"] = self.reaction_smarts
        return data


@dataclass(frozen=True)
class ModificationStep:
    modifier: str
    site: SiteSelector
    params: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ModificationStep":
        if not isinstance(data, Mapping):
            raise ModificationError("Each recipe step must be an object")
        unknown = set(data) - {"modifier", "site", "params"}
        if unknown:
            raise ModificationError(f"Unknown recipe step fields: {sorted(unknown)}")
        try:
            modifier = data["modifier"]
            site = SiteSelector.from_dict(data["site"])
        except KeyError as exc:
            raise ModificationError(f"Recipe step is missing {exc.args[0]!r}") from exc
        params = data.get("params", {})
        if not isinstance(params, Mapping):
            raise ModificationError("Recipe step params must be an object")
        return cls(modifier=modifier, site=site, params=dict(params))

    def to_dict(self) -> dict[str, Any]:
        data = {"modifier": self.modifier, "site": self.site.to_dict()}
        if self.params:
            data["params"] = dict(self.params)
        return data


@dataclass(frozen=True)
class Recipe:
    id: str
    steps: tuple[ModificationStep, ...]

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "Recipe":
        if not isinstance(data, Mapping) or set(data) - {"id", "steps"}:
            raise ModificationError("Recipe must contain only id and steps")
        try:
            recipe = cls(
                data["id"],
                tuple(ModificationStep.from_dict(step) for step in data["steps"]),
            )
        except (KeyError, TypeError) as exc:
            raise ModificationError(f"Malformed recipe: {exc}") from exc
        if not recipe.id or not recipe.steps:
            raise ModificationError("Recipe id and at least one step are required")
        return recipe

    def to_dict(self) -> dict[str, Any]:
        return {"id": self.id, "steps": [step.to_dict() for step in self.steps]}


class ModifierRegistry:
    def __init__(
        self,
        modifiers: Iterable[ModifierDefinition] = (),
        recipes: Iterable[Recipe] = (),
    ):
        self.modifiers: dict[str, ModifierDefinition] = {}
        self.recipes: dict[str, Recipe] = {}
        for modifier in modifiers:
            if modifier.id in self.modifiers:
                raise ModificationError(f"Duplicate modifier id: {modifier.id}")
            modifier.validate()
            self.modifiers[modifier.id] = modifier
        for recipe in recipes:
            if recipe.id in self.recipes:
                raise ModificationError(f"Duplicate recipe id: {recipe.id}")
            unknown = [
                step.modifier
                for step in recipe.steps
                if step.modifier not in self.modifiers
            ]
            if unknown:
                raise ModificationError(
                    f"Recipe {recipe.id} uses unknown modifiers: {unknown}"
                )
            self.recipes[recipe.id] = recipe

    def merge(self, overlay: "ModifierRegistry") -> "ModifierRegistry":
        return ModifierRegistry(
            [*self.modifiers.values(), *overlay.modifiers.values()],
            [*self.recipes.values(), *overlay.recipes.values()],
        )

    def to_document(self) -> dict[str, Any]:
        return {
            "schema_version": SCHEMA_VERSION,
            "kind": "p2smi_modifier_registry",
            "modifiers": [modifier.to_dict() for modifier in self.modifiers.values()],
            "recipes": [recipe.to_dict() for recipe in self.recipes.values()],
        }

    def write_json(self, file_path: Path) -> None:
        file_path.parent.mkdir(parents=True, exist_ok=True)
        file_path.write_text(
            json.dumps(self.to_document(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @classmethod
    def from_document(
        cls,
        document: Mapping[str, Any],
        base: Optional["ModifierRegistry"] = None,
    ) -> "ModifierRegistry":
        if not isinstance(document, Mapping):
            raise ModificationError("Modifier registry must be an object")
        unknown = set(document) - {"schema_version", "kind", "modifiers", "recipes"}
        if unknown:
            raise ModificationError(
                f"Unknown modifier registry fields: {sorted(unknown)}"
            )
        if document.get("schema_version") != SCHEMA_VERSION:
            raise ModificationError("Unsupported modifier registry schema version")
        if document.get("kind") != "p2smi_modifier_registry":
            raise ModificationError("Invalid modifier registry kind")
        modifiers = document.get("modifiers", [])
        recipes = document.get("recipes", [])
        if not isinstance(modifiers, list) or not isinstance(recipes, list):
            raise ModificationError("modifiers and recipes must be lists")
        parsed_modifiers = [ModifierDefinition.from_dict(item) for item in modifiers]
        parsed_recipes = [Recipe.from_dict(item) for item in recipes]
        if base is None:
            return cls(parsed_modifiers, parsed_recipes)
        return cls(
            [*base.modifiers.values(), *parsed_modifiers],
            [*base.recipes.values(), *parsed_recipes],
        )

    @classmethod
    def read_json(
        cls,
        file_path: Path,
        base: Optional["ModifierRegistry"] = None,
    ) -> "ModifierRegistry":
        try:
            document = json.loads(file_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ModificationError(
                f"Could not read modifier registry {file_path}: {exc}"
            ) from exc
        return cls.from_document(document, base=base)


def _validate_graft_fragment(fragment_smiles: Optional[str], modifier_id: str) -> None:
    if not isinstance(fragment_smiles, str):
        raise ModificationError(
            f"Graft modifier {modifier_id} requires fragment_smiles"
        )
    mol = Chem.MolFromSmiles(fragment_smiles)
    if mol is None:
        raise ModificationError(f"Modifier {modifier_id} has invalid fragment SMILES")
    dummies = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
    if len(dummies) != 1 or dummies[0].GetDegree() != 1:
        raise ModificationError(
            f"Modifier {modifier_id} requires one terminal '*' attachment atom"
        )


def _validate_reaction_smarts(reaction_smarts: Optional[str], modifier_id: str) -> None:
    if not isinstance(reaction_smarts, str) or not reaction_smarts.strip():
        raise ModificationError(
            f"Reaction modifier {modifier_id} requires reaction_smarts"
        )
    if len(reaction_smarts) > 4096:
        raise ModificationError(f"Reaction SMARTS for {modifier_id} is too long")
    try:
        reaction = rdChemReactions.ReactionFromSmarts(reaction_smarts)
    except Exception as exc:
        raise ModificationError(
            f"Modifier {modifier_id} has invalid reaction SMARTS"
        ) from exc
    if reaction is None:
        raise ModificationError(f"Modifier {modifier_id} has invalid reaction SMARTS")
    if (
        reaction.GetNumReactantTemplates() != 1
        or reaction.GetNumProductTemplates() != 1
    ):
        raise ModificationError(
            f"Reaction modifier {modifier_id} requires one reactant and one product"
        )
    reactant = reaction.GetReactantTemplate(0)
    product = reaction.GetProductTemplate(0)
    if reactant.GetNumAtoms() > 256 or product.GetNumAtoms() > 256:
        raise ModificationError(f"Reaction SMARTS for {modifier_id} is too large")
    reactant_maps = [
        atom.GetAtomMapNum() for atom in reactant.GetAtoms() if atom.GetAtomMapNum()
    ]
    product_maps = [
        atom.GetAtomMapNum() for atom in product.GetAtoms() if atom.GetAtomMapNum()
    ]
    if reactant_maps.count(1) != 1 or product_maps.count(1) != 1:
        raise ModificationError(
            f"Reaction modifier {modifier_id} requires atom map :1 once on each side"
        )
    if len(reactant_maps) != len(set(reactant_maps)) or len(product_maps) != len(
        set(product_maps)
    ):
        raise ModificationError(
            f"Reaction modifier {modifier_id} has duplicate atom maps"
        )
    if not set(product_maps).issubset(reactant_maps):
        raise ModificationError(
            f"Reaction modifier {modifier_id} has product atom maps absent "
            "from the reactant"
        )


def builtin_modifier_registry() -> ModifierRegistry:
    amines = ("n_terminal_amine", "sidechain_amine")
    definitions = [
        ModifierDefinition(
            "n_methyl",
            "Backbone N-methylation",
            "n_methyl",
            ("backbone_n", "backbone_amide_n"),
        ),
        ModifierDefinition(
            "n_acetyl",
            "N-terminal acetylation",
            "graft",
            ("n_terminal_amine",),
            "CC(=O)*",
        ),
        ModifierDefinition(
            "c_amide", "C-terminal amidation", "c_amide", ("c_terminal_carboxyl",)
        ),
        ModifierDefinition("mpeg_acetyl", "Methoxy-PEG acylation", "mpeg_acyl", amines),
        ModifierDefinition(
            "octanoyl", "Octanoylation", "graft", amines, "CCCCCCCC(=O)*"
        ),
        ModifierDefinition(
            "decanoyl", "Decanoylation", "graft", amines, "CCCCCCCCCC(=O)*"
        ),
        ModifierDefinition(
            "dodecanoyl", "Dodecanoylation", "graft", amines, "CCCCCCCCCCCC(=O)*"
        ),
        ModifierDefinition(
            "tetradecanoyl",
            "Tetradecanoylation",
            "graft",
            amines,
            "CCCCCCCCCCCCCC(=O)*",
        ),
        ModifierDefinition(
            "palmitoyl", "Palmitoylation", "graft", amines, "CCCCCCCCCCCCCCCC(=O)*"
        ),
        ModifierDefinition(
            "stearoyl", "Stearoylation", "graft", amines, "CCCCCCCCCCCCCCCCCC(=O)*"
        ),
        ModifierDefinition(
            "eicosanoyl",
            "Eicosanoylation",
            "graft",
            amines,
            "CCCCCCCCCCCCCCCCCCCC(=O)*",
        ),
        ModifierDefinition(
            "hexadecanedioyl",
            "Hexadecanedioic acid monoacylation",
            "graft",
            amines,
            "O=C(O)" + "C" * 14 + "C(=O)*",
        ),
        ModifierDefinition(
            "octadecanedioyl",
            "Octadecanedioic acid monoacylation",
            "graft",
            amines,
            "O=C(O)" + "C" * 16 + "C(=O)*",
        ),
        ModifierDefinition(
            "icosanedioyl",
            "Icosanedioic acid monoacylation",
            "graft",
            amines,
            "O=C(O)" + "C" * 18 + "C(=O)*",
        ),
        ModifierDefinition(
            "carbamidomethyl",
            "Thiol carbamidomethylation",
            "graft",
            ("sidechain_thiol",),
            "NC(=O)C*",
        ),
        ModifierDefinition(
            "phosphoryl",
            "Hydroxyl phosphorylation",
            "graft",
            ("sidechain_hydroxyl",),
            "P(=O)(O)(O)*",
        ),
        ModifierDefinition(
            "azidoacetyl",
            "Azidoacetyl click handle",
            "graft",
            amines,
            "[N-]=[N+]=NCC(=O)*",
        ),
        ModifierDefinition(
            "propiolyl", "Propiolyl click handle", "graft", amines, "C#CC(=O)*"
        ),
    ]
    return ModifierRegistry(definitions)


def _carboxyl_groups(mol: Chem.Mol) -> list[tuple[int, int]]:
    pattern = Chem.MolFromSmarts("[CX3:1](=[OX1])[OX2H1:2]")
    return [(match[0], match[2]) for match in mol.GetSubstructMatches(pattern)]


def _n_terminal_indices(mol: Chem.Mol) -> list[int]:
    candidates = []
    for match in mol.GetSubstructMatches(
        Chem.MolFromSmarts("[N;H1,H2;!$(N-C=O):1]-[C:2]-[C:3](=O)")
    ):
        candidates.append(match[0])
    return sorted(set(candidates))


def _c_terminal_indices(mol: Chem.Mol) -> list[int]:
    candidates = []
    for carbonyl_idx, _ in _carboxyl_groups(mol):
        carbonyl = mol.GetAtomWithIdx(carbonyl_idx)
        for alpha in carbonyl.GetNeighbors():
            if alpha.GetAtomicNum() != 6:
                continue
            if any(neighbor.GetAtomicNum() == 7 for neighbor in alpha.GetNeighbors()):
                candidates.append(carbonyl_idx)
                break
    return sorted(set(candidates))


def _role_indices(mol: Chem.Mol, role: str) -> list[int]:
    n_term = set(_n_terminal_indices(mol))
    if role == "n_terminal_amine":
        return sorted(n_term)
    if role == "c_terminal_carboxyl":
        return _c_terminal_indices(mol)
    if role == "backbone_n":
        return sorted(n_term)
    if role == "backbone_amide_n":
        pattern = Chem.MolFromSmarts("[N;H1;X3;$([N]-[C](=O));$([N]-[C]-[C](=O))]")
        return sorted({match[0] for match in mol.GetSubstructMatches(pattern)})
    if role == "sidechain_amine":
        pattern = Chem.MolFromSmarts("[N;H1,H2;!a;!$(N-C=O):1]")
        return [
            match[0]
            for match in mol.GetSubstructMatches(pattern)
            if match[0] not in n_term
        ]
    if role == "sidechain_hydroxyl":
        acid_oxygens = {oxygen for _, oxygen in _carboxyl_groups(mol)}
        pattern = Chem.MolFromSmarts("[O;H1;X2:1]")
        return [
            match[0]
            for match in mol.GetSubstructMatches(pattern)
            if match[0] not in acid_oxygens
        ]
    if role == "sidechain_thiol":
        pattern = Chem.MolFromSmarts("[S;H1:1]")
        return [match[0] for match in mol.GetSubstructMatches(pattern)]
    raise ModificationError(f"Unknown site role: {role}")


def resolve_site(mol: Chem.Mol, selector: SiteSelector) -> tuple[int, str]:
    selector.validate()
    if selector.residue is not None:
        raise ModificationError(
            "Residue-scoped selectors require sequence-aware assembly"
        )
    if selector.terminus:
        role = "n_terminal_amine" if selector.terminus == "N" else "c_terminal_carboxyl"
        indices = _role_indices(mol, role)
    elif selector.role:
        role = selector.role
        indices = _role_indices(mol, role)
    elif selector.atom_index is not None:
        if selector.atom_index >= mol.GetNumAtoms():
            raise ModificationError(
                f"atom_index {selector.atom_index} is outside the molecule"
            )
        return selector.atom_index, "atom_index"
    else:
        role = "smarts"
        pattern = Chem.MolFromSmarts(selector.smarts)
        matches = sorted(mol.GetSubstructMatches(pattern))
        if selector.match_index is not None:
            if selector.match_index >= len(matches):
                raise ModificationError("SMARTS match_index is outside the match list")
            return matches[selector.match_index][0], role
        indices = [match[0] for match in matches]

    if len(indices) != 1:
        raise ModificationError(
            f"Site selector resolved to {len(indices)} atoms; expected exactly one"
        )
    return indices[0], role


def _roles_for_atom(mol: Chem.Mol, atom_idx: int) -> set[str]:
    roles = {role for role in SITE_ROLES if atom_idx in _role_indices(mol, role)}
    if "backbone_n" in roles and atom_idx in _role_indices(mol, "backbone_amide_n"):
        roles.add("backbone_amide_n")
    return roles


def _mark_original_indices(mol: Chem.Mol) -> Chem.Mol:
    marked = Chem.Mol(mol)
    for atom in marked.GetAtoms():
        atom.SetIntProp("p2smi_original_index", atom.GetIdx())
    return marked


def _stable_selector(mol: Chem.Mol, selector: SiteSelector) -> SiteSelector:
    if selector.atom_index is None:
        return selector
    matches = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.HasProp("p2smi_original_index")
        and atom.GetIntProp("p2smi_original_index") == selector.atom_index
    ]
    if len(matches) != 1:
        raise ModificationError(
            f"Original atom_index {selector.atom_index} is no longer available"
        )
    return SiteSelector(atom_index=matches[0])


def _finalize_product(rw_mol: Chem.RWMol) -> Chem.Mol:
    product = rw_mol.GetMol()
    try:
        Chem.SanitizeMol(product)
    except Exception as exc:
        raise ModificationError(
            f"Modification produced invalid valence: {exc}"
        ) from exc
    if len(Chem.GetMolFrags(product)) != 1:
        raise ModificationError("Modification produced a disconnected molecule")
    for atom in product.GetAtoms():
        atom.SetAtomMapNum(0)
    return product


def _graft(mol: Chem.Mol, site_idx: int, fragment_smiles: str) -> Chem.Mol:
    fragment = Chem.MolFromSmiles(fragment_smiles)
    dummy_idx = next(
        atom.GetIdx() for atom in fragment.GetAtoms() if atom.GetAtomicNum() == 0
    )
    anchor_idx = fragment.GetAtomWithIdx(dummy_idx).GetNeighbors()[0].GetIdx()
    combined = Chem.CombineMols(mol, fragment)
    offset = mol.GetNumAtoms()
    rw_mol = Chem.RWMol(combined)
    rw_mol.AddBond(site_idx, offset + anchor_idx, Chem.BondType.SINGLE)
    rw_mol.RemoveAtom(offset + dummy_idx)
    return _finalize_product(rw_mol)


def _linked_acyl_fragment(fragment_smiles: str, linker_id: str) -> str:
    try:
        linker = Chem.MolFromSmiles(LIPID_LINKERS[linker_id])
    except KeyError as exc:
        raise ModificationError(f"Unknown lipid linker: {linker_id}") from exc
    cargo = Chem.MolFromSmiles(fragment_smiles)
    cargo_dummies = [atom for atom in cargo.GetAtoms() if atom.GetAtomicNum() == 0]
    if (
        len(cargo_dummies) != 1
        or cargo_dummies[0].GetDegree() != 1
        or cargo_dummies[0].GetNeighbors()[0].GetAtomicNum() != 6
        or not any(
            bond.GetBondType() == Chem.BondType.DOUBLE
            and bond.GetOtherAtom(cargo_dummies[0].GetNeighbors()[0]).GetAtomicNum()
            == 8
            for bond in cargo_dummies[0].GetNeighbors()[0].GetBonds()
        )
    ):
        raise ModificationError(
            "Lipid linkers require an acyl modifier with fragment R-C(=O)*"
        )

    linker_cargo_dummies = [
        atom
        for atom in linker.GetAtoms()
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == 1
    ]
    if len(linker_cargo_dummies) != 1:
        raise ModificationError(
            f"Lipid linker {linker_id} has invalid attachment chemistry"
        )
    linker_dummy = linker_cargo_dummies[0]
    linker_anchor = linker_dummy.GetNeighbors()[0].GetIdx()
    cargo_dummy = cargo_dummies[0]
    cargo_anchor = cargo_dummy.GetNeighbors()[0].GetIdx()

    combined = Chem.CombineMols(linker, cargo)
    offset = linker.GetNumAtoms()
    rw_mol = Chem.RWMol(combined)
    rw_mol.AddBond(linker_anchor, offset + cargo_anchor, Chem.BondType.SINGLE)
    for atom_idx in sorted(
        (linker_dummy.GetIdx(), offset + cargo_dummy.GetIdx()), reverse=True
    ):
        rw_mol.RemoveAtom(atom_idx)
    product = rw_mol.GetMol()
    try:
        Chem.SanitizeMol(product)
    except Exception as exc:
        raise ModificationError(
            f"Lipid linker produced invalid valence: {exc}"
        ) from exc
    linked_smiles = Chem.MolToSmiles(product, isomericSmiles=True)
    _validate_graft_fragment(linked_smiles, f"{linker_id} linked cargo")
    return linked_smiles


def _run_custom_reaction(
    mol: Chem.Mol,
    site_idx: int,
    reaction_smarts: str,
) -> Chem.Mol:
    reaction = rdChemReactions.ReactionFromSmarts(reaction_smarts)
    product_sets = reaction.RunReactants((mol,), maxProducts=257)
    if len(product_sets) >= 257:
        raise ModificationError("Custom reaction generated too many candidate products")

    source_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    products = {}
    for product_set in product_sets:
        product = product_set[0]
        selected_atoms = [
            atom
            for atom in product.GetAtoms()
            if atom.HasProp("old_mapno")
            and atom.GetIntProp("old_mapno") == 1
            and atom.HasProp("react_atom_idx")
            and atom.GetIntProp("react_atom_idx") == site_idx
        ]
        if len(selected_atoms) != 1:
            continue
        for atom in product.GetAtoms():
            if atom.HasProp("react_atom_idx"):
                source_atom = mol.GetAtomWithIdx(atom.GetIntProp("react_atom_idx"))
                if source_atom.HasProp("p2smi_original_index"):
                    atom.SetIntProp(
                        "p2smi_original_index",
                        source_atom.GetIntProp("p2smi_original_index"),
                    )
            for property_name in (
                "old_mapno",
                "react_atom_idx",
                "react_idx",
                "_ReactionDegreeChanged",
            ):
                if atom.HasProp(property_name):
                    atom.ClearProp(property_name)
        try:
            product = _finalize_product(Chem.RWMol(product))
        except ModificationError:
            continue
        if any(atom.GetAtomicNum() == 0 for atom in product.GetAtoms()):
            continue
        product_smiles = Chem.MolToSmiles(product, isomericSmiles=True)
        if product_smiles != source_smiles:
            products[product_smiles] = product

    if not products:
        raise ModificationError(
            "Custom reaction produced no valid product at the selected site"
        )
    if len(products) != 1:
        raise ModificationError(
            f"Custom reaction produced {len(products)} products at the selected site"
        )
    return next(iter(products.values()))


def apply_modifier(
    mol: Chem.Mol,
    definition: ModifierDefinition,
    selector: SiteSelector,
    params: Optional[Mapping[str, Any]] = None,
) -> Chem.Mol:
    params = {} if params is None else dict(params)
    site_idx, role = resolve_site(mol, selector)
    if role in {"atom_index", "smarts"}:
        actual_roles = _roles_for_atom(mol, site_idx)
        if not actual_roles.intersection(definition.allowed_roles):
            raise ModificationError(
                f"Modifier {definition.id} cannot target atom {site_idx}; "
                f"expected one of {definition.allowed_roles}"
            )
    elif role not in definition.allowed_roles:
        raise ModificationError(f"Modifier {definition.id} cannot target role {role}")

    if definition.operation == "graft":
        if not params:
            return _graft(mol, site_idx, definition.fragment_smiles)
        if set(params) != {"linker"} or not isinstance(params["linker"], str):
            raise ModificationError(
                f"Modifier {definition.id} accepts only a string linker parameter"
            )
        fragment = _linked_acyl_fragment(definition.fragment_smiles, params["linker"])
        return _graft(mol, site_idx, fragment)
    if definition.operation == "n_methyl":
        if params:
            raise ModificationError("n_methyl does not accept params")
        rw_mol = Chem.RWMol(mol)
        methyl_idx = rw_mol.AddAtom(Chem.Atom("C"))
        rw_mol.AddBond(site_idx, methyl_idx, Chem.BondType.SINGLE)
        return _finalize_product(rw_mol)
    if definition.operation == "c_amide":
        if params:
            raise ModificationError("c_amide does not accept params")
        groups = dict(_carboxyl_groups(mol))
        if site_idx not in groups:
            raise ModificationError(
                "C-terminal amidation requires a free carboxylic acid"
            )
        rw_mol = Chem.RWMol(mol)
        nitrogen_idx = rw_mol.AddAtom(Chem.Atom("N"))
        rw_mol.AddBond(site_idx, nitrogen_idx, Chem.BondType.SINGLE)
        rw_mol.RemoveAtom(groups[site_idx])
        return _finalize_product(rw_mol)
    if definition.operation == "mpeg_acyl":
        if (
            set(params) != {"units"}
            or not isinstance(params["units"], int)
            or params["units"] < 1
        ):
            raise ModificationError(
                "mpeg_acetyl requires a positive integer units parameter"
            )
        fragment = "CO" + "CCO" * params["units"] + "CC(=O)*"
        return _graft(mol, site_idx, fragment)
    if definition.operation == "reaction":
        if params:
            raise ModificationError(
                f"Reaction modifier {definition.id} does not accept params"
            )
        return _run_custom_reaction(mol, site_idx, definition.reaction_smarts)
    raise ModificationError(f"Unsupported modifier operation: {definition.operation}")


def apply_recipe(
    smiles: str,
    recipe: Recipe,
    registry: ModifierRegistry,
) -> tuple[str, list[str]]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ModificationError("Invalid input SMILES")
    mol = _mark_original_indices(mol)
    applied = []
    for step in recipe.steps:
        try:
            definition = registry.modifiers[step.modifier]
        except KeyError as exc:
            raise ModificationError(f"Unknown modifier: {step.modifier}") from exc
        mol = apply_modifier(
            mol,
            definition,
            _stable_selector(mol, step.site),
            step.params,
        )
        applied.append(definition.id)
    return Chem.MolToSmiles(mol, isomericSmiles=True), applied


def _unscoped_selector(selector: SiteSelector) -> SiteSelector:
    return SiteSelector(
        role=selector.role,
        atom_index=selector.atom_index,
        smarts=selector.smarts,
        match_index=selector.match_index,
    )


def build_modified_peptide(
    sequence,
    constraint: str,
    residue_registry,
    modifier_registry: ModifierRegistry,
    recipe: Recipe,
) -> tuple[str, list[str]]:
    """Apply residue steps before assembly and global steps to the peptide product."""
    from p2smi.registry import Residue, ResidueRegistry, parse_sequence_tokens
    from p2smi.utilities.smilesgen import constrained_peptide_smiles

    tokens = parse_sequence_tokens(sequence, residue_registry)
    seen_global_step = False
    residue_steps: dict[int, list[ModificationStep]] = {}
    global_steps = []
    for step in recipe.steps:
        if step.site.residue is None:
            seen_global_step = True
            global_steps.append(step)
        else:
            if seen_global_step:
                raise ModificationError(
                    "Residue-scoped recipe steps must precede global or terminal steps"
                )
            if step.site.residue > len(tokens):
                raise ModificationError(
                    f"Recipe residue {step.site.residue} is outside a "
                    f"{len(tokens)}-residue peptide"
                )
            residue_steps.setdefault(step.site.residue - 1, []).append(step)

    positional_residues = []
    positional_tokens = []
    applied = []
    for index, token in enumerate(tokens):
        source = residue_registry.resolve(token)
        mol = _mark_original_indices(Chem.MolFromSmiles(source.smiles))
        for step in residue_steps.get(index, []):
            try:
                definition = modifier_registry.modifiers[step.modifier]
            except KeyError as exc:
                raise ModificationError(f"Unknown modifier: {step.modifier}") from exc
            selector = _stable_selector(mol, _unscoped_selector(step.site))
            mol = apply_modifier(mol, definition, selector, step.params)
            applied.append(f"{definition.id}@{index + 1}")
        positional_id = f"POSITION_{index + 1}"
        positional_tokens.append(positional_id)
        positional_residues.append(
            Residue(
                id=positional_id,
                name=f"{source.name} at position {index + 1}",
                smiles=Chem.MolToSmiles(mol, isomericSmiles=True),
            )
        )

    positional_registry = ResidueRegistry(positional_residues)
    _, _, peptide_smiles = constrained_peptide_smiles(
        ",".join(positional_tokens),
        constraint,
        registry=positional_registry,
    )
    if global_steps:
        global_recipe = Recipe(f"{recipe.id}_global", tuple(global_steps))
        peptide_smiles, global_applied = apply_recipe(
            peptide_smiles, global_recipe, modifier_registry
        )
        applied.extend(global_applied)
    return peptide_smiles, applied
