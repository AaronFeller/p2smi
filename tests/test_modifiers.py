import json

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from p2smi.modifiers import (
    ModificationError,
    ModificationStep,
    ModifierDefinition,
    ModifierRegistry,
    Recipe,
    SiteSelector,
    apply_modifier,
    apply_recipe,
    build_modified_peptide,
    builtin_modifier_registry,
    resolve_site,
)
from p2smi.registry import Residue, ResidueRegistry


PEPTIDE = "NCC(=O)N[C@@H](C)C(=O)O"
LYSINE = "N[C@@H](CCCCN)C(=O)O"


def _apply(modifier_id, smiles, selector, params=None):
    registry = builtin_modifier_registry()
    return apply_modifier(
        Chem.MolFromSmiles(smiles),
        registry.modifiers[modifier_id],
        SiteSelector.from_dict(selector),
        params,
    )


def _smiles(mol):
    assert len(Chem.GetMolFrags(mol)) == 1
    assert not any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms())
    assert not any(atom.GetAtomMapNum() for atom in mol.GetAtoms())
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def test_terminal_caps_and_backbone_n_methylation_have_expected_products():
    acetylated = _apply("n_acetyl", PEPTIDE, {"terminus": "N"})
    amidated = _apply("c_amide", PEPTIDE, {"terminus": "C"})
    methylated = _apply("n_methyl", PEPTIDE, {"role": "backbone_amide_n"})

    assert acetylated.HasSubstructMatch(Chem.MolFromSmarts("CC(=O)NCC(=O)"))
    assert amidated.HasSubstructMatch(Chem.MolFromSmarts("C(=O)N"))
    assert methylated.HasSubstructMatch(Chem.MolFromSmarts("C(=O)N(C)[C@@H]"))


def test_c_terminal_amidation_does_not_consume_glutamate_sidechain_acid():
    peptide = "N[C@@H](CCC(=O)O)C(=O)N[C@@H](C)C(=O)O"
    product = _apply("c_amide", peptide, {"terminus": "C"})

    assert len(product.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[OH]"))) == 1
    assert product.HasSubstructMatch(Chem.MolFromSmarts("[C@@H](C)C(=O)N"))


@pytest.mark.parametrize(
    ("modifier_id", "minimum_carbons"),
    [
        ("octanoyl", 8),
        ("decanoyl", 10),
        ("dodecanoyl", 12),
        ("tetradecanoyl", 14),
        ("palmitoyl", 16),
        ("stearoyl", 18),
        ("eicosanoyl", 20),
    ],
)
def test_lipidation_builds_defined_lysine_amides(modifier_id, minimum_carbons):
    product = _apply(modifier_id, LYSINE, {"role": "sidechain_amine"})

    assert _smiles(product)
    assert (
        sum(atom.GetAtomicNum() == 6 for atom in product.GetAtoms())
        >= minimum_carbons + 6
    )
    assert product.HasSubstructMatch(Chem.MolFromSmarts("N-C(=O)-C"))


@pytest.mark.parametrize(
    "modifier_id",
    ["hexadecanedioyl", "octadecanedioyl", "icosanedioyl"],
)
def test_diacid_lipidation_retains_distal_acid(modifier_id):
    product = _apply(modifier_id, LYSINE, {"role": "sidechain_amine"})

    assert _smiles(product)
    assert len(product.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[OH]"))) == 2
    assert product.HasSubstructMatch(Chem.MolFromSmarts("N-C(=O)-C"))


@pytest.mark.parametrize(
    "linker",
    ["aeea", "aeea2", "gamma_glu", "gamma_glu_aeea", "gamma_glu_aeea2"],
)
def test_lipid_linkers_form_defined_amides(linker):
    product = _apply(
        "palmitoyl",
        LYSINE,
        {"role": "sidechain_amine"},
        {"linker": linker},
    )

    assert _smiles(product)
    assert product.HasSubstructMatch(Chem.MolFromSmarts("N-C(=O)"))
    assert product.HasSubstructMatch(Chem.MolFromSmarts("N-C(=O)-[CH2]"))


def test_gamma_glu_aeea2_c18_diacid_has_expected_linker_chemistry():
    product = _apply(
        "octadecanedioyl",
        LYSINE,
        {"role": "sidechain_amine"},
        {"linker": "gamma_glu_aeea2"},
    )

    assert _smiles(product)
    assert len(product.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[OH]"))) == 3
    ether_pattern = Chem.MolFromSmarts("[O;X2;H0]([#6])[#6]")
    assert len(product.GetSubstructMatches(ether_pattern)) == 4
    assert [label for _, label in Chem.FindMolChiralCenters(product)] == ["S", "S"]


def test_lipid_linker_rejects_non_acyl_grafts_and_unknown_linkers():
    with pytest.raises(ModificationError, match="acyl modifier"):
        _apply(
            "phosphoryl",
            "N[C@@H](CO)C(=O)O",
            {"role": "sidechain_hydroxyl"},
            {"linker": "aeea"},
        )

    with pytest.raises(ModificationError, match="Unknown lipid linker"):
        _apply(
            "palmitoyl",
            LYSINE,
            {"role": "sidechain_amine"},
            {"linker": "unknown"},
        )


def test_mpeg_requires_units_and_adds_exact_ether_count():
    with pytest.raises(ModificationError, match="units"):
        _apply("mpeg_acetyl", LYSINE, {"role": "sidechain_amine"})

    product = _apply(
        "mpeg_acetyl",
        LYSINE,
        {"role": "sidechain_amine"},
        {"units": 3},
    )
    assert _smiles(product)
    peg_chain = Chem.MolFromSmarts("COCCOCCOCCOCC(=O)N")
    assert product.HasSubstructMatch(peg_chain)


def test_thiol_phosphorylation_and_click_handle_products():
    carbamidomethyl = _apply(
        "carbamidomethyl", "N[C@@H](CS)C(=O)O", {"role": "sidechain_thiol"}
    )
    phosphorylated = _apply(
        "phosphoryl", "N[C@@H](CO)C(=O)O", {"role": "sidechain_hydroxyl"}
    )
    azide = _apply("azidoacetyl", LYSINE, {"role": "sidechain_amine"})
    alkyne = _apply("propiolyl", LYSINE, {"role": "sidechain_amine"})

    assert carbamidomethyl.HasSubstructMatch(Chem.MolFromSmarts("SCC(=O)N"))
    assert phosphorylated.HasSubstructMatch(Chem.MolFromSmarts("COP(=O)(O)O"))
    assert azide.HasSubstructMatch(Chem.MolFromSmarts("[N-]=[N+]=N"))
    assert alkyne.HasSubstructMatch(Chem.MolFromSmarts("C#CC(=O)N"))


def test_selectors_fail_closed_on_ambiguity_and_support_explicit_escape_hatches():
    dilysine = Chem.MolFromSmiles("N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)O")

    with pytest.raises(ModificationError, match="2 atoms"):
        resolve_site(dilysine, SiteSelector(role="sidechain_amine"))

    amines = sorted(
        match[0] for match in dilysine.GetSubstructMatches(Chem.MolFromSmarts("[N;H2]"))
    )
    index, role = resolve_site(dilysine, SiteSelector(atom_index=amines[-1]))
    assert index == amines[-1]
    assert role == "atom_index"

    index, role = resolve_site(
        dilysine,
        SiteSelector(smarts="[N;H2]", match_index=1),
    )
    assert index == amines[1]
    assert role == "smarts"

    with pytest.raises(ModificationError, match="cannot target atom"):
        _apply("n_acetyl", "CC", {"atom_index": 0})

    with pytest.raises(ModificationError, match="cannot target atom"):
        _apply("n_acetyl", "CCN", {"smarts": "[C]", "match_index": 0})


def test_atom_indices_refer_to_original_recipe_graph_after_atom_removal():
    registry = builtin_modifier_registry()
    source = Chem.MolFromSmiles(PEPTIDE)
    terminal_oxygen = max(
        atom.GetIdx()
        for atom in source.GetAtoms()
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() > 0
    )
    recipe = Recipe(
        "stable_indices",
        (
            ModificationStep("c_amide", SiteSelector(terminus="C")),
            ModificationStep("n_acetyl", SiteSelector(atom_index=terminal_oxygen)),
        ),
    )

    with pytest.raises(ModificationError, match="no longer available"):
        apply_recipe(PEPTIDE, recipe, registry)


def test_recipe_is_deterministic_and_input_is_unchanged_on_failure():
    registry = builtin_modifier_registry()
    recipe = Recipe(
        "caps",
        (
            ModificationStep("n_acetyl", SiteSelector(terminus="N")),
            ModificationStep("c_amide", SiteSelector(terminus="C")),
        ),
    )
    product, applied = apply_recipe(PEPTIDE, recipe, registry)

    assert Chem.MolFromSmiles(product) is not None
    assert applied == ["n_acetyl", "c_amide"]
    assert Chem.MolToSmiles(Chem.MolFromSmiles(PEPTIDE), isomericSmiles=True) != product

    failing = Recipe(
        "bad",
        (ModificationStep("n_acetyl", SiteSelector(role="sidechain_thiol")),),
    )
    with pytest.raises(ModificationError):
        apply_recipe(PEPTIDE, failing, registry)
    assert Chem.MolFromSmiles(PEPTIDE) is not None


def test_custom_modifier_registry_round_trip_and_validation(tmp_path):
    document = {
        "schema_version": 1,
        "kind": "p2smi_modifier_registry",
        "modifiers": [
            {
                "id": "fluoroacetyl",
                "name": "Fluoroacetylation",
                "operation": "graft",
                "allowed_roles": ["sidechain_amine"],
                "fragment_smiles": "FCC(=O)*",
            }
        ],
        "recipes": [
            {
                "id": "fluoro_lys",
                "steps": [
                    {
                        "modifier": "fluoroacetyl",
                        "site": {"role": "sidechain_amine"},
                    }
                ],
            }
        ],
    }
    registry_file = tmp_path / "modifiers.json"
    registry_file.write_text(json.dumps(document))
    registry = builtin_modifier_registry().merge(
        ModifierRegistry.read_json(registry_file)
    )

    product, applied = apply_recipe(LYSINE, registry.recipes["fluoro_lys"], registry)
    assert applied == ["fluoroacetyl"]
    assert Chem.MolFromSmiles(product).HasSubstructMatch(Chem.MolFromSmarts("FCC(=O)N"))

    document["modifiers"][0]["fragment_smiles"] = "CC"
    with pytest.raises(ModificationError, match="attachment atom"):
        ModifierRegistry.from_document(document)


def test_custom_reaction_smarts_targets_only_the_selected_atom():
    document = {
        "schema_version": 1,
        "kind": "p2smi_modifier_registry",
        "modifiers": [
            {
                "id": "custom_n_methyl",
                "name": "Custom amine methylation",
                "operation": "reaction",
                "allowed_roles": ["sidechain_amine"],
                "reaction_smarts": "[N;H1,H2:1]>>[N:1]C",
            }
        ],
        "recipes": [],
    }
    registry = ModifierRegistry.from_document(document)
    definition = registry.modifiers["custom_n_methyl"]

    product = apply_modifier(
        Chem.MolFromSmiles(LYSINE),
        definition,
        SiteSelector(role="sidechain_amine"),
    )

    assert product.HasSubstructMatch(Chem.MolFromSmarts("[N;H2][C@@H]"))
    assert product.HasSubstructMatch(Chem.MolFromSmarts("[CH2][N;H1][CH3]"))
    assert definition.to_dict() == document["modifiers"][0]


@pytest.mark.parametrize(
    ("reaction_smarts", "message"),
    [
        ("[N:1].[C:2]>>[N:1][C:2]", "one reactant"),
        ("[N:1]>>NC", "atom map :1"),
        ("[N:1]>>[N:1][C:2]", "absent from the reactant"),
        ("not reaction SMARTS", "invalid reaction SMARTS"),
    ],
)
def test_custom_reaction_smarts_validation(reaction_smarts, message):
    definition = {
        "id": "invalid_reaction",
        "name": "Invalid reaction",
        "operation": "reaction",
        "allowed_roles": ["sidechain_amine"],
        "reaction_smarts": reaction_smarts,
    }

    with pytest.raises(ModificationError, match=message):
        ModifierDefinition.from_dict(definition)


def test_custom_reaction_fails_when_selected_site_does_not_match_template():
    definition = ModifierDefinition.from_dict(
        {
            "id": "o_methyl",
            "name": "Hydroxyl methylation",
            "operation": "reaction",
            "allowed_roles": ["sidechain_amine"],
            "reaction_smarts": "[O;H1:1]>>[O:1]C",
        }
    )

    with pytest.raises(ModificationError, match="no valid product"):
        apply_modifier(
            Chem.MolFromSmiles(LYSINE),
            definition,
            SiteSelector(role="sidechain_amine"),
        )


def test_custom_reaction_preserves_original_indices_for_later_steps():
    reaction = ModifierDefinition.from_dict(
        {
            "id": "custom_n_methyl",
            "name": "Custom amine methylation",
            "operation": "reaction",
            "allowed_roles": ["sidechain_amine"],
            "reaction_smarts": "[N;H1,H2:1]>>[N:1]C",
        }
    )
    registry = builtin_modifier_registry().merge(ModifierRegistry([reaction]))
    source = Chem.MolFromSmiles(LYSINE)
    n_terminal_index = resolve_site(source, SiteSelector(terminus="N"))[0]
    recipe = Recipe(
        "reaction_then_cap",
        (
            ModificationStep(
                "custom_n_methyl",
                SiteSelector(role="sidechain_amine"),
            ),
            ModificationStep(
                "n_acetyl",
                SiteSelector(atom_index=n_terminal_index),
            ),
        ),
    )

    product, applied = apply_recipe(LYSINE, recipe, registry)

    assert Chem.MolFromSmiles(product).HasSubstructMatch(
        Chem.MolFromSmarts("CC(=O)N[C@@H]")
    )
    assert applied == ["custom_n_methyl", "n_acetyl"]


def test_lipid_linker_recipe_round_trips_through_registry_json():
    document = {
        "schema_version": 1,
        "kind": "p2smi_modifier_registry",
        "modifiers": [],
        "recipes": [
            {
                "id": "linked_lipid",
                "steps": [
                    {
                        "modifier": "octadecanedioyl",
                        "site": {"role": "sidechain_amine"},
                        "params": {"linker": "gamma_glu_aeea2"},
                    }
                ],
            }
        ],
    }
    registry = ModifierRegistry.from_document(
        document,
        base=builtin_modifier_registry(),
    )

    product, applied = apply_recipe(LYSINE, registry.recipes["linked_lipid"], registry)

    assert Chem.MolFromSmiles(product) is not None
    assert applied == ["octadecanedioyl"]
    assert registry.to_document()["recipes"][0] == document["recipes"][0]


def test_modifications_preserve_stereochemistry_and_expected_formula_delta():
    source = Chem.MolFromSmiles(LYSINE)
    product = _apply("n_acetyl", LYSINE, {"terminus": "N"})

    source_chiral = Chem.FindMolChiralCenters(source, includeUnassigned=True)
    product_chiral = Chem.FindMolChiralCenters(product, includeUnassigned=True)
    assert [label for _, label in source_chiral] == [
        label for _, label in product_chiral
    ]
    assert rdMolDescriptors.CalcMolFormula(source) != rdMolDescriptors.CalcMolFormula(
        product
    )


def test_sequence_recipe_targets_one_residue_before_peptide_assembly():
    residues = ResidueRegistry(
        [
            Residue(id="LYS", name="Lysine", smiles=LYSINE),
            Residue(id="ALA", name="Alanine", smiles="N[C@@H](C)C(=O)O"),
        ]
    )
    modifiers = builtin_modifier_registry()
    recipe = Recipe(
        "lipidated",
        (
            ModificationStep(
                "palmitoyl",
                SiteSelector(residue=2, role="sidechain_amine"),
            ),
            ModificationStep(
                "n_methyl",
                SiteSelector(residue=3, role="backbone_n"),
            ),
        ),
    )

    product, applied = build_modified_peptide(
        "ALA,LYS,ALA", "", residues, modifiers, recipe
    )

    mol = Chem.MolFromSmiles(product)
    assert mol is not None
    assert mol.HasSubstructMatch(Chem.MolFromSmarts("N-C(=O)-[CH2]-[CH2]"))
    assert mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)N(C)[C@@H]"))
    assert applied == ["palmitoyl@2", "n_methyl@3"]


def test_sequence_recipe_supports_residue_scoped_smarts_match_selection():
    residues = ResidueRegistry(
        [
            Residue(id="DIOL", name="Diol", smiles="N[C@@H](C(O)CO)C(=O)O"),
            Residue(id="ALA", name="Alanine", smiles="N[C@@H](C)C(=O)O"),
        ]
    )
    modifiers = builtin_modifier_registry()
    recipe = Recipe(
        "site_specific",
        (
            ModificationStep(
                "phosphoryl",
                SiteSelector(residue=1, smarts="[O;H1]", match_index=1),
            ),
        ),
    )

    product, applied = build_modified_peptide(
        "DIOL,ALA", "", residues, modifiers, recipe
    )

    assert Chem.MolFromSmiles(product).HasSubstructMatch(
        Chem.MolFromSmarts("OP(=O)(O)O")
    )
    assert applied == ["phosphoryl@1"]


def test_sequence_recipe_rejects_out_of_range_and_consumed_cyclization_site():
    residues = ResidueRegistry(
        [
            Residue(id="CYS", name="Cysteine", smiles="N[C@@H](CS)C(=O)O"),
            Residue(id="ALA", name="Alanine", smiles="N[C@@H](C)C(=O)O"),
        ]
    )
    modifiers = builtin_modifier_registry()

    out_of_range = Recipe(
        "outside",
        (
            ModificationStep(
                "carbamidomethyl",
                SiteSelector(residue=5, role="sidechain_thiol"),
            ),
        ),
    )
    with pytest.raises(ModificationError, match="outside"):
        build_modified_peptide("CYS,ALA", "", residues, modifiers, out_of_range)

    conflict = Recipe(
        "conflict",
        (
            ModificationStep(
                "carbamidomethyl",
                SiteSelector(residue=1, role="sidechain_thiol"),
            ),
        ),
    )
    with pytest.raises(Exception, match="reaction sites"):
        build_modified_peptide(
            "CYS,ALA,ALA,CYS", "SSCXXC", residues, modifiers, conflict
        )


def test_sequence_recipe_rejects_residue_step_after_global_step():
    residues = ResidueRegistry(
        [Residue(id="ALA", name="Alanine", smiles="N[C@@H](C)C(=O)O")]
    )
    recipe = Recipe(
        "wrong_order",
        (
            ModificationStep("c_amide", SiteSelector(terminus="C")),
            ModificationStep("n_methyl", SiteSelector(residue=1, role="backbone_n")),
        ),
    )

    with pytest.raises(ModificationError, match="must precede"):
        build_modified_peptide("ALA", "", residues, builtin_modifier_registry(), recipe)
