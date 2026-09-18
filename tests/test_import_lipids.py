import pytest
from rdkit import Chem

from p2smi.import_lipids import fatty_acid_to_acyl_fragment, import_lipid_modifiers
from p2smi.modifiers import (
    ModificationError,
    ModifierRegistry,
    SiteSelector,
    apply_modifier,
)


def test_import_lipid_maps_sdf_and_apply_to_lysine(tmp_path):
    lipid = Chem.MolFromSmiles("CCCCCCCC(=O)O")
    lipid.SetProp("LM_ID", "LMFA01010001")
    lipid.SetProp("NAME", "Octanoic acid")
    sdf_file = tmp_path / "lipids.sdf"
    writer = Chem.SDWriter(str(sdf_file))
    writer.write(lipid)
    writer.close()
    registry_file = tmp_path / "lipid-modifiers.json"

    count = import_lipid_modifiers(sdf_file, registry_file)
    registry = ModifierRegistry.read_json(registry_file)
    definition = registry.modifiers["LMFA01010001"]
    lysine = Chem.MolFromSmiles("N[C@@H](CCCCN)C(=O)O")
    product = apply_modifier(
        lysine,
        definition,
        SiteSelector(role="sidechain_amine"),
    )

    assert count == 1
    assert definition.name == "Octanoic acid"
    assert product.HasSubstructMatch(Chem.MolFromSmarts("N-C(=O)-C"))


def test_lmsd_import_skips_non_acids_unless_strict(tmp_path):
    non_acid = Chem.MolFromSmiles("CCCC")
    non_acid.SetProp("LM_ID", "LMFA00000001")
    non_acid.SetProp("NAME", "Not an acid")
    acid = Chem.MolFromSmiles("CCCCCCCC(=O)O")
    acid.SetProp("LM_ID", "LMFA01010001")
    acid.SetProp("NAME", "Octanoic acid")
    sdf_file = tmp_path / "structures.sdf"
    writer = Chem.SDWriter(str(sdf_file))
    writer.write(non_acid)
    writer.write(acid)
    writer.close()

    registry_file = tmp_path / "lipid-modifiers.json"
    assert import_lipid_modifiers(sdf_file, registry_file) == 1
    assert ModifierRegistry.read_json(registry_file).modifiers["LMFA01010001"].name == "Octanoic acid"

    with pytest.raises(ModificationError, match="one free carboxylic acid"):
        import_lipid_modifiers(sdf_file, tmp_path / "strict.json", strict=True)


def test_lmsd_import_skips_generic_r_group_acids(tmp_path):
    generic = Chem.MolFromSmiles("[*]C(=O)O")
    generic.SetProp("LM_ID", "LMFA01010000")
    generic.SetProp("NAME", "fatty acid")
    acid = Chem.MolFromSmiles("CCCC(=O)O")
    acid.SetProp("LM_ID", "LMFA01010004")
    acid.SetProp("NAME", "Butanoic acid")
    sdf_file = tmp_path / "structures.sdf"
    writer = Chem.SDWriter(str(sdf_file))
    writer.write(generic)
    writer.write(acid)
    writer.close()

    registry_file = tmp_path / "lipid-modifiers.json"
    assert import_lipid_modifiers(sdf_file, registry_file) == 1
    assert set(ModifierRegistry.read_json(registry_file).modifiers) == {"LMFA01010004"}

    with pytest.raises(ModificationError, match="terminal '\\*'"):
        import_lipid_modifiers(sdf_file, tmp_path / "strict.json", strict=True)


def test_lipid_import_rejects_ambiguous_or_non_acid_structures():
    with pytest.raises(ModificationError, match="one free carboxylic acid"):
        fatty_acid_to_acyl_fragment(Chem.MolFromSmiles("CCCC"))

    with pytest.raises(ModificationError, match="not symmetry-equivalent"):
        fatty_acid_to_acyl_fragment(
            Chem.MolFromSmiles("O=C(O)CCCC(C)CCC(=O)O")
        )

    with pytest.raises(ModificationError, match="connected component"):
        fatty_acid_to_acyl_fragment(Chem.MolFromSmiles("CCCC(=O)O.[Na+]"))


def test_lipid_import_accepts_symmetric_diacid_and_retains_distal_acid():
    fragment = fatty_acid_to_acyl_fragment(
        Chem.MolFromSmiles("O=C(O)CCCCCCCCCCCCCCCC(=O)O")
    )
    mol = Chem.MolFromSmiles(fragment)

    assert len(mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[OH]"))) == 1
    assert len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]) == 1


def test_lipid_import_filters_carbon_count_and_acid_mode(tmp_path):
    monoacid = Chem.MolFromSmiles("CCCCCCCC(=O)O")
    monoacid.SetProp("LM_ID", "C8_MONO")
    diacid = Chem.MolFromSmiles("O=C(O)CCCCCCCCCCCCCCCCCC(=O)O")
    diacid.SetProp("LM_ID", "C20_DIACID")
    too_long = Chem.MolFromSmiles("C" * 21 + "(=O)O")
    too_long.SetProp("LM_ID", "C21_MONO")
    sdf_file = tmp_path / "lipids.sdf"
    writer = Chem.SDWriter(str(sdf_file))
    for mol in (monoacid, diacid, too_long):
        writer.write(mol)
    writer.close()

    mono_file = tmp_path / "mono.json"
    assert import_lipid_modifiers(sdf_file, mono_file, acid_mode="monoacid") == 1
    assert set(ModifierRegistry.read_json(mono_file).modifiers) == {"C8_MONO"}

    diacid_file = tmp_path / "diacid.json"
    assert import_lipid_modifiers(sdf_file, diacid_file, acid_mode="diacid") == 1
    assert set(ModifierRegistry.read_json(diacid_file).modifiers) == {"C20_DIACID"}