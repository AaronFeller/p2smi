from rdkit import Chem

from p2smi.utilities.aminoacids import all_aminos


def test_amino_acid_library_entries_are_structurally_consistent():
    required_keys = {
        "Code",
        "Formula",
        "Letter",
        "MolWeight",
        "SMILES",
        "cterm",
        "disulphide",
        "ester",
        "nterm",
    }

    letters = []
    codes = []
    for name, props in all_aminos.items():
        assert required_keys.issubset(props), name

        base = Chem.MolFromSmiles(props["SMILES"])
        assert base is not None, name
        assert sum(atom.GetAtomicNum() == 0 for atom in base.GetAtoms()) == 0, name

        letters.append(props["Letter"])
        codes.append(props["Code"])

        for key in ["cterm", "disulphide", "ester", "nterm"]:
            value = props[key]
            if value is False:
                continue
            mol = Chem.MolFromSmiles(value)
            assert mol is not None, f"{name}:{key}"
            assert (
                sum(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()) == 1
            ), f"{name}:{key}"

    assert len(letters) == len(set(letters))
    assert len(codes) == len(set(codes))
