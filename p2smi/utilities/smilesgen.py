# Standard library imports
import itertools
import operator
import os
import os.path as path
import sys

# RDKit imports for chemical structure handling
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

# Import amino acid definitions
from p2smi.registry import parse_sequence_tokens
from p2smi.utilities.aminoacids import all_aminos

from functools import lru_cache

# build direct reverse maps once
LETTER2NAME = {props["Letter"]: name for name, props in all_aminos.items()}
CODE2NAME = {
    props.get("Code"): name for name, props in all_aminos.items() if "Code" in props
}

aminodata = all_aminos  # Current dictionary of amino acids

# Custom exceptions for specific error conditions


class CustomError(Exception):
    pass


class NoCysteineError(CustomError):
    pass


class BondSpecError(CustomError):
    pass


class FormatError(CustomError):
    pass


class UndefinedAminoError(CustomError):
    pass


class UndefinedPropertyError(CustomError):
    pass


class SmilesError(CustomError):
    pass


PEPTIDE_BOND_REACTION = AllChem.ReactionFromSmarts(
    "[N:5][C:4][C:1](=[O:2])[OH:3]."
    "[N;H1,H2:6][C:7][C:8](=[O:9])[OH:10]>>"
    "[N:5][C:4][C:1](=[O:2])[N:6][C:7][C:8](=[O:9])[OH:10]"
)


def add_amino(name):
    # Add an amino acid to aminodata if it exists in all_aminos and isn't already included
    if name in all_aminos and name not in aminodata:
        aminodata[name] = all_aminos[name]
        return True
    else:
        raise UndefinedAminoError(f"{name} not recognised as valid amino acid")


def remove_amino(name):
    # Remove an amino acid from aminodata if it exists
    if name in aminodata:
        del aminodata[name]
    else:
        raise UndefinedAminoError(f"{name} not found in amino acids")


def print_possible_aminos():
    # Return a list of all possible amino acid names
    return list(all_aminos.keys())


def print_included_aminos():
    # Return a list of currently included amino acid names
    return list(aminodata.keys())


def return_available_residues(out="Letter"):
    # Return a list of available residue properties (default: 'Letter')
    return [properties[out] for properties in aminodata.values()]


# precompute sets of residue *names* that satisfy each constraint
CONSTRAINT_RES_NAME_SETS = {
    key: frozenset([name for name, props in aminodata.items() if props.get(key)])
    for key in {"disulphide", "cterm", "nterm", "ester"}
}


def return_constraint_resis(constraint_type):
    # fast, no rebuild on each call
    return list(CONSTRAINT_RES_NAME_SETS[constraint_type])


@lru_cache(maxsize=None)
def property_to_name(prop, value):
    if prop == "Letter":
        try:
            return LETTER2NAME[value]
        except KeyError:
            raise UndefinedAminoError(f"{value} not found")
    if prop == "Code":
        try:
            return CODE2NAME[value]
        except KeyError:
            raise UndefinedAminoError(f"{value} not found")
    # fallback to old path for rare props
    for name, properties in aminodata.items():
        if properties.get(prop) == value:
            return name
    raise UndefinedAminoError(f"Amino-acid {value} for {prop} not found")


def gen_all_pos_peptides(pepliblen):
    # Generate all possible peptide sequences of a given length
    amino_keys = list(aminodata.keys())
    for pep in itertools.product(amino_keys, repeat=pepliblen):
        yield pep


def gen_all_matching_peptides(pattern):
    # Generate all peptide sequences matching a given pattern,
    # where "X" (or "x") is treated as a wildcard for any amino acid.
    pattern = (
        pattern.replace("x", "X")
        if isinstance(pattern, str)
        else ["X" if resi == "x" else resi for resi in pattern]
    )
    amino_keys = list(aminodata.keys())
    for pep in itertools.product(amino_keys, repeat=pattern.count("X")):
        pep = list(pep)
        outpep = []
        for resi in pattern:
            if resi != "X":
                # If residue is not a wildcard, use it directly or convert
                if resi in aminodata:
                    outpep.append(resi)
                else:
                    outpep.append(property_to_name("Letter", resi))
            else:
                outpep.append(pep.pop(0))
        yield outpep


_CONSTRAINT_LETTER_SETS = {
    "disulphide": frozenset(
        props["Letter"] for _, props in aminodata.items() if props.get("disulphide")
    ),
    "cterm": frozenset(
        props["Letter"] for _, props in aminodata.items() if props.get("cterm")
    ),
    "nterm": frozenset(
        props["Letter"] for _, props in aminodata.items() if props.get("nterm")
    ),
    "ester": frozenset(
        props["Letter"] for _, props in aminodata.items() if props.get("ester")
    ),
}


@lru_cache(maxsize=None)
def _is_valid_letter(letter: str) -> bool:
    # uses LETTER2NAME you already define elsewhere
    return letter in LETTER2NAME


def _normalize_seq_letters(seq):
    """Return the sequence as a list of one-letter codes; validate quickly."""
    letters = list(seq) if isinstance(seq, str) else list(seq)
    for r in letters:
        if not _is_valid_letter(r):
            raise UndefinedAminoError(f"{r} not recognised as amino acid letter")
    return letters


def _preserve_seq_type(orig, letters_list):
    """Return letters with the same container type convention your code expects."""
    if isinstance(orig, tuple):
        return tuple(letters_list)
    # For strings, downstream code often uses ','.join(seq), so keep as list
    return list(letters_list)


# -------------------------------------------------------------------


def _registry_sequence_sites(peptideseq, registry):
    tokens = parse_sequence_tokens(peptideseq, registry)
    sites = []
    for token in tokens:
        residue = registry.resolve(token)
        mol = _mol_from_smiles(residue.smiles, context=f"residue {residue.id}")
        candidates = _reactive_site_indices(mol)
        sites.append(
            {kind for kind, indices in candidates.items() if len(indices) == 1}
        )
    return tokens, sites


def can_ssbond(peptideseq, registry=None):
    """Disulphide: need at least two Cys-like residues;
    pick the pair with max separation (>=3 apart)."""
    if registry is not None:
        letters, sites = _registry_sequence_sites(peptideseq, registry)
        locs = [i for i, kinds in enumerate(sites) if "disulphide" in kinds]
    else:
        letters = _normalize_seq_letters(peptideseq)
        dis = _CONSTRAINT_LETTER_SETS["disulphide"]
        locs = [i for i, r in enumerate(letters) if r in dis]
    if len(locs) < 2:
        return False
    (a, b), sep = max(
        ((p, abs(p[0] - p[1])) for p in itertools.combinations(locs, 2)),
        key=operator.itemgetter(1),
    )
    if sep <= 2:
        return False
    pattern = "SS" + "".join("C" if i in (a, b) else "X" for i in range(len(letters)))
    return _preserve_seq_type(peptideseq, letters), pattern


def can_htbond(peptideseq, registry=None):
    """Your original heuristic: qualifies if len >= 5 or exactly 2."""
    letters = (
        parse_sequence_tokens(peptideseq, registry)
        if registry is not None
        else _normalize_seq_letters(peptideseq)
    )
    if len(letters) >= 5 or len(letters) == 2:
        return _preserve_seq_type(peptideseq, letters), "HT"
    return False


def can_scntbond(peptideseq, strict=False, registry=None):
    """Sidechain to N-terminus using a cterm-capable sidechain."""
    if registry is not None:
        letters, sites = _registry_sequence_sites(peptideseq, registry)
        locs = [i for i, kinds in enumerate(sites[3:], start=3) if "cterm" in kinds]
    else:
        letters = _normalize_seq_letters(peptideseq)
        cterm = _CONSTRAINT_LETTER_SETS["cterm"]
        locs = [i for i, r in enumerate(letters[3:], start=3) if r in cterm]
    if not locs or (len(locs) > 1 and strict):
        return False

    idx = locs[-1]  # keep your previous "last occurrence" behavior
    pattern = ["SC"] + ["Z" if i == idx else "X" for i in range(len(letters))]
    return _preserve_seq_type(peptideseq, letters), "".join(pattern)


def can_scctbond(peptideseq, strict=False, registry=None):
    """Sidechain to C-terminus using an nterm or ester-capable sidechain."""
    if registry is not None:
        letters, sites = _registry_sequence_sites(peptideseq, registry)
        locs = [(i, "N") for i, kinds in enumerate(sites[:-3]) if "nterm" in kinds]
        locs += [(i, "E") for i, kinds in enumerate(sites[:-3]) if "ester" in kinds]
    else:
        letters = _normalize_seq_letters(peptideseq)
        esters = _CONSTRAINT_LETTER_SETS["ester"]
        nterms = _CONSTRAINT_LETTER_SETS["nterm"]
        locs = [(i, "N") for i, r in enumerate(letters[:-3]) if r in nterms]
        locs += [(i, "E") for i, r in enumerate(letters[:-3]) if r in esters]
    if not locs or (len(locs) > 1 and strict):
        return False

    idx, code = locs[0]  # keep your previous "first eligible" behavior
    pattern = ["SC"] + [code if i == idx else "X" for i in range(len(letters))]
    return _preserve_seq_type(peptideseq, letters), "".join(pattern)


def can_scscbond(peptideseq, strict=False, registry=None):
    """Sidechain-to-sidechain: choose (cterm_pos, partner_pos) with max separation >= 2.
    Encode 'Z' at cterm_pos and 'N'/'E' at partner_pos depending on site set.
    """
    if registry is not None:
        letters, sites = _registry_sequence_sites(peptideseq, registry)
        locs_n = [i for i, kinds in enumerate(sites) if "nterm" in kinds]
        locs_c = [i for i, kinds in enumerate(sites) if "cterm" in kinds]
        locs_e = [i for i, kinds in enumerate(sites) if "ester" in kinds]
    else:
        letters = _normalize_seq_letters(peptideseq)
        nterms = _CONSTRAINT_LETTER_SETS["nterm"]
        cterms = _CONSTRAINT_LETTER_SETS["cterm"]
        esters = _CONSTRAINT_LETTER_SETS["ester"]
        locs_n = [i for i, r in enumerate(letters) if r in nterms]
        locs_c = [i for i, r in enumerate(letters) if r in cterms]
        locs_e = [i for i, r in enumerate(letters) if r in esters]

    if not locs_c or not (locs_n or locs_e):
        return False

    partners = [(j, "N") for j in locs_n] + [(j, "E") for j in locs_e]
    pairs = [
        ((i, j, code), abs(i - j))
        for i in locs_c
        for (j, code) in partners
        if abs(i - j) >= 2
    ]
    if not pairs:
        return False

    (ci, pj, code), _ = max(pairs, key=operator.itemgetter(1))
    pattern = "SC" + "".join(
        "Z" if k == ci else (code if k == pj else "X") for k in range(len(letters))
    )
    return _preserve_seq_type(peptideseq, letters), pattern


def what_constraints(peptideseq, registry=None):
    return [
        res
        for res in (
            can_ssbond(peptideseq, registry=registry),
            can_htbond(peptideseq, registry=registry),
            can_scctbond(peptideseq, registry=registry),
            can_scntbond(peptideseq, registry=registry),
            can_scscbond(peptideseq, registry=registry),
        )
        if res
    ]


def aaletter2aaname(aaletter):
    # Convert an amino acid letter to its full name
    for name, properties in all_aminos.items():
        if properties["Letter"] == aaletter:
            return name


def gen_library_strings(
    liblen,
    ssbond=False,
    htbond=False,
    scctbond=False,
    scntbond=False,
    scscbond=False,
    linear=False,
):
    # Generate a library of peptide strings based on specified bond constraints
    filterfuncs = []
    if ssbond:
        filterfuncs.append(can_ssbond)
    if htbond:
        filterfuncs.append(can_htbond)
    if scctbond:
        filterfuncs.append(can_scctbond)
    if scntbond:
        filterfuncs.append(can_scntbond)
    if scscbond:
        filterfuncs.append(can_scscbond)
    for sequence in gen_all_pos_peptides(liblen):
        for func in filterfuncs:
            if trialpeptide := func(sequence):
                yield trialpeptide
    if linear:
        for peptide in gen_all_pos_peptides(liblen):
            yield (peptide, "")


def gen_library_from_file(filepath, ignore_errors=False):
    # Generate peptide library entries from a file, ignoring commented/empty lines.
    with open(filepath) as peptides:
        for line in peptides:
            if line.startswith("#") or not line.strip():
                continue
            try:
                sequence, bond_def = map(str.strip, line.split(";"))
                if len(sequence.split(",")) == 1 and sequence not in all_aminos:
                    # Convert single-letter sequence to full names
                    sequence = [aaletter2aaname(letter) for letter in sequence]
                else:
                    sequence = sequence.split(",")
                yield constrained_peptide_smiles(sequence, bond_def)
            except Exception:
                if ignore_errors:
                    yield (None, None, None)
                else:
                    raise


def nmethylate_peptide_smiles(smiles):
    # N-methylate a peptide SMILES structure using substructure replacement
    mol = Chem.MolFromSmiles(smiles)
    n_pattern = Chem.MolFromSmarts("[$([Nh1](C)C(=O)),$([NH2]CC=O)]")
    methylated_pattern = Chem.MolFromSmarts("N(C)")
    rmol = AllChem.ReplaceSubstructs(
        mol, n_pattern, methylated_pattern, replaceAll=True
    )
    return Chem.MolToSmiles(rmol[0], isomericSmiles=True)


def nmethylate_peptides(structs):
    # Apply N-methylation to a sequence of peptide structures
    for struct in structs:
        seq, bond_def, smiles = struct
        if smiles:
            yield seq, bond_def, nmethylate_peptide_smiles(smiles)


@lru_cache(maxsize=None)
def return_smiles(resi):
    return return_constrained_smiles(resi, "SMILES")


@lru_cache(maxsize=None)
def return_constrained_smiles(resi, constraint):
    try:
        return aminodata[resi][constraint]
    except KeyError:
        try:
            return aminodata[property_to_name("Letter", resi)][constraint]
        except UndefinedAminoError:
            return aminodata[property_to_name("Code", resi)][constraint]


def _mol_from_smiles(smiles, *, context):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SmilesError(f"Could not parse {context} SMILES: {smiles}")
    return mol


def _run_peptide_bond_reaction(left_mol, right_mol):
    products = PEPTIDE_BOND_REACTION.RunReactants((left_mol, right_mol))
    if not products:
        raise SmilesError("SMARTS peptide-bond reaction did not produce a product")

    for product_set in products:
        try:
            product = Chem.Mol(product_set[0])
            Chem.SanitizeMol(product)
            return product
        except Exception:
            continue

    raise SmilesError("SMARTS peptide-bond reaction produced only invalid products")


def _has_double_bonded_oxygen(atom):
    return any(
        bond.GetBondType() == Chem.BondType.DOUBLE
        and bond.GetOtherAtom(atom).GetAtomicNum() == 8
        for bond in atom.GetBonds()
    )


def _find_terminal_carboxyl(mol):
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue

        double_oxygens = []
        hydroxyl_oxygens = []
        carbon_neighbors = []
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            if neighbor.GetAtomicNum() == 8:
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    double_oxygens.append(neighbor.GetIdx())
                elif (
                    bond.GetBondType() == Chem.BondType.SINGLE
                    and neighbor.GetDegree() == 1
                ):
                    hydroxyl_oxygens.append(neighbor.GetIdx())
            elif neighbor.GetAtomicNum() == 6:
                carbon_neighbors.append(neighbor)

        if len(double_oxygens) != 1 or len(hydroxyl_oxygens) != 1:
            continue

        if any(
            any(
                nn.GetAtomicNum() == 7
                for nn in carbon.GetNeighbors()
                if nn.GetIdx() != atom.GetIdx()
            )
            for carbon in carbon_neighbors
        ):
            matches.append((atom.GetIdx(), hydroxyl_oxygens[0]))

    if len(matches) != 1:
        raise SmilesError("Could not uniquely identify the peptide C-terminus")
    return matches[0]


def _find_n_terminal_amine(mol):
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue

        if any(
            neighbor.GetAtomicNum() == 6 and _has_double_bonded_oxygen(neighbor)
            for neighbor in atom.GetNeighbors()
        ):
            continue

        for carbon in atom.GetNeighbors():
            if carbon.GetAtomicNum() != 6:
                continue
            if any(
                neighbor.GetAtomicNum() == 6 and _has_double_bonded_oxygen(neighbor)
                for neighbor in carbon.GetNeighbors()
                if neighbor.GetIdx() != atom.GetIdx()
            ):
                matches.append(atom.GetIdx())
                break

    if len(matches) != 1:
        raise SmilesError("Could not uniquely identify the peptide N-terminus")
    return matches[0]


def _dummy_indices(mol):
    return [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]


def _dummy_details(mol, dummy_idx):
    atom = mol.GetAtomWithIdx(dummy_idx)
    neighbors = []
    oxygen_neighbor = None
    for bond in atom.GetBonds():
        neighbor = bond.GetOtherAtom(atom)
        if bond.GetBondType() == Chem.BondType.DOUBLE and neighbor.GetAtomicNum() == 8:
            oxygen_neighbor = neighbor.GetIdx()
        else:
            neighbors.append(neighbor.GetIdx())

    if oxygen_neighbor is not None and len(neighbors) == 1:
        return {"kind": "carbonyl", "atom": dummy_idx, "anchor": neighbors[0]}
    if oxygen_neighbor is None and len(neighbors) == 1:
        return {"kind": "terminal", "atom": dummy_idx, "anchor": neighbors[0]}
    raise SmilesError("Constraint placeholder has an unsupported bonding pattern")


def _finalize_mol(rw_mol):
    mol = rw_mol.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def _carboxyl_hydroxyl_idx(mol, carbonyl_idx):
    carbonyl = mol.GetAtomWithIdx(carbonyl_idx)
    hydroxyls = [
        bond.GetOtherAtom(carbonyl).GetIdx()
        for bond in carbonyl.GetBonds()
        if bond.GetBondType() == Chem.BondType.SINGLE
        and bond.GetOtherAtom(carbonyl).GetAtomicNum() == 8
        and bond.GetOtherAtom(carbonyl).GetDegree() == 1
    ]
    return (
        hydroxyls[0]
        if len(hydroxyls) == 1 and _has_double_bonded_oxygen(carbonyl)
        else None
    )


def _residue_backbone(mol):
    matches = []
    for alpha in mol.GetAtoms():
        if alpha.GetAtomicNum() != 6:
            continue
        amines = [atom for atom in alpha.GetNeighbors() if atom.GetAtomicNum() == 7]
        acids = [
            atom
            for atom in alpha.GetNeighbors()
            if _carboxyl_hydroxyl_idx(mol, atom.GetIdx()) is not None
        ]
        if len(amines) == len(acids) == 1:
            matches.append((alpha, amines[0], acids[0]))
    if len(matches) != 1:
        raise SmilesError("Could not uniquely identify an alpha-amino-acid backbone")
    return matches[0]


def _reactive_site_indices(mol):
    alpha, amine, carboxyl = _residue_backbone(mol)
    excluded = {alpha.GetIdx(), amine.GetIdx(), carboxyl.GetIdx()}
    excluded.update(atom.GetIdx() for atom in carboxyl.GetNeighbors())
    sites = {"cterm": [], "nterm": [], "ester": [], "disulphide": []}

    for atom in mol.GetAtoms():
        if atom.GetIdx() in excluded:
            continue
        if _carboxyl_hydroxyl_idx(mol, atom.GetIdx()) is not None:
            sites["cterm"].append(atom.GetIdx())
        if (
            atom.GetAtomicNum() == 16
            and atom.GetDegree() == 1
            and atom.GetTotalNumHs() > 0
        ):
            sites["disulphide"].append(atom.GetIdx())
        if (
            atom.GetAtomicNum() == 8
            and atom.GetDegree() == 1
            and atom.GetTotalNumHs() > 0
        ):
            if _carboxyl_hydroxyl_idx(mol, atom.GetNeighbors()[0].GetIdx()) is None:
                sites["ester"].append(atom.GetIdx())
        if (
            atom.GetAtomicNum() == 7
            and not atom.GetIsAromatic()
            and atom.GetTotalNumHs() > 0
        ):
            if not any(
                _carboxyl_hydroxyl_idx(mol, neighbor.GetIdx()) is not None
                for neighbor in atom.GetNeighbors()
            ):
                sites["nterm"].append(atom.GetIdx())
    return sites


def _annotate_registry_site(residue, residue_idx, kind):
    mol = _mol_from_smiles(residue.smiles, context=f"residue {residue.id}")
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    sites = _reactive_site_indices(mol)[kind]
    if len(sites) != 1:
        raise BondSpecError(
            f"Residue {residue.id} has {len(sites)} unambiguous {kind} reaction sites; "
            "selecting a site explicitly is not supported yet"
        )
    mol.GetAtomWithIdx(sites[0]).SetAtomMapNum(1000 + residue_idx)
    return mol


def _marked_site_idx(mol, residue_idx, kind):
    marker = 1000 + residue_idx
    matches = [
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() == marker
    ]
    if len(matches) != 1:
        raise SmilesError(
            f"Could not retain selected {kind} site for residue {residue_idx}"
        )
    return matches[0]


def _connect_acid_to_atom(mol, acid_idx, atom_idx):
    hydroxyl_idx = _carboxyl_hydroxyl_idx(mol, acid_idx)
    if hydroxyl_idx is None:
        raise SmilesError("Selected sidechain acid is not a free carboxylic acid")
    rw_mol = Chem.RWMol(mol)
    rw_mol.AddBond(acid_idx, atom_idx, Chem.BondType.SINGLE)
    rw_mol.RemoveAtom(hydroxyl_idx)
    return _finalize_mol(rw_mol)


def _connect_atom_to_c_terminus(mol, atom_idx):
    carbonyl_idx, hydroxyl_idx = _find_terminal_carboxyl(mol)
    rw_mol = Chem.RWMol(mol)
    rw_mol.AddBond(atom_idx, carbonyl_idx, Chem.BondType.SINGLE)
    rw_mol.RemoveAtom(hydroxyl_idx)
    return _finalize_mol(rw_mol)


def _connect_registry_constraint(mol, marked_sites, pattern):
    if pattern.startswith("SS"):
        (_, left_kind), (_, right_kind) = marked_sites
        left_idx = _marked_site_idx(mol, marked_sites[0][0], left_kind)
        right_idx = _marked_site_idx(mol, marked_sites[1][0], right_kind)
        rw_mol = Chem.RWMol(mol)
        rw_mol.AddBond(left_idx, right_idx, Chem.BondType.SINGLE)
        return _finalize_mol(rw_mol)

    if len(marked_sites) == 1:
        residue_idx, kind = marked_sites[0]
        site_idx = _marked_site_idx(mol, residue_idx, kind)
        if kind in {"nterm", "ester"}:
            return _connect_atom_to_c_terminus(mol, site_idx)
        return _connect_acid_to_atom(mol, site_idx, _find_n_terminal_amine(mol))

    acid_site = next((site for site in marked_sites if site[1] == "cterm"), None)
    donor_site = next(
        (site for site in marked_sites if site[1] in {"nterm", "ester"}), None
    )
    if acid_site is None or donor_site is None:
        raise BondSpecError(f"{pattern} does not specify a compatible sidechain pair")
    return _connect_acid_to_atom(
        mol,
        _marked_site_idx(mol, *acid_site),
        _marked_site_idx(mol, *donor_site),
    )


def _registry_constrained_peptide_smiles(peptideseq, pattern, registry):
    tokens = parse_sequence_tokens(peptideseq, registry)
    if not pattern:
        return peptideseq, "", linear_peptide_smiles(peptideseq, registry=registry)
    if pattern == "HT":
        mol = _build_peptide_mol(registry.resolve(token).smiles for token in tokens)
        return (
            peptideseq,
            pattern,
            Chem.MolToSmiles(_cyclize_head_to_tail(mol), isomericSmiles=True),
        )

    kind_for_code = {"C": "disulphide", "Z": "cterm", "N": "nterm", "E": "ester"}
    if pattern[:2] not in {"SS", "SC"} or len(pattern[2:]) != len(tokens):
        raise BondSpecError(
            f"{pattern} is not a valid constraint pattern for this sequence"
        )
    marked_sites = [
        (index, kind_for_code[code])
        for index, code in enumerate(pattern[2:])
        if code in kind_for_code
    ]
    if any(code not in {*kind_for_code, "X"} for code in pattern[2:]):
        raise BondSpecError(f"{pattern} contains an unknown constraint code")
    kinds = [kind for _, kind in marked_sites]
    valid_shapes = (
        (pattern.startswith("SS") and kinds == ["disulphide", "disulphide"])
        or (pattern.startswith("SC") and kinds in (["cterm"], ["nterm"], ["ester"]))
        or (
            pattern.startswith("SC")
            and len(kinds) == 2
            and "cterm" in kinds
            and any(kind in {"nterm", "ester"} for kind in kinds)
        )
    )
    if not valid_shapes:
        raise BondSpecError(f"{pattern} does not select the required reactive sites")

    fragments = []
    for index, token in enumerate(tokens):
        selected = next(
            (kind for site_index, kind in marked_sites if site_index == index), None
        )
        residue = registry.resolve(token)
        fragment = (
            _annotate_registry_site(residue, index, selected)
            if selected
            else _mol_from_smiles(residue.smiles, context=f"residue {residue.id}")
        )
        if not selected:
            for atom in fragment.GetAtoms():
                atom.SetAtomMapNum(0)
        fragments.append(fragment)
    mol = _build_peptide_mol(fragments)
    mol = _connect_registry_constraint(mol, marked_sites, pattern)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return peptideseq, pattern, Chem.MolToSmiles(mol, isomericSmiles=True)


def _connect_two_dummies(mol, bond_type=Chem.BondType.SINGLE):
    dummy_idxs = _dummy_indices(mol)
    if len(dummy_idxs) != 2:
        raise SmilesError(
            "Expected exactly two constraint placeholders for cyclization"
        )

    left = _dummy_details(mol, dummy_idxs[0])
    right = _dummy_details(mol, dummy_idxs[1])

    rw_mol = Chem.RWMol(mol)
    if left["kind"] == "carbonyl" and right["kind"] == "carbonyl":
        raise SmilesError("Cannot connect two carbonyl-center placeholders")

    if left["kind"] == "terminal" and right["kind"] == "terminal":
        rw_mol.AddBond(left["anchor"], right["anchor"], bond_type)
        remove_idxs = dummy_idxs
    else:
        carbonyl = left if left["kind"] == "carbonyl" else right
        terminal = right if left["kind"] == "carbonyl" else left
        rw_mol.GetAtomWithIdx(carbonyl["atom"]).SetAtomicNum(6)
        rw_mol.AddBond(terminal["anchor"], carbonyl["atom"], bond_type)
        remove_idxs = [terminal["atom"]]

    for idx in sorted(remove_idxs, reverse=True):
        rw_mol.RemoveAtom(idx)
    return _finalize_mol(rw_mol)


def _connect_dummy_to_atom(mol, atom_idx, bond_type=Chem.BondType.SINGLE):
    dummy_idxs = _dummy_indices(mol)
    if len(dummy_idxs) != 1:
        raise SmilesError("Expected exactly one constraint placeholder for cyclization")

    details = _dummy_details(mol, dummy_idxs[0])
    rw_mol = Chem.RWMol(mol)
    if details["kind"] == "carbonyl":
        rw_mol.GetAtomWithIdx(details["atom"]).SetAtomicNum(6)
        rw_mol.AddBond(details["atom"], atom_idx, bond_type)
    else:
        rw_mol.AddBond(details["anchor"], atom_idx, bond_type)
        rw_mol.RemoveAtom(details["atom"])
    return _finalize_mol(rw_mol)


def _connect_dummy_to_c_terminus(mol, bond_type=Chem.BondType.SINGLE):
    dummy_idxs = _dummy_indices(mol)
    if len(dummy_idxs) != 1:
        raise SmilesError(
            "Expected exactly one constraint placeholder for C-terminal cyclization"
        )

    details = _dummy_details(mol, dummy_idxs[0])
    if details["kind"] != "terminal":
        raise SmilesError("C-terminal cyclization requires a terminal placeholder")

    carbonyl_idx, hydroxyl_idx = _find_terminal_carboxyl(mol)

    rw_mol = Chem.RWMol(mol)
    rw_mol.AddBond(details["anchor"], carbonyl_idx, bond_type)
    for idx in sorted([details["atom"], hydroxyl_idx], reverse=True):
        rw_mol.RemoveAtom(idx)
    return _finalize_mol(rw_mol)


def _cyclize_head_to_tail(mol):
    n_term_idx = _find_n_terminal_amine(mol)
    carbonyl_idx, hydroxyl_idx = _find_terminal_carboxyl(mol)

    rw_mol = Chem.RWMol(mol)
    rw_mol.AddBond(n_term_idx, carbonyl_idx, Chem.BondType.SINGLE)
    rw_mol.RemoveAtom(hydroxyl_idx)
    return _finalize_mol(rw_mol)


def _build_peptide_mol(fragment_smiles):
    fragments = list(fragment_smiles)
    if not fragments:
        raise SmilesError("Cannot build a peptide from an empty sequence")

    def as_mol(fragment):
        if isinstance(fragment, Chem.Mol):
            return Chem.Mol(fragment)
        return _mol_from_smiles(fragment, context="residue")

    peptide = as_mol(fragments[0])
    for fragment in fragments[1:]:
        peptide = _run_peptide_bond_reaction(peptide, as_mol(fragment))
    return peptide


def _sequence_smiles(peptideseq, registry=None):
    if registry is None:
        return [return_smiles(resi) for resi in peptideseq]

    tokens = parse_sequence_tokens(peptideseq, registry)
    return [registry.resolve(token).smiles for token in tokens]


def linear_peptide_smiles(peptideseq, registry=None):
    """Build a linear peptide with RDKit SMARTS peptide-bond reactions."""
    if not peptideseq:
        return None

    mol = _build_peptide_mol(_sequence_smiles(peptideseq, registry=registry))
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def bond_counter(peptidesmiles):
    # Count and return the highest bond number found in the SMILES string
    return max([int(num) for num in peptidesmiles if num.isdigit()], default=0)


def pep_positions(linpepseq):
    # Calculate starting positions of residues in the linear peptide SMILES
    positions = []
    location = 0
    for resi in linpepseq:
        positions.append(location)
        location += len(return_smiles(resi)) - 1
    return positions


# Constrained peptide SMILES generator
def constrained_peptide_smiles(peptideseq, pattern, next_bond_id=None, registry=None):
    """
    Build constrained peptide SMILES with RDKit SMARTS assembly.

    The next_bond_id parameter is retained for API compatibility but is no longer
    used because cyclization is performed on molecular graphs rather than with
    SMILES ring indices.
    """
    if registry is not None:
        return _registry_constrained_peptide_smiles(peptideseq, pattern, registry)

    valid_codes = {"C": "disulphide", "Z": "cterm", "N": "nterm", "E": "ester", "X": ""}

    if not pattern:
        return peptideseq, "", linear_peptide_smiles(peptideseq)

    if pattern[:2] == "HT":
        mol = _build_peptide_mol(return_smiles(resi) for resi in peptideseq)
        mol = _cyclize_head_to_tail(mol)
        return peptideseq, pattern, Chem.MolToSmiles(mol, isomericSmiles=True)

    fragment_smiles = []
    for resi, code in zip(peptideseq, pattern[2:]):
        if code not in valid_codes:
            raise BondSpecError(f"{code} in pattern {pattern} not recognised")
        constraint = valid_codes[code]
        fragment_smiles.append(
            return_constrained_smiles(resi, constraint)
            if constraint
            else return_smiles(resi)
        )

    mol = _build_peptide_mol(fragment_smiles)

    pf = pattern.replace("X", "")
    if pf in {"SCN", "SCE"}:
        mol = _connect_dummy_to_c_terminus(mol)
    elif pf == "SCZ":
        mol = _connect_dummy_to_atom(mol, _find_n_terminal_amine(mol))
    elif pf.startswith("SC"):
        mol = _connect_two_dummies(mol)
    elif pf.startswith("SS"):
        mol = _connect_two_dummies(mol)
    else:
        raise BondSpecError(f"{pattern} not recognised as valid bond_def")

    return peptideseq, pattern, Chem.MolToSmiles(mol, isomericSmiles=True)


# generate structures from sequences with specified constraints
def gen_structs_from_seqs(
    sequences,
    ssbond=False,
    htbond=False,
    scctbond=False,
    scntbond=False,
    scscbond=False,
    linear=False,
):
    funcs = [
        (ssbond, can_ssbond),
        (htbond, can_htbond),
        (scctbond, can_scctbond),
        (scntbond, can_scntbond),
        (scscbond, can_scscbond),
    ]
    next_bond_id = 1

    for seq in sequences:
        emitted = False
        for check, func in funcs:
            if not check:
                continue
            result = func(seq)
            if not result:
                continue
            seq2, bonddef = result
            # patterns that consume a bond id
            needs_ring_id = bonddef.startswith("HT") or bonddef.startswith("SC")
            if needs_ring_id:
                out = constrained_peptide_smiles(seq2, bonddef, next_bond_id)
                next_bond_id += 1
            else:
                out = constrained_peptide_smiles(seq2, bonddef)
            yield out
            emitted = True

        if linear or not emitted:
            yield (seq, "", linear_peptide_smiles(seq))


def gen_library_structs(
    liblen,
    ssbond=False,
    htbond=False,
    scctbond=False,
    scntbond=False,
    scscbond=False,
    linear=False,
):
    # Generate peptide structures for library based on sequence length and constraints
    for peptideseq, bond_def in gen_library_strings(
        liblen, ssbond, htbond, scctbond, scntbond, scscbond, linear
    ):
        if bond_def == "":
            yield (peptideseq, "", linear_peptide_smiles(peptideseq))
        else:
            yield constrained_peptide_smiles(peptideseq, bond_def)


def filtered_output(output, filterfuncs, key=None):
    # Filter the output items based on provided functions
    for out_item in output:
        if key:
            if all(func(key(out_item)) for func in filterfuncs):
                yield out_item
        else:
            if all(func(out_item) for func in filterfuncs):
                yield out_item


def get_constraint_type(bond_def):
    # Determine the constraint type from a bond definition string
    type_id, defi = bond_def[:2], bond_def[2:]
    if defi == "":
        return "linear" if type_id == "" else type_id
    if type_id == "SS" and all(char in ["X", "C"] for char in defi):
        return "SS"
    if type_id == "SC":
        if all(char in ["X", "Z", "E", "N"] for char in defi):
            if defi.count("X") == len(defi) - 1:
                if "N" in defi or "E" in defi:
                    return "SCCT"
                if "Z" in defi:
                    return "SCNT"
            elif defi.count("X") == len(defi) - 2:
                return "SCSC"
    raise BondSpecError(f"{bond_def} not recognised as valid bond_def")


def count_constraint_types(inlist, ignore_errors=False):
    # Count the number of peptides for each constraint type
    count_dict = {
        "linear": 0,
        "SS": 0,
        "HT": 0,
        "SCSC": 0,
        "SCCT": 0,
        "SCNT": 0,
    }
    for pep in inlist:
        try:
            count_dict[get_constraint_type(pep[1])] += 1
        except Exception:
            if ignore_errors:
                continue
            else:
                raise
    return count_dict


def save_3Dmolecule(sequence, bond_def):
    # Generate and save a 3D structure file (SDF) for peptide with given bond definition
    fname = f"{''.join(sequence)}_{bond_def}.sdf"
    _, _, smiles = constrained_peptide_smiles(sequence, bond_def)
    mol = Chem.MolFromSmiles(smiles)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    writer = AllChem.SDWriter(fname)
    writer.write(mol)
    return fname


# --- update write_molecule to avoid 3D unless write=="structure" and minimise=True ---
def write_molecule(
    smiles,
    peptideseq,
    bond_def,
    outfldr,
    type="sdf",
    write="structure",
    return_struct=False,
    new_folder=True,
    minimise=False,
):
    twodfolder = threedfolder = outfldr
    if not return_struct and new_folder:
        twodfolder = path.join(outfldr, "2D-Files")
        threedfolder = path.join(outfldr, "3D-Files")

    bond_def = f"_{bond_def}" if bond_def else "_linear"
    try:
        name = peptideseq + bond_def
    except TypeError:
        try:
            name = (
                "".join([aminodata[resi]["Letter"] for resi in peptideseq]) + bond_def
            )
        except KeyError:
            name = ",".join(peptideseq) + bond_def

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise SmilesError(f"{smiles} returns None molecule")
    mol.SetProp("_Name", name)

    if write == "draw":
        if not path.exists(twodfolder):
            os.makedirs(twodfolder)
        AllChem.Compute2DCoords(mol)  # fast 2D
        Draw.MolToFile(mol, path.join(twodfolder, name + ".png"), size=(1000, 1000))
    elif write == "structure":
        if not path.exists(threedfolder):
            os.makedirs(threedfolder)
        if minimise:
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)
        else:
            # keep it 2D; most tools accept 2D SDF; much faster
            AllChem.Compute2DCoords(mol)
        block = Chem.MolToMolBlock(mol)
        if return_struct:
            return block
        else:
            with open(
                path.join(threedfolder, name + "." + type), "w", encoding="utf-8"
            ) as handle:
                handle.write(block)
    else:
        raise TypeError(f'"write" must be "draw" or "structure", got {write}')
    return True


def write_library(inputlist, outloc, write="text", minimise=False, write_to_file=False):
    # Write the peptide library output to file (text, drawn images, or structure files).
    count = 0
    if write == "text":
        with open(outloc, "w") as f:
            for peptide in inputlist:
                try:
                    seq, bond_def, smiles = peptide
                    bond_def = bond_def if bond_def else "linear"
                    f.write(f"{''.join(seq)}-{bond_def}: {smiles}\n")
                    count += 1
                except Exception as e:
                    print(e)
    # Handle drawing or structure writing
    elif write in {"draw", "structure"}:
        if write_to_file:
            with open(outloc, "w") as out:
                for peptide in inputlist:
                    peptideseq, bond_def, smiles = peptide
                    if not (peptideseq or bond_def or smiles):
                        continue
                    mol = Chem.MolFromSmiles(smiles)
                    if write == "structure":
                        # default to 2D unless minimise=True (caller controls)
                        AllChem.Compute2DCoords(mol)
                        name = ",".join(map(str, peptideseq)) + (
                            f"_{bond_def}" if bond_def else "_linear"
                        )
                        mol.SetProp("_Name", name)
                        out.write(Chem.MolToMolBlock(mol) + "\n$$$$\n")
                        count += 1
                    else:
                        # defer to write_molecule for PNGs
                        write_molecule(
                            smiles,
                            peptideseq,
                            bond_def,
                            path.dirname(outloc),
                            write="draw",
                        )
                        count += 1
        else:
            for peptide in inputlist:
                seq, bond_def, smiles = peptide
                try:
                    write_molecule(smiles, seq, bond_def, outloc, write=write)
                    count += 1
                except Exception as e:
                    print(e)
    else:
        raise TypeError(f'"write" must be set to "draw" or "structure", got {write}')
    return count


def main(pattern, out_file):
    # Main function: generate peptides matching a pattern and write to file.
    print(f"Writing all peptides for pattern {pattern}")
    out_f = f"{out_file}.sdf"
    peptides = gen_all_matching_peptides(pattern)
    structures = gen_structs_from_seqs(peptides, True, True, True, True, True, True)
    write_library(structures, out_f, "structure", False, True)


if __name__ == "__main__":
    # Execute main with command-line arguments
    main(*sys.argv[1:], sys.argv[0])
