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


def can_ssbond(peptideseq):
    """Disulphide: need at least two Cys-like residues; pick the pair with max separation (>=3 apart)."""
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


def can_htbond(peptideseq):
    """Your original heuristic: qualifies if len >= 5 or exactly 2."""
    letters = _normalize_seq_letters(peptideseq)
    if len(letters) >= 5 or len(letters) == 2:
        return _preserve_seq_type(peptideseq, letters), "HT"
    return False


def can_scntbond(peptideseq, strict=False):
    """Sidechain → C-terminal (via N-term constraint code 'Z' position)."""
    letters = _normalize_seq_letters(peptideseq)
    cterm = _CONSTRAINT_LETTER_SETS["cterm"]
    locs = [i for i, r in enumerate(letters[3:], start=3) if r in cterm]
    if not locs or (len(locs) > 1 and strict):
        return False
    idx = locs[-1]  # keep your previous "last occurrence" behavior
    pattern = ["SC"] + ["Z" if i == idx else "X" for i in range(len(letters))]
    return _preserve_seq_type(peptideseq, letters), "".join(pattern)


def can_scctbond(peptideseq, strict=False):
    """Sidechain ↔ C-term using N-term/ester site: encode 'N' or 'E' at the site."""
    letters = _normalize_seq_letters(peptideseq)
    esters = _CONSTRAINT_LETTER_SETS["ester"]
    nterms = _CONSTRAINT_LETTER_SETS["nterm"]

    locs = [(i, "N") for i, r in enumerate(letters[:-3]) if r in nterms]
    locs += [(i, "E") for i, r in enumerate(letters[:-3]) if r in esters]
    if not locs or (len(locs) > 1 and strict):
        return False

    i0, code = locs[0]  # match previous behavior (first eligible)
    pattern = ["SC"] + [code if i == i0 else "X" for i in range(len(letters))]
    return _preserve_seq_type(peptideseq, letters), "".join(pattern)


def can_scscbond(peptideseq, strict=False):
    """Sidechain-to-sidechain: choose (cterm_pos, partner_pos) with max separation >= 2.
    Encode 'Z' at cterm_pos and 'N'/'E' at partner_pos depending on site set.
    """
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


def what_constraints(peptideseq):
    return [
        res
        for res in (
            can_ssbond(peptideseq),
            can_htbond(peptideseq),
            can_scctbond(peptideseq),
            can_scntbond(peptideseq),
            can_scscbond(peptideseq),
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


def linear_peptide_smiles(peptideseq):
    """
    Build linear peptide SMILES by concatenating residue fragments.
    Start with 'O', then for each residue:
      (1) trim the last character of the current chain
      (2) append the full residue SMILES fragment
    Assumes each residue fragment in aminodata['SMILES'] ends with a connector
    atom that replaces the trimmed char from the chain (e.g., ...O).
    """
    if not peptideseq:
        return None

    parts = ["O"]  # starting oxygen atom (matches your existing convention)

    for resi in peptideseq:
        frag = return_smiles(resi)  # full fragment, unmodified
        # trim the LAST CHAR of the CURRENT chain, not the fragment
        parts[-1] = parts[-1][:-1]
        parts.append(frag)

    return "".join(parts)


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
def constrained_peptide_smiles(peptideseq, pattern, next_bond_id=None):
    """
    Build constrained peptide SMILES.
    Uses next_bond_id (int) for '*' placeholders when needed; returns (seq, pattern, smiles).
    """
    valid_codes = {"C": "disulphide", "Z": "cterm", "N": "nterm", "E": "ester", "X": ""}
    smiles = "O"

    if not pattern:
        return peptideseq, "", linear_peptide_smiles(peptideseq)

    if pattern[:2] == "HT":
        smi = linear_peptide_smiles(peptideseq)
        bid = (bond_counter(smi) + 1) if next_bond_id is None else next_bond_id
        sbid = str(bid)
        smi = smi[0] + sbid + smi[1:-5] + sbid + smi[-5:-1]
        return peptideseq, pattern, smi

    for resi, code in zip(peptideseq, pattern[2:]):
        smiles = smiles[:-1]
        if code in valid_codes:
            smiles += (
                return_constrained_smiles(resi, valid_codes[code])
                if valid_codes[code]
                else return_smiles(resi)
            )
        elif code == "X":
            smiles += return_smiles(resi)
        else:
            raise BondSpecError(f"{code} in pattern {pattern} not recognised")

    pf = pattern.replace("X", "")
    if pf in {"SCN", "SCE"}:
        smiles = smiles[:-5] + "*(=O)"
    elif pf == "SCZ":
        smiles = "N*" + smiles[1:]

    bid = (bond_counter(smiles) + 1) if next_bond_id is None else next_bond_id
    smiles = smiles.replace("*", str(bid))
    return peptideseq, pattern, smiles


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
            with open(path.join(threedfolder, name + "." + type), "wb") as handle:
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
