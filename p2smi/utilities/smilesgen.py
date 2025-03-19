import itertools
import operator
import os
import os.path as path
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from p2smi.utilities.aminoacids import all_aminos

aminodata = all_aminos


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
    if name in all_aminos and name not in aminodata:
        aminodata[name] = all_aminos[name]
        return True
    else:
        raise UndefinedAminoError(f"{name} not recognised as valid amino acid")


def remove_amino(name):
    if name in aminodata:
        del aminodata[name]
    else:
        raise UndefinedAminoError(f"{name} not found in amino acids")


def print_possible_aminos():
    return list(all_aminos.keys())


def print_included_aminos():
    return list(aminodata.keys())


def return_available_residues(out="Letter"):
    return [properties[out] for properties in aminodata.values()]


def return_constraint_resis(constraint_type):
    return [
        name
        for name, properties in aminodata.items()
        if properties.get(constraint_type)
    ]


def property_to_name(property, value):
    for name, properties in aminodata.items():
        if properties.get(property) == value:
            return name
    raise UndefinedAminoError(
        f"Amino-acid matching value {value} for {property} not found"
    )


def gen_all_pos_peptides(pepliblen):
    amino_keys = list(aminodata.keys())
    for pep in itertools.product(amino_keys, repeat=pepliblen):
        yield pep


def gen_all_matching_peptides(pattern):
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
                if resi in aminodata:
                    outpep.append(resi)
                else:
                    outpep.append(property_to_name("Letter", resi))
            else:
                outpep.append(pep.pop(0))
        yield outpep


def can_ssbond(peptideseq):
    disulphides = return_constraint_resis("disulphide")
    locs = [
        loc
        for loc, resi in enumerate(peptideseq)
        if resi in disulphides or property_to_name("Letter", resi) in disulphides
    ]
    if len(locs) < 2:
        return False
    pos_pairs = sorted(
        [(pair, abs(pair[0] - pair[1])) for pair in itertools.combinations(locs, 2)],
        key=operator.itemgetter(1),
        reverse=True,
    )
    best_pair = pos_pairs[0]
    if best_pair[1] <= 2:
        return False
    pattern = "SS" + "".join(
        ["C" if i in best_pair[0] else "X" for i in range(len(peptideseq))]
    )
    return peptideseq, pattern


def can_htbond(peptideseq):
    if len(peptideseq) >= 5 or len(peptideseq) == 2:
        return peptideseq, "HT"
    return False


def can_scntbond(peptideseq, strict=False):
    locs = [
        num + 3
        for num, resi in enumerate(peptideseq[3:])
        if resi in return_constraint_resis("cterm")
        or property_to_name("Letter", resi) in return_constraint_resis("cterm")
    ]
    if len(locs) == 0 or (len(locs) > 1 and strict):
        return False
    if not strict:
        pattern = ["SC"] + [
            "Z" if num == locs[-1] else "X" for num, _ in enumerate(peptideseq)
        ]
        return peptideseq, "".join(pattern)
    return False


def can_scctbond(peptideseq, strict=False):
    esters = return_constraint_resis("ester")
    nterms = return_constraint_resis("nterm")
    locs = [
        (num, "N")
        for num, resi in enumerate(peptideseq[:-3])
        if resi in nterms or property_to_name("Letter", resi) in nterms
    ]
    locs += [
        (num, "E")
        for num, resi in enumerate(peptideseq[:-3])
        if resi in esters or property_to_name("Letter", resi) in esters
    ]
    if len(locs) == 0 or (len(locs) > 1 and strict):
        return False
    if not strict:
        pattern = ["SC"] + [
            locs[0][1] if num == locs[0][0] else "X" for num, _ in enumerate(peptideseq)
        ]
        return peptideseq, "".join(pattern)
    return False


def can_scscbond(peptideseq, strict=False):
    nterms = return_constraint_resis("nterm")
    cterms = return_constraint_resis("cterm")
    esters = return_constraint_resis("ester")
    locs = {"nterms": [], "cterms": [], "esters": []}
    for loc, resi in enumerate(peptideseq):
        for bondtype, bondname in [
            (nterms, "nterms"),
            (cterms, "cterms"),
            (esters, "esters"),
        ]:
            if resi in bondtype or property_to_name("Letter", resi) in bondtype:
                locs[bondname].append(loc)
    if not locs["cterms"] or not (locs["nterms"] or locs["esters"]):
        return False
    possible_pairs = [
        (pair, abs(pair[0] - pair[1]))
        for pair in itertools.product(locs["cterms"], locs["nterms"] + locs["esters"])
        if abs(pair[0] - pair[1]) >= 2
    ]
    if not possible_pairs:
        return False
    best_pair = max(possible_pairs, key=operator.itemgetter(1))[0]
    pattern = "SC" + "".join(
        [
            (
                "Z"
                if num == best_pair[0]
                else (
                    "N"
                    if num == best_pair[1] and best_pair[1] in locs["nterms"]
                    else (
                        "E"
                        if num == best_pair[1] and best_pair[1] in locs["esters"]
                        else "X"
                    )
                )
            )
            for num, _ in enumerate(peptideseq)
        ]
    )
    return peptideseq, pattern


def what_constraints(peptideseq):
    return [
        result
        for func in [
            can_ssbond,
            can_htbond,
            can_scctbond,
            can_scntbond,
            can_scscbond,
        ]
        if (result := func(peptideseq))
    ]


def aaletter2aaname(aaletter):
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
    with open(filepath) as peptides:
        for line in peptides:
            if line.startswith("#") or not line.strip():
                continue
            try:
                sequence, bond_def = map(str.strip, line.split(";"))
                if len(sequence.split(",")) == 1 and sequence not in all_aminos:
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
    mol = Chem.MolFromSmiles(smiles)
    n_pattern = Chem.MolFromSmarts("[$([Nh1](C)C(=O)),$([NH2]CC=O)]")
    methylated_pattern = Chem.MolFromSmarts("N(C)")
    rmol = AllChem.ReplaceSubstructs(
        mol, n_pattern, methylated_pattern, replaceAll=True
    )
    return Chem.MolToSmiles(rmol[0], isomericSmiles=True)


def nmethylate_peptides(structs):
    for struct in structs:
        seq, bond_def, smiles = struct
        if smiles:
            yield seq, bond_def, nmethylate_peptide_smiles(smiles)


def return_smiles(resi):
    return return_constrained_smiles(resi, "SMILES")


def return_constrained_smiles(resi, constraint):
    try:
        return aminodata[resi][constraint]
    except KeyError:
        try:
            return aminodata[property_to_name("Letter", resi)][constraint]
        except UndefinedAminoError:
            try:
                return aminodata[property_to_name("Code", resi)][constraint]
            except UndefinedAminoError:
                raise UndefinedAminoError(f"{resi} not recognised as amino acid")


def linear_peptide_smiles(peptideseq):
    if not peptideseq:
        return None
    combsmiles = "O"
    for resi in peptideseq:
        combsmiles = combsmiles[:-1] + return_smiles(resi)
    return combsmiles


def bond_counter(peptidesmiles):
    return max([int(num) for num in peptidesmiles if num.isdigit()], default=0)


def pep_positions(linpepseq):
    positions = []
    location = 0
    for resi in linpepseq:
        positions.append(location)
        location += len(return_smiles(resi)) - 1
    return positions


def constrained_peptide_smiles(peptideseq, pattern):
    valid_codes = {
        "C": "disulphide",
        "Z": "cterm",
        "N": "nterm",
        "E": "ester",
        "X": "",
    }  # codes for type of constraint (X is no constraint)
    smiles = "O"  # start with O as the first oxygen atom

    if not pattern:
        return peptideseq, "", linear_peptide_smiles(peptideseq)

    if pattern[:2] == "HT":
        smiles = linear_peptide_smiles(peptideseq)
        bond_num = str(bond_counter(smiles) + 1)
        smiles = smiles[0] + bond_num + smiles[1:-5] + bond_num + smiles[-5:-1]
        return peptideseq, pattern, smiles

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

    # edit n-term or c-term depending on pattern
    # remove all X in pattern
    pattern_for_fixing = pattern.replace("X", "")
    # print(pattern_for_fixing)

    if pattern_for_fixing == "SCN":  # N acts as N term, which binds the CT
        smiles = (
            smiles[:-5] + "*(=O)"
        )  # removes (=O)O and adds *(=O) to cyclize to the C before
    elif pattern_for_fixing == "SCE":  # E is an ester, which binds the CT
        smiles = smiles[:-5] + "*(=O)"
    elif (
        pattern_for_fixing == "SCZ"
    ):  # Identifies sidechain to N term (Z is acting as C-term)
        smiles = "N*" + smiles[1:]  # replaces amino group (N) with N*

    bond_number = str(bond_counter(smiles) + 1)
    smiles = smiles.replace("*", bond_number)

    return peptideseq, pattern, smiles


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
    for seq in sequences:
        for check, func in funcs:
            if check and (result := func(seq)):
                seq, bonddef = result
                yield constrained_peptide_smiles(seq, bonddef)
        if linear:
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
    for peptideseq, bond_def in gen_library_strings(
        liblen, ssbond, htbond, scctbond, scntbond, scscbond, linear
    ):
        if bond_def == "":
            yield (peptideseq, "", linear_peptide_smiles(peptideseq))
        else:
            yield constrained_peptide_smiles(peptideseq, bond_def)


def filtered_output(output, filterfuncs, key=None):
    for out_item in output:
        if key:
            if all(func(key(out_item)) for func in filterfuncs):
                yield out_item
        else:
            if all(func(out_item) for func in filterfuncs):
                yield out_item


def get_constraint_type(bond_def):
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
    fname = f"{''.join(sequence)}_{bond_def}.sdf"
    _, _, smiles = constrained_peptide_smiles(sequence, bond_def)
    mol = Chem.MolFromSmiles(smiles)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    writer = AllChem.SDWriter(fname)
    writer.write(mol)
    return fname


def write_molecule(
    smiles,
    peptideseq,
    bond_def,
    outfldr,
    type="sdf",
    write="structure",
    return_struct=False,
    new_folder=True,
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

    mymol = Chem.MolFromSmiles(smiles)
    if not mymol:
        raise SmilesError(f"{smiles} returns None molecule")
    mymol.SetProp("_Name", name)

    if write == "draw":
        if not path.exists(twodfolder):
            os.makedirs(twodfolder)
        Draw.MolToFile(mymol, path.join(twodfolder, name + ".png"), size=(1000, 1000))
    elif write == "structure":
        if not path.exists(threedfolder):
            os.makedirs(threedfolder)
        AllChem.EmbedMolecule(mymol)
        AllChem.UFFOptimizeMolecule(mymol)
        if return_struct:
            return Chem.MolToMolBlock(mymol)
        else:
            with open(path.join(threedfolder, name + "." + type), "wb") as handle:
                handle.write(Chem.MolToMolBlock(mymol))
    else:
        raise TypeError(f'"write" must be set to "draw" or "structure", got {write}')
    return True


def write_library(inputlist, outloc, write="text", minimise=False, write_to_file=False):
    count = 0
    if write == "text":
        with open(outloc, "w") as f:
            for peptide in inputlist:
                try:
                    seq, bond_def, smiles = peptide
                    bond_def = bond_def if bond_def else "linear"
                    f.write(f"{','.join(seq)}-{bond_def}: {smiles}\n")
                    count += 1
                except Exception as e:
                    print(e)
    elif write in {"draw", "structure"}:
        if write_to_file:
            with open(outloc, "w") as out:
                for peptide in inputlist:
                    peptideseq, bond_def, smiles = peptide
                    if not peptideseq and not bond_def and not smiles:
                        continue
                    mol = Chem.MolFromSmiles(smiles)
                    AllChem.EmbedMolecule(mol)
                    AllChem.UFFOptimizeMolecule(mol)
                    try:
                        name = peptideseq + bond_def
                    except TypeError:
                        try:
                            name = (
                                "".join(
                                    [aminodata[resi]["Letter"] for resi in peptideseq]
                                )
                                + bond_def
                            )
                        except KeyError:
                            name = ",".join(map(str, peptideseq)) + bond_def
                    mol.SetProp("_Name", name)
                    molstr = Chem.MolToMolBlock(mol)
                    out.write(molstr + "\n$$$$\n")
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
    print(f"Writing all peptides for pattern {pattern}")
    out_f = f"{out_file}.sdf"
    peptides = gen_all_matching_peptides(pattern)
    structures = gen_structs_from_seqs(peptides, True, True, True, True, True, True)
    write_library(structures, out_f, "structure", False, True)


if __name__ == "__main__":
    main(*sys.argv[1:], sys.argv[0])
