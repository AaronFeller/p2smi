#!/usr/bin/env python
"""Module containing functions for combinatorial generation of lists of 
constrained and unconstrained peptides, for producing the approriate SMILES
string, and using openbabel to generate a 3D structure (pdb, sdf, etc) from the 
SMILES

All peptide sequences are handled as a string of 1 letter residue codes

Author:Fergal Duffy. email:fergaljd@gmail.com

Requires python 2.6 or better

Recommended to import with:
from PepLibGen.StructGen import StructGen as sg

"""

import itertools
import sys
import os
import operator
import os.path as path

# From my package
from aminoacids import aminos
from aminoacids import special_aminos
from aminoacids import d_aminos
from aminoacids import unnatural_aminos

# All aminos - a combined dictionary of all defined amino acids
all_aminos = {}
for amino_dict in [aminos, special_aminos, d_aminos, unnatural_aminos]:
    all_aminos.update(amino_dict)
# Aminodata - amino acids available to program
# Default - a copy of all_aminos
aminodata = dict(all_aminos)

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem import AllChem
except ImportError as e:
    print(e)


class NoCysteineError(Exception):
    pass


class BondSpecError(Exception):
    pass


class FormatError(Exception):
    pass


class UndefinedAminoError(Exception):
    pass


class UndefinedPropertyError(Exception):
    pass


class SmilesError(Exception):
    pass


def add_amino(name):
    """
    If an amino acid is defined in aminoacids.py, add it to aminodata list.

    Default aminodata is working amino list, only has natural aminos.
    """
    if name in list(all_aminos.keys()) and name not in list(aminodata.keys()):
        aminodata[name] = all_aminos[name]
        return True
    else:
        error_str = "%s not recognised as valid amino acid" % (name)
        raise UndefinedAminoError(error_str)


def remove_amino(name):
    """
    Remove an amino acid from the working amino acid list, aminodata
    """
    if name in list(aminodata.keys()):
        del aminodata[name]
    else:
        error_str = "%s not found in amino acids" % (name)
        raise UndefinedAminoError(error_str)


def print_possible_aminos():
    """
    Returns a list of names of all amino acids available to programme.

    Includes aminos in current working library, and not.
    """
    return list(all_aminos.keys())


def print_included_aminos():
    """
    Prints amino acids included in working library.

    Wraps around return_available_residues()
    Provides a function mirroring print_possible_aminos()
    """
    return list(aminodata.keys())


def return_available_residues(out="Letter"):
    """
    This function returns a list with all the 1-letter amino codes of
    amino acids which have been read in by the loadaminodata() function

    Can also return 3-letter codes - pass 'code' as input to function

    """
    return [properties[out] for (name, properties) in list(aminodata.items())]


def return_constraint_resis(constraint_type):
    """
    look in aminodata and returns a list of resi names that can bond in
    particular ways.
    Constraint_type can be 'nterm', 'cterm', 'ester'. 'nterm' for amino acids
    that can mimic a peptide n-terminal, 'cterm' the same for peptide
    c-terminals, 'ester' for residues that can ester bond to a peptide
    n-terminus
    """
    return [
        name
        for (name, properties) in list(aminodata.items())
        if properties[constraint_type]
    ]


def property_to_name(property, value):
    """
    Given a value of property for a amino acid, return the amino acid name
    """
    for name in aminodata:
        try:
            if aminodata[name][property] == value:
                return name
        except KeyError:
            errortext = "%s not a valid amino acid property" % property
            raise UndefinedPropertyError(errortext)

    raise UndefinedAminoError(
        "Amino-acid matching value %s for %s not found" % (property, value)
    )


def gen_all_pos_peptides(pepliblen):
    """This function will return a list of peptide sequences, containing all
    itertools.combinations of peptides of length pepliblen, using all all amino
    acids read
    in by the loadaminodata() function.

    """
    # Loop calculates cartesian of peptidelist with itself peplislen times.
    # Equivalent to all possible peptides of length peplislen
    # Rewritten as a generator for more efficient memory usage
    for pep in itertools.product(list(aminodata.keys()), repeat=pepliblen):
        yield pep


def gen_all_matching_peptides(pattern):
    """This function will generate all possible peptides matching a pattern.
    Pattern will a string of either amino-acids codes or 'X's.
    Will return a list of peptide strings with 'X's combinatorially replaced
    with all other peptides.
    For example, CXXXC will return all possible tri-peptides surrounded by two
    cysteine

    pattern can be an array of amino acid names and 'X's, or a string
    """
    # Procedure: find all 'X's, pull them out of the sequence, generate them as
    # linear peptide itertools.combinations, re-insert non 'X's, and add to output

    # Also allow lower case 'x' as placeholder
    try:
        pattern = pattern.replace("x", "X")
    except AttributeError:  # Input is a list
        pattern = ["X" if resi == "x" else resi for resi in pattern]

    for pep in itertools.product(list(aminodata.keys()), repeat=pattern.count("X")):
        pep = list(pep)
        outpep = []
        for resi in pattern:
            if resi != "X":
                if resi in aminodata:
                    outpep.append(resi)
                else:
                    # If letter code passed in, get name
                    outpep.append(property_to_name("Letter", resi))
            else:
                outpep.append(pep.pop(0))
        yield outpep


def can_ssbond(peptideseq):
    """Function will return True if the supplied peptide sequence is suitable
    for SS bonding. Will check for 2 Cysteines in a peptide with adequate
    spacing.
    Will return False if unsuitable.

    If peptide sequence is suitable for bonding, will return the given peptide
    sequnce with a constraint identifier 'SS' appended"""
    disulphides = return_constraint_resis("disulphide")
    locs = []
    pos_pairs = []
    for loc, resi in enumerate(peptideseq):
        try:
            if resi in disulphides or property_to_name("Letter", resi) in disulphides:
                locs.append(loc)
        except UndefinedAminoError:

            if resi in aminodata or property_to_name("Letter", resi) in aminodata:
                pass
            else:
                raise
    if len(locs) < 2:
        return False
    else:
        pos_pairs = sorted(
            [
                (pair, abs(pair[0] - pair[1]))
                for pair in itertools.combinations(locs, 2)
            ],
            key=operator.itemgetter(1),
            reverse=True,
        )
    best_pair = pos_pairs[0]
    if best_pair[1] <= 2:
        return False
    # Build pattern
    pattern = "SS"
    for num, resi in enumerate(peptideseq):
        if num in best_pair[0]:
            pattern += "C"
        else:
            pattern += "X"

    return peptideseq, pattern


def can_htbond(peptideseq):
    """Function to check if supplied peptide sequence is suitable for head-tail
    bonding. Will return False if unsuitable
    Criteria for head-tail binding is sequence length greater than 4 residues

    If peptide sequence is suitable for bonding, will return the given peptide
    sequnce with a constraint identifier '-HT' appended
    """
    if len(peptideseq) >= 5 or len(peptideseq) == 2:
        try:
            return peptideseq, "HT"
        except TypeError:  # Peptideseq could be  list of refs
            return peptideseq, "HT"
    else:
        return False


def can_scntbond(peptideseq, strict=False):
    """This function check if it there is a unique bond from a c-terminal
    mimicing residue side-chain to the supplied peptide sequence N-terminal.
    Will return False if cannot be bonded

    If bond is possible, the function will return the input sequence plus a
    constraint identifier '-SC' plus a bond identification string equal to the
    peptide sequence string with non side-chain side-chain bonding residues
    replace by 'X'

    If strict is True, will only return a peptide if there is one possible bond
    Otherwise, will return the longest possible bond.
    """

    locs = []
    # Pattern of X's, len(peptideseq) long, with the bonding residue indicated
    pattern = ["SC"]
    cterms = return_constraint_resis("cterm")

    # Check if there are any of the above residue, with a gap to the n-ter
    for num, resi in enumerate(peptideseq[3:]):
        try:
            if resi in cterms or property_to_name("Letter", resi) in cterms:
                locs.append(num + 3)
        except UndefinedAminoError:
            # Note, sequences under 4 residues will return False, not raise an error
            if resi in aminodata or property_to_name("Letter", resi) in aminodata:
                pass
            else:
                raise

    if len(locs) == 0:
        return False
    elif len(locs) > 1 and strict == True:
        return False
    elif strict == False:  # Return longest bond
        for num, resi in enumerate(peptideseq):
            #'Z' will signify a cterminus sidechain in the pattern
            if num == locs[-1]:
                pattern.append("Z")
            else:
                pattern.append("X")
        pattern = "".join(pattern)
        return peptideseq, pattern

    else:
        return False


def can_scctbond(peptideseq, strict=False):
    """This function will check if there is a unique bond from an n-terminal
    mimicing sidechain to
    the C terminus of the supplied peptide sequence.
    Will return False if cannot be bonded

    If bond is possible, the function will return the input sequence plus a
    constraint identifier '-SC' plus a bond identification string equal to the
    peptide sequence
    string with non side-chain side-chain bonding residues replace by 'X'

    If strict is True, will only return a peptide if there is one possible bond
    Otherwise, will return the longest possible bond.

    """
    esters = return_constraint_resis("ester")
    nterms = return_constraint_resis("nterm")
    locs = []
    # Pattern of X's, len(peptideseq) long, with the bonding residue indicated
    pattern = ["SC"]

    for num, resi in enumerate(peptideseq[:-3]):
        try:
            if resi in nterms or property_to_name("Letter", resi) in nterms:
                locs.append((num, "N"))
        except UndefinedAminoError:
            if resi in aminodata or property_to_name("Letter", resi) in aminodata:
                pass
            else:
                raise
        try:
            if resi in esters or property_to_name("Letter", resi) in esters:
                locs.append((num, "E"))
        except UndefinedAminoError:
            # Note, sequences under 4 residues will return False, not raise an error
            if resi in aminodata or property_to_name("Letter", resi) in aminodata:
                pass
            else:
                raise

    if len(locs) == 0:
        return False
    elif len(locs) > 1 and strict == True:
        return False
    elif strict == False:  # Return longest bond (first position in locs)
        for num, resi in enumerate(peptideseq):
            if num == locs[0][0]:
                pattern.append(locs[0][1])
            else:
                pattern.append("X")
        pattern = "".join(pattern)
        return peptideseq, pattern

    else:
        return False


def can_scscbond(peptideseq, strict=False):
    """This function will check if it can side chain - side chain bond
    Lys, Ser, Thr to Glu, Asp in the supplied peptide sequence
    Will return False if cannot be bonded

    If bond is possible, the function will return the input sequence plus a
    constraint identifier '-SC' plus a bond identification string equal to the
    peptide sequence
    string with non side-chain side-chain bonding residues replace by 'X'

    If more then one possible SC-SC bond and strict == True - then this funtion
    will return nothing.
    If strict == False, will return first to last"""

    # Find locations of cterm and nterm (and ester) bonds.
    nterms = return_constraint_resis("nterm")
    cterms = return_constraint_resis("cterm")
    esters = return_constraint_resis("ester")
    allresis = nterms + cterms + esters
    locs = {"nterms": [], "cterms": [], "esters": []}
    possible_pairs = []
    pattern = "SC"

    for loc, resi in enumerate(peptideseq):
        for bondtype in [(nterms, "nterms"), (cterms, "cterms"), (esters, "esters")]:
            try:
                if (
                    resi in bondtype[0]
                    or property_to_name("Letter", resi) in bondtype[0]
                ):
                    locs[bondtype[1]].append(loc)
            except UndefinedAminoError:
                if resi in aminodata or property_to_name("Letter", resi) in aminodata:
                    pass
                else:
                    raise
    # Check for possible bond pairs...
    if len(locs["cterms"]) == 0:
        return False
    elif len(locs["nterms"]) == 0 and len(locs["esters"]) == 0:
        return False
    else:
        for pair in itertools.product(
            locs["cterms"], (locs["nterms"] + locs["esters"])
        ):
            pair_gap = abs(pair[0] - pair[1])
            if pair_gap >= 2:
                possible_pairs.append((pair, pair_gap))
        possible_pairs.sort(key=operator.itemgetter(1), reverse=True)
        try:
            best_pair = possible_pairs[0][0]
        except IndexError:  # No possible pairs
            return False
    # Generate pattern
    for num, resi in enumerate(peptideseq):
        if num == best_pair[0]:
            pattern += "Z"
        elif num == best_pair[1]:
            if best_pair[1] in locs["nterms"]:
                pattern += "N"
            else:
                pattern += "E"
        else:
            pattern += "X"
    return peptideseq, pattern


def what_constraints(peptideseq):
    """
    Takes in a peptide sequence,
    and returns a list of what constraints are possible.
    """
    outlist = []
    for func in [can_ssbond, can_htbond, can_scctbond, can_scntbond, can_scscbond]:
        if func(peptideseq):
            peptide, pattern = func(peptideseq)
            outlist.append(pattern)
    return outlist


def aaletter2aaname(aaletter):
    """Translate a 1 letter amino acid string to the full name.
    >>>aaletter2aaname('A')
    >>>'L-Alanine'"""
    for name, amino_dict in list(all_aminos.items()):
        if amino_dict["Letter"] == aaletter:
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
    """This function will cycle through all possible peptides for a given
    library length and return peptide sequences with the desired constraints
    as defined in the input.

    If no constraints are set True, all non-constrained peptides are returned

    This is a generator function.
    """
    filterfuncs = []
    if ssbond == True:
        filterfuncs.append(can_ssbond)
    if htbond == True:
        filterfuncs.append(can_htbond)
    if scctbond == True:
        filterfuncs.append(can_scctbond)
    if scntbond == True:
        filterfuncs.append(can_scntbond)
    if scscbond == True:
        filterfuncs.append(can_scscbond)
    for sequence in gen_all_pos_peptides(liblen):
        for func in filterfuncs:
            trialpeptide = func(sequence)
            if trialpeptide != False:
                yield trialpeptide
    if linear == True:
        for peptide in gen_all_pos_peptides(liblen):
            yield (peptide, "")  # Blank string for bond_def


def gen_library_from_file(filepath, ignore_errors=False):
    """Generate peptides from inputfile descriptions.
    Inputfile should have lines of format <sequence, bond_def>
    and output will be generated as (sequence, bond_def, smiles)."""
    with open(filepath) as peptides:
        for line in peptides:
            if line.startswith("#") or not line.strip():
                continue
            try:
                # import pdb; pdb.set_trace()
                sequence, bond_def = [s.strip() for s in line.split(";")]
                # Sequence needs to be a tuple of amino-acid names
                # The raw string can be like 'L-Alanine,L-Alanine,..' or 'AA'
                # If it can't be split on commas and the string isn't just an
                # amino acid name, translate it from single letter codes.
                if (len(sequence.split(",")) == 1) and (
                    sequence not in list(all_aminos.keys())
                ):
                    sequence = [aaletter2aaname(l) for l in sequence]
                else:
                    sequence = sequence.split(",")
                yield constrained_peptide_smiles(sequence, bond_def)
            except Exception as e:
                if ignore_errors:
                    yield (None, None, None)
                else:
                    raise


def nmethylate_peptide_smiles(smiles):
    """N-Methylate the  of a peptide smiles string.
    Input a peptide smiles string.
    Returns a smiles string with all backbone nitrogens and the n-terminus
    methylated."""
    mol = Chem.MolFromSmiles(smiles)
    # Pattern to find N molecules in a peptide bond or at the N-terminus.
    # N molecule must have a free H, ignore Proline N-termini in
    # middle of sequences
    n_pattern = Chem.MolFromSmarts("[$([Nh1](C)C(=O)),$([NH2]CC=O)]")
    methylated_pattern = Chem.MolFromSmarts("N(C)")
    rmol = AllChem.ReplaceSubstructs(
        mol, n_pattern, methylated_pattern, replaceAll=True
    )
    return Chem.MolToSmiles(rmol[0], isomericSmiles=True)


def nmethylate_peptides(structs):
    """Wrapper for nmethylate_peptide_smiles to use when passing in peptide
    tuples of the format (seq, bond_def, smiles).
    Returns (seq, bond_def, nmethylated-smiles)."""
    for struct in structs:
        seq, bond_def, smiles = struct
        if smiles:
            yield seq, bond_def, nmethylate_peptide_smiles(smiles)


def return_smiles(resi):
    """
    Given amino-acid name, 3-letter code or letter, return SMILES

    Works from aminodata dictionary.
    """
    return return_constrained_smiles(resi, "SMILES")


def return_constrained_smiles(resi, constraint):
    """
    Returns the smiles of the amino-acid to be constrained.
    Constrained smiles differ from normal smiles in that they contain a '#'
    where the bond number is to be inserted
    """
    try:
        return aminodata[resi][constraint]
    except KeyError:
        try:
            return aminodata[property_to_name("Letter", resi)][constraint]
        except UndefinedAminoError:
            try:
                return aminodata[property_to_name("Code", resi)][constraint]
            except UndefinedAminoError:
                error_string = "%s not recognised as amino acid" % (resi)
                raise UndefinedAminoError(error_string)


def linear_peptide_smiles(peptideseq):
    """This function uses the return_smiles() function to determine
    the SMILES of a linear   peptide chain
    It mimics the peptide bond by concatentating the amino acid SMILES while
    dropping off the carbonly OH group"""
    # O here is placeholder to be removed - O has to be removed at start
    # of loop to keep final oxygen in linear peptide smiles
    if peptideseq == "":
        return None
    combsmiles = "O"
    for resi in peptideseq:
        # Mimics peptide bonding by slicing off final oxygen,
        # before adding next peptide
        combsmiles = combsmiles[:-1]
        smiles = return_smiles(resi)
        combsmiles = combsmiles + smiles
    return combsmiles  # , peptidestring


def bond_counter(peptidesmiles):
    """This is a helper function for functions which seek to add new loops to a
    SMILES structure.
    It counts the number of bond numbers used, and returns the highest
    Using this function conservatively (by adding one to result to get the
    next available bondnumber)
    limits the maximum total number of loops in a SMILES structure to 9"""
    bond_counter = 0
    for num in range(1, 10, 1):
        if peptidesmiles.find(str(num)) != -1:
            if num > bond_counter:
                bond_counter = num
    return bond_counter


def pep_positions(linpepseq):
    """Returns a list of locations of the first atom of all residues in a linear
    peptide smiles string, given the peptidesequence"""
    positions = []
    location = 0  # Where in peptidestring we are
    for i, resi in enumerate(linpepseq):
        positions.append(location)
        # Remove last 'O' molecule in peptide bond
        location += len(return_smiles(linpepseq[i])) - 1
    return positions


def constrained_peptide_smiles(peptideseq, pattern):
    """Return a SMILES string of a constrained peptide given an
    appropriate     peptide sequence and constraint spec pattern.

    The constraint spec is a string with the first 2 characters
    identifiying the constraint type, (HT, SS, SC),
    then one character per residue in the sequence, which is 'X' for
    residues not involved in constraints, 'C' for
    disulphide bonding residues, 'N' for n-terminal like sidechains, 'Z' for
    c-terminal like sidechains, and 'E' for ester bonding sidechains.

    Returns a tuple (peptideseq, pattern, SMILES)"""
    smiles = "O"
    usedcodes = ""
    codes = {"C": "disulphide", "Z": "cterm", "N": "nterm", "E": "ester"}

    if not pattern:
        return peptideseq, "", linear_peptide_smiles(peptideseq)

    if pattern[:2] == "HT":
        smiles = linear_peptide_smiles(peptideseq)
        bond_num = str(bond_counter(smiles) + 1)
        smiles = smiles[0] + bond_num + smiles[1:-5] + bond_num + smiles[-5:-1]
        return peptideseq, pattern, smiles

    # Put together smiles string with '*' placeholder
    for resi, code in zip(peptideseq, pattern[2:]):
        smiles = smiles[:-1]  # Remove final oxgen before each residue is added
        if code in list(codes.keys()):
            try:
                smiles += return_constrained_smiles(resi, codes[code])
            except TypeError:
                print(peptideseq, pattern)
                raise
            usedcodes += code
        elif code == "X":
            smiles += return_smiles(resi)
        else:
            raise BondSpecError(
                "%s in pattern %s not recognised" % code
            ).with_traceback(pattern)

            # Add in '*' for SCNT and SCCT bonds
    if len(usedcodes) == 1:
        if usedcodes == "N":  # SCCTbond
            smiles = smiles[:-1] + "*"
        elif usedcodes == "E":  # SCCT ester bond
            smiles = smiles[:-1] + "*"
        elif usedcodes == "Z":  # SCNT bond, insert # by Nterm
            smiles = smiles[0] + "*" + smiles[1:]

    # Make sure that only one constraint is in the pattern
    if smiles.count("*") > 2:
        raise BondSpecError(
            "Cannot handle more than one specified\
            constraint, found %s in sequence %s, \nsmiles %s"
            % (smiles.count("*"), peptideseq, smiles)
        )

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
    """
    Returns an iterator generating SMILES, given an iterator of sequences

    Will generate constrained SMILES according to what constraint are True
    and False
    """
    funcs = [
        (ssbond, can_ssbond),
        (htbond, can_htbond),
        (scctbond, can_scctbond),
        (scntbond, can_scntbond),
        (scscbond, can_scscbond),
    ]
    for seq in sequences:
        # print '\nSeq', seq
        for check, func in funcs:
            # print 'Check, func, func(seq)', check, func, func(seq)
            if check and func(seq):
                # print 'True', seq, func
                seq, bonddef = func(seq)
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
    """
    Given a desired peptide library length, and constraints, returns a list
    of tuples of peptide names (sequence-constraint descriptor) and SMILES
    strings

    If no constraints are set True, all non-constrained peptides are returned
    Can filter out peptides that break synthesis rules
    """
    for peptideseq, bond_def in gen_library_strings(
        liblen, ssbond, htbond, scctbond, scntbond, scscbond, linear
    ):
        if bond_def == "":
            yield (peptideseq, "", linear_peptide_smiles(peptideseq))
        else:
            yield constrained_peptide_smiles(peptideseq, bond_def)


def filtered_output(output, filterfuncs, key=None):
    """
    Applies functions to conditionally filter output (any iterable).

    filterfuncs must be return a True/False result for any output item from
    main func.
    Can set key to a function to return the element to be filtered by in each
    output item of main func
    """
    for out_item in output:
        if key:
            # If any test fails, do not yield item
            try:
                if [
                    func(key(out_item))
                    for func in filterfuncs
                    if func(key(out_item)) == False
                ]:
                    pass
                else:
                    yield out_item
            except TypeError:  # if filterfuncs not a list
                try:
                    if filterfuncs(key(out_item)):
                        yield out_item
                except Exception as e:
                    print("Out_item", out_item)
                    print("Key", key)
                    raise
        else:
            # If any test fails, do not yield item
            try:
                if [func(out_item) for func in filterfuncs if func(out_item) == False]:
                    pass
                else:
                    yield out_item
            except TypeError:  # if filterfuncs not a list
                if filterfuncs(out_item):
                    yield out_item


def get_constraint_type(bond_def):
    """
    pass in peptide bond_def, will return the type
    """
    type_id = bond_def[:2]
    defi = bond_def[2:]
    errortext = "%s not recognised as valid bond_def"

    if defi == "":  # Must be linear or HT
        if type_id == "":
            return "linear"
        elif type_id == "HT":
            return "HT"
        else:
            raise BondSpecError(errortext % bond_def)
    else:
        if type_id == "SS":
            if not [char for char in defi if char not in ["X", "C"]]:
                return "SS"
            else:
                raise BondSpecError("Invalid chars found in bond_def")
        elif type_id == "SC":
            if [char for char in defi if char not in ["X", "Z", "E", "N"]]:
                raise BondSpecError("Invalid chars found in bond_def")
            if defi.count("X") == len(defi) - 1:
                if "N" in defi or "E" in defi:
                    return "SCCT"
                elif "Z" in defi:
                    return "SCNT"
                else:
                    raise BondSpecError(errortext % bond_def)
            elif defi.count("X") == len(defi) - 2:
                return "SCSC"
            else:
                raise BondSpecError(errortext % bond_def)
        else:
            raise BondSpecError(errortext % bond_def)


def count_constraint_types(inlist, ignore_errors=False):
    """
    Returns the numbers of each type of constraints in inlist.


    List should be of the form (name, bond_def[, smiles])
    Works by figuring out the constraint type from the name
    returns a dictionary of type: number

    Can filter by synthesisabilty
    """
    # BondSpecError
    count_dict = {"linear": 0, "SS": 0, "HT": 0, "SCSC": 0, "SCCT": 0, "SCNT": 0}
    errortext = "pep {0} not recognised"
    for pep in inlist:
        try:
            count_dict[get_constraint_type(pep[1])] += 1
        except Exception as e:
            if ignore_errors:
                pass
            else:
                raise
    return count_dict


def save_3Dmolecule(sequence, bond_def):
    """Write a constrained peptide to an .sdf file.
    Will be saved to sequence_bonddef.sdf"""
    fname = "%s_%s.sdf" % ("".join(sequence), bond_def)
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
    """This function will take in the SMILES and namestring of the peptide,
    and write it to a either 3D PDBfile, or a 2D .png representation of the
    structure depending on whether write='draw' or write='structure'
    """

    if not return_struct:
        if new_folder:
            twodfolder = path.join(outfldr, "2D-Files")
            threedfolder = path.join(outfldr, "3D-Files")
        else:
            twodfolder = outfldr
            threedfolder = outfldr

    else:
        twodfolder = outfldr
        threedfolder = outfldr

    if bond_def:
        bond_def = "_" + bond_def
    else:
        bond_def = "_linear"
    try:
        name = peptideseq + bond_def
    except TypeError:  # If peptideseq is not a string
        try:
            name = "".join([aminos[resi]["Letter"] for resi in peptideseq]) + bond_def
        except KeyError:  # name without defined letters
            name = ",".join(peptideseq) + bond_def

    mymol = Chem.MolFromSmiles(smiles)
    if not mymol:
        raise SmilesError("%s returns None molecule" % smiles)
    mymol.SetProp("_Name", name)

    if write == "draw":
        if not path.exists(twodfolder):
            os.makedirs(twodfolder)
        # Writes 2D drawing for reference
        # try:
        # mymol.draw(False, path.join(twodfolder, name+'.png'))
        # except NameError:
        # if not oasaAvailable:
        # print 'OASA drawing library required to draw 2D molecules'
        # raise
        try:
            # print name
            Draw.MolToFile(
                mymol, path.join(twodfolder, name + ".png"), size=(1000, 1000)
            )
        except NameError:
            print("rdkit name_error")

    elif write == "structure":
        if not path.exists(threedfolder):
            os.makedirs(threedfolder)
        # Makes 3D by adding hydrogens, writing coordinates and doing

        AllChem.EmbedMolecule(mymol)
        AllChem.UFFOptimizeMolecule(mymol)
        # mymol = AllChem.AddHs(mymol)
        # #50 step local optimization usinsg MMFF94
        # mymol.make3D()
        # mymol.addh()
        # Overwrite = True
        if not return_struct:
            try:
                file = path.join(threedfolder, name + "." + type)
                handle = open(file, "wb")
            except IOError:
                # Messy, stupid bug.
                # PRN.XXX files are reserved as devices in windows -
                # so PRN.pdb is illegal
                # Will raise exception - so change name t work around
                file = path.join(threedfolder, name + "_." + type)
                handle = open(file, "wb")
            try:
                out = Chem.MolToMolBlock(mymol)
                handle.write(out)
            finally:
                handle.close()
        else:
            return Chem.MolToMolBlock(mymol)
    else:
        raise TypeError('"write" must be set to "draw" or "structure", got %s' % write)

    return True


def write_library(inputlist, outloc, write="text", minimise=False, write_to_file=False):
    """
    Takes in a list of tuples of (peptide_name, bond_def, SMILES) and writes
    them in the specified format.
    'text' is a text file of names, bond_defs and files to outloc
    here, outloc is a file
    'draw' will write 2D representations of all peptides to the outloc
    'structure' will write 3D pdb files to the outloc

    Returns number of peptides written
    """
    count = 0
    if write == "text":
        try:
            f = open(outloc, "w")
            for peptide in inputlist:
                try:
                    seq, bond_def, smiles = peptide
                    if bond_def == "":
                        bond_def = "linear"
                    # print("%s-%s: %s" % (",".join(seq), bond_def, smiles), file=f)
                    count += 1
                except Exception as e:
                    print(e)
        finally:
            f.close()
    elif write == "draw" or write == "structure":
        if write_to_file:  # Checking if we need to write to a file
            with open(outloc, "w") as out:  # Ensure this is "w" if you're writing text, "wb" if binary
                for peptide in inputlist:
                    peptideseq, bond_def, smiles = peptide
                    if not peptideseq and not bond_def and not smiles:
                        continue
                    mol = Chem.MolFromSmiles(smiles)
                    AllChem.EmbedMolecule(mol)
                    AllChem.UFFOptimizeMolecule(mol)
                    try:
                        name = peptideseq + bond_def
                    except TypeError:  # If peptideseq or bond_def is not a string
                        try:
                            name = "".join([aminos[resi]["Letter"] for resi in peptideseq]) + bond_def
                        except KeyError:  # name without defined letters
                            name = ",".join(map(str, peptideseq)) + bond_def
                    mol.SetProp("_Name", name)
                    molstr = Chem.MolToMolBlock(mol)
                    out.write(molstr + '\n$$$$\n')  # Writing directly to the file
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
        raise TypeError('"write" must be set to "draw" or "structure", got %s' % write)
    return count


def main(pattern, out_file):
    """
    If running module as a script, take the input as a pattern to generate
    combinations for, and ultimately to write to a single .sdf file
    Will use all possible constraints
    """
    print("Writing all peptides for pattern %s" % pattern)
    out_f = out_file+".sdf"
    peptides = gen_all_matching_peptides(pattern)
    structures = gen_structs_from_seqs(peptides, True, True, True, True, True, True)
    write_library(structures, out_f, "structure", False, True)

if __name__ == "__main__":
    main(*sys.argv[1:], sys.argv[0])
