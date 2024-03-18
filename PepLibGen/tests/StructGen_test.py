#!/usr/bin/env python
"""
test program for StructGen/StructGen module in PepLibGen package

"""
import sys
import unittest
import math
import random

from PepLibGen.StructGen import StructGen as sg
from PepLibGen.StructGen.aminoacids import aminos  # From my package
from PepLibGen.StructGen.aminoacids import special_aminos
from PepLibGen.StructGen.aminoacids import d_aminos
from PepLibGen.StructGen import synthrules
from PepLibGen.Analysis import chem_analysis


class StructGenTestCases(unittest.TestCase):
    def testAddAmino(self):
        lib_len = len(sg.print_included_aminos())
        for amino in ["D-Lysine", "L-Ornithine", "D-Glutamic_Acid"]:
            sg.add_amino(amino)
            lib_len += 1
            self.assertTrue(len(sg.print_included_aminos()) == lib_len)

        for not_amino in ["ASDDF", "lysine", "2"]:
            self.assertRaises(sg.UndefinedAminoError, sg.add_amino, not_amino)

    def testRemoveAmino(self):
        lib_len = len(sg.print_included_aminos())
        for amino in ["D-Lysine", "L-Ornithine", "D-Glutamic_Acid"]:
            sg.remove_amino(amino)
            lib_len -= 1
            self.assertTrue(len(sg.print_included_aminos()) == (lib_len))

        for amino in ["D-Glutamine", "D-Serine"]:
            self.assertRaises(sg.UndefinedAminoError, sg.remove_amino, amino)

    def testPropertytoName(self):
        conversions = {
            "L-Serine": ("Letter", "S"),
            "L-Glutamine": ("Code", "Gln"),
            "L-Isoleucine": ("SMILES", "N[C@@]([H])(C(CC)C)C(=O)O"),
        }

        bad_properties = [("Notreal", "Z"), ("Notapropery", 10000000)]

        bad_conversions = [("Letter", "Z"), ("Code", "QQQ")]

        for name, test in list(conversions.items()):
            self.assertTrue(sg.property_to_name(*test) == name)

        for property, value in bad_properties:
            self.assertRaises(
                sg.UndefinedPropertyError, sg.property_to_name, *(property, value)
            )

        for property, values in bad_conversions:
            self.assertRaises(
                sg.UndefinedAminoError, sg.property_to_name, *(property, value)
            )

    def testGenAllPosPeptides(self):
        for num in range(2, 5):
            library = list(sg.gen_all_pos_peptides(num))
            num_peps = len(library)
            num_no_dups = len(set(library))

            self.assertTrue(
                num_peps == len(sg.aminodata) ** num, "Incorrect combinations... "
            )
            # Check for duplicates - sets can't contain duplicates, so...
            self.assertTrue(num_no_dups == num_peps, "Duplicate Peptides")

        bad_iterator = sg.gen_all_pos_peptides("3a")
        self.assertRaises(TypeError, bad_iterator.__next__)

    def testGenAllMatchingPeptides(self):
        testpatterns = [
            "CXXXC",
            "SSSX",
            "DXDX",
            "YYXYXX",
            ["L-Cysteine", "X", "X"],
            ["L-Glutamic_Acid", "X", "X"],
        ]
        for pattern in testpatterns:
            library = list(sg.gen_all_matching_peptides(pattern))
            num_peps = len(library)

            self.assertTrue(
                num_peps == len(sg.aminodata) ** pattern.count("X"),
                "Incorrect combinations",
            )
            # Check for duplicates - sets can't contain duplicates, so...
            peps = []
            for pep in sg.gen_all_matching_peptides(pattern):
                # Lists not hashable, so create string representation
                peps.append("-".join([str(num) for num in pep]))
            self.assertTrue(len(set(peps)) == num_peps, "Duplicate peptides")

        bad_patterns = ["12345", "^$%"]
        bad_type = [23, 450]
        for pattern in bad_patterns:
            test = sg.gen_all_matching_peptides(pattern)
            self.assertRaises(sg.UndefinedAminoError, test.__next__)
        for pattern in bad_type:
            test = sg.gen_all_matching_peptides(pattern)
            self.assertRaises(TypeError, test.__next__)

    def testCanSSBond(self):
        truepatterns = [
            "CAAC",
            "CCCC",
            "HFHFCHCHHHC",
            ["L-Cysteine", "L-Alanine", "L-Glutamine", "L-Cysteine"],
            ["L-Cysteine", "L-Leucine", "L-Arginine", "L-Cysteine", "L-Alanine"],
        ]
        falsepatterns = [
            "CACA",
            "HDET",
            ["L-Cysteine", "L-Threonine", "L-Cysteine", "L-Tryptophan"],
            ["L-Asparagine", "Glycine", "L-Alanine", "L-Histidine"],
            [
                "L-Arginine",
                "L-Glutamic_Acid",
                "L-Proline",
                "L-Glutamic_Acid",
                "L-Proline",
            ],
        ]
        badpatterns = ["CBC", ["L-Notexistene", "D-Fakene"]]
        for pattern in truepatterns:
            self.assertTrue(
                sg.can_ssbond(pattern) != False and sg.can_ssbond(pattern) != None,
                "%s cannot ssbond" % pattern,
            )
        for pattern in falsepatterns:
            self.assertTrue(sg.can_ssbond(pattern) == False, "%s can ssbond" % pattern)
        for pattern in badpatterns:
            self.assertRaises(sg.UndefinedAminoError, sg.can_ssbond, pattern)

    def testCanHTBond(self):
        truepatterns = [
            "CAAFC",
            "FCCCC",
            "HFHFCHCHHHC",
            ["L-Cysteine", "L-Leucine", "L-Arginine", "L-Cysteine", "L-Alanine"],
            ["L-Cysteine", "L-Serine", "L-Arginine", "L-Cysteine", "L-Alanine"],
        ]
        falsepatterns = [
            "CACA",
            "HDET",
            "FG",
            ["L-Phenylalanine", "L-Glutamine"],
            ["L-Isoleucine", "L-Asparagine"],
        ]
        badpatterns = ["CBC", ["L-Notexistene", "D-Fakene"]]
        for pattern in truepatterns:
            self.assertTrue(
                sg.can_htbond(pattern) != False and sg.can_htbond(pattern) != None
            )
        for pattern in falsepatterns:
            self.assertTrue(sg.can_htbond(pattern) == False)
        # Will not test that aminos are valid

    def testCanSCNTBond(self):
        truepatterns = [
            "ASFD",
            "HRKD",
            "LMNE",
            "DEDE",
            "ADAE",
            ["L-Cysteine", "L-Leucine", "L-Arginine", "L-Glutamic_Acid", "L-Alanine"],
        ]
        falsepatterns = [
            "AAAA",
            "ADAA",
            "SYTA",
            "IEWPVIWPV",
            ["L-Cysteine", "L-Leucine", "L-Arginine", "L-Cysteine", "L-Alanine"],
        ]
        badpatterns = [
            "CBCB",
            ["L-Notexistene", "D-Fakene", "L-Notexistene", "D-Fakene"],
        ]
        for pattern in truepatterns:
            self.assertTrue(
                sg.can_scntbond(pattern) != False and sg.can_scntbond(pattern) != None,
                "%s cannot SCNT bond" % pattern,
            )
        for pattern in falsepatterns:
            self.assertTrue(
                sg.can_scntbond(pattern) == False, "%s can SCNT bond" % pattern
            )
        for pattern in badpatterns:
            self.assertRaises(sg.UndefinedAminoError, sg.can_scntbond, pattern)

    def testCanSCCTbond(self):
        # Can bond K,S,Y,T uniquely
        truepatterns = ["KHRK", "SMNS", "YQGC", "TYEN"]
        falsepatterns = ["FLWK", "GRCNS", "FLWANPY", "RFACGT"]
        badpatterns = [
            "CBCBC",
            ["L-Notexistene", "D-Fakene", "L-Notexistene", "D-Fakene"],
        ]
        for pattern in truepatterns:
            self.assertTrue(
                sg.can_scctbond(pattern) != False and sg.can_scctbond(pattern) != None,
                "%s cannot SCCT bond" % pattern,
            )
        for pattern in falsepatterns:
            self.assertTrue(
                sg.can_scctbond(pattern) == False, "%s can SCCT bond" % pattern
            )
        for pattern in badpatterns:
            self.assertRaises(sg.UndefinedAminoError, sg.can_scctbond, pattern)

    def testCanSCSCBond(self):
        # Can bond K,S,T to D,E
        truepatterns = ["KAAD", "DAAK", "SASGE"]
        falsepatterns = ["KQGQG", "A"]
        badpatterns = [
            "CBCB",
            ["L-Notexistene", "D-Fakene", "L-Notexistene", "D-Fakene"],
        ]
        for pattern in truepatterns:
            self.assertTrue(
                sg.can_scscbond(pattern) != False and sg.can_scscbond(pattern) != None,
                "%s cannot SCSC bond" % pattern,
            )
        for pattern in falsepatterns:
            self.assertTrue(
                sg.can_scscbond(pattern) == False, "%s can SCSC bond" % pattern
            )
        for pattern in badpatterns:
            self.assertRaises(sg.UndefinedAminoError, sg.can_scscbond, pattern)

    def testReturnSmiles(self):
        # test letter and 3 letter codes
        testdic = {
            "A": "N[C@@]([H])(C)C(=O)O",
            "W": "N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)O",
            "Val": "N[C@@]([H])(C(C)C)C(=O)O",
            "Gln": "N[C@@]([H])(CCC(=O)N)C(=O)O",
        }
        for code, smi in list(testdic.items()):
            self.assertTrue(sg.return_smiles(code) == smi)

    def testLinearPeptideSmiles(self):
        testdic = {
            "AAAA": "N[C@@]([H])(C)C(=O)N[C@@]([H])(C)C(=O)\
N[C@@]([H])(C)C(=O)N[C@@]([H])(C)C(=O)O",
            "AGAG": "N[C@@]([H])(C)C(=O)NCC(=O)N[C@@]([H])(C)C(=O)NCC(=O)O",
            "PQE": "N1[C@@]([H])(CCC1)C(=O)N[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])\
(CCC(=O)O)C(=O)O",
            "CNVR": "N[C@@]([H])(CS)C(=O)N[C@@]([H])(C(=O)N)C(=O)\
N[C@@]([H])(C(C)C)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)O",
        }
        for seq, smi in list(testdic.items()):
            self.assertTrue(
                sg.linear_peptide_smiles(seq) == smi,
                "Wrong smiles for %s: saw %s, should be %s"
                % (seq, sg.linear_peptide_smiles(seq), smi),
            )

    def testBondCounter(self):
        testdic = {
            "N1[C@@]([H])(C)C(=O)N[C@@]([H])(C)C(=O)N[C@@]([H])\
		(C)C(=O)N[C@@]([H])(C)C(=O)1": 1,
            "N[C@@]([H])(C)C(=O)NCC(=O)N[C@@]([H])(C)C(=O)NCC(=O)O": 0,
            "N9[C@@]([H])(CCC9)C(=O)N[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])\
		(CCC(=O)O)C(=O)O": 9,
            "N[C@@]([H])(CS)C(=O)N[C@@]([H])(C(=O)N)C(=O)N[C@@]([H])\
		(C(C)C)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)O": 0,
            "N[C@@]([H])(CS5)C(=O)N[C@@]([H])(C(=O)N)C(=O)N[C@@]([H])(C(C)C)C(=O)\
		N[C@@]([H])(CS5)C(=O)O": 5,
        }
        for smi, bondcount in list(testdic.items()):
            self.assertTrue(
                sg.bond_counter(smi) == bondcount,
                "Bondcounter: %s should be %s" % (sg.bond_counter(smi), bondcount),
            )

    def testConstrainedPeptide(self):
        testdic = {
            (
                "CAAC",
                "SSCXXC",
            ): "N[C@@]([H])(CS1)C(=O)N[C@@]([H])(C)C(=O)N[C@@]([H])(C)C(=O)N[C@@]([H])(CS1)C(=O)O",
            (
                "CACAC",
                "SSXXCXC",
            ): "N[C@@]([H])(CS)C(=O)N[C@@]([H])(C)C(=O)N[C@@]([H])(CS1)C(=O)N[C@@]([H])(C)C(=O)N[C@@]([H])(CS1)C(=O)O",
            (
                "CWCCWC",
                "SSXXCXXC",
            ): "N[C@@]([H])(CS)C(=O)N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)N[C@@]([H])(CS3)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)N[C@@]([H])(CS3)C(=O)O",
            (
                "CAAC",
                "HT",
            ): "N1[C@@]([H])(CS)C(=O)N[C@@]([H])(C)C(=O)N[C@@]([H])(C)C(=O)N[C@@]([H])(CS)C(=O)1",
            (
                "CACAC",
                "HT",
            ): "N1[C@@]([H])(CS)C(=O)N[C@@]([H])(C)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(C)C(=O)N[C@@]([H])(CS)C(=O)1",
            (
                "CWCCWC",
                "HT",
            ): "N3[C@@]([H])(CS)C(=O)N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)N[C@@]([H])(CS)C(=O)3",
            (
                "KPVEE",
                "SCNXXXZ",
            ): "N[C@@]([H])(CCCCN2)C(=O)N1[C@@]([H])(CCC1)C(=O)N[C@@]([H])(C(C)C)C(=O)N[C@@]([H])(CCC(=O)O)C(=O)N[C@@]([H])(CCC2(=O))C(=O)O",
            (
                "SHRDK",
                "SCEXXZX",
            ): "N[C@@]([H])(CO2)C(=O)N[C@@]([H])(CC1=CN=C-N1)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CC2(=O))C(=O)N[C@@]([H])(CCCCN)C(=O)O",
            (
                "AELMK",
                "SCXZXXN",
            ): "N[C@@]([H])(C)C(=O)N[C@@]([H])(CCC1(=O))C(=O)N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(CCSC)C(=O)N[C@@]([H])(CCCCN1)C(=O)O",
            (
                "SKYTE",
                "SCXNXXX",
            ): "N[C@@]([H])(CO)C(=O)N[C@@]([H])(CCCCN2)C(=O)N[C@@]([H])(Cc1ccc(O)cc1)C(=O)N[C@@]([H])(C(O)C)C(=O)N[C@@]([H])(CCC(=O)O)C(=O)2",
            (
                "QGCAE",
                "SCXXXXZ",
            ): "N1[C@@]([H])(CCC(=O)N)C(=O)NCC(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(C)C(=O)N[C@@]([H])(CCC1(=O))C(=O)O",
        }
        for (seq, bond_def), smi in list(testdic.items()):
            self.assertTrue(
                sg.constrained_peptide_smiles(seq, bond_def)[2] == smi,
                "Wrong smiles for %s: is %s should be %s"
                % (seq, sg.constrained_peptide_smiles(seq, bond_def)[2], smi),
            )

    def testWhatConstraints(self):
        tests = [
            ("CAAC", ["SSCXXC"]),
            ("AAGAA", ["HT"]),
            (
                ["L-Cysteine", "L-Serine", "L-Arginine", "L-Cysteine", "L-Alanine"],
                ["SSCXXCX", "HT", "SCXEXXX"],
            ),
        ]
        for seq, constraints in tests:

            self.assertTrue(sg.what_constraints(seq) == constraints)

    def testPepPositions(self):
        teststrings = []
        for k in range(3, 10):
            resis = random.sample(sg.return_available_residues(), k)
            random.shuffle(resis)
            teststrings.append("".join(resis))
        for seq in teststrings:
            locs = sg.pep_positions(seq)
            smi = sg.linear_peptide_smiles(seq)
            for loc in locs:
                self.assertTrue(
                    smi[loc] == "N",
                    "Seq %s loc %s in locs %s is %s" % (seq, loc, locs, smi[loc]),
                )

    # def testGenLibraryStrings(self):
    # #Generate library of linear peptides, and figure out how many can bind
    # #in what way.
    # #Then check that the library generated by sg.gen_library_strings
    # #matches
    # #Input = (liblen, ssbond=False,htbond=False,scctbond=False,
    # #scntbond=False,scscbond=False,linear=False)

    # for liblen in range(3,5):
    # allpeps = list(sg.gen_all_pos_peptides(liblen))
    # num_ssbonds = 0
    # num_scntbonds = 0
    # num_scctbonds = 0
    # num_scscbonds = 0
    # num_htbonds = 0
    # num_linear = len(allpeps)
    # for pep in allpeps:
    # if sg.can_ssbond(pep): num_ssbonds += 1
    # if sg.can_htbond(pep): num_htbonds += 1
    # if sg.can_scntbond(pep): num_scntbonds += 1
    # if sg.can_scctbond(pep): num_scctbonds += 1
    # if sg.can_scscbond(pep): num_scscbonds += 1
    # #Repeat 5 times, for random inputs
    # for num in range(0,5):
    # ssbond = random.randint(0,1)
    # htbond = random.randint(0,1)
    # scctbond = random.randint(0,1)
    # scntbond =random.randint(0,1)
    # scscbond =random.randint(0,1)
    # linear = random.randint(0,1)
    # testpeps = list(sg.gen_library_strings(liblen,ssbond, htbond, \
    # scctbond, scntbond, scscbond, linear))
    # expected_lib_len = num_ssbonds*ssbond + num_htbonds*htbond + \
    # num_scntbonds*scntbond \
    # + num_scctbonds*scctbond + num_scscbonds*scscbond \
    # + num_linear*linear
    # self.assertTrue(len(testpeps)==expected_lib_len, \
    # 'Peptide mismatch - got %s total, with ssbond %s, htbond %s, \
    # scctbond %s, scntbond %s, scscbond %s, should be %s with \
    # liblen %s ' % \
    # (len(testpeps), ssbond, htbond, scctbond, scntbond, \
    # scscbond, expected_lib_len, liblen))

    def testStructsFromSeqs(self):
        # gen_structs_from_seqs(sequences, ssbond=False,htbond=False,scctbond=False,
        # scntbond=False,scscbond=False,linear=False):
        seqs = [
            "ASDF",
            "QFDSA",
            "CRQSC",
            "DSCCK",
            ["L-Cysteine", "L-Leucine", "L-Arginine", "L-Cysteine", "L-Alanine"],
            ["L-Alanine", "L-Serine", "L-Arginine", "L-Cysteine", "L-Glutamic_Acid"],
        ]

        results = [
            (
                "ASDF",
                "",
                "N[C@@]([H])(C)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(Cc1ccccc1)C(=O)O",
            ),
            (
                "QFDSA",
                "HT",
                "N2[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(Cc1ccccc1)C(=O)N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(C)C(=O)2",
            ),
            (
                "QFDSA",
                "",
                "N[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(Cc1ccccc1)C(=O)N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(C)C(=O)O",
            ),
            (
                "CRQSC",
                "SSCXXXC",
                "N[C@@]([H])(CS1)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS1)C(=O)O",
            ),
            (
                "CRQSC",
                "HT",
                "N1[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS)C(=O)1",
            ),
            (
                "CRQSC",
                "",
                "N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS)C(=O)O",
            ),
            (
                "DSCCK",
                "HT",
                "N1[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCCN)C(=O)1",
            ),
            (
                "DSCCK",
                "SCXEXXX",
                "N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(CO1)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCCN)C(=O)1",
            ),
            (
                "DSCCK",
                "SCZXXXN",
                "N[C@@]([H])(CC1(=O))C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCCN1)C(=O)O",
            ),
            (
                "DSCCK",
                "",
                "N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCCN)C(=O)O",
            ),
            (
                ["L-Cysteine", "L-Leucine", "L-Arginine", "L-Cysteine", "L-Alanine"],
                "SSCXXCX",
                "N[C@@]([H])(CS1)C(=O)N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS1)C(=O)N[C@@]([H])(C)C(=O)O",
            ),
            (
                ["L-Cysteine", "L-Leucine", "L-Arginine", "L-Cysteine", "L-Alanine"],
                "HT",
                "N1[C@@]([H])(CS)C(=O)N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(C)C(=O)1",
            ),
            (
                ["L-Cysteine", "L-Leucine", "L-Arginine", "L-Cysteine", "L-Alanine"],
                "",
                "N[C@@]([H])(CS)C(=O)N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(C)C(=O)O",
            ),
            (
                [
                    "L-Alanine",
                    "L-Serine",
                    "L-Arginine",
                    "L-Cysteine",
                    "L-Glutamic_Acid",
                ],
                "HT",
                "N1[C@@]([H])(C)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCC(=O)O)C(=O)1",
            ),
            (
                [
                    "L-Alanine",
                    "L-Serine",
                    "L-Arginine",
                    "L-Cysteine",
                    "L-Glutamic_Acid",
                ],
                "SCXEXXX",
                "N[C@@]([H])(C)C(=O)N[C@@]([H])(CO1)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCC(=O)O)C(=O)1",
            ),
            (
                [
                    "L-Alanine",
                    "L-Serine",
                    "L-Arginine",
                    "L-Cysteine",
                    "L-Glutamic_Acid",
                ],
                "SCXXXXZ",
                "N1[C@@]([H])(C)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCC1(=O))C(=O)O",
            ),
            (
                [
                    "L-Alanine",
                    "L-Serine",
                    "L-Arginine",
                    "L-Cysteine",
                    "L-Glutamic_Acid",
                ],
                "SCXEXXZ",
                "N[C@@]([H])(C)C(=O)N[C@@]([H])(CO1)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCC1(=O))C(=O)O",
            ),
            (
                [
                    "L-Alanine",
                    "L-Serine",
                    "L-Arginine",
                    "L-Cysteine",
                    "L-Glutamic_Acid",
                ],
                "",
                "N[C@@]([H])(C)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCC(=O)O)C(=O)O",
            ),
        ]

        trial = sg.gen_structs_from_seqs(
            seqs,
            ssbond=True,
            htbond=True,
            scctbond=True,
            scntbond=True,
            scscbond=True,
            linear=True,
        )
        for a, b in zip(results, trial):
            self.assertTrue(a == b, " real %s trial %s" % (a, b))

    # def testGenLibraryStructs(self):
    # #Make sure that the function returns the correct amount of
    # #name,smi pairs
    # for liblen in range(3,5):
    # for num in range(0,5):
    # ssbond = random.randint(0,1)
    # htbond = random.randint(0,1)
    # scctbond = random.randint(0,1)
    # scntbond =random.randint(0,1)
    # scscbond =random.randint(0,1)
    # linear = random.randint(0,1)
    # expected_num = len(list(sg.gen_library_strings(liblen,ssbond, htbond,\
    # scctbond, scntbond, scscbond, linear)))
    # result = list(sg.gen_library_structs(liblen,ssbond, htbond,\
    # scctbond, scntbond, scscbond, linear))
    # self.assertTrue(expected_num == len(result),\
    # 'Returned %s structs, from %s inputs' \
    # % (expected_num, len(result)))
    # for struct in result:
    # seq, bond_def, smi = struct
    # self.assertTrue(isinstance(struct, tuple) and\
    # (isinstance(seq, str) or isinstance(seq, list) \
    # or isinstance(seq, tuple)) \
    # and isinstance(smi, str), \
    # 'struct %s, seq %s'\
    # % (struct, seq))

    def testFilteredOutput(self):
        input = [
            (
                "ASDF",
                "",
                "N[C@@]([H])(C)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(Cc1ccccc1)C(=O)O",
            ),
            (
                "QFDSA",
                "HT",
                "N2[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(Cc1ccccc1)C(=O)N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(C)C(=O)2",
            ),
            (
                "QFDSA",
                "",
                "N[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(Cc1ccccc1)C(=O)N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(C)C(=O)O",
            ),
            (
                "CRQSC",
                "SSCXXXC",
                "N[C@@]([H])(CS1)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS1)C(=O)O",
            ),
            (
                "CRQSC",
                "HT",
                "N1[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS)C(=O)1",
            ),
            (
                "CRQSC",
                "",
                "N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS)C(=O)O",
            ),
            (
                "DSCCK",
                "HT",
                "N1[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCCN)C(=O)1",
            ),
            (
                "DSCCK",
                "SCXEXXX",
                "N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(CO1)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCCN)C(=O)1",
            ),
            (
                "DSCCK",
                "SCZXXXN",
                "N[C@@]([H])(CC1(=O))C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCCN1)C(=O)O",
            ),
            (
                "DSCCK",
                "",
                "N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCCN)C(=O)O",
            ),
            (
                ["L-Cysteine", "L-Leucine", "L-Arginine", "L-Cysteine", "L-Alanine"],
                "SSCXXCX",
                "N[C@@]([H])(CS1)C(=O)N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS1)C(=O)N[C@@]([H])(C)C(=O)O",
            ),
            (
                ["L-Cysteine", "L-Leucine", "L-Arginine", "L-Cysteine", "L-Alanine"],
                "HT",
                "N1[C@@]([H])(CS)C(=O)N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(C)C(=O)1",
            ),
            (
                ["L-Cysteine", "L-Leucine", "L-Arginine", "L-Cysteine", "L-Alanine"],
                "",
                "N[C@@]([H])(CS)C(=O)N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(C)C(=O)O",
            ),
            (
                [
                    "L-Alanine",
                    "L-Serine",
                    "L-Arginine",
                    "L-Cysteine",
                    "L-Glutamic_Acid",
                ],
                "HT",
                "N1[C@@]([H])(C)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCC(=O)O)C(=O)1",
            ),
            (
                [
                    "L-Alanine",
                    "L-Serine",
                    "L-Arginine",
                    "L-Cysteine",
                    "L-Glutamic_Acid",
                ],
                "SCXEXXX",
                "N[C@@]([H])(C)C(=O)N[C@@]([H])(CO1)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCC(=O)O)C(=O)1",
            ),
            (
                [
                    "L-Alanine",
                    "L-Serine",
                    "L-Arginine",
                    "L-Cysteine",
                    "L-Glutamic_Acid",
                ],
                "SCXXXXZ",
                "N1[C@@]([H])(C)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCC1(=O))C(=O)O",
            ),
            (
                [
                    "L-Alanine",
                    "L-Serine",
                    "L-Arginine",
                    "L-Cysteine",
                    "L-Glutamic_Acid",
                ],
                "SCXEXXZ",
                "N[C@@]([H])(C)C(=O)N[C@@]([H])(CO1)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCC1(=O))C(=O)O",
            ),
            (
                [
                    "L-Alanine",
                    "L-Serine",
                    "L-Arginine",
                    "L-Cysteine",
                    "L-Glutamic_Acid",
                ],
                "",
                "N[C@@]([H])(C)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCC(=O)O)C(=O)O",
            ),
        ]

        over_four_resis = [
            (
                "ASDF",
                "",
                "N[C@@]([H])(C)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(Cc1ccccc1)C(=O)O",
            )
        ]

        def trialfunc_one(seq):
            if len(seq) > 4:
                return False
            else:
                return True

        def trialfunc_two(result):
            if len(result[0]) > 4:
                return False
            else:
                return True

        for func in [trialfunc_one, trialfunc_two]:
            if func is trialfunc_one:
                key = 0
            else:
                key = None
            filtered_results = sg.filtered_output(input, func, key)
            for output, test_case in zip(filtered_results, over_four_resis):
                self.assertTrue(
                    output == test_case, "output %s test %s" % (output, test_case)
                )

    def testGetConstraintType(self):
        tests = [
            ("HT", "HT"),
            ("SSCXXC", "SS"),
            ("SSXCXXC", "SS"),
            ("SCZXXE", "SCSC"),
            ("SCEXXXX", "SCCT"),
            ("SCXNXXX", "SCCT"),
            ("SCXXXZX", "SCNT"),
        ]

        errors = ["SCZXXD", "SSZXXD", "SCCXXC", "SCCXXZ"]

        for bond_def, constraint_type in tests:
            self.assertTrue(sg.get_constraint_type(bond_def) == constraint_type)

        for bond_def in errors:
            self.assertRaises(sg.BondSpecError, sg.get_constraint_type, bond_def)

    def testCountConstraintTypes(self):
        input = [
            (
                "ASDF",
                "",
                "N[C@@]([H])(C)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(Cc1ccccc1)C(=O)O",
            ),
            (
                "QFDSA",
                "HT",
                "N2[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(Cc1ccccc1)C(=O)N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(C)C(=O)2",
            ),
            (
                "QFDSA",
                "",
                "N[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(Cc1ccccc1)C(=O)N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(C)C(=O)O",
            ),
            (
                "CRQSC",
                "SSCXXXC",
                "N[C@@]([H])(CS1)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS1)C(=O)O",
            ),
            (
                "CRQSC",
                "HT",
                "N1[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS)C(=O)1",
            ),
            (
                "CRQSC",
                "",
                "N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS)C(=O)O",
            ),
            (
                "DSCCK",
                "HT",
                "N1[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCCN)C(=O)1",
            ),
            (
                "DSCCK",
                "SCXEXXX",
                "N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(CO1)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCCN)C(=O)1",
            ),
            (
                "DSCCK",
                "SCZXXXN",
                "N[C@@]([H])(CC1(=O))C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCCN1)C(=O)O",
            ),
            (
                "DSCCK",
                "",
                "N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCCCN)C(=O)O",
            ),
            (
                ["L-Cysteine", "L-Leucine", "L-Arginine", "L-Cysteine", "L-Alanine"],
                "SSCXXCX",
                "N[C@@]([H])(CS1)C(=O)N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS1)C(=O)N[C@@]([H])(C)C(=O)O",
            ),
            (
                ["L-Cysteine", "L-Leucine", "L-Arginine", "L-Cysteine", "L-Alanine"],
                "HT",
                "N1[C@@]([H])(CS)C(=O)N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(C)C(=O)1",
            ),
            (
                ["L-Cysteine", "L-Leucine", "L-Arginine", "L-Cysteine", "L-Alanine"],
                "",
                "N[C@@]([H])(CS)C(=O)N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(C)C(=O)O",
            ),
            (
                [
                    "L-Alanine",
                    "L-Serine",
                    "L-Arginine",
                    "L-Cysteine",
                    "L-Glutamic_Acid",
                ],
                "HT",
                "N1[C@@]([H])(C)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCC(=O)O)C(=O)1",
            ),
            (
                [
                    "L-Alanine",
                    "L-Serine",
                    "L-Arginine",
                    "L-Cysteine",
                    "L-Glutamic_Acid",
                ],
                "SCXEXXX",
                "N[C@@]([H])(C)C(=O)N[C@@]([H])(CO1)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCC(=O)O)C(=O)1",
            ),
            (
                [
                    "L-Alanine",
                    "L-Serine",
                    "L-Arginine",
                    "L-Cysteine",
                    "L-Glutamic_Acid",
                ],
                "SCXXXXZ",
                "N1[C@@]([H])(C)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCC1(=O))C(=O)O",
            ),
            (
                [
                    "L-Alanine",
                    "L-Serine",
                    "L-Arginine",
                    "L-Cysteine",
                    "L-Glutamic_Acid",
                ],
                "SCXEXXZ",
                "N[C@@]([H])(C)C(=O)N[C@@]([H])(CO1)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCC1(=O))C(=O)O",
            ),
            (
                [
                    "L-Alanine",
                    "L-Serine",
                    "L-Arginine",
                    "L-Cysteine",
                    "L-Glutamic_Acid",
                ],
                "",
                "N[C@@]([H])(C)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(CCC(=O)O)C(=O)O",
            ),
        ]

        count_dict = {"linear": 6, "SS": 2, "HT": 5, "SCSC": 2, "SCCT": 2, "SCNT": 1}

        test_count_dict = sg.count_constraint_types(input)
        self.assertTrue(
            count_dict == test_count_dict,
            "true %s test %s" % (count_dict, test_count_dict),
        )


if __name__ == "__main__":
    unittest.main()
