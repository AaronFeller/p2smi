#!/usr/bin/env python

import unittest
from PepLibGen.StructGen import synthrules


class StructGenTestCases(unittest.TestCase):
    bad_motifs = ["ASDPPP", "QNASD", "DGAS", "DPAS", "NDFG"]
    good_motifs = ["ASDF", "WEERT", "WSED"]

    mods = [
        (("SAAD", "EXXZ"), "GAAN"),
        (("KYYE", "NXXZ"), "LYYQ"),
        (("CDEC", "CXXC"), "MDEM"),
        (("YAAD", "SXXZ"), "FAAN"),
        (("TRKD", "EXXZ"), "VRKN"),
    ]

    hydrophobicities = {"ASDF": 0.73, "QWER": -3.31}

    pass_charges = ["ASDK", "QWREW", "AAAKAAA", "IWPV"]
    fail_charges = ["IWPVI", "NSYTI"]

    def test_motifs(self):
        for motif in self.bad_motifs:
            self.assertTrue(synthrules.test_motifs(motif))
        for motif in self.good_motifs:
            self.assertEqual(synthrules.test_motifs(motif), [])

    def test_modseq(self):
        for input, true_output in self.mods:
            self.assertEqual(synthrules.mod_seq(*input), true_output)

    def test_hydrophobicity(self):
        for sequence in self.hydrophobicities:
            self.assertAlmostEqual(
                synthrules.test_hydrophobicity(sequence),
                self.hydrophobicities[sequence],
            )

    def test_charge(self):
        for sequence in self.pass_charges:
            self.assertTrue(synthrules.test_charge(sequence))
        for sequence in self.fail_charges:
            self.assertFalse(synthrules.test_charge(sequence))


if __name__ == "__main__":
    unittest.main()
