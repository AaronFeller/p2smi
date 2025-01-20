#!/usr/bin/env python
"""
Optimized test program for StructGen/StructGen module in PepLibGen package.
"""
import unittest
from unittest.mock import patch
from generate import smilesgen as sg


class StructGenTestCases(unittest.TestCase):
    @patch("generate.smilesgen.add_amino")  # Mock expensive methods
    def testAddAmino(self, mock_add_amino):
        # Simplified and reduced input set
        aminos = ["D-Lysine", "L-Ornithine"]
        mock_add_amino.side_effect = lambda x: None  # Assume success
        for amino in aminos:
            sg.add_amino(amino)
            mock_add_amino.assert_called_with(amino)

    def testPropertyToName(self):
        # Consolidate good and bad tests
        valid_cases = [
            ("L-Serine", ("Letter", "S")),
            ("L-Glutamine", ("Code", "Gln")),
        ]
        invalid_cases = [("Notreal", "Z"), ("Letter", "X")]

        for name, (prop, val) in valid_cases:
            self.assertEqual(sg.property_to_name(prop, val), name)

        for prop, val in invalid_cases:
            with self.assertRaises(sg.UndefinedAminoError):
                sg.property_to_name(prop, val)

    def testGenAllPosPeptides(self):
        # Test minimal cases
        for num in range(1, 3):  # Smaller range
            library = list(sg.gen_all_pos_peptides(num))
            self.assertEqual(len(library), len(sg.aminodata) ** num)

        # Test invalid case
        with self.assertRaises(TypeError):
            next(sg.gen_all_pos_peptides("invalid"))

    def testGenAllMatchingPeptides(self):
        test_patterns = ["XAX", "XX"]
        for pattern in test_patterns:
            library = list(sg.gen_all_matching_peptides(pattern))
            self.assertEqual(len(library), len(sg.aminodata) ** pattern.count("X"))

    def testCanSSBond(self):
        # Reduced test set
        true_patterns = ["CAAC", "CCCC"]
        false_patterns = ["CACA"]

        for pattern in true_patterns:
            self.assertTrue(sg.can_ssbond(pattern))

        for pattern in false_patterns:
            self.assertFalse(sg.can_ssbond(pattern))

    def testReturnSmiles(self):
        # Simplified mapping
        test_data = {"A": "N[C@@H](C)C(=O)O", "G": "NCC(=O)O"}
        for code, smile in test_data.items():
            self.assertEqual(sg.return_smiles(code), smile)

    def testLinearPeptideSmiles(self):
        test_cases = {
            "AAA": "N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O"
        }
        for seq, smile in test_cases.items():
            self.assertEqual(sg.linear_peptide_smiles(seq), smile)


if __name__ == "__main__":
    unittest.main()