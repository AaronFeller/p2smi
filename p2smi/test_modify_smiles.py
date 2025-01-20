#!/usr/bin/env python
"""
Optimized test program for StructGen/StructGen module in PepLibGen package.
"""
import os
import unittest
from unittest.mock import patch
from generate import smilesgen as sg
from generate.modify_smiles import process_file as modify_smiles


class StructGenTestCases(unittest.TestCase):
    
    def testModifySmilesFile(self):
        # Create temporary input and output files
        input_file = "test_smiles_input.txt"
        output_file = "test_smiles_output.txt"
        with open(input_file, "w") as f:
            f.write("Seq1: N[C@@H](CN)C(=O)O\n")
            f.write("Seq2: N[C@@H](Cc1ccc(cc1)C)C(=O)N[C@@H](CC1=CN=C-N1)C(=O)O\n")

        try:
            # Call the modify_smiles function with both modifications
            modify_smiles(
                input_file=input_file,
                output_file=output_file,
                apply_methylation=True,
                apply_pegylation=True,
                methylation_percent=0.2,
            )

            # Validate the output file
            with open(output_file, "r") as f:
                lines = f.readlines()

            self.assertIn("N-methylation", lines[0], "Expected N-methylation annotation in header.")
            self.assertIn("PEGylation", lines[0], "Expected PEGylation annotation in header.")
            self.assertIn("N-methylation", lines[1], "Expected N-methylation annotation in header.")
            self.assertIn("PEGylation", lines[1], "Expected PEGylation annotation in header.")

        finally:
            # Clean up temporary files
            if os.path.exists(input_file):
                os.remove(input_file)
            if os.path.exists(output_file):
                os.remove(output_file)


if __name__ == "__main__":
    unittest.main()