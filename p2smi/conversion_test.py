import os
import unittest
from unittest.mock import MagicMock
from rdkit import Chem
from generate import fasta2smi
from generate.generate import generate_sequence


class TestSmilesValidation(unittest.TestCase):
    def setUp(self):
        # Create a directory for test files
        self.test_dir = "test_files"
        os.makedirs(self.test_dir, exist_ok=True)

        # Generate a test fasta file
        self.fasta_file = os.path.join(self.test_dir, "unit_test.fasta")
        self.smi_file = os.path.join(self.test_dir, "unit_test.smi")
        generate_sequence(5, 10, 0.2, 0.2, 2, None, outfile=self.fasta_file)

    def tearDown(self):
        # Remove test files and directory
        if os.path.exists(self.fasta_file):
            os.remove(self.fasta_file)
        if os.path.exists(self.smi_file):
            os.remove(self.smi_file)
        if os.path.exists(self.test_dir):
            os.rmdir(self.test_dir)

    def test_smiles_conversion_and_validation(self):
        # Mock the argument parser
        args = MagicMock()
        args.input_fasta = self.fasta_file
        args.out_file = self.smi_file

        # Call the main function
        fasta2smi.main(args)

        # Read the SMILES output
        with open(self.smi_file, "r") as f:
            smiles = f.read().strip().splitlines()

        # Validate SMILES
        for smi in smiles:
            smi = smi.split(': ')[-1]
            mol = Chem.MolFromSmiles(smi)
            self.assertIsNotNone(mol, f"Invalid SMILES: {smi}")


if __name__ == "__main__":
    unittest.main()