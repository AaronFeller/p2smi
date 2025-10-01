#!/bin/bash

#!/usr/bin/env bash

# record start
TOTAL_START=$(date +%s)

echo "Running generate-peptides..."
START=$(date +%s)
generate-peptides -min 100 --max 100 -ncaa 0.1 -d 0.1 -n 1000 -c all -o peptides.fasta
END=$(date +%s)
echo "generate-peptides took $((END - START)) seconds"

echo "Running fasta2smi..."
START=$(date +%s)
fasta2smi -i peptides.fasta -o peptides.p2smi
END=$(date +%s)
echo "fasta2smi took $((END - START)) seconds"

echo "Running modify-smiles..."
START=$(date +%s)
modify-smiles -i peptides.p2smi -o modified.p2smi --peg_rate 0.3 --nmeth_rate 0.2 --nmeth_residues 0.25
END=$(date +%s)
echo "modify-smiles took $((END - START)) seconds"

echo "Running smiles-props..."
START=$(date +%s)
smiles-props -i modified.p2smi -o properties.p2smi
END=$(date +%s)
echo "smiles-props took $((END - START)) seconds"

TOTAL_END=$(date +%s)
echo "Total pipeline time: $((TOTAL_END - TOTAL_START)) seconds"

# Clean up intermediate files
rm peptides.fasta peptides.p2smi modified.p2smi properties.p2smi