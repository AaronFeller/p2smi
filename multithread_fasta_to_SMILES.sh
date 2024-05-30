#!/bin/bash

# Directory containing the split FASTA files
input_dir="split_files"

# Output directory for the SMILES files
output_dir="output_smiles"
mkdir -p "$output_dir"

# Function call format
script_call="python3 PepLibGen/StructGen/fasta.py --out"

# Export the function call to be available in the parallel environment
export script_call

# Function to process each file
process_file() {
    file=$1
    base_name=$(basename "$file" .fasta)
    output_file="output_smiles/${base_name}_SMILES.txt"
    $script_call "$output_file" "$file"
}

export -f process_file

# Find all fasta files in the input directory and pass them to parallel
find "$input_dir" -name "*.fasta" | parallel -j 10 process_file {}
