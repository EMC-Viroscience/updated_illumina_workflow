#!/usr/bin/env python3

"""
Script: multiL_fasta_2singleL.py
Adapted by: Divyae K. Prasad (from Lee Bergstrand's script)

Description:
    - Converts multiline FASTA sequences to single-line format.
    - Converts all sequence characters to uppercase.
    - Removes duplicate sequences or sequence names using seqkit.

Usage:
    python multiL_fasta_2singleL.py <input_filename.fasta>

Example:
    python multiL_fasta_2singleL.py mySeqs.fasta

Output:
    - <input_filename>.fasta_SL (Single-line FASTA format)
"""

#===========================================================================================================
# Imports
import sys
import os
import re
from Bio import SeqIO

#===========================================================================================================

def check_args(expected_args):
    """
    Validates the number of command-line arguments.
    Exits with an error message if incorrect.

    :param expected_args: Expected number of arguments (including script name)
    """
    if len(sys.argv) != expected_args:
        print("Error: Incorrect number of arguments.")
        print("Usage: python " + sys.argv[0] + " <sequences.fasta>")
        sys.exit(1)  # Exit with error

#===========================================================================================================

def convert_fasta_to_single_line(input_file, output_file):
    """
    Reads a multiline FASTA file, converts sequences to uppercase,
    and writes them in a single-line FASTA format.

    :param input_file: Path to the input FASTA file
    :param output_file: Path to the output FASTA file
    """
    #print(f">> Reading and processing '{input_file}'...")

    try:
        with open(output_file, "w") as out_f:
            for record in SeqIO.parse(input_file, "fasta"):
                header = f">{record.id}\n"
                sequence = str(record.seq).upper() + "\n"  # Convert to uppercase & ensure single line
                out_f.write(header + sequence)

        # print(f">> Successfully converted {output_file} to single-line format.")

    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

#===========================================================================================================
# Main Execution
if __name__ == "__main__":
    check_args(2)  # Ensure correct number of arguments

    input_file = sys.argv[1]  # Input FASTA file
    output_file = f"{input_file}_SL"  # Output filename with "SL_" prefix

    convert_fasta_to_single_line(input_file, output_file)

    #print(">> Done!")
