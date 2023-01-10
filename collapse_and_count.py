import sys
import os
from Bio import SeqIO

""" This script removes duplicated sequences (same header) in a fasta file and
outputs the uniques sequences with the count in the header."""
""" Use: ./collapse_and_count.py input.fasta output.fasta """

input = sys.argv[1]
output = sys.argv[2]

counts = {}

for seq in SeqIO.parse(input, "fasta"):
    if seq.description not in counts:
        counts[seq.description] = [seq, 1]
    else:
        counts[seq.description][1] += 1

final_seq = []

for elt in counts:
    new_seq = counts[elt][0]
    new_seq.description = f"{new_seq.description}_count_{counts[elt][1]}"
    final_seq.append(new_seq)

with open(output, "w") as out:
    SeqIO.write(final_seq, out, "fasta")