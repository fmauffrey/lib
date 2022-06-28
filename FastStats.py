#!/usr/bin/python3

import sys
import statistics

try:
    from Bio import SeqIO
except ModuleNotFoundError:
    sys.exit("\nBiopython not installed\n")


def check_arg():
    if len(sys.argv) != 2:
        return False
    else:
        return True


def check_extension():
    extension = sys.argv[1].split(".")[-1]
    if extension in ["fa", "fna", "fasta"]:
        return "fasta"
    elif extension in ["fq", "fastq"]:
        return "fastq"
    else:
        sys.exit("\nIncorrect file format. File should end in .fq, .fastq, .fa, .fna or .fasta.\n")


def stats(file_ext):
    fast_sequences = SeqIO.parse(sys.argv[1], file_ext)
    seq_num = 0
    seq_length = []

    for seq in fast_sequences:
        seq_num += 1
        seq_length.append(len(seq.seq))

    total_length = sum(seq_length)
    max_length = max(seq_length)
    min_length = min(seq_length)
    mean_length = int(round(statistics.mean(seq_length), 0))
    median_length = int(statistics.median(seq_length))

    print(f"\nNumber of sequences = {seq_num}"
          f"\nTotal sequences length = {total_length} bp"
          f"\nMaximum sequence length = {max_length} bp"
          f"\nMinimum sequence length = {min_length} bp"
          f"\nAverage sequence length = {mean_length} bp"
          f"\nMedian sequence length = {median_length} bp")


if __name__ == "__main__":
    if check_arg():
        seq_type = check_extension()
        stats(seq_type)
    else:
        sys.exit("\n#####################################\n"
                 "########## FastStats.py ############\n"
                 "#####################################\n"
                 "\nExtract statistics from fasta/fastq file\n\n"
                 "Usage: FastStats.py [file]\n\n")
