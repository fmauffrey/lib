#!/usr/bin/python3

import sys

try:
    from Bio import SeqIO
except ModuleNotFoundError:
    sys.exit("\nBiopython not installed\n")

parameters = ["length_max", "contigs_max"]


def check_arg():
    if len(sys.argv) != 4:
        print("\nWrong number of arguments\n")
        return False
    elif sys.argv[1] not in parameters:
        print(f"\n{sys.argv[1]} parameter not valid\n")
        return False
    else:
        if sys.argv[1] in ["length_max", "contigs_max"]:
            try:
                int(sys.argv[2])
                return True
            except ValueError:
                print("\nlength_max, contigs_max: Value must be an integer\n")
                return False
        else:
            return True


def check_extension():
    extension = sys.argv[3].split(".")[-1]
    if extension in ["fa", "fna", "fasta"]:
        return "fasta"
    elif extension in ["fq", "fastq"]:
        return "fastq"
    else:
        sys.exit("\nIncorrect file format. File should end in .fq, .fastq, .fa, .fna or .fasta.\n")


def seq_split(file_ext):
    param = sys.argv[1]
    value = sys.argv[2]
    fast_sequences = SeqIO.parse(sys.argv[3], file_ext)

    if param == "contigs_max":
        seq_groups = []  # store the sequences to be saved together as lists
        seq_list = []
        for fasta in fast_sequences:
            if len(seq_list) < int(value):
                seq_list.append(fasta)
            else:
                seq_groups.append(seq_list)
                seq_list = [fasta]
        if seq_list:
            seq_groups.append(seq_list)

        for x in range(len(seq_groups)):
            with open(f"{sys.argv[3]}_{x}.fasta", "w") as output:
                SeqIO.write(seq_groups[x], output, file_ext)

    elif param == "length_max":
        seq_groups = []  # store the sequences to be saved together as lists
        seq_list = []
        length = 0
        for fasta in fast_sequences:
            if len(fasta.seq) > int(value):
                sys.exit(f"\n{fasta.description} is longer than indicated value\n"
                         "Exiting ...\n")
            elif length + len(fasta.seq) < int(value):
                seq_list.append(fasta)
                length += len(fasta.seq)
            else:
                seq_groups.append(seq_list)
                seq_list = [fasta]
                length = len(fasta.seq)
        seq_groups.append(seq_list)

        for x in range(len(seq_groups)):
            with open(f"{sys.argv[3]}_{x}.fasta", "w") as output:
                SeqIO.write(seq_groups[x], output, file_ext)


if __name__ == "__main__":
    if check_arg():
        seq_type = check_extension()
        seq_split(seq_type)
    else:
        sys.exit("\n#####################################\n"
                 "########## FastSplit.py ############\n"
                 "#####################################\n"
                 "\nSplit fasta/fastq file into multiple files based on chosen parameter\n\n"
                 "Usage: FastSplit.py [parameter] [value] [file]\n\n"
                 "Parameters are:\n"
                 "length_max = maximum cumulated length of contigs per file\n"
                 "contigs_max = maximum number of contigs per file\n\n")
