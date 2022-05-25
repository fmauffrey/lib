#!/usr/bin/python3

import sys
import math

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ModuleNotFoundError:
    sys.exit("\nBiopython not installed\n")

parameters = ["length_max", "contigs_max", "length_split"]


def check_arg():
    if len(sys.argv) != 4:
        print("\nWrong number of arguments\n")
        return False
    elif sys.argv[1] not in parameters:
        print(f"\n{sys.argv[1]} parameter not valid\n")
        return False
    else:
        try:
            int(sys.argv[2])
            return True
        except ValueError:
            print("\nParameter value must be an integer\n")
            return False


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
    value = int(sys.argv[2])
    fast_sequences = SeqIO.parse(sys.argv[3], file_ext)

    if param == "contigs_max":
        seq_groups = []  # store the sequences to be saved together as lists
        seq_list = []
        for fasta in fast_sequences:
            if len(seq_list) < value:
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
            if len(fasta.seq) > value:
                sys.exit(f"\n{fasta.description} is longer than indicated value\n"
                         "Exiting ...\n")
            elif length + len(fasta.seq) < value:
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

    elif param == "length_split":
        seq_list = []
        for seq in fast_sequences:
            if len(seq.seq) < value:
                seq_list.append(seq)
            else:
                nb = int(math.trunc(len(seq.seq) / value) + 1)
                for x in range(1, nb + 1):
                    record = SeqRecord(seq.seq[(x * value) - value:(x * value)],
                                       id=seq.id,
                                       name=seq.name,
                                       description=seq.description + "-" + str(x))
                    if len(record.seq) != 0:
                        seq_list.append(record)

        SeqIO.write(seq_list, f"{sys.argv[3]}_split{value}.fasta", "fasta")


if __name__ == "__main__":
    if check_arg():
        seq_type = check_extension()
        seq_split(seq_type)
    else:
        sys.exit("\n#####################################\n"
                 "########## FastSplit.py #############\n"
                 "#####################################\n"
                 "\nSplit fasta/fastq file into multiple files based on chosen parameter\n\n"
                 "Usage: FastSplit.py [parameter] [value] [file]\n\n"
                 "Parameters are:\n"
                 "length_max = maximum cumulated length of contigs per file\n"
                 "contigs_max = maximum number of contigs per file\n"
                 "length_split = split sequences at specific length (output = fasta)")
