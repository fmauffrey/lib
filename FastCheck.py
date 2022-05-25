#!/usr/bin/python3

import sys

try:
    from Bio import SeqIO
except ModuleNotFoundError:
    sys.exit("\nBiopython not installed\n")

parameters = ["filter_non_atgc", "bases_stat"]


def check_arg():
    if len(sys.argv) != 3:
        print("\nWrong number of arguments\n")
        return False
    elif sys.argv[1] not in parameters:
        print(f"\n{sys.argv[1]} parameter not valid\n")
        return False
    else:
        return True


def check_extension():
    extension = sys.argv[2].split(".")[-1]
    if extension in ["fa", "fna", "fasta"]:
        return "fasta"
    elif extension in ["fq", "fastq"]:
        return "fastq"
    else:
        sys.exit("\nIncorrect file format. File should end in .fq, .fastq, .fa, .fna or .fasta.\n")


def atgc(file_ext):
    param = sys.argv[1]
    fast_sequences = SeqIO.parse(sys.argv[2], file_ext)

    for seq in fast_sequences:
        comp = {"A": 0, "T": 0, "C": 0, "G": 0}
        for char in seq.seq:
            if char in comp:
                comp[char] += 1
            else:
                comp[char] = 1

        gc = round((comp["G"]+comp["C"])/(comp["G"]+comp["C"]+comp["A"]+comp["T"])*100, 2)

        if param == "bases_stat":
            string_comp = " ".join(f"{x}:{comp[x]}" for x in comp)
            print(f"{seq.id}\t{string_comp}\tGC:{gc} %")

        elif param == "filter_non_atgc":
            if len(comp) == 4:
                print(f">{seq.description}\n{seq.seq}")


if __name__ == "__main__":
    if check_arg():
        seq_type = check_extension()
        atgc(seq_type)
    else:
        sys.exit("\n#####################################\n"
                 "############ FastCheck.py ###########\n"
                 "#####################################\n"
                 "\nBases analyses of fasta/fastq file\n\n"
                 "Usage: FastCheck.py [parameter] [file]\n\n"
                 "Parameters are:\n"
                 "filter_non_atgc = print only sequences with ATGC bases\n"
                 "bases_stat = prints bases composition of sequences\n")
