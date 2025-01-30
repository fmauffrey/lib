#!/usr/bin/python3

import argparse
import sys
import statistics
from collections import Counter
try:
    from Bio import SeqIO
except ModuleNotFoundError:
    sys.exit("\nBiopython required\n")

def validate_file(fasta_file):
    # Check if file is in fasta format
    validation_fasta = SeqIO.parse(fasta_file, "fasta")
    return (any(validation_fasta))

def fasta_bases(input_fasta):
    # Bases statistics
    bases = {}
    for seq in SeqIO.parse(input_fasta, "fasta"):
        count = Counter(seq.seq)
        for base in count.keys():
            bases.setdefault(base, 0)
            bases[base] += count[base]

    print("\n".join(f"{base}:{count}" for base, count in bases.items()))

def fasta_filter(input_fasta, min_length, max_length, seq_present, seq_absent):
    # Filter contigs based on sequence length and motif presence/absence 
    remaining, discarded = 0, 0

    for seq in SeqIO.parse(input_fasta, "fasta"):
        if len(seq.seq) > min_length and len(seq.seq) < max_length and seq_present in seq.seq and seq_absent not in seq.seq:
            print(f">{seq.description}\n{seq.seq}")
            remaining += 1
        else:
            discarded += 1
    
    print(f"\033[94mRemaining sequences: {remaining}\nDiscarded sequences: {discarded}\033[0m", file = sys.stderr)

def fasta_split(input_fasta, parameter, value, output_folder):
    # Split fasta file based on number of contig or maximum length
    if parameter == "max_length":
        iteration = 1
        sequences = []
        cumulated_length = 0
        long_contigs = []
        for seq in SeqIO.parse(input_fasta, "fasta"):
            # If contigs larger than max size
            if len(seq.seq) > value:
                long_contigs.append(seq.id)
                with open(f"{output_folder}/{iteration}-{input_fasta}", "w") as out:
                    SeqIO.write([seq], out, "fasta")
                iteration += 1

            else:
                if cumulated_length + len(seq.seq) < value:
                    sequences.append(seq)
                    cumulated_length += len(seq.seq)
                else:
                    with open(f"{output_folder}/{iteration}-{input_fasta}", "w") as out:
                        SeqIO.write(sequences, out, "fasta")
                    iteration += 1
                    sequences = [seq]
                    cumulated_length = len(seq.seq)

        print(f"\033[94m{len(long_contigs)} contigs had a length > {value} bp\033[0m", file = sys.stderr)
    
    elif parameter == "max_contigs":
        iteration = 1
        sequences = []
        for seq in SeqIO.parse(input_fasta, "fasta"):
            if len(sequences) < value:
                sequences.append(seq)
            else:
                with open(f"{output_folder}/{iteration}-{input_fasta}", "w") as out:
                    SeqIO.write(sequences, out, "fasta")
                iteration += 1
                sequences = [seq]
        
        with open(f"{output_folder}/{iteration}-{input_fasta}", "w") as out:
            SeqIO.write(sequences, out, "fasta")

def fasta_stats(input_fasta):
    # Get fasta basic length stats
    seq_num = 0
    seq_length = []

    for seq in SeqIO.parse(input_fasta, "fasta"):
        seq_num += 1
        seq_length.append(len(seq.seq))

    total_length = sum(seq_length)
    max_length = max(seq_length)
    min_length = min(seq_length)
    mean_length = int(round(statistics.mean(seq_length), 0))
    median_length = int(statistics.median(seq_length))

    print(f"Number of sequences = {seq_num}"
          f"\nTotal sequences length = {total_length} bp"
          f"\nMaximum sequence length = {max_length} bp"
          f"\nMinimum sequence length = {min_length} bp"
          f"\nAverage sequence length = {mean_length} bp"
          f"\nMedian sequence length = {median_length} bp")

def fasta_collapse(input_fasta, output_fasta):
    # Remove sequences with identical header (keep the first one)
    counts = {}

    for seq in SeqIO.parse(input_fasta, "fasta"):
        if seq.description not in counts:
            counts[seq.description] = [seq, 1]
        else:
            counts[seq.description][1] += 1

    final_seq = []

    for elt in counts:
        new_seq = counts[elt][0]
        new_seq.description = f"{new_seq.description}_count_{counts[elt][1]}"
        final_seq.append(new_seq)

    with open(output_fasta, "w") as out:
        SeqIO.write(final_seq, out, "fasta")

if __name__ == "__main__":
    # General argument parser
    general_parser = argparse.ArgumentParser(description="Fasta files manipulation and stats.")

    # Add subparser to allow the different commands
    subparser = general_parser.add_subparsers(dest='command')
    bases_parser = subparser.add_parser('bases', description = "Get the occurence of each nucleotide in the assembly")
    filter_parser = subparser.add_parser('filter', description = "Filter contigs based on sequence length and motif presence/absence ")
    split_parser = subparser.add_parser('split', description = "Split fasta file based on number of contigs or maximum length")
    stats_parser = subparser.add_parser('stats', description = "Get basic length stats")
    collapse_parser = subparser.add_parser('collapse', description = "Remove sequences with identical header (keeping the first occurence)")

    # Check parser arguments
    bases_parser.add_argument('-i', type=str, metavar="FASTA", dest="input_fasta", help='Input fasta file', required=True)

    # Filter parser arguments
    filter_parser.add_argument('-i', type=str, metavar="FASTA", dest="input_fasta", help='Input fasta file', required=True)
    filter_parser.add_argument('--min', type=int, metavar="LENGTH", dest="min_length", help='Minimum contig length', default = 0)
    filter_parser.add_argument('--max', type=int, metavar="LENGTH", dest="max_length", help='Maximum contig length', default = 99999999999)
    filter_parser.add_argument('--present', type=str, metavar="SEQUENCE", dest="seq_present", help='Keep contigs containing this sequence', default = "")
    filter_parser.add_argument('--absent', type=str, metavar="SEQUENCE", dest="seq_absent", help='Keep contigs not containing this sequence', default = "ZZZZZZZZ")

    # Split argument parser
    split_parser.add_argument('-i', type=str, metavar="FASTA", dest="input_fasta", help='Input fasta file', required=True)
    split_parser.add_argument('-p', type=str, metavar="PARAMETER", dest="split_parameter", help="One of max_length (maximum cumulated length of contigs per file) or "
                 "max_contigs (maximum number of contigs per file) or length_split", required=True)
    split_parser.add_argument('-n', type=int, metavar="VALUE", dest="split_value", help='Value to use for splitting', required=True)
    split_parser.add_argument('-o', type=str, metavar="FOLDER", dest="output_folder", help='Output folder [\.]', default="./")

    # Stats argument parser
    stats_parser.add_argument('-i', type=str, metavar="FASTA", dest="input_fasta", help='Input fasta file', required=True)

    # Collapse argument parser
    collapse_parser.add_argument('-i', type=str, metavar="FASTA", dest="input_fasta", help='Input fasta file', required=True)
    collapse_parser.add_argument('-o', type=str, metavar="FASTA", dest="output_fasta", help='Output fasta file', required=True)
    
    args = general_parser.parse_args()
    
    if args.command not in ["bases", "filter", "split", "stats", "collapse"]:
        sys.exit("Commands available: bases, filter, split, stats, collapse")

    if not validate_file(args.input_fasta):
        sys.exit(f"{args.input_fasta} is not a fasta file")

    if args.command == "bases":
        fasta_bases(args.input_fasta)
    elif args.command == "filter":
        fasta_filter(args.input_fasta, int(args.min_length), int(args.max_length), args.seq_present, args.seq_absent)
    elif args.command == "split":
        if args.split_parameter in ["max_length", "max_contigs"]:
            fasta_split(args.input_fasta, args.split_parameter, int(args.split_value), args.output_folder)
        else:
            sys.exit("p = [max_length] or [max_contigs]")
    elif args.command == "stats":
        fasta_stats(args.input_fasta)
    elif args.command == "collapse":
        fasta_collapse(args.input_fasta, args.output_fasta)