#!/usr/bin/python3

import sys

try:
	from Bio import SeqIO
except ModuleNotFoundError:
	print("\nBiopython not installed\n")
	exit()

""" This script allows to filter contigs in fasta/fastq files """

parameters = ["length_max", "length_min", "seq_present", "seq_absent"]

def check_arg():
	if len(sys.argv) != 4:
		print("\nWrong number of arguments\n")
		return(False)
	elif sys.argv[1] not in parameters:
		print(f"\n{sys.argv[1]} parameter not valid\n")
		return(False)
	elif ".fa" not in sys.argv[3] and ".fasta" not in sys.argv[3]:
		print("\nFile name does not end in .fa or .fasta\n")
		return(False)
	else:
		if sys.argv[1] in ["length_max", "length_min"]:
			try:
				int(sys.argv[2])
				return(True)
			except ValueError:
				print("\nlength_max, length_min: Value must be an integer\n")
				return(False)
		else:
			return(True)

def filter():

	param = sys.argv[1]
	value = sys.argv[2]
	fasta_sequences = SeqIO.parse(sys.argv[3],'fasta')

	if param == "length_max":
		for fasta in fasta_sequences:
			if len(fasta.seq) <= int(value):
				print(f">{fasta.description}\n{fasta.seq}")

	elif param == "length_min":
		for fasta in fasta_sequences:
			if len(fasta.seq) >= int(value):
				print(f">{fasta.description}\n{fasta.seq}")

	elif param == "seq_present":
		for fasta in fasta_sequences:
			if value in fasta.seq:
				print(f">{fasta.description}\n{fasta.seq}")

	elif param == "seq_absent":
		for fasta in fasta_sequences:
			if value not in fasta.seq:
				print(f">{fasta.description}\n{fasta.seq}")


if __name__ == "__main__":
	if check_arg():
		filter()
	else:
		print("\n#####################################\n"
		"########## fasta_filt.py ############\n"
		"#####################################\n"
		"\nFilter fasta file based on chosen parameter\n\n"
		"Usage: fasta_filt.py [parameter] [value] [file.fasta]\n\n"
		"Parameters are:\n"
		"length_max = maximum length of contigs\n"
		"length_min = minimum length of contigs\n"
		"seq_present = contigs containing the sequence (or base)\n"
		"seq_absent = contigs not containing the sequence (or base)\n\n"
		"File must end in .fa or .fasta\n")
