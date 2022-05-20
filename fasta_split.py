#!/usr/bin/python3

import sys

try:
	from Bio import SeqIO
except ModuleNotFoundError:
	print("\nBiopython not installed\n")
	exit()

""" This script allows to split fasta/fastq files in smaller files with a maximum number of contigs per file or a
maximum length of all contigs per file """

parameters = ["length_max", "contigs_max"]


def check_arg():
	if len(sys.argv) != 4:
		print("\nWrong number of arguments\n")
		return False
	elif sys.argv[1] not in parameters:
		print(f"\n{sys.argv[1]} parameter not valid\n")
		return False
	elif ".fa" not in sys.argv[3] and ".fasta" not in sys.argv[3]:
		print("\nFile name does not end in .fa or .fasta\n")
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


def filter():

	param = sys.argv[1]
	value = sys.argv[2]
	fasta_sequences = SeqIO.parse(sys.argv[3],'fasta')
	x = 1  # Iteration count for naming new files

	if param == "contigs_max":
		seq_groups = []  # store the sequences to be saved together as lists
		seq_list = []
		for fasta in fasta_sequences:
			if len(seq_list) < int(value):
				seq_list.append(fasta)
			else:
				seq_groups.append(seq_list)
				seq_list = []
		if seq_list != []:
			seq_groups.append(seq_list)

		for x in range(len(seq_groups)):
			with open(f"{sys.argv[3]}_{x}.fasta", "w") as output:
				SeqIO.write(seq_groups[x], output, "fasta")

	elif param == "length_max":
		seq_groups = []  # store the sequences to be saved together as lists
		seq_list = []
		length = 0
		for fasta in fasta_sequences:
			if len(fasta.seq) > int(value):
				print(f"\n{fasta.description} is longer than indicated value\n"
				"Exiting ...\n")
				exit()
			elif length+len(fasta.seq) < int(value):
				seq_list.append(fasta)
				length += len(fasta.seq)
			else:
				seq_groups.append(seq_list)
				seq_list = [fasta]
				length = len(fasta.seq)
		seq_groups.append(seq_list)

		for x in range(len(seq_groups)):
			with open(f"{sys.argv[3]}_{x}.fasta", "w") as output:
				SeqIO.write(seq_groups[x], output, "fasta")


if __name__ == "__main__":
	if check_arg():
		filter()
	else:
		print("\n#####################################\n"
		"########## fasta_split.py ############\n"
		"#####################################\n"
		"\nSplit fasta file into multiple files based on chosen parameter\n\n"
		"Usage: fasta_split.py [parameter] [value] [file.fasta]\n\n"
		"Parameters are:\n"
		"length_max = maximum cumulated length of contigs per file\n"
		"contigs_max = maximum number of contigs per file\n\n"
		"File must end in .fa or .fasta\n")
