## FastManip

FastManip is a suite of scripts allowing diverse manipulation of fasta and fastq files.  
FastManip scripts are based on Biopython.

Requirement

Python >= 3.8  
Biopython >= 1.79  

### FastSplit

FastSplit can split a fasta/fastq file into multiple smaller files based on a filtering parameter.

FastSplit.py [parameter] [value] [file]  
Parameters:  
*length_max* = maximum cumulated length of contigs per file  
*contigs_max* = maximum number of contigs per file  
*length_split* = split sequences at specific length (output = fasta)

### FastFilt

FastFilt can filter out sequences from fasta/fastq files based on a filtering parameter.

FastSplit.py [parameter] [value] [file]  
Parameters:  
*length_max* = maximum length of contigs  
*length_min* = minimum length of contigs  
*seq_present* = contigs containing the sequence (or base)  
*seq_absent* = contigs not containing the sequence (or base)

### FastStats

FastStats can print different information of a fasta/fastq file.

FastStats.py [file]

### FastCheck

FastCheck can scan sequences string and return bases composition

FastCheck.py [parameter] [file]  
Parameters:  
*bases_stat* = sequence nucleotides composition   
*filter_non_atgc* = output only sequences composed of ATGC