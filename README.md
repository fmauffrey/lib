# FastManip

FastManip is a collection of scripts allowing diverse manipulation of fasta/fastq files. 
FastManip scripts are based on Biopython which must be installed in the local environment.

### Requirement

Python >= 3.8  
Biopython >= 1.79  

    pip instal biopython

### FastSplit

FastSplit splits a fasta/fastq file into multiple smaller files based on a maximum
number of sequence per file or a maximum cumulated sequences length per file. It can 
also split sequences into multiple sequence at a given length.
This can be useful for online analyses with file size limits (e.g. NCBI Blast). 

    FastSplit.py [parameters] [value] [file]  
 
- file: a fasta/fastq file
- parameter:
  - *length_max* -> maximum cumulated length of sequences per file  
  - *contigs_max* -> maximum number of sequences per file  
  - *length_split* -> split sequences at specific length (output = fasta)
- value: 

### FastFilt

FastFilt filters out sequences from fasta/fastq files based on a filtering parameter.

    FastSplit.py [parameter] [value] [file]  

- file: a fasta/fastq file
- parameter:
  - *length_max* -> maximum length of contigs  
  - *length_min* -> minimum length of contigs  
  - *seq_present* -> contigs containing the sequence (or base)  
  - *seq_absent* -> contigs not containing the sequence (or base)
- value: sequences length (for *length_max* and *length_min*) or a sequence STRING 
(for *seq_present* and *seq_absent*)


### FastStats

FastStats can output statistics of a fasta/fastq file (number of sequences, 
maximum/minimum sequences length, total length and sequences length mean/median.

    FastStats.py [file]

- file: a fasta/fastq file

### FastCheck

FastCheck scans sequences string and returns nucleotides composition or print sequences 
with only A/T/C/G bases.

    FastCheck.py [parameter] [file]  
 
- file: a fasta/fastq file
- parameters:
  - *bases_stat* -> sequence nucleotides composition 
  - *filter_non_atgc* -> prints sequences composed of A/T/G/C only.
