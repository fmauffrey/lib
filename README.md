Most of the scripts use Biopython for the analysis. It must be installed in your environment for using these scripts:

```
pip instal biopython
```

## FastSplit

FastSplit splits a fasta/fastq file into multiple smaller files based on a maximum
number of sequence per file or a maximum cumulated sequences length per file. It can 
also split sequences into multiple sequence at a given length.
This can be useful for online analyses with file size limits (e.g. NCBI Blast). 

```
FastSplit.py [parameters] [value] [file]
```
 
- file: a fasta/fastq file
- parameter:
  - *length_max* -> maximum cumulated length of sequences per file  
  - *contigs_max* -> maximum number of sequences per file  
  - *length_split* -> split sequences at specific length (output = fasta)
- value: 

## FastFilt

FastFilt filters out sequences from fasta/fastq files based on a filtering parameter.

```
FastFilt.py [parameter] [value] [file]
```

- file: a fasta/fastq file
- parameter:
  - *length_max* -> maximum length of contigs  
  - *length_min* -> minimum length of contigs  
  - *seq_present* -> contigs containing the sequence (or base)  
  - *seq_absent* -> contigs not containing the sequence (or base)
- value: sequences length (for *length_max* and *length_min*) or a sequence STRING 
(for *seq_present* and *seq_absent*)


## FastStats

FastStats can output statistics of a fasta/fastq file (number of sequences, 
maximum/minimum sequences length, total length and sequences length mean/median.

```
FastStats.py [file]
```

- file: a fasta/fastq file

## FastCheck

FastCheck scans sequences string and returns nucleotides composition or print sequences 
with only A/T/C/G bases.

```
FastCheck.py [parameter] [file]
```
 
- file: a fasta/fastq file
- parameters:
  - *bases_stat* -> sequence nucleotides composition 
  - *filter_non_atgc* -> prints sequences composed of A/T/G/C only.

## collapse_and_count

Remove duplicate sequences in a fasta file based on sequence description. The count of each
sequence is kept and added in the header.

```
collapse_and_count.py input output
```
 
- input: input fasta file (nucl or prot)
- output: output fasta file

## get_maps

Takes a list of Kegg KOs as input in a file and return a list of Kegg pathways with the number of KO including in each pathway.

```
get_maps.py -i FILE [-o FOLDER] [-l FOLDER]
```
 
-i: file with KOs
-o: output folder [.\]
-l: log folder [.\]
