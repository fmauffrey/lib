#!/usr/bin/python3

import argparse
import sys
import os
from collections import Counter
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

""" Takes a Genbank file as input and reorder it based on a specified CDS. """

def count_to_comment(count):
    # Convert features count dictionary into text for the comment section
    text = ""
    for value, feature in enumerate(count):
        if len(feature) > 3:
            text += f"{feature}\t\t{value}\n"
        else:
            text += f"{feature}\t\t\t{value}\n"
    
    return text

def check_files(folder, command):
    # Check if folder is empty or only one file is found
    if not os.listdir(folder):
        print("Provided folder is empty")
        sys.exit()
    elif len(os.listdir(folder)) == 1 and command in ["find", "order"]:
        print("Only one file in the folder. Pick any gene as reference")
        sys.exit()

    return(True)

def find_target(folder):
    # Find genes present in all files
    genes = {} # genes as key, count as variable

    for f in os.listdir(folder):
        gbk = SeqIO.parse(f"{folder}/{f}", "genbank")
        record = next(gbk)
        genes_list = []

        # Loop over all entries and count each gene. Remove duplicated genes.
        for entry in record.features[1:]:
            if entry.type == "gene":
                try:
                    gene_name = entry.qualifiers["gene"][0]
                    genes_list.append(gene_name)
                except KeyError:
                    pass

        # Remove duplicates
        nodup_genes_list = []
        [nodup_genes_list.append(x) for x in genes_list if genes_list.count(x) == 1 and x not in nodup_genes_list]

        # Count this gene in the global dictionary
        for gene_name in nodup_genes_list:
            if gene_name in genes:
                genes[gene_name] += 1
            else:
                genes[gene_name] = 1

    # Print all genes present in all samples
    print("The following genes can be used as reference:")
    for gene in genes:
        if genes[gene] == len(os.listdir(folder)):
            print(f"- {gene}")

def reorder(in_folder, out_folder, target):
    for f in os.listdir(in_folder):
        gbk = SeqIO.parse(f"{in_folder}/{f}", "genbank")
        record = next(gbk)
        new_features =[record.features[0]] # The first object is the total sequence which is unchanged
        total_length = record.features[0].location.end
        shift = 0
        new_origin = 0
        origin_strand = 1

        # First loop looking for the target gene and getting its position and its strand
        for entry in record.features[1:]:
            try:
                gene_name = entry.qualifiers["gene"][0]
                if gene_name == target:
                    new_origin = entry.location.start
                    origin_strand = entry.location.strand
                    shift = total_length - entry.location.end
                    break
            except KeyError:
                pass

        # Second loop changing position of features
        for entry in record.features[1:]:

            # If the origin gene is on strand 1, it should at the beginning of the sequence
            if origin_strand == 1:

                # If the gene spans over multiple region, it is certainly because it was cut by the original origin 
                if entry.location_operator == "join":
                    parts = entry.location.parts

                    current_start_location = max(parts[0].start, parts[1].start) # grab the good start position (always the higher number)
                    current_end_location = min(parts[0].end, parts[1].end) # grab the good end position (always the lowest number)
                    strand = parts[0].strand

                    new_start = current_start_location - new_origin
                    new_end = current_end_location - new_origin + total_length
                    new_location = FeatureLocation(new_start, new_end, strand)
                    entry.location = new_location

                else:
                    current_start_location = entry.location.start
                    current_end_location = entry.location.end
                    strand = entry.location.strand

                    # From the new origin, only substrate the target original location of the new origin 
                    if current_start_location >= new_origin:
                        new_start = current_start_location - new_origin
                        new_end = current_end_location - new_origin
                        new_location = FeatureLocation(new_start, new_end, strand)
                        entry.location = new_location
                    
                    # Before the new origin, substrate the target original location and add genome length 
                    else:
                        new_start = current_start_location - new_origin + total_length
                        new_end = current_end_location - new_origin + total_length
                        new_location = FeatureLocation(new_start, new_end, strand)
                        entry.location = new_location
                
                new_features.append(entry)

            # If the origin gene is on strand -1, it should at the end of the sequence
            if origin_strand == -1:

                # If the gene spans over multiple region, it is certainly because it was cut by the original origin 
                if entry.location_operator == "join":
                    parts = entry.location.parts

                    current_start_location = max(parts[0].start, parts[1].start) # grab the good start position (always the higher number)
                    current_end_location = min(parts[0].end, parts[1].end) # grab the good end position (always the lowest number)
                    strand = parts[0].strand

                    new_start = current_start_location + shift - total_length
                    new_end = current_end_location + shift
                    new_location = FeatureLocation(new_start, new_end, strand)
                    entry.location = new_location

                else:
                    current_start_location = entry.location.start
                    current_end_location = entry.location.end
                    strand = entry.location.strand

                    # Before the new origin, add the shift
                    if current_start_location <= new_origin:
                        new_start = current_start_location + shift 
                        new_end = current_end_location + shift
                        new_location = FeatureLocation(new_start, new_end, strand)
                        entry.location = new_location
                    
                    # After the new origin, add the shift and substrate the genome length 
                    else:
                        new_start = current_start_location + shift - total_length
                        new_end = current_end_location + shift - total_length
                        new_location = FeatureLocation(new_start, new_end, strand)
                        entry.location = new_location
                
                new_features.append(entry)

        # Replacing features with new position in the seqRecord object and save it
        record.features = new_features
        if not os.path.exists(out_folder):
            os.mkdir(out_folder)
        with open(f"{out_folder}/reordered_{f}", "w") as out:
            SeqIO.write(record, out, "genbank")

def extract_region(in_folder, out_folder, target, range):
    # Extract a specific region of all gbk around the given target
    for f in os.listdir(in_folder):
        gbk = SeqIO.parse(f"{in_folder}/{f}", "genbank")
        record = next(gbk)
        new_features =[record.features[0]] # The first object is the source which is unchanged
        range_start, range_end = 0, 0
        list_features = [] # list of feature type for the comment section

        # First loop looking for the target gene and getting its position and its strand
        for entry in record.features[1:]:
            try:
                gene_name = entry.qualifiers["gene"][0]
                if gene_name == target:
                    range_start = entry.location.start - int(range)
                    range_end = entry.location.end + int(range)
                    break
            except KeyError:
                pass
        
        # Second loop keeping only features in the position range
        for entry in record.features[1:]:
            if range_start <= entry.location.start <= range_end or range_start <= entry.location.end <= range_end:
                new_features.append(entry)
                list_features.append(entry.type)
        count_features = Counter(list_features)
        
        # Extract the position of the first and last features and adapting global positions
        first_start = new_features[1].location.start
        last_end = new_features[-1].location.end

        # Changing global position
        new_features[0].location = FeatureLocation(0, last_end - first_start)
        
        # Changing position of each feature
        for feature in new_features[1:]:
            index = new_features.index(feature)
            new_features[index].location = FeatureLocation(new_features[index].location.start - first_start,
                                                           new_features[index].location.end - first_start,
                                                           new_features[index].location.strand)
            
        # Adapting record information
        record.features = new_features
        record.description += f", {int(range)*2} bp sequence around {target}"
        record.seq = record.seq[first_start:last_end]
        record.annotations["comment"] = "Genbank file filtered with gbk_order.py\n\n" + count_to_comment(count_features)
        
        # Saving file
        if not os.path.exists(out_folder):
            os.mkdir(out_folder)
        with open(f"{out_folder}/extract_{f}", "w") as out:
            SeqIO.write(record, out, "genbank")

if __name__ == "__main__":
    # General argument parser
    parser = argparse.ArgumentParser(description="Genbank files functions.")

    # Add subparser to allow the different commands
    subparser = parser.add_subparsers(dest='command')
    find_parser = subparser.add_parser('find', description = "Find genes common to all Genbank files. Duplicated genes are ignored.")
    order_parser = subparser.add_parser('order', description = "Re-order Genbank files considering a gene in common as origin.")
    extract_parser = subparser.add_parser('extract', description = "Extract region surrounding a specified gene.")

    # Find parser arguments
    find_parser.add_argument('-i', type=str, metavar="FOLDER", dest="input_folder_find", help='Folder with all Genbank files', required=True)

    # Argparse argument parser
    order_parser.add_argument('-i', type=str, metavar="FOLDER", dest="input_folder_reorder", help='Folder with all Genbank files', required=True)
    order_parser.add_argument('-o', type=str, metavar="FOLDER", dest="output_folder", help='Output folder for all reordered Genbank files', required=True)
    order_parser.add_argument('-t', type=str, metavar="STR", dest="target", help='Target gene considered as origin', required=True)

    # Argparse argument parser
    extract_parser.add_argument('-i', type=str, metavar="FOLDER", dest="input_folder", help='Folder with all Genbank files', required=True)
    extract_parser.add_argument('-o', type=str, metavar="FOLDER", dest="output_folder", help='Output folder for all extracted regions', required=True)
    extract_parser.add_argument('-t', type=str, metavar="STR", dest="target", help='Target gene', required=True)
    extract_parser.add_argument('-r', type=str, metavar="NUM", dest="range", help='Sequence length to extract before and after the target', required=True)
    
    args = parser.parse_args()
    
    if args.command not in ["find", "order", "extract"]:
        sys.exit("Commands available: find, order, extract")

    if args.command == "find":
        if check_files(args.input_folder_find, "find"):
            find_target(args.input_folder_find)

    elif args.command == "order":
        if check_files(args.input_folder_reorder, "order"):
            reorder(args.input_folder_reorder, args.output_folder, args.target)
            print("Sequences reordered.")
    
    elif args.command == "extract":
        if check_files(args.input_folder, "extract"):
            extract_region(args.input_folder, args.output_folder, args.target, args.range)
            print("Regions extracted. For the moment, this script does not handle regions spanning over the end and the start of the sequence for circular sequences.")