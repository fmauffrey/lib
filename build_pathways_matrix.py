#!/usr/bin/python3

import argparse
import pandas as pd
import os
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

""" The script takes as input pathway counts file generated with get_maps.py and create a count matrix """

def create_dataframe(folder):
    # Create the dataframe with all pathways as column and all genomes as rows

    # Get files names
    files = os.listdir(folder)

    # Create dataframe
    df = pd.DataFrame()

    # Iterate over files and add the pathway column if not present in the dataframe
    for map_file in files:
        text = open(f"{folder}/{map_file}", "r").readlines()
        for line in text:
            pathway = line.split("\t")[0] # get pathway name
            if pathway not in list(df.columns):
                df[pathway] = [0]*len(files) # columns of 0
    
    # Change rows index
    df.index = files

    return df

def complete_df(dataframe, folder):
    # Iterate over files and add the values in the dataframe
    
    # Get files names
    files = os.listdir(folder)

    # Iterate over files and add the pathway column if not present in the dataframe
    for map_file in files:
        text = open(f"{folder}/{map_file}", "r").readlines()
        for line in text:
            pathway = line.split("\t")[0] # get pathway name
            count = line.split("\t")[1].replace("\n", "") # get pathway count
            dataframe.at[map_file, pathway] = count # update the corresponding cell value
    
    return dataframe

if __name__ == "__main__":

    # Argparse argument parser
    parser = argparse.ArgumentParser(description="Compute count matrix from pathway counts files generated with get_maps.py")
    parser.add_argument('-f', type=str, metavar="FOLDER", dest="input", help='folder containing pathways count files', required=True)
    parser.add_argument('-o', type=str, metavar="FILE", dest="output", help='output matrix name [KO_matrix.csv]', default="KO_matrix.csv")
    args = parser.parse_args()

    df = create_dataframe(args.input)
    completed_df = complete_df(df, args.input)
    completed_df.to_csv(args.output) # save dataframe