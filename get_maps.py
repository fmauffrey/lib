#!/usr/bin/python3

import sys
import requests
import re
import time
import argparse

""" The script takes a list of KO as input in a file and return a list of Kegg pathways with the number of KO including in each pathway. """


def kegg_request():
    # Send request to Kegg API for each KO to get corresponding maps and gather results for all KOs
    
    # Initialize variables
    KO_list = open(args.input, "r").read().splitlines() # list of KOs
    maps = {} # stores KO count for each map
    KO_total = len(KO_list)
    count = 0

    # Loop over KOs and send request. If failed, retry until success.
    with open(f"{args.logs}/log_{file_name}", "w") as log:
        for KO in KO_list:
            response = requests.get(f"https://rest.kegg.jp/link/pathway/{KO}") # link command of the Kegg API
            status = response.status_code
            content = ""

            # Success
            if status == 200:
                content = response.content.decode("utf-8")  # get content (bites) and convert into string

            # Fail due to too many requests at the same time
            elif status == 403:
                time.sleep(1)
                retry = True
                while retry:
                    print(f"API error ({str(status)}) for {KO}. Retrying ...")
                    log.write(f"API error ({str(status)}) for {KO}\n")
                    retry_response = requests.get(f"https://rest.kegg.jp/link/pathway/{KO}")
                    if retry_response.status_code == 200:
                        content = response.content.decode("utf-8")  # get content (bites) and convert into string
                        retry = False

            # Other errors
            else:
                print(f"API error ({str(status)}) for {KO}")
                log.write(f"API error ({str(status)}) for {KO}\n")
                count += 1
                print(f"{args.input}\t{count}/{KO_total}")
                pass

            # Parse the response content to get the modules names in a list
            modules = re.findall(r"map\d*", content)
            for elt in modules:
                if elt not in maps:
                    maps[elt] = 1
                else:
                    maps[elt] += 1
            
            count += 1
            print(f"{args.input}\t{count}/{KO_total}")
            log.write(f"{args.input}\t{count}/{KO_total}\n")
            
            time.sleep(0.5)

        return maps

def save_maps(maps_dict, name):
    # Save the maps count in a file
    with open(f"{args.output}/{name}", "w") as out:
        for elt in maps_dict:
            out.write(f"{elt}\t{str(maps_dict[elt])}\n")

if __name__ == "__main__":

    # Argparse argument parser
    parser = argparse.ArgumentParser(description="Get KO count per map from KO list")
    parser.add_argument('-i', type=str, metavar="FILE", dest="input", help='list of KO', required=True)
    parser.add_argument('-o', type=str, metavar="FOLDER", dest="output", help='output folder [./]', default="./")
    parser.add_argument('-l', type=str, metavar="FOLDER", dest="logs", help='logs folder [./]', default="./")
    args = parser.parse_args()

    file_name = args.input.split('/')[-1] # get the file name from the input path

    maps_count = kegg_request()
    save_maps(maps_count, file_name)