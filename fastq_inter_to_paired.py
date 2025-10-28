import argparse

"""
Convert interleaved FASTQ to paired FASTQ files.
"""

def convert(input_file, output_prefix):
    with open(input_file, 'r') as infile, \
         open(f"{output_prefix}_1.fastq", 'w') as outfile1, \
         open(f"{output_prefix}_2.fastq", 'w') as outfile2:
        while True:
            # Read 4 lines for read 1
            read1 = [infile.readline() for _ in range(4)]
            if not read1[0]:
                break  # End of file

            # Read 4 lines for read 2
            read2 = [infile.readline() for _ in range(4)]
            if not read2[0]:
                raise ValueError("Input FASTQ file has an odd number of reads.")

            # Write to respective output files
            outfile1.writelines(read1)
            outfile2.writelines(read2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert interleaved FASTQ to paired FASTQ files.")
    parser.add_argument("-i", type=str, dest="input", help="Input interleaved FASTQ file", required=True)
    parser.add_argument("-o", type=str, dest="output", help="Output FASTQ file prefix", required=True)
    args = parser.parse_args()

    convert(args.input, args.output)