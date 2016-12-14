#!/usr/bin/env python

import argparse
import sys
import os

READ_COUNT_IDENTIFIER="READ COUNT"
OUTPUT_ORDER=["raw","trimmed","decontaminated","final"]

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """

    parser = argparse.ArgumentParser(description='Create a table of read counts for all samples')
    parser.add_argument('--input',required=True,help='the input folder with kneaddata log files')
    parser.add_argument('--output',required=True,help='the output file to write')

    return parser.parse_args()

def get_read_count_type(line):
    """
    Return the count number and read type from the log line
    """

    data = line.rstrip().split(":")
    
    count = data[-1].strip()
    type = data[-3].lstrip().rstrip()
    
    return count, type

def get_reads(file, reads=None):
    """
    Get the read counts from the file 
    """

    if not reads:
        reads={}

    # get the sample name from the file
    sample=os.path.basename(file).split(".")[0]
    
    reads[sample]={}
    with open(file) as file_handle:
        for line in file_handle:
            if READ_COUNT_IDENTIFIER in line:
                count, type = get_read_count_type(line)
                reads[sample][type]=count

    return reads

def write_table(output, reads):
    """
    Write the table of read counts for all samples
    """

    with open(output, "w") as file_handle:
        # get the headers from all of the samples
        headers=set()
        for sample, counts in reads.items():
            headers.update(counts.keys())
        # order the headers
        header_order=[]
        for column in OUTPUT_ORDER:
            header_order+=sorted(list(filter(lambda x: x.startswith(column),headers)))
        
        file_handle.write("\t".join(["Sample"]+header_order)+"\n")
        for sample in sorted(reads.keys()):
            new_line=[sample]
            for column in header_order:
                try:
                    counts=reads[sample][column]
                except KeyError:
                    counts="NA"
                new_line.append(counts)
            file_handle.write("\t".join(new_line)+"\n")

def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)

    # Find the log files
    files=[os.path.join(args.input,file) for file in os.listdir(args.input)]
    logs=filter(lambda file: file.endswith(".log") and os.path.isfile(file), files)

    # Get the reads counts for all logs for all samples
    reads={}
    for file in logs:
        print("Reading log: " + file)
        reads=get_reads(file,reads)

    # Write the output table
    write_table(args.output, reads)
    print("Read count table written: " + args.output)

if __name__ == "__main__":
    main()
