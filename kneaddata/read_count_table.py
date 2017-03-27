#!/usr/bin/env python

import argparse
import sys
import os

try:
    from kneaddata import utilities
except ImportError:
    sys.exit("Please install kneaddata")

READ_COUNT_IDENTIFIER="READ COUNT"

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
    utilities.write_read_count_table(args.output, reads)
    print("Read count table written: " + args.output)

if __name__ == "__main__":
    main()
