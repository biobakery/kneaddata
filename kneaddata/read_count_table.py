#!/usr/bin/env python

import argparse
import sys
import os

INITIAL_DESC="Initial number of reads"
TRIMMED_DESC="Total reads after trimming"
FILTERED_DESC="Total reads after merging results from multiple databases"

FILE_EXTENSIONS={
    "raw pair1":"_R1_001.fastq",
    "raw pair2":"_R2_001.fastq",
    "trimmed pair1":"trimmed.1.fastq",
    "trimmed pair2":"trimmed.2.fastq",
    "trimmed orphan1":"trimmed.single.1.fastq",
    "trimmed orphan2":"trimmed.single.2.fastq",
    "decontaminated pair1":"paired_1.fastq",
    "decontaminated pair2":"paired_2.fastq",
    "decontaminated orphan1":"unmatched_1.fastq",
    "decontaminated orphan2":"unmatched_2.fastq"}

TABLE_COLUMNS=[
    "raw pair1",
    "raw pair2",
    "trimmed pair1",
    "trimmed pair2",
    "trimmed orphan1",
    "trimmed orphan2",
    "decontaminated pair1",
    "decontaminated pair2",
    "decontaminated orphan1",
    "decontaminated orphan2"]

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """

    parser = argparse.ArgumentParser(description='Create a table of read counts for all samples')
    parser.add_argument('logs',nargs='+',help='kneaddata log files')
    parser.add_argument('--output',required=True,help='the output file to write')

    return parser.parse_args()

def get_count(line):
    """
    Return the count number from the log line
    """

    return line.rstrip().split(":")[-1].strip()

def get_file_type(line):
    """
    Return the file name from the log line
    """

    file_name=line.split("(")[-1].split(")")[0].strip()

    file_type=None
    for type, ext in FILE_EXTENSIONS.iteritems():
        if file_name.endswith(ext):
            file_type=type
            break

    return file_type

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
            if INITIAL_DESC in line:
                reads[sample][get_file_type(line)]=get_count(line)
            elif TRIMMED_DESC in line:
                reads[sample][get_file_type(line)]=get_count(line)
            elif FILTERED_DESC in line:
                reads[sample][get_file_type(line)]=get_count(line)
    return reads

def write_table(output, reads):
    """
    Write the table of read counts for all samples
    """

    with open(output, "w") as file_handle:
        file_handle.write("\t".join(["Sample"]+TABLE_COLUMNS)+"\n")
        for sample in sorted(reads.keys()):
            new_line=[sample]
            for column in TABLE_COLUMNS:
                new_line.append(reads[sample][column])
            file_handle.write("\t".join(new_line)+"\n")

def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)

    # Get the reads counts for all logs for all samples
    reads={}
    for file in args.logs:
        reads=get_reads(file,reads)

    # Write the output table
    write_table(args.output, reads)

if __name__ == "__main__":
    main()
