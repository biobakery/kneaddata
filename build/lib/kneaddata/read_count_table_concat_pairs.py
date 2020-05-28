#!/usr/bin/env python

""" This script will count the paired/orphan reads from a folder of
    kneaddata output files. It will count raw, trimmed, contaminated (for all dbs), and final. 
    This script should be run on kneaddata output files where the input were paired end
    reads that were concatenated into a single file. """

import sys
import os

try:
    import argparse
except ImportError:
    sys.exit("Please upgrade to at least python v2.7")
    
try:
    from kneaddata import utilities
except ImportError:
    sys.exit("Please install kneaddata")

def read_file_n_lines(file,n):
    """ Read a file n lines at a time """
    
    line_set=[]
    with open(file) as file_handle:
        for line in file_handle:
            if len(line_set) == n:
                yield line_set
                line_set=[]
            line_set.append(line)
    
    # yield the last set
    if len(line_set) == n:
        yield line_set

def count_paired_orphan(file):
    """ Count the paired and orphan reads based on the read id """
    
    sequence_pair1=set()
    sequence_pair2=set()
    orphan=set()
    print("Reading file: " + file)

    # get the first line of the file to determine the read pair delimiter
    with open(file) as file_handle:
        line=file_handle.readline()
        if line.rstrip()[-2] == "/":
            delimiter="/"
        else:
            delimiter=" "
    
    for lines in read_file_n_lines(file, 4):
        # split the id to remove the paired identification
        # allow for both space and forward slash to differentiate read pairs
        if delimiter == " ":
            data = lines[0].rstrip().replace("@","",1).split(" ")
            try:
                pair_set = data[1].split(":")[0]
            except IndexError:
                pair_set = "NA"
        else:
            data = lines[0].rstrip().replace("@","",1).split("/")
            try:
                pair_set = data[1]
            except IndexError:
                pair_set = "NA"
   
        if pair_set == "1":
            sequence_pair1.add(data[0])
        elif pair_set == "2":
            sequence_pair2.add(data[0])
        else:
            orphan.add(data[0])
        
    # count pairs as those with two counts an orphan otherwise
    pair1=len(list(sequence_pair1.intersection(sequence_pair2)))
    pair2=pair1
    orphan1=len(list(sequence_pair1.difference(sequence_pair2)))+len(list(orphan))
    orphan2=len(list(sequence_pair2.difference(sequence_pair1)))
            
    print("Found counts of pair1, pair2, orphan1, orphan2: " + ",".join(str(i) for i in [pair1, pair2, orphan1, orphan2]))

    return pair1, pair2, orphan1, orphan2

def get_sample_name(file, sample_names):
    """ Figure out the sample name for the file based on the file name and sample names.
    This allows for multiple delimiters in the sample name. """
    
    # get the basename of the file
    file=os.path.basename(file)

    # get the names of all of the samples, sorted by length with the longest first
    for name in sorted(sample_names,key=len,reverse=True):
        if file.startswith(name):
            return name

    return "Unknown"

def parse_arguments(args):
    """ Parse the arguments from the user"""
    
    parser=argparse.ArgumentParser(
        description="Count paired/orphan reads from kneaddata run on a set of paired end reads concatenated into a single file.")
    parser.add_argument("--input", help="The folder of kneaddata output files for all samples", required=True)
    parser.add_argument("--output", help="The output file of reads to write.", required=True)
    
    return parser.parse_args()

def main():
    # parse arguments
    args = parse_arguments(sys.argv)
    
    # get all of the kneaddata fastq and log files from the folder
    files=[os.path.join(args.input, file) for file in os.listdir(args.input)]
    logs=filter(lambda x: x.endswith(".log"),files)
    trimmed_fastq=filter(lambda x: x.endswith(".fastq") and not "contam" in x and "trimm" in x,files)
    contam_fastq=filter(lambda x: x.endswith(".fastq") and "contam" in x and not "trimm" in x,files)
    final_fastq=filter(lambda x: x.endswith(".fastq") and not "contam" in x and not "trimm" in x,files)

    # get the raw single read counts from the log files, this is 2x the number of paired end reads
    reads={}
    for file in logs:
        for line in open(file):
            if "READ COUNT: raw single :" in line:
                read_count=line.rstrip().split()[-1]
                sample_name=os.path.basename(line.rstrip().split()[-3]).split(".")[0]
                if not sample_name in reads:
                    reads[sample_name]={}
                reads[sample_name]["raw pair1"]=int(read_count)/2
                reads[sample_name]["raw pair2"]=int(read_count)/2
    
    # get the counts for the trimmed files
    for file in trimmed_fastq:
        sample_name = get_sample_name(file, reads.keys())
        pair1, pair2, orphan1, orphan2 = count_paired_orphan(file)
        reads[sample_name]["trimmed pair1"]=pair1
        reads[sample_name]["trimmed pair2"]=pair2
        reads[sample_name]["trimmed orphan1"]=orphan1
        reads[sample_name]["trimmed orphan2"]=orphan2
    
    # get the counts for the contaminate files
    for file in contam_fastq:
        sample_name = get_sample_name(file, reads.keys())
        database_name = os.path.basename(file).replace(sample_name,"",1).replace("_bowtie2_contam.fastq","")
        pair1, pair2, orphan1, orphan2 = count_paired_orphan(file)
        tag = "decontaminated "+database_name
        reads[sample_name][tag+" pair1"]=reads[sample_name]["trimmed pair1"]-pair1
        reads[sample_name][tag+" pair2"]=reads[sample_name]["trimmed pair2"]-pair2
        reads[sample_name][tag+" orphan1"]=reads[sample_name]["trimmed orphan1"]-orphan1
        reads[sample_name][tag+" orphan2"]=reads[sample_name]["trimmed orphan2"]-orphan2   
        
    # get the counts for the final files
    for file in final_fastq:
        sample_name = get_sample_name(file, reads.keys())
        pair1, pair2, orphan1, orphan2 = count_paired_orphan(file)
        reads[sample_name]["final pair1"]=pair1
        reads[sample_name]["final pair2"]=pair2
        reads[sample_name]["final orphan1"]=orphan1
        reads[sample_name]["final orphan2"]=orphan2        
    
    # write the output file
    utilities.write_read_count_table(args.output,reads)
    
    
if __name__ == "__main__":
    main()
