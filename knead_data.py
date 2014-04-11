'''
pipeline.py
Author: Andy Shi

Pipeline for processing metagenomics sequencing data
'''
import argparse
import subprocess

# Global configuration options

# path to Trimmomatic executable
path_to_trim = "/n/sw/centos6/Trimmomatic-0.30/trimmomatic-0.30.jar"

def trim(infile, trimlen=60, prefix=None):
    '''
    input: 
        infile:     input fastq file
    optional input: 
        trimlen:    length to trim
        prefix:     output prefix
    output: At most 5 files
        (1) prefix.trimmed.1.fastq: trimmed first pair-end file
        (2) prefix.trimmed.2.fastq: trimmed second pair-end file
        (3) prefix.trimmed.single.1.fastq: trimmed sequences from the first file
        that lost their partner from the second file
        (4) prefix.trimmed.single.2.fastq: trimmed sequences from the second file
        that lost their partner from the first file
    Uses trimmomatic to trim reads to [arg] base pairs long
    '''
    if not prefix:
        prefix = infile[0]
        
    cmd = str("java -Xmx8g " + path_to_trim + " PE -phred33 " + infile[0] + 
            " " + infile[1] + " " + prefix + ".trimmed.1.fastq " + prefix + 
            ".trimmed.single.1.fastq " + prefix + ".trimmed.2.fastq " + 
            prefix + ".trimmed.single.2.fastq " + "MINLEN:" + str(trimlen) + 
            " &> " + prefix + ".log")
    print(cmd)

def tag():
    '''
    input:
    output:
    Uses BMTagger to tag and potentially remove unwanted reads
    '''
    pass

def main():
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input FASTQ file")
    parser.add_argument("--infile2", help="input FASTQ file mate")
    parser.add_argument("--trimlen", help="length to trim reads", default=60)

    args = parser.parse_args()
    print(args.trimlen)

    files = [args.infile]
    if args.infile2:
        files.append(args.infile2)

    trim(files, args.trimlen)
    tag()

if __name__ == '__main__':
    main()
