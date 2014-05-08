'''
pipeline.py
Author: Andy Shi

Pipeline for processing metagenomics sequencing data
'''
import argparse
import subprocess
import shlex

# Global configuration options

# path to Trimmomatic executable
path_to_trim = "/n/sw/centos6/Trimmomatic-0.32/trimmomatic-0.32.jar"
#path_to_trim = "/home/andy/bin/Trimmomatic-0.32/trimmomatic-0.32.jar"

def trim(infile, trimlen=60, prefix=None):
    '''
    input: 
        infile:     input fastq file list (either length 1 or length 2)
                    length 1: single end
                    length 2: paired end
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

    # check that we have the right number of input fastqs
    assert(len(infile) <= 2 and len(infile) >= 1)

    if not prefix:
        prefix = infile[0]
        
    trim_arg = ""
    if len(infile) == 1:
        trim_arg = str(path_to_trim + " SE -phred33 " + infile[0] + " " + prefix
                + ".trimmed.fastq " + "MINLEN:" + str(trimlen) + " &> " + prefix
                + ".log")
    else:
        trim_arg = str(path_to_trim + " PE -phred33 " + infile[0] + " " +
                infile[1] + " " + prefix + ".trimmed.1.fastq " + prefix +
                ".trimmed.single.1.fastq " + prefix + ".trimmed.2.fastq " +
                prefix + ".trimmed.single.2.fastq " + "MINLEN:" + str(trimlen))
        #trim_arg = str(path_to_trim + " PE -phred33 " + infile[0] + " " +
        #        infile[1] + " " + prefix + ".trimmed.1.fastq " + prefix +
        #        ".trimmed.single.1.fastq " + prefix + ".trimmed.2.fastq " +
        #        prefix + ".trimmed.single.2.fastq " + "MINLEN:" + str(trimlen) +
        #        " &> " + prefix + ".log")

    print("Trimmomatic arguments: " + trim_arg)
    cmd = 'java -Xmx500m -jar ' + trim_arg
    print(cmd)
    #subprocess.call(['java ', '-Xmx8g ', '-jar ', trim_arg], shell=True)
    args = shlex.split(cmd)
    print(args)
    subprocess.call(args)

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
    parser.add_argument("infile1", help="input FASTQ file")
    parser.add_argument("-2", "--infile2", help="input FASTQ file mate")
    parser.add_argument("--trimlen", help="length to trim reads", default=60)

    args = parser.parse_args()

    files = [args.infile1]
    if args.infile2:
        files.append(args.infile2)

    print("Running Trimmomatic...")
    trim(files, args.trimlen)
    tag()

if __name__ == '__main__':
    main()
