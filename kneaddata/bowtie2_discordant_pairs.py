#!/usr/bin/env python

"""
KneadData bowtie2 discordant pairs

This script wraps bowtie2 to allow for discordant paired alignments output to fastq.

Dependencies: Bowtie2

Copyright (c) 2017 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import sys
import os

try:
    import argparse
except ImportError:
    sys.exit("Please upgrade to python 2.7+")
    
import string
import tempfile
import subprocess

try:
    from kneaddata import utilities
except ImportError:
    sys.exit("Please install kneaddata")
    
def reverse_complement(sequence):
    table=string.maketrans("ATCG","TAGC")
    return sequence.translate(table)[::-1]

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
    description= "Kneaddata bowtie2 discordant pairs\n",
    formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument(
        "-1",
        dest="pair1",
        help="the fastq file of pair1 reads",
        required=True)
    parser.add_argument(
        "-2",
        dest="pair2",
        help="the fastq file of pair2 reads",
        required=True) 
    parser.add_argument(
        "-x",
        dest="index",
        help="the database index file",
        required=True)
    parser.add_argument(
        "--un-pair",
        help="the name of the output files for the paired reads without any alignments",
        required=True)
    parser.add_argument(
        "--al-pair",
        help="the name of the output files for the paired reads with concordant alignments",
        required=True)
    parser.add_argument(
        "--un-single",
        help="the name of the output files for the orphan reads without alignments",
        required=True)
    parser.add_argument(
        "--al-single",
        help="the name of the output files for the orphan reads with alignments",
        required=True)
    parser.add_argument(
        "-U",
        dest="orphan",
        help="the fastq files of orphan reads in comma-delimited list")
    parser.add_argument(
        "-S",
        dest="sam",
        help="the file to write the sam output")
    parser.add_argument(
        "--bowtie2",
        help="the path to the bowtie2 executable",
        default="bowtie2")
    parser.add_argument(
        "--threads",
        help="the number of threads to use",
        default=1)
    parser.add_argument(
        "--bowtie2-options",
        help="the bowtie2 options to apply")
    parser.add_argument(
        "--cat-pairs",
        action="store_true",
        help="concatenate pair files before aligning so reads are aligned as single end")
    parser.add_argument(
        "--reorder",
        action="store_true",
        help="print the sequences in the same order as the input files")
    
    return parser.parse_args()


def run_bowtie2(bowtie2_path,pair1,pair2,orphans,database,sam,threads,options,cat_pairs,reorder):
    """ Run bowtie2 with the options provided """

    command=[bowtie2_path]
    # if pairs are to be run as single end, then provide them all as orphans
    if cat_pairs:
        if orphans:
            orphans=",".join([pair1,pair2,orphans])
        else:
            orphans=",".join([pair1,pair2])
    else:
        command+=["-1",pair1,"-2",pair2]
    
    command+=["--threads",str(threads),"-x",database,"-S",sam,"--no-head"]
    if orphans:
        command+=["-U",orphans]
    if options:
        command+=utilities.format_options_to_list([options])
    if reorder:
        command+=["--reorder"]

    try:
        return_code=subprocess.check_call(command)
    except (EnvironmentError, subprocess.CalledProcessError) as e:
        message="Unable to run bowtie2: " +" ".join(command)
        if hasattr(e, 'output') and e.output:
            message+="\nError message returned from bowtie2:\n" + e.output
        sys.exit(message)
        
def organize_alignments_as_paired(sam,open_files,counts):
    """ Organize the alignments that were generated running as paired end reads """
    
    # read through the same file, writing reads to their corresponding files
    with open(sam) as file_handle:
        for line in file_handle:
            data=line.rstrip().split("\t")
            flag=int(data[1])
            
            # check for the pair information
            pair1=True
            if flag & 128:
                pair1=False
            elif data[0].endswith("2"):
                # the orphan reads do not have pair information
                # so need to check the read identifier for pair information
                pair1=False

            # check if the read aligned
            aligned=True
            if flag & 4:
                aligned=False

            # get the reverse complement of the sequence
            # if this aligns to the reverse reference strand
            if flag & 16:
                data[9]=reverse_complement(data[9])
                data[10]=data[10][::-1]

            # determine which file to write the sequence
            optional_fields=" ".join(data[11:])
            if aligned:
                # check if concordant alignment
                if flag & 2:
                    file_name = "pair1_aligned" if pair1 else "pair2_aligned"
                else:
                    file_name = "orphan1_aligned" if pair1 else "orphan2_aligned"
            else:
                # check if both in the pair do not have alignments
                if flag & 8:
                    file_name = "pair1_unaligned" if pair1 else "pair2_unaligned"
                else:
                    file_name = "orphan1_unaligned" if pair1 else "orphan2_unaligned"

            # write the read to the file
            open_files[file_name].write("\n".join(["@"+data[0],data[9],"+",data[10]])+"\n")
            # increase the count
            counts[file_name]=counts[file_name]+1  
    
def organize_alignments_single(sam,open_files,counts):
    """ Organize the alignments that were generated running the pairs as single end reads """
    
    # read through the sam file, organizing reads by those that aligned 
    aligned={}
    unaligned={}
    with open(sam) as file_handle:
        for line in file_handle:
            data=line.rstrip().split("\t")
            flag=int(data[1])
            
            query_id = data[0][:-1]
            pair1 = False if data[0][-1] == "2" else True
            
            # check if the read aligned
            if flag & 4:
                # this read did not align to the reference
                if not query_id in unaligned:
                    unaligned[query_id]=set()
                unaligned[query_id].add(pair1)
            else:
                # this read aligned to the reference
                if not query_id in aligned:
                    aligned[query_id]=set()
                aligned[query_id].add(pair1)
                
    # read through the sam file again to write the reads to the output files
    with open(sam) as file_handle:
        for line in file_handle:
            data=line.rstrip().split("\t")
            
            query_id = data[0][:-1]
            pair1 = False if data[0][-1] == "2" else True
            
            # check the alignment type of this query
            align_values = aligned.get(query_id,[])
            unalign_values = unaligned.get(query_id,[])
            
            if len(align_values) > 1:
                # both reads in the pair aligned to the reference
                file_name = "pair1_aligned" if pair1 else "pair2_aligned"
            elif len(unalign_values) > 1:
                # both reads did not align to the reference
                file_name = "pair1_unaligned" if pair1 else "pair2_unaligned"
            elif pair1 in align_values:
                # only this read from the pair aligned to the reference
                file_name = "orphan1_aligned" if pair1 else "orphan2_aligned"
            elif pair1 in unalign_values:
                # only this read from the pair did not align to the reference
                file_name = "orphan1_unaligned" if pair1 else "orphan2_unaligned"

            # write the read to the file
            open_files[file_name].write("\n".join(["@"+data[0],data[9],"+",data[10]])+"\n")
            # increase the count
            counts[file_name]=counts[file_name]+1                
    

def process_alignments(sam,aligned_pair,unaligned_pair,aligned_orphan,unaligned_orphan,cat_pairs):
    """ Read through the sam alignments and organize into the output files """
    
    # open the output files
    pair1_aligned=open(aligned_pair.replace("%","1"),"wt")
    pair2_aligned=open(aligned_pair.replace("%","2"),"wt")
    
    pair1_unaligned=open(unaligned_pair.replace("%","1"),"wt")
    pair2_unaligned=open(unaligned_pair.replace("%","2"),"wt")
    
    orphan1_aligned=open(aligned_orphan.replace("%","1"),"wt")
    orphan2_aligned=open(aligned_orphan.replace("%","2"),"wt")
    
    orphan1_unaligned=open(unaligned_orphan.replace("%","1"),"wt")
    orphan2_unaligned=open(unaligned_orphan.replace("%","2"),"wt")  
    
    open_files={"pair1_aligned":pair1_aligned,"pair2_aligned":pair2_aligned,
                "pair1_unaligned":pair1_unaligned,"pair2_unaligned":pair2_unaligned,
                "orphan1_aligned":orphan1_aligned,"orphan2_aligned":orphan2_aligned,
                "orphan1_unaligned":orphan1_unaligned,"orphan2_unaligned":orphan2_unaligned}
    counts={name:0 for name in open_files.keys()}  
    
    if cat_pairs:
        organize_alignments_single(sam,open_files,counts)
    else:
        organize_alignments_as_paired(sam,open_files,counts)
            
    # close all of the files
    for file_name, file_handle in open_files.items():
        file_handle.close()
        
    # write out the counts for each file
    for file_name,total in counts.items():
        print(file_name+" : "+str(total))
    

def main():
    # parse the command line arguments
    args = parse_arguments(sys.argv)
    
    # if no sam output is provided or it is set to dev/null, write to temp file
    output_dir=os.path.dirname(args.un_pair)
    temp_files=[]
    if args.sam is None or args.sam == os.devnull:
        file_out, args.sam = tempfile.mkstemp("kneaddata_", "_temp.sam", dir=output_dir)
        os.close(file_out)
        temp_files.append(args.sam)
    
    # run bowtie2
    run_bowtie2(args.bowtie2,args.pair1,args.pair2,args.orphan,args.index,args.sam,args.threads,args.bowtie2_options,args.cat_pairs,args.reorder)
    
    # write output files
    process_alignments(args.sam,args.al_pair,args.un_pair,args.al_single,args.un_single,args.cat_pairs)
    
    # remove the temp files
    for file in temp_files:
        utilities.remove_file(file)

if __name__ == "__main__":
    main()
    
