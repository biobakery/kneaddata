#!/usr/bin/env python

"""
KneadData parallel TRF

This script wraps bowtie2 to allow for discordant paired alignments output to fastq.

Dependencies: TRF

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

import argparse
    
import string
import tempfile
import subprocess

try:
    from kneaddata import utilities
except ImportError:
    sys.exit("Please install kneaddata")
    
def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
    description= "Kneaddata trf parallel\n",
    formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument(
        "--input",
        help="the fasta file of reads",
        required=True)
    parser.add_argument(
        "--output",
        help="the trf output file",
        required=True)
    parser.add_argument(
        "--trf-path",
        help="full path to the trf executable",
        required=True)
    parser.add_argument(
        "--trf-options",
        help="space delimited list of trf options",
        required=True)
    parser.add_argument(
        "--nproc",
        default=1,
        type=int,
        help="total number of processes to use")
    
    return parser.parse_args()


def run_trf(input,trf_path,trf_options,nproc,output,verbose=True):
    """ Run trf with the options provided """

    tempfile_list=[]
    datfile_list=[]
    commands=[]

    # check for one process and if so just run trf directly
    if nproc == 1:
        commands.append([[trf_path,input]+trf_options.split(" "),"trf",[input],[output],output])
        
        utilities.start_processes(commands,nproc,verbose)
    else:
        # get the total number of reads
        total_lines=0
        with open(input) as file_handle:
            for line in file_handle:
                total_lines+=1

        # split the input into multiple files and run in parallel
        for i in range(int(nproc)):
            file_out, new_file = tempfile.mkstemp(prefix=os.path.basename(output)+'_'+str(i)+'_temp_trf_output',dir=os.path.dirname(output))
            os.close(file_out)
            tempfile_list.append(new_file)
            datfile_list.append(new_file+".".join(trf_options.split("-")[0].split(" "))+"dat")

        # write the input file into all temp output files
        output_file_number=0
        lines_per_file = int(total_lines/int(nproc))
        lines_written=0
        file_handle_write=None
        for read_line in utilities.read_file_n_lines(input,2):
            if not file_handle_write:
                file_handle_write = open(tempfile_list[output_file_number],"wt")
            file_handle_write.write("".join(read_line))

            lines_written+=2
            if lines_written > lines_per_file:
                file_handle_write.close()
                lines_written=0
                output_file_number+=1
                file_handle_write = open(tempfile_list[output_file_number],"wt")
               
        file_handle_write.close() 

        # run commands
        for i, temp_in, temp_out in zip(range(len(tempfile_list)), tempfile_list, datfile_list):
            trf_command=[trf_path,temp_in]+trf_options.split(" ")
            commands.append([trf_command,"trf{}".format(i),[temp_in],[temp_out],temp_out])

        utilities.start_processes(commands,nproc,verbose)
    
        # merge all of the outputs to the final output file
        with open(output,"w") as file_write:
            for datfile in datfile_list:
                with open(datfile) as file_read:
                    for line in file_read:
                        file_write.write(line)
    
        # remove temp files
        for filename in tempfile_list+datfile_list:
            try:
                os.remove(filename)   
            except EnvironmentError:
                print("Unable to remove temp file: " + filename)         
    
def main():
    # parse the command line arguments
    args = parse_arguments(sys.argv)
   
    # remove extra quotes if included
    args.trf_options=args.trf_options.replace("'","")
 
    # add exe name if not included
    if not args.trf_path.endswith("trf"):
        if not args.trf_path.endswith("/"):
            args.trf_path+="/"
        args.trf_path+="trf"

    run_trf(args.input,args.trf_path,args.trf_options,args.nproc,args.output)
    
if __name__ == "__main__":
    main()
    
