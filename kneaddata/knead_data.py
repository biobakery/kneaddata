#!/usr/bin/env python

"""
KneadData

KneadData is a tool designed to perform quality control on metagenomic 
sequencing data, especially data from microbiome experiments. In these 
experiments, samples are typically taken from a host in hopes of learning 
something about the microbial community on the host. However, metagenomic 
sequencing data from such experiments will often contain a high ratio of host 
to bacterial reads. This tool aims to perform principled in silico separation 
of bacterial reads from these "contaminant" reads, be they from the host, 
from bacterial 16S sequences, or other user-defined sources.

Dependencies: Trimmomatic, Bowtie2 or BMTagger, and TRF (optional)

To Run: kneaddata -i <input.fastq> -o <output_dir>

Copyright (c) 2015 Harvard School of Public Health

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

# check for the required python version
required_python_version_major = 2
required_python_version_minor = 7
    
try:
    if (sys.version_info[0] < required_python_version_major or
        (sys.version_info[0] == required_python_version_major and
        sys.version_info[1] < required_python_version_minor)):
        sys.exit("CRITICAL ERROR: The python version found (version "+
            str(sys.version_info[0])+"."+str(sys.version_info[1])+") "+
            "does not match the version required (version "+
            str(required_python_version_major)+"."+
            str(required_python_version_minor)+"+)")
except (AttributeError,IndexError):
    sys.exit("CRITICAL ERROR: The python version found (version 1) " +
        "does not match the version required (version "+
        str(required_python_version_major)+"."+
        str(required_python_version_minor)+"+)")  

import os
import logging
import argparse
import re
import itertools

# Try to load one of the kneaddata modules to check the installation
try:
    from kneaddata import utilities
except ImportError:
    sys.exit("ERROR: Unable to find the kneaddata python package." +
        " Please check your install.")

from kneaddata import run
from kneaddata import config

VERSION="0.12.0"

# name global logging instance
logger=logging.getLogger(__name__)
# Global input files path list for FASTQC
original_input_files=[]

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "KneadData\n",
        formatter_class=argparse.RawTextHelpFormatter,
        prog="kneaddata")
    group1 = parser.add_argument_group("global options")
    group1.add_argument(
        "--version",
        action="version",
        version="%(prog)s v"+VERSION)
    parser.add_argument(
        "-v","--verbose",
        action="store_true",
        help="additional output is printed\n")
    group1.add_argument(
        "-i1", "--input1",
        help="Pair 1 input FASTQ file", 
        dest='input1')
    group1.add_argument(
        "-i2", "--input2",
        help="Pair 2 input FASTQ file", 
        dest='input2')
    group1.add_argument(
        "-un","--unpaired",
        help="unparied input FASTQ file", 
        dest='unpaired')
    group1.add_argument(
        "-o", "--output",
        dest='output_dir',
        help="directory to write output files",
        required=True)
    group1.add_argument(
        "-s", "--scratch",
        dest='scratch_dir',
        help="directory to write temp files",
        default="")
    group1.add_argument(
        "-db", "--reference-db",
        default=[], action="append",
        help="location of reference database (additional arguments add databases)")
    group1.add_argument(
        "--bypass-trim",
        action="store_true",
        help="bypass the trim step")
    group1.add_argument(
        "--output-prefix",
        help="prefix for all output files\n[ DEFAULT : $SAMPLE_kneaddata ]")
    group1.add_argument(
        "-t","--threads",
        type=int, 
        default=config.threads,
        metavar="<" + str(config.threads) + ">",  
        help="number of threads\n[ Default : "+str(config.threads)+" ]")
    group1.add_argument(
        "-p","--processes",
        type=int, 
        default=config.processes,
        metavar="<" + str(config.processes) + ">",  
        help="number of processes\n[ Default : "+str(config.processes)+" ]")
    group1.add_argument(
        "-q","--quality-scores",
        default=config.quality_scores,
        choices=config.quality_scores_options,
        dest='trimmomatic_quality_scores',
        help="quality scores\n[ DEFAULT : "+config.quality_scores+" ]")
    group1.add_argument(
        "--run-bmtagger",
        default=False,
        action="store_true",
        dest='bmtagger',
        help="run BMTagger instead of Bowtie2 to identify contaminant reads")
    group1.add_argument(
        "--bypass-trf",
        action="store_true",
        help="option to bypass the removal of tandem repeats")
    group1.add_argument(
        "--run-trf",
        action="store_true",
        help="legacy option to run the removal of tandem repeats (now run by default)")
    group1.add_argument(
        "--run-fastqc-start",
        default=False,
        dest='fastqc_start',
        action="store_true",
        help="run fastqc at the beginning of the workflow")
    group1.add_argument(
        "--run-fastqc-end",
        default=False,
        dest='fastqc_end',
        action="store_true",
        help="run fastqc at the end of the workflow")
    group1.add_argument(
        "--store-temp-output",
        action="store_true",
        help="store temp output files\n[ DEFAULT : temp output files are removed ]")
    group1.add_argument(
        "--remove-intermediate-output",
        action="store_true",
        help="remove intermediate output files\n[ DEFAULT : intermediate output files are stored ]")
    group1.add_argument(
        "--cat-final-output",
        action="store_true",
        help="concatenate all final output files\n[ DEFAULT : final output is not concatenated ]")
    group1.add_argument(
        "--log-level",
        default=config.log_level,
        choices=config.log_level_choices,
        help="level of log messages\n[ DEFAULT : "+config.log_level+" ]")
    group1.add_argument(
        "--log",
        help="log file\n[ DEFAULT : $OUTPUT_DIR/$SAMPLE_kneaddata.log ]")

    group2 = parser.add_argument_group("trimmomatic arguments")
    group2.add_argument(
        "--trimmomatic",
        dest='trimmomatic_path',
        help="path to trimmomatic\n[ DEFAULT : $PATH ]")
    group2.add_argument(
        "--run-trim-repetitive",
        default=False,
        dest='run_trim_repetitive',
        action="store_true",
        help="Trim fastqc generated overrepresented sequences\n")
    group2.add_argument(
        "--max-memory",
        default=config.trimmomatic_memory, 
        help="max amount of memory\n[ DEFAULT : "+config.trimmomatic_memory+" ]")
    group2.add_argument(
        "--trimmomatic-options",
        action="append",
        help="options for trimmomatic\n[ DEFAULT : "+" ".join(utilities.get_default_trimmomatic_options())+" ]\n"+\
             "MINLEN is set to "+str(config.trimmomatic_min_len_percent)+" percent of total input read length. The user can alternatively specify a length (in bases) for MINLEN.")
    group2.add_argument(
        "--sequencer-source",
        dest='sequencer_source',
        default=config.trimmomatic_provided_sequencer_default, 
        choices=config.trimmomatic_provided_sequencer_source,
        help="options for sequencer-source\n[ DEFAULT : "+config.trimmomatic_provided_sequencer_default+"]")

    group3 = parser.add_argument_group("bowtie2 arguments")
    group3.add_argument(
        "--bowtie2",
        dest='bowtie2_path',
        help="path to bowtie2\n[ DEFAULT : $PATH ]")
    group3.add_argument(
        "--bowtie2-options",
        action="append",
        help="options for bowtie2\n[ DEFAULT : "+ " ".join(config.bowtie2_options)+" ]")
    group3.add_argument(
        "--decontaminate-pairs",
        choices=["strict","lenient","unpaired"],
        default="strict",
        help="options for filtering of paired end reads (strict='remove both R1+R2 if either align', lenient='remove only if both R1+R2 align', unpaired='ignore pairing and remove as single end')\n"+\
             "[ DEFAULT : %(default)s ]")
    group3.add_argument(
        "--reorder",
        action="store_true",
        help="order the sequences in the same order as the input\n[ DEFAULT : Sequences are not ordered ]")
    group3.add_argument(
        "--serial",
        action="store_true",
        help="filter the input in serial for multiple databases so a subset of reads are processed in each database search")
        
    group4 = parser.add_argument_group("bmtagger arguments")
    group4.add_argument(
        "--bmtagger",
        dest='bmtagger_path',
        help="path to BMTagger\n[ DEFAULT : $PATH ]")

    group5 = parser.add_argument_group("trf arguments")
    group5.add_argument(
            "--trf",
            dest='trf_path',
            help="path to TRF\n[ DEFAULT : $PATH ]")
    group5.add_argument(
            "--match", 
            type=int,
            default=config.trf_match, 
            help="matching weight\n[ DEFAULT : "+str(config.trf_match)+" ]")
    group5.add_argument(
            "--mismatch",
            type=int,
            default=config.trf_mismatch, 
            help="mismatching penalty\n[ DEFAULT : "+str(config.trf_mismatch)+" ]")
    group5.add_argument(
            "--delta",
            type=int,
            default=config.trf_delta, 
            help="indel penalty\n[ DEFAULT : "+str(config.trf_delta)+" ]")
    group5.add_argument(
            "--pm",
            type=int,
            default=config.trf_match_probability,
            help="match probability\n[ DEFAULT : "+str(config.trf_match_probability)+" ]")
    group5.add_argument(
            "--pi",
            type=int,
            default=config.trf_pi,
            help="indel probability\n[ DEFAULT : "+str(config.trf_pi)+" ]")
    group5.add_argument(
            "--minscore",
            type=int,
            default=config.trf_minscore, 
            help="minimum alignment score to report\n[ DEFAULT : "+str(config.trf_minscore)+" ]")
    group5.add_argument(
            "--maxperiod",
            type=int,
            default=config.trf_maxperiod, 
            help="maximum period size to report\n[ DEFAULT : "+str(config.trf_maxperiod)+" ]")
    group6 = parser.add_argument_group("fastqc arguments")
    group6.add_argument(
            "--fastqc",
            dest='fastqc_path',
            help="path to fastqc\n[ DEFAULT : $PATH ]")

    return parser.parse_args()
    
def update_configuration(args):
    """ Update the run settings based on the arguments provided """

    # get the full path for the output directory
    args.output_dir = os.path.abspath(args.output_dir)
    if args.scratch_dir:
        args.scratch_dir = os.path.abspath(args.scratch_dir)    

    # set if temp output should be removed
    args.remove_temp_output = not args.store_temp_output
    
    # if intermediate output should be removed, then also remove temp output
    if args.remove_intermediate_output:
        args.remove_temp_output = True
    
    # check the input files are non-empty and readable
    args.input=[]
    if (args.input1 and args.input2):
        args.input.append(os.path.abspath(args.input1))
        args.input.append(os.path.abspath(args.input2))
    if (args.unpaired):
        args.input.append(os.path.abspath(args.unpaired))  
    utilities.is_file_readable(args.input[0],exit_on_error=True)
    if len(args.input) == 2:
        utilities.is_file_readable(args.input[1],exit_on_error=True)
    elif len(args.input) > 2:
        sys.exit("ERROR: Please provide at most 2 input files.")
    elif len(args.input) == 0:
        sys.exit("ERROR: Please provide --input1/--input2 or --unpaired (input) files.")
    
    #Store original file paths for FASTQC 
    for input in args.input:
        original_input_files.append(input)
    
    # create the output directory and scratch if needed
    utilities.create_directory(args.output_dir)
    if args.scratch_dir:
        utilities.create_directory(args.scratch_dir)
    
    # set bowtie2 options
    if args.bowtie2_options:
        # parse the options from the user into any array of options
        args.bowtie2_options=utilities.format_options_to_list(args.bowtie2_options)
    else:
        # if not set by user, then set to default options
        args.bowtie2_options = config.bowtie2_options
        
    # add the quality scores to the bowtie2 options
    args.bowtie2_options+=[config.bowtie2_flag_start+args.trimmomatic_quality_scores]    
  
    # set the mode for single end input file
    if len(args.input) == 1:
        args.decontaminate_pairs = "unpaired"
 
    # set the bowtie2 mode based on the pairs input
    args.discordant = False
    if args.decontaminate_pairs != "lenient" :
        args.discordant = True
 
    # update the quality score option into a flag for trimmomatic
    args.trimmomatic_quality_scores=config.trimmomatic_flag_start+args.trimmomatic_quality_scores
        
    # find the location of trimmomatic, trimmomatic does not need to be executable
    if not args.bypass_trim:
        args.trimmomatic_path=utilities.find_dependency(args.trimmomatic_path,config.trimmomatic_jar,"trimmomatic",
            "--trimmomatic", bypass_permissions_check=True)
    
    # find the location of bmtagger, if set to run
    if args.reference_db:
        if args.bmtagger:
            args.bmtagger_path=utilities.find_dependency(args.bmtagger_path,config.bmtagger_exe,"bmtagger",
                "--bmtagger", bypass_permissions_check=False)
            # add this folder to path, so as to be able to find other dependencies like bmfilter
            utilities.add_exe_to_path(os.path.dirname(args.bmtagger_path))
        else:
            # find the location of bowtie2, if not running with bmtagger
            args.bowtie2_path=utilities.find_dependency(args.bowtie2_path, config.bowtie2_exe, "bowtie2",
                "--bowtie2", bypass_permissions_check=False)        
    
    # find the location of trf, if set to run
    if not args.bypass_trf:
        args.trf_path=utilities.find_dependency(args.trf_path,config.trf_exe,"trf",
            "--trf", bypass_permissions_check=False)
        
    # if fastqc is set to be run, check if the executable can be found
    if args.fastqc_start or args.fastqc_end:
        args.fastqc_path=utilities.find_dependency(args.fastqc_path,config.fastqc_exe,"fastqc",
                                                   "--fastqc",bypass_permissions_check=False)

    # set the default output prefix 
    if args.output_prefix == None:
        if args.input[0].endswith(".gz") or args.input[0].endswith(".bz2"):
            # remove compression extension if present
            infile_base = os.path.splitext(os.path.splitext(os.path.basename(args.input[0]))[0])[0]
        else:
            infile_base = os.path.splitext(os.path.basename(args.input[0]))[0]
        args.output_prefix = infile_base + "_kneaddata"    

    # find the bowtie2 indexes for each of the reference databases
    # reference database inputs can be directories, indexes, or index files
    if args.reference_db:
        reference_indexes=[]
        database_type="bowtie2"
        if args.bmtagger:
            database_type="bmtagger"
        for directory in args.reference_db:
            reference_indexes.append(utilities.find_database_index(os.path.abspath(directory),database_type))
    
        args.reference_db=reference_indexes
    
    return args

def setup_logging(args):
    """ Set up the log file """
    
    if not args.log:
        args.log = os.path.join(args.output_dir,args.output_prefix+".log")

    # configure the logger
    logging.basicConfig(filename=args.log,format='%(asctime)s - %(name)s - %(levelname)s: %(message)s',
        level=getattr(logging,args.log_level), filemode='w', datefmt='%m/%d/%Y %I:%M:%S %p')
    
    # write the version of the software to the log
    logger.info("Running kneaddata v"+VERSION)
    
    # write the location of the output files to the log
    message="Output files will be written to: " + args.output_dir
    logger.info(message)
    
    # write out all of the argument settings
    message="Running with the following arguments: \n"
    for key,value in vars(args).items():
        if isinstance(value,list) or isinstance(value,tuple):
            value_string=" ".join([str(i) for i in value])
        else:
            value_string=str(value)
        
        message+=key+" = "+value_string+"\n"
    logger.debug(message)

def main():
    # Parse the arguments from the user
    args = parse_arguments(sys.argv)
    
    # Update the configuration
    args = update_configuration(args)
    
    # Start logging
    setup_logging(args)

    # set the prefix for the output files
    final_output_dir = args.output_dir
    if args.scratch_dir:
        full_path_output_prefix = os.path.join(args.scratch_dir, args.output_prefix)
        args.output_dir = args.scratch_dir
    else:
        full_path_output_prefix = os.path.join(args.output_dir, args.output_prefix)

    temp_output_files=[]
    # Check for compressed files, bam files, or sam files
    for index in range(len(args.input)):

        # check for gzipped/bz2 files
        if args.input[index].endswith(".gz") or args.input[index].endswith(".bz2"):
            args.input[index]=utilities.get_decompressed_file(args.input[index], args.output_dir, temp_output_files, args.input)
        elif args.input[index].endswith(".bam"):
            input_files_set=utilities.get_fastq_from_bam_file(args.input[index], args.output_dir, temp_output_files, args.input)
            if isinstance(input_files_set,list):
                args.input=input_files_set
            else:
                args.input[index]=input_files_set
        elif args.input[index].endswith(".sam"): 
            args.input[index]=utilities.get_fastq_from_sam_file(args.input[index], args.output_dir, temp_output_files, args.input)
        
    # Get the format of the first input file
    file_format=utilities.get_file_format(args.input[0])

    if file_format != "fastq":
        message="Your input file is of type: "+file_format+". Please provide an input file of fastq format."
        logger.critical(message)
        sys.exit(message)
    
    # if this is the new illumina identifier format, create temp files after reformatting the headers
    for index in range(len(args.input)):
        args.input[index]=utilities.get_reformatted_identifiers(args.input[index],index,args.output_dir, temp_output_files, args.input)
    
    # check for reads that are not ordered and order if needed (if trimmomatic is run)
    if not args.bypass_trim and len(args.input)==2:
        args.input=utilities.check_and_reorder_reads(args.input, args.output_dir, temp_output_files)
   
    # remove any temp files from decompress/reformat that are no longer needed
    utilities.update_temp_output_files(temp_output_files, [], args.input)
 
    # set trimmomatic options
    # this is done after the decompression and conversions from sam/bam
    # as the default requires the read length from the input sequences
    if args.trimmomatic_options:
        # parse the options from the user into an array of options
        args.trimmomatic_options = utilities.format_options_to_list(args.trimmomatic_options)
    else:
        # if trimmomatic options not set by user, then set to default options
        # use read length of input file for minlen
        args.trimmomatic_options = utilities.get_default_trimmomatic_options(utilities.get_read_length_fastq(args.input[0]),
            path=config.trimmomatic_adapter_folder,type="PE" if len(args.input) == 2 else "SE", sequencer_source=args.sequencer_source)
            

    # Get the number of reads initially
    utilities.log_read_count_for_files(args.input,"raw","Initial number of reads",args.verbose)
    
    # Run fastqc if set to run at start of workflow
    if args.fastqc_start or args.run_trim_repetitive:
        run.fastqc(args.fastqc_path, args.output_dir, original_input_files, args.threads, args.verbose)
        #Setting fastqc output zip and txt file path
        output_txt_files=[]
        for input_file_name in original_input_files:
            temp_file = os.path.splitext(input_file_name)[0]
            if (temp_file.count('fastq')>0 or temp_file.count('fq')>0 ):
                temp_file = os.path.splitext(temp_file)[0]
            output_txt_files.append(args.output_dir+"/fastqc/"+temp_file.split('/')[-1]+"_fastqc/fastqc_data.txt")

    if not args.bypass_trim:
        if args.run_trim_repetitive:
             # Get the Min Overrepresented Seq Length
            args.trimmomatic_options = utilities.get_updated_trimmomatic_parameters(output_txt_files, args.output_dir, args.trimmomatic_options)
        
        trimmomatic_output_files = run.trim(
        args.input, full_path_output_prefix, args.trimmomatic_path, 
        args.trimmomatic_quality_scores, args.max_memory, args.trimmomatic_options, 
        args.threads, args.verbose)
        
    else:
        message="Bypass trimming"	
        logger.info(message)	
        print(message)	
        trimmomatic_output_files=[args.input]
        
    # Get the number of reads after trimming
    utilities.log_read_count_for_files(trimmomatic_output_files,"trimmed","Total reads after trimming",args.verbose)
   
    # run TRF, if set
    if not args.bypass_trf:
        # run trf on all output files
        trf_output_files=run.tandem(trimmomatic_output_files, full_path_output_prefix, args.match,
                                      args.mismatch,args.delta,args.pm,args.pi,
                                      args.minscore,args.maxperiod,args.trf_path,
                                      args.processes,args.verbose,args.remove_temp_output,args.threads)
        # remove the aligment files, if intermediate output files should be removed
        if args.reference_db and args.remove_intermediate_output:
            temp_output_files+=utilities.resolve_sublists(trimmomatic_output_files)
    else:
        trf_output_files = trimmomatic_output_files
    # If a reference database is not provided, then bypass decontamination step
    if not args.reference_db:
        message="Bypass decontamination"
        logger.info(message)
        print(message)
        # resolve sub-lists if present
        final_output_files=trf_output_files
    else:
        final_output_files=run.decontaminate(args, full_path_output_prefix, trf_output_files)
        # remove trimmed output files, if set to remove intermediate outputx
        if not args.bypass_trim and args.remove_intermediate_output:
            temp_output_files+=utilities.resolve_sublists(trf_output_files)
        
    # If set, concat the final output files if there is more than one
    final_output_files = utilities.resolve_sublists(final_output_files)
    if args.cat_final_output and len(final_output_files) > 1:
        cat_output_file=full_path_output_prefix+config.fastq_file_extension
        utilities.cat_files(final_output_files,cat_output_file)
        
        # if removing intermediate output, then remove the files that were merged
        if args.remove_intermediate_output:
            temp_output_files+=final_output_files
            final_output_files=[cat_output_file]
        else:
            final_output_files.append(cat_output_file)
        
    # Remove any temp output files, if set
    if not args.store_temp_output:
        for file in temp_output_files:
            utilities.remove_file(file)
            
    # Run fastqc if set to run at end of workflow
    if args.fastqc_end:
        run.fastqc(args.fastqc_path, args.output_dir, final_output_files, args.threads, args.verbose)

    # If using scratch, then move final output files to output folder
    if args.scratch_dir:
        scratch_output_files=final_output_files
        final_output_files=[]
        for outfile in scratch_output_files:
            utilities.move_file(os.path.basename(outfile),args.output_dir,final_output_dir)
            final_output_files.append(os.path.join(final_output_dir,os.path.basename(outfile)))

    if len(final_output_files) > 1:
        message="\nFinal output files created: \n"
    else:
        message="\nFinal output file created: \n"
    
    message=message+ "\n".join(final_output_files) + "\n"
    logger.info(message)
    print(message)

if __name__ == '__main__':
    main()
