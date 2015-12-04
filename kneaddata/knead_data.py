#!/usr/bin/env python

'''
kneaddata.py
Author: Andy Shi

Pipeline for processing metagenomics sequencing data
'''

import sys

# check for the required python version
required_python_version_major = 2
required_python_version_minor = 7
    
try:
    if (sys.version_info[0] != required_python_version_major or
        sys.version_info[1] < required_python_version_minor):
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

# Try to load one of the kneaddata modules to check the installation
try:
    from . import utilities
except ImportError:
    sys.exit("ERROR: Unable to find the kneaddata python package." +
        " Please check your install.")

from . import storageheavy
from . import memoryheavy
from . import config

VERSION="0.4.6.1"

# name global logging instance
logger=logging.getLogger(__name__)

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
        "-i", "--input",
        help="input FASTQ file (add a second argument instance to run with paired input files)", 
        action="append",
        required=True)
    group1.add_argument(
        "-o", "--output",
        dest='output_dir',
        help="directory to write output files",
        required=True)
    group1.add_argument(
        "-db", "--reference-db",
        default=[], action="append",
        help="location of reference database (additional arguments add databases)")
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
        "-s", "--strategy",
        default=config.strategy, choices=config.strategy_choices,
        help="define operating strategy\n[ DEFAULT : "+config.strategy+" ]")
    group1.add_argument(
        "--run-bmtagger",
        default=False,
        action="store_true",
        dest='bmtagger',
        help="run BMTagger instead of Bowtie2 to identify contaminant reads")
    group1.add_argument(
        "--run-trf",
        default=False,
        dest='trf',
        action="store_true",
        help="run TRF to remove and/or mask tandem repeats")
    group1.add_argument(
        "--remove-temp-output",
        action="store_true",
        help="remove temp output files\n[ DEFAULT : temp output files are not removed ]")
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
        "--max-memory",
        default=config.trimmomatic_memory, 
        help="max amount of memory\n[ DEFAULT : "+config.trimmomatic_memory+" ]")
    group2.add_argument(
        "--trimmomatic-options",
        action="append",
        help="options for trimmomatic\n[ DEFAULT : "+" ".join(config.trimmomatic_options)+" ]")

    group3 = parser.add_argument_group("bowtie2 arguments")
    group3.add_argument(
        "--bowtie2",
        dest='bowtie2_path',
        help="path to bowtie2\n[ DEFAULT : $PATH ]")
    group3.add_argument(
        "--bowtie2-options",
        action="append",
        help="options for bowtie2\n[ DEFAULT : "+ " ".join(config.bowtie2_options)+" ]")
        
    group4 = parser.add_argument_group("bmtagger arguments")
    group4.add_argument(
        "--bmtagger",
        dest='bmtagger_path',
        help="path to BMTagger\n[ DEFAULT : $PATH ]")
    group4.add_argument(
        "--extract",
        default=False, action="store_true",
        help="output cleaned FASTQs without contaminant reads\n[ DEFAULT : output lists of contaminant reads ]")

    group5 = parser.add_argument_group("trf arguments")
    group5.add_argument(
            "--trf",
            default="trf",
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
    group5.add_argument(
            "--no-generate-fastq",
            default=True, 
            action="store_false", 
            help="don't generate fastq output from trf")
    group5.add_argument(
            "--mask",
            default=False,
            action="store_true",
            help="generate mask file from trf output")
    group5.add_argument(
            "--html",
            default=False, action="store_true",
            help="generate html file for trf output")

    return parser.parse_args()
    
def update_configuration(args):
    """ Update the run settings based on the arguments provided """

    # get the full path for the output directory
    args.output_dir = os.path.abspath(args.output_dir)
    
    # check the input files are non-empty and readable
    args.input[0] = os.path.abspath(args.input[0])
    utilities.is_file_readable(args.input[0],exit_on_error=True)
    
    if len(args.input) == 2:
        args.input[1] = os.path.abspath(args.input[1])
        utilities.is_file_readable(args.input[1],exit_on_error=True)
    elif len(args.input) > 2:
        sys.exit("ERROR: Please provide at most 2 input files.")
    
    # create the output directory if needed
    utilities.create_directory(args.output_dir)

    if (not args.no_generate_fastq) and (not args.mask) and args.trf:
        parser.error("\nYou cannot set the --no-generate-fastq flag without"
        " the --mask flag. Exiting...\n")

    # set trimmomatic options
    if args.trimmomatic_options:
        # parse the options from the user into an array of options
        args.trimmomatic_options=utilities.format_options_to_list(args.trimmomatic_options)
    else:
        # if not set by user, then set to default options
        args.trimmomatic_options = config.trimmomatic_options

    # set bowtie2 options
    if args.bowtie2_options:
        # parse the options from the user into any array of options
        args.bowtie2_options=utilities.format_options_to_list(args.bowtie2_options)
    else:
        # if not set by user, then set to default options
        args.bowtie2_options = config.bowtie2_options
        
    # add the quality scores to the bowtie2 options
    args.bowtie2_options+=[config.bowtie2_flag_start+args.trimmomatic_quality_scores]    
    
    # update the quality score option into a flag for trimmomatic
    args.trimmomatic_quality_scores=config.trimmomatic_flag_start+args.trimmomatic_quality_scores

    # set the default output prefix 
    if args.output_prefix == None:
        infile_base = os.path.splitext(os.path.basename(args.input[0]))[0]
        args.output_prefix = infile_base + "_kneaddata"
        
    # find the location of trimmomatic
    args.trimmomatic_path=utilities.find_dependency(args.trimmomatic_path,config.trimmomatic_jar,"trimmomatic",
        "--trimmomatic", True)
    
    # find the location of bowtie2
    args.bowtie2_path=utilities.find_dependency(args.bowtie2_path, config.bowtie2_exe, "bowtie2",
        "--bowtie2", False)
    
    # find the location of bmtagger, if set to run
    if args.bmtagger:
        args.bmtagger_path=utilities.find_dependency(args.bmtagger_path,config.bmtagger_exe,"bmtagger",
            "--bmtagger", True)
    
    # find the location of trf, if set to run
    if args.trf:
        args.trf_path=utilities.find_dependency(args.trf_path,config.trf_exe,"trf",
            "--trf", True)
    
    # find the bowtie2 indexes for each of the reference databases
    # reference database inputs can be directories, indexes, or index files
    reference_indexes=set()
    database_type="bowtie2"
    if args.bmtagger:
        database_type="bmtagger"
    for directory in args.reference_db:
        reference_indexes.add(utilities.find_database_index(os.path.abspath(directory),database_type))

    args.reference_db=list(reference_indexes)
    
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

    # Get the format of the first input file
    file_format=utilities.get_file_format(args.input[0])

    if file_format != "fastq":
        message="Your input file is of type: "+file_format+". Please provide an input file of fastq format."
        logger.critical(message)
        sys.exit(message)

    if args.strategy == 'memory':
        memoryheavy.memory_heavy(args)
    elif args.strategy == 'storage':
        storageheavy.storage_heavy(args)

if __name__ == '__main__':
    main()
