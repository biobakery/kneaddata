#!/usr/bin/env python

'''
knead_data.py
Author: Andy Shi

Pipeline for processing metagenomics sequencing data
'''
import os
import logging
import argparse
import gzip
import re
import sys

from knead_datalib import strategies, try_create_dir, parse_positive_int

VERSION="0.4.2"

def handle_cli():
    """parse command line arguments
    note: argparse converts dashes '-' in argument prefixes to
    underscores '_'
    """
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("global options")
    group1.add_argument(
        "-1", "--infile1",
        help="input FASTQ file", 
        required=True)
    group1.add_argument(
        "-2", "--infile2",
        help="input FASTQ file mate",
        default=None)
    group1.add_argument(
        "-db", "--reference-db",
        default=[], action="append", required=True,
        help=("prefix for reference databases used for either"
              " Bowtie2 or BMTagger"))
    group1.add_argument(
        "-o", "--output-prefix", default="knead_out",
        help="prefix for all output files")
    group1.add_argument(
        "-D", "--output-dir", default=os.getcwd(),
        help="where to put all output files")
    group1.add_argument(
        "--threads",
        type=parse_positive_int, default=None, 
        help=("Maximum number of processes to run."
              " Default uses all but one available CPUs"))
    group1.add_argument(
        "-s", "--strategy",
        default="storage", choices=['memory','storage'],
        help=("Define operating strategy: 'storage' for IO-heavy"
              " or 'memory' for memory-heavy"))
    group1.add_argument(
        '-l', '--logging',
        default="INFO",
        help=("Logging verbosity, options are debug, info, warning,"
              " and critical. If set to debug, temporary files are not"
              " removed"))
    group1.add_argument(
        '--logfile',
        default=None,
        help="Where to save logs")
    group1.add_argument(
        "--version",
        action="version",
        version="%(prog)s v"+VERSION)

    group2 = parser.add_argument_group("trimmomatic arguments")
    group2.add_argument(
        "-t", "--trim-path",
        help="path to Trimmomatic .jar executable")
    group2.add_argument(
        "--trimlen",
        type=parse_positive_int, default=60,
        help="minimum length for a trimmed read in Trimmomatic")
    group2.add_argument(
        "-m", "--max-mem",
        default="500m", 
        help=("Maximum amount of memory that will be used by "
              "Trimmomatic, as a string, ie 500m or 8g"))
    default_trimargs = ["SLIDINGWINDOW:4:20", "MINLEN:60"]
    group2.add_argument(
        "-a", "--trim-args",
        default=[], action="append",
        help=("Additional arguments for Trimmomatic, default: "
              +" ".join(default_trimargs)))

    group3 = parser.add_argument_group("bowtie2 arguments")
    group3.add_argument(
        "--bowtie2-path",
        default=None, help="path to bowtie2 if not found on $PATH")
    default_bowtie2_args = ["--very-sensitive"]
    group3.add_argument(
        "--bowtie2-args",
        default=[], action="append",
        help=("Additional arguments for Bowtie 2, default: " 
                + " ".join(default_bowtie2_args)))
        
    group4 = parser.add_argument_group("bmtagger arguments")
    group4.add_argument(
        "--bmtagger",
        default=False, action="store_true",
        help="If set, use BMTagger to identify contaminant reads")
    group4.add_argument(
        "--extract",
        default=False, action="store_true",
        help=("Only has an effect if --bmtagger is set. If this is set,"
              " knead_data outputs cleaned FASTQs, without contaminant reads."
              " Else, output a list or lists of contaminant reads."))
    group4.add_argument(
        "--bmtagger-path",
        default=None,
        help="path to BMTagger executable if not found in $PATH")

    group5 = parser.add_argument_group("trf arguments")
    group5.add_argument(
            "--trf",
            default=False, action="store_true",
            help="If set, use TRF to remove and/or mask tandem repeats")
    group5.add_argument(
            "--trf-path",
            default="trf",
            help="Path to TRF executable if not found in $PATH")
    group5.add_argument(
            "--match", type=parse_positive_int,
            default=2, 
            help="TRF matching weight. Default: 2")
    group5.add_argument(
            "--mismatch", type=parse_positive_int,
            default=7, 
            help="TRF mismatching penalty. Default: 7")
    group5.add_argument(
            "--delta", type=parse_positive_int,
            default=7, 
            help="TRF indel penalty. Default: 7")
    group5.add_argument(
            "--pm", type=parse_positive_int,
            default=80,
            help="TRF match probability (whole number). Default: 80")
    group5.add_argument(
            "--pi", type=parse_positive_int,
            default=10,
            help="TRF indel probability (whole number). Default: 10")
    group5.add_argument(
            "--minscore", type=parse_positive_int,
            default=50, 
            help="TRF minimum alignment score to report. Default: 50")
    group5.add_argument(
            "--maxperiod", type=parse_positive_int,
            default=500, 
            help="TRF maximum period size to report. Default: 500")
    group5.add_argument(
            "--no-generate-fastq",
            default=True, action="store_false", 
            help="If switched on, don't generate fastq output for trf")
    group5.add_argument(
            "--mask",
            default=False, action="store_true",
            help="If switched on, generate mask file for trf output")
    group5.add_argument(
            "--html",
            default=False, action="store_true",
            help="If switched on, generate html file for trf output")

    args = parser.parse_args()
    try_create_dir(args.output_dir)

    if (not args.no_generate_fastq) and (not args.mask) and args.trf:
        parser.error("\nYou cannot set the --no-generate-fastq flag without"
        " the --mask flag. Exiting...\n")

    # allow users to overwrite the defaults instead of appending to defaults
    if args.trim_args == []:
        args.trim_args = default_trimargs

    if args.bowtie2_args == []:
        args.bowtie2_args = default_bowtie2_args
        
    # find the location of trimmomatic
    trimmomatic_jar="trimmomatic-0.33.jar"
    args.trim_path=check_for_dependency(args.trim_path,trimmomatic_jar,"Trimmomatic",
        "--trim-path", True)
    # add the full path to the jar file
    args.trim_path=os.path.abspath(os.path.join(args.trim_path,trimmomatic_jar))
    
    # check for bowtie2
    bowtie2_exe="bowtie2"
    args.bowtie2_path=check_for_dependency(args.bowtie2_path, bowtie2_exe, "Bowtie2",
        "--bowtie2-path", False)
    # add the full path to bowtie2
    args.bowtie2_path=os.path.abspath(os.path.join(args.bowtie2_path,bowtie2_exe))

    return args

def check_for_dependency(path_provided,exe,name,path_option,bypass_permissions_check):
    """ 
    Check if the dependency can be found in the path provided or in $PATH
    Return the location of the dependency
    """

    found_path=""
    if path_provided:
        path_provided=os.path.abspath(path_provided)
        # check that the exe can be found
        try:
            files=os.listdir(path_provided)
        except EnvironmentError:
            sys.exit("ERROR: Unable to list files in "+name+" directory: "+ path_provided)
            
        if not exe in files:
            sys.exit("ERROR: The "+exe+" executable is not included in the directory: " + path_provided)
            
        found_path=path_provided
    else:
        # search for the exe
        exe_path=find_exe_in_path(exe, bypass_permissions_check)
        if not exe_path:
            sys.exit("ERROR: Unable to find "+name+". Please provide the "+
                "full path to "+name+" with "+path_option+".")
        else:
            found_path=exe_path  
        
    return found_path


def setup_logging(args):
    fmt = "%(asctime)s %(levelname)s: %(message)s"
    if not args.logfile:
        args.logfile = os.path.join(args.output_dir,
                                    args.output_prefix+".log")

    logger = logging.getLogger()
    logger.setLevel(getattr(logging, args.logging.upper()))
    logging.basicConfig(format=fmt)

    filehandler = logging.FileHandler(args.logfile)
    filehandler.setLevel(logging.DEBUG)
    filehandler.setFormatter(logging.Formatter(fmt=fmt))
    logger.addHandler(filehandler)
    
    return args


def get_file_format(file):
    """ Determine the format of the file """

    format="unknown"
    file_handle=None

    # check the file exists and is readable
    if not os.path.isfile(file):
        logging.critical("The input file selected is not a file: %s.",file)

    if not os.access(file, os.R_OK):
        logging.critical("The input file selected is not readable: %s.",file)

    try:
        # check for gzipped files
        if file.endswith(".gz"):
            file_handle = gzip.open(file, "r")
        else:
            file_handle = open(file, "r")

        first_line = file_handle.readline()
        second_line = file_handle.readline()
    except EnvironmentError:
        # if unable to open and read the file, return unknown
        return "unknown"
    finally:
        if file_handle:
            file_handle.close()

    # check that second line is only nucleotides or amino acids
    if re.search("^[A-Z|a-z]+$", second_line):
        # check first line to determine fasta or fastq format
        if re.search("^@",first_line):
            format="fastq"
        if re.search("^>",first_line):
            format="fasta"

    return format

def find_exe_in_path(exe, bypass_permissions_check=None):
    """
    Check that an executable exists in $PATH
    """
    
    paths = os.environ["PATH"].split(os.pathsep)
    for path in paths:
        fullexe = os.path.join(path,exe)
        if os.path.exists(fullexe):
            if bypass_permissions_check or os.access(fullexe,os.X_OK):
                return path
    return None

def main():
    args = handle_cli()
    # check args first

    args = setup_logging(args)

    logging.debug("Running knead_data.py with the following"
                  " arguments (from argparse): %s", str(args))

    # Get the format of the first input file
    file_format=get_file_format(args.infile1)

    if file_format != "fastq":
        logging.critical("Your input file is of type: %s . Please provide an input file of fastq format.",file_format)
        sys.exit(1)

    if args.strategy == 'memory':
        strategies.memory_heavy(args)
    elif args.strategy == 'storage':
        strategies.storage_heavy(args)
        
    logging.info("Done!")


if __name__ == '__main__':
    main()
