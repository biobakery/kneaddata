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
import gzip
import re

# Try to load one of the kneaddata modules to check the installation
try:
    from . import util
except ImportError:
    sys.exit("ERROR: Unable to find the kneaddata python package." +
        " Please check your install.")

from . import storageheavy
from . import memoryheavy

VERSION="0.4.6.1"

def handle_cli():
    """parse command line arguments
    note: argparse converts dashes '-' in argument prefixes to
    underscores '_'
    """
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("global options")
    group1.add_argument(
        "-i", "--input",
        help="input FASTQ file", 
        dest='infile1',
        required=True)
    group1.add_argument(
        "--input2",
        help="input FASTQ file mate",
        dest='infile2',
        default=None)
    group1.add_argument(
        "-db", "--reference-db",
        default=[], action="append",
        help=("prefix for reference databases used for either"
              " Bowtie2 or BMTagger"))
    group1.add_argument(
        "--output-prefix", default=None,
        help="prefix for all output files")
    group1.add_argument(
        "-o", "--output",
        dest='output_dir',
        help="where to write all output files",
        required=True)
    group1.add_argument(
        "--threads",
        type=util.parse_positive_int, default=1, 
        help=("Maximum number of processes to run."
              " Default is 1."))
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
              " kneaddata outputs cleaned FASTQs, without contaminant reads."
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
            "--match", type=util.parse_positive_int,
            default=2, 
            help="TRF matching weight. Default: 2")
    group5.add_argument(
            "--mismatch", type=util.parse_positive_int,
            default=7, 
            help="TRF mismatching penalty. Default: 7")
    group5.add_argument(
            "--delta", type=util.parse_positive_int,
            default=7, 
            help="TRF indel penalty. Default: 7")
    group5.add_argument(
            "--pm", type=util.parse_positive_int,
            default=80,
            help="TRF match probability (whole number). Default: 80")
    group5.add_argument(
            "--pi", type=util.parse_positive_int,
            default=10,
            help="TRF indel probability (whole number). Default: 10")
    group5.add_argument(
            "--minscore", type=util.parse_positive_int,
            default=50, 
            help="TRF minimum alignment score to report. Default: 50")
    group5.add_argument(
            "--maxperiod", type=util.parse_positive_int,
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
    
    # get the full path for the output directory
    args.output_dir = os.path.abspath(args.output_dir)
    
    # create the output directory if needed
    util.try_create_dir(args.output_dir)

    if (not args.no_generate_fastq) and (not args.mask) and args.trf:
        parser.error("\nYou cannot set the --no-generate-fastq flag without"
        " the --mask flag. Exiting...\n")

    # allow users to overwrite the defaults instead of appending to defaults
    if args.trim_args == []:
        args.trim_args = default_trimargs

    if args.bowtie2_args == []:
        args.bowtie2_args = default_bowtie2_args

    # set the default output prefix 
    if args.output_prefix == None:
        infile_base = os.path.splitext(os.path.basename(args.infile1))[0]
        args.output_prefix = infile_base + "_kneaddata"
        
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
    
    # find the bowtie2 indexes for each of the reference databases
    # reference database inputs can be directories, indexes, or index files
    reference_indexes=set()
    for directory in args.reference_db:
        reference_indexes.add(find_bowtie2_index(os.path.abspath(directory)))

    args.reference_db=list(reference_indexes)
    
    return args

def find_bowtie2_index(directory):
    """
    Search through the directory for the name of the bowtie2 index files
    Or if a file name is provided check it is a bowtie2 index
    """
    
    index=""
    # the extensions for standard bowtie2 index files
    bowtie2_index_ext_list=[".1.bt2",".2.bt2",".3.bt2",".4.bt2",
        ".rev.1.bt2",".rev.2.bt2"]
    # an extension for one of the index files for a large database
    bowtie2_large_index_ext=".1.bt2l"
    
    bowtie2_extensions=bowtie2_index_ext_list+[bowtie2_large_index_ext]
    
    if not os.path.isdir(directory):
        # check if this is the bowtie2 index file
        if os.path.isfile(directory):
            # check for the bowtie2 extension
            for ext in bowtie2_extensions:
                if re.search(ext+"$",directory):
                    index=directory.replace(ext,"")
                    break
        else:
            # check if this is the basename of the bowtie2 index files
            small_index=directory+bowtie2_index_ext_list[0]
            large_index=directory+bowtie2_large_index_ext
            if os.path.isfile(small_index) or os.path.isfile(large_index):
                index=directory
    else:
        # search through the files to find one with the bowtie2 extension
        for file in os.listdir(directory):
            # look for an extension for a standard and large bowtie2 indexed database
            for ext in [bowtie2_index_ext_list[-1],bowtie2_large_index_ext]:
                if re.search(ext+"$",file):
                    index=os.path.join(directory,file.replace(ext,""))
                    break
            if index:
                break
    
    if not index:
        sys.exit("ERROR: Unable to find bowtie2 index files in directory: " + directory)
    
    return index

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

    logging.debug("Running kneaddata with the following"
                  " arguments (from argparse): %s", str(args))

    # Get the format of the first input file
    file_format=get_file_format(args.infile1)

    if file_format != "fastq":
        logging.critical("Your input file is of type: %s . Please provide an input file of fastq format.",file_format)
        sys.exit(1)

    if args.strategy == 'memory':
        memoryheavy.memory_heavy(args)
    elif args.strategy == 'storage':
        storageheavy.storage_heavy(args)
        
    logging.info("Done!")


if __name__ == '__main__':
    main()
