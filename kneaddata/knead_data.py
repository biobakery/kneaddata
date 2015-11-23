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
    group1.add_argument(
        "-i", "--input",
        help="input FASTQ file", 
        dest='infile1',
        required=True)
    group1.add_argument(
        "--input2",
        help="input FASTQ file pair",
        dest='infile2',
        default=None)
    group1.add_argument(
        "-o", "--output",
        dest='output_dir',
        help="directory to write output files",
        required=True)
    group1.add_argument(
        "-db", "--reference-db",
        default=[], action="append",
        help="location of reference database")
    group1.add_argument(
        "--output-prefix",
        help="prefix for all output files\n[ DEFAULT : $SAMPLE_kneaddata ]")
    group1.add_argument(
        "--threads",
        type=int, 
        default=config.threads,
        metavar="<" + str(config.threads) + ">",  
        help="number of threads\n[ Default : "+str(config.threads)+" ]")
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
        "-t", "--trimmomatic",
        dest='trimmomatic_path',
        help="path to trimmomatic\n[ DEFAULT : $PATH ]")
    group2.add_argument(
        "-m", "--max-mem",
        default=config.trimmomatic_memory, 
        help="max amount of memory\n[ DEFAULT : "+config.trimmomatic_memory+" ]")
    group2.add_argument(
        "-a", "--trimmomatic-options",
        default=[], action="append",
        help="options for trimmomatic\n[ DEFAULT : "+" ".join(config.trimmomatic_options)+" ]")

    group3 = parser.add_argument_group("bowtie2 arguments")
    group3.add_argument(
        "--bowtie2",
        dest='bowtie2_path',
        help="path to bowtie2\n[ DEFAULT : $PATH ]")
    group3.add_argument(
        "--bowtie2-options",
        default=[], action="append",
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
    
    # create the output directory if needed
    utilities.try_create_dir(args.output_dir)

    if (not args.no_generate_fastq) and (not args.mask) and args.trf:
        parser.error("\nYou cannot set the --no-generate-fastq flag without"
        " the --mask flag. Exiting...\n")

    # allow users to overwrite the defaults instead of appending to defaults
    if args.trimmomatic_options == []:
        args.trimmomatic_options = config.trimmomatic_options

    if args.bowtie2_options == []:
        args.bowtie2_options = config.bowtie2_options

    # set the default output prefix 
    if args.output_prefix == None:
        infile_base = os.path.splitext(os.path.basename(args.infile1))[0]
        args.output_prefix = infile_base + "_kneaddata"
        
    # find the location of trimmomatic
    trimmomatic_jar="trimmomatic-0.33.jar"
    args.trimmomatic_path=check_for_dependency(args.trimmomatic_path,trimmomatic_jar,"Trimmomatic",
        "--trim-path", True)
    # add the full path to the jar file
    args.trimmomatic_path=os.path.abspath(os.path.join(args.trimmomatic_path,trimmomatic_jar))
    
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
    if not args.log:
        args.log = os.path.join(args.output_dir,
                                    args.output_prefix+".log")

    logger = logging.getLogger()
    logger.setLevel(getattr(logging, args.log_level))
    logging.basicConfig(format=fmt)

    filehandler = logging.FileHandler(args.log)
    filehandler.setLevel(logging.DEBUG)
    filehandler.setFormatter(logging.Formatter(fmt=fmt))
    logger.addHandler(filehandler)
    
    return args

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
    # Parse the arguments from the user
    args = parse_arguments(sys.argv)
    
    # Update the configuration
    args = update_configuration(args)

    args = setup_logging(args)

    logging.debug("Running kneaddata with the following"
                  " arguments (from argparse): %s", str(args))

    # Get the format of the first input file
    file_format=utilities.get_file_format(args.infile1)

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
