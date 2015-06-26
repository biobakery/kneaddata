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

from knead_datalib import strategies
from knead_datalib import try_create_dir

here = os.path.abspath(os.path.dirname(__file__))


def parse_positive_int(string):
    try:
        val = int(string)
    except ValueError:
        raise argparse.ArgumentTypeError("Unable to parse %s to int" %string) 
    if val <= 0:
        raise argparse.ArgumentTypeError("%s is not a positive integer" %string)
    return val


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

    group2 = parser.add_argument_group("trimmomatic arguments")
    group2.add_argument(
        "-t", "--trim-path",
        default=os.path.join(here, "Trimmomatic-0.33/trimmomatic-0.33.jar"),
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
        default=default_trimargs, action="append",
        help=("additional arguments for Trimmomatic, default: "
              +" ".join(default_trimargs)))

    group3 = parser.add_argument_group("bowtie2 arguments")
    group3.add_argument(
        "--bowtie2-path",
        default=None, help="path to bowtie2 if not found on $PATH")
    group3.add_argument(
        "--bowtie2-args",
        default=[], action="append",
        help="Additional arguments for Bowtie 2")
        
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
            help="TRF matching weight")
    group5.add_argument(
            "--mismatch", type=parse_positive_int,
            default=7, 
            help="TRF mismatching penalty")
    group5.add_argument(
            "--delta", type=parse_positive_int,
            default=7, 
            help="TRF indel penalty")
    group5.add_argument(
            "--pm", type=parse_positive_int,
            default=80,
            help="TRF match probability (whole number)")
    group5.add_argument(
            "--pi", type=parse_positive_int,
            default=10,
            help="TRF indel probability (whole number)")
    group5.add_argument(
            "--minscore", type=parse_positive_int,
            default=50, 
            help="TRF minimum alignment score to report")
    group5.add_argument(
            "--maxperiod", type=parse_positive_int,
            default=500, 
            help="TRF maximum period size to report")
    group5.add_argument(
            "--no-generate-fastq",
            default=True, action="store_false", 
            help="If switched on, don't generate fastq")
    group5.add_argument(
            "--mask",
            default=False, action="store_true",
            help="Generate mask file")
    group5.add_argument(
            "--html",
            default=False, action="store_true",
            help="Generate html file")

    args = parser.parse_args()
    try_create_dir(args.output_dir)

    if (not args.no_generate_fastq) and (not args.mask) and args.trf:
        parser.error("\nYou cannot set the --no-generate-fastq flag without"
        " the --mask flag. Exiting...\n")

    return args


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
