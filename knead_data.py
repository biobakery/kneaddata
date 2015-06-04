#!/usr/bin/env python

'''
knead_data.py
Author: Andy Shi

Pipeline for processing metagenomics sequencing data
'''
import os
import logging
import argparse

from knead_datalib import strategies

here = os.path.abspath(os.path.dirname(__file__))


def handle_cli():
    """parse command line arguments
    note: argparse converts dashes '-' in argument prefixes to
    underscores '_'
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-1", "--infile1",
        help="input FASTQ file", 
        required=True)
    parser.add_argument(
        "-2", "--infile2",
        help="input FASTQ file mate",
        default=None)
    parser.add_argument(
        "-o", "--output-prefix",
        help="prefix for all output files")
    parser.add_argument(
        "-D", "--output-dir", default=os.getcwd(),
        help="where to put all output files")
    parser.add_argument(
        "-db", "--reference-db",
        nargs = "+", default=[],
        help=("prefix for reference databases used for either"
              " Bowtie2 or BMTagger"))
    parser.add_argument(
        "-s", "--strategy",
        default="storage",
        help=("Define operating strategy: 'storage' for IO-heavy"
              " or 'memory' for memory-heavy"))
    # Consider using a params file
    parser.add_argument(
        "-t", "--trim-path",
        required=False,
        default=os.path.join(here, "Trimmomatic-0.33/trimmomatic-0.33.jar"),
        help="path to Trimmomatic .jar executable")
    parser.add_argument(
        "--bowtie2-path",
        default=None, help="path to bowtie2 if not found on $PATH")
    #parser.add_argument("-c", "--save-contaminants-to", help="File path",
    #                    default=False, action="store_true")
    parser.add_argument(
        "--trimlen",
        type=int, default=60,
        help="minimum length for a trimmed read in Trimmomatic")
    parser.add_argument(
        "-m", "--max-mem",
        default="500m", 
        help=("Maximum amount of memory that will be used by "
              "Trimmomatic, as a string, ie 500m or 8g"))
    parser.add_argument(
        "-a", "--trim-args",
        default="SLIDINGWINDOW:4:20",
        help="additional arguments for Trimmomatic")
    # don't know how to read in a "dictionary" for the additional arguments
    parser.add_argument(
        "--bowtie2-args",
        default="",
        help="Additional arguments for Bowtie 2")
    parser.add_argument(
        "--threads",
        type=int, default=None, 
        help="Maximum number of processes to run")
    parser.add_argument(
        "--bmtagger",
        default=False, action="store_true",
        help="If set, use BMTagger to identify contaminant reads")
    parser.add_argument(
        "--extract",
        default=False, action="store_true",
        help=("Only has an effect if --bmtagger is set. If this is set,"
              " knead_data outputs cleaned FASTQs, without contaminant reads."
              " Else, output a list or lists of contaminant reads."))
    parser.add_argument(
        "--bmtagger-path",
        default=None,
        help="path to BMTagger executable if not found in $PATH")
    parser.add_argument(
        '-l', '--logging',
        default="INFO",
        help=("Logging verbosity, options are debug, info, warning,"
              " and critical. If set to debug, temporary files are not"
              " removed"))
    return parser.parse_args()


def setup_logging(loglevel, logfile):
    fmt = "%(asctime)s %(levelname)s: %(message)s"

    logger = logging.getLogger()
    logger.setLevel(getattr(logging, loglevel.upper()))
    logging.basicConfig(format=fmt)

    filehandler = logging.FileHandler(logfile)
    filehandler.setLevel(logging.DEBUG)
    filehandler.setFormatter(logging.Formatter(fmt=fmt))
    logger.addHandler(filehandler)


def main():
    args = handle_cli()

    logfile = args.output_prefix + ".log"
    setup_logging(args.logging, logfile)

    if args.strategy == 'memory':
        strategies.memory_heavy(args)
    elif args.strategy == 'storage':
        strategies.storage_heavy(args)
    else:
        logging.critical("unrecognized strategy")
        
    logging.info("Done!")


if __name__ == '__main__':
    main()
