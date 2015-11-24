import os
import sys
import shlex
import logging
import tempfile
import gzip
import re
import logging

from math import floor
from functools import partial
from contextlib import contextmanager
from multiprocessing import cpu_count

# name global logging instance
logger=logging.getLogger(__name__)

def divvy_threads(args):
    avail_cpus = args.threads or cpu_count()-1
    n_consumers = len(args.reference_db)
    trim_threads = 1
    if n_consumers > 0:
        align_threads = max(1, floor(avail_cpus/float(n_consumers)))
    else:
        align_threads = 1
    return int(trim_threads), int(align_threads)
    

def try_create_dir(d):
    if not os.path.exists(d):
        logging.warning("Directory `%s' doesn't exist. Creating.", d)
        try:
            os.makedirs(d)
        except Exception as e:
            logging.crit("Unable to create directory `%s': %s", d, str(e))
            sys.exit(1)


@contextmanager
def mktempfifo(names=("a",)):
    tmpdir = tempfile.mkdtemp()
    names = map(partial(os.path.join, tmpdir), names)
    map(os.mkfifo, names)
    try:
        yield names
    finally:
        # still perform cleanup even if there were exceptions/errors in the
        # "with" block
        map(os.remove, names)
        os.rmdir(tmpdir)


@contextmanager
def mkfifo_here(names=("a",), mode=0600):
    for n in names:
        os.mkfifo(n, mode)
    try:
        yield names
    finally:
        for n in names:
            os.remove(n)


def process_return(name, retcode, stdout, stderr):
    if name:
        logging.debug("Finished running %s!" %name)
    if retcode:
        log = logging.critical
        log("%s exited with exit status %d", name, retcode)
    else:
        log = logging.debug
    if stdout:
        log("%s stdout:\n%s", name, stdout)
    if stderr:
        log("%s stderr:\n%s", name, stderr)
    if retcode:
        sys.exit(retcode)


def parse_positive_int(string):
    try:
        val = int(string)
    except ValueError:
        raise argparse.ArgumentTypeError("Unable to parse %s to int" %string) 
    if val <= 0:
        raise argparse.ArgumentTypeError("%s is not a positive integer" %string)
    return val


def _get_bowtie2_args(bowtie2_args):
    for arg in map(shlex.split, bowtie2_args):
        for a in arg:
            yield a
            
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

def is_file_fastq(file):
    """ Return true if the file is fastq """
    
    if get_file_format(file) == "fastq":
        return True
    else:
        return False


def log_run_and_arguments(executable, arguments, verbose):
    """ Log the run and arguments and print messages """
    
    message="Running "+executable+" ..."
    print(message)
    logger.info(message)
    # log the executable and arguments
    message=executable+" " + " ".join(arguments)
    if verbose:
        print(message)
    logger.debug(message)
