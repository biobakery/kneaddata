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

from . import config

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
    

def create_directory(directory):
    """ Try to create a directory if it does not exist """
    if not os.path.exists(directory):
        logger.debug("Creating output directory: "+directory)
        try:
            os.makedirs(directory)
        except EnvironmentError:
            message="Unable to create output directory: " + directory
            logger.critical(message)
            sys.exit(message)


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
            
def format_options_to_list(input_options):
    """ Take in a list of strings with each string containing one or more options
    Format into a list of options which can be appended to a command to be run as a subprocess
    """
    formatted_options_list=[]
    for option in input_options:
        formatted_options_list+=shlex.split(option)
        
    return formatted_options_list
            
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
    
def count_reads_in_fastq_file(file,verbose):
    """ Count the number of reads in a fastq file """
    
    total_lines=0
    try:
        # file is compressed based on extension
        if file.endswith(".gz"):
            file_handle=gzip.open(file)
        else:
            file_handle=open(file)
            
        # count the lines in the file
        for line in file_handle:
            total_lines+=1
            
        file_handle.close()
    except EnvironmentError:
        total_lines=0
        message="Unable to count reads in file: "+file
        if verbose:
            print(message)
        logger.debug(message)
        
    # divide the total line number to get the total number of reads
    total_reads=total_lines/4
    
    return total_reads

def log_read_count_for_files(files,message_base,verbose=None):
    """ Log the number of reads in the files """
        
    for file in files:
        total_reads=count_reads_in_fastq_file(file,verbose)
        message=message_base+" ( "+file+" ): " + str(total_reads)
        logger.info(message)
        print(message)
        
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
        
def find_dependency(path_provided,exe,name,path_option,bypass_permissions_check):
    """ 
    Check if the dependency can be found in the path provided or in $PATH
    Return the location of the dependency
    """

    if path_provided:
        path_provided=os.path.abspath(path_provided)
        # check that the exe can be found
        try:
            files=os.listdir(path_provided)
        except EnvironmentError:
            sys.exit("ERROR: Unable to list files in "+name+" directory: "+ path_provided)
            
        if not exe in files:
            sys.exit("ERROR: The "+exe+" executable is not included in the directory: " + path_provided)
        else:
            found_path=path_provided
    else:
        # search for the exe
        exe_path=find_exe_in_path(exe, bypass_permissions_check)
        if not exe_path:
            sys.exit("ERROR: Unable to find "+name+". Please provide the "+
                "full path to "+name+" with "+path_option+".")
        else:
            found_path=exe_path  
        
    return os.path.abspath(os.path.join(found_path,exe))

def find_database_index(directory, database_type):
    """
    Search through the directory for the name of the database index files
    Or if a file name is provided check it is a database index
    For bowtie2 and bmtagger databases
    """
    
    index=""
    if database_type == "bmtagger":
        all_extensions=config.bmtagger_db_endings
    else:
        all_extensions=config.bowtie2_db_endings+[config.bowtie2_large_index_ext]
        
    # sort the extensions with the longest first, to test the most specific first
    # to find the index
    all_extensions.sort(key=lambda x: len(x), reverse=True)
    
    if not os.path.isdir(directory):
        # check if this is the database index file
        if os.path.isfile(directory):
            # check for the database extension
            for extension in all_extensions:
                if re.search(extension+"$",directory):
                    index=directory.replace(extension,"")
                    break
        else:
            # check if this is the basename of the index files
            # only need to check the first three (to include bowtie2 large index)
            for extension in all_extensions[:3]:
                if os.path.isfile(directory+extension):
                    index=directory
                    break
    else:
        # search through the files to find one with the bowtie2 extension
        for file in os.listdir(directory):
            # look for an extension for a standard and large bowtie2 indexed database
            for extension in all_extensions:
                if re.search(extension+"$",file):
                    index=os.path.join(directory,file.replace(extension,""))
                    break
            if index:
                break
    
    if not index:
        sys.exit("ERROR: Unable to find "+database_type+" index files in directory: " + directory)
    
    return index

def is_file_nonempty_readable(file, exit_on_error=None):
    """ Check that the file exists, is readable, and is not empty """
    
    error_message=""
    # check the file exists
    if os.path.exists(file):
        # check for read access
        if os.access(file, os.R_OK):
            try:
                # try to get the file size
                if os.stat(file).st_size == 0:
                    error_message="ERROR: File is empty: " + file
            except EnvironmentError:
                error_message="ERROR: Unable to check size of file: " + file
        else:
            error_message="ERROR: File is not readable: " + file
    else:
        error_message="ERROR: File does not exist: " + file
        
    
    if exit_on_error and error_message:
        sys.exit(error_message)
        
    if not error_message:
        return True
    else:
        return False
    