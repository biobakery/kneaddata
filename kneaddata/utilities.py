import os
import sys
import shlex
import logging
import tempfile
import gzip
import re
import logging
import subprocess
import itertools
import multiprocessing
import datetime

from math import floor
from functools import partial
from contextlib import contextmanager

from . import config

# name global logging instance
logger=logging.getLogger(__name__)

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
        logger.debug("Finished running %s!" %name)
    if retcode:
        message="ERROR: " + name + " exited with code " + str(retcode)
        logger.critical(message)
        print(message)
    
    if stdout:
        logger.debug("%s stdout:\n%s", name, stdout)
    if stderr:
        logger.debug("%s stderr:\n%s", name, stderr)
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

def start_processes(commands,processes,verbose):
    """ Run the processes with the commands provided """
    
    # add verbose to command list
    commands = [i+[verbose] for i in commands]
    
    # create a pool of workers
    pool = multiprocessing.Pool(processes)
    returncodes = pool.map(run_command_returncode,commands)
    pool.close()
    pool.join()
    
    # exit if any subprocesses reported errors
    if sum(returncodes) > 0:
        sys.exit(1)
    
def run_command_returncode(args):
    """
    Convert the list of args to function arguments
    Catch errors and return code for subprocess
    """
    returncode=0
    try:
        run_command(*args, exit_on_error=False)
    except (EnvironmentError, subprocess.CalledProcessError, KeyboardInterrupt):
        returncode=1
        
    return returncode

def run_command(command,command_name,infiles,outfiles,verbose,exit_on_error):
    """ Run and log command """
    
    # convert any numbers in command to strings
    command = [str(i) for i in command]
    
    # check that the input files exist and are readable
    for file in infiles:
        logger.debug("Checking input file to "+command_name+" : "+file)
        is_file_readable(file, exit_on_error)
        
    message="Running " + command_name + " ... "
    print(message)
    logger.info(message)
    
    message=" ".join(command)
    logger.info("Execute command: " + message)
    if verbose:
        print("\n" + message + "\n")
        
    try:
        p_out = subprocess.check_output(command, stderr=subprocess.STDOUT)
        logger.debug(p_out)
    except (EnvironmentError, subprocess.CalledProcessError) as e:
        message="Error executing: " + " ".join(command) + "\n"
        if hasattr(e, 'output') and e.output:
            message+="\nError message returned from " + command_name + " :\n" + e.output
        logger.critical(message)
        log_system_status()
        if exit_on_error:
            sys.exit("CRITICAL ERROR: " + message)
        else:
            print(message)
            raise

    # check that the output files exist and are readable
    for file in outfiles:
        logger.debug("Checking output file from "+command_name+" : "+file)
        is_file_readable(file, exit_on_error) 
    
            
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
        
    # convert possible list of lists to list
    if files:
        if isinstance(files[0],list):
            files=itertools.chain.from_iterable(files)

    # count reads in each file
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

def file_size(file):
    """ Return the size of the file """
    
    try:
        size = os.stat(file).st_size
    except EnvironmentError:
        size = 0
        
    return size

def is_file_readable(file, exit_on_error=None):
    """ Check that the file exists and is readable """
    
    error_message=""
    # check the file exists
    if os.path.exists(file):
        # check for read access
        if not os.access(file, os.R_OK):
            error_message="ERROR: File is not readable: " + file
    else:
        error_message="ERROR: File does not exist: " + file
        
    if error_message:
        if exit_on_error:
            sys.exit(error_message)
        else:
            print(error_message)
            raise IOError
        
    if not error_message:
        return True
    else:
        return False

def byte_to_gigabyte(byte):
    """
    Convert byte value to gigabyte
    """
    
    return byte / (1024.0**3)
    
def log_system_status():
    """
    Print the status of the system
    """
    
    module_available=True
    try:
        import psutil
    except ImportError:
        module_available=False
        
    if module_available:
        try:
            # record the memory used
            memory = psutil.virtual_memory()
            logger.info("Total memory = " + str(byte_to_gigabyte(memory.total)) + " GB")
            logger.info("Available memory = " + str(byte_to_gigabyte(memory.available)) + " GB")
            logger.info("Free memory = " + str(byte_to_gigabyte(memory.free)) + " GB")
            logger.info("Percent memory used = " + str(memory.percent) + " %")
    
            # record the cpu info
            logger.info("CPU percent = " + str(psutil.cpu_percent()) + " %")
            logger.info("Total cores count = " + str(psutil.cpu_count()))
            
            # record the disk usage
            disk = psutil.disk_usage('/')
            logger.info("Total disk = " + str(byte_to_gigabyte(disk.total)) + " GB")
            logger.info("Used disk = "+ str(byte_to_gigabyte(disk.used)) + " GB")
            logger.info("Percent disk used = " + str(disk.percent) + " %")

            # record information about this current process
            process=psutil.Process()
            process_memory=process.memory_info()
            process_create_time=datetime.datetime.fromtimestamp(
                process.create_time()).strftime("%Y-%m-%d %H:%M:%S")
            process_cpu_times=process.cpu_times()
            # two calls required to cpu percent for non-blocking as per documentation
            process_cpu_percent=process.cpu_percent()
            process_cpu_percent=process.cpu_percent()
            
            logger.info("Process create time = " + process_create_time)
            logger.info("Process user time = " + str(process_cpu_times.user) + " seconds")
            logger.info("Process system time = " + str(process_cpu_times.system) + " seconds")
            logger.info("Process CPU percent = " + str(process_cpu_percent) + " %")
            logger.info("Process memory RSS = " + str(byte_to_gigabyte(process_memory.rss)) + " GB")
            logger.info("Process memory VMS = " + str(byte_to_gigabyte(process_memory.vms)) + " GB")
            logger.info("Process memory percent = " + str(process.memory_percent()) + " %")
            
        except (AttributeError, OSError, TypeError, psutil.Error):
            pass