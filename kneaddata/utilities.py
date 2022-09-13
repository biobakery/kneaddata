"""
KneadData: utilities module

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

import os
import fnmatch
import sys
import shlex
import logging
import tempfile
import gzip
import bz2
import re
import subprocess
import itertools
import multiprocessing
import datetime
import errno
import shutil
from zipfile import ZipFile

from kneaddata import config

# name global logging instance
logger=logging.getLogger(__name__)

def get_read_id_minus_pair(sequence_id_line):
    return sequence_id_line.rstrip()[:-1]

def update_temp_output_files(temp_output_files, new_file_list, all_input_files):
    """ Remove files on the list that are no longer needed and add the new file """

    for filename in temp_output_files:
        if not filename in all_input_files:
            remove_file(filename)
            temp_output_files.remove(filename)

    if not isinstance(new_file_list, list):
        new_file_list=[new_file_list]

    temp_output_files+=new_file_list

def check_and_reorder_reads(input_files, output_folder, temp_output_files):
    """ Check if reads are ordered and if not reorder """

    # read in the ids from the first pair (only check the first 100)
    ids = []
    for count, lines in zip(range(100),read_file_n_lines(input_files[0],4)):
        ids.append(get_read_id_minus_pair(lines[0]))

    mismatch=False
    for lines, pair_id in zip(read_file_n_lines(input_files[1],4), ids):
        if not get_read_id_minus_pair(lines[0]) == pair_id:
            mismatch=True
            break

    # reorder the pairs to match
    new_file_list = []
    if mismatch:
        message="Reordering read identifiers ..."
        print(message+"\n")
        logger.info(message)

        for index, infile in enumerate(input_files):
            file_out, new_file=tempfile.mkstemp(prefix="reordered_",
                suffix="_"+file_without_extension(infile), dir=output_folder)
            os.close(file_out)

            # read in all of the sequences then sort and write out
            ids={}
            for lines in read_file_n_lines(infile,4):
                id=get_read_id_minus_pair(lines[0])
                ids[id]=lines

            with open(new_file,"w") as file_handle:
                for id in sorted(ids.keys()):
                    file_handle.write("".join(ids[id]))

            # set the input file to the reordered temp file
            input_files[index]=new_file
            new_file_list.append(new_file)

        # add the temp file to the list and remove extra that are not needed
        update_temp_output_files(temp_output_files, new_file_list, input_files)

    return input_files

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
        print("Subprocess reported error. Please see log file for more details.")
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

def run_command(command,command_name,infiles,outfiles,stdout_file,verbose,exit_on_error,shell=False):
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
       
    if stdout_file:
        try:
            stdout=open(stdout_file,"w")
        except EnvironmentError:
            message="Unable to open file: " + stdout_file
            logger.critical(message)
            if exit_on_error:
                sys.exit("CRITICAL ERROR: " + message)
            else:
                raise EnvironmentError
    try:
        if stdout_file:
            if shell:
                p_out = subprocess.check_call(" ".join(command), stdout=stdout, shell=True)
            else:
                # run command, raise CalledProcessError if return code is non-zero
                p_out = subprocess.check_call(command, stdout=stdout)
        else:
            p_out = subprocess.check_output(command, stderr=subprocess.STDOUT)
        logger.debug(p_out)
    except (EnvironmentError, subprocess.CalledProcessError) as e:
        message="Error executing: " + " ".join(command) + "\n"
        if hasattr(e, 'output') and e.output:
            message+="\nError message returned from " + command_name + " :\n" + e.output.decode("utf-8")
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

def bunzip2_file(bz2_file, new_file):
    """ Return a new copy of the decompressed file """

    
    try:
        decompress_function=bz2.open
    except AttributeError:
        decompress_function=bz2.BZ2File

    message="Decompressing bzipped2 file ..."
    print(message+"\n")
    logger.info(message)    

    with open(new_file,"w") as file_write:
        with decompress_function(bz2_file, "rt") as file_read:
            for line in file_read:
                file_write.write(line)

    logger.info("Decompressed file created: " + new_file)

    return new_file

def gunzip_file(gzip_file, new_file):
    """
    Return a new copy of the file that is not gzipped
    """
    
    message="Decompressing gzipped file ..."
    print(message+"\n")
    logger.info(message)    
    
    try:
        file_handle_gzip=gzip.open(gzip_file,"rt")
        
        # write the gunzipped file
        file_handle=open(new_file,"w")
        shutil.copyfileobj(file_handle_gzip, file_handle)
        
    except EnvironmentError:
        sys.exit("Critical Error: Unable to gunzip input file: " + gzip_file)
    finally:
        file_handle.close()
        file_handle_gzip.close()
        
    logger.info("Decompressed file created: " + new_file)
        
    return new_file

def file_without_extension(file):
    """ Return the basename of the file without the extension """
    
    return os.path.splitext(os.path.basename(file))[0]

def get_decompressed_file(file, output_folder, temp_file_list, all_input_files):
    """ Check if a file is compressed, if so decompress """
    
    if file.endswith(".gz") or file.endswith(".bz2"):
        file_out, new_file=tempfile.mkstemp(prefix="decompressed_", 
            suffix="_"+file_without_extension(file), dir=output_folder)
        os.close(file_out)
        if file.endswith(".gz"):
            gunzip_file(file,new_file)
        else:
            bunzip2_file(file,new_file)
        update_temp_output_files(temp_file_list, new_file, all_input_files)
    else:
        new_file=file
        
    return new_file

def check_sequence_identifier_format(file):
    """ Check the fastq file to see if there are spaces in the identifier
        and the format of the id to see if this is the new illumina format """ 
    
    #checking first 100 (400/4) lines
    num_seq_to_check=100 
    num_lines_to_check=400
    
    #Fetching first and last 100 identifier sequences
    first_seq_identifiers_list=get_first_n_seq_identifiers(file,num_seq_to_check)
    last_seq_identifiers_list=get_last_n_seq_identifiers(file,num_lines_to_check)
    # Checking first and last 100 seq identifiers for spaces and new Illumina format
    for lines in first_seq_identifiers_list:
        new_format=sequence_identifier_format_conditions(lines)
    for lines in last_seq_identifiers_list:
        new_format=sequence_identifier_format_conditions(lines)
    return new_format
  
def sequence_identifier_format_conditions(identifier_seq):
    new_format=False
    if (" " in identifier_seq):
        new_format=True
    if not identifier_seq.endswith("/1\n") and not identifier_seq.endswith("/2\n"):
        new_format=True
    return new_format
    
def get_last_n_seq_identifiers(file, n):
    last_seq_identifiers=[]
    # Tail to find last lines
    try:
        process = subprocess.Popen(['tail', '-'+str(n), file], stdout=subprocess.PIPE)
    except subprocess.CalledProcessError:
        pass
    for i,line in enumerate(process.stdout.readlines()):
        if (i%4==0):
            last_seq_identifiers.append(line.decode("utf-8")) 
    return last_seq_identifiers
    
def get_first_n_seq_identifiers(file,n):
    count=0
    all_lines=read_file_n_lines(file,4)
    first_seq_identifiers=[]
    # Getting first nth seq identifier
    while(count<n):
        try:
            lines=next(all_lines)
        except StopIteration:
            # allow for less then n lines
            break
        first_seq_identifiers.append(lines[0])
        count+=1
    return first_seq_identifiers


        
def get_reformatted_identifiers(file, input_index, output_folder, temp_file_list, all_input_files):
    """ Reformat the sequence identifiers in the fastq file writing to a temp file """
    
    # check if the file needs to be reformatted
    reformat_file=check_sequence_identifier_format(file)
    
    if not reformat_file:
        return file
    
    message="Reformatting file sequence identifiers ..."
    print(message+"\n")
    logger.info(message)   
    
    file_out, new_file=tempfile.mkstemp(prefix="reformatted_identifiers",
        suffix="_"+file_without_extension(file), dir=output_folder)
    os.close(file_out)
    
    with open(new_file, "wt") as file_handle:
        for lines in read_file_n_lines(file,4):
            # reformat the identifier and write to temp file
            if " " in lines[0]:
                lines[0]=lines[0].replace(" ",".")
            if not lines[0].endswith("/1\n") and not lines[0].endswith("/2\n"):  
                if (input_index == 0):
                    lines[0]=lines[0].rstrip()+"#0/1\n"
                else:
                    lines[0]=lines[0].rstrip()+"#0/2\n"
                    
                    
            file_handle.write("".join(lines))
    
    # add the new file to the list of temp files
    update_temp_output_files(temp_file_list, new_file, all_input_files)
    
    return new_file

def bam_to_sam(bam_file, new_file):
    """
    Convert from a bam to sam file
    """

    exe=config.samtools_exe
    args=["view","-h",bam_file]
    
    message="Converting bam file to sam format ..."
    print(message)
    logger.info(message)
    
    run_command([exe]+args, exe, [bam_file], [new_file], new_file, verbose=False, exit_on_error=True)
    
    logger.info("Sam file created at: " + new_file)
    
def bam_to_fastq(bam_file, new_file, output_folder):
    """ Convert bam file to single or set of fastq files """

    exe=config.samtools_exe
    args=["bam2fq",bam_file]

    message="Converting bam file to fastq format ..."
    print(message)
    logger.info(message)
    
    run_command([exe]+args, exe, [bam_file], [new_file], new_file, verbose=False,exit_on_error=True)

    # check for paired-end output
    try:
        output=subprocess.check_output("grep '^@.*/2$' "+new_file,shell=True)
    except (EnvironmentError, subprocess.CalledProcessError):
        output=[]

    if output:
        # pairs found
        filename=file_without_extension(bam_file)
        pair1_file=os.path.join(output_folder,file_without_extension(filename)+"_decompressed_R1"+".fastq")
        pair2_file=os.path.join(output_folder,file_without_extension(filename)+"_decompressed_R2"+".fastq")

        run_command(["grep","'^@.*/1$'",new_file,"-A 3","--no-group-separator"],"bam to fastq R1",[new_file],[pair1_file],pair1_file,verbose=False,exit_on_error=True,shell=True)
        run_command(["grep","'^@.*/2$'",new_file,"-A 3","--no-group-separator"],"bam to fastq R2",[new_file],[pair2_file],pair2_file,verbose=False,exit_on_error=True,shell=True)

        remove_file(new_file)

        new_file=[pair1_file,pair2_file]

    return new_file


def get_fastq_from_bam_file(file, output_folder, temp_file_list, all_input_files):
    """ Look for paired input files from bam, requires samtools """
    if file.endswith(".bam"):
        # Check for the samtools software
        if not find_exe_in_path(config.samtools_exe):
            sys.exit("CRITICAL ERROR: The samtools executable can not be found. "
            "Please check the install or select another input format.")
        new_file=os.path.join(output_folder,file_without_extension(file)+"_decompressed"+".fastq")
        new_file=bam_to_fastq(file, new_file, output_folder) 
        update_temp_output_files(temp_file_list, new_file, all_input_files)
    else:
        new_file=file

    return new_file


def get_sam_from_bam_file(file, output_folder, temp_file_list, all_input_files):
    """ Check if a file is bam, if so create a sam file """
    
    if file.endswith(".bam"):
        # Check for the samtools software
        if not find_exe_in_path(config.samtools_exe):
            sys.exit("CRITICAL ERROR: The samtools executable can not be found. "
            "Please check the install or select another input format.")
        new_file=os.path.join(output_folder,file_without_extension(file)+"_decompressed"+".sam")
        bam_to_sam(file,new_file)
        update_temp_output_files(temp_file_list, new_file, all_input_files)
    else:
        new_file=file
        
    return new_file

def sam_to_fastq(file, new_file):
    """ Create a fastq file from the sam file """
    
    message="Converting sam file to fastq format ..."
    print(message+"\n")
    logger.info(message)    
    
    # First read in all of the ids
    # Read then write to not write duplicates
    # Also do not store all sequences and quality scores here to save space
    ids=set()
    with open(file) as file_handle:
        for line in file_handle:
            if not re.search("^@",line):
                data=line.rstrip().split(config.sam_delimiter)
                ids.add(data[config.sam_read_name_index])
        
    # Now read through the file again, writing out fastq sequences
    with open(file) as file_handle:
        with open(new_file, "w") as write_file_handle:
            for line in file_handle:
                if not re.search("^@",line):
                    data=line.rstrip().split(config.sam_delimiter)
                    read_id=data[config.sam_read_name_index]
                    if read_id in ids:
                        write_file_handle.write("@"+read_id+"\n")
                        write_file_handle.write(data[config.sam_read_index]+"\n")
                        write_file_handle.write("+\n")
                        write_file_handle.write(data[config.sam_read_quality]+"\n")
                        # remove the id so as to not write it more than once
                        ids.remove(read_id)

def get_fastq_from_sam_file(file, output_folder, temp_file_list, all_input_files):
    """ Check if a file is sam, if so create a fastq file """
    
    if file.endswith(".sam"):
        new_basename=file_without_extension(file)
        if not "decompressed" in new_basename:
            new_basename+="_decompressed"
        new_file=os.path.join(output_folder,new_basename+config.fastq_file_extension)
        sam_to_fastq(file,new_file)
        update_temp_output_files(temp_file_list, new_file, all_input_files)
    else:
        new_file=file
        
    return new_file
            
def get_file_format(file):
    """ Determine the format of the file """

    format="unknown"
    file_handle=None

    # check the file exists and is readable
    if not os.path.isfile(file):
        logger.critical("The input file selected is not a file: %s.",file)

    if not os.access(file, os.R_OK):
        logger.critical("The input file selected is not readable: %s.",file)

    try:
        # check for gzipped files
        if file.endswith(".gz"):
            file_handle = gzip.open(file, "rt")
        else:
            file_handle = open(file, "rt")

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

def resolve_sublists(lists):
    """ Resolve the sublists in list """
    
    if lists:
        if isinstance(lists[0], list):
            lists=list(itertools.chain.from_iterable(lists))
            
    return lists

def get_file_types(files,type,database_names):
    """ Get the types for each of the files based on the names and database"""
    
    # get the file types for raw files
    file_types=[]
    if type == "raw" and len(files) == 2:
        file_types=["pair1","pair2"]
    elif type == "trimmed":
        for file in files:
            for ending, name in config.trimmomatic_pe_names.items():
                if file.endswith(ending):
                    file_types.append(name)
                    break
    elif type == "decontaminated":
        for file, db_name in zip(files,database_names):
            basename=os.path.basename(file)
            if file.endswith("clean_1.fastq"):
                file_types.append(db_name+" pair1")
            elif file.endswith("clean_2.fastq"):
                file_types.append(db_name+" pair2")
            elif "_unmatched_1" in basename:
                file_types.append(db_name+" orphan1")
            elif "_unmatched_2" in basename:
                file_types.append(db_name+" orphan2")
            else:
                file_types.append(db_name+" single")
    elif type == "final":
        for file in files:
            for ending, name in config.final_file_types.items():
                if file.endswith(ending):
                    file_types.append(name)
                    break
    
    if not file_types:
        file_types=["single"]*len(files)  
        
    return file_types

def log_read_count_for_files(files,type,message_base,database_names=None,verbose=None):
    """ Log the number of reads in the files """
        
    # convert possible list of lists to list
    files=resolve_sublists(files)
    
    # get the types for the files for the log message
    file_types=get_file_types(files,type,database_names)

    # count reads in each file
    for file, file_type in zip(files,file_types):
        total_reads=count_reads_in_fastq_file(file,verbose)
        message=message_base+" ( "+file+" ): " + str(total_reads)
        logger.info("READ COUNT: "+type+" "+file_type+" : "+message)
        print(message)
        
def find_exe_in_path(exe, bypass_permissions_check=None, add_exe_to_path=None):
    """
    Check that an executable exists in $PATH
    """
    
    paths = os.environ["PATH"].split(os.pathsep)
    for path in paths:
        fullexe = os.path.join(path,exe)
        if os.path.exists(fullexe):
            if not bypass_permissions_check:
                check_file_executable(fullexe)
            if add_exe_to_path:
                path=fullexe
            return path
        elif os.path.isdir(path):
            # allow for filename filter matching
            exematch = fnmatch.filter(os.listdir(path),exe)
            if exematch and os.path.exists(os.path.join(path,exematch[0])):
                if not bypass_permissions_check:
                    check_file_executable(os.path.join(path,exematch[0]))
                if add_exe_to_path:
                    path=os.path.join(path,exematch[0])
                return path

    return None
        
def add_exe_to_path(exe_dir):
    """ 
    Add path to executable to $PATH
    """
    
    logger.debug("Add directory, %s, to path", exe_dir)
    
    os.environ["PATH"] = exe_dir + os.pathsep + os.environ["PATH"]        
    
def check_file_executable(exe):
    """
    Check the file can be executed
    """
    
    try:
        output=subprocess.check_output([exe,"--version"],stderr=subprocess.STDOUT)
    except EnvironmentError as error:
        if error.errno == errno.EACCES:
            sys.exit("ERROR: Unable to execute software: " + exe)
    except subprocess.CalledProcessError:
        pass

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
            
        found_paths=fnmatch.filter(files, exe)
        if not found_paths:
            sys.exit("ERROR: The "+exe+" executable is not included in the directory: " + path_provided)
        else:
            exe=found_paths[0]
            found_path=path_provided
            # check permissions
            if not bypass_permissions_check:
                check_file_executable(os.path.abspath(os.path.join(found_path,exe)))

        dependency_path=os.path.abspath(os.path.join(found_path,exe))

    else:
        # search for the exe
        dependency_path=find_exe_in_path(exe, bypass_permissions_check, add_exe_to_path=True)
        if not dependency_path:
            sys.exit("ERROR: Unable to find "+name+". Please provide the "+
                "full path to "+name+" with "+path_option+".")
        
    return dependency_path

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
   
def move_file(file, old_folder, new_folder):
    """ Try to move a file """

    old_file=os.path.join(old_folder,file)
    new_file=os.path.join(new_folder,file)

    remove_file(new_file)

    try:
        shutil.move(old_file,new_file)
    except EnvironmentError as e:
        logger.warning("Unable to move file "+old_file+" to "+new_file)
        logger.warning(e)
 
def remove_file(file):
    """ Try to remove the file """
    
    try:
        os.unlink(file)
    except EnvironmentError:
        logger.warning("Unable to remove file: " + file)

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
        
        
def read_file_n_lines(file,n):
    """ Read a file n lines at a time """
    
    line_set=[]
    with open(file) as file_handle:
        for line in file_handle:
            if len(line_set) == n:
                yield line_set
                line_set=[]
            line_set.append(line)
    
    # yield the last set
    if len(line_set) == n:
        yield line_set
        
def get_read_length_fastq(file):
    """ Get the read length from a fastq file """
    
    try:
        file_handle=open(file)
    except EnvironmentError:
        sys.exit("Unable to read file: " + file)
        
    sequence_id=file_handle.readline()
    sequence=file_handle.readline()
    
    file_handle.close()
    
    return len(sequence.rstrip())

def get_default_trimmomatic_options(read_length=None, path="", type="SE", sequencer_source=""):
    """ Get the default trimmomatic options """
    
    if not read_length:
        read_length=config.default_read_length
        
    # have the minlen equal to a percent of the read length
    minlen=int(read_length*(config.trimmomatic_min_len_percent/100.0))
    
    if (sequencer_source.lower()=="none"):
        return [config.trimmomatic_slidingwindow_option,
                config.trimmomatic_minlen_option_tag+config.trimmomatic_option_delimiter+str(minlen)]
    else:
        if type == "PE":
            adapter_settings = config.trimmomatic_trim_adapters_option_pe.replace("$PATH",path)
            if sequencer_source!=config.trimmomatic_provided_sequencer_default:
                adapter_settings = adapter_settings.replace(config.trimmomatic_provided_sequencer_default,sequencer_source)
        else:
            adapter_settings = config.trimmomatic_trim_adapters_option_se.replace("$PATH",path)
            if sequencer_source!=config.trimmomatic_provided_sequencer_default:
                adapter_settings = adapter_settings.replace(config.trimmomatic_provided_sequencer_default,sequencer_source)
                adapter_settings = adapter_settings.replace('PE','SE')
            
        return [config.trimmomatic_minlen_option_tag+config.trimmomatic_option_delimiter+str(config.min_init_read_length),
                adapter_settings,
                config.trimmomatic_slidingwindow_option,
                config.trimmomatic_minlen_option_tag+config.trimmomatic_option_delimiter+str(minlen)]
    
def cat_files(files,output_file):
    """ Cat the files to a single file """
    
    # check that the files exist
    file_list=list(filter(os.path.isfile,files))
    
    try:
        stdout=open(output_file,"w")
    except EnvironmentError:
        sys.exit("ERROR: Unable to open file: " + output_file)
    
    try:
        subprocess.check_call(["cat"]+file_list,stdout=stdout)
    except (subprocess.CalledProcessError,EnvironmentError):
        sys.exit("ERROR: Unable to cat files.")

def fastq_to_fasta(file, new_file):
    """
    Convert fastq file to fasta
    
    Fastq format short example:
    @SEQ_ID
    GATCTGG
    +
    !''****
    
    Fasta format short example:
    >SEQ_INFO
    GATCTGG
    """
    
    try:
        file_handle_read = open(file, "rt")
        line = file_handle_read.readline()
    except EnvironmentError:
        sys.exit("ERROR: Unable to read file: " + file)
    
    try:
        file_out=open(new_file,"w")
    except EnvironmentError:
        sys.exit("ERROR: Unable to write file: " + file)

    sequence=""
    sequence_id=""
    while line:
        if re.search("^@",line):
            # write previous sequence
            if sequence:
                file_out.write(sequence_id)
                file_out.write(sequence+"\n")
            
            sequence_id=line.replace("@",">",1)
            sequence=""
        elif re.search("^[A|a|T|t|G|g|C|c|N|n]+$", line):
            sequence+=line.rstrip()
        line=file_handle_read.readline()
        
    # write out the last sequence
    if sequence:
        file_out.write(sequence_id)
        file_out.write(sequence+"\n")    

    file_out.close()    
    file_handle_read.close()

    return new_file    

def write_read_count_table(output, reads):
    """
    Write the table of read counts for all samples
    """

    with open(output, "w") as file_handle:
        # get the headers from all of the samples
        headers=set()
        for sample, counts in reads.items():
            headers.update(counts.keys())
        # order the headers
        header_order=[]
        for column in ["raw","trimmed","decontaminated","final"]:
            order=sorted(list(filter(lambda x: x.startswith(column),headers)))
            
            # sort by paired then orphan, if present
            paired=list(filter(lambda x: "pair" in x, order))
            if paired:
                header_order+=paired
                orphan=list(filter(lambda x: "orphan" in x, order))
                header_order+=orphan
            else:
                header_order+=order
        
        file_handle.write("\t".join(["Sample"]+header_order)+"\n")
        for sample in sorted(reads.keys()):
            new_line=[sample]
            for column in header_order:
                try:
                    counts=reads[sample][column]
                except KeyError:
                    counts="NA"
                new_line.append(counts)
            file_handle.write("\t".join([str(i) for i in new_line])+"\n")

def get_updated_trimmomatic_parameters(input_list,output_dir,default_trimmomatic_options):
    adapter_dir_path=output_dir+"/adapters.fa"
    overreq_seq_length=0
    overreq_seq_length_list=[]
    seq_list=[]
    counter=0 
    fout = open(adapter_dir_path, "w")
    for input_file in input_list:
        try:
            f = open(input_file,"r")
            line = f.readline()
            while line: 
                if ">>Overrepresented sequences" in line: 
                    overrepresented_seq_status=line
                    #Logging Overrepresented sequences status
                    message = "\n>>Overrepresented sequences in "+input_file+"\n"
                    logger.info(message)
                    overreq_seq_list=[]
                    while not (">>END_MODULE" in line):
                        line = f.readline()
                        logger.info(line)
                        overreq_seq_list.append(line)
                                        
                    seq_list= overreq_seq_list[1:-1]
                    check_flag=False
                    for seq in seq_list: 
                        fout.write (">customAdapter"+str(counter)+"\n")
                        fout.write  (seq.split('\t')[0]+"\n")
                        check_flag=True
                        #Calculating length of the overrepresentted sequence
                        overreq_seq_length_list.append(len(seq.split('\t')[0]))
                        counter+=1
                    if not check_flag:
                        message = "\n>>No overrepresented sequences found in "+input_file+" Bypassing filtering for these sequences.\n"
                        logger.info(message)
                            
                #Logging Adapter Content status
                if ">>Adapter Content" in line:
                    message=line
                    logger.info(message)
                line = f.readline()
            f.close()
        except IOError:
            message = "Could not read FASTQC generated file: "+input_file
            logger.info(message)

    if overreq_seq_length_list:
        #Calculating the value of trimmomatic option based on overrepresented sequences
        for trimmomatic_options_count,trimmomatic_option in enumerate(default_trimmomatic_options):
            if trimmomatic_option.count("ILLUMINACLIP")>0: 
                sequence_adapter_path = trimmomatic_option.split(':')[1]
                trimmomatic_parameter = trimmomatic_option.split('.')
                with open(sequence_adapter_path) as sequence_adapter_file:
                    for loop_index,line in enumerate(sequence_adapter_file):
                        fout.write(line)
                        if loop_index%2!=0:
                            overreq_seq_length_list.append(len(line))
                            
                # Multiplying Minimun Overrepresented Seq Length with log function(0.6)
                min_overreq_seq_length =str(int(min(overreq_seq_length_list)*0.6))
                splited_trimmomatic_parameter = trimmomatic_parameter[-1].split(':')
                #Replace trimmomatic parameters with minimum overrepresented sequence length
                if len(input_list) == 2:
                    splited_trimmomatic_parameter[2] = min_overreq_seq_length
                else:
                    splited_trimmomatic_parameter[3] = min_overreq_seq_length
                temp_updated_parameter=':'.join(splited_trimmomatic_parameter)
                updated_parameter = ':'.join(temp_updated_parameter.split(':')[1:])
                #Updating the Global trimmomation_options value
                default_trimmomatic_options[trimmomatic_options_count]="ILLUMINACLIP:"+adapter_dir_path+":"+updated_parameter
    return default_trimmomatic_options
