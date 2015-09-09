import os
import re
import sys
import time
import shutil
import logging
import itertools
import subprocess
import collections
from functools import partial

from . import constants_knead_data as const
from . import divvy_threads, try_create_dir, mktempfifo, process_return, _get_bowtie2_args
import memoryheavy

def _generate_bowtie2_commands( infile_list, db_prefix_list,
                                bowtie2_path, output_prefix,
                                bowtie2_opts, tmp_dir, threads ):
    db_prefix_bases = _prefix_bases(db_prefix_list)
    for basename, fullpath in db_prefix_bases:
        is_paired = (len(infile_list) == 2)
        cmd = [bowtie2_path] 
        #cmd += list(dict_to_cmd_opts_iter(bowtie2_opts))
        cmd += [ "-x", fullpath ]
        cmd += [ "--threads", str(threads) ]

        output_str = output_prefix + "_" + basename
        outputs_to_combine = list()
        if is_paired:
            cmd += [ "-1", infile_list[0],
                     "-2", infile_list[1],
                     "--al-conc", output_str + "_contam_%.fastq",
                     "--un-conc", output_str + "_clean_%.fastq"]
            outputs_to_combine = [output_str + "_clean_1.fastq", 
                                  output_str + "_clean_2.fastq"]
        else:
            cmd += [ "-U", infile_list[0],
                     "--al", output_str + "_contam.fastq",
                     "--un", output_str + "_clean.fastq"]
            outputs_to_combine = [output_str + "_clean.fastq"]

        cmd += list(_get_bowtie2_args(bowtie2_opts))
        sam_out = os.path.join(tmp_dir, os.path.basename(output_str) + ".sam")
        cmd += [ "-S", sam_out ]
        yield (cmd, outputs_to_combine)


def _poll_workers(popen_list):
    """ Polls a list of processes initialized by subprocess.Popen. Returns the
    number of processes still running and a list of those processes. If one or
    more of the processes returned with non-zero exit code, raises an error. 

    popen_list: a list of objects from subprocess.Popen
    """
    failures = list()
    still_running = 0
    polls = [ (proc.poll(), proc, cmd) for proc, cmd in popen_list ]
    procs_still_running = list()
    for val, p, cmd in polls:
        if val is None:
            still_running += 1
            procs_still_running.append((p, cmd))
        elif val != 0:
            failures.append((val, cmd))
            
    if failures:
        msg = [ " ".join(cmd) + " Returned code " + val 
                for cmd, val in failures ]
        raise OSError("The following commands failed: "+msg)
    else:
        return (still_running, procs_still_running)


def align(infile_list, db_prefix_list, output_prefix, tmp_dir,
          bowtie2_path=None, n_procs=None, bowtie2_opts=list()):

    """Align a single-end sequence file or a paired-end duo of sequence
    files using bowtie2.

    :param infile_list: List; Input sequence files in fastq format. A
                        list of length 1 is interpreted as a single-end
                        sequence file, while a list of length 2 is 
                        interpreted as a paired-end sequence duo.
    :param db_prefix_list: List; Bowtie2 database prefixes. Multiple database
                           can be specified. Uses subprocesses to run bowtie2 
                           in parallel.
    :param output_prefix: String; Prefix of the output file
    :param tmp_dir: String; Path name to the temporary for holding bowtie2 
                    sam file output
    :keyword bowtie2_path: String; File path to the bowtie2 binary's location.
    :keyword n_procs: Int; Number of bowtie2 subprocesses to run.
    :keyword bowtie2_opts: List; List of additional arguments, as strings, to 
                            be passed to Bowtie2.
    """

    if not bowtie2_path:
        bowtie2_path = find_on_path("bowtie2")
        if not bowtie2_path:
            raise Exception("Could not find Bowtie2 path")

    if not n_procs:
        n_procs = len(db_prefix_list)

    commands_to_run = _generate_bowtie2_commands( infile_list,
                                                  db_prefix_list,
                                                  bowtie2_path,
                                                  output_prefix,
                                                  bowtie2_opts,
                                                  tmp_dir,
                                                  n_procs )
    commands_to_run = list(commands_to_run)
    
    procs_running = list()
    outputs = list()

    # poll to see if we can run more Bowtie2 instances
    while commands_to_run:
        cmd, output_to_combine = commands_to_run.pop()
        n_running, procs_running = _poll_workers(procs_running)
        if n_running >= n_procs:
            commands_to_run.append((cmd, output_to_combine))
            time.sleep(0.5)
        else:
            logging.debug("Running bowtie2 command: " + " ".join(cmd))
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            procs_running.append((proc, cmd))
            outputs.append(output_to_combine)

    # wait for everything to finish after we started running all the processes
    for p, cmd in procs_running:
        stdout, stderr = p.communicate()
        cmd = " ".join(cmd)
        logging.debug("bowtie2 command `%s' Complete", cmd)
        if stdout:
            logging.debug("bowtie2 command `%s' stdout: %s", cmd, stdout)
        if stderr:
            logging.debug("bowtie2 command `%s' stderr: %s", cmd, stderr)

    ret_codes = [p.returncode for p, _ in procs_running]

    # if Bowtie2 produced correct output, merge the files from multiple
    # databases
    combined_outs = []
    if (outputs != []):
        combined_outs = combine_tag(outputs, output_prefix)

    return(ret_codes, commands_to_run, combined_outs)


def tag(infile_list, db_prefix_list, temp_dir, prefix,
        bmtagger_path=None, n_procs=2, remove=False, debug=False):
    """
    Runs BMTagger on a single-end sequence file or a paired-end duo of sequence
    files. Returns a tuple (ret_codes, bmt_args). Both are lists, each which has
    the same number of elements as the number of databases. THe first contains
    the return codes for each of the BMTagger executions, and the second
    contains the command line calls for these executions.

    :param infile_list: List; a list of input files in FASTQ format. A list of
                        length 1 is interpreted as a single-end sequence file,
                        while a list of length 2 is interpreted as a paired-end
                        sequence duo
    :param db_prefix_list: List; Prefixes of the BMTagger databases. BMTagger
                           needs databases of the format db_prefix.bitmask, 
                           db_prefix.srprism.*, and db_prefix.{blastdb file
                           extensions} for any fixed db_prefix (see
                           const.BMTAGGER_DB_ENDINGS for full list of endings).
                           Multiple databases can be specified. Uses
                           subprocesses to run BMTagger in parallel
    :param temp_dir: String; Path name to temporary directory for BMTagger
    :param prefix: String; prefix for output files
    :keyword bmtagger_path: String; file path to the bmtagger.sh executable
    :keyword n_procs: Int; Max number of BMTagger subprocesses to run
    :keyword remove: Boolean; If True, output cleaned FASTQ files with reads
                     identified as contaminants removed. Otherwise, outputs
                     list(s) of reads identified as contaminants.
    :keyword debug: Boolean; If True, BMTagger does not delete temporary files
    """
    single_end = (len(infile_list) == 1)

    if not bmtagger_path:
        bmtagger_path = find_on_path("bmtagger.sh")
        if not bmtagger_path:
            raise Exception("Could not find BMTagger path!")

    logging.debug(db_prefix_list)
    db_list = list(_prefix_bases(db_prefix_list))
    bmt_args = [None for d in db_list]

    outputs = [None for d in db_list]
    # build arguments
    for (i, (basename, fullpath)) in enumerate(db_list):
        out_prefix = None
        if single_end:
            if remove:
                out_prefix = prefix + "_" + basename + "_clean"
                bmt_args[i] = [bmtagger_path, 
                              "-q", "1",
                              "-1", infile_list[0],
                              "-b", str(fullpath + ".bitmask"),
                              "-x", str(fullpath + ".srprism"),
                              "-T", temp_dir, 
                              "-o", out_prefix,
                              "--extract"] 
                outputs[i] = [out_prefix + const.BMTAGGER_SE_ENDING]

            else:
                out_prefix = prefix + "_" + basename + "_contam.out"
                bmt_args[i] = [bmtagger_path, 
                              "-q", "1",
                              "-1", infile_list[0],
                              "-b", str(fullpath + ".bitmask"),
                              "-x", str(fullpath + ".srprism"),
                              "-T", temp_dir, 
                              "-o", out_prefix]
                outputs[i] = [out_prefix]


        else:
            if remove:
                out_prefix = prefix + "_" + basename + "_clean"
                bmt_args[i] = [bmtagger_path, 
                              "-q", "1",
                              "-1", infile_list[0],
                              "-2", infile_list[1],
                              "-b", str(fullpath + ".bitmask"),
                              "-x", str(fullpath + ".srprism"),
                              "-T", temp_dir, 
                              "-o", out_prefix,
                              "--extract"]
                outputs[i] = [out_prefix + ending for ending in
                        const.BMTAGGER_PE_ENDINGS]
                                            
            else:
                out_prefix = prefix + "_" + basename + "_contam.out"
                bmt_args[i] = [bmtagger_path, 
                              "-q", "1",
                              "-1", infile_list[0],
                              "-2", infile_list[1],
                              "-b", str(fullpath + ".bitmask"),
                              "-x", str(fullpath + ".srprism"),
                              "-T", temp_dir, 
                              "-o", out_prefix]
                outputs[i] = [out_prefix]
        if debug:
            bmt_args[i] += ["--debug"]
    
    logging.debug("BMTagger commands to run:")
    logging.debug(bmt_args)

    if not bmt_args: # no databases specified
        return ([], [], [])

    # similar to what we did for Bowtie2
    # Poll to see if we can run more BMTagger instances. 
    procs_ran = list()
    procs_running = list()
    while bmt_args:
        cmd = bmt_args.pop() 
        n_running, procs_running = _poll_workers(procs_running)
        if n_running >= n_procs:
            bmt_args.append(cmd)
            time.sleep(0.5)
        else:
            logging.debug("Running BMTagger command: " + str(cmd))
            proc = subprocess.Popen(cmd, stderr=subprocess.PIPE,
                                    stdout=subprocess.PIPE)
            procs_running.append((proc, cmd))
            procs_ran.append(cmd)

    # wait for everything to finish after we started running all the processes
    for p, cmd in procs_running:
        stdout, stderr = p.communicate()
        logging.debug("BMTagger command `%s' Complete", cmd)
        if stdout:
            logging.debug("BMTagger command `%s' stdout: %s", cmd, stdout)
        if stderr:
            logging.debug("BMTagger command `%s' stderr: %s", cmd, stderr)

    ret_codes = [p.returncode for p, _ in procs_running]

    # if BMTagger produced correct output, merge the files from multiple
    # databases
    #if (outputs != []):
    set_retcodes = set(ret_codes)
    combined_outs = []
    if (len(set_retcodes) == 1) and (0 in set_retcodes):
        combined_outs = combine_tag(outputs, prefix)

    if not debug:
        for output_pair in outputs:
            for o in output_pair:
                try:
                    os.remove(o)
                except OSError as e:
                    logging.exception(e)
                    logging.warning("Could not remove file " + str(o))
                    

    logging.debug("processes run: %s", procs_ran)
    return (ret_codes, procs_ran, combined_outs)


# TODO: this is a pretty bad way of doing this. consider doing what Randall did
# in _trim_biopython
def check_fastq(strFname):
    '''
    Returns True if file strFname is a fastq file (based on file extension)
    '''
    isFastq = strFname.endswith('.fastq') or strFname.endswith('.fq')
    return (isFastq)

def intersect_fastq(lstrFiles, out_file):
    ''' 
    Intersects multiple fastq files with one another. Includes only the reads (4
    lines long each) that are common to all the files. Writes these reads to the
    output file specified in out_file. 

    Input:
    lstrFiles:  a list of fastq file names (as strings)
    out_file:   output file to write the results. 
    '''
    # optimize for the common case, where we are intersecting 1 file
    if len(lstrFiles) == 1:
        #shutil.copy(lstrFiles[0], out_file)
        os.rename(lstrFiles[0], out_file)
        return

    counter = collections.Counter()
    for fname in lstrFiles:
        with open(fname, "rU") as f:
            # nifty trick to read 4 lines at a time (fastq files have 4 lines
            # per read)
            for lines in itertools.izip_longest(*[f]*4):
                try:
                    read = ("".join(lines)).strip()
                except TypeError:
                    logging.critical("You probably passed in a bad fastq"
                                     " file, one that doesn't have 4 lines"
                                     " per read")
                    raise
                else:
                    counter[read] += 1

    num_files = len(lstrFiles)
    with open(out_file, "w") as f:
        for key in counter:
            # only includes reads that are in n or more files, for n input files
            if counter[key] >= num_files:
                f.write(key+"\n")
    return

def union_outfiles(lstrFiles, out_file):
    '''
    Union multiple .out files (lists of contaminant reads output by BMTagger).
    Includes 1 copy of every read name listed in any of the input files. 
    Input:
    lstrFiles:  a list of .out file names (as strings)
    out_file:   output file to write the results
    '''

    # optimize for the common case, where we are unioning 1 file
    if len(lstrFiles) == 1:
        #shutil.copy(lstrFiles[0], out_file)
        os.rename(lstrFiles[0], out_file)
        return

    counter = collections.Counter()
    for fname in lstrFiles:
        with open(fname, "rU") as f:
            for line in f:
                read = line.strip()
                counter[read] += 1
    with open(out_file, "w") as f:
        for key in counter:
            f.write(key+"\n")

    # TODO: Test which one (above or below) is faster

    '''
    # implementing union with command line tools
    strFiles = " ".join(lstrFiles)
    cmd = shlex.split("sort --unique " + strFiles)
    with open(out_file, "w") as f:
        subprocess.call(cmd, stdout=f)
    '''
    return 

def combine_tag(llstrFiles, out_prefix):
    '''
    Summary: Combines output if we run BMTagger/bowtie2 on multiple databases.
    Input:
        llstrFiles: a list of lists of BMTagger/bowtie2 outputs. The 'inner
        lists' are of length 1 or 2 
        out_prefix: Prefix for the output files. 

    Output:
        Returns a list of output files. Additionally, the log file is updated
        with read counts for the different input and output files.
    '''
    # get read counts
    msgs = [[str(f + ": " + str(get_num_reads(f))) for f in pair] for pair in
            llstrFiles]

    # flatten the list
    for msg in itertools.chain.from_iterable(msgs):
        logging.info("Read counts after tagging: %s", msg)

    single_end = True
    # Get the file names and concat them into strings (separated by spaces)
    fnames1 = [f[0] for f in llstrFiles]
    fnames2 = []
    if len(llstrFiles[0]) == 2:
        fnames2 = [f[1] for f in llstrFiles]
        single_end = False

    # Check if it was .fastq or not
    # Instead of this method, try passing in another parameter? 
    fIsFastq = check_fastq(fnames1[0])

    output_files = []
    if fIsFastq:
        # If BMTagger outputs fastq files, we only want the lines in common
        # after tagging against all the databases
        output_file = out_prefix + "_1.fastq"
        if single_end:
            output_file = out_prefix + ".fastq"

        intersect_fastq(fnames1, output_file)
        output_files.append(output_file)
        if fnames2 != []:
            if single_end:
                raise Exception("Cannot have single end file with two different .fastq's per database")
            output_file = out_prefix + "_2.fastq"
            # potential bug here. changed from fnames1 to fnames2
            intersect_fastq(fnames2, output_file)
            output_files.append(output_file)
    else:
        # Otherwise, if we generate .out files, we want the unique lines from
        # each, as to not have overlap when we produce the list of contaminant
        # reads
        union_outfiles(fnames1, out_prefix + ".out")
        output_files.append(out_prefix + ".out")
        if fnames2 != []:
            # This should not happen
            raise Exception("You have two different .out files for each database")

    # Get the read counts for the newly merged files
    logging.info("Read counts after merging from multiple databases: %s",
                 msg_num_reads(output_files))

    # remove output files from bowtie2/bmtagger if debug mode is not set
    if not logging.getLogger().isEnabledFor(logging.DEBUG):
        for group in [fnames1, fnames2]:
            if len(group) > 1:
                # if len(group) == 1, we renamed the files in intersect and
                # union
                for filename in group:
                    logging.info("Removing temporary file %s" %filename)
                    os.remove(filename)
    return output_files

def checkfile(fname, ftype="file", fail_hard=False):
    '''
    input:
        fnames:         One or more file names. Must be strings
        ftype:          A string describing what type of file we are checking. 
                        Used to clarify error messages.
        empty_failhard: Boolean, whether or not to raise an error when the file
                        is empty, or whether to return a value. 
        dne_failhard:   Boolean, whether or not to raise an error when the file
                        does not exist, or whether to return a value. 
    output:
        1:      the file exists
        0:      the file does not exist
        -1:     the file exists but is an empty file (size = 0 bytes)

        If fail_hard=True, then raises IOError on empty file and OSError on
        non-existent file instead of returning 0 and -1, respectively

    Summary: Helper function to test if a file exists and is nonempty, exists
    and is empty, or does not exist. Can fail loudly or silently
    '''
    try:
        if os.stat(fname).st_size > 0:
            return 1
        else:
            if fail_hard:
                raise IOError(str(fname + " is an empty " + ftype))
            else:
                return -1
    except OSError:
        if fail_hard:
            raise OSError(str("Could not find " + ftype + " " + fname))
        else:
            return 0

def checktrim_output(output_prefix, b_single_end):
    '''
    input:  output_prefix: a string containing the output prefix
            b_single_end:  True/False, single end or not

    output: a tuple (b, outputs, new_inputs)
            b:          a boolean. True if at least 1 Trimmomatic output file is
                        existing and nonempty
            outputs:    a list of the outputs we are checking against. 
            new_inputs: a list of lists [l_1, l_2, l_3, ...] where each l_i is
                        a (nonempty) input to BMTagger. Each l_i has length 1 or
                        length 2.
    Summary:    Checks if Trimmomatic output files exist, and if at least some
                of them do, assemble them into lists to pass to BMTagger.
    '''
    outputs = []
    ll_new_inputs = []
    if b_single_end:

        outputs.append(output_prefix + const.TRIM_SE_ENDING)
        checks = checkfile(outputs[0])
        if checks <= 0:
            return (False, outputs, ll_new_inputs)
        else:
            ll_new_inputs = [outputs]

    else:
        len_endings = len(const.TRIM_PE_ENDINGS)
        outputs = [output_prefix + ending for ending in const.TRIM_PE_ENDINGS]

        checks = [checkfile(out) for out in outputs]

        # check that all the files at least exist. If they don't exist this
        # means that something went wrong with Trimmomatic
        for i in range(len_endings):
            if checks[i] == 0:
                logging.critical("Could not find file %s",outputs[i])
                return(False, outputs, ll_new_inputs)
        
        if checks[0] == 1 and checks[1] == 1:
            ll_new_inputs.append([outputs[0], outputs[1]])
        elif checks[0] == 1:
            ll_new_inputs.append([outputs[0]])
        elif checks[1] == 1:
            ll_new_inputs.append([outputs[1]])
        for i in [2,3]:
            if checks[i] == 1:
                ll_new_inputs.append([outputs[i]])

        # If no valid inputs, return False
        if ll_new_inputs == []:
            return(False, outputs, ll_new_inputs)

    return (True, outputs, ll_new_inputs)


def get_num_reads(strFname):
    '''
    input: strFname: file name of the list of contaminant reads, or of the 
                    outputted fastq files from BMTagger.
           fIsFastq: is the file a fastq file or not.     
           output: the number of reads in the file, aka the number of lines in
                   the file divided by 4. If the command fails, return None
    Summary: Uses wc to find the number of reads in a file.
    '''
    pat = r'[0-9]+ '
    cmd = ["wc",  "--lines", strFname]
    fIsFastq = check_fastq(strFname)

    try:
        out = subprocess.check_output(cmd)
    except subprocess.CalledProcessError as e:
        logging.warning("Command %s failed with return code %i ",
                        str(e.cmd), e.returncode)
        return None
    else:
        # match to get the line numbers
        match = re.search(pat, out)
        if match:
            num_reads = int(match.group())
            if fIsFastq:
                # one read is 4 lines in the fastq file
                num_reads = num_reads/4

            return(num_reads)
        else:
            # should not happen
            logging.critical("This should not happen. "
                             "Can't find the regex in wc output!")
            return None
   
   
def msg_num_reads(lstrFiles):
    ''' 
    Joins the read counts for a list of files into a printable message, with one
    file per line.
    Input: lstrFiles: a list of either .fastq or .out (list of read headers)
                      files. 

    Returns: newline-separated string of the form 
        f: n
    f is the name of the file
    n is the number of reads
    '''
    return ("\n".join([f + ": " + str(get_num_reads(f)) for f in lstrFiles]))

def check_missing_files(args):
    """check inputs"""
    # deal with missing prefix
    if args.output_prefix == None:
        args.output_prefix = args.infile1 + "_output"

    # check for the existence of required files/paths
    paths = [args.infile1]
    if args.infile2 != None:
        paths.append(args.infile2)
    if args.trim_path:
        paths.append(args.trim_path)
    
    for p in paths:
        checkfile(p, fail_hard=True)

    for db_prefix in args.reference_db:
        dbs = None
        if args.bmtagger:
            dbs = map(lambda x: str(db_prefix + x), const.BMTAGGER_DB_ENDINGS)
        else:
            dbs = map(lambda x: str(db_prefix + x), const.BOWTIE2_DB_ENDINGS)
        logging.debug(dbs)
        checks = [ ( checkfile( db, ftype="reference database", 
                                fail_hard=False), db) 
                   for db in dbs ]
        logging.debug(checks)
        for check, db in checks:
            if check == 0:
                logging.critical("Could not find reference database file %s",
                                 db)
                sys.exit(1)
                

def dict_to_cmd_opts_iter(opts_dict, sep=" ", singlesep=" "):
    """sep separates long options and their values, singlesep separates
    short options and their values e.g. --long=foobar vs -M 2

    """

    for key, val in opts_dict.iteritems():
        key = key.replace('_', '-')
        if len(key) > 1:
            key = "--%s"% (key)
        else:
            key = "-%s"% (key)
            
        if val:
            yield key
            yield val
        else:
            yield key


def _prefix_bases(db_prefix_list):
    """From a list of absolute or relative paths, returns an iterator of the
    following tuple:
    (basename of all the files, full path of all the files)
    
    If more than one file as the same basename (but have different paths), the
    basenames get post-fixed with indices (starting from 0), based on how their
    basenames sort

    :param db_prefix_list: A list of databases' paths. Can be absolute or
    relative paths.
    """
    # sort by basename, but keep full path
    bases = sorted([ (os.path.basename(p), p) for p in db_prefix_list ], 
            key = lambda x: x[0])
    for name, group in itertools.groupby(bases, key=lambda x: x[0]):
        group = list(group)
        if len(group) > 1:
            for i, item in enumerate(group):
                yield ("%s_%i"%(item[0], i), item[1])
        else:
            yield (group[0][0], group[0][1])


def find_on_path(bin_str):
    """ Finds an executable living on the shells PATH variable.
    :param bin_str: String; executable to find

    :returns: Absolute path to `bin_str` or False if not found
    :rtype: str
    """

    for dir_ in os.environ['PATH'].split(':'):
        candidate = os.path.join(dir_, bin_str)
        if os.path.exists(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return False


def trim(infile, prefix, trimmomatic_path, 
         java_mem="500m", addl_args=list(), threads=1):
    '''
    Trim a sequence file using trimmomatic. 
        infile:     input fastq file list (either length 1 or length 2)
                    length 1: single end
                    length 2: paired end
        prefix:     prefix for outputs
        java_mem:   string, ie "500m" or "8g"; specifies how much memory
                    the Java VM is allowed to use
        addl_args:  list of string, additional arguments for Trimmomatic
        threads:    int, number of threads for Trimmomatic to use

    output: 4 files for paired end
        (1) prefix.trimmed.1.fastq: trimmed first pair-end file
        (2) prefix.trimmed.2.fastq: trimmed second pair-end file
        (3) prefix.trimmed.single.1.fastq: trimmed sequences from the first file
        that lost their partner from the second file
        (4) prefix.trimmed.single.2.fastq: trimmed sequences from the second
        file that lost their partner from the first file

            1 file for single end:
        prefix.trimmed.fastq

    returns:
        (res, cmd)              res is the error code of the command
                                cmd (string) is the command itself
    Summary: Calls Trimmomatic to trim reads based on quality
    '''

    single_end = (len(infile) == 1)
    trim_cmd = ["java", 
                "-Xmx" + java_mem, 
                "-d64",
                "-jar", trimmomatic_path]

    if single_end:
        trim_cmd += ["SE",
                    "-threads", str(threads), 
                    "-phred33", 
                    infile[0],
                    prefix + ".trimmed.fastq"] + addl_args
    else:
        trim_cmd += ["PE", 
                    "-threads", str(threads), 
                    "-phred33", 
                    infile[0], infile[1], 
                    prefix + ".trimmed.1.fastq", 
                    prefix + ".trimmed.single.1.fastq",
                    prefix + ".trimmed.2.fastq",
                    prefix + ".trimmed.single.2.fastq"] + addl_args

    logging.debug("Running trimmomatic with: " + " ".join(trim_cmd))
    proc = subprocess.Popen(trim_cmd,
                            stderr=subprocess.PIPE,
                            stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    retcode = proc.returncode
    process_return("Trimmomatic", retcode, stdout, stderr)
    return(retcode, trim_cmd)


def run_trf(fastqs, outs, match=2, mismatch=7, delta=7, pm=80, pi=10,
        minscore=50, maxperiod=500, generate_fastq=True, mask=False, html=False,
        trf_path="trf", n_procs=None):
    nfiles = len(fastqs)
    # 1-2 files being passed in
    assert(nfiles == 2 or nfiles == 1)

    procs = list()
    names = list()
    for (fastq, out) in zip(fastqs, outs):
        proc = memoryheavy.run_tandem(fastq, out, match, mismatch, delta, pm,
                pi, minscore, maxperiod, generate_fastq,mask,html, trf_path)
        procs.append(proc)
        names.append("tandem.py on " + fastq)

    for (proc, name) in zip(procs, names):
        stdout, stderr = proc.communicate()
        retcode = proc.returncode
        process_return(name, retcode, stdout, stderr)

def storage_heavy(args):
    check_missing_files(args)

    trim_threads, bowtie_threads = divvy_threads(args)
    output_prefix = os.path.join(args.output_dir, args.output_prefix)

    # determine single-ended or pair ends
    b_single_end = (args.infile2 is None)
    files = [args.infile1] if b_single_end else [args.infile1, args.infile2]

    # Get number of reads initially, then log
    msg = msg_num_reads(files)
    logging.info("Initial number of reads: "+msg)

    logging.info("Trimming...")
    trim(files, 
         threads          = trim_threads,
         trimmomatic_path = args.trim_path, 
         prefix           = output_prefix,
         java_mem         = args.max_mem, 
         addl_args        = args.trim_args)
    
    # TODO: run part of the pipeline (if you have the output, either overwrite
    # or do the next step)

    logging.info("Finished running Trimmomatic. "
                 "Checking output files exist... ")

    # check that Trimmomatic's output files exist
    b_continue, outputs, files_to_align = checktrim_output(output_prefix, 
                                                           b_single_end)

    msg_trim_body = msg_num_reads(outputs)
    logging.info("Number of reads after trimming: %s", msg_trim_body)

    if not b_continue:
        logging.critical("Trimmomatic just produced empty files.")
        logging.critical("Terminating the pipeline...")
        sys.exit(1)

    # make temporary directory for Bowtie2 files
    tempdir = output_prefix + "_temp"
    try_create_dir(tempdir)

    # TODO: Add parallelization. Use command line utility 'split' to split the
    # files, fork a thread for each file, and run BMTagger in each thread.

    # determine if logging level is debug or not:

    debug = logging.getLogger().isEnabledFor(logging.DEBUG)

    # Start aligning
    logging.info("Decontaminating")
    possible_orphan = (len(files_to_align) > 1)
    orphan_count = 1
    for files_list in files_to_align:
        prefix = output_prefix
        if possible_orphan and (len(files_list) == 1):
            prefix = output_prefix + "_se_" + str(orphan_count)
            orphan_count += 1
        elif len(files_list) == 2:
            prefix = output_prefix + "_pe"

        trf_out_base = None
        if args.trf:
            trf_out_base = prefix
            prefix = prefix + "_pre_tandem"

        if args.bmtagger:
            ret_codes, commands, c_outs = tag(infile_list    = files_list,
                                              db_prefix_list = args.reference_db,
                                              temp_dir       = tempdir,
                                              prefix         = prefix,
                                              bmtagger_path  = args.bmtagger_path,
                                              n_procs        = bowtie_threads,
                                              remove         = args.extract,
                                              debug          = debug)
        else:
            ret_codes, commands, c_outs = align(infile_list   = files_list,
                                               db_prefix_list = args.reference_db,
                                               output_prefix  = prefix,
                                               tmp_dir        = tempdir,
                                               bowtie2_path   = args.bowtie2_path,
                                               n_procs        = bowtie_threads,
                                               bowtie2_opts   = args.bowtie2_args)
        # run TRF (within loop)
        # iterate over all outputs from combining (there should either be 1 or
        # 2)
        if args.trf:
            trf_outs = [trf_out_base]
            if len(c_outs) == 2:
                trf_outs = [trf_out_base + str(i+1) for i in xrange(2)]
            run_trf(fastqs = c_outs, 
                    outs = trf_outs, 
                    match = args.match,
                    mismatch = args.mismatch, 
                    delta = args.delta, 
                    pm = args.pm,
                    pi = args.pi, 
                    minscore = args.minscore, 
                    maxperiod = args.maxperiod, 
                    generate_fastq = args.no_generate_fastq, 
                    mask = args.mask, 
                    html = args.html, 
                    trf_path = args.trf_path,
                    n_procs = bowtie_threads)
            if not debug:
                for c_out in c_outs:
                    os.remove(c_out)

    # check that everything returned correctly
    # gather non-zero return codes
    fails = [(i, ret_code) for i, ret_code in enumerate(ret_codes)
             if ret_code != 0]
    if len(fails) > 0:
        for i, ret_code in fails:
            logging.critical("The following command failed with return code %s: %s",
                         ret_code, commands[i])
        sys.exit(1)


    logging.info("Finished removing contaminants")

    #if not logging.getLogger().isEnabledFor(logging.DEBUG):
    if not debug:
        logging.debug("Removing temporary files...")
        # this removes trimmomatic outputs
        for output in outputs:
            os.remove(output)
        # this removes lots of temporary files generated by bowtie2/bmtagger
        shutil.rmtree(tempdir)

