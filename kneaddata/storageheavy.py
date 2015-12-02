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
import gzip
import tempfile

from . import utilities
from . import config
from . import memoryheavy

# name global logging instance
logger=logging.getLogger(__name__)

def _generate_bowtie2_commands( infile_list, db_prefix_list,
                                bowtie2_path, output_prefix,
                                bowtie2_opts, tmp_dir, threads, remove_temp_output ):
    db_prefix_bases = _prefix_bases(db_prefix_list)
    for basename, fullpath in db_prefix_bases:
        is_paired = (len(infile_list) == 2)
        cmd = [bowtie2_path]
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

        cmd += list(utilities._get_bowtie2_args(bowtie2_opts))
        if remove_temp_output:
            # if we are removing the temp output, then write the sam output to dev null to save space
            sam_out = os.devnull
        else:
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
        command_msg = "\n".join([ " ".join(cmd) + " Returned code " + str(val) 
                for val, cmd in failures ])
        message="The following command failed: " + command_msg
        logger.critical(message)
        sys.exit(message)
    else:
        return (still_running, procs_still_running)


def align(infile_list, db_prefix_list, output_prefix, tmp_dir, remove_temp_output,
          bowtie2_path, n_procs=None, bowtie2_opts=list(),verbose=None):

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
    :param remove_temp_output: Boolean; If set, remove temp output
    :keyword bowtie2_path: String; File path to the bowtie2 binary's location.
    :keyword n_procs: Int; Number of bowtie2 subprocesses to run.
    :keyword bowtie2_opts: List; List of additional arguments, as strings, to 
                            be passed to Bowtie2.
    """

    if not n_procs:
        n_procs = len(db_prefix_list)

    commands_to_run = _generate_bowtie2_commands( infile_list,
        db_prefix_list, bowtie2_path, output_prefix, bowtie2_opts,
        tmp_dir, n_procs, remove_temp_output)
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
            utilities.log_run_and_arguments("bowtie2",cmd[1:],verbose)
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            procs_running.append((proc, cmd))
            outputs.append(output_to_combine)

    # wait for everything to finish after we started running all the processes
    commands_with_errors=0
    ret_codes=[]
    commands_ran=[]
    for p, cmd in procs_running:
        stdout, stderr = p.communicate()
        cmd = " ".join(cmd)
        
        # check if the command finished without an error
        returncode = p.returncode
        ret_codes.append(returncode)
        commands_ran.append(cmd)
        if returncode == 0:
            message="bowtie2 run complete"
            if verbose:
                print(message)
            logger.debug(message)
        else:
            message="bowtie2 run reported error"
            print(message)
            logger.warning(message)
            
        if stdout:
            logger.debug("bowtie2 run stdout: " + stdout)
        if stderr:
            logger.debug("bowtie2 run stderr: " + stderr)

   
    # if Bowtie2 produced correct output, merge the files from multiple
    # databases
    combined_outs = []
    if (outputs != []):
        combined_outs = combine_tag(outputs, output_prefix, remove_temp_output)

    return(ret_codes, commands_ran, combined_outs)


def tag(infile_list, db_prefix_list, temp_dir, remove_temp_output, prefix,
        bmtagger_path, n_procs=2, remove=False, verbose=None):
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
                           config.bmtagger_db_endings for full list of endings).
                           Multiple databases can be specified. Uses
                           subprocesses to run BMTagger in parallel
    :param temp_dir: String; Path name to temporary directory for BMTagger
    :param remove_temp_output: Boolean; If set, remove temp output
    :param prefix: String; prefix for output files
    :keyword bmtagger_path: String; file path to the bmtagger.sh executable
    :keyword n_procs: Int; Max number of BMTagger subprocesses to run
    :keyword remove: Boolean; If True, output cleaned FASTQ files with reads
                     identified as contaminants removed. Otherwise, outputs
                     list(s) of reads identified as contaminants.
    """
    single_end = (len(infile_list) == 1)

    message="bmtagger prefix list "+db_prefix_list
    if verbose:
        print(message)
    logger.debug(message)
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
                outputs[i] = [out_prefix + config.bmtagger_se_ending]

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
                        config.bmtagger_pe_endings]
                                            
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
            utilities.log_run_and_command("bmtagger",cmd[1:],verbose)
            proc = subprocess.Popen(cmd, stderr=subprocess.PIPE,
                                    stdout=subprocess.PIPE)
            procs_running.append((proc, cmd))
            procs_ran.append(cmd)

    # wait for everything to finish after we started running all the processes
    for p, cmd in procs_running:
        stdout, stderr = p.communicate()
        message="bmtagger run complete"
        if verbose:
            print(message)
        logger.debug(message)
        if stdout:
            logger.debug("bmtagger run stdout: " + stdout)
        if stderr:
            logger.debug("bmtagger run stderr: " + stderr)

    ret_codes = [p.returncode for p, _ in procs_running]

    # if BMTagger produced correct output, merge the files from multiple
    # databases
    #if (outputs != []):
    set_retcodes = set(ret_codes)
    combined_outs = []
    if (len(set_retcodes) == 1) and (0 in set_retcodes):
        combined_outs = combine_tag(outputs, prefix, remove_temp_output)

    if remove_temp_output:
        for output_pair in outputs:
            for o in output_pair:
                try:
                    os.remove(o)
                except OSError as e:
                    message="Could not remove temp file: " + str(o)
                    if verbose:
                        print(message)
                    logger.debug(message)
                    
    return (ret_codes, procs_ran, combined_outs)

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
                    message="Fastq file is not correctly formatted"
                    print(message)
                    logger.critical(message)
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

def combine_tag(llstrFiles, out_prefix, remove_temp_output):
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
    
    # print out the reads for all files
    all_files=[]
    for pair in llstrFiles:
        all_files+=pair
    utilities.log_read_count_for_files(all_files,"Total reads after tagging")

    single_end = True
    # Get the file names and concat them into strings (separated by spaces)
    fnames1 = [f[0] for f in llstrFiles]
    fnames2 = []
    if len(llstrFiles[0]) == 2:
        fnames2 = [f[1] for f in llstrFiles]
        single_end = False

    # Check if it was .fastq or not
    # Instead of this method, try passing in another parameter? 
    fIsFastq = utilities.is_file_fastq(fnames1[0])

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
    utilities.log_read_count_for_files(output_files,"Total reads after merging results from multiple databases")

    # remove temp files if set
    if remove_temp_output:
        for group in [fnames1, fnames2]:
            if len(group) > 1:
                # if len(group) == 1, we renamed the files in intersect and
                # union
                for filename in group:
                    logger.debug("Removing temporary file %s" %filename)
                    os.remove(filename)
    return output_files

def checktrim_output(output_prefix, input_files):
    '''
    input:  output_prefix: a string containing the output prefix
            input_files: the list of input files

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
    
    # determine if paired input files
    b_single_end = True
    if len(input_files) == 2:
        b_single_end = False
    
    if b_single_end:
        outputs.append(output_prefix + config.trimomatic_se_ending)
        if not utilities.is_file_nonempty_readable(outputs[0]):
            return (False, outputs, ll_new_inputs)
        else:
            ll_new_inputs = [outputs]

    else:
        len_endings = len(config.trimomatic_pe_endings)
        outputs = [output_prefix + ending for ending in config.trimomatic_pe_endings]

        checks = [utilities.is_file_nonempty_readable(out) for out in outputs]

        # check that all the files at least exist. If they don't exist this
        # means that something went wrong with Trimmomatic
        for i in range(len_endings):
            if not checks[i]:
                message="File does not exist or is empty: " + outputs[i]
                print(message)
                logger.critical(message)
                return(False, outputs, ll_new_inputs)
        
        if checks[0] and checks[1]:
            ll_new_inputs.append([outputs[0], outputs[1]])
        elif checks[0]:
            ll_new_inputs.append([outputs[0]])
        elif checks[1]:
            ll_new_inputs.append([outputs[1]])
        for i in [2,3]:
            if checks[i]:
                ll_new_inputs.append([outputs[i]])

        # If no valid inputs, return False
        if ll_new_inputs == []:
            return(False, outputs, ll_new_inputs)

    return (True, outputs, ll_new_inputs)

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

def trim(infile, prefix, trimmomatic_path, quality_scores, 
         java_mem="500m", addl_args=list(), threads=1, verbose=None):
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
                    quality_scores, 
                    infile[0],
                    prefix + config.trimomatic_se_ending] + addl_args
    else:
        trim_cmd += ["PE", 
                    "-threads", str(threads),
                    quality_scores,
                    infile[0], infile[1], 
                    prefix + config.trimomatic_pe_endings[0], 
                    prefix + config.trimomatic_pe_endings[2],
                    prefix + config.trimomatic_pe_endings[1],
                    prefix + config.trimomatic_pe_endings[3]] + addl_args

    utilities.log_run_and_arguments("trimmomatic",trim_cmd,verbose)
    proc = subprocess.Popen(trim_cmd,
                            stderr=subprocess.PIPE,
                            stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    retcode = proc.returncode
    utilities.process_return("Trimmomatic", retcode, stdout, stderr)
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
        utilities.process_return(name, retcode, stdout, stderr)
        
def decontaminate(args, bowtie_threads, output_prefix, files_to_align):
    """
    Run bowtie2 or bmtagger then trf if set
    """

    # make temporary directory for temp output files
    if args.remove_temp_output:
        # create a temp folder if we are removing the temp output files
        tempdir=tempfile.mkdtemp(prefix=args.output_prefix+'_kneaddata_temp_',dir=args.output_dir)
    else:
        tempdir = output_prefix + "_temp"
        utilities.create_directory(tempdir)

    # Start aligning
    message="Decontaminating ..."
    print(message)
    logger.info(message)
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
            ret_codes, commands, c_outs = tag(infile_list=files_list,
                                              db_prefix_list=args.reference_db,
                                              temp_dir=tempdir,
                                              remove_temp_output=args.remove_temp_output,
                                              prefix=prefix,
                                              bmtagger_path=args.bmtagger_path,
                                              n_procs=bowtie_threads,
                                              remove=args.extract,
                                              verbose=args.verbose)
        else:
            ret_codes, commands, c_outs = align(infile_list=files_list,
                                                db_prefix_list=args.reference_db,
                                                output_prefix=prefix,
                                                tmp_dir=tempdir,
                                                remove_temp_output=args.remove_temp_output,
                                                bowtie2_path=args.bowtie2_path,
                                                n_procs=bowtie_threads,
                                                bowtie2_opts=args.bowtie2_options,
                                                verbose=args.verbose)
                
        # check that everything returned correctly
        # gather non-zero return codes
        fails = [(i, ret_code) for i, ret_code in enumerate(ret_codes) if ret_code != 0]
        if len(fails) > 0:
            for i, ret_code in fails:
                message="The following command failed: " + " ".join(commands[i])
                logger.critical(message)
                print(message)
            sys.exit(1)
        
        # run TRF (within loop)
        # iterate over all outputs from combining (there should either be 1 or
        # 2)
        if args.trf:
            trf_outs = [trf_out_base]
            if len(c_outs) == 2:
                trf_outs = [trf_out_base + str(i + 1) for i in xrange(2)]
            run_trf(fastqs=c_outs,
                    outs=trf_outs,
                    match=args.match,
                    mismatch=args.mismatch,
                    delta=args.delta,
                    pm=args.pm,
                    pi=args.pi,
                    minscore=args.minscore,
                    maxperiod=args.maxperiod,
                    generate_fastq=args.no_generate_fastq,
                    mask=args.mask,
                    html=args.html,
                    trf_path=args.trf_path,
                    n_procs=bowtie_threads)
            if args.remove_temp_output:
                for c_out in c_outs:
                    os.remove(c_out)
        message="Finished removing contaminants"
        logger.info(message)
        print(message)    

    if args.remove_temp_output:
        logger.debug("Removing temporary alignment files ...")
        # this removes lots of temporary files generated by bowtie2/bmtagger
        shutil.rmtree(tempdir)    
    

def storage_heavy(args):

    trim_threads, bowtie_threads = utilities.divvy_threads(args)
    output_prefix = os.path.join(args.output_dir, args.output_prefix)

    # Get number of reads initially, then log
    utilities.log_read_count_for_files(args.input,"Initial number of reads",args.verbose)

    message="Trimming ..."
    logger.info(message)
    print(message)
    trim(args.input, 
         threads          = trim_threads,
         trimmomatic_path = args.trimmomatic_path,
         quality_scores   = args.trimmomatic_quality_scores,
         prefix           = output_prefix,
         java_mem         = args.max_mem, 
         addl_args        = args.trimmomatic_options,
         verbose          = args.verbose)

    message="Finished running trimmomatic"
    if args.verbose:
        print(message)
    logger.info(message)

    # check that Trimmomatic's output files exist
    b_continue, outputs, files_to_align = checktrim_output(output_prefix, args.input)

    utilities.log_read_count_for_files(outputs,"Total reads after trimming",args.verbose)

    if not b_continue:
        message="Trimmomatic did not output any trimmed reads"
        logger.critical(message)
        sys.exit(message)

    # Start aligning
    if not args.reference_db:
        message="Bypass decontamination"
        logger.info(message)
        print(message)
    else:
        decontaminate(args, bowtie_threads, output_prefix, files_to_align)
        # remove temp trimmomatic files
        if args.remove_temp_output:
            logger.debug("Removing temporary trim files...")
            for output in outputs:
                os.remove(output)

