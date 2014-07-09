'''
knead_data.py
Author: Andy Shi

Pipeline for processing metagenomics sequencing data
'''
import argparse
import subprocess
import shlex
import os
import sys
import shutil
import re
import itertools
import collections
import constants_knead_data as const


def trim(infile, trimlen, trim_path, single_end, prefix, mem, addl_args):
    '''
    input: 
        infile:     input fastq file list (either length 1 or length 2)
                    length 1: single end
                    length 2: paired end
        trimlen:    length to trim
        trim_path:  Path to the Trimmomatic executable
        single_end: True/False, is it a single_end or paired end input
        prefix:     prefix for outputs
        mem:        string, ie "500m" or "8g", which specifies how much memory
                    the Java VM is allowed to use
        addl_args:  string of additional arguments for Trimmomatic
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

    # check that we have the right number of input fastqs
    assert((len(infile) == 2 and not single_end) or (single_end and len(infile)
        == 1))


    trim_arg = ""
    if single_end:
        trim_arg = str(trim_path + " SE -phred33 " + infile[0] + " " + prefix +
                ".trimmed.fastq " + "MINLEN:" + str(trimlen) + " " + addl_args)
    else:
        trim_arg = str(trim_path + " PE -phred33 " + infile[0] + " " +
                infile[1] + " " + prefix + ".trimmed.1.fastq " + prefix +
                ".trimmed.single.1.fastq " + prefix + ".trimmed.2.fastq " +
                prefix + ".trimmed.single.2.fastq " + "MINLEN:" + str(trimlen) +
                " " + addl_args)

    cmd = "java -Xmx" + mem + " -d64 -jar " + trim_arg
    print("Trimmomatic command that will be run: " + cmd)
    ret = subprocess.call(shlex.split(cmd))
    return(ret, cmd)
    #proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
    #stdout, stderr = proc.communicate()
    #return (proc.returncode, cmd)


def tag(infile, db_prefix, bmtagger_path, single_end, prefix, remove, debug,
        temp_dir, logfile, orphan=None):
    '''
    Runs multiple instances of BMTagger (in parallel), one for each input
    database. 
    input:
        infile:         a list of length 1 or 2 (for single and paired ends,
                        respectively) 
        db_prefix:      prefix of the BMTagger databases. BMTagger needs 
                        databases of the format [db_prefix].srprism.* and 
                        [db_prefix].[blastdb file extensions]
        bmtagger_path:  path to the bmtagger.sh executable
        single_end:     True/False, single end or paired ends
        prefix:         prefix for output files
        remove:         True/False, remove the reads from the input files 
                        (make a copy first) or not. 
        debug:          Boolean. If True, BMTagger does not delete temporary
                        files
        temp_dir:       Directory to put BMTagger's temporary files

    output:
        prefix.out:             a list of the contaminant reads
        prefix_depleted.fastq:  the input fastq with the unwanted reads removed

    returns:
        (ret_codes, bmt_args)   
        Both are lists, each which has the same number of elements as the number
        of databases. THe first contains the return codes for each of the
        BMTagger executions, and the second contains the command line calls for
        these executions
    Summary: Uses BMTagger to tag and potentially remove unwanted reads
    '''
    # quick check of inputs
    if single_end:
        if len(infile) != 1:
            print("Improper call to BMTagger!")
            return([],[])

    else:
        if len(infile) != 2:
            print("Improper call to BMTagger!")
            return([],[])


    db_len = len(db_prefix)
    bmt_args = [None for i in range(db_len)]

    # correctly manage the output prefix
    outputs = [None for i in range(db_len)]
    # build arguments
    for i in range(db_len):
        db = db_prefix[i]
        out_prefix = None
        if single_end:
            if remove:
                out_prefix = prefix + "_db" + str(i)
                if orphan != None:
                    out_prefix = prefix + "_db" + str(i) + "_se_" + str(orphan)
                bmt_args[i] = str(bmtagger_path + " -q 1 -1 " + infile[0] + 
                        " -b " + db + ".bitmask -x " + db + ".srprism -T " + 
                        temp_dir + " -o " + out_prefix + " --extract")
                outputs[i] = [out_prefix + const.BMTAGGER_SE_ENDING]

            else:
                out_prefix = prefix + "_db" + str(i) + ".out"
                if orphan != None:
                    out_prefix = str(prefix + "_db" + str(i) + "_se_" + str(orphan)
                            + ".out")
                bmt_args[i] = str(bmtagger_path + " -q 1 -1 " + infile[0] + 
                        " -b " + db + ".bitmask -x " + db + ".srprism -T " + 
                        temp_dir + " -o " + out_prefix)
                outputs[i] = [out_prefix]


        else:
            if remove:
                out_prefix = prefix + "_db" + str(i) + "_pe"
                bmt_args[i] = str(bmtagger_path + " -q 1 -1 " + infile[0] + 
                        " -2 " + infile[1] + " -b " + db + ".bitmask -x " + db +
                        ".srprism -T " + temp_dir + " -o " + out_prefix + 
                        " --extract")
                outputs[i] = [out_prefix + ending for ending in
                        const.BMTAGGER_PE_ENDINGS]
            
            else:
                out_prefix = prefix + "_db" + str(i) + "_pe" + ".out"
                bmt_args[i] = str(bmtagger_path + " -q 1 -1 " + infile[0] + 
                        " -2 " + infile[1] + " -b " + db + ".bitmask -x " + db +
                        ".srprism -T " + temp_dir + " -o " + out_prefix)
                outputs[i] = [out_prefix]
        if debug:
            bmt_args[i] = bmt_args[i] + " --debug"
    
    print("BMTagger commands to run:")
    print(bmt_args)
    procs = [subprocess.Popen(shlex.split(arg)) for arg in bmt_args]

    ret_codes = []
    # wait for all subprocesses to terminate
    for p in procs:
        p.wait()
        ret_codes.append(p.returncode)

    # if BMTagger produced correct output, merge the files from multiple
    # databases
    if (outputs != []) and (all(ret == 0 for ret in ret_codes)):
        combined_prefix = prefix
        if orphan != None:
            combined_prefix = prefix + "_se_" + str(orphan)
        combine_tag(outputs, logfile, combined_prefix, single_end)

    if not debug:
        for output_pair in outputs:
            for o in output_pair:
                try:
                    os.remove(o)
                except OSError as e:
                    print("Could not remove file " + str(o))
                    print("OS Error {0}: {1}".format(e.errno, e.strerror))

    return (ret_codes, bmt_args)


def check_fastq(strFname):
    '''
    Returns True if file strFname is a fastq file (based on file extension)
    '''
    return (strFname[-6:] == '.fastq')

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
        shutil.copy(lstrFiles[0], out_file)
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
                    print("Error encountered!!")
                    print("You probably passed in a bad fastq file, one that doesn't have 4 lines per read")
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
        shutil.copy(lstrFiles[0], out_file)
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
    return 

# TODO: Test which one (above or below) is faster

'''
# Dead code, for now. Tried to implement union using command line tools. 
def union_outfiles(lstrFiles, out_file):
    strFiles = " ".join(lstrFiles)
    # subprocess.call does not directly support redirection. Try to use
    # subprocess.Popen, and manually pipe the stdout to a file
    print(shlex.split("sort -u " + strFiles + " > " + out_file))
    #subprocess.call(shlex.split("sort -u " + strFiles + " > " + out_file))
'''

def combine_tag(llstrFiles, logfile, out_prefix, single_end):
    '''
    Summary: Combines output if we run BMTagger on multiple databases.
    Input:
        llstrFiles: a list of lists of BMTagger outputs. This is passed in from
        the tag function. The 'inner lists' are of length 1 or 2
        logfile: the file used to log statistics about the files

    Output:
        the log file is updated with read counts
    '''
    # get read counts
    msgs = [[str(f + ": " + str(get_num_reads(f))) for f in pair] for pair in llstrFiles]

    # flatten the list
    msg_flattened = itertools.chain.from_iterable(msgs)
    msg_to_print = "\n".join(msg_flattened)
    print("Read counts after tagging:")
    print(msg_to_print)
    with open(logfile, "a") as f:
        f.write("\nRead counts after tagging:\n")
        f.write(msg_to_print)

    # Get the file names and concat them into strings (separated by spaces)
    fnames1 = [f[0] for f in llstrFiles]
    fnames2 = []
    if len(llstrFiles[0]) == 2:
        fnames2 = [f[1] for f in llstrFiles]

    # Check if it was .fastq or not
    fIsFastq = check_fastq(fnames1[0])

    output_files = []
    # TODO: If it's only 1 file (common case), just copy it over and skip the
    # merging
    if fIsFastq:
        # If BMTagger outputs fastq files, we only want the lines in common
        # after tagging against all the databases
        output_file = out_prefix + "_pe_1.fastq"
        if single_end:
            output_file = out_prefix + ".fastq"

        intersect_fastq(fnames1, output_file)
        output_files.append(output_file)
        if fnames2 != []:
            if single_end:
                raise Exception("Cannot have single end file with two different .fastq's per database")
            output_file = out_prefix + "_pe_2.fastq"
            intersect_fastq(fnames1, output_file)
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
    joined_read_msgs_print = msg_num_reads(output_files)
    print("Read counts after merging from multiple databases:")
    print(joined_read_msgs_print)
    with open(logfile, "a") as f:
        f.write("\nRead counts after merging from multiple databases:\n")
        f.write(joined_read_msgs_print)
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

def is_single_end(file1, file2):
    '''
    Determine from two input fastq's if we are running single-end or paired-end.
    Outputs a list of file(s) that should be processsed. 
    Input:      file1:      an input fastq. INVARIANT: this argument is never 
                            None.
                file2:      possibly another input fastq. Could be None.

    Output:     (b_single_end, out_files)
                b_single_end:   True when single-end, False otherwise. 
                out_files:      A list of file(s) that will be processed by the
                                rest of the pipeline. If single-end, should just
                                be [file1]. If paired end, should be 
                                [file1, file2]
    If single-end, file2 should be None. Else, file2 should not be None. 
    '''
    b_single_end = True
    out_files = [file1]
    assert(file1)
    if file2:
        out_files.append(file2)
        b_single_end = False

    return (b_single_end, out_files)
    
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
                print("Could not find file " + outputs[i])
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
    cmd = "wc --lines " + strFname
    fIsFastq = check_fastq(strFname)

    try:
        out = subprocess.check_output(shlex.split(cmd))
    except subprocess.CalledProcessError as e:
        print("Command " + str(e.cmd) + " failed with return code " +
                str(e.returncode))
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
            print("This should not happen. Can't find the regex in wc output!")
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


def main():
    # parse command line arguments
    # note: argparse converts dashes '-' in argument prefixes to underscores '_' 
    parser = argparse.ArgumentParser()
    parser.add_argument("-1", "--infile1", help="input FASTQ file", required =
            True)
    parser.add_argument("-2", "--infile2", help="input FASTQ file mate")
    parser.add_argument("--trimlen", type=int, help="length to trim reads",
            default=60)
    parser.add_argument("-o", "--output-prefix",
            help="prefix for all output files")
    parser.add_argument("-db", "--reference-db", nargs = "+", default=[],
            help="prefix for reference databases used in BMTagger")
    # Consider using a params file
    parser.add_argument("-t", "--trim-path", help="path to Trimmomatic",
            required = True)
    parser.add_argument("-b", "--bmtagger-path", help="path to BMTagger",
            required = True)

    parser.add_argument("-x", "--extract", help="Remove contaminant reads",
            default=False, action="store_true")
    parser.add_argument("-m", "--max-mem", default="500m", 
            help="Maximum amount of memory that will be used by Trimmomatic, as a string, ie 500m or 8g")
    parser.add_argument("-a", "--trim-args", default="",
            help="additional arguments for Trimmomatic")
    parser.add_argument("-d", "--debug", default=False, action="store_true",
            help="If set, temporary files are not removed")
    # parser.add_argument("-S", "--slurm", help="Running in a slurm environment",
    #        action = "store_true")

    args = parser.parse_args()

    # check inputs
    # deal with missing prefix
    if args.output_prefix == None:
        args.output_prefix = args.infile1 + "_output"

    # check for the existence of required files/paths
    paths = [args.infile1, args.trim_path, args.bmtagger_path]
    if args.infile2 != None:
        paths.append(args.infile2)
    [checkfile(p, fail_hard=True) for p in paths]

    for db_prefix in args.reference_db:
        dbs = map(lambda x: str(db_prefix + x), const.DB_ENDINGS)
        [checkfile(db, ftype="BMTagger database", fail_hard=True) for db in dbs]

    # determine single-ended or pair ends
    b_single_end, files = is_single_end(args.infile1, args.infile2)

    # Get number of reads initially, then log
    # num_reads_init = [get_num_reads(f) for f in files]
    lenfiles = len(files)
    msg_init = "\nInitial number of reads:\n"
    msg = msg_num_reads(files)
    print(msg_init)
    print(msg)

    # Log file. Create a new one the first time, then keep appending to it.
    logfile = args.output_prefix + ".log"
    with open(logfile, "w") as f:
        f.write("Running knead_data.py with the following arguments (from argparse):\n"
                + str(args))
        f.write(msg_init)
        f.write(msg)

    print("Running Trimmomatic...")

    trim(files, trimlen = args.trimlen, trim_path = args.trim_path, single_end =
            b_single_end, prefix = args.output_prefix, mem = args.max_mem,
            addl_args = args.trim_args)
    
    # TODO: run part of the pipeline (if you have the output, either overwrite
    # or do the next step)

    print("Finished running Trimmomatic. Checking output files exist... ")

    # check that Trimmomatic's output files exist
    b_continue, outputs, bmt_inputs = checktrim_output(args.output_prefix, 
            b_single_end)

    msg_trim_init = "\nNumber of reads after trimming:\n"
    msg_trim_body = msg_num_reads(files)
    print(msg_trim_init)
    print(msg_trim_body)

    with open(logfile, "a") as f:
        f.write(msg_trim_init)
        f.write(msg_trim_body)

    if not b_continue:
        print("Trimmomatic produced no non-empty files.")
        print("Terminating the pipeline...")
        sys.exit(1)

    # make temporary directory for BMTagger files
    tempdir = args.output_prefix + "_temp"
    try:
        os.makedirs(tempdir)
    except OSError:
        if os.path.isdir(tempdir):
            print("BMTagger temporary directory already exists! Using it...")
        else:
            print("Cannot make the BMTagger temporary directory")
            raise

    # TODO: Add parallelization. Use command line utility 'split' to split the
    # files, fork a thread for each file, and run BMTagger in each thread.

    # Start tagging
    out_files = []

    if b_single_end:
        tag(infile = bmt_inputs[0], db_prefix = args.reference_db, 
                bmtagger_path = args.bmtagger_path, single_end = True, prefix =
                args.output_prefix, remove = args.extract, temp_dir = tempdir,
                debug = args.debug, logfile = logfile)
    else:
        orphan_counter = 1
        for inp in bmt_inputs:
            if len(inp) == 2:
                # Run paired end BMTagger
                tag(infile = inp, db_prefix = args.reference_db, bmtagger_path =
                        args.bmtagger_path, single_end = False, prefix =
                        args.output_prefix, remove = args.extract, temp_dir = tempdir, 
                        debug = args.debug, logfile = logfile)
            else:
                # Run single end BMTagger
                tag(infile = inp, db_prefix = args.reference_db, bmtagger_path =
                        args.bmtagger_path, single_end = True, prefix =
                        args.output_prefix, remove = args.extract, temp_dir = tempdir,
                        debug = args.debug, orphan = orphan_counter, logfile =
                        logfile)
                orphan_counter += 1

    print("Finished running BMTagger.")

    if not args.debug:
        print("Removing temporary files...")
        for output in outputs:
            os.remove(output)
        shutil.rmtree(tempdir)
    print("Done!")

if __name__ == '__main__':
    main()
