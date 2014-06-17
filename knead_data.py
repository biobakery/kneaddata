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

# Constants

# File endings for BMTagger's required database files
DB_ENDINGS =    [".bitmask", ".srprism.amp", 
                ".srprism.idx", ".srprism.imp",
                ".srprism.map", ".srprism.pmp", 
                ".srprism.rmp", ".srprism.ss",
                ".srprism.ssa", ".srprism.ssd", 
                ".nhr", ".nin", ".nsq"]

# Trimmomatic file endings for single end and paired end, respectively
TRIM_SE_ENDING = ".trimmed.fastq"

TRIM_PE_ENDINGS =   [".trimmed.1.fastq", 
                    ".trimmed.2.fastq", 
                    ".trimmed.single.1.fastq", 
                    ".trimmed.single.2.fastq"]

# BMTagger file endings if you choose to remove the contaminant reads. For
# single end and paired end reads, respectively.
BMTAGGER_SE_ENDING = ".fastq"
BMTAGGER_PE_ENDINGS =   ["_1.fastq",
                        "_2.fastq"]

def trim(infile, trimlen, trim_path, single_end, prefix, mem):
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
        trim_arg = str(trim_path + " SE -phred33 " + infile[0] + " " + prefix
                + ".trimmed.fastq " + "MINLEN:" + str(trimlen))
    else:
        trim_arg = str(trim_path + " PE -phred33 " + infile[0] + " " +
                infile[1] + " " + prefix + ".trimmed.1.fastq " + prefix +
                ".trimmed.single.1.fastq " + prefix + ".trimmed.2.fastq " +
                prefix + ".trimmed.single.2.fastq " + "MINLEN:" + str(trimlen))

    cmd = "java -Xmx" + mem + " -d64 -jar " + trim_arg
    print("Trimmomatic command that will be run: " + cmd)
    ret = subprocess.call(shlex.split(cmd))
    return(ret, cmd)
    #proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
    #stdout, stderr = proc.communicate()
    #return (proc.returncode, cmd)


def tag(infile, db_prefix, bmtagger_path, single_end, prefix, remove, temp_dir):
    '''
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
        temp_dir:       Directory to put BMTagger's temporary files

    output:
        prefix.out:             a list of the contaminant reads
        prefix_depleted.fastq:  the input fastq with the unwanted reads removed

    returns:
        (res, cmd)              res is the error code of the command
                                cmd (string) is the command itself
    Summary: Uses BMTagger to tag and potentially remove unwanted reads
    '''
    # check inputs
    db_len = len(db_prefix)
    assert (db_len > 0)

    bmt_args = ["" for i in range(db_len)]
    # build arguments
    for i in range(db_len):
        db = db_prefix[i]
        if single_end:
            bmt_args[i] = str(bmtagger_path + " -q 1 -1 " + infile[0] + " -b " +
                    db + ".bitmask -x " + db + ".srprism -T " + temp_dir + 
                    " -o " + prefix) 
        else:
            bmt_args[i] = str(bmtagger_path + " -q 1 -1 " + infile[0] + 
                    " -2 " + infile[1] + " -b " + db + ".bitmask -x " + db + 
                    ".srprism -T " + temp_dir + " -o " + prefix)
        if remove:
            # remove the contaminant reads
            bmt_args[i] = bmt_args[i] + " --extract"
    
    for arg in bmt_args:
        # Run all the BMTagger instances 
        print("BMTagger command to be run: " + arg)
        res = subprocess.call(shlex.split(arg))

    return (res,bmt_args)


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

        outputs.append(output_prefix + TRIM_SE_ENDING)
        checks = checkfile(outputs[0])
        if checks <= 0:
            return (False, outputs, ll_new_inputs)
        else:
            ll_new_inputs = [outputs]

    else:
        len_endings = len(TRIM_PE_ENDINGS)
        outputs = [output_prefix + TRIM_PE_ENDINGS[i] for i in
                        range(len_endings)]

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


def get_num_reads(fname):
    '''
    input:  fname: fastq file name, as a string
    output: the number of reads in the file, aka the number of lines in the file 
            divided by 4. If the command fails, return None
    Summary: Uses wc to find the number of reads in a file.
    '''
    pat = r'[0-9]+ '
    cmd = "wc --lines " + fname
    out = ""
    try:
        out = subprocess.check_output(shlex.split(cmd))
    except subprocess.CalledProcessError as e:
        print("Command " + str(e.cmd) + " failed with return code " +
                str(e.returncode))
        return None
    
    # match to get the line numbers
    match = re.search(pat, out)
    if match:
        # one read is 4 lines in the fastq file
        num_reads = int(match.group())/4
        return(num_reads)
    else:
        # should not happen
        print("This should not happen. Can't find the regex in wc output!")
        return None


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
    parser.add_argument("-db", "--reference-db", nargs = "+", 
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
    parser.add_argument("-S", "--slurm", help="Running in a slurm environment",
            action = "store_true")

    args = parser.parse_args()

    # check inputs
    # deal with missing prefix
    if args.output_prefix == None:
        args.output_prefix = args.infile1 + "_out"

    # if not extracting, append .out to the prefix
    if not args.extract:
        args.output_prefix = args.output_prefix + ".out"


    # check for the existence of required files/paths
    paths = [args.infile1, args.trim_path, args.bmtagger_path]
    if args.infile2 != None:
        paths.append(args.infile2)
    [checkfile(p, fail_hard=True) for p in paths]

    for db_prefix in args.reference_db:
        dbs = map(lambda x: str(db_prefix + x), DB_ENDINGS)
        [checkfile(db, ftype="BMTagger database", fail_hard=True) for db in dbs]

    # determine single-ended or pair ends
    b_single_end, files = is_single_end(args.infile1, args.infile2)

    # Get number of reads initially, then log
    num_reads_init = map(get_num_reads, files)
    lenfiles = len(files)
    msg = "Initial number of reads:\n"
    for i in range(len(files)):
        msg = msg + files[i] + ": " + str(num_reads_init[i]) + "\n"
    print(msg)

    # Log file. Create a new one the first time, then keep appending to it.
    logfile = args.output_prefix + ".log"
    with open(logfile, "w") as f:
        f.write(msg)

    print("Running Trimmomatic...")

    trim(files, trimlen = args.trimlen, trim_path = args.trim_path, single_end =
            b_single_end, prefix = args.output_prefix, mem = args.max_mem)
    
    # TODO: run part of the pipeline (if you have the output, either overwrite
    # or do the next step)

    print("Finished running Trimmomatic. Checking output files exist... ")

    # check that Trimmomatic's output files exist
    b_continue, outputs, bmt_inputs = checktrim_output(args.output_prefix, 
            b_single_end)

    msg = "Number of reads after trimming:\n"
    for output in outputs:
        msg = msg + output + ": " + str(get_num_reads(output)) + "\n"
    print(msg)

    with open(logfile, "a") as f:
        f.write(msg)

    if not b_continue:
        print("Trimmomatic produced no non-empty files.")
        print("Terminating the pipeline...")
        sys.exit(1)

    # make temporary directory for BMTagger files
    tempdir = args.output_prefix + "_temp"
    try:
        os.makedirs(tempdir)
    except OSError:
        if not os.path.isdir(tempdir):
            print("Cannot make the BMTagger temporary directory")
            raise

    # TODO: Add parallelization. Use command line utility 'split' to split the
    # files, fork a thread for each file, and run BMTagger in each thread.

    # Start tagging
    out_files = []

    if b_single_end:
        tag(infile = bmt_inputs[0], db_prefix = args.reference_db, 
                bmtagger_path = args.bmtagger_path, single_end = True, prefix =
                args.output_prefix, remove = args.extract, temp_dir = tempdir)

        # Get the proper output file name for logging purposes
        if args.extract:
            out_files.append(args.output_prefix + BMTAGGER_SE_ENDING)
        else:
            out_files.append(args.output_prefix)

    else:
        for inp in bmt_inputs:
            if len(inp) == 2:
                # Run paired end BMTagger
                if not args.extract:
                    out_prefix = args.output_prefix + "_pe"
                tag(infile = inp, db_prefix = args.reference_db, bmtagger_path =
                        args.bmtagger_path, single_end = False, prefix =
                        out_prefix, remove = args.extract, temp_dir = tempdir)

                # Get the proper output file name for logging purposes
                if args.extract:
                    for ending in BMTAGGER_PE_ENDINGS:
                        out_files.append(out_prefix + ending)
                else:
                    out_files.append(args.output_prefix)

            else:
                # Run single end BMTagger
                if not args.extract:
                    out_prefix = args.output_prefix + inp[0] + "_se"
                tag(infile = inp, db_prefix = args.reference_db, bmtagger_path =
                        args.bmtagger_path, single_end = True, prefix =
                        out_prefix, remove = args.extract, temp_dir = tempdir)

                if args.extract:
                    out_files.append(out_prefix + BMTAGGER_SE_ENDING)
                else:
                    out_files.append(out_prefix)

    print("Finished running BMTagger.")
    print(out_files)

    msg = "Number of reads after tagging:\n"
    # Calculate the number of reads remaining
    percent_reads_left = None
    for i in range(len(out_files)):
        num_reads_orig = get_num_reads(out_files[i])
        try:
            num_reads = float(num_reads_orig)
            if b_single_end:
                percent_reads_left = num_reads/num_reads_init[0]

            # if --extract is set, look for one of the paired end files.
            # Otherwise, look for the file containing the list of reads
            else:
                if args.extract and (out_files[i] == args.output_prefix + "_pe" +
                        BMTAGGER_PE_ENDINGS[0]):
                    percent_reads_left = num_reads/num_reads_init[0]
                elif (not args.extract) and (out_files[i] == args.output_prefix):
                    percent_reads_left = (num_reads_init[0] -
                            (num_reads*4))/num_reads_init[0]
        except TypeError:
            pass

        msg = msg + out_files[i] + ": " + str(num_reads_orig) + "\n"

    msg2 = "Proportion of reads that survived: " + str(percent_reads_left)
    print(msg)
    print(msg2)

    with open(logfile, "a") as f:
        f.write(msg)
        f.write(msg2)

    print("Removing temporary files...")
    for output in outputs:
        os.remove(output)
    shutil.rmtree(tempdir)
    print("Done!")

if __name__ == '__main__':
    main()
