'''
kneaddata.py
Author: Andy Shi

Pipeline for processing metagenomics sequencing data
'''
import argparse
import subprocess
import shlex
import os
import sys

def trim(infile, trimlen, trim_path, single_end, prefix):
    '''
    input: 
        infile:     input fastq file list (either length 1 or length 2)
                    length 1: single end
                    length 2: paired end
    optional input: 
        trimlen:    length to trim
        prefix:     prefix for outputs
    output: At most 5 files
        (1) prefix.trimmed.1.fastq: trimmed first pair-end file
        (2) prefix.trimmed.2.fastq: trimmed second pair-end file
        (3) prefix.trimmed.single.1.fastq: trimmed sequences from the first file
        that lost their partner from the second file
        (4) prefix.trimmed.single.2.fastq: trimmed sequences from the second
        file that lost their partner from the first file
        (5) prefix.log: log file with stdout output from Trimmomatic
    Uses trimmomatic to trim reads to [arg] base pairs long
    '''

    # check that we have the right number of input fastqs
    assert(len(infile) == 2 or len(infile) == 1)

    trim_arg = ""
    if single_end:
        trim_arg = str(trim_path + " SE -phred33 " + infile[0] + " " + prefix
                + ".trimmed.fastq " + "MINLEN:" + str(trimlen) + " &> " + prefix
                + ".log")
    else:
        trim_arg = str(trim_path + " PE -phred33 " + infile[0] + " " +
                infile[1] + " " + prefix + ".trimmed.1.fastq " + prefix +
                ".trimmed.single.1.fastq " + prefix + ".trimmed.2.fastq " +
                prefix + ".trimmed.single.2.fastq " + "MINLEN:" + str(trimlen))

    cmd = "java -Xmx500m -jar " + trim_arg
    print("Trimmomatic command that will be run: " + cmd)
    call = shlex.split(cmd)
    subprocess.call(call)
    return 0


def tag(infile, db_prefix, bmtagger_path, single_end, prefix):
    '''
    input:
        infile: a list of length 1 or 2 (for single and paired ends,
        respectively) 
    output:
    Uses BMTagger to tag and potentially remove unwanted reads
    '''
    # check inputs
    db_len = len(db_prefix)
    assert (db_len > 0)

    bmt_args = ["" for i in xrange(db_len)]
    # build arguments
    for i in xrange(db_len):
        db = db_prefix[i]
        if single_end:
            bmt_args[i] = str(bmtagger_path + " -q 1 -1 " + infile[0] + 
                    " -b " + db + ".bitmask -x " + db + 
                    ".srprism -T ./temp_dir -o " + prefix + ".out") 
        else:
            bmt_args[i] = str(bmtagger_path + " -q 1 -1 " + infile[0] + 
                    " -2 " + infile[1] + " -b " + db + ".bitmask -x " + db + 
                    ".srprism -T ./temp_dir -o " + prefix + ".out")
    
    for arg in bmt_args:
        # Run all the BMTagger instances 
        print("BMTagger command to be run: " + arg)
        call = shlex.split(arg)
        print(call)
        subprocess.call(call)
    return 0


def checkfile(fname):
    '''
    input:
        fname: the file name. A string.
    output:
        1: the file exists
        0: the file does not exist
        -1: the file exists but is an empty file (size = 0 bytes)
    '''
    try:
        if os.stat(fname).st_size > 0:
            return 1
        else:
            return -1
    except OSError:
        return 0

def main():
    # parse command line arguments
    # note: argparse converts dashes '-' in arguement prefixes to underscores
    # '_' 
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
    parser.add_argument("-t", "--trim-path", help="path to Trimmomatic",
            required = True)
    parser.add_argument("-b", "--bmtagger-path", help="path to BMTagger",
            required = True)
    parser.add_argument("-S", "--slurm", help="Running in a slurm environment",
            action = "store_true")

    args = parser.parse_args()

    # check inputs
    # deal with missing prefix
    if not args.output_prefix:
        args.output_prefix = args.infile1

    # check for the existence of required files/paths
    paths = [args.infile1, args.infile2, args.trim_path, args.bmtagger_path]
    for path in paths:
        if path != None and not os.path.exists(path):
            print("Could not find file " + str(path))
            print("Aborting...")
            sys.exit(2)
    for db in args.reference_db:
        if not os.path.exists(db):
            print("Could not find BMTagger database " + db)
            print("Aborting...")
            sys.exit(2)

    # determine single-ended or pair ends
    b_single_end = True
    files = [args.infile1]
    if args.infile2:
        files.append(args.infile2)
        b_single_end = False

    print("Running Trimmomatic...")

    trim(files, trimlen = args.trimlen, trim_path = args.trim_path, single_end =
            b_single_end, prefix = args.output_prefix)
    
    print("Finished running Trimmomatic. Checking output files exist... ")

    # check that Trimmomatic's output files exist

    # TODO: Take better care to make the temporary directory for BMTagger files
    outputs = []
    bmt_inputs = []
    if b_single_end:
        outputs.append(str(args.output_prefix + "trimmed.fastq"))
        if not os.path.exists(outputs[0]):
            print("Could not find file " + output)
            print("Trimmomatic failed. Exiting...")
            sys.exit(1)
        bmt_inputs = [outputs]
    else:
        outputs = [args.output_prefix + ".trimmed." for i in xrange(4)]
        outputs[0] = outputs[0] + "1.fastq"
        outputs[1] = outputs[1] + "2.fastq"
        outputs[2] = outputs[2] + "single.1.fastq"
        outputs[3] = outputs[3] + "single.2.fastq"
        checks = map(checkfile, outputs)
        for i in xrange(4):
            if checks[i] == 0:
                print("Could not find file " + output[i])
                print("Trimmomatic failed. Exiting...")

        if checks[0] == 1 and checks[1] == 1:
            bmt_inputs.append([outputs[0], outputs[1]])
        elif checks[0] == 1:
            bmt_inputs.append([outputs[0]])
        elif checks[1] == 1:
            bmt_inputs.append([outputs[1]])
        for i in [2,3]:
            if checks[i] == 1:
                bmt_inputs.append([outputs[i]])

    if b_single_end:
        tag(infile = bmt_inputs, db_prefix = args.reference_db, bmtagger_path =
                args.bmtagger_path, single_end = True, prefix =
                args.output_prefix)
    else:
        for input in bmt_inputs:
            if len(input) == 2:
                tag(infile = input, db_prefix = args.reference_db, bmtagger_path
                    = args.bmtagger_path, single_end = False, prefix =
                    args.output_prefix + "_pe")
            else:
                tag(infile = input, db_prefix = args.reference_db,
                    bmtagger_path = args.bmtagger_path, single_end = True,
                    prefix = args.output_prefix + "_se_" + input)

    print("Finished running BMTagger.")

if __name__ == '__main__':
    main()
