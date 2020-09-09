#!/usr/bin/env python

'''
generate_db.py
Author: Andy Shi

Helper script to generate databases for the knead_data.py pipeline
Requires the BMTagger suite and NCBI BLAST.
'''

import argparse
import subprocess
import shlex
import os
import sys
import threading

def exists(s_cmd):
    '''
    Prints out error message with appropriate command name substituted in
    '''
    out = str("The output files for " + s_cmd + " already exist! Moving on...")
    print(out)
    return  


class run_proc_thread(threading.Thread):
    def __init__(self, command, name, logdir="log"):
        threading.Thread.__init__(self)

        # command to call
        self.command = command

        # name of the command (ie bmtool or srprism)
        self.name = name

        # directory to put log files
        self.logdir = logdir

    def run(self):
        print(self.command)
        print(shlex.split(self.command))
        p = subprocess.Popen(shlex.split(self.command), stdout =
                subprocess.PIPE, stderr = subprocess.PIPE)
        # returns (stdout, stderr)
        p_outputs = p.communicate()

        endings = [".out", ".err"]
        logfiles = [os.path.join(self.logdir, self.name + ending) for ending in
                endings]

        for (logfile,output) in zip(logfiles, p_outputs):
            with open(logfile, "w") as f:
                f.write(output.decode("utf-8"))
        

def run_proc(s_cmd, name, logdir):
    '''
    Helper function to run a command in a new thread.
    input:     s_cmd (string), which is the command to be run
    returns:   a (pointer to) the thread which is running the command
    '''
    print("Command to be run:")
    print(s_cmd)
    t = run_proc_thread(s_cmd, name, logdir)
    t.start()
    return t


def main():
    # parse command line arguments
    # note: argparse converts dashes '-' in argument prefixes to underscores '_'
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help="input FASTA file")
    parser.add_argument("-o", "--output-prefix", 
        help="prefix for all output files")

    # assumes that bmtool, srprism, and makeblastdb are all in the user's $PATH.
    # Otherwise you can specify their locations.
    parser.add_argument("-b", "--bmtool-path", default="bmtool",
            help="path to bmtool executable")
    parser.add_argument("-s", "--srprism-path", default="srprism",
            help="path to srprism executable")
    parser.add_argument("-m", "--makeblastdb-path", default="makeblastdb",
            help="path to makeblastdb executable")
    parser.add_argument("-l", "--logdir", default="log",
            help = "location to store log files")

    args = parser.parse_args()

    # Take care of missing prefix
    if not args.output_prefix:
        args.output_prefix = args.fasta

    # TODO: Consider factoring out the checking into a util.py file. 
    # TODO: If the user supplies a path, check it. Manually set the defaults

    # check if input files exist
    #inputs = [args.fasta, args.bmtool_path, args.srprism_path,
    #        args.makeblastdb_path]
    inputs = [args.fasta]
    for inp in inputs:
        if not os.path.exists(inp):
            print("Could not find file " + inp)
            print("Aborting... ")
            sys.exit(1)

    # make the log directory
    try:
        os.makedirs(args.logdir)
    except OSError:
        if os.path.isdir(args.logdir):
            print("Log directory exists! Overwriting its contents...")
        else:
            raise

    # before running each command, check that output files don't already exist.
    
    threads = []

    # bmtool
    if not os.path.exists(args.output_prefix + ".bitmask"):
        print("Running bmtool!")
        cmd = str(args.bmtool_path + " -d " + args.fasta + " -o " +
            args.output_prefix + ".bitmask -A 0 -z -w 18")
        print("Command to be run:")
        print(cmd)
        t = run_proc(cmd, "bmtool", args.logdir)
        threads.append(t)
    else:
        exists("bmtool")

    # srprism 
    srprism_ext = [".amp", ".idx", ".imp", ".map", ".pmp", ".rmp", ".ss",
            ".ssa", ".ssd"]
    srprism_files = map(lambda x: str(args.output_prefix + ".srprism" + x),
            srprism_ext)
    print("Checking for the following files:")
    print(srprism_files)

    # run if >= 1 srprism output file does not exist
    b_srprism_run = False

    for f in srprism_files:
        if not os.path.exists(f):
            b_srprism_run = True
            break

    if b_srprism_run:
        print("Running srprism!")
        cmd = str(args.srprism_path + " mkindex -i " + args.fasta + " -o " +
            args.output_prefix + ".srprism -M 7168")
        print("Command to be run:")
        print(cmd)
        t = run_proc(cmd, "srprism", args.logdir)
        threads.append(t)
    else:
        exists("srprism")


    # makeblastdb
    blastdb_ext = [".nhr", ".nin", ".nsq"]

    blastdb_files = map(lambda x: str(args.output_prefix + x), blastdb_ext)
    print("Checking for the following files:")
    print(blastdb_files)

    # run if >= 1 makeblastdb output file does not exist
    b_blastdb_run = False

    for f in blastdb_files:
        if not os.path.exists(f):
            b_blastdb_run = True
            break

    if b_blastdb_run:
        print("Running makeblastdb!")
        cmd = str(args.makeblastdb_path + " -in " + args.fasta + 
            " -dbtype nucl -out " + args.output_prefix)
        print("Command to be run:")
        print(cmd)
        t = run_proc(cmd, "makeblastdb", args.logdir)
        threads.append(t)
    else:
        exists("makeblastdb")
    
    # wait for all threads to finish
    for t in threads:
        t.join()

if __name__ == '__main__':
    main()
