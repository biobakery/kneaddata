import os
import re
import sys
import logging
import tempfile
import subprocess
from glob import glob
from functools import partial
from contextlib import contextmanager
from itertools import tee, izip_longest

from . import divvy_threads


def rmext(name_str, all=False):
    """removes file extensions 

    :keyword all: Boolean; removes all extensions if True, else just
    the outside one

    """

    _match = lambda name_str: re.match(r'(.+)(\..*)', name_str)
    path, name_str = os.path.split(name_str)
    match = _match(name_str)
    while match:
        name_str = match.group(1)
        match = _match(name_str)
        if not all:
            break

    return os.path.join(path, name_str)
    

def sliding_window(it, l, fill=None):
    args = tee(it, l)
    # advance each iterator as many steps as its rank-1
    # thus, advance the 1st iterator none, 2nd iterator once. etc.
    for i in range(len(args)-1, 0, -1):
        for iter_ in args[i:]:
            next(iter_, None)
    return izip_longest(*args, fillvalue=fill)


def trimmomatic(fastq_in, fastq_out, filter_args_list, jar_path,
                maxmem="500m", threads=1):
    args = ["java", "-Xmx"+maxmem, "-d64",
            "-jar", jar_path,
            "SE", "-threads", str(threads),
            fastq_in, fastq_out]
    args += filter_args_list
    logging.debug("Running trimmomatic with arguments %s", args)
    return subprocess.Popen(args, stderr=subprocess.PIPE,
                            stdout=subprocess.PIPE)


def bowtie2(index_str, input_fastq, output_clean_fastq,
            output_con_fastq, threads):
    args = ["bowtie2", "--very-sensitive",
            "-x", index_str,
            "-U", input_fastq,
            "--un", output_clean_fastq,
            "--al", output_con_fastq,
            "-S", os.devnull,
            "--threads", str(threads),
            '--quiet']
    logging.debug("Running bowtie2 with arguments %s", args)
    return subprocess.Popen(args, stderr=subprocess.PIPE,
                            stdout=subprocess.PIPE)



@contextmanager
def mktempfifo(names=("a",)):
    tmpdir = tempfile.mkdtemp()
    names = map(partial(os.path.join, tmpdir), names)
    map(os.mkfifo, names)
    yield names
    map(os.remove, names)
    os.rmdir(tmpdir)


def decontaminate_reads(in_fname, index_strs, output_prefix,
                        output_dir, filter_args_list, filter_jar_path,
                        trim_threads, bowtie_threads):
    tmpfilebases = ['filter']+map(os.path.basename, index_strs[:-1])

    with mktempfifo(tmpfilebases) as filenames:
        clean_file = os.path.join(output_dir, output_prefix+".fastq")
        filter_proc = trimmomatic(in_fname, filenames[0],
                                  filter_args_list, filter_jar_path,
                                  threads=trim_threads)
        staggered = sliding_window(filenames, 2)

        def _procs():
            for (in_cur, in_next), index_str in zip(staggered, index_strs):
                contam_name = "{}_{}_contam.fastq".format(
                    output_prefix, os.path.basename(index_str))
                contam_file = os.path.join(output_dir, contam_name)
                in_next = in_next or clean_file
                yield bowtie2(index_str, in_cur, in_next, contam_file,
                              bowtie_threads)

        procs = [filter_proc]+list(_procs())
        names = ["trimmomatic"] + [ "bowtie2 (%s)"%(idx) for idx in index_strs ]
        for proc, name in zip(procs, names):
            stdout, stderr = proc.communicate()
            retcode = proc.returncode
            if retcode:
                log = logging.critical
                log("%s exited with exit status %d", name, retcode)
            else:
                log = logging.debug
            if stdout:
                log("%s stdout: %s", name, stdout)
            if stderr:
                log("%s stderr: %s", name, stderr)

            if retcode:
                sys.exit(retcode)


def check_args(args):
    if not args.output_prefix:
        args.output_prefix = rmext(os.path.basename(args.infile1), all=True)

    if args.infile2:
        logging.critical("memory heavy strategy only supports single-end reads")
        sys.exit(1)

    if not os.path.exists(args.output_dir):
        logging.warning("Output dir %s doesn't exist; creating",
                        args.output_dir)
        os.mkdir(args.output_dir)

    for db_base in args.reference_db:
        if not glob(db_base+"*"):
            logging.critical("Unable to find database `%s'", db_base)
            sys.exit(1)

    return args


def memory_heavy(args):
    args = check_args(args)
    trim_threads, bowtie_threads = divvy_threads(args)
    decontaminate_reads(args.infile1, args.reference_db,
                        args.output_prefix, args.output_dir,
                        args.trim_args, args.trim_path, trim_threads,
                        bowtie_threads)
    

    
