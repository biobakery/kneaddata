import os
import re
import sys
import logging
import subprocess
from glob import glob
from itertools import tee, izip_longest

from . import divvy_threads, mktempfifo, process_return
here = os.path.dirname(os.path.realpath(__file__))


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
            output_con_fastq, threads, bowtie2_args, bowtie2_path="bowtie2"):
    args = [bowtie2_path, 
            "-x", index_str,
            "-U", input_fastq,
            "--un", output_clean_fastq,
            "--al", output_con_fastq,
            "-S", os.devnull,
            "--threads", str(threads)] + bowtie2_args
    logging.debug("Running bowtie2 with arguments %s", args)
    return subprocess.Popen(args, stderr=subprocess.PIPE,
                            stdout=subprocess.PIPE)


def run_tandem(in_fastq, out, match=2, mismatch=7, delta=7, pm=80, pi=10,
        minscore=50, maxperiod=500, generate_fastq=True, mask=False, html=False,
        trf_path="trf"):
    tandem_cmd = ["python", os.path.join(here, "tandem.py"),
                  in_fastq, out, 
                  "--match", str(match),
                  "--mismatch", str(mismatch),
                  "--delta", str(delta),
                  "--pm", str(pm),
                  "--pi", str(pi),
                  "--minscore", str(minscore),
                  "--maxperiod", str(maxperiod),
                  "--trf-path", trf_path]

    if mask:
        tandem_cmd += ["--mask"]
    if not generate_fastq:
        tandem_cmd += ["--no-generate-fastq"]
    if html:
        tandem_cmd += ["--html"]

    logging.debug("Running tandem.py with: %s" %(" ".join(tandem_cmd)))
    return subprocess.Popen(tandem_cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)


def decontaminate_reads(in_fname, index_strs, output_prefix,
                        output_dir, filter_args_list, filter_jar_path,
                        trim_threads, bowtie_threads, bowtie2_args, 
                        bowtie2_path, trf, match, mismatch, delta, pm, pi,
                        minscore, maxperiod, generate_fastq, mask, html,
                        trf_path):
    
    tmpfilebases = ['filter']+map(os.path.basename, index_strs[:-1])
    if trf:
        tmpfilebases += [os.path.basename(index_strs[-1])]

    with mktempfifo(tmpfilebases) as filenames:
        clean_file = os.path.join(output_dir, output_prefix+".fastq")
        filter_proc = trimmomatic(in_fname, filenames[0],
                                  filter_args_list, filter_jar_path,
                                  threads=trim_threads)
        staggered = sliding_window(filenames, 2)

        def _procs():
            for (in_cur, in_next), index_str in izip_longest(staggered,
                                                             index_strs):
                if index_str != None:
                    contam_name = "{}_{}_contam.fastq".format(
                        output_prefix, os.path.basename(index_str))
                    contam_file = os.path.join(output_dir, contam_name)
                if trf:
                    if in_next != None:
                        yield bowtie2(index_str, in_cur, in_next, contam_file,
                                    bowtie_threads, bowtie2_args, 
                                    bowtie2_path=bowtie2_path)
                    else:
                        assert(index_str == None)
                        yield run_tandem(in_cur, clean_file, match, mismatch,
                                delta, pm, pi, minscore, maxperiod,
                                generate_fastq, mask, html, trf_path)
                else:
                    in_next = in_next or clean_file
                    yield bowtie2(index_str, in_cur, in_next, contam_file,
                                bowtie_threads, bowtie2_args, 
                                bowtie2_path=bowtie2_path)


        procs = [filter_proc]+list(_procs())
        names = ["trimmomatic"] + ["bowtie2 (%s)"%(idx) for idx in index_strs]

        if trf:
            names += ["tandem.py"]

        assert(len(procs) == len(names))

        for proc, name in zip(procs, names):
            stdout, stderr = proc.communicate()
            retcode = proc.returncode
            process_return(name, retcode, stdout, stderr)


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

    if not args.bowtie2_path:
        args.bowtie2_path = "bowtie2"

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
                        bowtie_threads, bowtie2_args=args.bowtie2_args, 
                        bowtie2_path=args.bowtie2_path,
                        trf=args.trf, match=args.match, mismatch=args.mismatch,
                        delta=args.delta, pm=args.pm, pi=args.pi,
                        minscore=args.minscore, maxperiod=args.maxperiod,
                        generate_fastq=args.no_generate_fastq, mask=args.mask,
                        html=args.html, trf_path=args.trf_path)
