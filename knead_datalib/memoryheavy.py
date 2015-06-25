import os
import re
import sys
import logging
import subprocess
from glob import glob
from itertools import tee, izip_longest

import tandem
from . import divvy_threads, mktempfifo, process_return


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
            output_con_fastq, threads, bowtie2_path="bowtie2"):
    #args = [bowtie2_path, "--very-sensitive",
    #        "-x", index_str,
    #        "-U", input_fastq,
    #        "--un", output_clean_fastq,
    #        "--al", output_con_fastq,
    #        "-S", os.devnull,
    #        "--threads", str(threads),
    #        '--quiet']
    args = [bowtie2_path, "--very-sensitive",
            "-x", index_str,
            "-U", input_fastq,
            "--un", output_clean_fastq,
            "--al", output_con_fastq,
            "-S", os.devnull,
            "--threads", str(threads)]
    logging.debug("Running bowtie2 with arguments %s", args)
    return subprocess.Popen(args, stderr=subprocess.PIPE,
                            stdout=subprocess.PIPE)


def decontaminate_reads(in_fname, index_strs, output_prefix,
                        output_dir, filter_args_list, filter_jar_path,
                        trim_threads, bowtie_threads, bowtie2_path="bowtie2",
                        trf=True, match=2, mismatch=7, delta=7, pm=80, pi=10,
                        minscore=5, maxperiod=500, generate_fastq=True,
                        dat=False, mask=False, html=False, trf_path="trf"):

    trf_args = map(str, [match, mismatch, delta, pm, pi, minscore, maxperiod])
    fasta_fname = os.path.basename(in_fname) + ".fasta"
    trf_out_fname = fasta_fname + "." + ".".join(trf_args) + ".stream.dat"
    mask_fname = None
    if mask:
        mask_fname = fasta_fname + "." + ".".join(trf_args) + ".mask"

    tmpfilebases = (['filter']+map(os.path.basename, index_strs) +
            [fasta_fname, trf_out_fname])

    with mktempfifo(tmpfilebases) as filenames:
        # final bowtie2 output, if trf
        last_bowtie2 = filenames[-3]

        fasta_file = filenames[-2]
        # trf output, to be converted to final output
        trf_out_file = filenames[-1]
        # final output
        clean_file = os.path.join(output_dir, output_prefix+".fastq")
        filter_proc = trimmomatic(in_fname, filenames[0],
                                  filter_args_list, filter_jar_path,
                                  threads=trim_threads)
        staggered = sliding_window(filenames[:-2], 2)
        #print(list(staggered))
        #sys.exit(1)

        def _bowtie2_procs():
            for (in_cur, in_next), index_str in zip(staggered, index_strs):
                contam_name = "{}_{}_contam.fastq".format(
                    output_prefix, os.path.basename(index_str))
                contam_file = os.path.join(output_dir, contam_name)
                if in_next == last_bowtie2:
                    if not trf: 
                        in_next = clean_file
                #print("\nPRINTED STUFF HERE\n")
                #print(in_cur, in_next)
                yield bowtie2(index_str, in_cur, in_next, contam_file,
                              bowtie_threads, bowtie2_path=bowtie2_path)

        procs = [filter_proc] + list(_bowtie2_procs())
        names = ["trimmomatic"] + [ "bowtie2 (%s)"%(idx) for idx in index_strs ] 
        if not trf:
            print("RETURNING SOME STUFF")
            for (proc, name) in zip(procs, names):
                stdout, stderr = proc.communicate()
                retcode = proc.returncode
                process_return(name, retcode, stdout, stderr)
            return

        fastq_to_fasta_proc = tandem._fastq_to_fasta(last_bowtie2, fasta_file)
        # need to make fifo otherwise converting won't work
        #if mask:
        #    os.mkfifo(mask_fname)

        converter_procs = tandem._convert(last_bowtie2, trf_out_file,
                clean_file, mask_fname, generate_fastq, mask)

        trf_out_fp = open(trf_out_file, "w")
        try:
            trf_proc = tandem._trf(fasta_file, trf_out_fp, match, mismatch,
                    delta, pm, pi, minscore, maxperiod, dat, mask, html,
                    trf_path)

            procs = procs + [fastq_to_fasta_proc, trf_proc]
            names = names + ["fastq_to_fasta", "trf"]
            print(procs)
            for (proc, name) in zip(procs, names):
                stdout, stderr = proc.communicate()
                retcode = proc.returncode
                process_return(name, retcode, stdout, stderr)
        finally: 
            trf_out_fp.close()

        # wait for converter procs
        converter_names = ["converting %s -> %s" %(trf_out_file, [clean_file,
            clean_file + ".mask"][i]) for (i,_) in enumerate(converter_procs)]
        for proc, name in zip(converter_procs, converter_names):
            stdout, stderr = proc.communicate()
            retcode = proc.returncode
            process_return(name, retcode, stdout, stderr)

        # remove mask file
        #if mask:
        #    os.remove(mask_fname)


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
                        bowtie_threads, bowtie2_path=args.bowtie2_path,
                        trf = args.trf, match = args.match, mismatch =
                        args.mismatch, delta = args.delta, pm = args.pm, pi =
                        args.pi, minscore = args.minscore, maxperiod =
                        args.maxperiod, generate_fastq = args.no_generate_fastq,
                        mask = args.mask, html = args.html,
                        trf_path = args.trf_path)
