import os
import sys
import logging
import argparse
import itertools
import subprocess

# import stuff from the __init__.py in this directory
from __init__ import parse_positive_int, mktempfifo, try_create_dir, process_return
here = os.path.dirname(os.path.realpath(__file__))

def _fastq_to_fasta(fastq, fasta):
    '''
    Converts fastq to fasta using fastq_to_fasta.py
    fastq: Input fastq filename
    fasta: Output fasta filename
    '''
    fastq_to_fasta_cmd = ["python", 
                          os.path.join(here, "fastq_to_fasta.py"),
                          fastq, 
                          "--fasta", fasta]
    logging.debug("Running fastq_to_fasta with %s", fastq_to_fasta_cmd)
    fastq_to_fasta = subprocess.Popen(fastq_to_fasta_cmd,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return(fastq_to_fasta)


def _trf(fasta, outfp, match=2, mismatch=7, delta=7, pm=80, pi=10, minscore=50,
        maxperiod=500, mask=False, html=False, trf_path="trf"):
    '''
    Generate the TRF command and run it using subprocess. TRF is run such that
    the reads identified as containing tandem repeats are piped from stdout to
    another program

    fasta: input fasta filename for TRF
    outfp: A file pointer that specifies where to write the TRF stdout
    trf_path: path to trf. If not specified, defaults to "trf"

    match, mismatch, delta, pm, pi, minscore, and maxperiod are all TRF
    parameters. 

    mask: Tells TRF whether to generate a mask file.
    html: Tells TRF whether to generate a HTML output
    '''

    trf_args = map(str, [match, mismatch, delta, pm, pi, minscore, maxperiod])
    trf_cmd = [trf_path, fasta] + trf_args + ["-ngs"]

    if mask:
        trf_cmd += ["-m"]
    if not html:
        trf_cmd += ["-h"]

    logging.debug("Running trf with %s", trf_cmd)
    trf = subprocess.Popen(trf_cmd, stdout=outfp, stderr=subprocess.PIPE)
    return(trf)


def _generate_fastq(fastq, trf_output, out):
    cmd = ["python", os.path.join(here, "tandem_process.py"), 
            fastq, trf_output, out]
    logging.debug("Running generate_fastq with %s", cmd)
    proc = subprocess.Popen(cmd)
    return(proc)


def _generate_mask(fastq, trf_output, out):
    cmd = ["python", os.path.join(here, "tandem_process.py"), 
            fastq, trf_output, out, "--mask"]
    logging.debug("Running generate_mask with %s", cmd)
    proc = subprocess.Popen(cmd)
    return(proc)
    

def _convert(fastq, trf_output, out, mask_fname, generate_fastq, mask):
    converter_procs = []
    converter_names = []

    if generate_fastq: 
        converter_procs.append(_generate_fastq(fastq, trf_output, out))

    if mask:
        mask_out = out
        if generate_fastq:
            mask_out = out + ".mask"
        assert(mask_fname != None)
        converter_procs.append(_generate_mask(fastq, mask_fname,  mask_out))

    return(converter_procs)


def handle_cli():
    parser = argparse.ArgumentParser(description="Remove tandem repeats.")
    parser.add_argument("fastqs", nargs="+", help="input fastqs")
    parser.add_argument("-o", "--outputs", nargs="+", required=True, 
            help="output prefixes")

    parser.add_argument(
            "--trf-path",
            default="trf",
            help="Path to TRF executable if not found in $PATH")
    parser.add_argument(
            "--match", type=parse_positive_int,
            default=2, 
            help="TRF matching weight")
    parser.add_argument(
            "--mismatch", type=parse_positive_int,
            default=7, 
            help="TRF mismatching penalty")
    parser.add_argument(
            "--delta", type=parse_positive_int,
            default=7, 
            help="TRF indel penalty")
    parser.add_argument(
            "--pm", type=parse_positive_int,
            default=80,
            help="TRF match probability (whole number)")
    parser.add_argument(
            "--pi", type=parse_positive_int,
            default=10,
            help="TRF indel probability (whole number)")
    parser.add_argument(
            "--minscore", type=parse_positive_int,
            default=50, 
            help="TRF minimum alignment score to report")
    parser.add_argument(
            "--maxperiod", type=parse_positive_int,
            default=500, 
            help="TRF maximum period size to report")
    parser.add_argument(
            "--no-generate-fastq",
            default=True, action="store_false", 
            help="If switched on, don't generate fastq")
    parser.add_argument(
            "--mask",
            default=False, action="store_true",
            help="Generate mask file")
    parser.add_argument(
            "--html",
            default=False, action="store_true",
            help="Generate html file")
    
    args = parser.parse_args()

    if (not args.no_generate_fastq) and (not args.mask) and args.trf:
        parser.error("\nYou cannot set the --no-generate-fastq flag without"
        " the --mask flag. Exiting...\n")

    assert(len(args.fastqs) == len(args.outputs))

    return args


def main():
    args = handle_cli()

    '''
    # TEMP LOGGING
    fmt = "%(asctime)s %(levelname)s: %(message)s"
    logger = logging.getLogger()
    logger.setLevel(getattr(logging, "DEBUG"))
    logging.basicConfig(format=fmt)
    '''

    nfiles = len(args.fastqs)

    fasta_fnames = [os.path.basename(fastq) + ".fasta" for fastq in args.fastqs]
    # file names for trf mask (trf default provided, can't change it)
    trf_args = map(str, [args.match, args.mismatch, args.delta, args.pm,
        args.pi, args.minscore, args.maxperiod])
    trf_out_fnames = [f + "." + ".".join(trf_args) + ".stream.dat" for f in
            fasta_fnames]
    mask_fnames = [None] * nfiles
    if args.mask:
        mask_fnames = [f + "." + ".".join(trf_args) + ".mask" for f in
                fasta_fnames]

    with mktempfifo(fasta_fnames + trf_out_fnames) as filenames:
        fasta_outs = filenames[:nfiles]
        trf_outs = filenames[nfiles:]
        if not args.no_generate_fastq:
            trf_outs = [os.devnull for f in args.fastqs]

        # process names
        fastq_to_fasta_names = []
        converter_names = []
        trf_names = []

        # subprocess.Popen instances
        fastq_to_fasta_procs = []
        converter_procs = []
        trf_procs = []

        for (fastq_in, fasta_out) in zip(args.fastqs, fasta_outs):
            proc = _fastq_to_fasta(fastq_in, fasta_out)
            fastq_to_fasta_names.append("fastq_to_fasta %s -> %s" %(fastq_in,
                                                                    fasta_out))
            fastq_to_fasta_procs.append(proc)

        # must start converter_procs before opening the file handle, otherwise
        # will hang
        # must make a fifo for the mask output otherwise converter will just
        # read through the whole file
        if args.mask:
            for mask_fname in mask_fnames:
                logging.debug("Making fifo %s" %mask_fname)
                os.mkfifo(mask_fname)

        try:
            for (fastq, trf_out, out, mask_fname) in zip(args.fastqs, trf_outs,
                    args.outputs, mask_fnames):
                converter_proc_grp = _convert(fastq, trf_out, out, mask_fname,
                        args.no_generate_fastq, args.mask)
                for (i, c) in enumerate(converter_proc_grp):
                    converter_procs.append(c)
                    cur_out = [out, out + ".mask"][i]
                    converter_names.append("converting %s -> %s" %(trf_out,
                                                                   cur_out))

            trf_out_fps = [open(t, "w") for t in trf_outs]
            # steps: 
            # 1. wait for the fastq_to_fasta and trf procs. 
            # 2. close the trf_out_fp file handle
            # 3. wait for the converter_procs
            # 4. if (not debug) and mask: remove the mask_fname
            try:
                for (fasta, trf_out_fp) in zip(fasta_outs, trf_out_fps):
                    trf_proc = _trf(fasta, trf_out_fp, args.match,
                            args.mismatch, args.delta, args.pm, args.pi,
                            args.minscore, args.maxperiod, args.mask, args.html,
                            args.trf_path)
                    trf_names.append("trf on %s" %fasta)
                    trf_procs.append(trf_proc)
                for (name, proc) in zip(itertools.chain(fastq_to_fasta_names,
                    trf_names), itertools.chain(fastq_to_fasta_procs,
                        trf_procs)):
                    stdout, stderr = proc.communicate()
                    retcode = proc.returncode
                    process_return(name, retcode, stdout, stderr)
            finally:
                for fp in trf_out_fps:
                    fp.close()

            for (name, proc) in zip(converter_names, converter_procs):
                stdout, stderr = proc.communicate()
                retcode = proc.returncode
                process_return(name, retcode, stdout, stderr)

        finally:
            # remove mask files
            if args.mask:
                for mask_fname in mask_fnames:
                    logging.debug("Removing %s" %mask_fname)
                    os.remove(mask_fname)

    return 100
        

if __name__ == "__main__":
    main()
