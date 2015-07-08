import os
import sys
import logging
import argparse
import itertools
import subprocess

# import stuff from the __init__.py in this directory
from __init__ import parse_positive_int, mktempfifo, try_create_dir, process_return
here = os.path.dirname(os.path.realpath(__file__))

############## constants ######################
# number of lines per entry in trf output
TANDEM_LENGTH = 2

# number of lines per entry in fastq 
FASTQ_LENGTH = 4
###############################################


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
    #trf_cmd = [trf_path, fasta] + trf_args + ["-ngs"]
    trf_cmd = [trf_path, fasta] + trf_args 

    if mask:
        trf_cmd += ["-m"]
    if not html:
        trf_cmd += ["-h"]

    logging.debug("Running trf with %s", trf_cmd)
    #trf = subprocess.Popen(trf_cmd, stdout=outfp, stderr=subprocess.PIPE)
    trf = subprocess.Popen(trf_cmd, stdout=outfp,
            stderr=subprocess.PIPE, bufsize=1)
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


def run_tandem_process(fastq, fasta, output, tandem, maskfile):
    cmd = ["python", os.path.join(here, "tandem_process.py"), fastq, fasta,
            output]

    if tandem != None:
        cmd += ["--tandem", tandem]
    if maskfile != None:
        cmd += ["--maskfile", maskfile]

    logging.debug("Running tandem_process.py with %s" %str(cmd))
    #return(subprocess.Popen(cmd, stdout=subprocess.PIPE,
    #    stderr=subprocess.PIPE))
    return(subprocess.Popen(cmd))


def run_tandem(fastq, output, match=2, mismatch=7, delta=7, pm=80, pi=10,
        minscore=50, maxperiod=500, generate_fastq=True, mask=False, html=False,
        trf_path="trf"):
    fasta_fname = os.path.basename(fastq) + ".fasta"
    trf_args = map(str, [match, mismatch, delta, pm, pi, minscore, maxperiod])
    trf_out_fname = fasta_fname + "." + ".".join(trf_args) + ".stream.dat"

    mask_fname = None
    if mask:
        mask_fname = fasta_fname + "." + ".".join(trf_args) + ".mask"

    #with mktempfifo((fasta_fname, trf_out_fname)) as filenames:
    with mktempfifo((fasta_fname,)) as filenames:
        if not generate_fastq:
            trf_out = os.devnull
            tandem = None
        else:
            trf_out = trf_out_fname
            tandem = trf_out
            os.mkfifo(tandem)

        trf_out_fp = None
        if mask:
            os.mkfifo(mask_fname)
        try:
            tandem_process_proc = run_tandem_process(fastq, fasta_fname, output,
                    tandem, mask_fname)
            trf_out_fp = open(trf_out, "w")

            trf_proc = _trf(filenames[0], trf_out_fp, match, mismatch, delta,
                    pm, pi, minscore, maxperiod, mask, html, trf_path)
            procs = [trf_proc, tandem_process_proc]
            names = ["trf", "tandem_process"]
            for (proc, name) in zip(procs, names):
                print(name)
                stdout, stderr = proc.communicate()
                retcode = proc.returncode
                process_return(name, retcode, stdout, stderr)

        finally:
            if trf_out_fp != None:
                trf_out_fp.close()
            if generate_fastq:
                os.remove(tandem)
            if mask:
                os.remove(mask_fname)
    

def handle_cli():
    parser = argparse.ArgumentParser(description="Remove tandem repeats.")
    parser.add_argument("fastq", help="input fastqs")
    parser.add_argument("output", help="output prefixes")

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

    #assert(len(args.fastqs) == len(args.outputs))

    return args


def main():
    args = handle_cli()

    # TEMP LOGGING
    fmt = "%(asctime)s %(levelname)s: %(message)s"
    logger = logging.getLogger()
    logger.setLevel(getattr(logging, "DEBUG"))
    logging.basicConfig(format=fmt)

    run_tandem(args.fastq, args.output, args.match, args.mismatch, args.delta,
            args.pm, args.pi, args.minscore, args.maxperiod,
            args.no_generate_fastq, args.mask, args.html, args.trf_path)
    return


if __name__ == "__main__":
    main()
