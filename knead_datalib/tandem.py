import os
import logging
import subprocess

# import stuff from the __init__.py in this directory
from . import divvy_threads, mktempfifo, process_return
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
    logging.debug("Running fastq_to_fasta with:\n" + str(fastq_to_fasta_cmd))
    fastq_to_fasta = subprocess.Popen(fastq_to_fasta_cmd, 
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.PIPE)
    return(fastq_to_fasta)


def _trf(fasta, outfp, match=2, mismatch=7, delta=7, pm=80, pi=10, minscore=50,
        maxperiod=500, dat=False, mask=False, html=False, trf_path="trf"):
    '''
    Generate the TRF command and run it using subprocess. TRF is run such that
    the reads identified as containing tandem repeats are piped from stdout to
    another program

    fasta: input fasta filename for TRF
    outfp: A file pointer that specifies where to write the TRF stdout
    trf_path: path to trf. If not specified, defaults to "trf"

    match, mismatch, delta, pm, pi, minscore, and maxperiod are all TRF
    parameters. 

    dat: Tells TRF whether to generate a .dat file
    mask: Tells TRF whether to generate a mask file.
    html: Tells TRF whether to generate a HTML output
    '''

    trf_args = map(str, [match, mismatch, delta, pm, pi, minscore, maxperiod])
    trf_cmd = [trf_path, fasta] + trf_args + ["-ngs"]

    if mask:
        trf_cmd += ["-m"]
    if not html:
        trf_cmd += ["-h"]

    logging.debug("Running trf with:\n" + str(trf_cmd))
    trf = subprocess.Popen(trf_cmd, stdout=outfp, stderr=subprocess.PIPE)
    return(trf)


def _generate_fastq(fastq, trf_output, out):
    cmd = ["python", os.path.join(here, "tandem_process.py"), 
            fastq, trf_output, out]
    logging.debug("Running tandem_process.py generate_fastq with:\n" + str(cmd))
    proc = subprocess.Popen(cmd)
    return(proc)


def _generate_mask(fastq, trf_output, out):
    cmd = ["python", os.path.join(here, "tandem_process.py"), 
            fastq, trf_output, out, "--mask"]
    logging.debug("Running tandem_process.py generate_mask with:\n" + str(cmd))
    proc = subprocess.Popen(cmd)
    return(proc)
    

def _convert(fastq, trf_output, out, mask_fname, generate_fastq, mask):
    converter_procs = []
    converter_names = []

    if generate_fastq: 
        converter_procs.append(_generate_fastq(fastq, trf_output, out))
        converter_names.append("tandem_process.py generate_fastq")

    if mask:
        mask_out = out
        if generate_fastq:
            mask_out = out + ".mask"
        assert(mask_fname != None)
        converter_procs.append(_generate_mask(fastq, mask_fname,  mask_out))
        converter_names.append("tandem_process.py generate_mask")

    return(converter_procs, converter_names)


def trf(fastq, out, match=2, mismatch=7, delta=7, pm=80, pi=10, minscore=50,
        maxperiod=500, generate_fastq=True, dat=False, mask=False, html=False,
        trf_path="trf"):

    fasta_fname = os.path.basename(fastq) + ".fasta"
    # file name for trf mask (trf default provided, can't change it)
    trf_args = map(str, [match, mismatch, delta, pm, pi, minscore, maxperiod])
    trf_out_fname = fasta_fname + "." + ".".join(trf_args) + "stream.dat"
    mask_fname = None
    if mask:
        mask_fname = fasta_fname + "." + ".".join(trf_args) + ".mask"

    with mktempfifo((fasta_fname, trf_out_fname)) as filenames:

        fastq_to_fasta = _fastq_to_fasta(fastq, filenames[0])
        
        trf_out = filenames[1]
        if dat:
            trf_out = fasta_fname + "." + ".".join(trf_args) + ".dat"

        converter_procs, converter_names = _convert(fastq, trf_out, out,
                                                    mask_fname, generate_fastq, 
                                                    mask)
        trf_out_fp = open(trf_out, "w")
        trf = _trf(filenames[0], trf_out_fp, match, mismatch, delta, pm, pi,
                minscore, maxperiod, dat, mask, html,trf_path)
        try:
            procs = [fastq_to_fasta, trf]
            proc_names = ["fastq_to_fasta " + fastq, "trf " + filenames[0]]
            for (proc, name) in zip(procs, proc_names):
                stdout, stderr = proc.communicate()
                retcode = proc.returncode
                process_return(name, retcode, stdout, stderr)
        finally:
            trf_out_fp.close()

        for (proc, name) in zip(converter_procs, converter_names):
            stdout, stderr = proc.communicate()
            retcode = proc.returncode
            process_return(name, retcode, stdout, stderr)

        if not logging.getLogger().isEnabledFor(logging.DEBUG) and mask:
            os.remove(mask_fname)


