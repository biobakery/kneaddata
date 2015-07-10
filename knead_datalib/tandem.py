import os
import sys
import Queue
import logging
import argparse
import threading
import itertools
import subprocess

import time

# import stuff from the __init__.py in this directory
from __init__ import parse_positive_int, mkfifo_here, try_create_dir, process_return
here = os.path.dirname(os.path.realpath(__file__))

############## constants ######################
# number of lines per entry in fastq 
FASTQ_LENGTH = 4

# sentinel to indicate we've reached the end of the queue
SENTINEL_DONE = "DONE"
SENTINEL_CONTINUE = "CONTINUE"
###############################################


def _trf(fasta, trf_out_queue, match=2, mismatch=7, delta=7, pm=80,
        pi=10, minscore=50, maxperiod=500, html=False, trf_path="trf"):
    '''
    Generate the TRF command and run it using subprocess. TRF is run such that
    the reads identified as containing tandem repeat are written to a fifo.
    Then, read TRF's output and write a tuple (header, TRF output string) to a
    Queue.Queue object.

    fasta: input fasta filename for TRF
    trf_out_queue: A Queue.Queue object where TRF's output will be written
    trf_path: path to trf. If not specified, defaults to "trf"

    match, mismatch, delta, pm, pi, minscore, and maxperiod are all TRF
    parameters. 

    html: Tells TRF whether to generate a HTML output
    '''

    trf_args = map(str, [match, mismatch, delta, pm, pi, minscore, maxperiod])
    trf_cmd = [trf_path, fasta] + trf_args + ["-ngs"]

    if not html:
        trf_cmd += ["-h"]

    logging.debug("Running trf with %s", trf_cmd)
    trf_proc = subprocess.Popen(trf_cmd, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, bufsize=1)

    line = trf_proc.stdout.readline().rstrip()
    header = None
    output_lines = []
    while line:
        if line[0] == "@":
            if header != None:
                assert(len(output_lines) > 0)
                trf_out_queue.put( (SENTINEL_CONTINUE, header, output_lines) )
                output_lines = []
            header = line
        else:
            output_lines.append(line)
        line = trf_proc.stdout.readline().rstrip()
    # flush last one
    trf_out_queue.put( (SENTINEL_CONTINUE, header, output_lines) )

    trf_proc.stdout.close()
    trf_out_queue.put( (SENTINEL_DONE, None, None) )


def _fastq_to_fasta(in_fastq, fasta_fifo, qual_queue):
    fasta_fifo_fp = None
    in_fastq_fp = open(in_fastq, "r")
    try:
        fasta_fifo_fp = open(fasta_fifo, "w")
        fastq_iter = itertools.izip_longest(*[in_fastq_fp]*FASTQ_LENGTH)
        for (fastq_header, fastq_seq, fastq_plus, fastq_quals) in fastq_iter:
            # checking
            for f in [fastq_header, fastq_seq, fastq_plus, fastq_quals]:
                assert(f != None)
            assert(len(fastq_seq) == len(fastq_quals))
            assert(fastq_plus.rstrip() == "+")

            print("THREAD WRITING TO FASTA")
            fasta_header = fastq_header.replace("@", ">", 1).rstrip()
            fasta_fifo_fp.write(fasta_header + "\n")
            fasta_fifo_fp.write(fastq_seq.rstrip() + "\n")

            print("THREAD PUTTING INTO QUEUE")
            qual_queue.put( (fastq_header.rstrip(), fastq_seq.rstrip(),
                fastq_quals.rstrip()) )
    finally:
        for fp in [in_fastq_fp, fasta_fifo_fp]:
            if fp != None:
                fp.close()


def _convert(qual_queue, trf_out_queue, fastq_out, mask_out):
    fastq_out_fp = None
    mask_out_fp = None

    if fastq_out != None:
        fastq_out_fp = open(fastq_out, "w")
    if mask_out != None:
        mask_out_fp = open(mask_out, "w")

    # assumes everything has been stripped of trailing newline
    def _write(fastq_header, trf_header, fastq_seq, trf_output, fastq_qual):
        assert(len(fastq_seq) == len(fastq_qual))
        if (fastq_out != None) and (fastq_header != trf_header):
            _write_fastq(fastq_header, fastq_seq, fastq_qual)
        if mask_out != None:
            _write_masked(fastq_header, fastq_seq, trf_output, fastq_qual)

    def _write_masked(fastq_header, fastq_seq, trf_output, fastq_qual):
        masked_seq = fastq_seq
        for out_line in trf_output:
            # use TRF's first two outputs to tell where the beginning and end of
            # the tandem repeated section is. TRF uses 1-indexing, so we must
            # subtract 1 to be compatible with Python
            spl = out_line.split(" ")
            lower = int(spl[0]) - 1
            upper = int(spl[1]) - 1
            len_rep = upper - lower
            masked_seq = masked_seq[:lower] + 'N' * len_rep + masked_seq[upper:]

        mask_out_fp.write(fastq_header + "\n")
        mask_out_fp.write(masked_seq + "\n")
        mask_out_fp.write("+\n")
        mask_out_fp.write(fastq_qual + "\n")

    def _write_fastq(fastq_header, fastq_seq, fastq_qual):
        fastq_out_fp.write(fastq_header + "\n")
        fastq_out_fp.write(fastq_seq + "\n")
        fastq_out_fp.write("+\n")
        fastq_out_fp.write(fastq_qual + "\n")

    is_more = True
    trf_header = None
    trf_output = []

    while True:
        fastq_header, fastq_seq, fastq_qual = qual_queue.get()
        print(fastq_header, fastq_seq, fastq_qual)
        if (trf_header == None) and (len(trf_output) == 0):
            # try to get more TRF output
            if is_more:
                (sentinel, trf_header, trf_output) = trf_out_queue.get()
                if sentinel == SENTINEL_DONE:
                    is_more = False
                    trf_header = None
                    trf_output = []
            else:
                _write(fastq_header, fastq_seq, trf_output, fastq_qual)
                continue

        _write(fastq_header, trf_header, fastq_seq, trf_output, fastq_qual)
        trf_header = None
        trf_output = []
        qual_queue.task_done()


def run_tandem(fastq, output, match=2, mismatch=7, delta=7, pm=80, pi=10,
        minscore=50, maxperiod=500, generate_fastq=True, mask=False, html=False,
        trf_path="trf"):

    fasta_fname = os.path.basename(fastq) + ".fasta"
    trf_args = map(str, [match, mismatch, delta, pm, pi, minscore, maxperiod])
    #trf_out_fname = fasta_fname + "." + ".".join(trf_args) + ".stream.dat"

    # set output properly
    fastq_out = None
    mask_out = None
    if generate_fastq:
        fastq_out = output + ".fastq"
    if mask:
        if generate_fastq:
            mask_out = output + ".fastq.mask"
        else:
            mask_out = output + ".mask"

    #with mkfifo_here((fasta_fname, trf_out_fname)) as filenames:
    with mkfifo_here((fasta_fname, )) as filenames:
        print(filenames)
        qual_queue = Queue.Queue()
        trf_out_queue = Queue.Queue()

        thread_fastq_fasta = threading.Thread(target=_fastq_to_fasta,
                args=(fastq, filenames[0], qual_queue))
        thread_fastq_fasta.daemon = True
        thread_fastq_fasta.start()

        thread_trf = threading.Thread(target=_trf, args=(filenames[0],
            trf_out_queue, match, mismatch, delta, pm, pi, minscore, maxperiod,
            html, trf_path))
        thread_trf.daemon = True
        thread_trf.start()

        thread_convert = threading.Thread(target=_convert, args=(qual_queue,
            trf_out_queue, fastq_out, mask_out))
        thread_convert.daemon = True
        thread_convert.start()

        qual_queue.join()


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

    #qual_queue = Queue.Queue()
    #qual_queue.put(("@read1","ATCG", "AAAA"))

    #trf_out_queue = Queue.Queue()
    #trf_out_queue.put((SENTINEL_CONTINUE, "@read1", ["1 2 other stuff"]))
    #fastq_to_fasta(args.fastq, "fasta", q)
    #_trf("fasta", q)
    #print(q.get())
    #thread_convert = threading.Thread(target=_convert, args=(qual_queue,
    #    trf_out_queue, "fastq_out", "mask_out"))
    #thread_convert.daemon = True
    #thread_convert.start()
    #_convert(qual_queue, trf_out_queue, "fastq_out", "mask_out")
    #qual_queue.join()
    #return

    run_tandem(args.fastq, args.output, args.match, args.mismatch, args.delta,
            args.pm, args.pi, args.minscore, args.maxperiod,
            args.no_generate_fastq, args.mask, args.html, args.trf_path)


if __name__ == "__main__":
    main()
