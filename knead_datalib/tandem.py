import os
import sys
import logging
import argparse
import itertools
import subprocess

import Queue
#import multiprocessing
import threading

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

    return(trf_proc.returncode)


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

            #logging.debug("THREAD WRITING TO FASTA\n")
            fasta_header = fastq_header.replace("@", ">", 1).rstrip()
            fasta_fifo_fp.write(fasta_header + "\n")
            fasta_fifo_fp.write(fastq_seq.rstrip() + "\n")

            #logging.debug("THREAD PUTTING INTO QUEUE\n")
            qual_queue.put( (SENTINEL_CONTINUE, fastq_header.rstrip(),
                fastq_seq.rstrip(), fastq_quals.rstrip()) )
    finally:
        qual_queue.put( (SENTINEL_DONE, None, None, None) )
        for fp in [in_fastq_fp, fasta_fifo_fp]:
            if fp != None:
                fp.close()


def _convert(qual_queue, trf_out_queue, fastq_out, mask_out):
    # read fastq quality scores from qual_queue, trf output from trf_out_queue,
    # and write a fastq without tandem repeats to fastq_out, and a mask file
    # with tandem repeats masked by N's to mask_out
    fastq_out_fp = None
    mask_out_fp = None

    if fastq_out != None:
        fastq_out_fp = open(fastq_out, "w")
    if mask_out != None:
        mask_out_fp = open(mask_out, "w")

    # assumes everything has been stripped of trailing newline
    def _write_and_set_next(fastq_header, trf_header, fastq_seq, trf_output,
            fastq_qual):
        assert(len(fastq_seq) == len(fastq_qual))
        new_trf_header = None
        new_trf_output = []
        if (fastq_header != trf_header):
            new_trf_header = trf_header
            new_trf_output = trf_output
            if (fastq_out_fp != None):
                _write_fastq(fastq_out_fp, fastq_header, fastq_seq, fastq_qual)
            if (mask_out_fp != None):
                _write_fastq(mask_out_fp, fastq_header, fastq_seq, fastq_qual)
        elif mask_out_fp != None:
            _write_masked(fastq_header, fastq_seq, trf_output, fastq_qual)
        else:
            print("SEQUENCE NOT MATCH")

        return((new_trf_header, new_trf_output))

    def _write_masked(fastq_header, fastq_seq, trf_output, fastq_qual):
        masked_seq = fastq_seq
        print(trf_output)
        for out_line in trf_output:
            # use TRF's first two outputs to tell where the beginning and end of
            # the tandem repeated section is. TRF uses 1-indexing, so we must
            # subtract 1 to be compatible with Python
            spl = out_line.split(" ")
            lower = int(spl[0]) - 1
            upper = int(spl[1]) - 1
            # indices are INCLUSIVE
            len_rep = upper - lower + 1
            masked_seq = (masked_seq[:lower] + ('N' * len_rep) +
                    masked_seq[(upper+1):])

        logging.debug(masked_seq)
        mask_out_fp.write(fastq_header + "\n")
        mask_out_fp.write(masked_seq + "\n")
        mask_out_fp.write("+\n")
        mask_out_fp.write(fastq_qual + "\n")

    def _write_fastq(fastq_fp, fastq_header, fastq_seq, fastq_qual):
        fastq_fp.write(fastq_header + "\n")
        fastq_fp.write(fastq_seq + "\n")
        fastq_fp.write("+\n")
        fastq_fp.write(fastq_qual + "\n")

    is_more = True
    trf_header = None
    trf_output = []

    while True:
        fastq_sentinel, fastq_header, fastq_seq, fastq_qual = qual_queue.get()
        logging.debug(str((fastq_header, fastq_seq, fastq_qual)))
        if fastq_sentinel == SENTINEL_DONE:
            logging.debug("BREAKING")
            break
        if (trf_header == None) and (len(trf_output) == 0):
            # try to get more TRF output
            if is_more:
                (sentinel, trf_header, trf_output) = trf_out_queue.get()
                if sentinel == SENTINEL_DONE:
                    is_more = False
                    trf_header = None
                    trf_output = []
            else:
                # didn't get any TRF output
                (trf_header, trf_output) = _write_and_set_next(fastq_header,
                        trf_header, fastq_seq, trf_output, fastq_qual)
                assert(trf_header == None)
                assert(trf_output == [])
                continue
        # figure out what to write and update trf_header and trf_output
        (trf_header, trf_output) = _write_and_set_next(fastq_header, trf_header,
                fastq_seq, trf_output, fastq_qual)

    logging.debug("TRYING TO CLOSE FDS")
    for fp in [fastq_out_fp, mask_out_fp]:
        logging.debug(fp)
        if fp != None:
            logging.debug("CLOSING " + str(fp))
            fp.close()
            logging.debug("CLOSED " + str(fp))
    logging.debug("TRYING TO RETURN")
    return


def run_tandem(fastq, output, match=2, mismatch=7, delta=7, pm=80, pi=10,
        minscore=50, maxperiod=500, generate_fastq=True, mask=False, html=False,
        trf_path="trf"):

    fasta_fname = os.path.basename(fastq) + ".fasta"
    trf_args = map(str, [match, mismatch, delta, pm, pi, minscore, maxperiod])

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

    with mkfifo_here((fasta_fname, )) as filenames:
        #logging.debug(filenames)
        #qual_queue = multiprocessing.JoinableQueue()
        #trf_out_queue = multiprocessing.JoinableQueue()
        qual_queue = Queue.Queue()
        trf_out_queue = Queue.Queue()

        #thread_fastq_fasta = multiprocessing.Process(target=_fastq_to_fasta,
        #        args=(fastq, filenames[0], qual_queue))
        thread_fastq_fasta = threading.Thread(target=_fastq_to_fasta,
                args=(fastq, filenames[0], qual_queue))
        thread_fastq_fasta.daemon = True
        thread_fastq_fasta.start()

        #thread_trf = multiprocessing.Process(target=_trf, args=(filenames[0],
        #    trf_out_queue, match, mismatch, delta, pm, pi, minscore, maxperiod,
        #    html, trf_path))
        thread_trf = threading.Thread(target=_trf, args=(filenames[0],
            trf_out_queue, match, mismatch, delta, pm, pi, minscore, maxperiod,
            html, trf_path))
        thread_trf.start()

        #thread_convert = multiprocessing.Process(target=_convert,
        #        args=(qual_queue, trf_out_queue, fastq_out, mask_out))
        #thread_convert = threading.Thread(target=_convert, args=(qual_queue,
        #    trf_out_queue, fastq_out, mask_out))
        #thread_convert.daemon = True
        #thread_convert.start()

        _convert(qual_queue, trf_out_queue, fastq_out, mask_out)


def handle_cli():
    parser = argparse.ArgumentParser(description="Remove tandem repeats.")
    parser.add_argument("fastq", help="input fastq")
    parser.add_argument("output", help="output prefix")

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

    return args


def main():
    args = handle_cli()

    # TEMP LOGGING
    #fmt = "%(asctime)s %(levelname)s: %(message)s"
    #logger = logging.getLogger()
    #logger.setLevel(getattr(logging, "DEBUG"))
    #logging.basicConfig(format=fmt)

    run_tandem(args.fastq, args.output, args.match, args.mismatch, args.delta,
            args.pm, args.pi, args.minscore, args.maxperiod,
            args.no_generate_fastq, args.mask, args.html, args.trf_path)


if __name__ == "__main__":
    main()
