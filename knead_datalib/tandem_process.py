import re
import os
import sys
import time
import logging
import argparse
from itertools import izip_longest
from __init__ import mktempfifo

import threading
import Queue

############## constants ######################
# number of lines per entry in trf output
TANDEM_LENGTH = 2

# number of lines per entry in fastq 
FASTQ_LENGTH = 4
###############################################


# depreciated
def generate_fastq(fastq_file, tandem_out_file, out_file):
    tandem_out = None
    fastq = None
    out = None

    def _write_fastq(fastq_lines):
        try:
            out.write("".join(fastq_lines))
        except TypeError:
            # TypeError can happen when you do 
            # "".join([str, str, str, None]). This means that the
            # input fastq didn't have 4 lines per read
            logging.critical("You probably passed in a bad fastq file, %s, one"
                    " that doesn't have 4 lines" " per read" % fastq_file)
            raise
    try:
        tandem_out = open(tandem_out_file, "r")
        fastq = open(fastq_file, "r")
        out = open(out_file, "w")

        # read FASTQ_LENGTH lines at a time
        fastq_iter = izip_longest(*[fastq]*FASTQ_LENGTH)
        # read TANDEM_LENGTH lines at a time
        tandem_iter = izip_longest(*[tandem_out]*TANDEM_LENGTH)
        for (fastq_lines, tandem_out_lines) in izip_longest(fastq_iter,
                tandem_iter):
            if tandem_out_lines == None:
                assert(fastq_lines != None)
                # exhausted tandem_iter
                _write_fastq(fastq_lines)
            elif fastq_lines[0] != tandem_out_lines[0]:
                # fastq header is NOT a tandem repeat, so keep it
                _write_fastq(fastq_lines)
    finally:
        # close everything, even if there was an exception above
        for outfile in [tandem_out, fastq, out]:
            if outfile != None:
                outfile.close()


# depreciated
def generate_masked(fastq_file, masked_file, out_file):
    fasta_head = re.compile(r'^>')
    masked = None
    fastq = None
    out = None

    def _write_output(out_fp, fastq_header, seq, plus, quals):
        seq = seq + "\n"
        assert(len(seq) == len(quals))
        assert(header[1:] == fastq_header[1:])
        out_fp.write("".join([fastq_header, seq, plus, quals]))

    try:
        masked = open(masked_file, "r")
        fastq = open(fastq_file, "r")
        out = open(out_file, "w")

        mask_line = masked.readline()

        header = None
        fastq_header = None
        plus = None
        quals = None
        seq = ""
        for fastq_lines in izip_longest(*[fastq]*FASTQ_LENGTH):
            seq = ""
            fastq_header = fastq_lines[0]
            plus = fastq_lines[2]
            quals = fastq_lines[3]

            assert(plus.strip() == '+')

            while mask_line:
                if re.search(fasta_head, mask_line):
                    if header != None:
                        _write_output(out, fastq_header, seq, plus, quals)
                        header = mask_line
                        mask_line = masked.readline()
                        break
                    header = mask_line
                else:
                    seq += mask_line.rstrip()
                mask_line = masked.readline()

        # flush the last guy
        _write_output(out, fastq_header, seq, plus, quals)

    finally:
        # close everything if there was an exception above
        for outfile in [masked, fastq, out]:
            if outfile != None:
                outfile.close()


def process_tandem(fastq, fasta, tandem, output, maskfile):
    print("PROCESS_TANDEM")
    mask = (maskfile != None)
    generate_fastq = (tandem != None)
    assert(mask or generate_fastq)

    def _write_masked(out_fp, fastq_header, mask_header, seq, plus, quals):
        seq = seq + "\n"
        assert(len(seq) == len(quals))
        assert(mask_header[1:] == fastq_header[1:])
        out_fp.write("".join([fastq_header, seq, plus, quals]))

    fastq_fp = None
    fasta_fp = None
    tandem_fp = None
    maskfile_fp = None
    fastqout_fp = None
    maskout_fp = None

    fasta_head = re.compile(r'^>')

    try:
        fastq_fp = open(fastq, "r")
        fasta_fp = open(fasta, "w")
        if generate_fastq:
            print("OPENING " + tandem)
            tandem_fp = open(tandem, "r")
            print("DONE OPENING " + tandem)
            fastqout_fp = open(output + ".fastq", "w")
        if mask:
            print("OPENING " + maskfile)
            maskfile_fp = open(maskfile, "r")
            print("DONE OPENING MASK")
            if generate_fastq:
                maskout_fp = open(output + ".fastq.mask", "w")
            else:
                maskout_fp = open(output + ".fastq", "w")

        if generate_fastq:
            tandem_out_iter = izip_longest(*[tandem_fp]*TANDEM_LENGTH)
            tandem_header = None

        mask_header = None
        fastq_iter = izip_longest(*[fastq_fp]*FASTQ_LENGTH)
        fastq_lines = None
        for fastq_lines in fastq_iter:

            fastq_header = fastq_lines[0]
            fastq_seq = fastq_lines[1]
            fastq_plus = fastq_lines[2]
            fastq_quals = fastq_lines[3]
            print(fastq_lines)
            for f in [fastq_header, fastq_seq, fastq_plus, fastq_quals]:
                assert(f != None)
            assert(len(fastq_seq) == len(fastq_quals))
            assert(fastq_plus.rstrip() == "+")

            # write to fasta file
            fasta_header = fastq_header.replace("@", ">", 1).rstrip()
            fasta_fp.write(fasta_header + "\n")
            fasta_fp.write(fastq_seq.rstrip() + "\n")
            print(fasta_header)
            print(fastq_seq.rstrip())
            print("WROTE TO FASTA")

            tandem_header = None

            if generate_fastq:
                if tandem_header == None:
                    try:
                        # HANGING HERE
                        print("TRYING TO READ TANDEM_OUT")
                        tandem_out_lines = tandem_out_iter.next()
                        tandem_header = tandem_out_lines[0]
                    except StopIteration:
                        tandem_header = "StopIteration"
                print(tandem_header)
                if fastq_header != tandem_header:
                    # write when the headers don't match up
                    fastqout_fp.write("".join([fastq_header, fastq_seq,
                        fastq_plus, fastq_quals]))
                else:
                    # reset when they are equal, so we can read another header
                    tandem_header = None
                    
            if mask:
                print("READING MASK FILE")
                mask_line = maskfile_fp.readline()
                seq = ""
                while mask_line:
                    if re.search(fasta_head, mask_line):
                        # a new header is found
                        if mask_header != None:
                            # if not the first header, write
                            _write_masked(maskout_fp, fastq_header, mask_header,
                                    seq, fastq_plus, fastq_quals)
                            mask_header = mask_line
                            break
                        # first iteration, when we have no information
                        # set mask_header to the first line in mask file
                        mask_header = mask_line
                    else:
                        # this is a sequence
                        seq += mask_line.rstrip()
                    mask_line = maskfile_fp.readline()
        # out of the while loop
        if mask:
            # flush the last guy
            _write_masked(maskout_fp, fastq_header, mask_header, seq,
                    fastq_plus, fastq_quals)
    finally:
        for outfile in [fastq_fp, fasta_fp, tandem_fp, maskfile_fp, fastqout_fp,
            maskout_fp]:
            if outfile != None:
                outfile.close()


def handle_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq", help="Input FASTQ file")
    parser.add_argument("fasta", help="Temporary fasta file")
    parser.add_argument("output", help="Output file")
    parser.add_argument("--tandem", help=("Input file from TRF -ngs, "
                                    "containing tandem repeats"))
    parser.add_argument("--maskfile", help="TRF mask file")

    args = parser.parse_args()

    if (args.tandem == None) and (args.maskfile == None):
        parser.error(("\nYou must set at least one of --tandem or --maskfile."
                    " Exiting..."))
    return args


def main():
    args = handle_cli()
    logging.debug("Running tandem_process.py with: " + str(args))
    process_tandem(args.fastq, args.fasta, args.tandem, args.output,
            args.maskfile)

if __name__ == "__main__":
    main()
