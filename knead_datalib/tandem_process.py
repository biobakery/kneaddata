import re
import sys
import logging
import argparse
from itertools import izip_longest


############## constants ######################
# number of lines per entry in trf output
TANDEM_LENGTH = 2

# number of lines per entry in fastq 
FASTQ_LENGTH = 4
###############################################


def generate_fastq(fastq_file, tandem_out_file, out_file):
    tandem_out = open(tandem_out_file, "r")
    fastq = open(fastq_file, "r")
    out = open(out_file, "w")

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
        # if the tandem_out file is empty, we still want to write the fastq file
        wrote = False
        # read FASTQ_LENGTH lines at a time
        for fastq_lines in izip_longest(*[fastq]*FASTQ_LENGTH):
            # read TANDEM_LENGTH lines at a time
            for tandem_out_lines in izip_longest(*[tandem_out]*TANDEM_LENGTH):
                # fastq header is NOT a tandem repeat, so keep it
                if fastq_lines[0] != tandem_out_lines[0]:
                    _write_fastq(fastq_lines)
                    wrote = True
            if wrote == False:
                _write_fastq(fastq_lines)
    finally:
        # close everything, even if there was an exception above
        tandem_out.close()
        fastq.close()
        out.close()


def generate_masked(fastq_file, masked_file, out_file):
    masked = open(masked_file, "r")
    fastq = open(fastq_file, "r")
    out = open(out_file, "w")
    fasta_head = re.compile(r'^>')

    def _write_output(out_fp, fastq_header, seq, plus, quals):
        seq = seq + "\n"
        assert(len(seq) == len(quals))
        assert(header[1:] == fastq_header[1:])
        out_fp.write("".join([fastq_header, seq, plus, quals]))

    try:
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
        masked.close()
        fastq.close()
        out.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq", help="Input FASTQ file")
    parser.add_argument("tandem", help="Input file from TRF")
    parser.add_argument("output", help="Output file")
    parser.add_argument("--mask", default=False, action="store_true", 
            help=("If set, specifies that the input is a masked FASTA. A "
                 "masked FASTQ will be generated"))

    args = parser.parse_args()

    if args.mask:
        generate_masked(args.fastq, args.tandem, args.output)
    else:
        generate_fastq(args.fastq, args.tandem, args.output)

    return 0


if __name__ == "__main__":
    main()
