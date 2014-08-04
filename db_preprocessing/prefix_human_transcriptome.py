'''
Ad Hoc script to filter Human Transcriptome database FASTAs. Just prepends
something to each FASTA header so we can distinguish it more easily when
post-processing data
'''

import argparse
import re

def filter_file(infile, outfile):
    f_in = open(infile, "r")
    f_out = open(outfile, "w")

    for line in f_in:
        new_line = None
        if line[0] == '>':
            new_line = '>' + "000001_Human_mRNA_" + line[1:]
        else:
            new_line = line
        f_out.write(new_line)

    f_in.close()
    f_out.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input fasta containing RNA data")
    parser.add_argument("outfile", help="output file")

    args = parser.parse_args()

    filter_file(args.infile, args.outfile)

if __name__ == '__main__':
    main()
