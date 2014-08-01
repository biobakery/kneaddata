import resultParser
import argparse
import numpy as np
import itertools

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bowtie2", nargs="+", default=[],
            help="Bowtie2 output files")
    parser.add_argument("-w", "--bwa", nargs="+", default=[],
            help="BWA output files")
    parser.add_argument("-t", "--bmtagger", nargs="+", default=[],
            help="BMTagger output files")
    parser.add_argument("-o", "--orig", nargs="+", default=[],
            help="Original input files")

    args = parser.parse_args()

    regex = r'^([A-Za-z_]+):'

    katCounterSams = resultParser.SamCounter(pattern=regex)
    katCounterFastq = resultParser.FastqCounter(pattern=regex)
    katCounterBMT = resultParser.BMTOutCounter(pattern=regex)
    
    for f in args.orig:
        print("Original input file: " + f)
        print(katCounterFastq.count(f, bSingleEnd=False))

    for f in args.bowtie2:
        print("Bowtie2 output file: " + f)
        print(katCounterSams.count(f, bSingleEnd=False))

    for f in args.bwa:
        print("BWA output file: " + f)
        print(katCounterSams.count(f, bSingleEnd=False))

    for f in args.bmtagger:
        print("BMTagger output file: " + f)
        print(katCounterBMT.count(f, bSingleEnd=False))

    return 0


if __name__ == '__main__':
    main()
