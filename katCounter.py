import resultParser
import argparse
import numpy as np
import itertools

def combine(match):
    return (str(match.group(1)))

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

    regex = r'([A-Za-z_]+):'

    for f in args.orig:
        print("Original input file: " + f)
        katCounterFastq = resultParser.FastqCounter(pattern=regex,
                combineName=combine)
        katCounterFastq.count(f, bSingleEnd=False)
        print(katCounterFastq.get())

    for f in args.bowtie2:
        katCounterSams = resultParser.SamCounter(pattern=regex,
                combineName=combine)
        print("Bowtie2 output file: " + f)
        katCounterSams.count(f, bSingleEnd=True)
        print(katCounterSams.get())

    for f in args.bwa:
        katCounterSams = resultParser.SamCounter(pattern=regex,
                combineName=combine)
        print("BWA output file: " + f)
        katCounterSams.count(f, bSingleEnd=True)
        print(katCounterSams.get())

    for f in args.bmtagger:
        katCounterBMT = resultParser.BMTOutCounter(pattern=regex,
                combineName=combine)
        print("BMTagger output file: " + f)
        katCounterBMT.count(f, bSingleEnd=False)
        print(katCounterBMT.get())

    return 0


if __name__ == '__main__':
    main()
