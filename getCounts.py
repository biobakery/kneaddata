import re
import argparse
from collections import Counter
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", 
        help="Input file name") 

    args = parser.parse_args()

    regex = r'\d\d\d\d\d\d_([A-Z][a-z]*)_([a-z]*)_'
    counterCount = Counter()
    iTotalReads = 0
    iLineCounter = 0

    iLineMod = 1
    if args.filename[-6:] == '.fastq':
        # Only consider every 4th line for .fastq
        iLineMod = 4

    with open(args.filename, "rU") as f:
        for line in f:
            if (iLineCounter % iLineMod) == 0:
                match = re.search(regex, line)
                if match:
                    counterCount[str(match.group(1) + " " + match.group(2))] += 1
                    iTotalReads += 1
                else:
                    print("DID NOT MATCH REGEX")
                    print(line)
                    sys.exit(1)
            iLineCounter += 1

    print("Total reads: " + str(iTotalReads))
    print(counterCount)

if __name__ == '__main__':
    main()
