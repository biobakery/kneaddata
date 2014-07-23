'''
Quick Python script to parse the output sequences from synmetap before and after
they are run through the knead_data.py pipeline. Accepts .fastq and .out (output
from knead_data.py) files, and reports total number of reads, number of human
reads, and read counts for each species present.
'''

import re
import argparse
from collections import Counter
import sys

def getCount(strFname):
    regex = r'\d\d\d\d\d\d_([A-Z][a-z]*)_([A-Za-z-]*)_'
    counterCount = Counter()

    iTotalReads = 0
    iLineCounter = 0

    iLineMod = 1
    if strFname[-6:] == '.fastq':
        # Only consider every 4th line for .fastq
        iLineMod = 4

    try:
        f = open(strFname, "rU")
        for line in f:
            if (iLineCounter % iLineMod) == 0:
                # keep this number between 1-4 so taking mods is easy
                iLineCounter = 0    

                # use regex to get the genus, species name
                match = re.search(regex, line)
                if match:
                    counterCount[str(match.group(1) + " " + match.group(2))] += 1
                    iTotalReads += 1
                else:
                    # If there is no regex match, print the line so we can see
                    # what is going wrong 
                    pass
                    #print("DID NOT MATCH REGEX")
                    #print(line)
                    # sys.exit(1)
            iLineCounter += 1
    except IOError:
        print("Could not find file " + strFname)
        return None

    return (iTotalReads, counterCount['Homo sapiens'], counterCount)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", 
        help="Input file name") 

    args = parser.parse_args()

    iTotalReads, iHumanReads, counter = getCount(args.filename)
    print("Total reads: " + str(iTotalReads))
    print("Human reads: " + str(iHumanReads))
    print(counter)

if __name__ == '__main__':
    main()
