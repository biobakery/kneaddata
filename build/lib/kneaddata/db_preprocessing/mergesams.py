import argparse

def merge(infiles, outfile):

    setReads = set()
    for infile in infiles:
        with open(infile, "r") as fileIn:
            for strLine in fileIn:
                if strLine.startswith('@'):
                    continue
                strSplit = strLine.split("\t")
                if strSplit[2] != '*':
                    setReads.add(strSplit[0])

    with open(outfile, "w") as fileOut:
        for strRead in setReads:
            fileOut.write(strRead + "\n")
    return(len(setReads))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infiles", nargs="+", 
            help="sam files you wish to merge")
    parser.add_argument("outfile", help="output file")

    args = parser.parse_args()

    iNumUniqReads = merge(args.infiles, args.outfile)
    print("Number of merged reads: " + str(iNumUniqReads))
    return 0

if __name__ == '__main__':
    main()
