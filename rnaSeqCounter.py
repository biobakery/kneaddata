import argparse
import resultParser
import mergesams

def combine(match):
    return (str(match.group(1) + " " + match.group(2)))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bowtie2", nargs="+", default=[],
            help="Bowtie2 output files")
    parser.add_argument("--bwa", nargs="+", default=[],
            help="BWA output files")
    parser.add_argument("--list", nargs="+", default=[],
            help="List of contaminant reads")
    parser.add_argument("--orig", default=[],
            help="Original input files")
    parser.add_argument("-t", "--trim", nargs="+", default=[], 
            help="Trimmed files")
    parser.add_argument("-p", "--parent-dir", 
            help="Parent directory of the files")
    parser.add_argument("-o", "--output", help="Output .csv file")

    args = parser.parse_args()

    # regex to match the FASTA headers
    regex = r'\d\d\d\d\d\d_([A-Z][a-z]*)_([A-Za-z-]*)_'

    lstrAligners = ["bowtie2-trim", 
                    "bowtie2-notrim", 
                    "bwa-trim", 
                    "bwa-notrim",
                    "bmtagger"]
    iLenAligners = len(lstrAligners)

    lstrFieldNames = ["OrigTotal",
                      "OrigHuman",
                      "OrigSilva",
                      "MergedTotal",
                      "MergedHuman",
                      "MergedSilva",
                      "HumanTotal",
                      "HumanHuman",
                      "SilvaTotal",
                      "SilvaSilva"]
    iLenFieldNames = len(lstrFieldNames)

    arrResults = np.zeros((iLenAligners,iLenFieldNames))
    
    # count the `original' files
    for i in range(iLenAligners):
        arrTemp = np.zeros(iLenFieldNames)
        lstrFnames = None
        if (i == 1) or (i == 3):
            # count the original input files
            lstrFnames = [os.path.join(args.parent_dir, args.orig)]
            #counterResult = counterFastq.count(strFname)
            #arrTemp[0] = np.sum(counterResult.values())
            #arrTemp[1] = counter["Silva rRNA"]
            #arrTemp[2] = counter["Human mRNA"]
        else:
            # count the trimmed files
            lstrFnames = [os.path.join(args.parent_dir, f) for f in args.trim]

        lstCounters = [resultParser.FastqCounter(pattern=regex,
            combineName=combine) for i in len(lstrFnames)]
        lstResult = [c.count(f) for f in lstrFnames]
        arrTemp[0] = np.sum([np.sum(c.values()) for c in lstResult])
        arrTemp[1] = np.sum(c["Silva rRNA"] for c in lstResult])
        arrTemp[1] = np.sum(c["Human mRNA"] for c in lstResult])

        # merge the SAM files
        lstrNewBowtieFnames = [f + ".uniq" for f in args.bowtie2] 
        lstrNewBwaFnames = [f + ".uniq" for f in args.bwa] 
        for strOldFname, strNewFname in zip(args.bowtie2, lstrNewBowtieFnames):
            merge(strOldFname, strNewFname)

        for strOldFname, strNewFname in zip(args.bowtie2, lstrNewBwaFnames):
            merge(strOldFname, strNewFname)

        # get merged counts



if __name__ == '__main__'():
    main()
