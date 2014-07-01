import getCounts
import constants_knead_data as const
import argparse
import os
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("datasets", nargs="+", 
        help="A list of the different data sets you want to parse")
    parser.add_argument("parent_dir", 
            help="parent directory of the datasets")
    parser.add_argument("savefile", help="Output file for numpy array (as csv)")

    args = parser.parse_args()
    
    # assuming all of our data is paired end, and we are using BMTagger to
    # output the putative human reads in a .out file, instead of producing a
    # .fastq file with the human reads removed.

    if os.path.isfile(os.path.join(args.parent_dir, args.savefile)):
        print("Output file already exists! Code did not execute")
        return 0

    ORIG_ENDINGS = ['.fir.fastq', '.sec.fastq']
    iLenDatasets = len(args.datasets)

    lFileNames = ['' for i in range(iLenDatasets)]
    iFileNameCounter = 0
        
    lstrFieldNames = ['NumReadsOrig', 'NumReadsTrimmedPE',
            'NumReadsTrimmedSE1', 'NumReadsTrimmedSE2', 'HumanOrig',
            'HumanTrimmedPE', 'HumanTrimmedSE1', 'HumanTrimmedSE2',
            'HumanRemovedPE', 'HumanRemovedSE1', 'HumanRemovedSE2',
            'NonHumanRemovedPE', 'NonHumanRemovedSE1', 'NonHumanRemovedSE2']

    iLenField = len(lstrFieldNames)

    lRes = np.array([[-1 for i in range(iLenField)] for i in
        range(iLenDatasets)], dtype="float_")

    for iDataInd in range(iLenDatasets):
        dataset = args.datasets[iDataInd]
        strPathToDataSet = os.path.abspath(
                os.path.join(args.parent_dir, dataset))
        
        lFileNames[iFileNameCounter] = dataset
        iFileNameCounter += 1

        res = np.array([0 for i in range(iLenField)], dtype="float_")
        
        # count the input files
        strFile = os.path.join(strPathToDataSet, dataset + ORIG_ENDINGS[0])
        print(strFile)
        try:
            iTotalReads, iHumanReads, counter = getCounts.getCount(strFile)
            res[0] = iTotalReads
            res[4] = iHumanReads
        except TypeError:
            print("Could not get counts for " + strFile)

        # count the trimmed files
        endings = const.TRIM_PE_ENDINGS[1:] # only include one of the PE files
        for i in range(len(endings)):
            strFile = os.path.join(strPathToDataSet,
                    dataset + ".out" + endings[i])
            print(strFile)
            try:
                iTotalReads, iHumanReads, counter = getCounts.getCount(strFile)
                res[i+1] = iTotalReads
                res[i+5] = iHumanReads
            except TypeError:
                print("Could not get counts for " + strFile)

        # count the BMTagger output
        for i in range(len(const.BMTAGGER_OUTPUT_ENDINGS)):
            strFile = os.path.join(strPathToDataSet,
                    dataset + const.BMTAGGER_OUTPUT_ENDINGS[i])
            print(strFile)
            try:
                iTotalReads, iHumanReads, counter = getCounts.getCount(strFile)
                res[i+8] = iHumanReads
                res[i+11] = iTotalReads - iHumanReads
            except TypeError:
                print("Could not get counts for " + strFile)
        lRes[iDataInd] = res
    print(lFileNames)
    print(lstrFieldNames)
    print(lRes)

    np.savetxt(args.savefile, lRes, fmt="%d", delimiter=",", newline="\n",
            header=",".join(lstrFieldNames), comments="")


if __name__ == '__main__':
    main()
