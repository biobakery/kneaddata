import argparse
import resultParser
import mergesams


'''
Input file specifications:
Aligner name, original, combined, silva, human [newline]
Aligner name, original, combined, silva, human [newline]
...
'''

# constants
# regex to match the FASTA headers
REGEX = r'\d\d\d\d\d\d_([A-Z][a-z]*)_([A-Za-z-]*)_'

# things we are trying to count
FIELDNAMES = ["OrigTotal",
             "OrigHuman",
             "OrigSilva",
             "MergedTotal",
             "MergedHuman",
             "MergedSilva",
             "HumanTotal",
             "HumanHuman",
             "SilvaTotal",
             "SilvaSilva"]

def combine(match):
    return (str(match.group(1) + " " + match.group(2)))


def read_input(strInfile, iNumAligners, fir_delim=",", sec_delim=None):
    ''' 
    Given an input file, returns a list of lists of string lists. Each entry in
    the big list are all the files for each aligner. Each entry in each sublist
    is a list of the files for original, silva, human, and combined. Some
    entries may be empty in the 3rd level list. 
    '''
    # for sanity checks
    iLineCounter = 0
    strAlignerName = None
    lllInFiles = [None for i in range(iNumAligners)]
    with open(strInfile, "r") as f:
        for strLine in f:
            strLineSplit = strLine.split(fir_delim)
            strAlignerName = strLineSplit[0]
            try:
                lllInFiles[iLineCounter] = [s.split(sec_delim) for s in
                        strLineSplit[1:]]
            except IndexError:
                print("You seem to have entered an invalid input file or an invalid number of aligners. Please double-check your inputs.")
                raise
            iLineCounter += 1

    if iLineCounter != iNumAligners:
        raise IOError("You seem to have entered an invalid input file or an invalid number of aligners. Please double-check your inputs.")

    print("Input files:")
    print(lllInFiles)
    return(lllInFiles)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--aligners", type=int, required=True,
            help="Number of aligners")
    parser.add_argument("-p", "--parent-dir", 
            help="Parent directory of the files")
    parser.add_argument("-o", "--output", default="out.csv",
            help="Output .csv file")
    parser.add_argument("-i", "--input", required=True, 
            help="Input file. See documentation for specifics")
    parser.add_argument("--delimiter", default=",", 
            help="Delimiter in input file")

    args = parser.parse_args()

    if args.aligners <= 0: 
        raise IOError("Number of aligners must be greater than 0")

    if not args.parent_dir:
        args.parent_dir = os.getcwd()

    lllInputs = read_input(strInfile=args.input, iNumAligners=args.aligners,
            fir_delim=args.delimiter)

    iLenFieldNames = len(FIELDNAMES)

    arrResults = np.zeros((args.aligners,iLenFieldNames))
    
    # count the `original' files
    for input_index in range(args.aligners):
        arrTempResult = np.zeros(iLenFieldNames)
        # lInput = [original, combined, silva, human]
        for i in range(len(lInput)):
            if i == 1 # and is sam and does not exist
                # merge files, make a BMTOutCounter
            elif isSam:
                # make a SamCounter for each of the files
            elif isOut:
                # make a BMTOutCounter for each of the files
            else:
                # you done fcked up, A-ARON

            # map and count everything

            # store answer
        arrResults[input_index] = arrTempResult


        
if __name__ == '__main__'():
    main()
