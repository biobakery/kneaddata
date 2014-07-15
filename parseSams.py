import getCounts
import argparse
import os
import numpy as np
import constants_knead_data as const
import re

ORIG_ENDINGS = ['.fir.fastq']
TRIM_ENDINGS = ['_pe.sam', '_se_1.sam', '_se_2.sam']
NOTRIM_ENDING = '_pe_notrim.sam'

def getFiles(strDataset, strParentDir, strAligner, fTrimmed):
    ''' Get file names we need to count '''
    orig_fnames = [os.path.join(strParentDir, strDataset, strDataset + ending)
            for ending in ORIG_ENDINGS]    

    trimmed_fnames = []
    sam_fnames = []
    if fTrimmed:
        trimmed_fnames = [os.path.join(strParentDir, strDataset, strDataset +
                ".out" + ending) for ending in const.TRIM_PE_ENDINGS[1:]]
        sam_fnames = [os.path.join(strParentDir, strDataset, strDataset + "-" +
                strAligner + ending) for ending in TRIM_ENDINGS]
    else:
        sam_fnames = [os.path.join(strParentDir, strDataset, strDataset + "-" +
                strAligner + NOTRIM_ENDING)]

    return (orig_fnames, trimmed_fnames, sam_fnames)


def countSam(strSamfile):
    ''' Count the total number of aligned reads and the total number of aligned
    human reads in the sam file
    '''
    iTotalAligned = 0
    iAlignedHuman = 0
    regex = r'\d\d\d\d\d\d_([A-Z][a-z]*)_([a-z]*)_'
    
    # determine single end or pair end based on file name
    single_end = False
    strLastLetters = strSamfile[-9:]
    if strLastLetters == TRIM_ENDINGS[1] or strLastLetters == TRIM_ENDINGS[2]:
        single_end = True

    with open(strSamfile, "r") as f:
        for line in f:
            # comment line, skip
            if line[0] == '@':
                pass
            else:
                split_line = line.split("\t")
                if split_line[2].strip() != "*":
                    iTotalAligned += 1
                    match = re.search(regex, line)
                    if match:
                        if str(match.group(1) + " " + match.group(2)) == 'Homo sapiens': 
                            iAlignedHuman += 1
    # .sam files for paired ends have one entry for every paired end read
    if not single_end:
        iTotalAligned = iTotalAligned/2
        iAlignedHuman = iAlignedHuman/2

    return (iTotalAligned, iAlignedHuman)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("datasets", nargs="+",
        help="A list of the different data sets you want to parse")
    parser.add_argument("--parent-dir", default = None,
            help="parent directory of the datasets")
    parser.add_argument("--savefile", help="Output file for numpy array (as csv)")
    parser.add_argument("--aligner", required = True, help="Which aligner used")
    parser.add_argument("--trimmed", default=False, action="store_true",
            help="Running other aligners against trimmed output")

    args = parser.parse_args()
    if args.parent_dir == None:
        args.parent_dir = os.getcwd()

    if os.path.isfile(os.path.join(args.parent_dir, args.savefile)):
        print("Output file already exists! Code did not execute")
        return 1

    lstrFieldNames = ['OrigTotal', 'OrigHuman', 'TrimTotal', 'TrimHuman',
        'OutTotal', 'OutHuman']
    iLenField = len(lstrFieldNames)
    iLenDatasets = len(args.datasets)

    arrResult = np.array([[-1 for i in range(iLenField)] for i in
        range(iLenDatasets)], dtype = "float_")

    fnames_all = [getFiles(dataset, args.parent_dir, args.aligner, args.trimmed)
        for dataset in args.datasets]

    for i in range(iLenDatasets):
        (lstr_orig_fnames, lstr_trimmed_fnames, lstr_sam_fnames) = fnames_all[i]

        temp_result = np.array([-1 for j in range(iLenField)], dtype="float_")

        # get counts for original .fastq files
        orig_counts = [getCounts.getCount(fname) for fname in lstr_orig_fnames]
        temp_result[0] = orig_counts[0][0]
        temp_result[1] = orig_counts[0][1]

        # get counts for trimmed files
        if lstr_trimmed_fnames == []:
            temp_result[2] = temp_result[0]
            temp_result[3] = temp_result[1]
        else:
            print(lstr_trimmed_fnames)
            trim_counts = [getCounts.getCount(fname) for fname in
                lstr_trimmed_fnames]
            print(trim_counts)
            temp_result[2] = np.sum([f[0] for f in trim_counts])
            temp_result[3] = np.sum([f[1] for f in trim_counts])

        # get counts for output 
        out_counts = [countSam(fname) for fname in lstr_sam_fnames]
        print(out_counts)
        temp_result[4] = np.sum([count[0] for count in out_counts])
        temp_result[5] = np.sum([count[1] for count in out_counts])

        arrResult[i] = temp_result

    print(lstrFieldNames)
    print(arrResult)
    # save the file as a csv
    np.savetxt(args.savefile, arrResult, fmt="%d", delimiter=",", newline="\n",
            header=",".join(lstrFieldNames), comments="")

if __name__ == '__main__':
    main()
