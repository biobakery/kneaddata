'''
Code to downsample a ribosomal RNA (rRNA) database like Silva. It will only
include a small percentage of all the rRNAs, chosen randomly from the database.
'''
import argparse
import random

def downsample(infile, outfile, fltProp):
    f_in = open(infile, "r")
    f_out = open(outfile, "w")

    random.seed(0)
    iOrigReads = 0
    iDownsampledReads = 0
    fUse= False

    for line in f_in:
        if line[0] == '>':
            # only take fltProp of the reads
            iOrigReads += 1
            if random.random() < fltProp:
                fUse = True
                iDownsampledReads += 1
            else:
                fUse = False
        if fUse:
            f_out.write(line)

    f_in.close()
    f_out.close()

    return(iOrigReads, iDownsampledReads)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input fasta containing RNA data")
    parser.add_argument("outfile", help="output file")
    parser.add_argument("probability", type=float, 
            help="probability to keep a specific read")

    args = parser.parse_args()

    if args.probability > 1.0 or args.probability < 0.0:
        print("Invalid input!")
        raise IOError("Probability must be between 0 and 1")


    iOrig, iDown= downsample(args.infile, args.outfile, args.probability)

    print("Number of original entries: " + str(iOrig))
    print("Number of downsampled entries: " + str(iDown))
    print("Proportion of entries kept: " + str(iDown / (iOrig + 0.0)))

if __name__ == '__main__':
    main()
