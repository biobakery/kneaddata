import random
import argparse

def sample(lElts, iKeep):
    ''' 
    Implements reservoir sampling. Keeps iKeep number of things from an lElts
    iterable.
    '''
    if iKeep == 0:
        return []
    if iKeep < 0:
        raise IOError("Must sample a nonnegative number of elements")

    lRes = []
    iCount = 0
    for elt in lElts:
        if iCount < iKeep:
            lRes.append(elt)
        else:
            iRandInd = random.randint(0, iCount-1)
            if iRandInd < iKeep:
                lRes[iRandInd] = elt
        iCount += 1
    return lRes

def fastaReader(filename):
    '''
    Reads a FASTA file. Returns an iterable, whose elements are each entry
    (header + sequence) in the FASTA file.
    '''
    currHeader = None
    currRead = None
    with open(filename, "r") as f:
        currHeader = f.next()
        currRead = []
        for line in f:
            if line[0] == '>':
                yield (currHeader, "".join(currRead))
                currHeader = line
                currRead = []
            else:
                currRead.append(line)
        # flush the last read
        yield (currHeader, "".join(currRead))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fastas", nargs="+", required=True, 
            help="Input FASTA files to sample from")
    parser.add_argument("-o", "--outfiles", nargs="+", required=True,
            help="Locations to write output")
    parser.add_argument("-n", "--num-samples", type=int, required=True,
            help="Number of samples to take")

    args = parser.parse_args()

    if args.num_samples <= 0:
        raise IOError("Number of samples must be >= 1")

    if len(args.fastas) != len(args.outfiles):
        raise IOError("Number of input FASTAs must be the same as number of output files")

    for fasta, outfile in zip(args.fastas, args.outfiles):
        lResult = sample(lElts = fastaReader(fasta), iKeep = args.num_samples)
        with open(outfile, "w") as f:
            for (header, seq) in lResult:
                f.write(header)
                f.write(seq)

if __name__ == '__main__':
    main()
