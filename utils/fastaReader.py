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
