'''
Ad Hoc script to filter Silva database FASTAs. First, it changes U's to T's.
Then, it removes spaces in the sequence entries. Finally, it removes sequences
that are not from bacteria or archaea.
'''

import argparse
import re

def filter_silva(infile, outfile):
    '''
    Converts all the U to T from the input FASTA, infile. Writes the output to
    the output FASTA specified by outfile
    '''
    f_in = open(infile, "r")
    f_out = open(outfile, "w")

    regex = r'\d+ ([A-Za-z]+);'

    fIsBacteria = False
    for line in f_in:
        if line[0] == '>':
            # Only keep those that are bacteria or archaea
            match = re.search(regex, line)
            if match:
                if match.group(1) == 'Bacteria' or match.group(1) == 'Archaea':
                    fIsBacteria = True
                    new_line = '>' + "000000_Silva_rRNA_" + line[1:]
                    f_out.write(new_line)
                    continue
                else:
                    fIsBacteria = False
                    continue
            else:
                print("No match found in the following line: " + line)

        if fIsBacteria:
            # preserve upper and lower case
            new_line = line.replace('U', 'T')
            new_line = new_line.replace('u', 't')
            new_line = new_line.replace(' ', '')
            f_out.write(new_line)

    f_in.close()
    f_out.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input fasta containing RNA data")
    parser.add_argument("outfile", help="output file")

    args = parser.parse_args()

    filter_silva(args.infile, args.outfile)

if __name__ == '__main__':
    main()
