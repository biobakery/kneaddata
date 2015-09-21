'''
Converts an RNA FASTA file to the corresponding cDNA file. Basically replaces
all the Uracils in the sequence with Thymines. Removes extraneous spaces in the
sequences, too
'''
import argparse

def convert(infile, outfile):
    '''
    Converts all the U to T from the input FASTA, infile. Writes the output to
    the output FASTA specified by outfile
    '''
    f_in = open(infile, "r")
    f_out = open(outfile, "w")

    for line in f_in:
        # don't convert FASTA headers
        if line[0] == '>':
            f_out.write(line)
            continue
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

    convert(args.infile, args.outfile)

if __name__ == '__main__':
    main()
