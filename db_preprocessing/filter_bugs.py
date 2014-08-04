import argparse
import re
import sets

def filter_bugs(infile, outfile):
    '''
    Picks out specific bugs from the silva database.
    '''
    f_in = open(infile, "r")
    f_out = open(outfile, "w")

    regex = r';([A-Za-z]+) ([A-Za-z]+).*\n'

    bugs_list = ["Fusobacterium nucleatum",
                 "Yersinia pestis",
                 "Trichodesmium erythraeum",
                 "Treponema pallidum",
                 "Streptococcus sanguinis",
                 "Prevotella melaninogenica"]
    bugs = sets.ImmutableSet(bugs_list)
    print(bugs)

    fKeep = False
    for line in f_in:
        if line[0] == '>':
            # Only keep those that are in our set
            match = re.search(regex, line)
            if match:
                strBugName = match.group(1) + " " + match.group(2)
                if strBugName in bugs:
                    fKeep = True
                    print(line)
                else:
                    fKeep = False
            else:
                fKeep = False

        if fKeep:
            f_out.write(line)
            # write to own separate file
            with open("_".join(strBugName.split()) + ".fna", "a") as f:
                f.write(line)

    f_in.close()
    f_out.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="input fasta containing RNA data")
    parser.add_argument("outfile", help="output file")

    args = parser.parse_args()

    filter_bugs(args.infile, args.outfile)

if __name__ == '__main__':
    main()
