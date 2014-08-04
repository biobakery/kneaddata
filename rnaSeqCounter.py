import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bowtie2", nargs="+", default=[],
            help="Bowtie2 output files")
    parser.add_argument("--bwa", nargs="+", default=[],
            help="BWA output files")
    parser.add_argument("--list", nargs="+", default=[],
            help="List of contaminant reads")
    parser.add_argument("--orig", nargs="+", default=[],
            help="Original input files")
    parser.add_argument("-t", "--trim", nargs="+", default=[], 
            help="Trimmed files")
    parser.add_argument("-p", "--parent-dir", 
            help="Parent directory of the files")
    parser.add_argument("-o", "--output", help="Output .csv file")

    args = parser.parse_args()

    regex = r'\d\d\d\d\d\d_([A-Z][a-z]*)_([A-Za-z-]*)_'

if __name__ == '__main__'():
    main()
