import re
import argparse

def fastq_to_fasta(fastq_in, fasta_out=None):
    # idea copied from Humann2

    with open(fastq_in, "r") as fp_fastq_in:
        fp_fasta_out = None

        if fasta_out:
            fp_fasta_out = open(fasta_out, "w")

        line = fp_fastq_in.readline()
        while line:
            if re.search(r'^@', line):
                sequence_id = line.replace("@", ">", 1).rstrip()
                line = fp_fastq_in.readline()
                sequence=""
                while line:
                    if re.search(r'^\+', line):
                        if not fp_fasta_out:
                            print(sequence_id)
                            print(sequence)
                        else:
                            fp_fasta_out.write(sequence_id+"\n")
                            fp_fasta_out.write(sequence+"\n")
                        break
                    else:
                        sequence += line.rstrip()
                    line = fp_fastq_in.readline()
            line = fp_fastq_in.readline()
        
        if fp_fasta_out:
            fp_fasta_out.close()
    return


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq", help="input fastq")
    parser.add_argument("--fasta", default=None, help=("output fasta. If not"
        " specified, the output will be printed to stdout"))
    args = parser.parse_args()

    fastq_to_fasta(args.fastq, args.fasta)


if __name__ == "__main__":
    main()
