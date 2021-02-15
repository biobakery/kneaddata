#!/usr/bin/env python

"""
KneadData bowtie2 discordant pairs

This script wraps bowtie2 to allow for discordant paired alignments output to fastq.

Dependencies: Bowtie2

Copyright (c) 2017 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import sys
import os
import logging
logger=logging.getLogger(__name__)

try:
    import argparse
except ImportError:
    sys.exit("Please upgrade to python 2.7+")
    
import string
import tempfile
import subprocess

try:
    from kneaddata import utilities
except ImportError:
    sys.exit("Please install kneaddata")
    
def reverse_complement(sequence):
    try:
        table=string.maketrans("ATCG","TAGC")
    except AttributeError:
        # allow for python3 in which maketrans is from str class
        table=str.maketrans("ATCG","TAGC")
    return sequence.translate(table)[::-1]

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
    description= "Kneaddata bowtie2 discordant pairs\n",
    formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument(
        "-1",
        dest="pair1",
        help="the fastq file of pair1 reads",
        required=True)
    parser.add_argument(
        "-2",
        dest="pair2",
        help="the fastq file of pair2 reads",
        required=True) 
    parser.add_argument(
        "-x",
        dest="index",
        help="the database index file",
        required=True)
    parser.add_argument(
        "--un-pair",
        help="the name of the output files for the paired reads without any alignments",
        required=True)
    parser.add_argument(
        "--al-pair",
        help="the name of the output files for the paired reads with concordant alignments",
        required=True)
    parser.add_argument(
        "--un-single",
        help="the name of the output files for the orphan reads without alignments",
        required=True)
    parser.add_argument(
        "--al-single",
        help="the name of the output files for the orphan reads with alignments",
        required=True)
    parser.add_argument(
        "-U",
        dest="orphan",
        help="the fastq files of orphan reads in comma-delimited list")
    parser.add_argument(
        "-S",
        dest="sam",
        help="the file to write the sam output")
    parser.add_argument(
        "--bowtie2",
        help="the path to the bowtie2 executable",
        default="bowtie2")
    parser.add_argument(
        "--threads",
        help="the number of threads to use",
        default=1)
    parser.add_argument(
        "--bowtie2-options",
        help="the bowtie2 options to apply")
    parser.add_argument(
        "--mode",
        choices=["strict","unpaired"],
        help="the run mode")
    parser.add_argument(
        "--cat-pairs",
        action="store_true",
        help="concatenate pair files before aligning so reads are aligned as single end")
    parser.add_argument(
        "--reorder",
        action="store_true",
        help="print the sequences in the same order as the input files")
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="print diagnostics")
    
    return parser.parse_args()

def bowtie2_unpaired_command(args, input_fastq_path, output_sam_path, maintain_input_ordering):

    # Bowtie2 from input as unpaired reads to output, using the provided index
    command=[args.bowtie2, "-U", input_fastq_path, "-S", output_sam_path, "-x", args.index]
    

    # Runtime options
    command+=["--threads",str(args.threads)]
    if args.bowtie2_options:
        command+=utilities.format_options_to_list([args.bowtie2_options])

    if not args.verbose:
        command+=["--quiet"]
   
    # Maintain ordering of the fastq, if required
    if maintain_input_ordering:
        command+=["--reorder"]

    return command

def run_command(command, **kwargs):
    try:
        return_code=subprocess.check_call(command, **kwargs)
    except (EnvironmentError, subprocess.CalledProcessError) as e:
        message="Unable to run: " +" ".join(command)
        if hasattr(e, 'output') and e.output:
            message+="\nError message returned:\n" + e.output
        sys.exit(message)
    

def read_sam_line(line):
    data=line.rstrip().split("\t")

    query_id = data[0][:-1]
    mate = data[0][-1]
    is_aligned = not (int(data[1]) & 4)
    read="\n".join(["@"+data[0], data[9], "+", data[10], ""])
    return (query_id, mate, is_aligned, read)

def process_alignments(pair1_sam, pair2_sam, orphan_sam, aligned_pair, unaligned_pair, aligned_orphan, unaligned_orphan, treat_pair_as_aligned_if_either_read_aligned):
    """ Read through the paired sam alignments and organize into the output files """

    open_files={
       1: {
         'both_aligned': open(aligned_pair.replace("%","1"),"wt"),
         'both_unaligned': open(unaligned_pair.replace("%","1"),"wt"),
         'only_this_aligned': open(aligned_orphan.replace("%","1"),"wt"),
         'only_this_unaligned': open(unaligned_orphan.replace("%","1"),"wt")
       },
       2: {
         'both_aligned': open(aligned_pair.replace("%","2"),"wt"),
         'both_unaligned': open(unaligned_pair.replace("%","2"),"wt"),
         'only_this_aligned': open(aligned_orphan.replace("%","2"),"wt"),
         'only_this_unaligned': open(unaligned_orphan.replace("%","2"),"wt")
       }
    }
    counts={
         'both_aligned': 0,
         'both_unaligned': 0,
         'only_this_aligned': 0,
         'only_this_unaligned': 0
    }


    with open(pair1_sam) as fh_1, open(pair2_sam) as fh_2:
        line_1=fh_1.readline()
        while line_1.startswith("@"):
            line_1=fh_1.readline()
        line_2=fh_2.readline()
        while line_2.startswith("@"):
            line_2=fh_2.readline()
        line_count = 1
        while line_1 and line_2:
            query_id_1, mate_1, is_aligned_1, read_1 = read_sam_line(line_1)
            query_id_2, mate_2, is_aligned_2, read_2 = read_sam_line(line_2)
            if not (query_id_1 == query_id_2 and mate_1 != mate_2):
                raise ValueError(
                    "sam files do not match on line {0}: found IDs {1}{2} and {3}{4}"
                        .format(line_count, query_id_1, mate_1, query_id_2, mate_2)
                )
            aa = is_aligned_1 and is_aligned_2 or (treat_pair_as_aligned_if_either_read_aligned and (is_aligned_1 or is_aligned_2))

            x1 = 'both_aligned' if aa else 'only_this_aligned' if is_aligned_1 else 'only_this_unaligned' if is_aligned_2 else 'both_unaligned'
            x2 = 'both_aligned' if aa else 'only_this_aligned' if is_aligned_2 else 'only_this_unaligned' if is_aligned_1 else 'both_unaligned'
                
            counts[x1]+=1
            open_files[1][x1].write(read_1)
            open_files[2][x2].write(read_2)
            line_1=fh_1.readline()
            line_2=fh_2.readline()
            line_count+=1

        while line_1:
            query_id_1, mate_1, is_aligned_1, read_1 = read_sam_line(line_1)
            x1 = 'only_this_aligned' if is_aligned_1 else 'only_this_unaligned'
            counts[x1]+=1
            open_files[1][x1].write(read_1)
            line_1=fh_1.readline()

        while line_2:
            query_id_2, mate_2, is_aligned_2, read_2 = read_sam_line(line_2)
            x2 = 'only_this_aligned' if is_aligned_2 else 'only_this_unaligned'
            counts[x2]+=2
            open_files[2][x2].write(read_2)
            line_2=fh_2.readline()

    logger.info("Paired reads: {0} both aligned, {1} both unaligned, {2} only pair1 aligned, {3} only pair2 aligned"
     .format(counts["both_aligned"], counts["both_unaligned"], counts["only_this_aligned"], counts["only_this_unaligned"])
    )
    if orphan_sam:
        orphan_counts = {'only_this_aligned' : 0, 'only_this_unaligned': 0}
        with open(orphan_sam) as fh:
            line = fh.readline()
            while line.startswith("@"):
                line=fh.readline()
            while line:
                query_id, mate, is_aligned, read = read_sam_line(line)
                m = 2 if mate == '2' else 1
                x = 'only_this_aligned' if is_aligned else 'only_this_unaligned'
                orphan_counts[x]+=1
                open_files[m][x].write(read)
                line = fh.readline()
        logger.info("Orphan reads: {0} aligned, {1} unaligned"
          .format(orphan_counts["only_this_aligned"], orphan_counts["only_this_unaligned"])
        )

    for fh in open_files[1].values():
       fh.close()
    for fh in open_files[2].values():
       fh.close()

def temp_sam_path(args, name):
    output_dir=os.path.dirname(args.un_pair)
    file_out, sam = tempfile.mkstemp("kneaddata_"+name+".sam", dir=output_dir)
    os.close(file_out)
    return sam

def main():
    # parse the command line arguments
    args = parse_arguments(sys.argv)
    temp_files = []

    ch = logging.StreamHandler()
    ch.setFormatter(logging.Formatter('%(asctime)s kneaddata_bowtie2_discordant_pairs - %(message)s'))
    logger.addHandler(ch)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    logger.debug("info")
    pair1_sam=temp_sam_path(args, "pair1")
    temp_files.append(pair1_sam)
    command=bowtie2_unpaired_command(args, input_fastq_path = args.pair1, output_sam_path = pair1_sam, maintain_input_ordering = True)
    logger.debug("Aligning pair1 reads: "+ " ".join(command))
    run_command(command)

    pair2_sam=temp_sam_path(args, "pair2")
    temp_files.append(pair2_sam)
    command=bowtie2_unpaired_command(args, input_fastq_path = args.pair2, output_sam_path = pair2_sam, maintain_input_ordering = True)
    logger.debug("Aligning pair2 reads: "+ " ".join(command))
    run_command(command)

    orphan_sam = None
    if args.orphan:
        orphan_sam=temp_sam_path(args, "orphan")
        temp_files.append(orphan_sam)
        command=bowtie2_unpaired_command(args, input_fastq_path = args.orphan, output_sam_path = orphan_sam, maintain_input_ordering = args.reorder)
        logger.debug("Aligning orphan reads: "+ " ".join(command))
        run_command(command)

    logger.debug("Processing the alignments")
    process_alignments(pair1_sam, pair2_sam, orphan_sam,
      aligned_pair = args.al_pair,
      unaligned_pair = args.un_pair,
      aligned_orphan = args.al_single,
      unaligned_orphan = args.un_single,
      treat_pair_as_aligned_if_either_read_aligned = (args.mode == "strict")
    )

    if args.sam and args.sam != os.devnull:
        logger.debug("Aggregating .sam output to " + args.sam)
        command = ["cat"]
        command.extend(temp_files)
        with open(args.sam, 'w') as fh: 
            run_command(command, stdout = fh)

    logger.debug("Removing temporary files")
    for file in temp_files:
        utilities.remove_file(file)

if __name__ == "__main__":
    main()
    
