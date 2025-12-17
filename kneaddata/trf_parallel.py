#!/usr/bin/env python

"""
KneadData parallel TRF

This script wraps bowtie2 to allow for discordant paired alignments output to fastq.

Dependencies: TRF

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
import argparse
import tempfile
import subprocess

try:
    from kneaddata import utilities
except ImportError:
    sys.exit("Please install kneaddata")


def parse_arguments(args):
    """
    Parse the arguments from the user
    """

    parser = argparse.ArgumentParser(
        description="Kneaddata trf parallel\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "--input",
        help="the fasta file of reads",
        required=True)
    parser.add_argument(
        "--output",
        help="the trf output file (.dat merged from chunks)",
        required=True)
    parser.add_argument(
        "--trf-path",
        help="full path to the trf executable (or directory containing it)",
        required=True)
    parser.add_argument(
        "--trf-options",
        help="space delimited list of trf options (e.g. '2 7 7 80 10 50 500 -h -ngs')",
        required=True)
    parser.add_argument(
        "--nproc",
        default=1,
        type=int,
        help="total number of processes to use")

    return parser.parse_args()


def run_trf(input, trf_path, trf_options, nproc, output, verbose=True):
    """Run TRF with the options provided, optionally in parallel.

    When nproc > 1:
      - The input FASTA is split into nproc chunks (2 lines per read).
      - Each chunk is written to a short-named temp file.
      - TRF is run on each temp file in parallel.
      - The resulting .dat files are merged into the final output.

    To avoid TRF segfaults due to extremely long paths, this version:
      - Uses short prefixes for temporary files.
      - Allows overriding the temp directory via KNEADDATA_TRF_TMP.
    """

    tempfile_list = []
    tempfile_written_list = []
    datfile_list = []
    datfile_to_write_list = []
    commands = []

    # Choose a directory for TRF temp files.
    # Allow override via env var to keep paths short (recommended).
    tmp_dir = os.environ.get("KNEADDATA_TRF_TMP")
    if tmp_dir:
        try:
            os.makedirs(tmp_dir, exist_ok=True)
        except OSError:
            # If mkdir fails, fall back to output dir
            tmp_dir = os.path.dirname(output) or "."
    else:
        tmp_dir = os.path.dirname(output) or "."

    # Helper: numeric options part for TRF .dat naming
    # Original code used:
    #   new_file + ".".join(trf_options.split("-")[0].split(" ")) + "dat"
    # which corresponds to <input>.<numeric_opts>.dat
    numeric_opts_raw = trf_options.split("-")[0].split()
    numeric_opts = ".".join(numeric_opts_raw) if numeric_opts_raw else ""
    if numeric_opts:
        numeric_opts = "." + numeric_opts  # leading dot, e.g. ".2.7.7.80.10.50.500"

    # Single-process mode: just run TRF directly
    if nproc == 1:
        # NOTE: If input path itself is extremely long and TRF still chokes,
        # you may want to symlink input to a short path before invoking this.
        commands.append([
            [trf_path, input] + trf_options.split(" "),
            "trf",
            [input],
            [output],
            output
        ])

        utilities.start_processes(commands, nproc, verbose)
        return

    # Multi-process mode
    # Count total lines in input
    total_lines = 0
    with open(input) as file_handle:
        for line in file_handle:
            total_lines += 1

    # Split the input into multiple files and run in parallel
    for i in range(int(nproc)):
        # Use a short, stable prefix instead of the full output basename
        fd, new_file = tempfile.mkstemp(
            prefix="kd_trf_{:03d}_".format(i),
            suffix=".fa",
            dir=tmp_dir
        )
        os.close(fd)

        tempfile_list.append(new_file)
        # Expected TRF .dat output for this temp input
        datfile_name = new_file + numeric_opts + ".dat"
        datfile_list.append(datfile_name)

    # Write the input file into all temp output files (FASTA: 2 lines per read)
    output_file_number = 0
    lines_per_file = int(total_lines / int(nproc)) if nproc > 0 else total_lines
    lines_written = 0
    file_handle_write = None

    for read_line in utilities.read_file_n_lines(input, 2):
        if not file_handle_write:
            file_handle_write = open(tempfile_list[output_file_number], "wt")
            tempfile_written_list.append(tempfile_list[output_file_number])
            datfile_to_write_list.append(datfile_list[output_file_number])

        file_handle_write.write("".join(read_line))

        lines_written += 2
        # Only roll over to a new file if we haven't reached the last one
        if (lines_written > lines_per_file and
                output_file_number < (len(tempfile_list) - 1)):
            file_handle_write.close()
            lines_written = 0
            output_file_number += 1
            file_handle_write = None

    if file_handle_write is not None:
        file_handle_write.close()

    # Run TRF on each temporary chunk
    for i, temp_in, temp_out in zip(
        range(len(tempfile_written_list)),
        tempfile_written_list,
        datfile_to_write_list,
    ):
        trf_command = [trf_path, temp_in] + trf_options.split(" ")
        commands.append([
            trf_command,
            "trf{}".format(i),
            [temp_in],
            [temp_out],
            temp_out
        ])

    utilities.start_processes(commands, nproc, verbose)

    # Merge all of the outputs to the final output file
    with open(output, "w") as file_write:
        for datfile in datfile_to_write_list:
            with open(datfile) as file_read:
                for line in file_read:
                    file_write.write(line)

    # Remove temp files
    for filename in tempfile_list + datfile_to_write_list:
        try:
            os.remove(filename)
        except EnvironmentError:
            print("Unable to remove temp file: " + filename)


def main():
    # parse the command line arguments
    args = parse_arguments(sys.argv)

    # remove extra quotes if included
    args.trf_options = args.trf_options.replace("'", "")

    # add exe name if not included
    if not args.trf_path.endswith("trf"):
        if not args.trf_path.endswith("/"):
            args.trf_path += "/"
        args.trf_path += "trf"

    run_trf(args.input, args.trf_path, args.trf_options, args.nproc, args.output)


if __name__ == "__main__":
    main()

