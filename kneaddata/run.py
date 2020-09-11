"""
KneadData: run module

Copyright (c) 2015 Harvard School of Public Health

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

import os
import re
import sys
import time
import shutil
import logging
import itertools
import subprocess
import gzip
import tempfile
from kneaddata import utilities
from kneaddata import config

# name global logging instance
logger=logging.getLogger(__name__)

def fastqc(fastqc_path, output_dir, input_files, threads, verbose):
    """ Run fastq on the input files, placing output in directory provided """
    # write output to a subfolder
    fastqc_output_dir=os.path.join(output_dir, "fastqc")
    # create the directory if it does not already exist
    utilities.create_directory(fastqc_output_dir)
    
    command=[fastqc_path]+input_files+["--threads",str(threads)]+["--outdir",fastqc_output_dir]+["--extract"]
    # run fastqc command
    utilities.run_command(command,"fastqc",input_files,[],None,verbose,exit_on_error=True)


def align(infile_list, db_prefix_list, output_prefix, remove_temp_output,
          bowtie2_path, threads, processors, bowtie2_opts, verbose, 
          discordant=None, reorder=None, serial=None, decontaminate_pairs=None):
    """ Runs bowtie2 on a single-end sequence file or a paired-end set of files. 
    For each input file set and database provided, a bowtie2 command is generated and run."""

    # determine if the input are paired reads
    is_paired = (len(infile_list) == 2)

    # create the bowtie2 commands
    commands = []
    all_outputs_to_combine = []
    database_names = []
    if discordant:
        all_outputs_to_combine = [[],[],[]]
        database_names = [[],[],[]]
    all_contaminated_outputs = []
    bowtie2_command = [bowtie2_path, "--threads", str(threads)] + bowtie2_opts
    
    for basename, fullpath in _prefix_bases(db_prefix_list):
        output_str = output_prefix + "_" + basename + "_bowtie2"
        cmd = bowtie2_command + ["-x", fullpath]
        if discordant:
            # if running in serial mode, use the last set of outputs as input
            if serial and all_outputs_to_combine[0]:
                current_infile_list=[all_outputs_to_combine[0][-1][-2],all_outputs_to_combine[0][-1][-1]]
                # check for orphans from the prior filter run
                if all_outputs_to_combine[1]:
                    current_infile_list.append(all_outputs_to_combine[1][-1][0])
                    
                if all_outputs_to_combine[2]:
                    current_infile_list.append(all_outputs_to_combine[2][-1][0])
            else:
                current_infile_list=infile_list
            
            # run the pairs allowing for all alignments (including those generating orphans)
            cmd=["kneaddata_bowtie2_discordant_pairs","--bowtie2",bowtie2_path,"--threads", str(threads),"-x",fullpath,"--mode",decontaminate_pairs]
            
            if bowtie2_opts:
                cmd+=["--bowtie2-options","\""+" ".join(bowtie2_opts)+"\""]
            
            # add the input and output files for the pairs
            pair_output_str = output_str + "_paired"
            cmd += ["-1", current_infile_list[0], "-2", current_infile_list[1],
                    "--un-pair", pair_output_str + "_clean_%" + config.fastq_file_extension]
            cmd+=["--al-pair", pair_output_str + "_contam_%" + config.fastq_file_extension]
            
            all_contaminated_outputs.append(pair_output_str + "_contam_1" + config.fastq_file_extension)
            all_contaminated_outputs.append(pair_output_str + "_contam_2" + config.fastq_file_extension)
                
            outputs_to_combine= [pair_output_str + "_clean_1" + config.fastq_file_extension, 
                                  pair_output_str + "_clean_2" + config.fastq_file_extension]
            all_outputs_to_combine[0].append(outputs_to_combine)
            
            # add the orphan input and output files
            single_output_str = output_str + "_unmatched_%"
            if len(current_infile_list) > 2:
                cmd+=["-U", ",".join(current_infile_list[2:])]
            cmd+=["--un-single", single_output_str + "_clean" + config.fastq_file_extension]
            cmd+=["--al-single", single_output_str + "_contam" + config.fastq_file_extension]
            
            if reorder:
                cmd+=["--reorder"]
            
            all_contaminated_outputs.append(output_str + "_unmatched_1_contam" + config.fastq_file_extension)
            all_contaminated_outputs.append(output_str + "_unmatched_2_contam" + config.fastq_file_extension)
            all_outputs_to_combine[1].append([output_str + "_unmatched_1_clean" + config.fastq_file_extension])
            all_outputs_to_combine[2].append([output_str + "_unmatched_2_clean" + config.fastq_file_extension])     

            database_names[0]+=[basename,basename]
            database_names[1]+=[basename]
            database_names[2]+=[basename]
                            
        elif is_paired:
            if serial and all_outputs_to_combine:
            # if running in serial mode, take the last set of output files as inputs
                cmd += ["-1", all_outputs_to_combine[-1][-2], "-2", all_outputs_to_combine[-1][-1]]
            else:
                cmd += ["-1", infile_list[0], "-2", infile_list[1]]
            cmd+= ["--un-conc", output_str + "_clean_%" + config.fastq_file_extension]
            cmd+=["--al-conc", output_str + "_contam_%" + config.fastq_file_extension]
            all_contaminated_outputs.append(output_str + "_contam_1" + config.fastq_file_extension)
            all_contaminated_outputs.append(output_str + "_contam_2" + config.fastq_file_extension)
            outputs_to_combine = [output_str + "_clean_1" + config.fastq_file_extension, 
                                  output_str + "_clean_2" + config.fastq_file_extension]
            all_outputs_to_combine.append(outputs_to_combine)
            database_names+=[basename,basename]

        else:
            # if running in serial mode, take the last output file as input
            if serial and all_outputs_to_combine:
                cmd += ["-U", all_outputs_to_combine[-1][0]]
            else:
                cmd += ["-U", infile_list[0]]
            cmd += ["--un", output_str + "_clean" + config.fastq_file_extension]
            cmd+=["--al", output_str + "_contam" + config.fastq_file_extension]
            all_contaminated_outputs.append(output_str + "_contam" + config.fastq_file_extension)
            outputs_to_combine = [output_str + "_clean" + config.fastq_file_extension]
            all_outputs_to_combine.append(outputs_to_combine)
            database_names+=[basename]

        if remove_temp_output:
            # if we are removing the temp output, then write the sam output to dev null to save space
            sam_out = os.devnull
        else:
            sam_out = output_str + ".sam"
        cmd += [ "-S", sam_out ]
        
        commands.append([cmd,"bowtie2",infile_list,outputs_to_combine,None])

    # run the bowtie2 commands with the number of processes specified
    utilities.start_processes(commands,processors,verbose)

    # write out total number of contaminated reads found
    for file in all_contaminated_outputs:
        total_contaminates=utilities.count_reads_in_fastq_file(file, verbose)
        message="Total contaminate sequences in file ( " + file + " ) : " + str(total_contaminates)
        logger.info(message)
        if verbose:
            print(message)   

    # if bowtie2 produced output, merge the files from multiple databases
    combined_outs = []
    if all_outputs_to_combine:
        if discordant:
            combined_outs1 = combine_fastq_output_files(all_outputs_to_combine[0], output_prefix + "_paired", remove_temp_output, database_names[0])
            combined_outs2 = combine_fastq_output_files(all_outputs_to_combine[1], output_prefix + "_unmatched_1", remove_temp_output, database_names[1])
            combined_outs3 = combine_fastq_output_files(all_outputs_to_combine[2], output_prefix + "_unmatched_2", remove_temp_output, database_names[2])
            combined_outs = [combined_outs1,combined_outs2,combined_outs3]
        else:
            combined_outs = combine_fastq_output_files(all_outputs_to_combine, output_prefix, remove_temp_output, database_names)

    return combined_outs

def write_tagged_sequences_from_fastq(input_fastq, bmtagger_output, output_fastq, verbose):
    """ Find the sequences bmtagger has tagged as contaminates from the extract output file """
    
    # store all of the sequences bmtagger has not tagged as contaminates
    untagged_sequences=set()
    for lines in utilities.read_file_n_lines(bmtagger_output,4):
        untagged_sequences.add(lines[0])
                        
    try:
        file_handle_write=open(output_fastq,"w")
    except EnvironmentError:
        sys.exit("ERROR: Unable to open file: " + output_fastq)
        
    tagged_sequences=0
    for lines in utilities.read_file_n_lines(input_fastq,4):
        # check if the sequence was identified by bmtagger
        if not lines[0] in untagged_sequences:
            tagged_sequences+=1
            file_handle_write.write("".join(lines))
        
    # log the number of sequences
    message="Total contaminate sequences in file ( " + output_fastq + " ): " + str(tagged_sequences)
    logger.info(message)
    if verbose:
        print(message)

def tag(infile_list, db_prefix_list, remove_temp_output, output_prefix,
        bmtagger_path, processes, verbose):
    """ Runs BMTagger on a single-end sequence file or a paired-end set of files. 
    For each input file set and database provided, a bmtagger command is generated and run."""

    # determine if the input are paired reads
    is_paired = (len(infile_list) == 2)

    # create a temp directory for bmtagger
    tempdir=tempfile.mkdtemp(prefix=os.path.basename(output_prefix)+'_temp_',dir=os.path.dirname(output_prefix))

    # create the bmtagger commands
    commands = []
    all_outputs_to_combine = []
    contaminated_outputs = []
    bmtagger_command = [bmtagger_path, "-q", "1", "-1", infile_list[0], 
                        "-T", tempdir,"--extract"]

    # build arguments
    database_names=[]
    for (basename, fullpath) in _prefix_bases(db_prefix_list):
        prefix = output_prefix + "_" + basename + "_bmtagger"
        cmd = bmtagger_command + ["-b", str(fullpath + ".bitmask"),
                                  "-x", str(fullpath + ".srprism"),
                                  "-o", prefix]
        if is_paired:
            cmd += ["-2", infile_list[1]]
            outputs_to_combine = [prefix + "_1" + config.fastq_file_extension,
                                  prefix + "_2" + config.fastq_file_extension]
            # name the corresponding contaminated output files
            contaminated_outputs.append([prefix + "_contam_1" + config.fastq_file_extension,
                                         prefix + "_contam_2" + config.fastq_file_extension])
            database_names+=[basename,basename]
        else:
            outputs_to_combine = [prefix + config.fastq_file_extension]
            # name the corresponding contaminated output files
            contaminated_outputs.append([prefix + "_contam" + config.fastq_file_extension])
            database_names+=[basename]

        commands.append([cmd,"bmtagger",infile_list,outputs_to_combine,None])
        all_outputs_to_combine.append(outputs_to_combine)
        
    # run the bmtagger commands with the number of processes specified
    utilities.start_processes(commands,processes,verbose)
    
    # write the files of contaminate sequences
    for index, outputs in enumerate(all_outputs_to_combine):
        for input_fastq, bmtagger_output, contam_output_fastq in zip(infile_list, outputs, contaminated_outputs[index]):
            write_tagged_sequences_from_fastq(input_fastq, bmtagger_output, contam_output_fastq, verbose)

    # remove the temp directory
    try:
        shutil.rmtree(tempdir)
    except EnvironmentError:
        logger.debug("Unable to remove temp directory: " +tempdir)

    # merge the output files from multiple databases
    combined_outs = []
    if all_outputs_to_combine:
        combined_outs = combine_fastq_output_files(all_outputs_to_combine, output_prefix, remove_temp_output,database_names)
                    
    return combined_outs

def intersect_fastq(fastq_files, out_file, remove_temp_output=None):
    """ Intersects multiple fastq files with one another. Includes only the reads (4
    lines long each) that are common to all the files. Writes these reads to the
    output file specified in out_file. 
    """
    
    # optimize for the common case, where we are intersecting 1 file
    if len(fastq_files) == 1:
        if remove_temp_output:
            shutil.move(fastq_files[0], out_file)
        else:
            shutil.copyfile(fastq_files[0], out_file)
    else:
        # store the number of files that contain each sequence
        sequence_count={}
        for fname in fastq_files:
            for lines in utilities.read_file_n_lines(fname, 4):
                sequence_count[lines[0]]=sequence_count.get(lines[0],0)+1
    
        num_files = len(fastq_files)
        with open(out_file, "w") as file_handle:
            # read through one of the files, writing out each sequence that 
            # is found in all of the files
            for lines in utilities.read_file_n_lines(fastq_files[0],4):
                if sequence_count.get(lines[0],0) >= num_files:
                    file_handle.write("".join(lines))

def combine_fastq_output_files(files_to_combine, out_prefix, remove_temp_output, database_names):
    """ Combines fastq output created by BMTagger/bowtie2 on multiple databases and 
    returns a list of output files. Also updates the log file with read counts for the 
    input and output files.
    """
    
    # print out the reads for all files
    utilities.log_read_count_for_files(files_to_combine,"decontaminated",
        "Total reads after removing those found in reference database", database_names)

    # create lists of all of the output files for pair 1 and for pair 2
    files_for_pair1 = [f[0] for f in files_to_combine]
    try:
        files_for_pair2 = [f[1] for f in files_to_combine]
    except IndexError:
        files_for_pair2 = []

    # select an output prefix based on if the outputs are paired or not
    output_file = out_prefix + "_1" + config.fastq_file_extension
    if not files_for_pair2:
        output_file = out_prefix + config.fastq_file_extension

    # create intersect file from all output files for pair 1
    intersect_fastq(files_for_pair1, output_file, remove_temp_output)
    output_files=[output_file]
    
    # create an intersect file from all output files for pair 2
    if files_for_pair2:
        output_file = out_prefix + "_2" + config.fastq_file_extension
        intersect_fastq(files_for_pair2, output_file, remove_temp_output)
        output_files.append(output_file)

    # Get the read counts for the newly merged files
    utilities.log_read_count_for_files(output_files,"final","Total reads after merging results from multiple databases")

    # remove temp files if set
    if remove_temp_output:
        for group in [files_for_pair1, files_for_pair2]:
            for filename in group:
                utilities.remove_file(filename)
                
    return output_files

def _prefix_bases(db_prefix_list):
    """From a list of absolute or relative paths, returns an iterator of the
    following tuple: (basename of all the files, full path of all the files)
    
    If more than one file as the same basename (but have different paths), the
    basenames get post-fixed with indices (starting from 0), based on how their
    basenames sort.
    """
    # sort by basename, but keep full path
    bases = sorted([ (os.path.basename(p), p) for p in db_prefix_list ], 
        key = lambda x: x[0])
    
    processed=set()
    # check for databases with the same basename
    for name, group in itertools.groupby(bases, key=lambda x: x[0]):
        group = list(group)
        if len(group) > 1:
            for i, item in enumerate(group):
                yield ("%s_%i"%(item[0], i), item[1])
                processed.add(item[1])
    
    # return all databases with unique basenames in the original order
    for file in db_prefix_list:
        if not file in processed:
            yield os.path.basename(file), file

def trim(infiles, outfiles_prefix, trimmomatic_path, quality_scores, 
         java_memory, additional_options, threads, verbose):
    """ Creates and runs trimmomatic commands based on input files and options. 
    Returns a list of the output files.
    """

    command = ["java", "-Xmx" + java_memory, "-jar", trimmomatic_path]

    # determine if paired end input files
    paired_end=False
    if len(infiles) == 2:
        paired_end = True
        
    if paired_end:
        # set options for paired end input files
        mode = "PE"
        outfiles = [outfiles_prefix + config.trimomatic_pe_endings[0], 
                    outfiles_prefix + config.trimomatic_pe_endings[2],
                    outfiles_prefix + config.trimomatic_pe_endings[1],
                    outfiles_prefix + config.trimomatic_pe_endings[3]]
    else:
        # set options for single input file
        mode = "SE"
        outfiles = [outfiles_prefix + config.trimomatic_se_ending]
    
    # add positional arguments to command
    command += [mode, "-threads", str(threads), quality_scores] + infiles + outfiles
    # add optional arguments to command
    command += additional_options

    # run trimmomatic command
    utilities.run_command(command,"Trimmomatic",infiles,outfiles,None,verbose,exit_on_error=True)
    
    # now check all of the output files to find which are non-empty and return as 
    # sets for running the alignment steps
    nonempty_outfiles=[]
    outfile_size = [utilities.file_size(file) for file in outfiles]
    if paired_end:
        # if paired fastq files remain after trimming, preserve pairing
        if outfile_size[0] > 0 and outfile_size[2] > 0:
            nonempty_outfiles.append([outfiles[0],outfiles[2]])
        elif outfile_size[0] > 0:
            nonempty_outfiles.append([outfiles[0]])
            # remove the second paired file if empty
            utilities.remove_file(outfiles[2])
        elif outfile_size[2] > 0:
            nonempty_outfiles.append([outfiles[2]])
            # remove the second paired file if empty
            utilities.remove_file(outfiles[0])
        
        # add sequences without pairs, if present
        if outfile_size[1] > 0:
            nonempty_outfiles.append([outfiles[1]])
        else:
            # remove the file if empty
            utilities.remove_file(outfiles[1])
            
        if outfile_size[3] > 0:
            nonempty_outfiles.append([outfiles[3]])
        else:
            # remove the file if empty
            utilities.remove_file(outfiles[3])
        
    else:
        if outfile_size[0] > 0:
            nonempty_outfiles=[[outfiles[0]]]
        else:
            # remove the file if empty
            utilities.remove_file(outfiles[0])
        
    if not nonempty_outfiles:
        sys.exit("ERROR: Trimmomatic created empty output files.")
        
    return nonempty_outfiles
        
def remove_repeats_from_fastq(input_fastq, trf_output, output_fastq):
    """ Remove the sequences from TRF that contain repeats from the output files """
    
    sequences_with_repeats=set()
    try:
        with open(trf_output) as file_handle:
            for line in file_handle:
                # sequences start with "@"
                if line[0] == "@":
                    sequences_with_repeats.add(line)
    except EnvironmentError:
        pass
                
    try:
        file_handle_write=open(output_fastq,"w")
    except EnvironmentError:
        sys.exit("ERROR: Unable to open file: " + output_fastq)
        
    removed_sequences=0
    for lines in utilities.read_file_n_lines(input_fastq,4):
        # check if the sequence was identified by TRF
        if lines[0] in sequences_with_repeats:
            removed_sequences+=1
        else:
            file_handle_write.write("".join(lines))
        
    # log the number of sequences removed for repeats
    logger.info("Total number of sequences with repeats removed from file ( " + 
                input_fastq + " ): " + str(removed_sequences))
        
def tandem(input_files, output_prefix, match, mismatch, delta, pm, pi, minscore,
               maxperiod, trf_path, processors, verbose, remove_temp_output, threads):
    """ Run TRF on all input files """

    # Convert all arguments to strings    
    trf_args = list(map(str, [match, mismatch, delta, pm, pi, minscore, maxperiod]))

    output_files=[]
    pairs=False
    unmatched=1
    output_prefix+=".repeats.removed"
    for input_fastq_files in input_files:
        # Get the names for the output files
        if len(input_fastq_files) > 1:
            pairs=True
            output_fastq_files = [output_prefix + "." + str(i) + config.fastq_file_extension for i in range(1,len(input_fastq_files)+1)]
        elif pairs:
            output_fastq_files = [output_prefix + ".unmatched." + str(unmatched) + config.fastq_file_extension]
            unmatched+=1
        else:
            output_fastq_files = [output_prefix + config.fastq_file_extension]
        
        commands=[]
        temp_fasta_files=[]
        trf_output_files=[]
        for input_fastq in input_fastq_files:
            # create a temp fasta file from the fastq file
            input_fasta = input_fastq.replace(os.path.splitext(input_fastq)[-1],config.fasta_file_extension)
            utilities.fastq_to_fasta(input_fastq, input_fasta)
            temp_fasta_files.append(input_fasta)
            
            trf_output_file = input_fasta+".trf.parameters."+".".join(trf_args)+".dat"
            trf_output_files.append(trf_output_file)
            
            # suppress html output and write reduced data file to standard output
            trf_command=["kneaddata_trf_parallel","--input",input_fasta,"--output",trf_output_file,"--trf-path",
                    trf_path,"--trf-options","'"+" ".join(trf_args+["-h","-ngs"])+"'","--nproc",str(threads)]

            # only run trf if the fasta file is not empty
            if os.path.getsize(input_fasta) > 0:
                commands.append([trf_command,"trf",[input_fasta],[trf_output_file],trf_output_file])
            
        # run the trf commands with the number of processes specified
        utilities.start_processes(commands,processors,verbose)
        
        # remove all fasta files when complete
        for file in temp_fasta_files:
            utilities.remove_file(file)
    
        # use the trf output to print the final fastq output files
        for i in range(len(input_fastq_files)):
            remove_repeats_from_fastq(input_fastq_files[i], trf_output_files[i], output_fastq_files[i])
            
        # remove trf output if remove temp output is set
        if remove_temp_output:
            for file in trf_output_files:
                utilities.remove_file(file)
        
        # sets for running the alignment steps
        if pairs:
            if (len(output_fastq_files)==2):
                output_files.append([output_fastq_files[0],output_fastq_files[1]])
            else:
                output_files.append([output_fastq_files[0]])
        else: 
            output_files.append([output_fastq_files[0]])
            
    return output_files
        
def decontaminate(args, output_prefix, files_to_align):
    """
    Run bowtie2 or bmtagger then trf if set
    """

    # Start aligning
    message="Decontaminating ..."
    print(message)
    logger.info(message)
    possible_orphan = (len(files_to_align) > 1)
    orphan_count = 1
    output_files=[]
    
    # if running bowtie2 with discordant and pairs, run all reads at once
    if not args.bmtagger and args.discordant and isinstance(files_to_align[0], list) and len(files_to_align[0]) == 2:
        alignment_output_files = align([files_to_align[0][0],files_to_align[0][1]]+utilities.resolve_sublists(files_to_align[1:]), 
            args.reference_db, output_prefix, args.remove_temp_output, args.bowtie2_path, args.threads,
            args.processes, args.bowtie2_options, args.verbose, discordant=args.discordant, 
            reorder=args.reorder, serial=args.serial, decontaminate_pairs=args.decontaminate_pairs)
        output_files=alignment_output_files
    else:
        for files_list in files_to_align:
            prefix = output_prefix
            if possible_orphan and (len(files_list) == 1):
                prefix = output_prefix + "_unmatched_" + str(orphan_count)
                orphan_count += 1
            elif len(files_list) == 2:
                prefix = output_prefix + "_paired"
        
            if args.bmtagger:
                alignment_output_files = tag(files_list, args.reference_db,
                             args.remove_temp_output, prefix, args.bmtagger_path,
                             args.processes, args.verbose)
            else:
                alignment_output_files = align(files_list, args.reference_db, prefix, 
                               args.remove_temp_output, args.bowtie2_path, args.threads,
                               args.processes, args.bowtie2_options, args.verbose, serial=args.serial)
             
            output_files.append(alignment_output_files)   
            
    return output_files
        
