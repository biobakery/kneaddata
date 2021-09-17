"""
KneadData: config module

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

# Default settings for command line arguments
threads=1
processes=1

log_level_choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"]
log_level=log_level_choices[0]

bmtagger_exe="bmtagger.sh"

trimmomatic_provided_sequencer_source=["NexteraPE","TruSeq2","TruSeq3","none"]
trimmomatic_provided_sequencer_default="NexteraPE"
trimmomatic_adapter_folder=os.path.join(os.path.dirname(os.path.abspath(__file__)),"adapters")
trimmomatic_jar="trimmomatic*"
trimmomatic_memory="500m"
trimmomatic_flag_start="-"
trimmomatic_trim_adapters_option_se="ILLUMINACLIP:$PATH/"+trimmomatic_provided_sequencer_default+"-PE.fa:2:30:10"
trimmomatic_trim_adapters_option_pe="ILLUMINACLIP:$PATH/"+trimmomatic_provided_sequencer_default+"-PE.fa:2:30:10:8:TRUE"
trimmomatic_slidingwindow_option="SLIDINGWINDOW:4:20"
trimmomatic_minlen_option_tag="MINLEN"
trimmomatic_option_delimiter=":"
trimmomatic_min_len_percent=50
default_read_length=100
min_init_read_length=60

# quality score flags for trimmomatic and bowtie2
quality_scores_options=["phred33","phred64"]
quality_scores=quality_scores_options[0]

bowtie2_exe="bowtie2"
bowtie2_flag_start="--"
bowtie2_options=["--very-sensitive-local"]

trf_exe="trf"
trf_match=2
trf_mismatch=7
trf_delta=7
trf_match_probability=80
trf_pi=10
trf_minscore=50
trf_maxperiod=500

samtools_exe="samtools"

#set the locations of data in the sam file
sam_read_name_index=0
sam_flag_index=1
sam_reference_index=2
sam_mapq_index=4
sam_cigar_index=5
sam_read_index=9
sam_start_optional_index=11

sam_read_quality=10
sam_unmapped_flag=0x4
sam_delimiter="\t"

# File endings for BMTagger's required database files
bowtie2_db_endings = [
    ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
bowtie2_large_index_ext = ".rev.1.bt2l"
bmtagger_db_endings = [
    ".bitmask", ".srprism.amp", ".srprism.idx", ".srprism.imp",
    ".srprism.map", ".srprism.pmp", ".srprism.rmp", ".srprism.ss",
    ".srprism.ssa", ".srprism.ssd", ".nhr", ".nin", ".nsq"]
    
    
# File extensions
fastq_file_extension=".fastq"
fasta_file_extension=".fasta"

# Trimmomatic file endings for single end and paired end, respectively
trimomatic_se_ending = ".trimmed.fastq"


# Trimmomatic in paired end mode writes four output files
# The first two files are the trimmed fastq files for input file 1 and 2
# The second two files are the trimmed sequences that lost their pair
# For example the third file is the trimmed sequences from input file 1
# that do not have a pair from input file 2 after trimming
trimomatic_pe_endings = [
    ".trimmed.1.fastq", ".trimmed.2.fastq", 
    ".trimmed.single.1.fastq", ".trimmed.single.2.fastq"]

trimmomatic_pe_names = {
    trimomatic_pe_endings[0]:"pair1", trimomatic_pe_endings[1]:"pair2",
    trimomatic_pe_endings[2]:"orphan1", trimomatic_pe_endings[3]:"orphan2"}

# Set the names of the file types for alignment files
final_file_types = { "paired_1.fastq": "pair1", "paired_2.fastq": "pair2",
    "unmatched_1.fastq":"orphan1", "unmatched_2.fastq":"orphan2"}



fastqc_exe="fastqc"

