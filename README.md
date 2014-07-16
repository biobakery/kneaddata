KneadData User Guide v0.1
========================

Last updated on July 16 2014.

Authors: Andy Shi and Aleksandar Kostic  
Huttenhower Lab, Harvard School of Public Health,  
Boston, MA

You can access this repository with SSH or with HTTPS.

Table of Contents
-----------------
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Quick Start Guide](#quick-start-guide)

    1. [Data Locations](#data-locations)
    2. [Indexing](#indexing)
    3. [How to Run](#how-to-run)

4. [Detailed Documentation](#detailed-documentation)

# 1. Introduction

KneadData is a tool designed to perform quality control on metagenomic
sequencing data, especially data from microbiome experiments. In these
experiments, samples are typically taken from a host in hopes of learning
something about the microbial community on the host. However, metagenomic
sequencing data from such experiments will often contain a high ratio of host to
bacterial reads. This tool aims to perform principled  *in silico* removal of
these "contaminant" reads, be they from the host, from bacterial 16S sequences,
or other user-defined sources.

# 2. Installation

To download the latest stable release, use one of the links below and extract
the files.

+ [ZIP](https://bitbucket.org/biobakery/kneaddata/get/v0.1.zip)
+ [GZ](https://bitbucket.org/biobakery/kneaddata/get/v0.1.tar.gz)
+ [BZ](https://bitbucket.org/biobakery/kneaddata/get/v0.1.tar.bz2)

Currently, KneadData is only supported on Linux and Macs. 

KneadData requires
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic),
[BMTagger](ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/), and [NCBI
BLAST](http://www.ncbi.nlm.nih.gov/books/NBK1762/). Please see these respective
project websites for download and installation instructions. 

# 3. Quick Start Guide

### 3.1. Data Locations

KneadData requires reference sequences for the contamination you are trying to
remove. If you wish to remove reads from a particular "host" (broadly defined,
the host can be an organism, or a set of organisms, or just a set of sequences),
you simply must provide KneadData with a
[FASTA](http://en.wikipedia.org/wiki/FASTA_format) file containing these
sequences. Usually, researchers want to remove reads from the human genome, the
human transcriptome, or ribosomal RNA. You can access some of these FASTA files
using the resources below:

+ Ribosomal RNA: [Silva](http://www.arb-silva.de/) provides a comprehensive
  database for ribosomal RNA sequences spanning all three domains of life
  (*Bacteria*, *Archaea*, and *Eukarya*). 

+ Human Genome & Transcriptome: Information about the newest assembly of human
  genomic data can be found at the [NCBI project
  page](http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/). USCS
  provides a convenient
  [website](http://hgdownload.cse.ucsc.edu/downloads.html#human) to download
  this data. 

### 3.2. Generating KneadData Databases

KneadData requires that your reference sequences (FASTA files) be indexed to
form KneadData databases beforehand. This only needs to be done once per
reference sequence. KneadData includes `generate_db.py`, a Python script that
will automatically generate these databases. Simply run

`python generate_db.py reference.fasta`

By default, this will generate the reference databases, whose names are prefixed
with `reference.fasta`. 

A note on PATH: The above command will fail if the tools in the BMTagger suite
(specifically, bmtool and srprism) and the NCBI BLAST executables are not in
your PATH. If this is the case, you can specify a path to these tools. Run 

`python generate_db.py -h` 


##### Example

Let's say you want to remove human reads from your metagenomic sequencing data.
You've downloaded the human genome in a file called `Homo_sapiens.fasta`. You
can generate the KneadData database by executing

`python generate_db.py Homo_sapiens.fasta -o Homo_sapiens_db`

All of the required KneadData database files will have file names prefixed by
`Homo_sapiens_db` and have various file extensions.

or see the [Detailed Documentation](#detailed-documentation) for more information.

### 3.3. How to Run

After generating your database file, we can start to remove contaminant reads.
As input, KneadData requires FASTQ files. It supports both single end and paired
end reads. 

A note on outputs: KneadData by default outputs a *list* of contaminant reads.
It can also output new FASTQ files containing only the non-contaminant reads,
with the contaminant reads removed. This feature can be enabled by using the
`-x` or the `--extract` flag. 

#### Single End

To run KneadData in single end mode, run

`python knead_data.py -1 seq.fastq -db DB_NAME -t TRIM_PATH -b BMTAGGER_PATH`

By default, this will create a file called `seq.fastq_output.out` containing a
list of contaminant reads found in `seq.fastq`. 

+ `seq.fastq`: Your input FASTQ file
+ `DB_NAME`: Prefix for the KneadData database. 
+ `TRIM_PATH`: Path to the Trimmomatic executable.
+ `BMTAGGER_PATH`: Path to the BMTagger executable.


##### Example

Continuing from the previous example, let's say we want to remove human reads
from a file called `seq.fastq` using the *Homo sapiens* database we generated
earlier. Additionally, suppose the Trimmomatic executable file was located at
`~/bin/Trimmomatic/trimmomatic-0.32.jar` and the BMTagger executable was located
at `~/bin/bmtagger.sh`. Let's say you want your contaminant reads to be stored
in a file called `seq_contams.out`. You would then run

`python knead_data.py -1 seq.fastq -db Homo_sapiens_db -t
~/bin/Trimmomatic/trimmomatic-0.32.jar -b ~/bin/bmtagger.sh -o seq_contams` 

Let's say that, instead of outputting your contaminant reads in a separate file,
you just want a "cleaned" FASTQ file that contains no contaminant reads. If you
execute

`python knead_data.py -1 seq.fastq -db Homo_sapiens_db -t
~/bin/Trimmomatic/trimmomatic-0.32.jar -b ~/bin/bmtagger.sh -o seq_clean` 

you will get a file `seq_clean.fastq` which contains all the non-contaminant
reads--in this case, the ones that weren't human reads. 

#### Paired End
To run KneadData in paired end mode, run

`python knead_data.py -1 seq1.fastq -2 seq2.fastq -db DB_NAME -t TRIM_PATH -b BMTAGGER_PATH -o seq_output`

+ `seq1.fastq`: Your input FASTQ file, first mate
+ `seq2.fastq`: Your input FASTQ file, second mate
+ `DB_NAME`: Prefix for the KneadData database. 
+ `TRIM_PATH`: Path to the Trimmomatic executable.
+ `BMTAGGER_PATH`: Path to the BMTagger executable.
+ `seq_output`: Optional, but recommended. Prefix for the output files. 

Outputs: There can be 1-3 outputs if you run KneadData in paired end mode and
output the list of contaminant reads, and 2-4 outputs if you run KneadData in
paired end mode and want cleaned FASTQ files without the contaminant reads. This
all depends on what happens during the quality filtering and trimming part of
the pipeline. 

When performing quality filtering and trimming for paired end files, 3 things
can happen:

1. Both reads in the pair pass. 
2. The read in the first mate passes, and the one in the second doesn't. 
3. The read in the second mate passes, and the one in the first doesn't. 

The different number of outputs are a function of the read quality. See the
example below for more concreteness.

##### Example

Instead of single end reads, let's say you have paired end reads and you again
want to identify all the human reads using the database generated earlier. If your
files are `seq1.fastq` and `seq2.fastq`, you would execute

`python knead_data.py -1 seq1.fastq -2 seq2.fastq -db Homo_sapiens_db -t
~/bin/Trimmomatic/trimmomatic-0.32.jar -b ~/bin/bmtagger.sh -o seq_contams` 

This will output:

+ `seq_contams.out`: A list of reads that were identified as human from situation
1 above. 

Depending on the input, one or more of the files below may not be output. 

+ `seq_contams_se_1.out`: A list of reads that were identified as human from
situation 2 above. 
+ `seq_contams_se_2.out`: A list of reads that were identified as human from
situation 3 above. 

If instead you want cleaned FASTQ files that have the human reads removed, you
would execute 

`python knead_data.py -1 seq1.fastq -2 seq2.fastq -db Homo_sapiens_db -t
~/bin/Trimmomatic/trimmomatic-0.32.jar -b ~/bin/bmtagger.sh -o seq_clean
--extract`

This will output:

+ `seq_clean_pe_1.fastq` & `seq_clean_pe_2.fastq`: These files correspond to
situation 1 described above. These files contain the cleaned reads where both
mates passed the quality filtering process. 

Depending on the input, one or more of the files below may not be output. 

+ `seq_clean_se_1.fastq`: This file corresponds to situation 2 described above.
  It contains cleaned reads where only the first mate passed the filtering
  process.
  
+ `seq_clean_se_2.fastq`: This file corresponds to situation 3 described above.
  It contains cleaned reads where only the second mate passed the filtering
  process.


### A Note on Memory

KneadData requires quite a bit of memory, around 8-9 gigabytes, to run. 


# 4. Detailed Documentation

    usage: knead_data.py [-h] -1 INFILE1 [-2 INFILE2] [--trimlen TRIMLEN]
                        [-o OUTPUT_PREFIX] [-db REFERENCE_DB [REFERENCE_DB ...]]
                        -t TRIM_PATH -b BMTAGGER_PATH [-x] [-m MAX_MEM]
                        [-a TRIM_ARGS] [-d]

    optional arguments:
    -h, --help            show this help message and exit
    -1 INFILE1, --infile1 INFILE1
                            input FASTQ file
    -2 INFILE2, --infile2 INFILE2
                            input FASTQ file mate
    --trimlen TRIMLEN     length to trim reads
    -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                            prefix for all output files
    -db REFERENCE_DB [REFERENCE_DB ...], --reference-db REFERENCE_DB [REFERENCE_DB ...]
                            prefix for reference databases used in BMTagger
    -t TRIM_PATH, --trim-path TRIM_PATH
                            path to Trimmomatic
    -b BMTAGGER_PATH, --bmtagger-path BMTAGGER_PATH
                            path to BMTagger
    -x, --extract         Remove contaminant reads
    -m MAX_MEM, --max-mem MAX_MEM
                            Maximum amount of memory that will be used by
                            Trimmomatic, as a string, ie 500m or 8g
    -a TRIM_ARGS, --trim-args TRIM_ARGS
                            additional arguments for Trimmomatic
    -d, --debug           If set, temporary files are not removed

