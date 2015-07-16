# KneadData User Guide v0.35

Last updated on July 16, 2015.

Authors: Andy Shi, Aleksandar Kostic, Randall Schwager, Lauren McIver, Curtis
Huttenhower
Huttenhower Lab, Harvard School of Public Health,  
Boston, MA

You can access this repository with SSH or with HTTPS.

If you use this software, please cite our paper: (TBA)


## Table of Contents
- [Introduction](#markdown-header-introduction)
- [Tutorial and Demo](#markdown-header-tutorial-and-demo) 
- [Installation](#markdown-header-installation)
- [Quick Start Guide](#markdown-header-quick-start-guide)

    1. Data Locations
    2. Indexing
    3. How to Run

- Detailed Documentation


## Introduction

KneadData is a tool designed to perform quality control on metagenomic and
metatranscriptomic sequencing data, especially data from microbiome experiments.
In these experiments, samples are typically taken from a host in hopes of
learning something about the microbial community on the host. However,
sequencing data from such experiments will often contain a high ratio of host to
bacterial reads. This tool aims to perform principled  *in silico* separation of
bacterial reads from these "contaminant" reads, be they from the host, from
bacterial 16S sequences, or other user-defined sources. Additionally, KneadData
can be used for other filtering tasks. For example, if one is trying to clean
data derived from a human sequencing experiment, KneadData can be used to
separate the human and the non-human reads. 

----------------------------------------------

## Tutorial and Demo

For a brief tutorial and demo, please visit the [KneadData
homepage](http://huttenhower.org/kneaddata)

----------------------------------------------

## Installation

To download the latest stable release, use one of the links below and extract
the files.

- [ZIP](https://bitbucket.org/biobakery/kneaddata/get/v0.35.zip)
- [GZ](https://bitbucket.org/biobakery/kneaddata/get/v0.35.tar.gz)
- [BZ](https://bitbucket.org/biobakery/kneaddata/get/v0.35.tar.bz2)

Or

- `pip install -e 'hg+https://bitbucket.org/biobakery/kneaddata@master#egg=knead_datalib-dev'`

Currently, KneadData is only supported on Linux and Macs. 

KneadData requires
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), and 
[Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). 

If installing via pip, Trimmomatic is optional.

The following dependencies are optional:
[BMTagger](ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/) and [NCBI
BLAST](http://www.ncbi.nlm.nih.gov/books/NBK1762/). 

Please see these respective project websites for download and installation
instructions. 

### Configuration
The Bowtie 2 and/or BMTagger programs need to be able to find their executables:

Bowtie2 executable:
- `bowtie2`

BMTagger executables:
- `srprism`
- `bmfilter`
- `extract_fullseq`
- `blastn` (included in NCBI blast)

There are three ways to do this:

1. Specify an argument (Bowtie2 ONLY). You can specify the path to the Bowtie2
   executable using the `--bowtie2-path` argument. 

2. Update your PATH. If the executables are in your `~/bin` folder, you should
   run the following in your shell:

        PATH=$PATH:~/bin
        export PATH

You may want to consider putting this in your `~/.bashrc` or `~/.bash_profile`
file to make it permanent. 

2. Use a `bmtagger.conf` file (BMTagger ONLY). This file is included with Knead
   Data. You can set the locations of the above executables using the variables
   `SRPRISM`, `BMFILTER`, `EXTRACT_FA`, and `BLASTN` respectively. Make sure
   this file is in your working directory when you run Knead Data, i.e. if you
   are running Knead Data from your home directory, copy this `bmtagger.conf`
   file to your home directory also. 

### Bowtie2 vs. BMTagger

KneadData supports two methods for detecting contaminant reads, either
[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) or
[BMTagger](ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/). Bowtie2 is the
default. In our tests, it uses much less memory compared to BMTagger. Both tools
are similarly accurate, but BMTagger can run faster on larger databases. By
default, all runs are performed with Bowtie2, but it is possible to use BMTagger
as well. 

----------------------------------------------

## Quick Start Guide

### Data Locations

KneadData requires reference sequences for the contamination you are trying to
remove. Let's say you wish to filter reads from a particular "host." Broadly
defined, the host can be an organism, or a set of organisms, or just a set of
sequences. Then, you simply must provide KneadData with a
[FASTA](http://en.wikipedia.org/wiki/FASTA_format) file containing these
sequences. Usually, researchers want to remove reads from the human genome, the
human transcriptome, or ribosomal RNA. You can access some of these FASTA files
using the resources below:

- Ribosomal RNA: [Silva](http://www.arb-silva.de/) provides a comprehensive
  database for ribosomal RNA sequences spanning all three domains of life
  (*Bacteria*, *Archaea*, and *Eukarya*). 

- Human Genome & Transcriptome: Information about the newest assembly of human
  genomic data can be found at the [NCBI project
  page](http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/). USCS
  provides a convenient
  [website](http://hgdownload.cse.ucsc.edu/downloads.html#human) to download
  this data. 


### Generating KneadData Databases

KneadData requires that your reference sequences (FASTA files) be indexed to
form KneadData databases beforehand. This only needs to be done once per
reference sequence. 

For certain common databases, we provide indexed files. If you use these, you
can skip the steps below. 

Links to common databases:

- TBD

#### Bowtie 2

Simply run the `bowtie2-build` indexer included with Bowtie 2 as follows:

`bowtie2-build <reference> <db-name>`

Where `<reference>` is the reference FASTA file, and `<db-name>` is the name you
wish to call your Bowtie2 database. The `bowtie2-build` command must be in your
PATH, or you should provide a full path to the commad. For more details, refer
to the [bowtie2-build
documentation](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer)

#### BMTagger

KneadData includes `generate_db.py`, a Python script that
will automatically generate these databases for BMTagger. Simply run

`python generate_db.py reference.fasta`

By default, this will generate the reference databases, whose names are prefixed
with `reference.fasta`. 

A note on PATH: The above command will fail if the tools in the BMTagger suite
(specifically, bmtool and srprism) and the NCBI BLAST executables are not in
your PATH. If this is the case, you can specify a path to these tools using the
`-b`, `-s`, and `-m` options. Run 

`python generate_db.py -h` 

for more details.


##### Example

Let's say you want to remove human reads from your metagenomic sequencing data.
You've downloaded the human genome in a file called `Homo_sapiens.fasta`. You
can generate the KneadData database by executing

`bowtie2-build Homo_sapiens.fasta -o Homo_sapiens_db`

for Bowtie2, or

`python generate_db.py Homo_sapiens.fasta -o Homo_sapiens_db`

for BMTagger. 

All of the required KneadData database files will have file names prefixed by
`Homo_sapiens_db` and have various file extensions.


### How to Run

After generating your database file, we can start to remove contaminant reads.
As input, KneadData requires FASTQ files. It supports both single end and paired
end reads. KneadData uses either Bowtie 2 (default) or BMTagger to identify the
contaminant reads. 

A note on outputs: KneadData by default outputs new FASTQ files containing only
non-contaminant reads, and new FASTQ files containing only the contaminant
reads. 

When using BMTagger, the default output is a *list* of contaminant reads.
It can also output new FASTQ files containing only the non-contaminant reads,
with the contaminant reads removed. This feature can be enabled by using the
`-x` or the `--extract` flag. 

#### Single End

To run KneadData in single end mode, run 

    ./knead_data.py -1 seq.fastq -db DB_NAME -t TRIM_PATH

This will create files called 

+ `seq.fastq_DB_NAME_clean.fastq`: FASTQ file containing reads that were not
  identified as contaminants. 
+ `seq.fastq_DB_NAME_contam.fastq`: FASTQ file containing reads that were
  identified as contaminants.
+ `seq.fastq_output.fastq`: Same contents as `seq.fastq_DB_NAME_clean.fastq`.
  KneadData with Bowtie2 will produce one clean/contaminated FASTQ output pair
  per reference database. Then, the intersection of the reads in all the clean
  FASTQs will be output in this file. This represents reads that were not in any
  of the reference databases. Since we only searched against one reference
  database, this file is a copy of `seq.fastq_DB_NAME_clean.fastq`.

To run KneadData in single end mode with BMTagger, run

    ./knead_data.py -1 seq.fastq -db DB_NAME -t TRIM_PATH --bmtagger --bmtagger-path BMTAGGER_PATH

By default, this will create a file called `seq.fastq_output.out` containing a
list of contaminant reads found in `seq.fastq`. If instead you pass the
`--extract` option to `knead_data.py`, you will get a file called
`seq.fastq_output.fastq` which contains all the non-contaminant reads found in
`seq.fastq`.

+ `seq.fastq`: Your input FASTQ file
+ `DB_NAME`: Prefix for the KneadData database. 
+ `TRIM_PATH`: Path to the Trimmomatic executable.
+ `BMTAGGER_PATH`: Path to the BMTagger executable. If not specified, the
  program will attempt to find `bmtagger.sh` in your `$PATH`. 


##### Example

Continuing from the previous example, let's say we want to remove human reads
from a file called `seq.fastq` using the *Homo sapiens* database we generated
earlier. Additionally, suppose the Trimmomatic executable file was located at
`~/bin/Trimmomatic/trimmomatic-0.32.jar`. To run with Bowtie2:

    ./knead_data.py -1 seq.fastq -db Homo_sapiens_db -t
    ~/bin/Trimmomatic/trimmomatic-0.32.jar -o seq_output

This will create files called 

+ `seq_output_Homo_sapiens_db_clean.fastq`: FASTQ file containing reads that
  were not identified as being human.
+ `seq_output_Homo_sapiens_db_contam.fastq`: FASTQ file containing reads that
  were identified as human reads. 
+ `seq_output.fastq`: Same contents as `seq.fastq_Homo_sapiens_db_clean.fastq`.
  Contains reads that were not identified as being human. 

If you wanted to use BMTagger, suppose the BMTagger executable was located
at `~/bin/bmtagger.sh`. Let's say you want your contaminant reads to be stored
in a file called `seq_contams.out`. You would then run

    ./knead_data.py -1 seq.fastq -db DB_NAME -t TRIM_PATH --bmtagger --bmtagger-path BMTAGGER_PATH -o seq_output

Let's say that, instead of outputting your contaminant reads in a separate file,
you just want a "cleaned" FASTQ file that contains no contaminant reads. If you
execute

    ./knead_data.py -1 seq.fastq -db Homo_sapiens_db -t ~/bin/Trimmomatic/trimmomatic-0.32.jar --bmtagger --bmtagger-path ~/bin/bmtagger.sh -o seq_clean

you will get a file `seq_clean.fastq` which contains all the non-contaminant
reads, the ones that were not identified as human reads. 


#### Paired End

To run KneadData in paired end mode with Bowtie2, run

`./knead_data.py -1 seq1.fastq -2 seq2.fastq -db DB_NAME -t TRIM_PATH -o
seq_output`

To run KneadData in paired end mode with BMTagger, run

    python knead_data.py -1 seq1.fastq -2 seq2.fastq -db DB_NAME -t TRIM_PATH --bmtagger --bmtagger-path BMTAGGER_PATH -o seq_output

+ `seq1.fastq`: Your input FASTQ file, first mate
+ `seq2.fastq`: Your input FASTQ file, second mate
+ `DB_NAME`: Prefix for the KneadData database. 
+ `TRIM_PATH`: Path to the Trimmomatic executable.
+ `BMTAGGER_PATH`: Path to the BMTagger executable.
+ `seq_output`: Optional, but recommended. Prefix for the output files. 

The outputs depend on what happens during the quality filtering and trimming
part of the pipeline. 

When performing quality filtering and trimming for paired end files, 3 things
can happen:

1. Both reads in the pair pass. 
2. The read in the first mate passes, and the one in the second doesn't. 
3. The read in the second mate passes, and the one in the first doesn't. 

The number of outputs are a function of the read quality. 

KneadData + Bowtie2 Outputs: There can be up to 8 outputs per reference
database, plus up to 5 aggregate outputs. 

KneadData + BMTagger Outputs: There can be 1-3 outputs per reference database if
you run KneadData in paired end mode and output the list of contaminant reads,
and 2-4 outputs per reference database if you run KneadData in paired end mode
and want cleaned FASTQ files without the contaminant reads. 


##### Example

Instead of single end reads, let's say you have paired end reads and you want to
separate the reads that came from bacterial mRNA, bacterial rRNA, and human RNA.
You have two databases, one prefixed `bact_rrna_db` and the other prefixed
`human_rna_db`, and your sequence files are `seq1.fastq` and `seq2.fastq`. To
run with Bowtie2, execute

    ./knead_data.py -1 seq1.fastq -2 seq2.fastq -db bact_rrna_db human_rna_db -t ~/bin/Trimmomatic/trimmomatic-0.32.jar -o seq_out

This will output:

Files for just the `bact_rrna_db` database:

+ `seq_out_pe_bact_rrna_db_contam_1.fastq`: Reads from the first mate in
  situation (1) above that were identified as belonging to the `bact_rrna_db`
  database. 
+ `seq_out_pe_bact_rrna_db_contam_2.fastq`: Reads from the second mate in
  situation (1) above that were identified as belonging to the `bact_rrna_db`
  database. 
+ `seq_out_pe_bact_rrna_db_clean_1.fastq`: Reads from the first mate in
  situation (1) above that were identified as NOT belonging to the
  `bact_rrna_db` database. 
+ `seq_out_pe_bact_rrna_db_clean_2.fastq`: Reads from the second mate in
  situation (1) above that were identified as NOT belonging to the
  `bact_rrna_db` database. 

Depending on the input FASTQ, one or more of the following may be output:

+ `seq_out_se_1_bact_rrna_db_contam.fastq`: Reads from the first mate in
  situation (2) above that were identified as belonging to the `bact_rrna_db`
  database. 
+ `seq_out_se_1_bact_rrna_db_clean.fastq`: Reads from the first mate in
  situation (2) above that were identified as NOT belonging to the
  `bact_rrna_db` database. 
+ `seq_out_se_2_bact_rrna_db_contam.fastq`: Reads from the second mate in
  situation (3) above that were identified as belonging to the `bact_rrna_db`
  database. 
+ `seq_out_se_2_bact_rrna_db_clean.fastq`: Reads from the second mate in
  situation (3) above that were identified as NOT belonging to the
  `bact_rrna_db` database. 

Files for just the `human_rna_db` database:

+ `seq_out_pe_human_rna_db_contam_1.fastq`: Reads from the first mate in
  situation (1) above that were identified as belonging to the `human_rna_db`
  database. 
+ `seq_out_pe_human_rna_db_contam_2.fastq`: Reads from the second mate in
  situation (1) above that were identified as belonging to the `human_rna_db`
  database. 
+ `seq_out_pe_human_rna_db_clean_1.fastq`: Reads from the first mate in
  situation (1) above that were identified as NOT belonging to the
  `human_rna_db` database. 
+ `seq_out_pe_human_rna_db_clean_2.fastq`: Reads from the second mate in
  situation (1) above that were identified as NOT belonging to the
  `human_rna_db` database. 

Depending on the input FASTQ, one or more of the following may be output:

+ `seq_out_se_1_human_rna_db_contam.fastq`: Reads from the first mate in
  situation (2) above that were identified as belonging to the `human_rna_db`
  database. 
+ `seq_out_se_1_human_rna_db_clean.fastq`: Reads from the first mate in
  situation (2) above that were identified as NOT belonging to the
  `human_rna_db` database. 
+ `seq_out_se_2_human_rna_db_contam.fastq`: Reads from the second mate in
  situation (2) above that were identified as belonging to the `human_rna_db`
  database. 
+ `seq_out_se_2_human_rna_db_clean.fastq`: Reads from the second mate in
  situation (2) above that were identified as NOT belonging to the
  `human_rna_db` database. 

Aggregated files:

+ `seq_out.log`: Log file containing statistics about the run. 
+ `seq_out_pe_1.fastq`: Reads from the first mate in situation (1) identified as
  NOT belonging to any of the reference databases. 
+ `seq_out_pe_2.fastq`: Reads from the second mate in situation (1) identified as
  NOT belonging to any of the reference databases. 
+ `seq_out_se_1.fastq`: Reads from the first mate in situation (2) identified as
  NOT belonging to any of the reference databases. 
+ `seq_out_se_2.fastq`: Reads from the second mate in situation (3) identified as
  NOT belonging to any of the reference databases. 


To run with BMTagger, execute

    ./knead_data.py -1 seq1.fastq -2 seq2.fastq -db bact_rrna_db human_rna_db -t ~/bin/Trimmomatic/trimmomatic-0.32.jar --bmtagger --bmtagger-path ~/bin/bmtagger.sh -o seq_contams

---------------------------------

## Detailed Documentation

This documentation can be accessed via `./knead_data.py -h`.

```
#!text
usage: knead_data.py [-h] -1 INFILE1 [-2 INFILE2] -db REFERENCE_DB
                     [-o OUTPUT_PREFIX] [-D OUTPUT_DIR] [--threads THREADS]
                     [-s {memory,storage}] [-l LOGGING] [--logfile LOGFILE]
                     [--version] [-t TRIM_PATH] [--trimlen TRIMLEN]
                     [-m MAX_MEM] [-a TRIM_ARGS] [--bowtie2-path BOWTIE2_PATH]
                     [--bowtie2-args BOWTIE2_ARGS] [--bmtagger] [--extract]
                     [--bmtagger-path BMTAGGER_PATH] [--trf]
                     [--trf-path TRF_PATH] [--match MATCH]
                     [--mismatch MISMATCH] [--delta DELTA] [--pm PM] [--pi PI]
                     [--minscore MINSCORE] [--maxperiod MAXPERIOD]
                     [--no-generate-fastq] [--mask] [--html]

optional arguments:
  -h, --help            show this help message and exit

global options:
  -1 INFILE1, --infile1 INFILE1
                        input FASTQ file
  -2 INFILE2, --infile2 INFILE2
                        input FASTQ file mate
  -db REFERENCE_DB, --reference-db REFERENCE_DB
                        prefix for reference databases used for either Bowtie2
                        or BMTagger
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        prefix for all output files
  -D OUTPUT_DIR, --output-dir OUTPUT_DIR
                        where to put all output files
  --threads THREADS     Maximum number of processes to run. Default uses all
                        but one available CPUs
  -s {memory,storage}, --strategy {memory,storage}
                        Define operating strategy: 'storage' for IO-heavy or
                        'memory' for memory-heavy
  -l LOGGING, --logging LOGGING
                        Logging verbosity, options are debug, info, warning,
                        and critical. If set to debug, temporary files are not
                        removed
  --logfile LOGFILE     Where to save logs
  --version             Print version and exit

trimmomatic arguments:
  -t TRIM_PATH, --trim-path TRIM_PATH
                        path to Trimmomatic .jar executable
  --trimlen TRIMLEN     minimum length for a trimmed read in Trimmomatic
  -m MAX_MEM, --max-mem MAX_MEM
                        Maximum amount of memory that will be used by
                        Trimmomatic, as a string, ie 500m or 8g
  -a TRIM_ARGS, --trim-args TRIM_ARGS
                        Additional arguments for Trimmomatic, default:
                        SLIDINGWINDOW:4:20 MINLEN:60

bowtie2 arguments:
  --bowtie2-path BOWTIE2_PATH
                        path to bowtie2 if not found on $PATH
  --bowtie2-args BOWTIE2_ARGS
                        Additional arguments for Bowtie 2, default: --very-
                        sensitive

bmtagger arguments:
  --bmtagger            If set, use BMTagger to identify contaminant reads
  --extract             Only has an effect if --bmtagger is set. If this is
                        set, knead_data outputs cleaned FASTQs, without
                        contaminant reads. Else, output a list or lists of
                        contaminant reads.
  --bmtagger-path BMTAGGER_PATH
                        path to BMTagger executable if not found in $PATH

trf arguments:
  --trf                 If set, use TRF to remove and/or mask tandem repeats
  --trf-path TRF_PATH   Path to TRF executable if not found in $PATH
  --match MATCH         TRF matching weight. Default: 2
  --mismatch MISMATCH   TRF mismatching penalty. Default: 7
  --delta DELTA         TRF indel penalty. Default: 7
  --pm PM               TRF match probability (whole number). Default: 80
  --pi PI               TRF indel probability (whole number). Default: 10
  --minscore MINSCORE   TRF minimum alignment score to report. Default: 50
  --maxperiod MAXPERIOD
                        TRF maximum period size to report. Default: 500
  --no-generate-fastq   If switched on, don't generate fastq output for trf
  --mask                If switched on, generate mask file for trf output
  --html                If switched on, generate html file for trf output
```

----------------------------------

## Usage Notes

### Specifying Additional Arguments

If you want to specify additional arguments for Bowtie2 using the
`--bowtie2-args` flag, you will need to use the equals sign. For each additional
argument, you should specify the `--bowtie2-args` flag. 

For example:

`./knead_data.py ... --bowtie2-args=--very-fast --bowtie2-args=-p 2`

A similar approach is used to specify additional arguments for Trimmomatic:

`./knead_data.py ... --trim-args "LEADING:3" --trim-args "TRAILING:3"`

*NOTE*: Manually specifying additional arguments will completely override the
defaults. 

### Memory

KneadData requires quite a bit of memory when running with BMTagger, around 8-9
gigabytes. When running with Bowtie2, KneadData only requires around 2-4
gigabytes.
