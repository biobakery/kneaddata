# KneadData User Guide v0.4.6

Last updated on September 10, 2015.

Authors: Andy Shi, Aleksandar Kostic, Randall Schwager, Lauren McIver, Curtis
Huttenhower
Huttenhower Lab, Harvard School of Public Health,  
Boston, MA

You can access this repository with SSH or with HTTPS.

If you use this software, please cite our paper: (TBA)


## Table of Contents
- [Introduction](#markdown-header-introduction)
- [Tutorial and Demo](#markdown-header-tutorial-and-demo) 
- [Requirements](#markdown-header-requirements)
- [Installation](#markdown-header-installation)
- [How to Run](#markdown-header-how-to-run)
- [Complete Option List](#markdown-header-complete-option-list)


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

## Requirements

1.  [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (version == 0.33) (automatically installed)
2.  [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (version >= 2.1) (automatically installed)
3.  [Python](http://www.python.org/) (version >= 2.7)
4.  [Java Runtime Environment](http://www.oracle.com/technetwork/java/javase/downloads/jre7-downloads-1880261.html)
5.  Memory (>= 4 Gb if using Bowtie2, >= 8 Gb if using BMTagger)
6.  Operating system (Linux or Mac)

Optionally, [BMTagger](ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/) can be used instead of Bowtie2. Bowtie2 is run by default. In our tests, we found bowtie2 uses less memory and BMTagger runs faster on larger databases.

The executables for the required software packages should be installed in your $PATH. Alternatively, you can provide the location of the Bowtie2 install ($BOWTIE2_DIR) with the following KneadData option “--bowtie2 $BOWTIE2_DIR”.

## Installation

### Download KneadData

You can download the latest KneadData release or the development version.

Option 1: Latest Release (Recommended)

* [Download](https://bitbucket.org/biobakery/kneaddata/downloads/kneaddata-v0.4.6.1.tar.gz) and unpack the latest release of KneadData.

Option 2: Development Version

* Create a clone of the repository:

    `` $ hg clone https://bitbucket.org/biobakery/kneaddata ``

    Note: Creating a clone of the repository requires [Mercurial](http://mercurial.selenic.com/) to be installed. Once the clone is created you can always update to the latest version of the repository with `` $ hg pull --update ``.

### Install KneadData

Before installing KneadData, please install the Java Runtime Environment (JRE). First [download](http://www.oracle.com/technetwork/java/javase/downloads/jre7-downloads-1880261.html) the JRE for your platform. Then follow the instructions for your platform: [Linux 64-bit](http://docs.oracle.com/javase/8/docs/technotes/guides/install/linux_jre.html#CFHIEGAA) or [Mac OS](http://docs.oracle.com/javase/8/docs/technotes/guides/install/mac_jre.html#jre_8u40_osx). At the end of the installation, add the location of the java executable to your $PATH.

1. Download and unpack the KneadData software
    * Download the software: [kneaddata.tar.gz](https://bitbucket.org/biobakery/kneaddata/downloads/kneaddata-v0.4.6.1.tar.gz)
    * `` $ tar zxvf kneaddata.tar.gz ``
    * `` $ cd kneaddata ``
2. From the KneadData directory, install KneadData
    * `` $ python setup.py install ``
    * This command will automatically install Trimmomatic and Bowtie2. To bypass the install of dependencies, add the option "--bypass-dependencies-install".
    * If you do not have write permissions to '/usr/lib/', then add the option "--user" to the install command. This will install the python package into subdirectories of '~/.local'. Please note when using the "--user" install option on some platforms, you might need to add '~/.local/bin/' to your $PATH as it might not be included by default. You will know if it needs to be added if you see the following message ``kneaddata: command not found`` when trying to run KneadData after installing with the "--user" option.
3. Download the human reference database to $DIR
    * `` $ kneaddata_database --download human bowtie2 $DIR ``
    * When running this command, $DIR should be replaced with the full path to the directory you have selected to store the database.

### Creating a Custom Reference Database

A human reference database can be downloaded to use when running KneadData. Alternatively, you can create your own custom reference database.

#### Selecting Reference Sequences

First you must select reference sequences for the contamination you are trying to
remove. Say you wish to filter reads from a particular "host." Broadly
defined, the host can be an organism, or a set of organisms, or just a set of
sequences. Then, you simply must generate a reference database for KneadData from a
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


#### Generating KneadData Databases

KneadData requires that your reference sequences (FASTA files) be indexed to
form KneadData databases beforehand. This only needs to be done once per
reference sequence. 

For certain common databases, we provide indexed files. If you use these, you
can skip the manual build steps below. 

To download the indexed human reference database, run the following command:
    * `` $ kneaddata_database --download human bowtie2 $DIR ``
    * When running this command, $DIR should be replaced with the full path to the directory you have selected to store the database.

##### Creating a Bowtie2 Database

Simply run the `bowtie2-build` indexer included with Bowtie2 as follows:

`$ bowtie2-build <reference> <db-name>`

Where `<reference>` is the reference FASTA file, and `<db-name>` is the name you
wish to call your Bowtie2 database. For more details, refer
to the [bowtie2-build
documentation](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer)

##### Creating a BMTagger Database

KneadData includes `kneaddata_build_database`, an executable that
will automatically generate these databases for BMTagger. Simply run

`$ kneaddata_build_database reference.fasta`

By default, this will generate the reference databases, whose names are prefixed
with `reference.fasta`. 

A note on PATH: The above command will fail if the tools in the BMTagger suite
(specifically, bmtool and srprism) and the NCBI BLAST executables are not in
your PATH. If this is the case, you can specify a path to these tools using the
`-b`, `-s`, and `-m` options. Run 

`$ kneaddata_build_database --help` 

for more details.


##### Example Custom Database

Say you want to remove human reads from your metagenomic sequencing data.
You downloaded the human genome in a file called `Homo_sapiens.fasta`. You
can generate the KneadData database by executing

`$ bowtie2-build Homo_sapiens.fasta -o Homo_sapiens_db`

for Bowtie2, or

`$ kneaddata_build_database Homo_sapiens.fasta -o Homo_sapiens_db`

for BMTagger. 

All of the required KneadData database files will have file names prefixed by
`Homo_sapiens_db` and have various file extensions.


### How to Run

After downloading or generating your database file, you can start to remove contaminant reads.
As input, KneadData requires FASTQ files. It supports both single end and paired
end reads. KneadData uses either Bowtie2 (default) or BMTagger to identify the
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

    ` $ kneaddata --input seq.fastq --reference-db $DATABASE --output kneaddata_output `

This will create files in the folder `kneaddata_output` named 

+ `seq_kneaddata_$DATABASE_contam.fastq`: FASTQ file containing reads that were
  identified as contaminants from the database (named $DATABASE).
+ `seq_kneaddata.fastq`: This represents reads that were not in the reference databases.
+ `seq_kneaddata.log`

To run KneadData in single end mode with BMTagger, run

    ` $ kneaddata --input seq.fastq --reference-db $DATABASE --run-bmtagger`

By default, this will create a file called `seq_kneaddata_output.out` containing a
list of contaminant reads found in `seq.fastq`. If instead you pass the
`--extract` option to `kneaddata`, you will get a file called
`seq_kneaddata.fastq` which contains all the non-contaminant reads found in
`seq.fastq`.

+ `seq.fastq`: Your input FASTQ file
+ `$DATABASE`: Prefix for the KneadData database. 

##### Example

Continuing from the previous example, say we want to remove human reads from a
file called `seq.fastq` using the *Homo sapiens* database we generated earlier.
To run with Bowtie2:

    ` $ kneaddata --input seq.fastq --reference-db Homo_sapiens_db --output seq_output `

This will create files in the folder `seq_output` named:

+ `seq_kneaddata_Homo_sapiens_db_contam.fastq`: FASTQ file containing reads that
  were identified as human reads. 
+ `seq_kneaddata.fastq`: Contains reads that were not identified as being human. 
+ `seq_kneaddata.log`

If you wanted to use BMTagger, suppose the BMTagger executable was located at
`~/bin/bmtagger.sh`. Say you want your contaminant reads to be stored in a file
called `seq_contams.out` in the folder `kneaddata_output`. You would then run

    ` $ kneaddata --input seq.fastq -db Homo_sapiens_db --run-bmtagger --bmtagger ~/bin/bmtagger.sh --output-prefix seq_output --output kneaddata_output `

Say that, instead of outputting your contaminant reads in a separate file, you
just want a "cleaned" FASTQ file that contains no contaminant reads. If you
execute

    ` $ kneaddata --input seq.fastq -db Homo_sapiens_db --run-bmtagger --bmtagger ~/bin/bmtagger.sh --output kneaddata_output`

you will get a file `seq_kneaddata.fastq` which contains all the non-contaminant
reads, the ones that were not identified as human reads. 


#### Paired End

To run KneadData in paired end mode with Bowtie2, run

    ` $ kneaddata --input seq1.fastq --input seq2.fastq -db $DATABASE --output kneaddata_output`

To run KneadData in paired end mode with BMTagger, run

    ` $ kneaddata --input seq1.fastq --input seq2.fastq -db $DATABASE --run-bmtagger --output kneaddata_output `

+ `seq1.fastq`: Your input FASTQ file, first mate
+ `seq2.fastq`: Your input FASTQ file, second mate
+ `$DATABASE`: Prefix for the KneadData database. 
+ `kneaddata_output`: The folder to write the output files. 

The outputs depend on what happens during the quality filtering and trimming
part of the pipeline. 

When performing quality filtering and trimming for paired end files, 3 things
can happen:

1. Both reads in the pair pass. 
2. The read in the first mate passes, and the one in the second does not pass. 
3. The read in the second mate passes, and the one in the first does not pass. 

The number of outputs are a function of the read quality. 

KneadData + Bowtie2 Outputs: There can be up to 8 outputs per reference
database, plus up to 5 aggregate outputs. 

KneadData + BMTagger Outputs: There can be 1-3 outputs per reference database if
you run KneadData in paired end mode and output the list of contaminant reads,
and 2-4 outputs per reference database if you run KneadData in paired end mode
and want cleaned FASTQ files without the contaminant reads. 


##### Example

Instead of single end reads, say you have paired end reads and you want to
separate the reads that came from bacterial mRNA, bacterial rRNA, and human RNA.
You have two databases, one prefixed `bact_rrna_db` and the other prefixed
`human_rna_db`, and your sequence files are `seq1.fastq` and `seq2.fastq`. To
run with Bowtie2, execute

    `$ kneaddata --input seq1.fastq --input seq2.fastq -db bact_rrna_db -db human_rna_db --output seq_out `

This will output files in the folder `seq_out` named:

Files for just the `bact_rrna_db` database:

+ `seq_kneaddata_pe_bact_rrna_db_contam_1.fastq`: Reads from the first mate in
  situation (1) above that were identified as belonging to the `bact_rrna_db`
  database. 
+ `seq_kneaddata_pe_bact_rrna_db_contam_2.fastq`: Reads from the second mate in
  situation (1) above that were identified as belonging to the `bact_rrna_db`
  database. 
+ `seq_kneaddata_pe_bact_rrna_db_clean_1.fastq`: Reads from the first mate in
  situation (1) above that were identified as NOT belonging to the
  `bact_rrna_db` database. 
+ `seq_kneaddata_pe_bact_rrna_db_clean_2.fastq`: Reads from the second mate in
  situation (1) above that were identified as NOT belonging to the
  `bact_rrna_db` database. 

Depending on the input FASTQ, one or more of the following may be output:

+ `seq_kneaddata_se_1_bact_rrna_db_contam.fastq`: Reads from the first mate in
  situation (2) above that were identified as belonging to the `bact_rrna_db`
  database. 
+ `seq_kneaddata_se_1_bact_rrna_db_clean.fastq`: Reads from the first mate in
  situation (2) above that were identified as NOT belonging to the
  `bact_rrna_db` database. 
+ `seq_kneaddata_se_2_bact_rrna_db_contam.fastq`: Reads from the second mate in
  situation (3) above that were identified as belonging to the `bact_rrna_db`
  database. 
+ `seq_kneaddata_se_2_bact_rrna_db_clean.fastq`: Reads from the second mate in
  situation (3) above that were identified as NOT belonging to the
  `bact_rrna_db` database. 

Files for just the `human_rna_db` database:

+ `seq_kneaddata_pe_human_rna_db_contam_1.fastq`: Reads from the first mate in
  situation (1) above that were identified as belonging to the `human_rna_db`
  database. 
+ `seq_kneaddata_pe_human_rna_db_contam_2.fastq`: Reads from the second mate in
  situation (1) above that were identified as belonging to the `human_rna_db`
  database. 
+ `seq_kneaddata_pe_human_rna_db_clean_1.fastq`: Reads from the first mate in
  situation (1) above that were identified as NOT belonging to the
  `human_rna_db` database. 
+ `seq_kneaddata_pe_human_rna_db_clean_2.fastq`: Reads from the second mate in
  situation (1) above that were identified as NOT belonging to the
  `human_rna_db` database. 

Depending on the input FASTQ, one or more of the following may be output:

+ `seq_kneaddata_se_1_human_rna_db_contam.fastq`: Reads from the first mate in
  situation (2) above that were identified as belonging to the `human_rna_db`
  database. 
+ `seq_kneaddata_se_1_human_rna_db_clean.fastq`: Reads from the first mate in
  situation (2) above that were identified as NOT belonging to the
  `human_rna_db` database. 
+ `seq_kneaddata_se_2_human_rna_db_contam.fastq`: Reads from the second mate in
  situation (2) above that were identified as belonging to the `human_rna_db`
  database. 
+ `seq_kneaddata_se_2_human_rna_db_clean.fastq`: Reads from the second mate in
  situation (2) above that were identified as NOT belonging to the
  `human_rna_db` database. 

Aggregated files:

+ `seq_kneaddata.log`: Log file containing statistics about the run. 
+ `seq_kneaddata_pe_1.fastq`: Reads from the first mate in situation (1) identified as
  NOT belonging to any of the reference databases. 
+ `seq_kneaddata_pe_2.fastq`: Reads from the second mate in situation (1) identified as
  NOT belonging to any of the reference databases. 
+ `seq_kneaddata_se_1.fastq`: Reads from the first mate in situation (2) identified as
  NOT belonging to any of the reference databases. 
+ `seq_kneaddata_se_2.fastq`: Reads from the second mate in situation (3) identified as
  NOT belonging to any of the reference databases. 

---------------------------------

## Complete Option List

All options can be accessed with `$ kneaddata --help`.

```
usage: kneaddata [-h] [--version] [-v] -i INPUT -o OUTPUT_DIR
                 [-db REFERENCE_DB] [--output-prefix OUTPUT_PREFIX] [-t <1>]
                 [-p <1>] [-q {phred33,phred64}] [-s {memory,storage}]
                 [--run-bmtagger] [--run-trf] [--remove-temp-output]
                 [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--log LOG]
                 [--trimmomatic TRIMMOMATIC_PATH] [--max-memory MAX_MEMORY]
                 [--trimmomatic-options TRIMMOMATIC_OPTIONS]
                 [--bowtie2 BOWTIE2_PATH] [--bowtie2-options BOWTIE2_OPTIONS]
                 [--bmtagger BMTAGGER_PATH] [--trf TRF_PATH] [--match MATCH]
                 [--mismatch MISMATCH] [--delta DELTA] [--pm PM] [--pi PI]
                 [--minscore MINSCORE] [--maxperiod MAXPERIOD]
                 [--no-generate-fastq] [--mask] [--html]

KneadData

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         additional output is printed

global options:
  --version             show program's version number and exit
  -i INPUT, --input INPUT
                        input FASTQ file (add a second argument instance to run with paired input files)
  -o OUTPUT_DIR, --output OUTPUT_DIR
                        directory to write output files
  -db REFERENCE_DB, --reference-db REFERENCE_DB
                        location of reference database (additional arguments add databases)
  --output-prefix OUTPUT_PREFIX
                        prefix for all output files
                        [ DEFAULT : $SAMPLE_kneaddata ]
  -t <1>, --threads <1>
                        number of threads
                        [ Default : 1 ]
  -p <1>, --processes <1>
                        number of processes
                        [ Default : 1 ]
  -q {phred33,phred64}, --quality-scores {phred33,phred64}
                        quality scores
                        [ DEFAULT : phred33 ]
  -s {memory,storage}, --strategy {memory,storage}
                        define operating strategy
                        [ DEFAULT : storage ]
  --run-bmtagger        run BMTagger instead of Bowtie2 to identify contaminant reads
  --run-trf             run TRF to remove and/or mask tandem repeats
  --remove-temp-output  remove temp output files
                        [ DEFAULT : temp output files are not removed ]
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        level of log messages
                        [ DEFAULT : DEBUG ]
  --log LOG             log file
                        [ DEFAULT : $OUTPUT_DIR/$SAMPLE_kneaddata.log ]

trimmomatic arguments:
  --trimmomatic TRIMMOMATIC_PATH
                        path to trimmomatic
                        [ DEFAULT : $PATH ]
  --max-memory MAX_MEMORY
                        max amount of memory
                        [ DEFAULT : 500m ]
  --trimmomatic-options TRIMMOMATIC_OPTIONS
                        options for trimmomatic
                        [ DEFAULT : SLIDINGWINDOW:4:20 MINLEN:60 ]

bowtie2 arguments:
  --bowtie2 BOWTIE2_PATH
                        path to bowtie2
                        [ DEFAULT : $PATH ]
  --bowtie2-options BOWTIE2_OPTIONS
                        options for bowtie2
                        [ DEFAULT : --very-sensitive ]

bmtagger arguments:
  --bmtagger BMTAGGER_PATH
                        path to BMTagger
                        [ DEFAULT : $PATH ]

trf arguments:
  --trf TRF_PATH        path to TRF
                        [ DEFAULT : $PATH ]
  --match MATCH         matching weight
                        [ DEFAULT : 2 ]
  --mismatch MISMATCH   mismatching penalty
                        [ DEFAULT : 7 ]
  --delta DELTA         indel penalty
                        [ DEFAULT : 7 ]
  --pm PM               match probability
                        [ DEFAULT : 80 ]
  --pi PI               indel probability
                        [ DEFAULT : 10 ]
  --minscore MINSCORE   minimum alignment score to report
                        [ DEFAULT : 50 ]
  --maxperiod MAXPERIOD
                        maximum period size to report
                        [ DEFAULT : 500 ]
  --no-generate-fastq   don't generate fastq output from trf
  --mask                generate mask file from trf output
  --html                generate html file for trf output
```

----------------------------------

## Usage Notes

### Specifying Additional Arguments

If you want to specify additional arguments for Bowtie2 using the
`--bowtie2-options` flag, you will need to use the equals sign along with quotes. Add additional flags for each option.

For example:

`$ kneaddata --input demo.fastq --output kneaddata_output --reference-db database_folder --bowtie2-options="--very-fast" --bowtie2-options="-p 2"`

A similar approach is used to specify additional arguments for Trimmomatic:

`$ kneaddata --input demo.fastq --output kneaddata_output --reference-db database_folder --trimmomatic-options="LEADING:3" --trimmomatic-options="TRAILING:3"`

*NOTE*: Manually specifying additional arguments will completely override the
defaults. 

