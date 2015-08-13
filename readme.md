# KneadData #

KneadData is a tool designed to perform quality control on metagenomic sequencing data, especially data from microbiome experiments. In these experiments, samples are typically taken from a host in hopes of learning something about the microbial community on the host. However, metagenomic sequencing data from such experiments will often contain a high ratio of host to bacterial reads. This tool aims to perform principled in silico separation of bacterial reads from these "contaminant" reads, be they from the host, from bacterial 16S sequences, or other user-defined sources.

**If you use the KneadData software, please cite our manuscript: TBD**

For additional information, please see the [KneadData User Manual](http://huttenhower.sph.harvard.edu/kneaddata/manual).

## Contents ##

* [Requirements](#markdown-header-requirements)
* [Installation](#markdown-header-installation)
* [How to run](#markdown-header-how-to-run)
    * [Basic usage](#markdown-header-basic-usage)
    * [Demo runs](#markdown-header-demo-runs)


## Requirements ##

1.  [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
2.  [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (version >= 2.1)
3.  [Python](http://www.python.org/) (version >= 2.7)
4.  Operating system (Linux or Mac)

Optionally, [BMTagger](ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/) can be used instead of Bowtie2.

The executables for the required software packages should be installed in your $PATH. Alternatively, you can provide the location of the Bowtie2 install with the following KneadData option “--bowtie2-path”. 

## Installation ##

1. Download and unpack the [KneadData software](https://bitbucket.org/biobakery/kneaddata/get/tip.tar.gz)
2. From the KneadData directory, install KneadData
 
    `` $ python setup.py install ``
    

3. Download the reference database to $INSTALL_LOCATION

    `` $ kneaddata_database --download human bowtie2 $INSTALL_LOCATION ``


## How to Run ##

### Basic usage ###

`` $ kneaddata --infile1 $INPUT --reference-db $DATABASE --trim-path $TRIM_PATH ``

```
$INPUT = a single end fastq file
$DATABASE = the index of the KneadData database
$TRIM_PATH = the full path to the Trimmomatic Java archive
```

For paired end reads, add the option “--infile2 $INPUT2” to provide the second set of reads. Also please note that more than one reference database can be provided.

Three types of output files will be created:

1. `` $INPUT_$DATABASE_clean.fastq ``
2. `` $INPUT_$DATABASE_contam.fastq ``
3. `` $INPUT_output.fastq ``

Files of type #1 and #2 will be created for each reference database provided.


### Demo run ###

The examples folder contains a demo input file. This file is a single read, fastq format.

`` $ kneaddata --infile1 examples/demo.fastq --reference-db $DATABASE --trim-path $TRIM_PATH ``


