# KneadData #

----

 * Download the KneadData software ( [kneaddata.tar.gz](https://bitbucket.org/biobakery/kneaddata/downloads/kneaddata-v0.5.0.tar.gz) ) then follow the [steps to install and run](#markdown-header-getting-started-with-kneaddata).

 * If you use the KneadData software, please cite our manuscript: TBD

 * If you use this tool, sign up for the [KneadData Users Google Group](https://groups.google.com/d/forum/kneaddata-users) and pass along any issues or feedback.

----

KneadData is a tool designed to perform quality control on metagenomic sequencing data, especially data from microbiome experiments. In these experiments, samples are typically taken from a host in hopes of learning something about the microbial community on the host. However, metagenomic sequencing data from such experiments will often contain a high ratio of host to bacterial reads. This tool aims to perform principled in silico separation of bacterial reads from these "contaminant" reads, be they from the host, from bacterial 16S sequences, or other user-defined sources.


## Getting Started with KneadData ##

### Requirements ###

1.  [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (version == 0.33) (automatically installed)
2.  [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (version >= 2.1) (automatically installed)
3.  [Python](http://www.python.org/) (version >= 2.7)
4.  [Java Runtime Environment](http://www.oracle.com/technetwork/java/javase/downloads/jre7-downloads-1880261.html)
5.  Memory (>= 4 Gb if using Bowtie2, >= 8 Gb if using BMTagger)
6.  Operating system (Linux or Mac)

Optionally, [BMTagger](ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/) can be used instead of Bowtie2.

The executables for the required software packages should be installed in your $PATH. Alternatively, you can provide the location of the Bowtie2 install ($BOWTIE2_DIR) with the following KneadData option “--bowtie2 $BOWTIE2_DIR”. 

### Installation ###

Before installing KneadData, please install the Java Runtime Environment (JRE). First [download](http://www.oracle.com/technetwork/java/javase/downloads/jre7-downloads-1880261.html) the JRE for your platform. Then follow the instructions for your platform: [Linux 64-bit](http://docs.oracle.com/javase/8/docs/technotes/guides/install/linux_jre.html#CFHIEGAA) or [Mac OS](http://docs.oracle.com/javase/8/docs/technotes/guides/install/mac_jre.html#jre_8u40_osx). At the end of the installation, add the location of the java executable to your $PATH.

1. Download and unpack the KneadData software
    * Download the software: [kneaddata.tar.gz](https://bitbucket.org/biobakery/kneaddata/downloads/kneaddata-v0.5.0.tar.gz)
    * `` $ tar zxvf kneaddata.tar.gz ``
    * `` $ cd kneaddata ``
2. From the KneadData directory, install KneadData
    * `` $ python setup.py install ``
    * This command will automatically install Trimmomatic and Bowtie2. To bypass the install of dependencies, add the option "--bypass-dependencies-install".
    * If you do not have write permissions to '/usr/lib/', then add the option "--user" to the install command. This will install the python package into subdirectories of '~/.local'. Please note when using the "--user" install option on some platforms, you might need to add '~/.local/bin/' to your $PATH as it might not be included by default. You will know if it needs to be added if you see the following message ``kneaddata: command not found`` when trying to run KneadData after installing with the "--user" option.
3. Download the human reference database to $DIR
    * `` $ kneaddata_database --download human bowtie2 $DIR ``
    * When running this command, $DIR should be replaced with the full path to the directory you have selected to store the database.


### How to Run ###

#### Basic usage ####

`` $ kneaddata --input $INPUT --reference-db $DATABASE --output $OUTPUT_DIR ``

```
$INPUT = a single end fastq file
$DATABASE = the index of the KneadData database
$OUTPUT_DIR = the output directory
```

For paired end reads, add a second input argument “--input $INPUT2” (with $INPUT2 replaced with the second input file). Also please note that more than one reference database can be provided in the same manner by using multiple database options (for example, "--reference-db $DATABASE1 --reference-db $DATABASE2").

Three types of output files will be created (where $INPUTNAME is the basename of $INPUT):

1. The final file of filtered sequences after trimming
    * `` $OUTPUT_DIR/$INPUTNAME_kneaddata.fastq ``

2. The contaminant sequences from testing against a database (with this database name replacing $DATABASE)
    * `` $OUTPUT_DIR/$INPUTNAME_kneaddata_$DATABASE_contam.fastq ``

3. The log file from the run
    * `` $OUTPUT_DIR/$INPUTNAME_kneaddata.log ``

4. The fastq file of trimmed sequences
    * `` $OUTPUT_DIR/$INPUTNAME_kneaddata.trimmed.fastq ``
    * Trimmomatic is run with the following arguments by default "SLIDINGWINDOW:4:20 MINLEN:60". To change the Trimmomatic arguments, use the option "--trimmomatic-options".


If there is more than one reference database, there will be a fourth output file type. Files of this type will be named `` $OUTPUT_DIR/$INPUTNAME_kneaddata_$DATABASE_clean.fastq `` and will contain the filtered sequences after testing against a specific database (with this database name replacing $DATABASE in the file name). The file `` $OUTPUT_DIR/$INPUTNAME_kneaddata.fastq `` is the set of all sequences contained in these filtered files.

If running with two input files, each type of fastq output file will be created for each one of the pairs of the input files. If running with the TRF step, an additional set of files with repeats removed will be written.

#### Demo run ####

The examples folder contains a demo input file. This file is a single read, fastq format.

`` $ kneaddata --input examples/demo.fastq --reference-db examples/demo_db --output kneaddata_demo_output ``

This will create three output files:

1. `` kneaddata_demo_output/demo_kneaddata.fastq ``
2. `` kneaddata_demo_output/demo_kneaddata_demo_db_contam.fastq ``
3. `` kneaddata_demo_output/demo_kneaddata.log ``
3. `` kneaddata_demo_output/demo_kneaddata.trimmed.log ``

## User Manual ##
For the full user manual and more advanced usage, see the [User Manual](doc/UserManual.md).
