
# KneadData History #

## v0.7.1 10-15-2018

* Relax criteria for detecting new Illumina sequence identifier to include ids that are truncated to at least 7 tokens of colon-delimited information.

## v0.7.0 11-14-2017

* Adapter trimming has been added to the Trimmomatic step by default.
* The paired default mode of concordant plus discordant alignment has been removed and replaced with discordant alignment. An evaluation of these two modes found they are relatively the same (reads filtered differed by at most 1%) and the concordant plus discordant mode is not supported in newer versions of bowtie2. Concordant only alignments is still an option available for paired end reads.
* The default min length percent has been reduced from 70% to 50% for the Trimmomatic step.
* Newer functions (like reformat identifiers and discordant pairs) have been modified to fix error messages when running with python3.

## v0.6.1 06-14-2017

* A new option "--remove-intermediate-output" allows for the removal of output files that were generated during the workflow. This option is helpful in reducing the total size of products generated which can be large for many samples.

## v0.6.0 05-24-2017

* A new option "--serial" allows for alignments to be processed in a sequential manner. The clean reads from aligning to the first reference database will be used as input to the second database. The input/outputs will chain until all databases have been processed. Databases will be run in the same order as they are provided on the command line. This differs from the default mode which will run the same input files through each database, in parallel if the "--processes N" option is set to more than one.
* The log messages were updated so single end runs with multiple reference databases will have more descriptive column names, each filtered column will include the database name, when read counts are compiled from the logs with the utility script.

## v0.5.5 05-01-2017 ##

* By default alignments for paired input files now include discordant alignments. To restrict these alignments and only allow concordant alignments (which was the default in prior versions) add the option "--no-discordant".
* A new option "--cat-pairs" allows for alignments of paired end reads as single end (ignoring pair identifiers). This differs from the default alignment mode as the default mode will first align as pairs and then as single end if no concordant alignments are identified.
* A new option "--reorder" which can be used with discordant alignments will reorder the clean fastq files to the same order as the original input files.
* If new illumina sequence identifiers are found, they are replaced with identifiers similar to the original format to maintain pair identification with alignments.
* A mouse bowtie2 database is now available for download.

## v0.5.4 12-15-2016 ##

### New Features ###

* A human metatranscriptome bowtie2 database is now available for download.
* A ribosomal RNA bowtie2 database is now available for download.
* A new option "--cat-final-output" will create a single output file from all final output files for paired-end inputs.
* A new script kneaddata_read_count_table will create a table of read counts (raw, trimmed, decontaminated, final) for all logs provided.

## v0.5.3 12-09-2016 ##

### Bug Fixes ###

* If the output prefix selected is the same as the basename of the input file and the input file is gzipped, the final output file for single-end reads would be deleted as the name was the same as the temp gunzipped input file. A temp file is now used to store the gunzipped input file to prevent overlap.

## v0.5.2 10-05-2016 ##

### New Features ###

* KneadData is now python3 compatible.

## v0.5.1 02-02-2016 ##

### New Features ###

* New options were added to run fastqc at the beginning and/or end of the workflow. Please note these options require fastqc be installed.
* Input files can now be gzipped. They can also be of the SAM or BAM format. Please note with BAM input files SAMTools must be installed.
* The default Trimmomatic options now have the minimum length set to 70% of the length of the reads in the input file.

## v0.5.0 01-22-2016 ##

### New Features ###

* Bypass and run options have been added so a workflow can include any number of tools. There are nine possible combinations to run Trimmomatic, Bowtie2/BMTagger, and TRF.
* Functional and unit tests were added. To run use $ kneaddata_test. Please note this requires external dependencies be installed like Bowtie2.

