
# KneadData History #

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

