
# KneadData History #

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

