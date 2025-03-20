
# KneadData History #
## v0.12.3 TBD
* Fixed auto-close issue template
* Fixed hg39 T2T Human genome download issue
* Added cat database

## v0.12.2 02-24-2025
* Fixed trimmomatic conda executable 
* Added the hg39 T2T Human genome reference database
* Fixed TRF parallel run bug 
* TRF broken download links updated
  
## v0.12.1 01-06-2025
* Added dog reference genome database
* PyPi + conda release updated

## v0.12.0 09-13-2022
* Allow for tracking of paired-end reads from bam input files
* Add an option to use a scratch directory

## v0.11.0 06-03-2022
* Change the input format to --input1 --input2 --unpaired
* Trim the spaces in the identifiers sequences of the input files
* Add the identifier sequence number if not present in the input files
* Allow for files with less then N lines (in checks for sequence identifiers)

## v0.10.0 02-08-2021

* Remove the file sort from the function that counts the reads so the file name and database name are always in sync. Resolves issue with workflows where read counts are out of sync for database names.
* Change trimmomatic min read length to 60 nt.
* Add back in legacy option "--run-trf" for backwards compatibility with workflows.

## v0.9.0 11-24-2020

* Changed bowtie2 settings from --very-sensitive to --very-sensitive-local
* Checking 1st and last 100 seq and removing spaces from seq identifiers if present
* Add min len to start of trimmomatic run to try to resolve trimmomatic hanging during adapter trimming with short reads

## v0.8.0 07-31-2020
* Add another bowtie2 option for filtering of paired end reads. This "strict" mode, the new default, removes both reads if either pair aligns. Modify existing pair alignment options into a single option.

## v0.7.10 07-29-2020

* Remove trim-repetitive sequences option (prior set as default for filtering to remove adapters and overrepresented sequences).
* Add more adapters plus option to select adapter type.

## v0.7.9 07-14-2020

* Modified check for ordered pairs to only check the first 100 sequence ids to reduce memory usage.
* Move resolve sublists to allow for multiple lists in final output files (resolves error with running with multiple databases).

## v0.7.8 07-08-2020

* Remove temp files from decompression and reformatting after they are no longer needed to save space.
* Update SILVA_128_LSUParc_SSUParc_ribosomal_RNA to Version 0.2 
* Change to the read mode for the decompress function for bz2 files for python 3 compatibility.

## v0.7.7-alpha 05-27-2020
* Changes in the workflow order of TRF and bowtie2
* Add trim-repetitive option flag to remove adapters and overrepresented sequences.

## v0.7.6 05-05-2020

* Add trf install to setup.

## v0.7.5 03-30-2020

* Add option for user to provide the database location (local or remote) with version check.

## v0.7.4 12-04-2019

* Check for order of pairs if trimmomatic is run and reorder if needed (trimmomatic requires ordered pairs).
* Allow for bz2 input files.

## v0.7.3 10-16-2019

* Fix merge final file with run trf and concat pairs as prior just included orphans (used with workflow).
* Update to download new human database (hg37).

## v0.7.2 12-04-2018

* Remove extra java option in trimmomatic call allowing for compatibility with the latest java versions.

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

