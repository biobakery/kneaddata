import unittest
import subprocess
import tempfile
import shutil
import os

import cfg

def run_kneaddata(command):
    """ Run the kneaddata command """
    
    try:
        stdout=subprocess.check_output(command)
    except (EnvironmentError, subprocess.CalledProcessError):
        print("Warning: Unable to execute kneaddata in test.\n"+" ".join(command))    
        
def remove_temp_folder(tempdir):
    """ Remove the temp folder """
    
    try:
        shutil.rmtree(tempdir)
    except EnvironmentError:
        print("Warning: Unable to remove temp directory: " + tempdir)

def check_output(output_files_expected,output_folder):
    """ Check the output folder has the expected file and they are all non-zero """
    
    for file in output_files_expected:
        expected_file = os.path.join(output_folder,file)
        # check the file exists
        yield (os.path.isfile(os.path.join(expected_file)), "File does not exist: " + file)
        
        # check the file is not empty
        yield (os.stat(expected_file).st_size > 0, "File is empty: " + file)
    
def file_basename(file):
    """ Get the basename for the file without the extension """
    
    return os.path.splitext(os.path.basename(file))[0]

def database_basename(folder):
    """ Get the basename for the database in the folder """
    
    # get a single database file from the folder
    database_file=os.listdir(folder)[0]
    
    # some files have multiple extensions, so use split instead of splitext
    return database_file.split(".")[0]

def get_filtered_file_basename(basename,db_folder,software):
    """ Get the basename of the files created by bowtie2 and bmtagger """
    
    database=database_basename(db_folder)
        
    return basename+cfg.name_delimiter+database+"_"+software
    

class TestFunctionalKneadData(unittest.TestCase):
    """
    Test KneadData workflows
    """
    
    def test_trimmomatic_only_no_reference_database_single_end(self):
        """
        Test running the default flow of trimmomatic on single end input as no
        reference database is provided
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,
                   "--output",tempdir]
        run_kneaddata(command)
        
        # get the basename of the input file
        basename=file_basename(cfg.fastq_file)
        
        expected_output_files=[basename+cfg.log_extension,
                               basename+cfg.single_trim_extension]
        
        # check the output files are as expected
        for expression, message in check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        remove_temp_folder(tempdir)
        
    def test_trimmomatic_only_no_reference_database_paired_end(self):
        """
        Test running the default flow of trimmomatic on paired end input as no
        reference database is provided
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,"--input",cfg.fastq_file,
                   "--output",tempdir]
        run_kneaddata(command)
        
        # get the basename of the input file
        basename=file_basename(cfg.fastq_file)
        
        expected_output_files=[basename+cfg.log_extension,
                               basename+cfg.paired_trim_extensions[0],
                               basename+cfg.paired_trim_extensions[1]]
        
        # check the output files are as expected
        for expression, message in check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        remove_temp_folder(tempdir)
        
    def test_trimmomatic_bowtie2_database_single_end(self):
        """
        Test running the default flow of trimmomatic on single end input with
        bowtie2 database provided
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,
                   "--output",tempdir,"--reference-db",cfg.bowtie2_db_folder]
        run_kneaddata(command)
        
        # get the basename of the input file
        basename=file_basename(cfg.fastq_file)
        filtered_file_basename=get_filtered_file_basename(basename,cfg.bowtie2_db_folder,"bowtie2")
        
        expected_output_files=[basename+cfg.log_extension,
                               basename+cfg.single_trim_extension,
                               filtered_file_basename+cfg.clean_extension,
                               filtered_file_basename+cfg.contaminated_extension,
                               filtered_file_basename+cfg.sam_extension,
                               basename+cfg.final_extension]
        
        # check the output files are as expected
        for expression, message in check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        remove_temp_folder(tempdir)
        
    def test_trimmomatic_bowtie2_database_paired_end_remove_intermedite_temp_output(self):
        """
        Test running the default flow of trimmomatic on paired end input with a
        bowtie2 database provided
        Test running with remove intermediate temp output files
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,"--input",cfg.fastq_file,
                   "--output",tempdir,"--reference-db",cfg.bowtie2_db_folder,
                   "--remove-temp-output"]
        run_kneaddata(command)
        
        # get the basename of the input file
        basename=basename=file_basename(cfg.fastq_file)
        filtered_file_basename=get_filtered_file_basename(basename,cfg.bowtie2_db_folder,"bowtie2")
        
        expected_output_files=[basename+cfg.log_extension,
                               basename+cfg.paired_trim_extensions[0],
                               basename+cfg.paired_trim_extensions[1],
                               filtered_file_basename+cfg.clean_extension,
                               filtered_file_basename+cfg.contaminated_extension,
                               filtered_file_basename+cfg.sam_extension,
                               basename+cfg.final_extensions_paired[0],
                               basename+cfg.final_extensions_paired[1]]
        
        # check the output files are as expected
        for expression, message in check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        remove_temp_folder(tempdir)

