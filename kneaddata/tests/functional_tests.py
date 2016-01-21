import unittest
import tempfile
import os

import cfg
import utils


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
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=utils.file_basename(cfg.fastq_file)
        
        expected_output_files=[basename+cfg.log_extension,
                               basename+cfg.single_trim_extension]
        
        # check the output files are as expected
        for expression, message in utils.check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
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
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=utils.file_basename(cfg.fastq_file)
        
        expected_output_files=[basename+cfg.log_extension,
                               basename+cfg.paired_trim_extensions[0],
                               basename+cfg.paired_trim_extensions[1]]
        
        # check the output files are as expected
        for expression, message in utils.check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    def test_trimmomatic_bowtie2_database_single_end(self):
        """
        Test running the default flow of trimmomatic on single end input with
        bowtie2 database provided
        Test with keeping temp files
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,
                   "--output",tempdir,"--reference-db",cfg.bowtie2_db_folder,
                   "--store-temp-output"]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=utils.file_basename(cfg.fastq_file)
        filtered_file_basename=utils.get_filtered_file_basename(basename,cfg.bowtie2_db_folder,"bowtie2")
        
        expected_output_files=[basename+cfg.log_extension,
                               basename+cfg.single_trim_extension,
                               filtered_file_basename+cfg.clean_extension,
                               filtered_file_basename+cfg.contaminated_extension,
                               filtered_file_basename+cfg.sam_extension,
                               basename+cfg.final_extension]
        
        # check the output files are as expected
        for expression, message in utils.check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
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
                   "--output",tempdir,"--reference-db",cfg.bowtie2_db_folder]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=basename=utils.file_basename(cfg.fastq_file)
        filtered_file_basename=utils.get_filtered_file_basename(basename,cfg.bowtie2_db_folder,"bowtie2",True)
        
        expected_non_empty_output_files=[basename+cfg.log_extension,
                               basename+cfg.paired_trim_extensions[0],
                               basename+cfg.paired_trim_extensions[1],
                               basename+cfg.final_extensions_paired[0],
                               basename+cfg.final_extensions_paired[1]]
        
        # check the output files are as expected
        for expression, message in utils.check_output(expected_non_empty_output_files, tempdir):
            self.assertTrue(expression,message)

        # add the expected output files which can be empty
        expected_output_files=expected_non_empty_output_files
        expected_output_files+=[filtered_file_basename+cfg.paired_contaminated_extension[0],
                               filtered_file_basename+cfg.paired_contaminated_extension[1]]
            
        # check there are only three files in the output folder
        actual_output_files=os.listdir(tempdir)
        self.assertEqual(len(actual_output_files), len(expected_output_files))

        # remove the temp directory
        utils.remove_temp_folder(tempdir)

    def test_trimmomatic_bowtie2_two_databases_paired_end_remove_intermedite_temp_output(self):
        """
        Test running the default flow of trimmomatic on paired end input with two
        bowtie2 database provideded (both with the same name)
        Test running with remove intermediate temp output files
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,"--input",cfg.fastq_file,
                   "--output",tempdir,"--reference-db",cfg.bowtie2_db_folder,
                   "--reference-db",cfg.bowtie2_db_folder]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=basename=utils.file_basename(cfg.fastq_file)
        filtered_file_basename=utils.get_filtered_file_basename(basename,cfg.bowtie2_db_folder,"bowtie2")
        
        expected_non_empty_output_files=[basename+cfg.log_extension,
                               basename+cfg.paired_trim_extensions[0],
                               basename+cfg.paired_trim_extensions[1],
                               basename+cfg.final_extensions_paired[0],
                               basename+cfg.final_extensions_paired[1]]
        
        # check the output files are as expected
        for expression, message in utils.check_output(expected_non_empty_output_files, tempdir):
            self.assertTrue(expression,message)
            
        # add the expected output files which can be empty
        expected_output_files=expected_non_empty_output_files
        expected_output_files+=[filtered_file_basename+cfg.paired_contaminated_extension[0],
                               filtered_file_basename+cfg.paired_contaminated_extension[1]]
          
        # check there are only three files in the output folder
        actual_output_files=os.listdir(tempdir)
        self.assertEqual(len(actual_output_files), len(expected_output_files))

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
