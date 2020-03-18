import unittest
import tempfile
import os

from kneaddata import utilities
from kneaddata import config

import cfg
import utils


def skipIfExeNotFound(exe):
    if isinstance(exe,str):
        exe=[exe]
    if all([utilities.find_exe_in_path(requires, bypass_permissions_check=True) for requires in exe]):
        return lambda func: func
    return unittest.skip("{} is not installed so test is skipped".format(",".join(exe)))

class TestFunctionalKneadData(unittest.TestCase):
    """
    Test KneadData workflows
    """
    
    @skipIfExeNotFound(config.trimmomatic_jar)
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
        
    @skipIfExeNotFound(config.trimmomatic_jar)
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
        
    @skipIfExeNotFound([config.trimmomatic_jar, config.bowtie2_exe])
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
        
    @skipIfExeNotFound([config.trimmomatic_jar, config.bowtie2_exe])
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
                   "--output",tempdir,"--reference-db",cfg.bowtie2_db_folder,"--no-discordant"]
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

    @skipIfExeNotFound([config.trimmomatic_jar, config.bowtie2_exe])
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
                   "--reference-db",cfg.bowtie2_db_folder,"--no-discordant"]
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
          
        # check there are at least the main expected files in the output folder
        actual_output_files=os.listdir(tempdir)
        self.assertGreater(len(actual_output_files), len(expected_output_files))

        # remove the temp directory
        utils.remove_temp_folder(tempdir)

    @skipIfExeNotFound([config.trimmomatic_jar, config.bowtie2_exe])
    def test_trimmomatic_bowtie2_two_databases_paired_end_serial(self):
        """
        Test running the default flow of trimmomatic on paired end input with two
        bowtie2 database provided (both with the same name)
        Test running in serial alignment mode
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,"--input",cfg.fastq_file,
                   "--output",tempdir,"--reference-db",cfg.bowtie2_db_folder,
                   "--reference-db",cfg.bowtie2_db_folder,"--no-discordant","--serial"]
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
          
        # check there are at least the main expected files in the output folder
        actual_output_files=os.listdir(tempdir)
        self.assertGreater(len(actual_output_files), len(expected_output_files))

        # remove the temp directory
        utils.remove_temp_folder(tempdir)

    @skipIfExeNotFound([config.trimmomatic_jar, config.bowtie2_exe])
    def test_trimmomatic_bowtie2_paired_end_remove_intermedite_temp_output_discordant(self):
        """
        Test running the default flow of trimmomatic on paired end input with one
        bowtie2 database provided
        Test running with remove intermediate temp output files
        Test with discordant alignments
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,"--input",cfg.fastq_pair_file,
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
          
        # check there are the expected number of files in the output folder
        actual_output_files=list(filter(os.path.getsize,[os.path.join(tempdir,file) for file in os.listdir(tempdir)]))
        self.assertEqual(len(actual_output_files), len(expected_non_empty_output_files))

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    @skipIfExeNotFound([config.trimmomatic_jar, config.bowtie2_exe, config.trf_exe])
    def test_trimmomatic_bowtie2_paired_end_remove_intermedite_temp_output_discordant_trf(self):
        """
        Test running the default flow of trimmomatic on paired end input with one
        bowtie2 database provided
        Test running with remove intermediate temp output files
        Test with discordant alignments
        Test with TRF
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,"--input",cfg.fastq_pair_file,
                   "--output",tempdir,"--reference-db",cfg.bowtie2_db_folder,
                   "--reference-db",cfg.bowtie2_db_folder, "--run-trf"]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=basename=utils.file_basename(cfg.fastq_file)
        filtered_file_basename=utils.get_filtered_file_basename(basename,cfg.bowtie2_db_folder,"bowtie2")
        
        expected_non_empty_output_files=[basename+cfg.log_extension,
                               basename+cfg.paired_trim_extensions[0],
                               basename+cfg.paired_trim_extensions[1],
                               basename+cfg.final_extensions_paired[0],
                               basename+cfg.final_extensions_paired[1],
                               basename+cfg.paired_repeats_removed_extensions[0],
                               basename+cfg.paired_repeats_removed_extensions[1]]
        
        # check the output files are as expected
        for expression, message in utils.check_output(expected_non_empty_output_files, tempdir):
            self.assertTrue(expression,message)
          
        # check there are the expected number of files in the output folder
        actual_output_files=list(filter(os.path.getsize,[os.path.join(tempdir,file) for file in os.listdir(tempdir)]))
        self.assertEqual(len(actual_output_files), len(expected_non_empty_output_files))

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    @skipIfExeNotFound([config.trimmomatic_jar, config.bowtie2_exe, config.trf_exe])
    def test_trimmomatic_bowtie2_database_and_trf_single_end(self):
        """
        Test running the default flow of trimmomatic on single end input with
        bowtie2 database provided
        Test with keeping temp files
        Test with TRF
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,
                   "--output",tempdir,"--reference-db",cfg.bowtie2_db_folder,
                   "--store-temp-output", "--run-trf"]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=utils.file_basename(cfg.fastq_file)
        filtered_file_basename=utils.get_filtered_file_basename(basename,cfg.bowtie2_db_folder,"bowtie2")
        
        expected_output_files=[basename+cfg.log_extension,
                               basename+cfg.single_trim_extension,
                               filtered_file_basename+cfg.clean_extension,
                               filtered_file_basename+cfg.contaminated_extension,
                               filtered_file_basename+cfg.sam_extension,
                               basename+cfg.final_extension,
                               basename+cfg.repeats_removed_extension]
        
        # check the output files are as expected
        for expression, message in utils.check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    @skipIfExeNotFound([config.trimmomatic_jar, config.bowtie2_exe, config.trf_exe])
    def test_trimmomatic_bowtie2_database_and_trf_paired_end_remove_intermedite_temp_output(self):
        """
        Test running the default flow of trimmomatic on paired end input with a
        bowtie2 database provided
        Test running with remove intermediate temp output files
        Test running trf
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,"--input",cfg.fastq_file,
                   "--output",tempdir,"--reference-db",cfg.bowtie2_db_folder,"--run-trf","--no-discordant"]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=basename=utils.file_basename(cfg.fastq_file)
        filtered_file_basename=utils.get_filtered_file_basename(basename,cfg.bowtie2_db_folder,"bowtie2",True)
        
        expected_non_empty_output_files=[basename+cfg.log_extension,
                               basename+cfg.paired_trim_extensions[0],
                               basename+cfg.paired_trim_extensions[1],
                               basename+cfg.final_extensions_paired[0],
                               basename+cfg.final_extensions_paired[1],
                               basename+cfg.paired_repeats_removed_extensions[0],
                               basename+cfg.paired_repeats_removed_extensions[1]]
        
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
        
    @skipIfExeNotFound([config.trimmomatic_jar, config.trf_exe])
    def test_trimmomatic_and_trf_no_reference_database_single_end(self):
        """
        Test running the default flow of trimmomatic on single end input as no
        reference database is provided
        Test with also running trf
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,
                   "--output",tempdir,"--run-trf"]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=utils.file_basename(cfg.fastq_file)
        
        expected_output_files=[basename+cfg.log_extension,
                               basename+cfg.single_trim_extension,
                               basename+cfg.repeats_removed_extension]
        
        # check the output files are as expected
        for expression, message in utils.check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    @skipIfExeNotFound([config.trimmomatic_jar, config.trf_exe])
    def test_trimmomatic_and_trf_no_reference_database_paired_end(self):
        """
        Test running the default flow of trimmomatic on paired end input as no
        reference database is provided
        Test with also running trf
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,"--input",cfg.fastq_file,
                   "--output",tempdir,"--run-trf"]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=utils.file_basename(cfg.fastq_file)
        
        expected_output_files=[basename+cfg.log_extension,
                               basename+cfg.paired_trim_extensions[0],
                               basename+cfg.paired_trim_extensions[1],
                               basename+cfg.paired_repeats_removed_extensions[0],
                               basename+cfg.paired_repeats_removed_extensions[1]]
        
        # check the output files are as expected
        for expression, message in utils.check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)

    @skipIfExeNotFound(config.bowtie2_exe)
    def test_bowtie2_only_single_end(self):
        """
        Test on single end input with bowtie2 database provided
        Test with keeping temp files
        Test bypassing trim step
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,
                   "--output",tempdir,"--reference-db",cfg.bowtie2_db_folder,
                   "--store-temp-output","--bypass-trim"]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=utils.file_basename(cfg.fastq_file)
        filtered_file_basename=utils.get_filtered_file_basename(basename,cfg.bowtie2_db_folder,"bowtie2")
        
        expected_output_files=[basename+cfg.log_extension,
                               filtered_file_basename+cfg.clean_extension,
                               filtered_file_basename+cfg.contaminated_extension,
                               filtered_file_basename+cfg.sam_extension,
                               basename+cfg.final_extension]
        
        # check the output files are as expected
        for expression, message in utils.check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    @skipIfExeNotFound(config.bowtie2_exe)
    def test_bowtie2_only_paired_end_remove_intermedite_temp_output(self):
        """
        Test running the default flow of trimmomatic on paired end input with a
        bowtie2 database provided
        Test running with remove intermediate temp output files
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,"--input",cfg.fastq_file,
                   "--output",tempdir,"--reference-db",cfg.bowtie2_db_folder,"--bypass-trim","--no-discordant"]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=basename=utils.file_basename(cfg.fastq_file)
        filtered_file_basename=utils.get_filtered_file_basename(basename,cfg.bowtie2_db_folder,"bowtie2",True)
        
        expected_non_empty_output_files=[basename+cfg.log_extension,
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

    @skipIfExeNotFound(config.trf_exe)
    def test_trf_only_single_end(self):
        """
        Test running trf only on single end input
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,
                   "--output",tempdir,"--run-trf","--bypass-trim"]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=utils.file_basename(cfg.fastq_file)
        
        expected_output_files=[basename+cfg.log_extension,
                               basename+cfg.repeats_removed_extension]
        
        # check the output files are as expected
        for expression, message in utils.check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    @skipIfExeNotFound(config.trf_exe)
    def test_trf_only_paired_end(self):
        """
        Test running only trf on paired end input
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,"--input",cfg.fastq_file,
                   "--output",tempdir,"--run-trf","--bypass-trim"]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=utils.file_basename(cfg.fastq_file)
        
        expected_output_files=[basename+cfg.log_extension,
                               basename+cfg.paired_repeats_removed_extensions[0],
                               basename+cfg.paired_repeats_removed_extensions[1]]
        
        # check the output files are as expected
        for expression, message in utils.check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    @skipIfExeNotFound([config.trimmomatic_jar, config.bowtie2_exe, config.trf_exe])
    def test_trimmomatic_bowtie2_database_and_trf_single_end_gzipped_input(self):
        """
        Test running the default flow of trimmomatic on single end input with
        bowtie2 database provided
        Test with keeping temp files
        Test with TRF
        Test with gzipped input fastq file
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file_gzipped,
                   "--output",tempdir,"--reference-db",cfg.bowtie2_db_folder,
                   "--store-temp-output", "--run-trf"]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        fastq_file_basename=utils.file_basename(cfg.fastq_file_gzipped)
        basename=utils.file_basename(cfg.fastq_file)
        filtered_file_basename=utils.get_filtered_file_basename(basename,cfg.bowtie2_db_folder,"bowtie2")
        
        expected_output_files=[basename+cfg.log_extension,
                               basename+cfg.single_trim_extension,
                               filtered_file_basename+cfg.clean_extension,
                               filtered_file_basename+cfg.contaminated_extension,
                               filtered_file_basename+cfg.sam_extension,
                               basename+cfg.final_extension,
                               basename+cfg.repeats_removed_extension]
	
        # check the output files are as expected
        for expression, message in utils.check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    @skipIfExeNotFound([config.trimmomatic_jar, config.bowtie2_exe, config.trf_exe, config.samtools_exe])
    def test_trimmomatic_bowtie2_database_and_trf_single_end_bam_input(self):
        """
        Test running the default flow of trimmomatic on single end input with
        bowtie2 database provided
        Test with keeping temp files
        Test with TRF
        Test with bam input fastq file
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.file_bam,
                   "--output",tempdir,"--reference-db",cfg.bowtie2_db_folder,
                   "--store-temp-output", "--run-trf"]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=utils.file_basename(cfg.file_bam)
        filtered_file_basename=utils.get_filtered_file_basename(basename,cfg.bowtie2_db_folder,"bowtie2")
        
        expected_output_files=[basename+"_decompressed.fastq",
                               basename+"_decompressed.sam",
                               basename+cfg.log_extension,
                               basename+cfg.single_trim_extension,
                               filtered_file_basename+cfg.clean_extension,
                               filtered_file_basename+cfg.contaminated_extension,
                               filtered_file_basename+cfg.sam_extension,
                               basename+cfg.final_extension,
                               basename+cfg.repeats_removed_extension]
        
        # check the output files are as expected
        for expression, message in utils.check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    @skipIfExeNotFound([config.trimmomatic_jar, config.bowtie2_exe, config.trf_exe])
    def test_trimmomatic_bowtie2_database_and_trf_single_end_sam_input(self):
        """
        Test running the default flow of trimmomatic on single end input with
        bowtie2 database provided
        Test with keeping temp files
        Test with TRF
        Test with sam input fastq file
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.file_sam,
                   "--output",tempdir,"--reference-db",cfg.bowtie2_db_folder,
                   "--store-temp-output", "--run-trf"]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=utils.file_basename(cfg.file_sam)
        filtered_file_basename=utils.get_filtered_file_basename(basename,cfg.bowtie2_db_folder,"bowtie2")
        
        expected_output_files=[basename+"_decompressed.fastq",
                               basename+cfg.log_extension,
                               basename+cfg.single_trim_extension,
                               filtered_file_basename+cfg.clean_extension,
                               filtered_file_basename+cfg.contaminated_extension,
                               filtered_file_basename+cfg.sam_extension,
                               basename+cfg.final_extension,
                               basename+cfg.repeats_removed_extension]
        
        # check the output files are as expected
        for expression, message in utils.check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    @skipIfExeNotFound([config.trimmomatic_jar, config.fastqc_exe])
    def test_trimmomatic_fastqc_start_no_reference_database_single_end(self):
        """
        Test running the default flow of trimmomatic on single end input as no
        reference database is provided
        Test running fastqc at the beginning of the workflow
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,
                   "--output",tempdir,"--run-fastqc-start"]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=utils.file_basename(cfg.fastq_file)
        
        expected_output_files=[os.path.join("fastqc",basename+cfg.fastqc_extensions[0]),
                               os.path.join("fastqc",basename+cfg.fastqc_extensions[1]),
                               basename+cfg.log_extension,
                               basename+cfg.single_trim_extension]
        
        # check the output files are as expected
        for expression, message in utils.check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    @skipIfExeNotFound([config.trimmomatic_jar, config.fastqc_exe])
    def test_trimmomatic_fastqc_start_no_reference_database_paired_end(self):
        """
        Test running the default flow of trimmomatic on paired end input as no
        reference database is provided
        Test running fastqc at the beginning of the workflow
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,"--input",cfg.fastq_file,
                   "--output",tempdir,"--run-fastqc-start"]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=utils.file_basename(cfg.fastq_file)
        
        expected_output_files=[os.path.join("fastqc",basename+cfg.fastqc_extensions[0]),
                               os.path.join("fastqc",basename+cfg.fastqc_extensions[1]),
                               basename+cfg.log_extension,
                               basename+cfg.paired_trim_extensions[0],
                               basename+cfg.paired_trim_extensions[1]]
        
        # check the output files are as expected
        for expression, message in utils.check_output(expected_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    @skipIfExeNotFound([config.trimmomatic_jar, config.bowtie2_exe, config.fastqc_exe])
    def test_trimmomatic_bowtie2_database_fastqc_end_single_end(self):
        """
        Test running the default flow of trimmomatic on single end input with
        bowtie2 database provided
        Test with keeping temp files
        Test running fastqc at the end of the workflow
        """
        
        # create a temp directory for output
        tempdir = tempfile.mkdtemp(suffix="test_kneaddata_")
        
        # run kneaddata test
        command = ["kneaddata","--input",cfg.fastq_file,
                   "--output",tempdir,"--reference-db",cfg.bowtie2_db_folder,
                   "--store-temp-output", "--run-fastqc-end"]
        utils.run_kneaddata(command)
        
        # get the basename of the input file
        basename=utils.file_basename(cfg.fastq_file)
        final_basename=utils.file_basename(basename+cfg.final_extension)
        filtered_file_basename=utils.get_filtered_file_basename(basename,cfg.bowtie2_db_folder,"bowtie2")
        
        expected_output_files=[os.path.join("fastqc",final_basename+cfg.fastqc_extensions[0]),
                               os.path.join("fastqc",final_basename+cfg.fastqc_extensions[1]),
                               basename+cfg.log_extension,
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
	
