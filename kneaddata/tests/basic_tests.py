import unittest
import tempfile
import os
import logging
import filecmp

import cfg
import utils

from kneaddata import run
from kneaddata import utilities

class TestHumann2Functions(unittest.TestCase):
    """
    Test the functions found in kneaddata
    """
    
    def setUp(self):
        # set up nullhandler for logger
        logging.getLogger('kneaddata').addHandler(logging.NullHandler())
        
    def test_read_in_file_n_lines(self):
        """
        Test the function that reads in a file n lines at a time
        """
        
        # Get the sequences from the file, removing end-lines and the fastq "@" character
        sequences=[lines[0].rstrip()[1:] for lines in utilities.read_file_n_lines(cfg.merge_files[0], 4)]
        
        self.assertEqual(sorted(sequences), sorted(cfg.merge_files_1_sequences))

    def test_intersect_fastq(self):
        """
        Test the intersect_fastq function
        """
        
        file_handle, temp_output_file=tempfile.mkstemp(prefix="kneaddata_test")
        
        run.intersect_fastq(cfg.merge_files, temp_output_file)
        
        # read in the sequences from the fastq output file
        # Get the sequences from the file, removing end-lines and the fastq "@" character
        sequences=[lines[0].rstrip()[1:] for lines in utilities.read_file_n_lines(temp_output_file, 4)]
        
        # remove the temp output file
        utils.remove_temp_file(temp_output_file)
        
        self.assertEqual(sorted(sequences), sorted(cfg.merge_files_sequences_intersect))
        
    def test_count_reads_in_fastq_file(self):
        """
        Test the count reads in fastq file function 
        """
        
        read_count=utilities.count_reads_in_fastq_file(cfg.merge_files[0],False)

        self.assertEqual(read_count, len(cfg.merge_files_1_sequences))
        

    def test_is_file_fastq(self):
        """
        Test the is file fastq function and also the get file format function
        """
        
        self.assertTrue(utilities.is_file_fastq(cfg.merge_files[0]))
        
    def test_find_database_index_folder(self):
        """
        Test the find database index function with a folder as input
        """
        
        db_index=utilities.find_database_index(cfg.bowtie2_db_folder, "bowtie2")
        
        self.assertEqual(db_index, cfg.bowtie2_db_index)
        
    def test_find_database_index_file(self):
        """
        Test the find database index function with a file as input
        """
        
        db_index=utilities.find_database_index(cfg.bowtie2_db_file, "bowtie2")
        
        self.assertEqual(db_index, cfg.bowtie2_db_index)
        
    def test_find_database_index_index(self):
        """
        Test the find database index function with an index as input
        """
        
        db_index=utilities.find_database_index(cfg.bowtie2_db_index, "bowtie2")
        
        self.assertEqual(db_index, cfg.bowtie2_db_index)
        
    def test_sam_to_fastq(self):
        """
        Test the sam to fastq function
        Test sam file contains one read with two mappings (to test it is only
        written once to the fastq output file)
        """
        
        file_handle, temp_output_file=tempfile.mkstemp(prefix="kneaddata_test")
        
        utilities.sam_to_fastq(cfg.file_sam, temp_output_file)
        
        self.assertTrue(filecmp.cmp(temp_output_file, cfg.fastq_file_matches_sam_and_bam,
                                     shallow=False))
        
        utils.remove_temp_file(temp_output_file)
        

