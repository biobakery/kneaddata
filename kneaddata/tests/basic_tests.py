import unittest
import tempfile
import os
import logging

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
        
        self.assertEqual(sorted(sequences), sorted(cfg.merges_files_1_sequences))

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

