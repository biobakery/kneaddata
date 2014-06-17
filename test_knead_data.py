'''
test_knead_data.py: Python script for running unit tests for knead_data.py. You
might want to run this from its own folder because it creates a lot of temporary
files which may interfere with yours. 

Author: Andy Shi
'''

import knead_data as kd
import unittest
import os
import subprocess
import shlex

class test_checkfile(unittest.TestCase):
    def setUp(self): 
        with open("empty.trimmed.fastq", "w") as f:
            pass
        with open("non-empty.trimmed.1.fastq", "w") as f:
            f.write("Hello world!")
        self.files = ["empty.trimmed.fastq", "non-empty.trimmed.1.fastq",
                "dne.txt"]

        self.bad_file1 = (1,2,3)
        self.bad_file2 = 1

    def test_checkfile(self):
        files_res = map(kd.checkfile, self.files)
        self.assertEqual(files_res, [-1, 1, 0])

        self.assertRaises(TypeError, kd.checkfile, self.bad_file1, 
                "passing tuple to checkfile")
        self.assertRaises(TypeError, kd.checkfile, self.bad_file2,
                "Passing in int to checkfile")
        self.assertRaises(IOError, kd.checkfile, self.files[0], fail_hard=True)
        self.assertRaises(OSError, kd.checkfile, self.files[2], fail_hard=True)

    def tearDown(self):
        os.remove("empty.trimmed.fastq")
        os.remove("non-empty.trimmed.1.fastq")



class test_is_single_end(unittest.TestCase):
    def setUp(self):
        self.good1 = "abc.fastq"
        self.good2 = "efg.fastq"
        self.bad = None

    def test_is_single_end(self):
        bool1, out1 = kd.is_single_end(self.good1, self.good2)
        bool2, out2 = kd.is_single_end(self.good1, self.bad)
        self.assertFalse(bool1, "b_single_end, pass in 2 fastq's")
        self.assertEqual(out1, [self.good1, self.good2])
        self.assertTrue(bool2, "b_single_end, pass in 1 fastq")
        self.assertEqual(out2, [self.good1])
        


class test_checktrim_output(unittest.TestCase):
    def setUp(self):
        self.out = "out"
        self.files =    ["out.trimmed.fastq", 
                        "out.trimmed.1.fastq",
                        "out.trimmed.2.fastq", 
                        "out.trimmed.single.1.fastq",
                        "out.trimmed.single.2.fastq"]
        for fname in self.files:
            with open(fname, "w") as f:
                pass

    def test_badInput(self):
        ''' passing in random things instead of strings, looking to get type
        errors '''
        self.assertRaises(TypeError, kd.checktrim_output, 1, True)
        self.assertRaises(TypeError, kd.checktrim_output, [1,2,3], True)

    def test_file_dne(self):
        ''' >= 1 of the files is not there '''
        b, o, inp = kd.checktrim_output("dne", True)
        self.assertFalse(b)

    def test_se_empty(self):
        ''' single end, file is empty '''
        b1, o1, in1 = kd.checktrim_output(self.out, True)
        self.assertFalse(b1, "single end with empty input")
        self.assertEqual(o1, ["out.trimmed.fastq"], 
            "single end correct output")
        self.assertEqual(in1, [], "input for failed checking should be []")

    def test_se_nonempty(self):
        ''' single end, file is non-empty '''
        with open("out.trimmed.fastq", "w") as f:
            f.write("hello world")
        b, out, inp = kd.checktrim_output(self.out, True)
        self.assertTrue(b, "single end with nonempty input")
        self.assertEqual(out, ["out.trimmed.fastq"], 
            "single end correct output")
        self.assertEqual(inp, [["out.trimmed.fastq"]], "single end correct input")

    def test_pe_allEmpty(self):
        ''' paired end, all files are empty '''
        b, o, inp = kd.checktrim_output(self.out, False)
        self.assertFalse(b, "paired end with empty input should be false")
        self.assertEqual(o, ["out.trimmed.1.fastq", "out.trimmed.2.fastq",
        "out.trimmed.single.1.fastq", "out.trimmed.single.2.fastq"],
        "paired end correct output")
        self.assertEqual(inp, [], "paired end empty; should be empty input")
    
    def test_pe_peEmpty(self):
        ''' paired end, one of the paired files is empty '''
        with open("out.trimmed.1.fastq", "w") as f:
            f.write("hello world")
        b1, o1, inp1 = kd.checktrim_output(self.out, False)
        self.assertTrue(b1, "paired end with non-empty input")
        self.assertTrue(["out.trimmed.1.fastq"] in inp1)

        # switch the roles of out.trimmed.1.fastq and out.trimmed.2.fastq

        os.remove("out.trimmed.1.fastq")
        with open("out.trimmed.1.fastq", "w") as f:
            pass
        with open("out.trimmed.2.fastq", "w") as f:
            f.write("non-empty")

        b2, o2, inp2 = kd.checktrim_output(self.out, False)
        self.assertTrue(b2)
        self.assertTrue(["out.trimmed.2.fastq"] in inp2)

    def test_pe_bothPENonempty(self):
        ''' paired end, both paired files non-empty '''
        with open("out.trimmed.1.fastq", "w") as f:
            f.write("non-empty")
        with open("out.trimmed.2.fastq", "w") as f:
            f.write("non-empty")

        b, o, inp = kd.checktrim_output(self.out, False)
        self.assertTrue(b)
        self.assertTrue(["out.trimmed.1.fastq", "out.trimmed.2.fastq"] in inp)


    def test_pe_seEmpty(self):
        ''' paired end, one of the single end files is empty '''
        with open("out.trimmed.single.2.fastq", "w") as f:
            f.write("non-empty")
        
        b, o, inp = kd.checktrim_output(self.out, False)
        self.assertTrue(b)
        self.assertTrue(["out.trimmed.single.2.fastq"] in inp)

    def tearDown(self):
        for fname in self.files:
            os.remove(fname)


class test_get_num_reads(unittest.TestCase):
    def setUp(self):
        # create files
        self.testfile = "words.txt"
        cmd = "head --lines 4000 /usr/share/dict/words"

        out = subprocess.check_output(shlex.split(cmd))

        with open(self.testfile, "w") as f:
            f.write(out)

        with open("empty.txt", "w") as f:
            pass

    def test_notExist(self):
        ''' Nonexistent file '''
        self.assertEqual(kd.get_num_reads("dne.txt", False), None, 
                "Should return None for nonexistent file")

    def test_empty(self):
        ''' Empty file '''
        self.assertEqual(kd.get_num_reads("empty.txt", False), 0, 
                "Should return 0 for empty file")
    
    def test_nonEmpty_fastq(self):
        ''' Normal, nonempty fastq file'''
        self.assertEqual(kd.get_num_reads(self.testfile, True), 1000, 
                "Should return the number of lines in the file divided by 4")

    def test_nonEmpty_out(self):
        '''Normal, nonempty LIST of hits'''
        self.assertEqual(kd.get_num_reads(self.testfile, False), 4000, 
                "Should return the number of lines in the file")

    def tearDown(self):
        os.remove(self.testfile)
        os.remove("empty.txt")

if __name__ == '__main__':
    unittest.main()
