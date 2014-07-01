'''
test_knead_data.py: Python script for running unit tests for knead_data.py. You
might want to run this from its own folder because it creates a lot of temporary
files which may interfere with yours. 

Author: Andy Shi
'''

# Note to self: all tests must begin with test

import knead_data as kd
import unittest
import os
import subprocess
import shlex

class Checkfile(unittest.TestCase):
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



class IsSingleEnd(unittest.TestCase):
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
        


class ChecktrimOutput(unittest.TestCase):
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


class GetNumReads(unittest.TestCase):
    def setUp(self):
        # create files
        self.testfile1 = "words.fastq"
        cmd = "head --lines 4000 /usr/share/dict/words"
        self.testfile2 = "words.out"

        out = subprocess.check_output(shlex.split(cmd))

        with open(self.testfile1, "w") as f:
            f.write(out)

        with open(self.testfile2, "w") as f:
            f.write(out)

        with open("empty.txt", "w") as f:
            pass

    def test_notExist(self):
        ''' Nonexistent file '''
        self.assertEqual(kd.get_num_reads("dne.txt"), None, 
                "Should return None for nonexistent file")

    def test_empty(self):
        ''' Empty file '''
        self.assertEqual(kd.get_num_reads("empty.txt"), 0, 
                "Should return 0 for empty file")
    
    def test_nonEmpty_fastq(self):
        ''' Normal, nonempty fastq file'''
        self.assertEqual(kd.get_num_reads(self.testfile1), 1000, 
                "Should return the number of lines in the file divided by 4")

    def test_nonEmpty_out(self):
        '''Normal, nonempty LIST of hits'''
        self.assertEqual(kd.get_num_reads(self.testfile2), 4000, 
                "Should return the number of lines in the file")

    def tearDown(self):
        os.remove(self.testfile1)
        os.remove(self.testfile2)
        os.remove("empty.txt")

class Tag(unittest.TestCase):
    def setUp(self):
        # create a fake BMTagger path
        self.testScript = "test.sh"
        with open(self.testScript, "w") as f:
            f.write("#!/bin/bash\necho $@")
        os.chmod(self.testScript, 0744)
        self.testFile1 = "test1.fastq"
        self.testFile2 = "test2.fastq"

    def testNone(self):
        ''' BMTagger is called with 0 databases'''
        self.assertEqual(kd.tag(
            infile = [self.testFile1], db_prefix = [], bmtagger_path =
            "./test.sh", single_end = True, prefix = "test",
            remove = True, debug = True, temp_dir = "."), ([], []))
    
    def testBadScript(self):
        '''Should return nonzero exit code'''
        self.assertNotEqual(kd.tag(
            infile = [self.testFile1], db_prefix = [], bmtagger_path =
            "ls", single_end = True, prefix = "test",
            remove = True, debug = True, temp_dir = ".")[0], [0])
        
    def testBuildsCorrectly(self):
        ''' Check the BMTagger call, provided correct inputs '''

        self.assertEqual(kd.tag(infile = [self.testFile1], db_prefix =
            ["test_db"], bmtagger_path = "./test.sh", single_end = True, prefix =
            "test", remove = False, debug = False, temp_dir = "temp"), ([0],
            ["./test.sh -q 1 -1 test1.fastq -b test_db.bitmask -x test_db.srprism -T temp -o test_db0.out"])) 
        self.assertEqual(kd.tag(infile = [self.testFile1], db_prefix =
            ["test_db"], bmtagger_path = "./test.sh", single_end = True, prefix =
            "test", remove = True, debug = False, temp_dir = "temp"), ([0],
            ["./test.sh -q 1 -1 test1.fastq -b test_db.bitmask -x test_db.srprism -T temp -o test_db0 --extract"])) 

        self.assertEqual(kd.tag(infile = [self.testFile1], db_prefix =
            ["test_db"], bmtagger_path = "./test.sh", single_end = True, prefix =
            "test", remove = True, debug = True, temp_dir = "temp"), ([0],
            ["./test.sh -q 1 -1 test1.fastq -b test_db.bitmask -x test_db.srprism -T temp -o test_db0 --extract --debug"])) 

        self.assertEqual(kd.tag(infile = [self.testFile1], db_prefix =
            ["test_db"], bmtagger_path = "./test.sh", single_end = True, prefix =
            "test", remove = False, debug = True, temp_dir = "temp"), ([0],
            ["./test.sh -q 1 -1 test1.fastq -b test_db.bitmask -x test_db.srprism -T temp -o test_db0.out --debug"])) 

        self.assertEqual(kd.tag(infile = [self.testFile1], db_prefix =
            ["test_db"], bmtagger_path = "./test.sh", single_end = True, prefix =
            "test", remove = False, debug = True, temp_dir = "temp", orphan=1), ([0],
            ["./test.sh -q 1 -1 test1.fastq -b test_db.bitmask -x test_db.srprism -T temp -o test_db0_se_1.out --debug"])) 

        self.assertEqual(kd.tag(infile = [self.testFile1], db_prefix =
            ["test_db"], bmtagger_path = "./test.sh", single_end = True, prefix =
            "test", remove = False, debug = True, temp_dir = "temp", orphan=2), ([0],
            ["./test.sh -q 1 -1 test1.fastq -b test_db.bitmask -x test_db.srprism -T temp -o test_db0_se_2.out --debug"])) 

        self.assertEqual(kd.tag(infile = [self.testFile1, self.testFile2],
            db_prefix = ["test_db"], bmtagger_path = "./test.sh", single_end =
            False, prefix = "test", remove = False, debug = True, temp_dir =
            "temp"), ([0], 
            ["./test.sh -q 1 -1 test1.fastq -2 test2.fastq -b test_db.bitmask -x test_db.srprism -T temp -o test_db0_pe.out --debug"])) 

        self.assertEqual(kd.tag(infile = [self.testFile1, self.testFile2],
            db_prefix = ["test_db"], bmtagger_path = "./test.sh", single_end =
            False, prefix = "test", remove = True, debug = True, temp_dir =
            "temp"), ([0], 
            ["./test.sh -q 1 -1 test1.fastq -2 test2.fastq -b test_db.bitmask -x test_db.srprism -T temp -o test_db0_pe --extract --debug"])) 

        self.assertEqual(kd.tag(infile = [self.testFile1, self.testFile2],
            db_prefix = ["test_db"], bmtagger_path = "./test.sh", single_end =
            False, prefix = "test", remove = False, debug = False, temp_dir =
            "temp"), ([0], 
            ["./test.sh -q 1 -1 test1.fastq -2 test2.fastq -b test_db.bitmask -x test_db.srprism -T temp -o test_db0_pe.out"])) 

    def testMultipleDB(self):
        self.assertEqual(kd.tag(infile = [self.testFile1, self.testFile2],
            db_prefix = ["test_db1", "test_db2"], bmtagger_path = "./test.sh",
            single_end = False, prefix = "test", remove = False, debug = False,
            temp_dir = "temp"), 
            ([0,0], ["./test.sh -q 1 -1 test1.fastq -2 test2.fastq -b test_db1.bitmask -x test_db1.srprism -T temp -o test_db0_pe.out", "./test.sh -q 1 -1 test1.fastq -2 test2.fastq -b test_db2.bitmask -x test_db2.srprism -T temp -o test_db1_pe.out"])) 
        self.assertEqual(kd.tag(infile = [self.testFile1, self.testFile2],
            db_prefix = ["test_db1", "test_db2"], bmtagger_path = "./test.sh",
            single_end = False, prefix = "test", remove = True, debug = False,
            temp_dir = "temp"), 
            ([0,0], ["./test.sh -q 1 -1 test1.fastq -2 test2.fastq -b test_db1.bitmask -x test_db1.srprism -T temp -o test_db0_pe --extract", "./test.sh -q 1 -1 test1.fastq -2 test2.fastq -b test_db2.bitmask -x test_db2.srprism -T temp -o test_db1_pe --extract"])) 


    def tearDown(self):
        os.remove(self.testScript)

class CheckFastq(unittest.TestCase):
    def tests(self):
        ''' True if file name ends in .fastq, else false '''
        self.assertTrue(kd.check_fastq("test.fastq"))
        self.assertFalse(kd.check_fastq("test.out"))
        self.assertFalse(kd.check_fastq(""))

class CombineTagBase(unittest.TestCase):
    def setUp(self):
        # create log file
        self.logfile = "logfile.log"
        with open(self.logfile, "w") as f:
            pass
        # create files to join
        self.db1_pe1_fastq = "pre_db1_pe_1.fastq"
        self.db1_pe2_fastq = "pre_db1_pe_2.fastq"
        self.db2_pe1_fastq = "pre_db2_pe_1.fastq"
        self.db2_pe2_fastq = "pre_db2_pe_2.fastq"
        self.db1_pe1_out = "pre_db1_pe.out"
        self.db2_pe1_out = "pre_db2_pe.out"
        self.db1_se1_fastq = "pre_db1_se_1.fastq"
        self.db2_se1_fastq = "pre_db2_se_1.fastq"
        self.files = [self.db1_pe1_fastq, self.db1_pe2_fastq,
                self.db2_pe1_fastq, self.db2_pe2_fastq, self.db1_pe1_out,
                self.db2_pe1_out, self.db1_se1_fastq, self.db2_se1_fastq]
        for fname in self.files:
            with open(fname, "w") as f:
                for i in range(4):
                    f.write("Have this in common\n")
                for i in range(4):
                    f.write(fname + "\n")

    '''
    def tearDown(self):
        os.remove("logfile.log")
        for fname in self.files:
            os.remove(fname)
    '''

class CombineTag(CombineTagBase, unittest.TestCase):
    def test_combine_tag(self):
        kd.combine_tag([[self.db1_pe1_fastq,
            self.db1_pe2_fastq],[self.db2_pe1_fastq, self.db2_pe2_fastq]],
            logfile = self.logfile, out_prefix = "pre", single_end=False)

class IntersectFastq(CombineTagBase, unittest.TestCase):
    def test_intersect_fastq(self):
        ''' aaaaaaa '''
        kd.intersect_fastq([self.db1_pe1_fastq,self.db2_pe1_fastq], "")

if __name__ == '__main__':
    unittest.main()
