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
import random
import collections
import itertools
import shutil

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



'''
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
'''
        


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

class PrefixBases(unittest.TestCase):
    def testNoRepeats(self):
        l = ['/home/a.db', '/dir/b.db']
        true_ans = [('a.db', '/home/a.db'), ('b.db', '/dir/b.db')]
        ans = kd._prefix_bases(l)
        for (t, a) in zip(true_ans, ans):
            self.assertEqual(t, a)

class newDir(unittest.TestCase):
    '''Makes a new temporary directory to put temp files. Nukes this directory
    when finished. '''

    def setUp(self):
        # save current directory, make new directory for temporary files
        self.curwd = os.getcwd()
        self.newdir = "temp"
        try:
            os.mkdir(os.path.join(self.curwd, self.newdir))
        except OSError:
            if not os.path.isdir(os.path.join(self.curwd, self.newdir)):
                raise

        os.chdir(os.path.join(self.curwd, self.newdir))

    def tearDown(self):
        os.chdir(self.curwd)
        shutil.rmtree(os.path.join(self.curwd,self.newdir))


class Tag(newDir, unittest.TestCase):
    def setUp(self):

        newDir.setUp(self)

        # create a fake BMTagger path
        self.testScript = "test.sh"

        with open(self.testScript, "w") as f:
            f.write("#!/bin/bash\necho $@")

        os.chmod(self.testScript, 0744)
        self.testFile1 = "test1.fastq"
        self.testFile2 = "test2.fastq"

        self.logfile = "test.log"
        with open(self.logfile, "w") as f:
            pass
        self.fakeOutputs = ["test_db_contam.out",
                            "test_db_clean.fastq",
                            "test_db_clean_1.fastq",
                            "test_db_clean_2.fastq",
                            "test_db1_contam.out",
                            "test_db2_contam.out",
                            "test_db1_clean_1.fastq",
                            "test_db1_clean_2.fastq",
                            "test_db2_clean_1.fastq",
                            "test_db2_clean_2.fastq"]

        for fakeOutput in self.fakeOutputs:
            with open(fakeOutput, "w") as f:
                f.write("Read 1\nGATTACA\n+\nKTACRIA\n")


    def testNone(self):
        ''' BMTagger is called with 0 databases'''
        self.assertEqual(kd.tag(
                infile_list = [self.testFile1], db_prefix_list = [], 
                bmtagger_path = "./test.sh", prefix = "test", remove = True, 
                debug = True, temp_dir = ".", logfile = self.logfile), ([], []))
    
    def testBadScript(self):
        '''Should return nonzero exit code'''
        self.assertNotEqual(kd.tag( infile_list = [self.testFile1],
                db_prefix_list = ['a.db'], bmtagger_path = "ls", prefix = "test",
                remove = True, debug = True, temp_dir = ".", 
                logfile = self.logfile)[0], [0])
        
    def testBuildsCorrectly(self):
        ''' Check the BMTagger call, provided correct inputs '''

        # be very careful of the order of these tests. Run the ones that have
        # the --debug flag set first, because otherwise you will remove the
        # *_db[number]* files and the subsequent tests that rely on them will
        # fail. 

        # debug on/off, can correctly choose between clean vs. contam outputs
        self.assertEqual(kd.tag(infile_list = [self.testFile1], db_prefix_list =
            ["db"], bmtagger_path = "./test.sh", prefix =
            "test", remove = False, debug = True, temp_dir = "temp", logfile =
            self.logfile), ([0],
            [shlex.split("./test.sh -q 1 -1 test1.fastq -b db.bitmask -x db.srprism -T temp -o test_db_contam.out --debug")])) 

        self.assertEqual(kd.tag(infile_list = [self.testFile1], db_prefix_list =
            ["db"], bmtagger_path = "./test.sh", prefix =
            "test", remove = True, debug = True, temp_dir = "temp", logfile =
            self.logfile), ([0],
            [shlex.split("./test.sh -q 1 -1 test1.fastq -b db.bitmask -x db.srprism -T temp -o test_db_clean --extract --debug")]))

        # simulate paired end inputs
        self.assertEqual(kd.tag(infile_list = [self.testFile1, self.testFile2],
            db_prefix_list = ["db"], bmtagger_path = "./test.sh", prefix =
            "test", remove = False, debug = True, temp_dir = "temp", logfile =
            self.logfile), ([0], 
            [shlex.split("./test.sh -q 1 -1 test1.fastq -2 test2.fastq -b db.bitmask -x db.srprism -T temp -o test_db_contam.out --debug")])) 

        self.assertEqual(kd.tag(infile_list = [self.testFile1, self.testFile2],
            db_prefix_list = ["db"], bmtagger_path = "./test.sh", prefix =
            "test", remove = True, debug = True, temp_dir = "temp", logfile =
            self.logfile), ([0], 
            [shlex.split("./test.sh -q 1 -1 test1.fastq -2 test2.fastq -b db.bitmask -x db.srprism -T temp -o test_db_clean --extract --debug")])) 

        # testing with debug off
        self.assertEqual(kd.tag(infile_list = [self.testFile1], db_prefix_list =
            ["db"], bmtagger_path = "./test.sh", prefix =
            "test", remove = False, debug = False, temp_dir = "temp", logfile =
            self.logfile), ([0],
            [shlex.split("./test.sh -q 1 -1 test1.fastq -b db.bitmask -x db.srprism -T temp -o test_db_contam.out")])) 

        self.assertEqual(kd.tag(infile_list = [self.testFile1], db_prefix_list =
            ["db"], bmtagger_path = "./test.sh", prefix =
            "test", remove = True, debug = False, temp_dir = "temp", logfile =
            self.logfile), ([0],
            [shlex.split("./test.sh -q 1 -1 test1.fastq -b db.bitmask -x db.srprism -T temp -o test_db_clean --extract")])) 



    def testMultipleDBContam(self):
        ''' test output of contaminated reads for multiple DBs '''
        true_output = map(shlex.split, ["./test.sh -q 1 -1 test1.fastq -2 test2.fastq -b db1.bitmask -x db1.srprism -T temp -o test_db1_contam.out --debug", "./test.sh -q 1 -1 test1.fastq -2 test2.fastq -b db2.bitmask -x db2.srprism -T temp -o test_db2_contam.out --debug"])

        output = kd.tag(infile_list = [self.testFile1, self.testFile2],
                db_prefix_list = ["db1", "db2"], bmtagger_path = "./test.sh", 
                prefix = "test", remove = False, debug = True, 
                temp_dir = "temp", logfile = self.logfile) 
        self.assertEqual(output[0], [0,0])
        self.assertTrue(output[1][0] in true_output)
        self.assertTrue(output[1][1] in true_output)

    def testMultipleDBClean(self):
        ''' test output of clean reads for multiple DBs '''
        true_output = map(shlex.split, ["./test.sh -q 1 -1 test1.fastq -2 test2.fastq -b db1.bitmask -x db1.srprism -T temp -o test_db1_clean --extract --debug", "./test.sh -q 1 -1 test1.fastq -2 test2.fastq -b db2.bitmask -x db2.srprism -T temp -o test_db2_clean --extract --debug"])

        output = kd.tag(infile_list = [self.testFile1, self.testFile2],
                db_prefix_list = ["db1", "db2"], 
                bmtagger_path = "./test.sh", prefix = "test", remove = True, 
                debug = True, temp_dir = "temp", logfile = self.logfile)
        self.assertEqual(output[0], [0,0])
        self.assertTrue(output[1][0] in true_output)
        self.assertTrue(output[1][1] in true_output)

    def tearDown(self):
        os.remove(self.testScript)
        os.remove(self.logfile)

        newDir.tearDown(self)


class CheckFastq(unittest.TestCase):
    def tests(self):
        ''' True if file name ends in .fastq, else false '''
        self.assertTrue(kd.check_fastq("test.fastq"))
        self.assertFalse(kd.check_fastq("test.out"))
        self.assertFalse(kd.check_fastq(""))

class CombineTagBase(newDir):
    def setUp(self):

        newDir.setUp(self)

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
        self.fastqs = [self.db1_pe1_fastq, self.db1_pe2_fastq,
                self.db2_pe1_fastq, self.db2_pe2_fastq, self.db1_se1_fastq,
                self.db2_se1_fastq]
        self.outs = [self.db1_pe1_out, self.db2_pe1_out]


        self.iterations = 10
        
        self.common_fastq = 0
        self.common_out = 0
        self.notcommon_fastq = 0
        self.notcommon_out = 0
        
        self.sequences = ["Read in common1\nsequence1\n+\nscore1", 
                        "Read in common2\nsequence2\n+\nscore2",
                        "Read in common3\nsequence3\n+\nscore3"]
        self.headers = ["Read header in common1",
                        "Read header in common2",
                        "Read header in common3"]

        for fname in self.fastqs:
            with open(fname, "w") as f:
                for i in range(self.iterations):
                    floatRand = random.random()
                    if floatRand < 0.25:
                        f.write(self.sequences[0]+"\n")
                        self.common_fastq += 1
                    elif floatRand < 0.5:
                        f.write(self.sequences[1]+"\n")
                        self.common_fastq += 1
                    elif floatRand < 0.75:
                        f.write(self.sequences[2]+"\n")
                        self.common_fastq += 1
                    else:
                        f.write(fname + str(self.notcommon_fastq) + 
                                "\nsequence\n+\nscore\n")
                        self.notcommon_fastq += 1

        for fname in self.outs:
            with open(fname, "w") as f:
                for i in range(self.iterations):
                    floatRand = random.random()
                    if floatRand < 0.25:
                        f.write(self.headers[0]+"\n")
                        self.common_out += 1
                    if floatRand < 0.5:
                        f.write(self.headers[1]+"\n")
                        self.common_out += 1
                    if floatRand < 0.75:
                        f.write(self.headers[2]+"\n")
                        self.common_out += 1
                    else:
                        f.write(fname + str(self.notcommon_out) + "\n")
                        self.notcommon_out += 1

    def tearDown(self):
        newDir.tearDown(self)

# Inherit the base class that does setUp and tearDown BEFORE inheriting
# unittest.TestCase
class IntersectFastq(CombineTagBase, unittest.TestCase):
    def test_intersect_fastq(self):
        ''' Testing intersection for fastqs''' 
        kd.intersect_fastq([self.db1_se1_fastq, self.db2_se1_fastq], "pre_se.fastq")
        lines = []
        counter = collections.Counter()
        with open("pre_se.fastq", "rU") as f:
            for lines in itertools.izip_longest(*[f]*4):
                read = ("".join(lines)).strip()
                counter[read] += 1
                self.assertTrue(read in self.sequences)
        os.remove("pre_se.fastq")
        for key in counter:
            self.assertTrue(key in self.sequences)
            self.assertEqual(counter[key], 1)


class UnionOutfiles(CombineTagBase, unittest.TestCase):
    def test_union_outfiles(self):
        ''' Testing unions for .out files'''
        kd.union_outfiles(self.outs, "pre.out")
        lines = []
        counter = collections.Counter()
        with open("pre.out", "rU") as f:
            for line in f:
                counter[line.strip()] += 1
        for key in counter:
            # test that every noncommon line appears once
            self.assertEqual(counter[key], 1)

            # test that the keys of the counter are legit
            if key[:4] != "Read":
                self.assertTrue(key[:-1] in self.outs)
            else:
                self.assertTrue(key in self.headers)
        os.remove("pre.out")


class CombineTag(CombineTagBase, unittest.TestCase):
    def setUp(self):
        CombineTagBase.setUp(self)
        # create log file
        self.logfile = "logfile.log"
        with open(self.logfile, "w") as f:
            pass

    # since we tested union and intersect before, we just need to check that we
    # got the output filenames and the logging correct.
    def test_combine_tag_outputs(self):
        # good paired end input 
        outputs = kd.combine_tag([[self.db1_pe1_fastq,
            self.db1_pe2_fastq],[self.db2_pe1_fastq, self.db2_pe2_fastq]],
            logfile = self.logfile, out_prefix = "pre")
        self.assertEqual(outputs, ["pre_1.fastq", "pre_2.fastq"])

        # bad outfile input. can't have >1 outfile per database
        self.assertRaises(Exception, kd.combine_tag, llstrFiles =
                [[self.db1_pe1_out, self.db1_pe1_out],[self.db2_pe1_out,
                    self.db2_pe1_out]], logfile = self.logfile, out_prefix =
                "pre")

        # aardvark figure out what's happening here
        # bad fastq input
        self.assertRaises(Exception, kd.combine_tag, llstrFiles =
                [[self.db1_pe1_fastq, self.db1_pe2_fastq],
                    [self.db2_pe1_fastq]], logfile = self.logfile, out_prefix =
                "pre")

        outputs = kd.combine_tag([[self.db1_pe1_out],[self.db2_pe1_out]],
                logfile = self.logfile, out_prefix = "pre")
        self.assertEqual(outputs, ["pre.out"])

        outputs = kd.combine_tag([[self.db1_pe1_out],[self.db2_pe1_out]],
                logfile = self.logfile, out_prefix = "pre")
        self.assertEqual(outputs, ["pre.out"])

        outputs = kd.combine_tag([[self.db1_se1_fastq],[self.db2_se1_fastq]],
                logfile = self.logfile, out_prefix = "pre")
        self.assertEqual(outputs, ["pre.fastq"])

    def test_combine_tag_logfile(self):
        kd.combine_tag([[self.db1_pe1_fastq,
            self.db1_pe2_fastq],[self.db2_pe1_fastq, self.db2_pe2_fastq]],
            logfile = self.logfile, out_prefix = "pre")

        file_contents = None
        with open(self.logfile, "r") as f:
            file_contents = f.read()
            self.assertTrue(len(file_contents) > 0)

    def tearDown(self):
        CombineTagBase.tearDown(self)

class MsgNumReads(CombineTagBase, unittest.TestCase):
    def test_msg_num_reads(self):
        msg = kd.msg_num_reads(self.files[:4])
        true_msg = "pre_db1_pe_1.fastq: 10\npre_db1_pe_2.fastq: 10\npre_db2_pe_1.fastq: 10\npre_db2_pe_2.fastq: 10" 
        self.assertEqual(msg, true_msg)


if __name__ == '__main__':
    unittest.main()
