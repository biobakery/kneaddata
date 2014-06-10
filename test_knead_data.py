import knead_data as kd
import unittest
import os

class test_checkfile(unittest.TestCase):
    def setUp(self): 
        with open("empty.txt", "w") as f:
            pass
        with open("non-empty.txt", "w") as f:
            f.write("Hello world!")
        self.files = ["empty.txt", "non-empty.txt", "dne.txt"]
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
        os.remove("empty.txt")
        os.remove("non-empty.txt")


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


if __name__ == '__main__':
    unittest.main()
