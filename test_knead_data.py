import knead_data as kd
import unittest
import os

class testCheckFiles(unittest.TestCase):
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

        self.assertEqual(0, kd.checkfile(self.bad_file1), 
            "Passing in tuple in checkfile")
        self.assertEqual(kd.checkfile(self.bad_file2), 0,
                "Passing in int to checkfile")

    def test_checkexists(self):
        self.assertTrue(kd.checkexists(self.files[:2]),
                "passing existing files to checkexists")
        self.assertRaises(IOError, kd.checkexists, self.files[2],
                "passing nonexistent file to checkexists")
        self.assertRaises(TypeError, kd.checkexists, self.bad_file1,
                "passing tuple to checkexists")
        self.assertRaises(TypeError, kd.checkexists, self.bad_file2,
                "passing int to checkexists")

    def tearDown(self):
        os.remove("empty.txt")
        os.remove("non-empty.txt")


if __name__ == '__main__':
    unittest.main()
