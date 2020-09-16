#!/usr/bin/env python

"""
This software is used to test the KneadData pipeline.
Dependencies: KneadData (and all KneadData dependencies)
"""

import os
import sys
import unittest

# Try to load the kneaddata package to check the installation
try:
    from kneaddata import utilities
    from kneaddata import config
except ImportError:
    sys.exit("CRITICAL ERROR: Unable to find the kneaddata python package." +
        " Please check your install.") 

# required python versions (2.7+ or 3.0+)
required_python_version_major = [2,3]
required_python_version_minor = [7,0]

# check for either of the required versions
pass_check=False
try:
    for major, minor in zip(required_python_version_major, required_python_version_minor):
        if (sys.version_info[0] == major and sys.version_info[1] >= minor):
            pass_check=True
except (AttributeError,IndexError):
    sys.exit("CRITICAL ERROR: The python version found (version 1) " +
        "does not match the version required (version "+
        str(required_python_version_major)+"."+
        str(required_python_version_minor)+"+)")

if not pass_check:
    sys.exit("CRITICAL ERROR: The python version found (version "+
        str(sys.version_info[0])+"."+str(sys.version_info[1])+") "+
        "does not match the version required (version "+
        str(required_python_version_major)+"."+
        str(required_python_version_minor)+"+)")

def check_dependency(exe,bypass_permissions_check=None):
    """ Check the dependency can be found """
    if not utilities.find_exe_in_path(exe, bypass_permissions_check):
        print("Warning: Unable to find "+exe+". Tests requiring this tool will be skipped.")

import argparse

def parse_arguments(args):
    """
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "KneadData Test\n",
        formatter_class=argparse.RawTextHelpFormatter,
        prog="kneaddata_test")
    parser.add_argument(
        "--bypass-functional-tests", 
        help="do not run the kneaddata end to end functional tests\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "--bypass-unit-tests", 
        help="do not run the unit tests\n", 
        action="store_true",
        default=False)

    return parser.parse_args()

def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    # Check for dependencies
    if not args.bypass_functional_tests:
        check_dependency(config.trimmomatic_jar,bypass_permissions_check=True)
        check_dependency(config.bowtie2_exe)
        check_dependency(config.trf_exe)
        check_dependency(config.samtools_exe)
        check_dependency(config.fastqc_exe)
    
    # Get the unittests
    directory_of_tests=os.path.dirname(os.path.abspath(__file__))
    
    tests_to_run=[]
    if not args.bypass_functional_tests:
        tests_to_run.append(unittest.TestLoader().discover(directory_of_tests,pattern='functional_tests*.py'))
    if not args.bypass_unit_tests:
        tests_to_run.append(unittest.TestLoader().discover(directory_of_tests,pattern='basic_tests*.py'))

    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(tests_to_run))

if __name__ == '__main__':
    main()
