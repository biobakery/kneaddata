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
    import kneaddata
except ImportError:
    sys.exit("CRITICAL ERROR: Unable to find the kneaddata python package." +
        " Please check your install.") 

# check for the required python version
required_python_version_major = 2
required_python_version_minor = 7
    
try:
    if (sys.version_info[0] != required_python_version_major or
        sys.version_info[1] < required_python_version_minor):
        sys.exit("CRITICAL ERROR: The python version found (version "+
            str(sys.version_info[0])+"."+str(sys.version_info[1])+") "+
            "does not match the version required (version "+
            str(required_python_version_major)+"."+
            str(required_python_version_minor)+"+)")
except (AttributeError,IndexError):
    sys.exit("CRITICAL ERROR: The python version found (version 1) " +
        "does not match the version required (version "+
        str(required_python_version_major)+"."+
        str(required_python_version_minor)+"+)")  

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
        "-r","--run-functional-tests", 
        help="also run the functional tests\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "-b","--bypass-unit-tests", 
        help="do not run the unit tests\n", 
        action="store_true",
        default=False)
    
    return parser.parse_args()

def get_testdirectory():
    """ Return the location of all of the tests """
    
    return os.path.dirname(os.path.abspath(__file__))

def get_funtionaltests():
    """ Return all of the functional tests """
    
    directory_of_tests=get_testdirectory()
    
    functional_suite = unittest.TestLoader().discover(directory_of_tests,pattern='functional_tests*.py')
    
    return functional_suite

def get_unittests():
    """ Return all of the unit tests """
    
    directory_of_tests=get_testdirectory()
    
    basic_suite = unittest.TestLoader().discover(directory_of_tests,pattern='basic_tests*.py')
    advanced_suite = unittest.TestLoader().discover(directory_of_tests, pattern='advanced_tests*.py')    
    
    return [basic_suite, advanced_suite]

def unittests_suite_only():
    """ Return a TestSuite of just the unit tests """
    
    return unittest.TestSuite(get_unittests())

def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    # Get the unittests
    test_suites=[]
    if not args.bypass_unit_tests:
        test_suites=get_unittests()
    
    # Get the functional tests if requested
    if args.run_functional_tests:
        test_suites+=get_funtionaltests()

    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(test_suites))