#!/usr/bin/env python2.7

import unittest
import doctest

import imp
utils = imp.load_source('utils', '/home/hugh.rand/mnt/biob/svn/Biostats/rand/snppipeline/snppipeline/utils.py')

class Test(unittest.TestCase):
    """Unit tests for utils."""

    def test_doctests(self):
        """Run utils doctests"""
        doctest.testmod(utils)

if __name__ == "__main__":
    unittest.main()
    
