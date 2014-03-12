#!/usr/bin/env python2.7

import unittest
import doctest

import imp
utilsnew = imp.load_source('utilsnew', '/home/hugh.rand/projects/snppipeline/snppipeline/utilsnew.py')

class Test(unittest.TestCase):
    """Unit tests for utilsnew."""

    def test_doctests(self):
        """Run utilsnew doctests"""
        doctest.testmod(utilsnew)

if __name__ == "__main__":
    unittest.main()
    
