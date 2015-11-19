#!/usr/bin/env python

import sys
if sys.version_info < (2,7,):
    import unittest2 as unittest
else:
    import unittest
import doctest

from snppipeline import pileup

def load_tests(loader, tests, pattern):
    """Run pileup module doctests"""
    tests.addTests(doctest.DocTestSuite(pileup))
    return tests

if __name__ == "__main__":
    unittest.main()
