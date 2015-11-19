#!/usr/bin/env python

import sys
if sys.version_info < (2,7,):
    import unittest2 as unittest
else:
    import unittest
import doctest

from snppipeline import utils

def load_tests(loader, tests, pattern):
    """Run utils doctests"""
    tests.addTests(doctest.DocTestSuite(utils))
    return tests

if __name__ == "__main__":
    unittest.main()
