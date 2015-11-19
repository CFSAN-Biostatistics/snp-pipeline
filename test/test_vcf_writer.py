#!/usr/bin/env python

import sys
if sys.version_info < (2,7,):
    import unittest2 as unittest
else:
    import unittest
import doctest

from snppipeline import vcf_writer

def load_tests(loader, tests, pattern):
    """Run vcf_writer module doctests"""
    tests.addTests(doctest.DocTestSuite(vcf_writer))
    return tests

if __name__ == "__main__":
    unittest.main()
