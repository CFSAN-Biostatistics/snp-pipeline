#!/usr/bin/env python

"""
This module runs doctests unit tests on all modules with doctests.
"""

import sys
if sys.version_info < (2,7,):
    import unittest2 as unittest
else:
    import unittest
import doctest

from snppipeline import command
from snppipeline import fastq
from snppipeline import pileup
from snppipeline import utils
from snppipeline import vcf_writer

def load_tests(loader, tests, pattern):
    """Run doctests for all modules"""
    tests.addTests(doctest.DocTestSuite(command))
    tests.addTests(doctest.DocTestSuite(fastq))
    tests.addTests(doctest.DocTestSuite(pileup))
    tests.addTests(doctest.DocTestSuite(utils))
    tests.addTests(doctest.DocTestSuite(vcf_writer))
    return tests

if __name__ == "__main__":
    unittest.main()
