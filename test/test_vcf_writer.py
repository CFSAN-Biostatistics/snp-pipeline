#!/usr/bin/env python

import unittest
import doctest

from snppipeline import vcf_writer

class Test(unittest.TestCase):
    """Unit tests for vcf_writer module."""

    def test_doctests(self):
        """Run vcf_writer module doctests"""
        doctest.testmod(vcf_writer, verbose=False)

if __name__ == "__main__":
    unittest.main()
