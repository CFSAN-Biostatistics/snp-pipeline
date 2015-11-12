#!/usr/bin/env python

import unittest
import doctest

from snppipeline import pileup

class Test(unittest.TestCase):
    """Unit tests for pileup module."""

    def test_doctests(self):
        """Run pileup module doctests"""
        doctest.testmod(pileup, verbose=False)

if __name__ == "__main__":
    unittest.main()
