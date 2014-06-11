#!/usr/bin/env python2.7

import unittest
import doctest

from snppipeline import utils

class Test(unittest.TestCase):
    """Unit tests for utils."""

    def test_doctests(self):
        """Run utils doctests"""
        doctest.testmod(utils, verbose=True)

if __name__ == "__main__":
    unittest.main()
