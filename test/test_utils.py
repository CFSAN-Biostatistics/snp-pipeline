#!/usr/local/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 08:42:02 2014

@author: hugh.rand
"""

# -*- coding: utf-8 -*-
import unittest
import doctest

import googlemaps

class Test(unittest.TestCase):
    """Unit tests for googlemaps."""

    def test_doctests(self):
        """Run googlemaps doctests"""
        doctest.testmod(googlemaps)

if __name__ == "__main__":
    unittest.main()