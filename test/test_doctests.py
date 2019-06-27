#!/usr/bin/env python

"""
This module runs doctests unit tests on all modules with doctests.
"""

from __future__ import print_function

import doctest
import glob
import os
import sys
if sys.version_info < (2, 7,):
    import unittest2 as unittest
else:
    import unittest

import snppipeline as source_package


def load_module_by_path(path):
    """Load a python module from its path.

    Parameters
    ----------
    path : str
        Path to the module source file.

    Returns
    -------
    mod : module
        Loaded module.
    """
    import imp

    module_file_basename = os.path.basename(path)
    module_name, ext = os.path.splitext(module_file_basename)
    mod = imp.load_source(module_name, path)
    return mod


def file_contains_doctests(path):
    """Scan a python source file to determine if it contains any doctest examples.

    Parameters
    ----------
    path : str
        Path to the module source file.

    Returns
    -------
    flag : bool
        True if the module source code contains doctest examples.
    """
    with open(path) as f:
        for line in f:
            if ">>>" in line:
                return True
    return False


def load_tests(loader, tests, pattern):
    """Run doctests for all modules"""
    source_dir = os.path.dirname(source_package.__path__[0])
    python_source_glob = os.path.join(source_dir, source_package.__name__, "*.py")
    python_source_files = glob.glob(python_source_glob)
    for python_source_file in python_source_files:
        if not file_contains_doctests(python_source_file):
            continue
        module = load_module_by_path(python_source_file)
        tests.addTests(doctest.DocTestSuite(module))
    return tests


if __name__ == "__main__":
    unittest.main()
