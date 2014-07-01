===============================
CFSAN SNP Pipeline
===============================

.. Image showing the PyPI version badge - links to PyPI
.. .. image:: https://badge.fury.io/py/snp-pipeline.png
..    :target: http://badge.fury.io/py/snp-pipeline
    
.. Image showing the Travis Continuous Integration test status
.. .. image:: https://travis-ci.org/CFSAN-Biostatistics/snp-pipeline.png?branch=master
..        :target: https://travis-ci.org/CFSAN-Biostatistics/snp-pipeline

.. Image showing the PyPi download count
.. .. image:: https://pypip.in/d/snp-pipeline/badge.png
..        :target: https://pypi.python.org/pypi/snp-pipeline

The CFSAN SNP Pipeline is a Python-based system for the production of SNP 
matrices from sequence data used in the phylogenetic analysis of pathogenic 
organisms sequenced from samples of interest to food safety.

The SNP Pipeline was developed by the United States Food 
and Drug Administration, Center for Food Safety and Applied Nutrition.

* Free software: See license below. 
* Documentation: http://snp-pipeline.rtfd.org
* Source Code: https://github.com/CFSAN-Biostatistics/snp-pipeline
* PyPI Distribution: https://pypi.python.org/pypi/snp-pipeline

Features
--------

The CFSAN SNP Pipeline uses reference-based alignments to create a matrix of
SNPs for a given set of samples. The process generally starts off by finding
a reference that is appropriate for the samples of interest, and collecting
the sample sequence data into an appropriate directory structure. The SNP
pipeline can then be used to perform the alignment of the samples to the
reference. Once the sample sequences are aligned, a list of SNP positions is
generated. The list of SNP positions is then used in combination with
alignments of the samples to the reference sequence to call SNPs. The SNP
calls are organized into a matrix containing (just) the SNP calls for all
of the sequences.

The CFSAN SNP Pipeline is written in a combination of bash and python. The
code (including the bash scripts) is designed to be straighforward to
install. A script is provided that allows the running of the Python code
from the command line. All of the parameters used are available as command
line options. In situations where additional customization is desired, the
code is not highly complex and should be easy to modify as necessary.

Two examples of using the code are provided. These examples serve as both
unit tests, and as examples that can be modified to work on other data
sets of interest.

Citing SNP Pipeline
-------------------

To cite SNP Pipeline, please reference the SNP Pipeline GitHub repository:

    https://github.com/CFSAN-Biostatistics/snp-pipeline

and cite the associated paper:

    An evaluation of alternative methods for constructing
    phylogenies from whole genome sequence data: A case
    study with *Salmonella*


License
-------

This project constitutes a work of the United States Government and is not subject to domestic copyright protection under 17 USC § 105.

This program is free software: you can redistribute it and/or modify it under the terms of the included License.

This program is distributed in the hope that it will be useful. Responsibility
for the use of the system and interpretation of documentation and results lies
solely with the user. In no event shall CFSAN be liable for direct, indirect,
special, incidental, or consequential damages resulting from the use, misuse,
or inability to use the system and accompanying documentation. Third parties'
use of or acknowledgment of the sustem does not in any way represent that
CFSAN endorses such third parties or expresses any opinion with respect to
their statements. 
